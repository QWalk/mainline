/*

An attempt on a general tool to process log files (without need to
know what quantities we are actually dealing with). The only thing
we need to be present in a block is 'totweight'. Right now, the
whole thing is geared towards the local moments and charges, but it
might have a potential to become more general.

As of now
 - 'label' is not respected (we do not distinguish between
   VMC, DMC, RMC & RMC_cen ...)
 - does not handle array valued quantities
 - no dropping of thermalization blocks (difficult to do automatically,
   as different quantities are likely to have different thermalization
   time); easy fix would be a manual switch with a number of first
   blocks to drop 

*/

#include "Qmc_std.h"
#include "qmc_io.h"

int main(int argc, char **argv) {

  cout.precision(10);

  int reblock=1;
  string filename;
  if(!strcmp(argv[1], "-reblock") && argc == 4) {
      reblock=atoi(argv[2]);
      filename=argv[3];
  } else if ( argc==2 ) {
    filename=argv[1];
  } else {
    cout << "Usage: locmom [-reblock <number>] filename" << endl;
  }
    
  // read log file into words
  ifstream input(filename.c_str());
  vector < string >  words;
  parsefile(input, words);
  input.close();
  
  // divide words into blocks
  vector < vector < string> > blocktext;
  vector <string> tmp;
  unsigned int pos=0;
  while(readsection(words, pos, tmp, "BLOCK"))
    blocktext.push_back(tmp);
  int num_blocks=blocktext.size();

  // extract the section names (i.e., quantities) from the first block
  vector <string> quantity;
  int level=0;
  for(pos=0; pos< blocktext[0].size(); pos++) {
    if(blocktext[0][pos]==startsec) {
      level++;
    }
    else if(blocktext[0][pos]==endsec) {
      level--;
    }
    // for now I don't consider array-valued quantities
    if(level==0 && blocktext[0][pos]!=endsec && blocktext[0][pos+1]==startsec
       && blocktext[0][pos]!="z_pol" && blocktext[0][pos]!="autocorr_energy") {
      //cout << blocktext[0][pos] << endl;
      quantity.push_back(blocktext[0][pos]);
    }
  }
  int num_quantities=quantity.size();

  // extract data from blocks, last element of raw_block_data is totweight
  Array2 <doublevar> raw_block_data(num_blocks,num_quantities+1);
  for (int b=0; b<num_blocks; b++) {
    string value;
    for (int j=0; j<num_quantities; j++) {
      vector<string> section;
      //readvalue(blocktext[b],pos=0,value,quantity[j].c_str());
      readsection(blocktext[b],pos=0,section,quantity[j].c_str());
      raw_block_data(b,j)=atof(section[0].c_str());
    }
    readvalue(blocktext[b],pos=0,value,"totweight");
    raw_block_data(b,num_quantities)=atof(value.c_str());
  }

  // here will be reblocking, for now just copy raw_block_data
  /*
  int num_reblocks=num_blocks;
  Array2 <doublevar> reblock_data(num_reblocks,num_quantities+1);
  for (int b=0; b<num_reblocks; b++) {
    for (int j=0; j<num_quantities+1; j++) {
      reblock_data(b,j)=raw_block_data(b,j);
    }
  }
  */
  
  // reblocking
  int num_reblocks=num_blocks/reblock;
  cout << num_blocks << " total blocks reblocked into " << num_reblocks << endl;
  int reblock_start=num_blocks-reblock*num_reblocks;
  if ( reblock_start > 0 ) {
    cout << "First " << reblock_start
	 << " original blocks dropped (not enough to build another block)." << endl;
  }
  Array2 <doublevar> reblock_data(num_reblocks,num_quantities+1);
  reblock_data=0.0;
  while(reblock_start<num_blocks) {
    doublevar sum_weights=0.0;
    for (int b=reblock_start; b<reblock_start+reblock; b++)
      sum_weights+=raw_block_data(b,num_quantities);
    int nb=reblock_start/reblock;
    for (int b=reblock_start; b<reblock_start+reblock; b++) {
      for (int j=0; j<num_quantities; j++)
	reblock_data(nb,j)+=raw_block_data(b,j)*raw_block_data(b,num_quantities)/sum_weights;
      reblock_data(nb,num_quantities)=sum_weights;
    }
    reblock_start+=reblock;
  }

  // calculate averages and statistical errors
  Array1 <doublevar> average(num_quantities);
  average=0.0;
  Array1 <doublevar> error(num_quantities);
  error=0.0;
  doublevar sum_weights=0.0;
  for (int b=0; b<num_reblocks; b++) { 
    sum_weights+=reblock_data(b,num_quantities);
  }
  for (int j=0; j<num_quantities; j++) {
    for (int b=0; b<num_reblocks; b++) {
      average(j)+=reblock_data(b,j)*reblock_data(b,num_quantities)/sum_weights;
    }
    for (int b=0; b<num_reblocks; b++) {
      error(j)+=(reblock_data(b,j)-average(j))*(reblock_data(b,j)-average(j))
	*reblock_data(b,num_quantities)/sum_weights/(num_reblocks);
      // Lucas uses just /num_reblocks instead of /(num_reblocks-1)
    }
    // print out
    int field=cout.precision()+5;
    cout << setw(20)    << left  << quantity[j]
	 << setw(field) << right << average(j) << " +/- "
	 << setw(field) << sqrt(error(j)) << endl;
  } 

}

  
