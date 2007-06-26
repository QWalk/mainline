/*
 
Copyright (C) 2007 Lucas K. Wagner

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License along
with this program; if not, write to the Free Software Foundation, Inc.,
51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 
*/
#include "Qmc_std.h"
#include "Properties_block.h"
#include "qmc_io.h"
#include "Properties_average.h"
#include <algorithm>
#include <deque>

using namespace Properties_types;
//----------------------------------------------------------------------

void getBlocks(vector <string> & words, vector <string> & labels, 
               Array1 < Array1 < Properties_block> > & allblocks) {

  vector < vector < string> > blocktext;
  vector <string> tmp;
  unsigned int pos=0;
  while(readsection(words, pos, tmp, "BLOCK"))
    blocktext.push_back(tmp);

  int nblock=blocktext.size();
  //cout << blocktext.size() << " blocks \n";

  //vector <string> labels;
  vector <int> belong_to;
  vector <int> lab_nb; //number of blocks that belong to each label
  string label;
  for(int b=0; b< nblock; b++) {

    if(!readvalue(blocktext[b], pos=0, label, "LABEL"))
      error("didn't find label in the file");
    vector<string>::iterator place=find(labels.begin(), labels.end(), label);
    belong_to.push_back(place-labels.begin());
    if(place==labels.end()) {
      labels.push_back(label);
      lab_nb.push_back(1);
    }
    else lab_nb[belong_to[b]]++;
  }
  
  int nlabels=labels.size();
  //cout << nlabels << " labels" << endl;
  //for(int i=0; i< nlabels; i++) {
  //  cout << labels[i] << "  " << lab_nb[i] << endl;
  //}

  allblocks.Resize(nlabels);
  for(int i=0; i< nlabels; i++) {
    allblocks(i).Resize(lab_nb[i]);
  }
  Array1 <int> nb_read(nlabels, 0);
  for(int b=0; b< nblock; b++) {
    int lab=belong_to[b];
    int i=nb_read(lab);
    nb_read(lab)++;
    allblocks(lab)(i).restoreFromLog(blocktext[b]);
  }
}

//----------------------------------------------------------------------

void output_trace(vector <string> & labels, 
                  Array1 < Array1 <Properties_block> > & allblocks) {
  int nlabels=labels.size();
  for(int i=0; i< nlabels; i++) {
    string nm=labels[i]+".trace";
    ofstream tr(nm.c_str());
    tr.precision(15);
    for(int b=0; b< allblocks(i).GetDim(0); b++) {
      tr << allblocks(i)(b).avg(total_energy,0) << endl;
    }
  }
}
                 
//----------------------------------------------------------------------

void output_trace_force(vector <string> & labels,
                        Array1 < Array1 <Properties_block> > & allblocks) {
  int nlabels=labels.size();

  for(int i=0; i< nlabels; i++) {
    int naux=allblocks(i)(0).aux_energy.GetDim(0);
    int n_cvg=allblocks(i)(0).aux_energy.GetDim(1);
    for(int a=0; a< naux; a++) {
      for(int n=0; n< n_cvg; n++) {
        string nm=labels[i]+"f";
        append_number(nm,a);
        nm+="-"; append_number(nm,n);
        nm+=".trace";
        ofstream tr(nm.c_str());
        tr.precision(15);
        for(int b=0; b< allblocks(i).GetDim(0); b++) {
          tr << (allblocks(i)(b).aux_energy(a,n)-allblocks(i)(b).avg(total_energy,0))
          /allblocks(i)(b).aux_size(a) << endl;
      }
      }
    }
  }
  
}
//----------------------------------------------------------------------

int reblock_average(Array1 <Properties_block> & orig_block, int reblock, int equil,
            Properties_final_average & avg) {
  int nblock=orig_block.GetDim(0);
  if(reblock > nblock) 
    return 0;
  
  int neffblock=0;
  Array1<Properties_block> reblocks;

  reblocks.Resize(nblock/reblock);
  int start=nblock%reblock;
  int count=0;
  for(int i=start; i< nblock; i+=reblock) {
    int top=min(i+reblock, nblock);
    reblocks(count++).reduceBlocks(orig_block, i, top);
  }
  neffblock=nblock/reblock;
  avg.blockReduce(reblocks, 0, nblock/reblock, equil);

  return neffblock;
}

//----------------------------------------------------------------------

struct Gosling_options {
  string label;
  int tot_energy;
  int trace;
  int equil;
  int show_autocorr;
  int two_pt_forces;
  int trace_force;
  int reblock;
  vector<string> filenames;
  Gosling_options() { 
    tot_energy=0;
    trace=0;
    equil=1;
    show_autocorr=0;
    two_pt_forces=0;
    trace_force=0;
    reblock=1;
  }
};

//----------------------------------------------------------------------

void get_options(int argc, char ** argv, 
		 Gosling_options & options) { 


  for(int i=1; i< argc; i++) {
    if(!strcmp(argv[i], "-label") && argc > i+1) {
      options.label=argv[++i];
    }
    else if(!strcmp(argv[i], "-tot_energy"))
      options.tot_energy=1;
    else if(!strcmp(argv[i], "-trace"))
      options.trace=1;
    else if(!strcmp(argv[i], "-no_equil"))
      options.equil=0;
    else if(!strcmp(argv[i], "-show_autocorr"))
      options.show_autocorr=1;
    else if(!strcmp(argv[i], "-two_point_forces"))
      options.two_pt_forces=1;
    else if(!strcmp(argv[i], "-trace_force"))
      options.trace_force=1;
    else if(!strcmp(argv[i], "-reblock") && argc > i+1)
      options.reblock=atoi(argv[++i]);
    else if(!strcmp(argv[i], "-h")) {
      cout << "usage: gosling <options> <log file>" << endl
           << "-label <label>    : only print information for the given label" << endl
           << "-trace            : print out energy traces for each label, into files named label.trace" << endl
           << "-trace_force      : print out displacement traces for each label, into files named labelf#.trace" << endl
	   << "-reblock <number> : reblock the blocks into groups of <number>.  If they don't evenly divide, gosling" << endl
	   << "                    will throw out the first blocks until they do." << endl
           << "-no_equil         : don't try to find the equilibration steps; just average over everything" << endl
           << "-show_autocorr    : show estimates for the autocorrelation as a function of step for each label" << endl
           << "-two_point_forces : attempt to collate plus and minus differences into forces and print them out" << endl
           << "                    this will only work if the forces are in order + force, - force ..." << endl
           << "-h                : print this help message " << endl; 
    }
    else if(argv[i][0]=='-') 
      error("Unknown option: ",argv[i]);
    else 
      options.filenames.push_back(argv[i]);
  }

  //options.filenames.push_back(string(argv[argc-1]));

}


//----------------------------------------------------------------------

struct Label_list { 
  string label;
  Array1 <Properties_final_average> avg;
  int navg;
  Label_list(int max) { 
    avg.Resize(max);
    navg=0;
  }
};

int main(int argc, char ** argv) {
  cout.precision(10);

  if(argc <= 1) error("Need an argument");
  Gosling_options options;
  get_options(argc, argv, options);

 
  vector <Label_list> all_averages;
 
  int nfiles=options.filenames.size();
  for(int file=0; file < nfiles; file++) {
  

    ifstream input(options.filenames[file].c_str());
    
    vector < string >  words;
    parsefile(input, words);
    input.close();
    
    vector <string> labels;
    Array1 < Array1 <Properties_block > > allblocks;
    getBlocks(words, labels, allblocks);
    int nlabels=labels.size();
    
    if(options.trace) output_trace(labels, allblocks);
    if(options.trace_force) output_trace_force(labels, allblocks);
    
    
    //cout << "Last updated: $Date: 2006/11/13 18:23:07 $ GMT\n";
    //Get the blocks from each file; if only one
    //file, go ahead and print it out
    Properties_final_average avg;
    avg.showAutocorr(options.show_autocorr);
    for(int lab=0; lab < nlabels; lab++) {
      if(labels[lab]==options.label || options.label=="") {
	int neffblock;
	neffblock=reblock_average(allblocks(lab), options.reblock, 
				  options.equil, avg);
	if(neffblock==0) {
	  cout << "Skipping reblocking on " << labels[lab]
	       << " because there aren't enough blocks " << endl;
	  neffblock=reblock_average(allblocks(lab), 1, options.equil, avg);
	}

	if(nfiles < 2) { 
	  cout << "#####################" << endl;
	  cout << labels[lab] << ":  "<< allblocks(lab).GetDim(0) 
	       << " total blocks reblocked into " << neffblock << endl;
	  cout << "#####################" << endl;
	  
	  avg.showSummary(cout);
	}

	int found=0;
	for(vector<Label_list>::iterator i=all_averages.begin();
	    i!= all_averages.end(); i++) {
	  if(labels[lab] == i->label) {
	    i->avg(i->navg++)=avg;
	    found=1;
	  }
	}
	if(!found) { 
	  Label_list nlabel(nfiles);
	  nlabel.avg=avg;
	  nlabel.navg=1;
	  nlabel.label=labels[lab];
	  all_averages.push_back(nlabel);
	}

      }
    }
  }
  

  //Print out averages accumulated from 
  //all files

  if(nfiles > 1) { 
    Properties_final_average avg;
    avg.showAutocorr(options.show_autocorr);
    
    for(vector<Label_list> :: iterator lab=all_averages.begin();
	lab!= all_averages.end(); lab++) { 
      avg.averageReduce(lab->avg, 0, lab->navg);
      cout << "#########################" << endl;
      cout << lab->label << " reaccumulated \n";
      cout << "#########################" << endl;
      avg.showSummary(cout);
    }
  }

    /*
 
  if(options.two_pt_forces) {
    for(int lab=0; lab < nlabels; lab++) {
      if(labels[lab]==options.label || (options.label=="" && lab==0)) {
        avg.blockReduce(allblocks(lab), 0, allblocks(lab).GetDim(0), options.equil);
        Array2 <doublevar> forces;
        avg.twoPointForces(forces);
        for(int i=0; i < forces.GetDim(0); i++) {
          cout << "force" << i << "   " << forces(i,0) 
              << "  +/- " << forces(i,1) << endl;
        }
      }
    }
  }
    */

}
