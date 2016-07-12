/*
 
Copyright (C) 2007 Michal Bajdich

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
#include "converter.h"
#include <iostream>
#include <iomanip>
#include <cassert>
#include <string>
#include <vector>
#include <fstream>
#include <cmath>
#include <cstring>
#include <cstdlib>

using namespace std;

void reorder_orbs(vector <int> & orbs, double & det_weights){
  // sort by ascending value
  int n=orbs.size();
  int tmp;
  for(int i=0; i < n; i++) {
    if(abs(orbs[i])>99)
      cout <<"Excitation with orbital > 99. GAMESS formated file is not fully usable. Inspect your result! M.B."<<endl;
    for(int j=0; j < n; j++) {
      if (orbs[i]>orbs[j]){
	tmp=orbs[j];
	orbs[j]=orbs[i];
	orbs[i]=tmp;
	det_weights=-det_weights;
      }
    }
  }
  /*
  cout <<"reorder_orbs"<<endl;
  for(int i=0; i < n; i++) {
    cout <<orbs[i]<<"  ";
  }
  cout <<endl;
  */
     
}


void normal_reorder_orbs(int & nup, int & ndown, int & ncore, vector <int> & orbs, double & det_weights){
  // put to normal order 
  int n=orbs.size();
  int tmp;

  /*
  cout <<"before normal_reorder_orbs"<<endl;
  for(int i=0; i < n; i++) {
    cout <<orbs[i]<<"  ";
  }
  cout <<endl;
  */

  for(int i=0; i < n; i++) {
    if (orbs[i]>0) {
      if(orbs[i]!=nup-i){
	for(int j=ncore; j < nup; j++){
	  if(orbs[i]==j+1 && nup-j-1!=i ){
	    tmp=orbs[nup-j-1];
	    orbs[nup-j-1]=orbs[i];
	    orbs[i]=tmp;
	    det_weights=-det_weights;  
	  }
	}
      }
    }
    else{
      if(-orbs[i]!=(i+1-(nup))){
	//cout <<"orbs[i] "<<orbs[i]<<"i+1-(nup) "<<i+1-(nup)<<endl;
	for(int j=ncore; j < ndown; j++){	
	  if(-orbs[i]==j+1 && nup-2*ncore+j!=i ){
	    //cout <<" swap pos "<<nup-2*ncore+j<<" with pos "<<i<<endl;
	    tmp=orbs[nup-2*ncore+j];
	    orbs[nup-2*ncore+j]=orbs[i];
	    orbs[i]=tmp;
	    det_weights=-det_weights;  
	  }
	}
      }
    }
  } 
  
  /*
  cout <<"after normal_reorder_orbs"<<endl;
  for(int i=0; i < n; i++) {
    cout <<orbs[i]<<"  ";
  }
  cout <<endl;
  */
  
}


void spin_separate(vector <int> & vals, vector <int> & spinup, vector <int> & spindown, int & tmp_up, int & tmp_down){
  int n=vals.size();
  tmp_up=tmp_down=0;
  for (int i=0; i < n; i++){
    if(vals[i]>tmp_up)
      tmp_up=vals[i];
    if(vals[i]<tmp_down)
      tmp_down=vals[i];
    if(vals[i]>0)
      spinup.push_back(vals[i]);
    else
      spindown.push_back(-vals[i]);
  }
  tmp_down=-tmp_down;
  //cout <<"The largest spin-up orbital is "<<tmp_up<<" and spin-down orbital is "<<-tmp_down<<endl;
}


void sort_abs_largest_first(vector <double> & vals, vector <int> & list){
  int n=vals.size();
  list.resize(n); 
  for (int i=0; i < n; i++) 
    list[i] = i;
  
  for (int i=1; i < n; i++) {
    const double temp = vals[i];
    const double abstemp =  fabs(vals[i]);
    int j;
    for (j=i-1; j>=0 && fabs(vals[j])<abstemp; j--) {
      vals[j+1] = vals[j];
      list[j+1] = list[j];
    }
    vals[j+1] = temp;
    list[j+1] = i;
  }
  //for (int i=0; i < n; i++) 
  //cout << list[i]<<endl;
}


void usage(const char * name) {
  cout << "usage: " << name <<   " <options> <output> " << endl;
  cout << "Where options can be: \n";
  cout << "-csf      write CSFs rather than separate determinants\n";
  cout << "-o        desired output file\n";
  cout << "-wthresh  threshold weight value (Default 0.01)\n";
  cout << "-state    write determinants for selected state \n";
  cout << "-norder  orbitals in determinants follow natural order\n";
  cout << "-print_level (Default 1)\n";
  exit(1);
}

int main(int argc, char ** argv) {
  string infilename;
  string outputname;
  int computed_states=1;
  int iroot=1;
  int nstate=1;
  double wtresh=0.01;
  int symmetry=0;
  int reorder=0;
  int printout=1;
  
  for(int i=1; i< argc-1; i++) {
    if(!strcmp(argv[i], "-wthresh") && argc > i+1) {
            wtresh=double(atof(argv[++i]));
    }
    else if(!strcmp(argv[i], "-o") && argc > i+1) {
      outputname=argv[++i];
    }
    else if(!strcmp(argv[i], "-csf")) {
      symmetry=1;
    }
    else if (!strcmp(argv[i], "-state")) {
      nstate=atoi(argv[++i]);
    }
    else if (!strcmp(argv[i], "-norder")) {
      reorder=1;
    }
    else if(!strcmp(argv[i], "-print_level") && argc > i+1) {
            printout=int(atoi(argv[++i]));
    }
    else {
    cout << "Didn't understand option " << argv[i] << endl;
    usage(argv[0]);
    }
  }

  if(argc >= 2) {
    infilename=argv[argc-1];
  }
  else { usage(argv[0]); }

  if(outputname=="") {
    outputname="gamessci2qmc.out";
  }

  cout << "Starting converter"<<endl;
  cout << "Using weight treshhold "<<wtresh<<endl;

  ifstream is(infilename.c_str());
  if(!is) {
    cout << "Couldn't open" <<infilename<< endl;
    exit(1);
  }
  vector <string> words;
  string line;
  string space=" ";
  int using_guga=0;
  int using_determinants=0;
  int csfmax=0; //how many CSFs was there total
  int number_of_core_orbitals=0;
  vector <int> nelectrons(2);
  while(getline(is, line)){
    words.clear();
    split(line, space, words);
    //find GUGA DISTINCT ROW TABLE
    if(words[0]=="GUGA" && words[1]=="DISTINCT" && words[2]=="ROW" && words[3]=="TABLE"){
      //cout <<"GUGA CITYP using CSF's"<<endl;
      using_guga=1;
      while(getline(is,line)) {
	words.clear();
	split(line, space, words);
	//first read in CORE ALPHA BETA
	if(words[0]=="NUMBER" && words[2]=="CORE" && words[3]=="MOLECULAR" && words[4]=="ORBITALS"){
	  number_of_core_orbitals=atoi(words[6].c_str());
	  cout <<number_of_core_orbitals<<" core orbitals"<<endl<<endl;
	}
	if(words[0]=="NUMBER" && words[2]=="ALPHA" && words[3]=="ELECTRONS"){
	  nelectrons[0]=atoi(words[5].c_str());
	  cout <<nelectrons[0]<<" alpha electrons"<<endl;
	}
	if(words[0]=="NUMBER" && words[2]=="BETA" && words[3]=="ELECTRONS"){
	  nelectrons[1]=atoi(words[5].c_str());
	  cout <<nelectrons[1]<<" beta electrons"<<endl;
	}


	//find COMPUTING THE HAMILTONIAN FOR THE    ????? CSF-S
	if(words[0]=="COMPUTING" && words[2]=="HAMILTONIAN" && words[6]=="CSF-S..."){
	  csfmax=atoi(words[5].c_str());
	  cout <<"GUGA CITYP using "<<csfmax<<" CSF's"<<endl;
	}

        // find  NUMBER OF REQUESTED STATES
        if(words[0]=="DAVIDSON" && words[1]=="METHOD" && words[3]=="DIAGONALIZATION" ){
          while(getline(is,line)) {
            words.clear();
            split(line, space, words);
            if(words[0]=="NUMBER" && words[2]=="STATES" && words[3]=="REQUESTED"){
              computed_states=atoi(words[5].c_str());
              cout << "Number of states computed in CI calculation : "<< computed_states << endl;
              break;
            }
          }
        }

       // find IROOT value
        if(words[0]=="PROPERTIES" && words[3]=="COMPUTED" && words[5]=="STATE" && words[6]=="-IROOT-" ){
          iroot=atoi(words[7].c_str()); 
          cout <<"IROOT value in CI calculation : "<< iroot <<endl;
          break;
        }
     }
     break;
    }

    //find AMES LABORATORY DETERMINANTAL FULL CI
    if(words[0]=="AMES" && words[1]=="LABORATORY" && words[2]=="DETERMINANTAL" && words[3]=="FULL" && words[4]=="CI" ){
      cout <<"Found ALDET CITYP using determinantal CI"<<endl;
      using_determinants=1;
      while(getline(is,line)) {
        words.clear();
        split(line, space, words);
        if(words[0]=="NUMBER" && words[3]=="STATES" && words[4]=="REQUESTED" && words[5]=="="){
          computed_states=atoi(words[6].c_str());
          cout <<computed_states<<" states computed in CI calculation"<<endl;
        }
        if(words[0]=="CI" && words[1]=="PROPERTIES" && words[6]=="ROOT" && words[7]=="NUMBER"){
          iroot=atoi(words[8].c_str());
          cout <<"Properties of state "<<iroot<<" were requested in CI"<<endl;
          break;
        }
      }
      break;
    }
  }
  is.close();

  if (nstate < 1) {
    cout << "Wrong input for -state . Setting STATE value to default (STATE=IROOT).\n" ;
    nstate=iroot;
  }
  else if (nstate == 0) {
    cout << "Setting STATE value to default (STATE=IROOT).\n" ;
    nstate=iroot;
  }
  else if (nstate>computed_states) {
    cout << "STATE value is higher than number of computed states in CI.\n" ;
    cout << "Setting STATE value to default (STATE=IROOT).\n" ;
    nstate=iroot;
  }


  if(using_guga){
    if(symmetry)
      cout << "Keeping the whole CSFs" <<endl;
    vector <int> csf_full;
    vector <int> csf;
    vector <double> csf_weights;
    vector <string> csf_occupations_str;
    
    
    int line_position;
    is.open(infilename.c_str());
    while(getline(is, line)) {
      words.clear();
      split(line, space, words);
      if(words[0]=="DAVIDSON" && words[1]=="METHOD"){
        if(printout>1) cout <<"found DAVIDSON METHOD CI-MATRIX DIAGONALIZATION"<<endl;
        //read in info on  DAVIDSON METHOD CI-MATRIX DIAGONALIZATION   
        while(getline(is, line)) {
          words.clear();
          split(line, space, words);
          if(words[0]=="ITER." && words[2]=="IMPROVED") {
            while(getline(is, line)) {
              if(printout>1) cout << line <<endl;
              if(line=="") break;
            }
          }//---done energy info
          //cout << line<<endl;
          words.clear();
          split(line, space, words);
          //read in contributing CSFs
          if (words[0]=="STATE"&& words[1]=="#" && words[3]=="ENERGY" && atoi(words[2].c_str())==nstate){
            cout << "Printing CSFs for state : " << atoi(words[2].c_str()) << endl;
            int i=0;
            double smallest_csf_weight=1.0;
            while(getline(is,line)) {
              words.clear();
              split(line, space, words);
              //cout << "in state : " << line<<endl;
              if(words.size() >0 and (words[0]=="......" || words[0]=="STATE")) break;

              if(words.size()>1 && atoi(words[0].c_str())>0  && !(words[0]=="---") && !(words[0]=="CSF")){
                csf.push_back(atoi(words[0].c_str())-1);
                csf_weights.push_back(atof(words[1].c_str()));
                csf_occupations_str.push_back(words[2]);
                if(smallest_csf_weight > abs(csf_weights[i]))
                  smallest_csf_weight=abs(csf_weights[i]);
                i++;
              }
              else if (words.size()==1){
                csf_occupations_str.back()+=words[0];
              }
            }
            if(csf.size()>0){
              cout <<"Found "<<csf.size()<<" contributing CSFs with smallest weight of "<<smallest_csf_weight<<endl;
              if(printout>1)
                cout <<"MB: to increase printout set PRTTOL smaller than "<<smallest_csf_weight<<" in $GUGDIA"<<endl;
            }
            else{
              cout <<"did not find any CSF COEF OCCUPANCY array in the DAVIDSON METHOD CI-MATRIX DIAGONALIZATION printout"<<endl;
              exit(1);
            }
          }
        }
      }//if DAVIDSON METHOD CI-MATRIX DIAGONALIZATION
    }
    is.close();

    vector <vector <double> > det_weights(csf.size());
    vector <vector <string> > det_str(csf.size());



    int counter_csf=0;
    is.open(infilename.c_str());
    while(getline(is, line)) {
      words.clear();
      split(line, space, words);
      //here where all CSFs printout starts
      if (words[0]=="DETERMINANT" && words[1]=="CONTRIBUTION") {
	if(printout>1) cout <<"looping to read the relevant CSFs"<<endl;
	int kk=0;
	for(int k=0;k<csf.size();k++){
	  if(printout>1) cout <<"looking for CSF with index "<<csf[k]+1<<" .........";
	  while(getline(is,line)) {
	    words.clear();
	    split(line, space, words);
	    if(words[0]=="......" && words[1]=="END" && words[2]=="OF" && words[3]=="-DRT-" && words[4]=="GENERATION"){
	      cout <<"\n found the end of DRT GENERATION"<<endl;
	      if(kk!=csfmax){
		cout <<"did not found enought CSF'S!, did you run gamess with NPRT=2 in $CIDRT?"<<endl;
		exit(1);
	      }
	      break;
	    }
	    if(words[0]=="CASE" && words[1]=="VECTOR" && words[2]== "=" && words.size()==4){
	      if(atoi(words[3].c_str()) > csf[k]+1){
	     	cout <<"error, the CSF with index "<<csf[k]+1<<" could not be found"<<endl;
	      	exit(1);
	      }
	      if(atoi(words[3].c_str())-1==csf[k]){
		if(printout>1)  cout <<"found!"<<endl;
		counter_csf++;
		csf_full.push_back(atoi(words[4].c_str())-1);


		//read in the determinants in:
		double sum_of_weights=0.0;
		int tries=0;
		int max_tries=1000;
		
		while(getline(is,line)) {
		  words.clear();
		  split(line, space, words);
		  /*
		  the old way
		  if ( line.size() > 33  ){
		    if (line.substr(13,2)=="C(" ){
		      if(printout>2) cout <<line<<endl;
		      det_weights[k].push_back(atof(line.substr(20,10).c_str()));
		      det_str[k].push_back(line.substr(33));
		      line_position=is.tellg();
		      //if there is a line which is too long this while has to be there
		      while(getline(is,line)) {
			words.clear();
			split(line, space, words);
			if ((words[0]=="C(" ) || (words[0]=="CASE" && words[1]=="VECTOR") || line.size()==0){
			  is.seekg(line_position);
			  break;
			}
			else{
			  det_str[k].back() += " ";
			  det_str[k].back() +=line.substr(33);
			}
		      }
		    }//if
		  }//if
		  */
		  //new way
		  if ( line.size() > 33  ){
		    if (line.substr(13,2)=="C(" ){
		      if(printout>2) cout <<line<<endl;
		      det_weights[k].push_back(atof(line.substr(20,10).c_str()));
		      det_str[k].push_back(line.substr(33));
		      //MB: new way of doing things: if there is a line which is too long this while has to be there
		      int sum_of_string_leghts=det_str[k].back().size();
		      while(sum_of_string_leghts <  (nelectrons[0]+nelectrons[1]-2*number_of_core_orbitals)*3){
			getline(is,line);
			split(line, space, words);
			string remainder=line.substr(33);
			det_str[k].back() += " ";
			det_str[k].back() +=remainder;
			sum_of_string_leghts+=remainder.size();
		      }
		      sum_of_weights+=det_weights[k].back()*det_weights[k].back();
		      if(printout>3) cout <<"w: "<<det_weights[k].back()<<" at line number "<<line_position<<" sum of weights: "<<sum_of_weights<<endl;
		      if(abs(sqrt(sum_of_weights)-1.0)<1e-5)
			break;
		      //M.B.: this was working before, but tellg() was giving negative number for files > 2Gb
		      //if there is a line which is too long this while has to be there
		      /*
		      line_position=is.tellg();
		      while(getline(is,line)) {
			words.clear();
			split(line, space, words);
			tries++;
			if ((words[0]=="C(" ) || (words[0]=="CASE" && words[1]=="VECTOR") || line.size()==0){
			  is.seekg(line_position);
			  cout <<"going back to line position" <<line_position<<endl;
			  break;
			}
			else{
			  if(printout>2) cout <<line<<endl;
			  det_str[k].back() += " ";
			  det_str[k].back() +=line.substr(33);
			}
		      }//while
		      */
		      
		    }
		    if(tries>max_tries){
		      cout<<"could not fully read CSF "<<csf[k]+1<<" ...exiting"<<endl; 
		      exit(1); 
		    }
		    tries++;
		  }//if(line.size() > 33)
		}//end of while
		if(!(abs(sqrt(sum_of_weights)-1.0)<1e-5)){
		  cout<<"could not fully read CSF "<<csf[k]+1<<" ...exiting"<<endl; 
		  exit(1); 
		}
		//need to break in order to find another CSF
		break;
	      }//end of if CASE VECTOR = csf[k]
	      if(kk>csfmax){
		cout<<"Number of CSF's is greater than designed limit of "<<csfmax<<"CSF's"<<endl;
		exit(1);
	      }
	      kk++;
	    }
	  }//while(getline(is,line))
	}//int k loop
      }//end of reading all CSFs
    }
    is.close();

    if(counter_csf!=csf.size()){
      cout <<"counter_csf "<<"and "<<csf.size()<<"dont match "<<endl;
      exit(1);
    }
    cout << "done readout"<<endl;
    //end of read out
    //storing determinants and reodering orbs
    det_weights.resize(csf.size());
    det_str.resize(csf.size());
    vector < vector < vector <int> > > det(csf.size());
    
    for(int i=0;i<det_weights.size();i++){
      if(!det_weights[i].size()){
	cout<<" CSF with index "<<csf[i]+1<<" was zero determinants in!, exitting"<<endl;
	exit(1);
      }
      
      det[i].resize(det_weights[i].size());
      for(int j=0;j<det_weights[i].size();j++){
	for(int k=0; k < det_str[i][j].size(); k=k+3){
	  det[i][j].push_back(atoi(det_str[i][j].substr(k,3).c_str()));
	}
	
	reorder_orbs(det[i][j],det_weights[i][j]);
	if(reorder)
	  normal_reorder_orbs(nelectrons[0],nelectrons[1], number_of_core_orbitals, det[i][j],det_weights[i][j]);

	//if(i>0)
	//  reorder3(det[i][j],det_weights[i][j],det[0][0]);
	// det_spinup[i].resize(det_weights[i].size());
	// det_spindown[i].resize(det_weights[i].size());
	// spin_separate(det[i][j],det_spinup[i][j],det_spindown[i][j] );
	//cout << "Weight: "<<det_weights[i][j]<<"   state: ";
	//for(int k=0; k < det[i][j].size();k++){
	//  cout << det[i][j][k]<< " ";
	// }
	//cout <<endl;
      }
    }
    
    cout << "done storing determinants and reodering orbs"<<endl;
    
    //storing occupation array  
    vector < vector <int> > csf_occupation(csf.size());
    for(int i=0;i< csf.size();i++){
      for(int j=0;j< csf_occupations_str[i].size();j++){
	csf_occupation[i].push_back(atoi(csf_occupations_str[i].substr(j,1).c_str()));
      }
    }
    cout << "done storing occupation array"<<endl;
    
    if(symmetry==0){
      // assign full weight to each determimant
      for(int i=0;i<csf.size();i++){
	for(int j=0;j<det_weights[i].size();j++){
	  //if(csf[i]==847){
	  //cout << csf_weights[i]<< " "<<det_weights[csf[i]-1][j]<<endl;
	  //}
	  det_weights[i][j]*=csf_weights[i];
	}
      }
      cout << "done assign full weight to each determimant"<<endl;
      
      //find the same determimants and add their contributions
      for(int i=0;i<csf.size();i++){
	//cout << csf[i]<<" "<< csf_occupations_str[i]<<endl;
	for(int j=i+1;j<csf.size();j++){
	  if(csf_occupations_str[i]==csf_occupations_str[j]){
	    //cout << "CSF "<<csf[i]<<" and "<<csf[j]<<  " are the same"<<endl;
	    for(int k=0;k<det[i].size();k++)
	      for(int l=0;l<det[j].size();l++)
		if(det[i][k]==det[j][l]){
		  det_weights[i][k]+=det_weights[j][l];
		  det_weights[j][l]=0.0;
		  //cout << "determimant"<<k<<" and "<<l<<" are the same"<<endl;
		  //for(int h=0;h<det[csf[i]-1][k].size();h++)
		  //  cout <<det[csf[i]-1][k][h]<<" ";
		  //cout <<endl;
		  //for(int h=0;h<det[csf[j]-1][k].size();h++)
		  //  cout <<det[csf[j]-1][l][h]<<" ";
		  //cout <<endl;
		}
	  }
	}
      }
      cout << "done finding the same determimants and adding their contributions"<<endl;
    }
  

    //use weight treshold and store for printout
    vector <double> det_weights_printing;
    //int ircounter=0;
    //vector <int> det_weights_bonds;
    vector < vector <int> > det_printing;
    vector < vector <int> > det_printing_old;
    vector <double> csf_weights_printing;
    vector <int> csf_printing;
    vector < vector <double> >  det_weights_printing_symmetry;
    for(int i=0;i<csf.size();i++){
      if(symmetry){
	if(abs(csf_weights[i])> wtresh){
	  csf_weights_printing.push_back(csf_weights[i]);
	  csf_printing.push_back(i);
	  //det_weights_printing_symmetry.push_back(det_weights[csf[]]);
	  //for(int k=0;k<det[csf[i]].size();k++){
	  // det_weights_printing.push_back(det_weights[csf[i]][k]);
	  //det_printing.push_back(det[csf[i]][k]);
	    //det_weights_bonds.push_back(ircounter);
	  //}
	  //ircounter++;
	}
      }
      else{
	for(int k=0;k<det[i].size();k++){
	  if(abs(det_weights[i][k])> wtresh){
	    // cout <<"det_weights "<< det_weights[csf[i]][k]<<endl;; 
	    det_weights_printing.push_back(det_weights[i][k]);
	    det_printing_old.push_back(det[i][k]);
	  }
	}
      }
    }
    
    //sorting
    vector <int> list;
    if(symmetry){
      sort_abs_largest_first(csf_weights_printing, list);
      for(int i=0;i<csf_weights_printing.size();i++){
	det_weights_printing_symmetry.push_back(det_weights[csf_printing[list[i]]]);
	for(int k=0;k<det[csf_printing[list[i]]].size();k++){
	  det_weights_printing.push_back(det_weights[csf_printing[list[i]]][k]);
	  det_printing.push_back(det[csf_printing[list[i]]][k]);
	}
      }
    }
    else{
      sort_abs_largest_first(det_weights_printing, list);
      for (int i=0;i<det_weights_printing.size();i++){
	det_printing.push_back(det_printing_old[list[i]]);
      }
      
    }


    if(symmetry){
      cout << "After cutoff applied: found "<<det_weights_printing.size()<<" determinats with weights"<<endl;
      cout << "and number of independent weights "<<csf_weights_printing.size()<<endl;
    }
    else
      cout << "After cutoff applied: found "<<det_weights_printing.size()<<" unique determinats with weights"<<endl; 
    
    
    vector < vector <int> > det_up(det_weights_printing.size());
    vector < vector <int> > det_down(det_weights_printing.size());
    int max_up, max_down, tmp_up,tmp_down;
    max_up=max_down=0;
    for (int i=0;i<det_weights_printing.size();i++){
      spin_separate(det_printing[i], det_up[i], det_down[i], tmp_up, tmp_down);
      if(tmp_up>max_up)
	max_up=tmp_up;
      if(tmp_down>max_down)
	max_down=tmp_down;
    }
    cout <<"the largest spin-up orbital is "<<max_up<<" and spin-down orbital is "<<max_down<<endl;
    

    

    //FINAL OUTPUT 
    ofstream output(outputname.c_str());
    string gap="   ";
    
    if(symmetry){
      output <<gap<< "# NCSF = "<<csf_weights_printing.size()<<endl;
      for (int i=0;i<csf_weights_printing.size();i++){
	output.precision(7);
	output.width(10);
	output <<gap<<"CSF {  "<< csf_weights_printing[i] <<"  ";
	for(int j=0;j<det_weights_printing_symmetry[i].size();j++){
	  output.precision(7);
	  output.width(10);
	  output << det_weights_printing_symmetry[i][j]<<"  ";	
	}
	output << "}"<<endl;
      }
    }
    else {
      output <<gap<< "DETWT {"<<endl;
      output <<gap<< "# NDET = "<<det_weights_printing.size()<<endl;
      output <<gap;
      for (int i=0;i<det_weights_printing.size();i++){
	output.precision(7);
	output.width(10);
	output << det_weights_printing[i] <<"  ";
	if ((i+1)%6==0)
	  output <<endl<<gap;
      }
      output << "}"<<endl;
    }
    output <<gap<<"STATES {"<<endl<<gap;
    if(symmetry){
      int counter=0;
      for (int i=0;i<csf_weights_printing.size();i++){
	output <<"#  CSF "<<i+1<<": weight: "<<csf_weights_printing[i]<<endl<<gap;
	for(int j=0;j<det_weights_printing_symmetry[i].size();j++){
	  output <<"#  Determinant "<<j+1<<": weight: "<<det_weights_printing_symmetry[i][j]<<endl<<gap;
	  for(int m=0;m<number_of_core_orbitals;m++)
	    output << m+1 <<"  ";
	  for(int k=det_up[counter].size();k>0;k--)
	    output << det_up[counter][k-1] <<"  ";
	  output <<endl<<gap;
	  for(int m=0;m<number_of_core_orbitals;m++)
	    output << m+1 <<"  ";
	  for(int k=0;k<det_down[counter].size();k++)
	    output << det_down[counter][k] <<"  ";
	  output <<endl<<gap;
	  counter++;
	}
      }
    }
    else {
      for (int i=0;i<det_weights_printing.size();i++){
	output <<"#  Determinant "<<i+1<<": weight: "<<det_weights_printing[i]<<endl<<gap;
	for(int j=0;j<number_of_core_orbitals;j++)
	  output << j+1 <<"  ";
	for(int k=det_up[i].size();k>0;k--)
	  output << det_up[i][k-1] <<"  ";
	output <<endl<<gap;
	for(int j=0;j<number_of_core_orbitals;j++)
	  output << j+1 <<"  ";
	for(int k=0;k<det_down[i].size();k++)
	  output << det_down[i][k] <<"  ";
	output <<endl<<gap;
      }
    }
    output << "}"<<endl;
    output.close();
    cout <<"done"<<endl;
  }//end using guga
  else if(using_determinants){
    vector <double> det_weights;
    vector < vector <string> > det_occupations_str(2);
    int ncorr_orbs;
    int nacc_orbs;
    vector <int> nelectrons(2);
    int orbs_total;
    int ndet;
    int cycle=1;
    if(symmetry)
      cout << "ignoring option -csf" <<endl;
     is.open(infilename.c_str());
     while(getline(is, line) && cycle<2 ) {
       words.clear();
       split(line, space, words);
       
       if(words[0]=="NUMBER" && words[2]=="CORE" && words[3]=="ORBITALS"){
	 ncorr_orbs=atoi(words[5].c_str());
	 cout <<ncorr_orbs<<" core orbitals"<<endl;
       }
       if(words[0]=="NUMBER" && words[2]=="ACTIVE" && words[3]=="ORBITALS"){
	 nacc_orbs=atoi(words[5].c_str());
	 cout <<nacc_orbs<<" active orbitals"<<endl;
       }
       
       if(words[0]=="NUMBER" && words[2]=="ALPHA" && words[3]=="ELECTRONS"){
	 nelectrons[0]=atoi(words[5].c_str());
	 cout <<nelectrons[0]<<" alpha electrons"<<endl;
       }
       if(words[0]=="NUMBER" && words[2]=="BETA" && words[3]=="ELECTRONS"){
	 nelectrons[1]=atoi(words[5].c_str());
	 cout <<nelectrons[1]<<" beta electrons"<<endl;
       }
       if(words[0]=="NUMBER" && words[2]=="OCCUPIED" && words[3]=="ORBITALS" && words[4]=="="){
	 orbs_total=atoi(words[5].c_str());
	 cout <<orbs_total<<" total number orbitals"<<endl;
       }
       
       if(words[0]=="ITERATION" && words[1]=="ENERGY" && words[2]=="GRADIENT") {
	 while(getline(is, line)) {
	   words.clear();
	   split(line, space, words);
	   cout << line <<endl;
	   if(words[0]=="ALL" && words[1]=="STATES" && words[2]=="CONVERGED.") break;
	 }
       }//---done energy info
       //read in actual determinamnts

      if (words[0]=="STATE" && words[2]=="ENERGY=" && atoi(words[1].c_str())==nstate){
        while(getline(is, line)) {
          words.clear();
          split(line, space, words);
          ndet=0;
          if(words[0]=="ALPHA" && words[2]=="BETA" && words[4]=="COEFFICIENT") {
            while(getline(is, line)) {
              words.clear();
              split(line, space, words);
              if(words[1]=="|" && words[3]=="|"){
//                 cout << words[0]<<" "<<words[2]<<" "<<words[4]<<endl;
                det_weights.push_back(atof(words[4].c_str()));
                det_occupations_str[0].push_back(words[0]);
                det_occupations_str[1].push_back(words[2]);
                ndet++;
              }
              //..... DONE WITH DETERMINANT CI COMPUTATION .....
              //find the end
              if(words[0]=="....." || line.size()==0) {
                cout <<"found "<<ndet<<" determinants"<<endl;
                cycle++;
                break;
              }
            }
           break;
          }
        }
      }
    }
     is.close();
     cout <<"done readout"<<endl;
     
     vector < vector < vector <int> > > det_occupations(2);
    
     for(int spin=0;spin<2;spin++){
       det_occupations[spin].resize(ndet);
       for(int det=0;det<ndet;det++){
	 for(int j=0;j<ncorr_orbs;j++){
	   det_occupations[spin][det].push_back(j+1);
	 }
	 for(int j=0;j<nacc_orbs;j++){
	   if(atoi(det_occupations_str[spin][det].substr(j,1).c_str())==1)
	     det_occupations[spin][det].push_back(ncorr_orbs+j+1);
	 }
       }//ndets
     }//spin
     

     
     
     //select determinants
     vector <double> det_weights_print;
     vector < vector < vector <int> > > det_occupations_print(2);
     for(int det=0;det<ndet;det++){
       if(fabs(det_weights[det])>wtresh){
	 det_weights_print.push_back(det_weights[det]);
	 for(int spin=0;spin<2;spin++){
	   det_occupations_print[spin].push_back(det_occupations[spin][det]);
	 }
       }
       else{
	 //cout <<" det "<<det<<endl;
       }
     }
     cout <<"done selecting determinants"<<endl; 
     cout <<det_weights_print.size()<<" determinants has weight larger than "<<wtresh<<endl;
          
     vector <int> list;
     sort_abs_largest_first(det_weights_print, list);

     //final printout
     ofstream output(outputname.c_str());
     string gap="   ";
     output <<gap<< "DETWT {"<<endl;
     output <<gap<< "# NDET = "<<det_weights_print.size()<<endl;
     output <<gap;
     for (int i=0;i<det_weights_print.size();i++){
	output.precision(7);
	output.width(10);
	output << det_weights_print[i] <<"  ";
	if ((i+1)%6==0)
	  output <<endl<<gap;
     }
     output << "}"<<endl;
     output <<gap<<"STATES {"<<endl<<gap;
     for (int i=0;i<det_weights_print.size();i++){
       output <<"#  Determinant "<<i+1<<": weight: "<<det_weights_print[i]<<endl<<gap;
       for(int spin=0;spin<2;spin++){
	 for(int j=0;j<det_occupations_print[spin][list[i]].size();j++)
	   output <<det_occupations_print[spin][list[i]][j]<<"  ";
	 output <<endl<<gap;
       }
     }
     output << "}"<<endl;


     output.close();
     cout <<"done"<<endl;

     


  }
  else{
    cout <<"Could not find GUGA nor ALDET type CI";
    exit(1);
  }
}
