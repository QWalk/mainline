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
using namespace std;

void reorder(vector <int> & orbs, double & det_weights){
  // sort by ascending value
  int n=orbs.size();
  int tmp;
  for(int i=0; i < n; i++) {
      for(int j=0; j < n; j++) {
        if (orbs[i]>orbs[j]){
          tmp=orbs[j];
          orbs[j]=orbs[i];
          orbs[i]=tmp;
          det_weights=-det_weights;
        }
      }
  }
}

void spin_separate(vector <int> & vals, vector <int> & spinup, vector <int> & spindown){
  int n=vals.size();
  for (int i=0; i < n; i++){
    if(vals[i]>0)
      spinup.push_back(vals[i]);
    else
      spindown.push_back(-vals[i]);
  }
}

int main(int argc, char ** argv) {
  string infilename;
  string outputname;
  double wtresh=0.01;
  int symmetry=0;
  
  for(int i=1; i< argc; i++) {
    if(!strcmp(argv[i], "-i") && argc > i+1) {
      infilename=argv[++i];
    }
    else if(!strcmp(argv[i], "-wthresh") && argc > i+1) {
            wtresh=double(atof(argv[++i]));
    }
    else if(!strcmp(argv[i], "-o") && argc > i+1) {
      outputname=argv[++i];
    }
    else if(!strcmp(argv[i], "-keep_symmetry")) {
      symmetry=1;
    }
    else {
    cout << "Didn't understand option " << argv[i]
         << endl;
    cout << "Use either -i or -wthresh"<<endl;
    }
  }
  
  if(infilename=="") {
    infilename="gam.out"; 
  }
  
  if(outputname=="") {
    outputname=infilename+".o";
  }

  cout << "Starting converter"<<endl;

  int csfmax=1000000;
  vector <int> csf_full;
  vector <int> csf;
  vector <double> csf_weights;
  vector <string> csf_occupations_str;
  vector <vector <double> > det_weights(csfmax);
  vector <vector <string> > det_str(csfmax);
  
  
  cout << "Using weight treshhold "<<wtresh<<endl;
  if(symmetry)
    cout << "Keeping symmetry in determinant weights" <<endl<<endl;
  
  ifstream is(infilename.c_str());
  if(!is) {
    cout << "Couldn't open" <<infilename<< endl;
    exit(1);
  }

  string line;
  string space=" ";
  vector <string> words;
  int i,k;
  int counter, counter_csf;
  while(getline(is, line)) {
    split(line, space, words);
    if (words[0]=="DETERMINANT" && words[1]=="CONTRIBUTION") {
      k=0;
      while(getline(is,line)) {
        words.clear();
        split(line, space, words);
        //cout << line << endl;
        if(words[0]=="......" && words[1]=="END" && words[2]=="OF" && words[3]=="-DRT-" && words[4]=="GENERATION"){
          //cout <<"found the end of DRT GENERATION"<<endl;
          break;
        }
   
        if(words[0]=="CASE" && words[1]=="VECTOR"){
          csf_full.push_back(atoi(words[4].c_str())-1);
          k++;
        }
        
        if ( line.size() ){
          if ( (words[0]=="CSF" && words[2]=="C(" )|| (words[0]=="C(" )){
            det_weights[k-1].push_back(atof(line.substr(20,10).c_str()));
            det_str[k-1].push_back(line.substr(33));
          }
        }
      }
    }
    counter_csf=k;
    words.clear();
    
    split(line, space, words);
    if(words[0]=="ITER." && words[2]=="IMPROVED") {
      while(getline(is, line)) {
        cout << line <<endl;
        if(line=="") break;
      }
    }//---done energy info
    words.clear();
    
    split(line, space, words);
    if (words[0]=="CSF" && words[1]=="COEF" && words[2]=="OCCUPANCY") {
      i=0;
      while(getline(is,line)) {
        words.clear();
        split(line, space, words);
        if(words[0]=="......") break;
       
        if(words.size()>1 && atoi(words[0].c_str())>0  && !(words[0]=="---")){
          csf.push_back(atoi(words[0].c_str())-1);
          csf_weights.push_back(atof(words[1].c_str()));
          csf_occupations_str.push_back(words[2]);
          i++;
        }
        else if (words.size()==1){
          csf_occupations_str.back()+=words[0];
        }
        
      }
      counter=i;
    }
        
  }
  is.close();
  cout << "done readout"<<endl;
  //end of read out

  //storing determinants and reodering orbs
  det_weights.resize(counter_csf);
  det_str.resize(counter_csf);
  vector < vector < vector <int> > > det(counter_csf);
  //vector < vector < vector <int> > > det_spinup(counter_csf);
  //vector < vector < vector <int> > > det_spindown(counter_csf);
  
  for(int i=0;i<det_weights.size();i++){
    det[i].resize(det_weights[i].size());
    for(int j=0;j<det_weights[i].size();j++){
      for(int k=0; k < det_str[i][j].size(); k=k+3){
        det[i][j].push_back(atoi(det_str[i][j].substr(k,3).c_str()));
      }
      reorder(det[i][j],det_weights[i][j]);
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
  vector < vector <int> > csf_occupation(counter);
  for(int i=0;i< csf_occupations_str.size();i++){
    for(int j=0;j< csf_occupations_str[i].size();j++){
      csf_occupation[i].push_back(atoi(csf_occupations_str[i].substr(j,1).c_str()));
    }
  }
  cout << "done storing occupation array"<<endl;

  // assign full weight to each determimant
  for(int i=0;i<csf.size();i++){
    for(int j=0;j<det_weights[csf[i]].size();j++){
      //if(csf[i]==847){
      //cout << csf_weights[i]<< " "<<det_weights[csf[i]-1][j]<<endl;
      //}
      det_weights[csf[i]][j]*=csf_weights[i];
    }
  }
  
  cout << "done assign full weight to each determimant"<<endl;

  //find the same determimants and add their contributions
  if(symmetry==0){
    for(int i=0;i<csf.size();i++){
      //cout << csf[i]<<" "<< csf_occupations_str[i]<<endl;
      for(int j=i+1;j<csf.size();j++){
        if(csf_occupations_str[i]==csf_occupations_str[j]){
          //cout << "CSF "<<csf[i]<<" and "<<csf[j]<<  " are the same"<<endl;
          for(int k=0;k<det[csf[i]].size();k++)
            for(int l=0;l<det[csf[j]].size();l++)
              if(det[csf[i]][k]==det[csf[j]][l]){
                det_weights[csf[i]][k]+=det_weights[csf[j]][l];
                det_weights[csf[j]][l]=0.0;
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
  int ircounter=0;
  vector <int> det_weights_bonds;
  vector < vector <int> > det_printing;
  for(int i=0;i<csf.size();i++){
    if(symmetry){
      if(abs(csf_weights[i])> wtresh){
        for(int k=0;k<det[csf[i]].size();k++){
          det_weights_printing.push_back(det_weights[csf[i]][k]);
          det_printing.push_back(det[csf[i]][k]);
          det_weights_bonds.push_back(ircounter);
        }
        ircounter++;
      }
    }
    else{
      for(int k=0;k<det[csf[i]].size();k++){
        if(abs(det_weights[csf[i]][k])> wtresh){
          //      cout <<"det_weights "<< det_weights[csf[i]][k]<<endl;; 
          det_weights_printing.push_back(det_weights[csf[i]][k]);
          det_printing.push_back(det[csf[i]][k]);
        }
      }
    }
  }
  if(symmetry){
    cout << "found "<<det_weights_printing.size()<<" determinats with weights"<<endl;
    cout << "number of independent weights "<<ircounter<<endl;
  }
  else
    cout << "found "<<det_weights_printing.size()<<" unique determinats with weights"<<endl; 
  

  vector < vector <int> > det_up(det_weights_printing.size());
  vector < vector <int> > det_down(det_weights_printing.size());
  for (int i=0;i<det_weights_printing.size();i++)
    spin_separate(det_printing[i], det_up[i], det_down[i]);

   
  //FINAL OUTPUT 
  ofstream output(outputname.c_str());
  string gap="   ";
  
  output <<gap<< "DETWT {"<<endl;
  output <<gap<< "# NDET = "<<det_weights_printing.size()<<endl;
  output <<gap;
  for (int i=0;i<det_weights_printing.size();i++){
    output.precision(6);
    output.width(10);
    output << det_weights_printing[i] <<"  ";
    if ((i+1)%6==0)
      output <<endl<<gap;
  }
  output << "}"<<endl;

  if(symmetry){
   output <<gap<< "DETWT_BONDS {"<<endl; 
   output <<gap;
   for (int i=0;i<det_weights_bonds.size();i++){
     output.precision(3);
     output.width(3);
     output << det_weights_bonds[i] <<"  ";
     if ((i+1)%16==0)
       output <<endl<<gap;
   }
   output << "}"<<endl;
  }

  output <<gap<<"STATES {"<<endl<<gap;
  for (int i=0;i<det_weights_printing.size();i++){
    for(int k=det_up[i].size();k>0;k--)
       output << det_up[i][k-1] <<"  ";
    output <<endl<<gap;
    for(int k=0;k<det_down[i].size();k++)
       output << det_down[i][k] <<"  ";
    output <<endl<<gap;
  }
  output << "}"<<endl;
  output.close();
  cout <<"done"<<endl;
}
