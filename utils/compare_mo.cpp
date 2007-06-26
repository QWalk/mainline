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
#include "Array.h"


void read_mo(string & orbin, int nmo, int nfunc, Array2 <double> & moCoeff)
{
  ifstream in(orbin.c_str());
  if(!in) error("Couldn't open ", orbin);
  string dummy;
  in >> dummy;
  while(dummy != "COEFFICIENTS")
    in >> dummy;

  moCoeff.Resize(nmo, nfunc);
  for(int m=0; m< nmo; m++) {
    for(int f=0; f< nfunc; f++) {
      in >> moCoeff(m,f);
     }
  }
}

#include <iomanip>
/*!
Compare the sets of molecular orbitals, and
return true if they're the same, and false if not
Currently it doesn't use the lcao_overlap matrix, but
should be modified to do that; it's approximate otherwise.
*/
bool compare_mo(Array2 <double> & oldMOCoeff,
                Array2 <double> & newMOCoeff,
                Array2 <double> & lcao_overlap,
                Array1 <int> compare_list) {

 assert(oldMOCoeff.GetDim(1)==newMOCoeff.GetDim(1));

 int nfunctions=oldMOCoeff.GetDim(1);

 int ncompare=compare_list.GetDim(0);
 if(newMOCoeff.GetDim(0) < ncompare)
   error("the punch file doesn't have enough molecular orbitals to compare.");
 if(oldMOCoeff.GetDim(0) < ncompare)
   error("the old punch file doesn't have enough MO's to compare.");

  vector <int> unresolved_mos;

  //First check to see if the mo's are in the same place
  //(most should be)
  cout.precision(15);
  for(int i=0; i< ncompare; i++) {
    int mo=compare_list(i);
    double dot=0, mag_old=0, mag_new=0;
    //cout << "-----------------------mo " << mo << endl;
    for(int f=0; f< nfunctions; f++) {
      //cout << setw(10) << "function " << setw(4) << f
      //     << " new " << setw(15) << newMOCoeff(mo,f)
      //     << " old " << setw(15) << oldMOCoeff(mo,f)
      //     << "  diff "  
      //     << setw(15) << newMOCoeff(mo,f)-oldMOCoeff(mo,f)
      //     << endl;
      dot+=newMOCoeff(mo,f)*oldMOCoeff(mo,f);
      mag_old+=oldMOCoeff(mo,f)*oldMOCoeff(mo,f);
      mag_new+=newMOCoeff(mo,f)*newMOCoeff(mo,f);
    }
    dot /= sqrt(mag_old*mag_new);
    dot =fabs(dot);
    cout << "mo " << mo << "  dot " << dot << endl;
    if(fabs(dot-1) > .01) {
      unresolved_mos.push_back(mo);
    }
  }

  int nunresolved=unresolved_mos.size();
  for(int i=0; i< nunresolved; i++) {
    cout << "not matched: " << unresolved_mos[i] << endl;
  }
  bool are_same=true;
  //See if any just swapped..
  for(int i=0; i< nunresolved; i++) {
    int mo1=unresolved_mos[i];
    bool resolved_swapping=false;

    for(int j=0; j< nunresolved; j++) {
      int mo2=unresolved_mos[j];

      double dot=0, mag_old=0, mag_new=0;
      for(int f=0; f< nfunctions; f++) {
        dot+=newMOCoeff(mo1,f)*oldMOCoeff(mo2,f);
        mag_old+=oldMOCoeff(mo2,f)*oldMOCoeff(mo2,f);
        mag_new+=newMOCoeff(mo1,f)*newMOCoeff(mo1,f);
      }
      dot /= sqrt(mag_old*mag_new);
      dot=fabs(dot);
      if(fabs(dot-1) < .01) {
        cout << "switched orbital: mo " << mo2 << " went to " << mo1
             << " dot product " << dot
             << endl;
        resolved_swapping=true;
      }
    }

    if(!resolved_swapping) {
      cout << "Unresolvable change in mo " << mo1 << endl;
      are_same=false;
    }

  }

  return are_same;
}

//----------------------------------------------------------------------

int main(int argc, char * argv[]) {

  if(argc <5) {
    cout << "usage: compare_mo <nfunc> <nmo> <file1> <file2>" << endl;
    exit(1);
  }
  cout.precision(15);

  string file1=argv[3];
  string file2=argv[4];
  int nmo=atoi(argv[2]);
  int nfunc=atoi(argv[1]);
  Array2 <doublevar> mo1(nmo, nfunc), mo2(nmo, nfunc);
  Array2 <doublevar> lcao_overlap(nfunc, nfunc);
  lcao_overlap=0.0;
  for(int i=0; i< nfunc; i++)
    lcao_overlap(i,i)=1.0;  
  
  Array1<int> compare_list(nmo);
  for(int m=0; m< nmo; m++) compare_list(m)=m;

  read_mo(file1, nmo, nfunc, mo1);
  read_mo(file2, nmo, nfunc, mo2);
  compare_mo(mo1, mo2, lcao_overlap, compare_list);
}
