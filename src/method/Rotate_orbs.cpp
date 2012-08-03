/*
 
Copyright (C) 2007 Zachary Helms
 with further modifications by Lucas K. Wagner

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

#include "Program_options.h"
#include "Rotate_orbs.h"
#include "qmc_io.h"
#include "System.h"
#include "MatrixAlgebra.h"
#include "ulec.h"
/*!
Read the "words" from the method section in the input file
via doinput() parsing, and store section information in private
variables orbs, resolution, and minmax.
Set up MO_matrix and Sample_point objects for wavefunction from input
*/
void Rotate_orbs_method::read(vector <string> words,
                       unsigned int & pos,
                       Program_options & options)
{
  cout << "rrreeead " << endl;
  
  sys=NULL;
  allocate(options.systemtext[0],  sys);

  vector <vector < string> > orbgroup_txt;
  pos=0;
  vector < string> dummy;
  while( readsection(words,pos,dummy,"ORB_GROUP")) { 
    orbgroup_txt.push_back(dummy);
  }
  cout << "aaaahree " << endl;
  orbital_groups.Resize(orbgroup_txt.size());
  for(int i=0; i< orbital_groups.GetDim(0); i++) { 
    int norb_this=orbgroup_txt[i].size();
    orbital_groups[i].Resize(norb_this);
    for(int j=0; j< norb_this; j++) { 
      orbital_groups[i][j]=atoi(orbgroup_txt[i][j].c_str())-1;
    }
  }

  cout << "here " << endl;


  vector <string> orbtext;
  if(!readsection(words, pos=0, orbtext, "ORBITALS")){
      error("Need ORBITALS");
  }
  allocate(orbtext, sys, mymomat);

  mywalker=NULL;
  sys->generateSample(mywalker);

}

//----------------------------------------------------------------------
int Rotate_orbs_method::showinfo(ostream & os) { 
  os << "Rotate_orbs method " << endl;
  return 1;
}


//----------------------------------------------------------------------
void Rotate_orbs_method::run(Program_options & options, ostream & output) {

  int norb=orbital_groups[0].GetDim(0);
  cout << "norb " << norb << endl;
  ifstream rot_inp("rotation");
  Array2 <doublevar> Rtot(norb,norb,0.0);
  for(int i=0; i< norb; i++) { 
    for(int j=0; j< norb; j++) { 
      rot_inp >> Rtot(i,j);
    }
  }
  ofstream testorb("test.orb");
  mymomat->writeorb(testorb, Rtot,orbital_groups[0]);
  testorb.close();
  
}
