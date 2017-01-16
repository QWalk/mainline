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
//----------------------------------------------------------------------
//utils/create_mesh_orb.cpp

#include <iostream>
#include <cstdlib>

using namespace std;

int main(int argc, char ** argv) {
  if(argc < 2) {
    cerr << "usage: " << argv[0] << " <number of molecular orbitals>"
	 << endl;
    exit(1);
  }

  int nmo=atoi(argv[1]);

  //We output for one center the mo's
  //There are two coefficients: 1 is 0.0 and 2 is 1.0
  for(int m=0; m < nmo; m++) {
    for(int m2=0; m2 < nmo; m2++) {
      int coeff=1;
      if(m==m2) coeff=2;
      cout << m+1 << "    " << m2+1 << "   1    " << coeff << endl;
    }
  }

  cout << "COEFFICIENTS" << endl;
  cout << "0.0      1.0 " << endl;
  

}


//----------------------------------------------------------------------
