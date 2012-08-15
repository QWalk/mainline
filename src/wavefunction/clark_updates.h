/*
 
Copyright (C) 2012 Lucas K. Wagner

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
#ifndef CLARK_UPDATES_H_INCLUDED
#define CLARK_UPDATES_H_INCLUDED
#include "Array.h"

//
//
//----------------------------------------------------------------------
//Excitation management code
struct Excitation { 
  Array1< Array1 < int> > g; //remove these orbitals, indices (s,orb#)
  Array1< Array1 <int> >e; //replace with these orbitals, indices (s,orb#)
  Array1 <int> sign;  //Sign of the permutation compared to the input orbital ordering
};

void build_excitation_list(Array3 <Array1 <int> > & occupation,int f,//(function,det,spin) (orb #)
    Array1 <Excitation> & ex);
//---------

template <class T> void clark_updates(Array2 <T> & ginv, Array2 <T> & M,
    Array1 <Excitation> & excitations, int s, Array1 <T> & ratios) { 
  int norb=M.GetDim(1);
  int ne=M.GetDim(0);
  int nex=excitations.GetDim(0);
  /*
  cout << "ne " << ne << " norb " << norb << endl;

  cout << "[ ";
  for(int e=0; e< ne; e++) { 
    cout << "[ ";
    for(int o=0; o< norb; o++) { 
      cout << M(e,o)  <<  ", ";
    }
    cout << "],";
    cout << endl;
  }
  cout << "]";

  cout << "inverse " << endl;
  cout << "[ ";
  for(int e=0; e< ne; e++) { 
    cout << "[ ";
    for(int o=0; o< ne; o++) { 
      cout << ginv(e,o)  <<  ", ";
    }
    cout << "],";
    cout << endl;
  }
  cout << "]";
*/

  Array2 <T>  tmat,detmat;
  vector <int> allg, alle;
  for(int e=0; e< nex; e++) { 
    for(int i=0; i< excitations(e).g(s).GetDim(0); i++) { 
      int g=excitations(e).g(s)(i);
      bool found=false;
      for(vector<int>::iterator gi=allg.begin(); gi!=allg.end(); gi++) {
        if(g==*gi) {
          found=true;
          break;
        }
      }
      if(!found) allg.push_back(g);
      int ex=excitations(e).e(s)(i);
      found=false;
      for(vector<int>::iterator ei=alle.begin(); ei!=alle.end(); ei++) {
        if(ex==*ei) {
          found=true;
          break;
        }
      }
      if(!found) alle.push_back(ex);
    }
  }
//  for(vector<int>::iterator gi=allg.begin(); gi!=allg.end(); gi++)
//    cout << "g " << *gi << endl;
//  for(vector <int>::iterator ei=alle.begin(); ei!=alle.end(); ei++) 
//    cout << "e " << *ei << endl;


  tmat.Resize(allg.size(),alle.size());
  int counti=0;
  for(vector<int>::iterator gi=allg.begin(); gi!=allg.end(); gi++) { 
    int countj=0;
    for(vector <int>::iterator ei=alle.begin(); ei!=alle.end(); ei++) { 
      int i=*gi;
      int j=*ei;
      T dot=0.0;
      for(int e=0; e< ne; e++) {
        dot+=ginv(e,i)*M(e,j);
      }
//      cout << "tmat " << i << " " << j << " : " << dot << endl;
      tmat(counti,countj)=dot;
      countj++;
    }
    counti++;
  }
    
  Array1 <Excitation> remap(nex);
  for(int e=0; e< nex; e++) { 
    int n=excitations(e).g(s).GetDim(0);
    remap(e).g.Resize(2);
    remap(e).e.Resize(2);
    remap(e).g(s).Resize(n);
    remap(e).e(s).Resize(n);

    int ng=allg.size();
    int ne=alle.size();
    for(int i=0; i< n; i++) { 
      for(int j=0; j< ng; j++) { 
        if(allg[j]==excitations(e).g(s)(i)) {
          remap(e).g(s)(i)=j;
          break;
        }
      }
      for(int j=0; j< ne; j++) { 
        if(alle[j]==excitations(e).e(s)(i)) { 
          remap(e).e(s)(i)=j;
          break;
        }
      }
    }
  }

  ratios.Resize(nex);
  ratios=T(1.0);
  for(int ex=1; ex < nex; ex++) { 
    int n=excitations(ex).g(s).GetDim(0);
    detmat.Resize(n,n);
    for(int i=0; i < n; i++ ) { 
      for(int j=0; j < n; j++) { 
        detmat(i,j)=tmat(remap(ex).g(s)(i),
                         remap(ex).e(s)(j));
      }
    }
    ratios(ex)=Determinant(detmat,n)*T(excitations(ex).sign(s));
  }

}

//--------
#endif //CLARK_UPDATES_H_INCLUDED
