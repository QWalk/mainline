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

class Excitation_list { 
  public:
     void build_excitation_list(Array3 <Array1 <int> > & occupation,int f);//(function,det,spin) (orb #) )
     template <class T> void clark_updates(Array2 <T> & ginv, Array2 <T> & M,
         int s, Array1 <T> & ratios);
  private:
     Array1 <Excitation> excitations;
     Array1 <Excitation> remap;
     vector <vector <int> > allg,alle;
};

//---------

template <class T> void Excitation_list::clark_updates(Array2 <T> & ginv, Array2 <T> & M,
     int s, Array1 <T> & ratios) { 
  int norb=M.GetDim(1);
  int ne=M.GetDim(0);
  int nex=excitations.GetDim(0);

  Array2 <T>  tmat,detmat;

  tmat.Resize(allg[s].size(),alle[s].size());
  int counti=0;
  for(vector<int>::iterator gi=allg[s].begin(); gi!=allg[s].end(); gi++) { 
    int countj=0;
    for(vector <int>::iterator ei=alle[s].begin(); ei!=alle[s].end(); ei++) { 
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
    
  ratios.Resize(nex);
  ratios=T(1.0);
  
  for(int ex=1; ex < nex; ex++) { 
    Array1 <int> & g=remap(ex).g(s);
    Array1 <int> & e=remap(ex).e(s);
    int n=excitations(ex).g(s).GetDim(0);
    switch(n) { 
      case 0:
        break;
      case 1:
        ratios(ex)=T(excitations(ex).sign(s))*tmat(g(0),e(0));
        break;
      case 2:
        ratios(ex)=T(excitations(ex).sign(s))* (tmat(g(0),e(0))*tmat(g(1),e(1))
              - tmat(g(1),e(0))*tmat(g(0),e(1)));
        break;
      default:
        detmat.Resize(n,n);
        for(int i=0; i < n; i++ ) { 
          for(int j=0; j < n; j++) { 
            detmat(i,j)=tmat(g(i),e(j));
          }
        }
        ratios(ex)=Determinant(detmat,n)*T(excitations(ex).sign(s));
        break;
    }
    //if(n==0) { }  //do nothing
    //else if(n==1) { 
    //  ratios(ex)=T(excitations(ex).sign(s))*tmat(g(0),e(0));
    //}
    //else  { 
    //  detmat.Resize(n,n);
    //  for(int i=0; i < n; i++ ) { 
    //    for(int j=0; j < n; j++) { 
    //      detmat(i,j)=tmat(g(i),e(j));
    //    }
    //  }
    //  ratios(ex)=Determinant(detmat,n)*T(excitations(ex).sign(s));
    //}
  }

}

//--------
#endif //CLARK_UPDATES_H_INCLUDED
