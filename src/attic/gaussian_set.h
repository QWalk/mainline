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
#ifndef GAUSSIAN_SET_H_INCLUDED
#define GAUSSIAN_SET_H_INCLUDED


#include "Qmc_std.h"
//----------------------------------------------------------------------

class Contracted_gaussian {
public:
  Array1 <doublevar> alpha; //exponent
  Array1 <doublevar> coeff; //coefficients for each exponent
  Array1 <int> lvals; // symmetry in x^l(0) * y^l(1)  * z^l(2)
  string center_name;
  Array1 < int > center;

  Contracted_gaussian() {
    lvals.Resize(3);
  }

  Contracted_gaussian operator=(Contracted_gaussian & gauss) {
    center_name=gauss.center_name;
    for(int j=0; j< 3; j++) lvals(j)=gauss.lvals(j);
    int nexpand=gauss.alpha.GetDim(0);
    alpha.Resize(nexpand);
    coeff.Resize(nexpand);
    for(int j=0; j< nexpand; j++) {
      alpha(j)=gauss.alpha(j);
      coeff(j)=gauss.coeff(j);
    }
    return *this;
  }
};


//----------------------------------------------------------------------

class Center {
public:
  Array1 <doublevar> pos;
  string label;
  Array1 <int> basis;
  Center() {
    pos.Resize(3);
  }
};

void read_centerpos(string & filename, Array2 <doublevar> & position, vector <string> & labels);
void read_basis(string & basisin, Array1 <Contracted_gaussian> & basis);
void create_local_basis(string & centerin, string & basisin,
                        Array1 <Center> & centers,
                        Array1 <Contracted_gaussian >  & basis );
                        
//----------------------------------------------------------------------
/*!
A class that holds some convienent lookup tables.
*/
class Gaussian_lookups {
public:
  Array1 <int> totbasis2cen;
  Array1 <int> totbasis2bas;
  Array1 <int> center_start, center_end, equivalent_center;
  void set_lookup_tables(Array2 <doublevar> & latvec,
                         Array1 <doublevar> & origin,
                         Array1 <Center> & centers);
};

#endif //GAUSSIAN_SET_H_INCLUDED
