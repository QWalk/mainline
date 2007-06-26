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

#ifndef SPACE_WARPER_H_INCLUDED
#define SPACE_WARPER_H_INCLUDED

#include "Qmc_std.h"
#include "Basis_function.h"
class Sample_point;

class Space_warper {
public:
  Space_warper();

  ~Space_warper() { 
    //if(weight_basis) delete weight_basis; 
  }

  void set_warp(int warp) {
    assert(warp==0 || warp==1);
    warp_on=warp;
  }

  void read(vector <string> & words);
  int space_warp(Sample_point * refsample,
                 Sample_point * sample,
                 int e, Array1 <doublevar> & R_old,
                 Array1 <doublevar> & R_new, 
                 doublevar & jacobian );
private:
  Basis_function * weight_basis;
  int warp_on;
  int ex; //!< exponent of warping transformation
};

#endif //SPACE_WARPER_H_INCLUDED

//----------------------------------------------------------------------
