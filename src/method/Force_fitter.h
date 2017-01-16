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

#ifndef FORCE_FITTER_H_INCLUDED
#define FORCE_FITTER_H_INCLUDED


#include "Qmc_std.h"

class Force_fitter {
 public:
  void setup(doublevar  cut, int nexpansion );

  doublevar cutoff() { return cut;}
  void fit_force(Array1 <doublevar> & r, Array1 <doublevar> & bare_force,
                 Array1 <doublevar> & fit_force);
  


 private:
  int m;//the base expansion
  doublevar cut;
  int nexp;
  Array1 <doublevar> coeff;

};


#endif //FORCE_FITTER_H_INCLUDED
