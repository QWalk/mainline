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
#ifndef PROPERTIES_POINT_INCLUDED
#define PROPERTIES_POINT_INCLUDED

#include "Qmc_std.h"
#include "Wavefunction.h"
#include "Average_generator.h"
/*!
  Universal quantities at one point.  This includes the wave function value,
  any averaging variables, and the various energy components of the configuration.
 */
struct Properties_point {

  Properties_point() {
    reset();
  }
  
  void setSize(int nwf);
  void reset() {
    weight=0;
    count=0;
  }
  void mpiSend(int node);
  void mpiReceive(int node);
  void write(string & indent, ostream & os);
  void read(istream & is);

  const doublevar  energy(int w) const {
    return kinetic(w)+potential(w)+nonlocal(w);
  }
  //Add the values contained in pt to this point, respecting the weight
  void weighted_add(const Properties_point & pt);
  //Add the values contained in pt, without respecting the weight, and 
  //with an optional prefactor
  void unweighted_add(const Properties_point & pt,doublevar pre=1.0);
 
  int count; //whether to count this point or not

  //Properties to track
  Array1 <doublevar> kinetic;
  Array1 <doublevar> potential;
  Array1 <doublevar> nonlocal;
  Array1 <doublevar> weight; //!< averaging weight
  Wf_return wf_val; //!< wavefunction value  
  Array2 <Average_return> avgrets; 
  //general accumulations from Average_generators.  First index is the wf/sys number
  // and the second one is the Average_generator number
  
  private:
};




#endif //PROPERTIES_POINT_INCLUDED
