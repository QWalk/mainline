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
  All relevant quantities at one point
 */
struct Properties_point {

  Properties_point() {
    //There are some possible problems with hard-coding maxchildren,
    //but the only branching algorithm is DMC, and the way it is at
    //the time of writing, one can only branch once per step, so
    //it should really be safe at maxchildren=2.  Maybe we should
    //use a vector instead; need to check memory usage.
    maxchildren=3;
    children.Resize(maxchildren);
    reset();
  }
  
  void setSize(int nwf);
  void reset() {
    nchildren=0;
    parent=-1;
    weight=0;
    count=0;
    
  }


  void mpiSend(int node);

  void mpiReceive(int node);

  void write(string & indent, ostream & os);
  void read(istream & is);


  doublevar energy(int w) {
    return kinetic(w)+potential(w)+nonlocal(w);
  }

 
  int nchildren;
  int parent;
  Array1 <int> children;

  int count; //whether to count this point or not

  //Properties to track
  Array1 <doublevar> kinetic;
  Array1 <doublevar> potential;
  Array1 <doublevar> nonlocal;

  Array1 <doublevar> weight; //!< averaging weight
  Wf_return wf_val; //!< wavefunction value  
  Array1 <Average_return> avgrets; //general accumulations from Average_generators
  
  private:
  int maxchildren;
};




#endif //PROPERTIES_POINT_INCLUDED
