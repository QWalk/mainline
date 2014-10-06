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


#ifndef WANNIER_METHOD_H_INCLUDED
#define WANNIER_METHOD_H_INCLUDED

#include <vector>
#include <fstream>
#include <math.h>
#include <stdlib.h>
#include "Qmc_method.h"
#include "qmc_io.h"
#include "Array.h"
#include "MO_matrix.h"
#include "Sample_point.h"
#include "Array45.h"
class System;
class Wavefunction;
class Program_options;

/*!
 */
class Wannier_method : public Qmc_method
{
public:
  void read(vector <string> words,
            unsigned int & pos,
            Program_options & options);
  void run(Program_options & options, ostream & output);
  int showinfo(ostream & os);
  Wannier_method():mymomat(NULL), mywalker(NULL),sys(NULL) {}
  ~Wannier_method()
  {
    if(mymomat!=NULL) delete mymomat;
    if(mywalker != NULL) delete mywalker;
    if(sys != NULL) delete sys;
  }

private:
  void calculate_overlap(Array1 <int> &  orb_list, Array3 <dcomplex> & eikr, Array2 <doublevar> & phi2phi2);
  void optimize_rotation(Array3 <dcomplex> &  eikr,Array2 <doublevar> & R ); 
  
  doublevar evaluate_local(const Array3 <dcomplex> & eikr,
    Array2 <doublevar> & Rgen, Array2 <doublevar> & R)  ;
  doublevar eval_tstep(Array3 <dcomplex> & eikr, Array2 <doublevar> & Rgen,
    Array2 <doublevar> & Rgen_save, Array2 <doublevar> & deriv, doublevar tstep,
    Array2 <doublevar> & R);
  
  Array1 < Array1  <int> > orbital_groups;	//groups of orbitals to localize together
  doublevar resolution;	//grid coarsness: 10=coarser 0.1=finer
  MO_matrix * mymomat; //will hold MO information
  Sample_point * mywalker; //a single configuration/walker
  System * sys;
  Array2 <doublevar> mymovals; //(i,j) where i=MO#, j=0 default (for now)
  Array2 <doublevar> LatticeVec;
  Array1 <doublevar> origin;
  doublevar shake;
  string out_orbs;


};

#endif //WANNIER_METHOD_H_INCLUDED
//------------------------------------------------------------------------
