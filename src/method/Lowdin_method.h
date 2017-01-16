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


#ifndef LOWDIN_METHOD_H_INCLUDED
#define LOWDIN_METHOD_H_INCLUDED

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
class Lowdin_method : public Qmc_method
{
public:
  void read(vector <string> words,
            unsigned int & pos,
            Program_options & options);
  void run(Program_options & options, ostream & output);
  int showinfo(ostream & os);
  Lowdin_method():mymomat(NULL), mywalker(NULL),sys(NULL) {}
  ~Lowdin_method()
  {
    if(mymomat!=NULL) delete mymomat;
    if(mywalker != NULL) delete mywalker;
    if(sys != NULL) delete sys;
  }

private:
  void calculate_overlap(Array2 <doublevar> & S);
  doublevar resolution;	//grid coarsness: 10=coarser 0.1=finer
  MO_matrix * mymomat; //will hold MO information
  Sample_point * mywalker; //a single configuration/walker
  System * sys;
  Array2 <doublevar> LatticeVec;
  Array1 <doublevar> origin;
  string out_orbs;
};

#endif //LOWDIN_METHOD_H_INCLUDED
//------------------------------------------------------------------------
