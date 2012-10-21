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


#ifndef ENERGY_SLICE_H_INCLUDED
#define ENERGY_SLICE_H_INCLUDED

#include <vector>
#include <fstream>
#include <math.h>
#include <stdlib.h>
#include "Qmc_method.h"
#include "qmc_io.h"
#include "Array.h"
#include "MO_matrix.h"
#include "Sample_point.h"
#include "Average_generator.h"
class System;
class Wavefunction;
class Program_options;

/*!
 */
class Energy_slice : public Qmc_method
{
public:
  void read(vector <string> words,
            unsigned int & pos,
            Program_options & options);
  void run(Program_options & options, ostream & output);
  int showinfo(ostream & os);
  Energy_slice():sample(NULL),sys(NULL) {}
  ~Energy_slice()
  {
    if(sample != NULL) delete sample;
    if(sys != NULL) delete sys;
    if(wf !=NULL) delete wf;
    if(pseudo!=NULL) delete pseudo;
  }

private:
  Sample_point * sample; //a single configuration/walker
  System * sys;
  Wavefunction * wf;
  Wavefunction_data * wfdata;
  Pseudopotential * pseudo;
  //Array1 <Local_density_accumulator * > densplt;
  Array1 <Average_generator * > average_var;
  vector <vector <string > > avg_words;
  int nconfig;
  int nbin;
};

#endif //ENERGY_SLICE_H_INCLUDED
//------------------------------------------------------------------------
