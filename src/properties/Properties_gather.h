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
//----------------------------------------------------------------------
//include/Properties_gather.h

#ifndef PROPERTIES_GATHER_H_INCLUDED
#define PROPERTIES_GATHER_H_INCLUDED

#include "Space_warper.h"
#include "Wavefunction.h"
#include "System.h"
#include "Sample_point.h"
#include "Properties_point.h"
#include "Pseudopotential.h"
#include "Guiding_function.h"
class Dynamics_generator;
class Dynamics_info;

class Properties_gather {
 public:

  Properties_gather() {
    //square_weight=0;
  }

  ~Properties_gather();
  void read(vector <string> & words);

  void gatherData(Properties_point &,
                  Pseudopotential *, System *, Wavefunction_data *, 
                  Wavefunction *, Sample_point *, Guiding_function *, int n_converge=1,
                  int aux_updated=0); 
  /*
   void extendedGather(Properties_point & myprop,
                                   Pseudopotential * psp, 
                                   System * sys, 
                                   Wavefunction_data * wfdata,
                                   Wavefunction * wf, 
                                   Sample_point * sample, 
                                   Guiding_function * guide, int n_converge,
				   int aux_updated,
				   Array2 <doublevar> & drift,
				   Array1 < Array2 <doublevar> > & aux_drift,
				       Array1 <Array2 <doublevar> > & aux_positions);
*/
/*
  void updateAuxFunctions(Sample_point * sample);

 
  int nAux() {
    return aux_sys.GetDim(0);
  }

  void squareWeight(int s) {
    assert(s==0 || s==1);
    square_weight=s;
  }  
*/
 private:
/*
  int square_weight;

  //Auxillary wave function data
  Space_warper warper;
  Array1 <System * > aux_sys;
  Array1 <Wavefunction_data *> aux_wfdata;
  Array1 <Wavefunction * > aux_wf;
  Array1 <Sample_point * > aux_sample;
  */

};



#endif //PROPERTIES_GATHER_H_INCLUDED

//----------------------------------------------------------------------
