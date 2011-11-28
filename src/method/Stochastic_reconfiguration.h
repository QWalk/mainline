/*
 
Copyright (C) 2011 Lucas K. Wagner

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

#ifndef STOCHASTIC_RECONFIGURATION_H_INCLUDED
#define STOCHASTIC_RECONFIGURATION_H_INCLUDED

#include "Qmc_std.h"
#include "Qmc_method.h"
#include "Wavefunction.h"
#include "Wavefunction_data.h"
#include "Sample_point.h"
#include "Guiding_function.h"
#include "System.h"
#include "Pseudopotential.h"
#include "Program_options.h"

class Stochastic_reconfiguration_method : public Qmc_method
{
public:
  void read(vector <string> words,
            unsigned int & pos,
            Program_options & options);
  void run(Program_options & options, ostream & output);
  ~Stochastic_reconfiguration_method()
  {
    if(pseudo) delete pseudo;
    if(sys) delete sys;
    if(wfdata) delete wfdata;
  }

  Stochastic_reconfiguration_method() { 
    wfdata=NULL; sys=NULL;
    pseudo=NULL;
  }
  int showinfo(ostream & os);
private:
  doublevar tau;
  int iterations;
  int vmc_nstep;   
  Wavefunction_data * wfdata;
  Pseudopotential * pseudo;
  System * sys;
  Primary guide_wf; 
  string wfoutputfile;
  Program_options options;

  void wavefunction_derivative(Array1<doublevar> & energies,Array2 <doublevar> & S,Array1<doublevar> & en);

};

#endif //STOCHASTIC_RECONFIGURATION_H_INCLUDED

