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


#ifndef REPTATION_METHOD_H_INCLUDED
#define REPTATION_METHOD_H_INCLUDED

#include "Qmc_std.h"
#include "Qmc_method.h"
#include "Wavefunction.h"
#include "Wavefunction_data.h"
#include "Sample_point.h"
#include "Guiding_function.h"
#include "System.h"
#include "Pseudopotential.h"
#include "Split_sample.h"
#include "Space_warper.h"
class Program_options;
#include "Properties.h"
#include <deque>

struct Reptile_point;

/*!
\brief
Evaluates the expectation value \f$ <\Psi|H|\Psi>/<\Psi|\Psi> \f$
stochastically.  Keyword: REPTATION 
*/
class Reptation_method : public Qmc_avg_method
{
public:

  Reptation_method() {
    have_read_options=0;
    have_generated_variables=0;
    have_attached_variables=0;
    sys=NULL;
    mywfdata=NULL;
    guidewf=NULL;
    pseudo=NULL;
    wf=NULL;
  }
  void read(vector <string> words,
            unsigned int & pos,
            Program_options & options);
  void run(Program_options & options, ostream & output);
  
  int showinfo(ostream & os);

  int generateVariables(Program_options & options);

  virtual void runWithVariables(Properties_manager & prop, 
                                System * sys,
                                Wavefunction_data * wfdata,
                                Pseudopotential * psp,
                                ostream & output);

  ~Reptation_method()
  {
    
    if(have_generated_variables) {
      if(sys) delete sys;
      if(mywfdata) delete mywfdata;
      if(pseudo) delete pseudo;
    }

    if(guidewf) delete guidewf;

    for(int i=0; i< densplt.GetDim(0); i++) {
      if(densplt(i)) delete densplt(i);
    }
    
   }

private:
  int allocateIntermediateVariables(System * , Wavefunction_data *);
  int deallocateIntermediateVariables();
  
  //return 1 if we read an entire reptile, 0 if not
  int readcheck(string & , int & direction, deque <Reptile_point> & reptile);
  void storecheck(int direction, deque <Reptile_point> & reptile,string &);

  doublevar getAcceptance(deque <Reptile_point> & reptile,
                          Reptile_point & pt, int direction);
  doublevar slither(int direction, 
                    deque <Reptile_point> & reptile,
                    Properties_gather & mygather, 
                    Reptile_point & pt,
                    doublevar & main_diffusion);

  void get_avg(deque <Reptile_point> & reptile, 
	       Properties_point & pt);
  void get_center_avg(deque <Reptile_point> & reptile, 
		      Properties_point & pt);
  int have_read_options;
  int have_generated_variables;
  int have_attached_variables;

  string center_trace;
  int trace_wait;

  
  doublevar energy_cutoff; //!< when to cut off the branching term
  doublevar eref; //!< reference energy



  int nblock;
  int nstep;
  string readconfig;
  string log_label;
  string storeconfig;
  doublevar timestep;
  int reptile_length;

  Properties_gather mygather;
  Dynamics_generator * sampler;
  Guiding_function * guidewf;
  Pseudopotential * pseudo;
  System * sys;
  Sample_point *  sample; 
  Wavefunction * wf; 
  Wavefunction_data * mywfdata;


  Array1 < Local_density_accumulator *> densplt;
  vector <vector <string> > dens_words;
  Array1 < Average_generator * > average_var;
  vector <vector <string> > avg_words;

  
};

#endif //REPTATION_METHOD_H_INCLUDED
//------------------------------------------------------------------------
