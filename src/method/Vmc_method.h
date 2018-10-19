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

#ifndef VMC_METHOD_H_INCLUDED
#define VMC_METHOD_H_INCLUDED

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

/*!
\brief
Evaluates the expectation value \f$ <\Psi|H|\Psi>/<\Psi|\Psi> \f$
stochastically.  Keyword: VMC 
*/
class Vmc_method : public Qmc_avg_method
{
public:

  Vmc_method() {
    have_read_options=0;
    have_generated_variables=0;
    have_attached_variables=0;
    sysprop=NULL;
    mywfdata=NULL;
    guidewf=NULL;
    pseudo=NULL;
    wf=NULL;
    sample=NULL;
    sampler=NULL;
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

  ~Vmc_method()
  {
    
    if(have_generated_variables) {
      if(sysprop) delete sysprop;
      if(mywfdata) delete mywfdata;
      if(pseudo) delete pseudo;
    }

    if(guidewf) delete guidewf;
    if(sampler) delete sampler;

    for(int i=0; i< densplt.GetDim(0); i++) {
      if(densplt(i)) delete densplt(i);
    }

    
   }

private:
  int allocateIntermediateVariables(System * , Wavefunction_data *);
  int deallocateIntermediateVariables();
  void readcheck(string & );
  void storecheck(string &, int append=0);
  int have_read_options;
  int have_generated_variables;
  int have_attached_variables;

  Dynamics_generator * sampler;

  int nblock;
  int nstep;
  int nconfig;
  int ndecorr;
  bool auto_timestep;
  string storeconfig;
  string readconfig;
  string log_label;
  string jsonfile;
  string dump_file; //!< dump the energy and electron positions

  doublevar timestep;
  int nelectrons;

  int ndim;
  int print_wf_vals; //!< whether or not to put the values of the wavefunction to cout
  int low_io; //!< Only output the configurations and any Density objects at the end of the run, not every block
  string config_trace; //!< where to lay down a trace of configurations

  string guidetype;

  Guiding_function * guidewf;

  Pseudopotential * pseudo;
  System * sysprop;
  Sample_point * sample;
  Array1 <Config_save_point> config_pos;

  Wavefunction *   wf;  
  Wavefunction_data * mywfdata;

  Properties_manager myprop;

  Properties_gather mygather;

  Array1 < Local_density_accumulator *> densplt;
  vector <vector <string> > dens_words;
  Array1 < Nonlocal_density_accumulator *> nldensplt;
  vector <vector <string> > nldens_words;
  Array1 < Average_generator * > average_var;
  vector <vector <string> > avg_words;

};

#endif //VMC_METHOD_H_INCLUDED
//------------------------------------------------------------------------
