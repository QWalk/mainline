/*
 
Copyright (C) 2009 Michal Bajdich and Fernando A. Reboredo

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


#ifndef SHDMC_METHOD_H_INCLUDED
#define SHDMC_METHOD_H_INCLUDED

#include "Qmc_std.h"
#include "Qmc_method.h"
#include "Wavefunction.h"
#include "Wavefunction_data.h"
#include "Sample_point.h"
#include "Guiding_function.h"
#include "System.h"
#include "Pseudopotential.h"
#include "Properties.h"

class Program_options;

/*!
\brief
Keyword: SHDMC
Optimize the wave function lambda or dlambda parameters in DMC method.
At the end, the wave function will be the best (maximal overlap) FN wave function 
for given parameter space. See the SHDMC paper of MB and FAR, PR? 
*/

class Shdmc_method : public Qmc_method
{
public:
  void read(vector <string> words,
            unsigned int & pos,
            Program_options & options);

  void readcheck(string & readconfig);

  void run(Program_options & options, ostream & output);
  ~Shdmc_method()
  {
    if(pseudo) delete pseudo;
    if(sysprop) delete sysprop;

    deallocate(wfdata);
    deallocate(wf);
    delete sample;
  }

  int showinfo(ostream & os);
private:
  void get_averages(Array1 <doublevar> & parms, 
		    int & iter,
		    ostream & output, 
		    Program_options & options,
		    doublevar & vmc_function,
		    doublevar & vmc_energy,
		    doublevar & vmc_energy_err,
		    doublevar & vmc_variance,
		    doublevar & weight_variance,
		    Array1 <doublevar> & new_parms);

  void wf_printout(Array1 <doublevar> & parms, int iter, 
		   doublevar & value, doublevar & energy, doublevar & energy_err, doublevar & variance, doublevar & weight_variance, 
		   ostream & output);
  
  doublevar eref; //!< reference energy for variance minimization
  int eref_exists; //!< whether the reference energy for variance minimization was supplied
  int nconfig;   //!< Number of configurations(walkers)
  int vmc_nblocks; //!< Number of blocks in vmc for distribution adjustment
  int iterations; //!< Number of times that the variance can be called.  The variance routine is called about once per parameter per iteration
  int nfunctions; //!< Number of wavefunctions we have.
  enum min_function_type { min_variance, min_energy, min_mixed, min_weight_variance } min_function;
  int use_weights; //!< Whether to use weights in the correlated sampling
  doublevar mixing; //!< mixing weight for energy component in mixed minimization (0.95 default)

  string readconfig;
  string readconfig_non_cannonical;
  string storeconfig;
  string storeconfig_non_cannonical;
  vector < vector <string> > mc_words;

  int nparms;
  Array1 <doublevar> calcpot;
  string wfoutputfile;//!< Where to put the wavefunction output
  Primary guide_wf; //!< Guiding function
    
  int use_extended_output; //!< Whether to print out wavefunction at every iter. of min
  System * sysprop;
  Sample_point * sample;
  Array1 <Config_save_point> config_pos;
  Array1 <doublevar> dmc_weight;
  Wavefunction * wf;  
  Wavefunction_data * wfdata;
  Pseudopotential * pseudo;
  Qmc_avg_method * avg_method;
  Properties_manager myprop;
  int nblocks_max;
  doublevar block_increase_factor;
  doublevar scale_step;
  Array1 <doublevar> delta_parms;
  Array1 <doublevar> old_delta_parms;
};

#endif //NEWTON_OPT_METHOD_H_INCLUDED
//------------------------------------------------------------------------

