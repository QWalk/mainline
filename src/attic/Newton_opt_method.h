/*
 
Copyright (C) 2007 Michal Bajdich

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


#ifndef NEWTON_OPT_METHOD_H_INCLUDED
#define NEWTON_OPT_METHOD_H_INCLUDED

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
Optimize the wavefunction parameters with a given configuration.
Keyword: NEWTON_OPT
 
*/
class Newton_opt_method : public Qmc_method
{
public:
  void read(vector <string> words,
            unsigned int & pos,
            Program_options & options);

  void readcheck(string & readconfig);

  void run(Program_options & options, ostream & output);
  ~Newton_opt_method()
  {
    if(pseudo) delete pseudo;
    if(sysprop) delete sysprop;

    deallocate(wfdata);
    deallocate(wf);
    delete sample;
  }

  int showinfo(ostream & os);
private:
  void generate_all_quanities (Array1 <doublevar> & parms,
			       Array1 <doublevar> & delta, 
			       Array1 <doublevar> & local_energy,
			       Array1 < Array1 <doublevar> > & local_energy_gradient,
			       Array1 < Array1 <doublevar> > & wf_gradient,
			       Array1 < Array1 <doublevar> > & wf_hessian
			       );

 
  void calculate_first_averages(Array1 <doublevar> & parms,
				Array1 <doublevar> & local_energy, 
				Array1 < Array1 <doublevar> > local_energy_gradient,
				Array1 < Array1 <doublevar> > wf_gradient,
				doublevar & function_mean,
				doublevar & function_mean_est,
				doublevar & energy_mean_est, 
				doublevar & variance_mean_est,
				Array1 <doublevar> & wf_gradient_mean,
				Array1 <doublevar> & local_energy_gradient_mean,
				Array1 <doublevar> & energy_gradient_mean
				);
  
  void calculate_second_averages(Array1 <doublevar> & parms, 
				 Array1 <doublevar> & local_energy, 
				 Array1 < Array1 <doublevar> > local_energy_gradient,
				 Array1 < Array1 <doublevar> > wf_gradient,
				 Array1 < Array1 <doublevar> > wf_hessian,
				 doublevar & energy_mean,
				 Array1 <doublevar> & wf_gradient_mean,
				 Array1 <doublevar> & local_energy_gradient_mean,
				 Array1 <doublevar> & energy_gradient_mean,
				 Array1 <doublevar> & variance_grad_mean,
				 Array2 <doublevar> & hess1_eng_mean,
				 Array2 <doublevar> & hess2_eng_mean,
				 Array2 <doublevar> & hess3_1eng_mean,
				 Array2 <doublevar> & hess3_2eng_mean,
				 Array2 <doublevar> & hess1_var_mean
				 );


  void build_gradient(Array1 <doublevar> & energy_gradient_mean,
		      Array1 <doublevar> & variance_grad_mean,
		      Array1 <doublevar> & grad_var
		      );

  void build_hessian(Array1 <doublevar> & wf_gradient_mean,
		     Array1 <doublevar> & local_energy_gradient_mean,
		     Array2 <doublevar> & hess1_eng_mean,
		     Array2 <doublevar> & hess2_eng_mean,
		     Array2 <doublevar> & hess3_1eng_mean,
		     Array2 <doublevar> & hess3_2eng_mean,
		     Array2 <doublevar> & hess1_var_mean,
		     Array2 <doublevar> & hessian
		     );
  
  void Get_correlated_energies_and_variances(Array1 < Array1 <doublevar> > & parms,
					     int & iter,
					     Array1 <doublevar> & y, 
					     Array1 <doublevar> & energies,
					     Array1 <doublevar> & variances,
					     ostream & output,
					     Program_options & options);
  void adjust_distribution(Array1 <doublevar> & parms, 
			   int & iter,
			   string log_label,
			   ostream & output, 
			   Program_options & options,
			   doublevar & vmc_function,
			   doublevar & vmc_energy,
			   doublevar & vmc_energy_err,
			   doublevar & vmc_variance,
			   doublevar & weight_variance
			   );

  void wf_printout(Array1 <doublevar> & parms, int iter, doublevar & value, doublevar & energy, doublevar & energy_err, doublevar & variance, doublevar & weight_variance, int nconfig, doublevar mu, 
		   ostream & output);

  int Levenberg_marquad(Array1 <doublevar> & gradient, 
					   Array2 <doublevar> & hessian, 
					   int & iter, 
					   Array1 <doublevar> & parms, 
					   doublevar & function_mean,
					   doublevar & energy_mean,
					   doublevar & energy_mean_err,
					   doublevar & variance_mean,
					   int & nu, 
					   doublevar & damping,
					   ostream & vmcoutput,
					   ostream & output,
					   Program_options & options);
  
  doublevar eref; //!< reference energy for variance minimization
  int eref_exists; //!< whether the reference energy for variance minimization was supplied
  int nconfig;   //!< Number of configurations(walkers)
  int vmc_nblocks; //!< Number of blocks in vmc for distribution adjustment
  int iterations; //!< Number of times that the variance can be called.  The variance routine is called about once per parameter per iteration
  int nfunctions; //!< Number of wavefunctions we have.
  enum min_function_type { min_variance, min_energy, min_mixed, min_weight_variance } min_function;
  int use_weights; //!< Whether to use weights in the correlated sampling
  doublevar mixing; //!< mixing weight for energy component in mixed minimization (0.95 default)
  int plus_version_of_hessian_of_energy;
  int use_correlated_sampling;
  int store_hessian_to_file;

  string readconfig;
  string readconfig_non_cannonical;
  string storeconfig;
  string storeconfig_non_cannonical;
  vector < vector <string> > mc_words;

  int nparms;
  Array1 <doublevar> calcpot;
  string hessianstore; //!< Where to put the temporary pseudo file
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
  int  calculate_projectors;
  int nblocks_max;
  Array1 <doublevar> delta_parms;
  Array1 <doublevar> old_delta_parms;



};

#endif //NEWTON_OPT_METHOD_H_INCLUDED
//------------------------------------------------------------------------
