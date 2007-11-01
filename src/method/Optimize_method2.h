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


#ifndef OPTIMIZE_METHOD2_H_INCLUDED
#define OPTIMIZE_METHOD2_H_INCLUDED

#include "Qmc_std.h"
#include "Qmc_method.h"
#include "Wavefunction.h"
#include "Wavefunction_data.h"
#include "Sample_point.h"
#include "Guiding_function.h"
#include "System.h"
#include "Pseudopotential.h"
class Program_options;

/*!
\brief
Optimize the wavefunction parameters with a given configuration.
Keyword: OPTIMIZE
 
*/
class Optimize_method2 : public Qmc_method
{
public:
  void read(vector <string> words,
            unsigned int & pos,
            Program_options & options);
  void run(Program_options & options, ostream & output);
  ~Optimize_method2()
  {
    if(pseudo) delete pseudo;
    if(sysprop) delete sysprop;

    deallocate(wfdata);
    for(int i=0; i< nconfig; i++)
    {
      deallocate(wf(i));
      delete electrons(i);
    }

  }

  int showinfo(ostream & os);
private:
  void func_val(int n, const Array1 <double> & parms, double & val, double & energy, double & variance, 
                int & min_nconfig, ostream & output);

  void energy_grad(Array1 <double> & parms, int nparms_start, int nparms_end, Array1 <double> & grad, doublevar & ereff, 
                                     Array1 <double> & delta, int & min_nconfig, ostream & output);

  void func_hessian_rev(Array1 <double> & parms, int nparms_start, int nparms_end, Array2 <double> & hessian,  
                    Array1 <doublevar>  & grad_var_tmp,  doublevar & ereff,
		    Array1 <double> & delta, int & min_nconfig, ostream & output);

  void energy_grad_analytical(Array1 <double> & parms, int nparms_start, int nparms_end, Array1 <double> & grad,  
						Array1 <double> & grad_wf, doublevar & energy_mean,
						Array1 <double> & delta, int & min_nconfig, ostream & output);
  void func_hessian_rev_analytical(Array1 <double> & parms, int nparms_start, int nparms_end, Array2 <double> & hessian,
                                    Array1 <doublevar>  & grad_var, doublevar & energy_mean, 
						     Array1 <double> & delta, int & min_nconfig, ostream & output);
  int LEVMAR_DER(Array1 <double> & parms, int nparms_start, int nparms_end, Array1 <double> & delta,
                 int & min_nconfig, int iter_min, 
                 const int itmax, int & iter, ostream & output);

  void wf_printout(int iter, doublevar value, doublevar energy, doublevar variance, int nconfig, doublevar mu, 
				     ostream & output);
  doublevar eref; //!< reference energy for variance minimization
  int eref_exists; //!< whether the reference energy for variance minimization was supplied
  int nconfig;   //!< Number of configurations(walkers)
  int iterations; //!< Number of times that the variance can be called.  The variance routine is called about once per parameter per iteration
  int nfunctions; //!< Number of wavefunctions we have.
  enum min_function_type { min_variance, min_energy, min_mixed } min_function;
  int use_weights; //!< Whether to use weights in the correlated sampling
  
  Array1 <Wf_return > orig_vals; //!< Original wave function values before opt
  doublevar ln_norm_orig_vals;
  string pseudostore; //!< Where to put the temporary pseudo file
  Pseudo_buffer psp_buff;
  string wfoutputfile;//!< Where to put the wavefunction output
  Primary guide_wf; //!< Guiding function
  Array1 <doublevar> local_energy; //!< local energy(that doesn't change when we optimize
  Array1 <doublevar> lastparms;
  doublevar mixing; //!< mixing weight for energy component in mixed minimization (0.95 default)
  int use_extended_output; //!< Whether to print out wavefunction at every iter. of min
  // int use_analytical_grad_and_hess; //!< Whether to use analytical gradient and hessian of WF with parms
  int iter_min_read;//!< number of iterations before the dumping is being decreased
  int min_nconfig_read;//!< starting number of configuration used
  int multiply;//!< by home much do I increase the nconfig
  int maxnparmsatonce; //!< maximum parameters optimized at once
  Array1 <doublevar> weight; //weights for each walker in optimization
  System * sysprop;
  Array1 <Sample_point *> electrons;
  Array1 <Wavefunction * > wf;
  Wavefunction_data * wfdata;
  Pseudopotential * pseudo;
};

#endif //OPTIMIZE_METHOD2_H_INCLUDED
//------------------------------------------------------------------------
