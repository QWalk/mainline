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

#ifndef OPTIMIZE_METHOD_H_INCLUDED
#define OPTIMIZE_METHOD_H_INCLUDED

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
class Optimize_method : public Qmc_method
{
public:
  void read(vector <string> words,
            unsigned int & pos,
            Program_options & options);
  void run(Program_options & options, ostream & output);
  ~Optimize_method()
  {
    if(pseudo) delete pseudo;
    if(sys) delete sys;
    if(wfdata) delete wfdata;
    if(wf) delete wf;
    if(sample) delete sample;
  }

  Optimize_method() { 
    sample=NULL; wf=NULL;
    wfdata=NULL; sys=NULL;
    pseudo=NULL;
  }
  int showinfo(ostream & os);
  doublevar variance(int n, Array1 <double> & parms, double & val, 
                     int check=1 //!< =1: make sure weights and variance make sense, =0: don't
                     );
  doublevar derivatives(int n, Array1 <double> & parms, Array1 <double> & deriv, 
                                         double & val, int check);
  void iteration_print(double f, double gg, double tol,  int itn, ostream & output);

private:

  doublevar eref; //!< reference energy
  int nconfig;   //!< Number of configurations(walkers)
  int iterations; //!< Number of optimization steps to take
  int nfunctions; //!< Number of wavefunctions we have.
  enum min_function_type { min_variance, min_abs, min_lorentz, min_energy, min_mixed } min_function;
  int use_weights; //!< Whether to use weights in the correlated sampling
  int guess_eref; //!< whether to guess eref from the incoming walker distribution
  int update_psp; //!< Whether to update the pseudopotential every iteration or not
  Array1 <Wf_return > orig_vals; //!< Original wave function values before opt
  string pseudostore; //!< Where to put the temporary pseudo file
  Pseudo_buffer psp_buff;
  string wfoutputfile;//!< Where to put the wavefunction output
  Primary guide_wf; //!< Guiding function
  Array1 <doublevar> local_energy; //!< local energy(that doesn't change when we optimize.  Includes the pseudopotential when update_psp is off.
  Array1 < Array1 <doublevar> >psp_test; //!< For random evaluation of the psp--we want to use the same random numbers
  Array1 <doublevar> lastparms;
  doublevar mixing; //!< mixing weight for energy component in mixed minimization (0.95 default)
  int use_extended_output; //!< Whether to print out wavefunction at every iter. of min
  System * sys;
  Array1 <Config_save_point> config_pos;
  Sample_point * sample;
  Wavefunction * wf;
  Wavefunction_data * wfdata;
  Pseudopotential * pseudo;
  
  Array1 <doublevar> last_ens;
};


#include "macopt.h"
//Interface to the macopt functions
class Optimize_fn:public Macopt {
public:
  Optimize_fn(int _n,int _verbose, double _tol,
              int _itmax,int _rich):Macopt(_n,_verbose,_tol,_itmax,_rich) { 
                param_n=_n; opt_method=NULL; output=NULL;
              } 
  ~Optimize_fn(){ };
  double func(double * _p);
  double dfunc(double * _p, double * _g);
  int param_n;
  void iteration_print(double f, double gg, double tol,  int itn) {
    assert(opt_method!=NULL);
    assert(output != NULL);
    opt_method->iteration_print(f,gg,tol,itn,*output);
  }

  Optimize_method * opt_method;
  ostream * output; //link to the output stream we should be using.
};


#endif //OPTIMIZE_METHOD_H_INCLUDED
//------------------------------------------------------------------------
