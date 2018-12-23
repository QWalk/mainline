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

#ifndef LINEAR_OPTIMIZATION_H_INCLUDED
#define LINEAR_OPTIMIZATION_H_INCLUDED

#include "Qmc_std.h"
#include "Qmc_method.h"
#include "Wavefunction.h"
#include "Wavefunction_data.h"
#include "Sample_point.h"
#include "Guiding_function.h"
#include "System.h"
#include "Pseudopotential.h"
#include "Program_options.h"
struct Average_return;
class Linear_optimization_method : public Qmc_method
{
public:
  void read(vector <string> words,
            unsigned int & pos,
            Program_options & options);
  void run(Program_options & options, ostream & output);
  ~Linear_optimization_method()
  {
    if(pseudo) delete pseudo;
    if(sys) delete sys;
    if(wfdata) delete wfdata;
  }

  Linear_optimization_method() { 
    wfdata=NULL; sys=NULL;
    pseudo=NULL;
  }
  int showinfo(ostream & os);
private:
  int iterations;
  int vmc_nstep;   
  int nconfig_eval;
  int max_nconfig_eval;
  int max_vmc_nstep;
  int max_zero_iterations;
  double svd_tolerance;
  double nodal_cutoff; // Cutoff for the parameter derivatives
  bool do_uncorrelated_evaluation;
  bool pseudopotential_derivatives;
  doublevar sig_H_threshold;
  doublevar en_convergence;
  doublevar minimum_psi0;
  Wavefunction_data * wfdata;
  Pseudopotential * pseudo;
  System * sys;
  Primary guide_wf; 
  string wfoutputfile;
  Program_options options;

  void wavefunction_derivative(Array2<doublevar> & H,Array2 <doublevar> & S,Array1<doublevar> & en);
  void wavefunction_energy(Array1 <doublevar> & en);
  
  double fit_stabil(Array1 <doublevar> & stabil, Array2 <doublevar> & energies,int nfit=5);
  
  //returns the estimated energy change
  double line_minimization(Array2 <doublevar> & S, 
                           Array2 <doublevar> & Sinv, 
                           Array2 <doublevar> & H,
                           Array1 <doublevar> & alpha,
                           ostream & os);

  void correlated_evaluation(Array1 <Array1 <doublevar> > & alphas,
                             int ref_alpha,
                             Array2 <doublevar> & energies);

  void uncorrelated_evaluation(Array1 <Array1 <doublevar> > & alphas,
                               Array2 <doublevar> & energies);
  
  bool deriv_is_significant(Average_return & avg, Average_return & err,int n);
  
};




#endif //LINEAR_OPTIMIZATION_H_INCLUDED

