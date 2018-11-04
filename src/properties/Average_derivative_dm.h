/*
 
 Copyright (C) 2017 Lucas K. Wagner
 
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

#ifndef AVERAGE_DERIVATIVE_DM_H_INCLUDED
#define AVERAGE_DERIVATIVE_DM_H_INCLUDED

#include "Average_generator.h"
#include "Average_density_matrix.h"

/*!
 \brief
 Evaluate the correlation of moments of the Hamiltonian with Average_generator objects.
 */
class Average_derivative_dm:public Average_generator {
 public:
  virtual void evaluate(Wavefunction_data * wfdata, Wavefunction * wf,
                        System * sys, Sample_point * sample, Average_return & ) { 
    error("need to use Properties_point version of evaluate");
  }
  virtual void evaluate(Wavefunction_data * wfdata, Wavefunction * wf,
                        System * sys, Sample_point * sample, Properties_point & pt,
                        Average_return &avg);
  
  virtual void read(System * sys, Wavefunction_data * wfdata, vector
                    <string> & words);
  virtual void write_init(string & indent, ostream & os);
  virtual void randomize(Wavefunction_data * wfdata, Wavefunction * wf,
                        System * sys, Sample_point * sample) { 
    dm_eval.randomize(wfdata,wf,sys,sample);
  } 
  virtual void read(vector <string> & words);
  virtual void write_summary(Average_return &,Average_return &, ostream & os);
  virtual void jsonOutput(Average_return &,Average_return &, ostream & os);
 private:
  Average_tbdm_basis dm_eval;
  int nparms,nmo,ndm_elements;
  //doublevar cutoff; //Cutoff for evaluating wf gradient wrt det coeff: distance per electron
  Array1<doublevar> cutoff_list; //List of cutoffs to use when evaling gradients; distance per electron
};




#endif //AVERAGE_DERIVATIVE_DM_H_INCLUDED
