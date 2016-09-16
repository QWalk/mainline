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

#ifndef AVERAGE_DENSITY_MATRIX_H_INCLUDED
#define AVERAGE_DENSITY_MATRIX_H_INCLUDED

#include "Average_generator.h"
#include "MO_matrix.h"

//Evaluate the two-body density matrix on a basis.
class Average_tbdm_basis:public Average_generator { 
 public:
  virtual void evaluate(Wavefunction_data * wfdata, Wavefunction * wf,
                        System * sys, Sample_point * sample, Average_return & );
  virtual void read(System * sys, Wavefunction_data * wfdata, vector
                    <string> & words);
  virtual void write_init(string & indent, ostream & os);
  virtual void randomize(Wavefunction_data * wfdata, Wavefunction * wf,
                        System * sys, Sample_point * sample);
  virtual void read(vector <string> & words);
  virtual void write_summary(Average_return &,Average_return &, ostream & os);
  virtual void jsonOutput(Average_return &,Average_return &, ostream & os);
  virtual ~Average_tbdm_basis() { 
    if(momat) delete momat;
    if(cmomat) delete cmomat;
  }
  Average_tbdm_basis() { 
    momat=NULL;
    cmomat=NULL;
  }
 private:
  MO_matrix * momat;
  Complex_MO_matrix * cmomat;
  int nmo;
  int npoints_eval; //Number of integration points to use
  int nstep_sample; //Number of steps to use when sampling the auxilliary points
  doublevar gen_sample(int nstep, doublevar  tstep, int e, Array2 <dcomplex> & movals, Sample_point * sample) ;
  void calc_mos(Sample_point *, int e, Array2 <dcomplex> & movals);
  void evaluate_obdm(Wavefunction_data * wfdata, Wavefunction * wf,
                        System * sys, Sample_point * sample, Average_return & avg);

  //  void evaluate_ekt();
  void evaluate_tbdm(Wavefunction_data * wfdata, Wavefunction * wf,
                        System * sys, Sample_point * sample, Average_return & avg);
  void evaluate_old(Wavefunction_data * wfdata, Wavefunction * wf,
                        System * sys, Sample_point * sample, Average_return & avg);
  
  

  enum tbdm_t { tbdm_uu,tbdm_ud, tbdm_du,tbdm_dd } ;
  int tbdm_index(tbdm_t typ, int i, int j, int k, int l) { 
    return nmo+4*nmo*nmo+2*(typ*nmo*nmo*nmo*nmo+i*nmo*nmo*nmo+j*nmo*nmo+k*nmo+l);
  }

  void jsontbdmvalHelper(Average_return & avg,tbdm_t tbdm,ostream & os);
  void jsontbdmerrHelper(Average_return & avg,Average_return & err,tbdm_t tbdm,ostream & os);
  void jobdmoutput(Array2 <dcomplex> &obdm, ostream & os);
  void jdgoutput(string s, int indx,  Average_return &avg1, Average_return &avg2, ostream & os);
  
  Array1 <Array1 <int> > occupations;
  Array1 < Array1 <doublevar> > saved_r;
  Array1 < int > rk; //the electron number that each position corresponds to
  bool complex_orbitals; 
  bool eval_tbdm;
  bool eval_old;
  bool tbdm_diagonal;
  
};

#endif //AVERAGE_DENSITY_MATRIX_H_INCLUDED
