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
#ifndef AVERAGE_EKT_H_INCLUDED
#define AVERAGE_EKT_H_INCLUDED
#include "Average_generator.h"
#include "MO_matrix.h"
#include "Pseudopotential.h"
#include <iostream>
using namespace std; 
//Evaluate extended Koopmans' theorem.
class Average_ekt:public Average_generator { 
 public:
  virtual void evaluate(Wavefunction_data * wfdata, Wavefunction * wf,
                        System * sys, Sample_point * sample, Average_return & ) {
    cout << "need psp object" << endl;
  };
  virtual void evaluate(Wavefunction_data * wfdata, Wavefunction * wf,
                        System * sys, Pseudopotential *psp, Sample_point * sample, Average_return & );
  virtual void evaluate(Wavefunction_data * wfdata, Wavefunction * wf,
		System * sys, Pseudopotential *psp, Sample_point * sample, Properties_point & pt, Average_return & );
  virtual void read(System * sys, Wavefunction_data * wfdata, vector
                    <string> & words);
  virtual void write_init(string & indent, ostream & os);
  virtual void randomize(Wavefunction_data * wfdata, Wavefunction * wf,
			 System * sys, Sample_point * sample);
  virtual void read(vector <string> & words);
  virtual void write_summary(Average_return &,Average_return &, ostream & os);
  virtual void joutput(Array2 <dcomplex> &obdm, ostream & os);
  virtual void jsonOutput(Average_return &,Average_return &, ostream & os);
  virtual ~Average_ekt() { 
    if(momat) delete momat;
    if(cmomat) delete cmomat;
  }
  Average_ekt() { 
    momat=NULL;
    cmomat=NULL;
  }

  //		  bool do_tmoves,vector <Tmove> & tmoves); 
 private:
  MO_matrix * momat;
  Complex_MO_matrix * cmomat;
  int nmo;
  int npoints_eval; //Number of integration points to use
  int nstep_sample; //Number of steps to use when sampling the auxilliary points
  int deterministic_psp; //Whether to use deterministic evaluation of the pseudopotential
  doublevar gen_sample(int nstep, doublevar  tstep, int e, Array2 <dcomplex> & movals, Sample_point * sample) ;
  void calc_mos(Sample_point *, int e, Array2 <dcomplex> & movals);
  void calc_mosLap(Sample_point * sample, int e, Array2 <dcomplex> & molaps); 

  void evaluate_valence(Wavefunction_data * wfdata, Wavefunction * wf,
                        System * sys, Pseudopotential *psp, 
			Sample_point * sample, Average_return & avg);
  void evaluate_obdm(Wavefunction_data * wfdata, Wavefunction * wf,
		     System * sys, Sample_point * sample, 
		     Average_return & avg);
  //  void evaluate_ekt();
  void evaluate_conduction(Wavefunction_data * wfdata, Wavefunction * wf,
			   System * sys, Pseudopotential *psp, 
			   Sample_point * sample, Average_return & avg); 
  void calcPseudoMo(System * sys,
		    Sample_point * sample,
		    Pseudopotential *psp, 
		    const Array1 <doublevar> & accept_var,
		    Array1 <dcomplex> & totalv);//, 
  Array1 <Array1<int> > occupations;
  Array1 < Array1 <doublevar> > saved_r;
  Array1 < int > rk; //the electron number that each position corresponds to
  bool complex_orbitals; 
  bool eval_conduction;
  bool eval_valence;
  bool eval_obdm; 
  int totnelectrons;
  bool dump_data;
  doublevar ekt_cutoff; //!< The cutoff for EKT values; if |ekt| > cutoff, value is set to sign(ekt)*cutoff
};

#endif //AVERAGE_EKT_H_INCLUDED
