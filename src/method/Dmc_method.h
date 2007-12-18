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


#ifndef DMC_METHOD_H_INCLUDED
#define DMC_METHOD_H_INCLUDED

#include "Qmc_std.h"
#include "Qmc_method.h"
#include "Wavefunction.h"
#include "Wavefunction_data.h"
#include "Sample_point.h"
#include "Guiding_function.h"
#include "Pseudopotential.h"
#include "System.h"
#include "Split_sample.h"
#include "Properties.h"
#include <deque>

class Program_options;

struct Dmc_history { 
  doublevar main_en;
  Array2 <doublevar> aux_en;
  void mpiSend(int node);
  
  void mpiReceive(int node);

};

struct Dmc_point { 
  Properties_point prop;
  deque <Dmc_history> past_energies;
  doublevar weight;
  int ignore_walker;
  Config_save_point config_pos;  
  Array1 <doublevar> age;  //!< age of each electron
  Dmc_point() { 
    weight=1;
    ignore_walker=0;
  }
  void mpiSend(int node);
  void mpiReceive(int node);
  
};


class Dmc_method : public Qmc_avg_method
{
public:


  void read(vector <string> words,
            unsigned int & pos,
            Program_options & options);
  int generateVariables(Program_options & options);
  void run(Program_options & options, ostream & output);

  virtual void runWithVariables(Properties_manager & prop, 
                                System * sys,
                                Wavefunction_data * wfdata,
                                Pseudopotential * psp,
                                ostream & output);


  Dmc_method() {
    have_read_options=0;
    have_allocated_variables=0;

    wf=NULL;
    mypseudo=NULL;
    mysys=NULL;
    mywfdata=NULL;
    guidingwf=NULL;
    dyngen=NULL;
    sample=NULL;
    
  }
  ~Dmc_method()
  {
    if(have_allocated_variables) {
      if(mypseudo) delete mypseudo;
      if(mysys) delete mysys;
      deallocate(mywfdata);
      
    }
    if(guidingwf) delete guidingwf;
    if(dyngen) delete dyngen;
    for(int i=0; i< densplt.GetDim(0); i++) {
      if(densplt(i)) delete densplt(i);
    }
  }


  int showinfo(ostream & os);
 private:

  Properties_manager myprop;
  Properties_gather mygather;



  int allocateIntermediateVariables(System * , Wavefunction_data *);
  void deallocateIntermediateVariables() {
    if(sample) delete sample;
    sample=NULL;
    deallocate(wf);
    wf=NULL;
    

  }

  void savecheckpoint(string & filename, Sample_point *);
  void restorecheckpoint(string & filename, System * sys,
			 Wavefunction_data * wfdata,Pseudopotential * pseudo);
  void cdmcReWeight(Array2 <doublevar> & energy_temp, 
                    Array1 < Wf_return> & value_temp);

  
  doublevar getWeight(Dmc_point & pt,
                      doublevar teff, doublevar etr);

  void getAuxWeight(Dmc_point & pt,
                    doublevar teff, 
                    Array1 <doublevar>&  aux_teff, 
                    Array2 <doublevar> & aux_weight);
  
  int calcBranch();
  void find_cutoffs();
  void loadbalance();
  void updateEtrial(doublevar feedback);
  

  int have_allocated_variables;
  int have_read_options;
  int do_cdmc;
  int low_io; //!< write out configs and densities only at the end.
  int aux_converge; //!< Number of projection lengths at which to evaluate correlated sampling
  int tmoves; //!< whether to do Casula's t-moves

  int nblock, nstep;
  doublevar eref; //!< reference energy-best guess at dmc energy
  doublevar etrial; //!< current trial energy
  doublevar timestep;
  string readconfig, storeconfig;
  string log_label;

  int nconfig;

  int nelectrons;
  int ndim;
  int nwf;
  
  doublevar feedback;
  int feedback_interval;
  doublevar start_feedback;

  int nhist; //!< amount of history to go back when doing correlated sampling

  Dynamics_generator * dyngen;

  Dmc_guiding_function * guidingwf;
  Array1 <doublevar> sum_proportions;
  string guidetype;
  
  
  //If a walker happens to have a very low energy, we will cut
  //it off smoothly with these points
  doublevar branchcut_start, branchcut_stop;

  //these times the standard deviation are branchcut_start, branchcut_stop
  doublevar branch_start_cutoff, branch_stop_cutoff;

  System * mysys;
  Pseudopotential * mypseudo;
  Sample_point * sample;
  Wavefunction * wf;
  Wavefunction_data * mywfdata;

  Array1 <Dmc_point> pts;

  Array1 < Local_density_accumulator *> densplt;
  vector <vector <string> > dens_words;
  Array1 < Nonlocal_density_accumulator *> nldensplt;
  vector <vector <string> > nldens_words;

};




#endif //DMC_METHOD_H_INCLUDED
//------------------------------------------------------------------------
