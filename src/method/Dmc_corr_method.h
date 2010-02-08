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


#ifndef DMC_CORR_METHOD_H_INCLUDED
#define DMC_CORR_METHOD_H_INCLUDED

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

struct Dmc_corr_history { 
  Array1 <doublevar> main_en;
  void mpiSend(int node);
  
  void mpiReceive(int node);

};



struct Dmc_corr_point { 
  Properties_point prop;
  deque <Dmc_corr_history> past_energies;
  doublevar weight;
  int system; // the system that sampled this point.
  Config_save_point config_pos;  
  Array1 <doublevar> age;  //!< age of each electron
  Array1 <doublevar> jacobian;
  Dmc_corr_point() { 
    weight=1;
  }
  void mpiSend(int node);
  void mpiReceive(int node);
  
};


class Dmc_corr_method : public Qmc_method
{
public:


  void read(vector <string> words,
            unsigned int & pos,
            Program_options & options);
  int generateVariables(Program_options & options);
  void run(Program_options & options, ostream & output);

  Dmc_corr_method() {
    have_read_options=0;
    have_allocated_variables=0;
  }
  ~Dmc_corr_method()
  {
    /*
    if(have_allocated_variables) {
      if(mypseudo) delete mypseudo;
      if(mysys) delete mysys;
      deallocate(mywfdata);
      
    }
    if(guidingwf) delete guidingwf;
    if(dyngen) delete dyngen;
    for(int i=0; i< densplt.GetDim(0); i++) {
      if(densplt(i)) delete densplt(i);
    }*/
    for(int i=0; i< mywfdata.GetDim(0); i++)  if(mywfdata(i)) delete mywfdata(i);
    for(int i=0; i< mysys.GetDim(0); i++)  if(mysys(i)) delete mysys(i);
    if(guidingwf) delete guidingwf;
    if(dyngen) delete dyngen;
    if(mypseudo) delete mypseudo;

  }


  int showinfo(ostream & os);
 private:




  void deallocateIntermediateVariables() {
    for(int i=0; i< wf.GetDim(0); i++)  if(wf(i)) delete wf(i);
    for(int i=0; i< sample.GetDim(0); i++) if(sample(i)) delete sample(i);
  }

  void savecheckpoint(string & filename, Sample_point *);
  void restorecheckpoint(string & filename, System * sys,
			 Wavefunction_data * wfdata,Pseudopotential * pseudo);

  doublevar propagate_walker(int walker);
  void add_point(int walker); //calculate energies and other values for a new point

  
  
  int calcBranch();
  void find_cutoffs();
  void updateEtrial();
  int have_allocated_variables;
  int have_read_options;
  

  int nblock, nstep;
  doublevar timestep;
  string readconfig, storeconfig;
  string log_label;
  int npsteps;
  int nconfig;
  int nhist; //!< amount of history to go back when doing correlated sampling

  Dynamics_generator * dyngen;

  Dmc_guiding_function * guidingwf;
  string guidetype;
  //If a walker happens to have a very low energy, we will cut
  //it off smoothly with these points
  //doublevar branchcut_start, branchcut_stop;

  //these times the standard deviation are branchcut_start, branchcut_stop
  //doublevar branch_start_cutoff, branch_stop_cutoff;

  Array1 <System *> mysys;
  Pseudopotential * mypseudo;
  Array1 <Sample_point *> sample;
  Array1 <Wavefunction *> wf;
  Array1 <Wavefunction_data *> mywfdata;
  Array1 <doublevar> eref;
  Array1 <doublevar> etrial;

  Array1 <Dmc_corr_point> pts;

  Space_warper warper;
  Properties_manager myprop;
  Properties_gather mygather;

};




#endif //DMC_CORR_METHOD_H_INCLUDED
//------------------------------------------------------------------------
