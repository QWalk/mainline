/*
 
Copyright (C) 2010 Lucas K. Wagner

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


#ifndef RMC_CORR_METHOD_H_INCLUDED
#define RMC_CORR_METHOD_H_INCLUDED

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

//All the important variables in a given slice
struct Rmc_corr_history { 
  Array1 <doublevar> main_en;
  Array2 <Wf_return> wfs; //system number, electron number
  Array1 <Config_save_point> configs;
  Array1 <doublevar> jacobian;
  Properties_point prop;

  void mpiSend(int node);
  void mpiReceive(int node);
  void read(istream & is);
  void write(ostream & os);
};


//This is a point in the dNs space, where d is dimension,
//N is the number of electrons, and s is the number of slices.
struct Rmc_corr_point { 
  //Properties_point prop;
  deque <Rmc_corr_history> path;
  int system; // the system that sampled this point.
  Array1 <doublevar> weight;
  Array1 <doublevar> branching;
  Array1 <doublevar> totgf;
  int recalc_gf;
  int direction;
  Rmc_corr_point() { 
    //weight=1;
    recalc_gf=1;
    direction=1;
  }
  
  //Careful when using these; I've implemented them to send
  //only the most basic parts, so a lot of stuff needs
  //to be recalculated.  Works fine for reading from the config
  //files.
  void mpiSend(int node);
  void mpiReceive(int node);
  void read(istream & is);
  void write(ostream & os);
  
};


class Rmc_corr_method : public Qmc_method
{
public:


  void read(vector <string> words,
            unsigned int & pos,
            Program_options & options);
  int generateVariables(Program_options & options);
  void run(Program_options & options, ostream & output);

  Rmc_corr_method() {
    have_read_options=0;
    have_allocated_variables=0;
  }
  ~Rmc_corr_method()
  {
    for(int i=0; i< mywfdata.GetDim(0); i++)  if(mywfdata(i)) delete mywfdata(i);
    for(int i=0; i< mysys.GetDim(0); i++)  if(mysys(i)) delete mysys(i);
    if(guidingwf) delete guidingwf;
    if(dyngen) delete dyngen;
    if(mypseudo) delete mypseudo;

  }


  int showinfo(ostream & os);
 private:


  void readcheck();
  void storecheck();
  

  void deallocateIntermediateVariables() {
    for(int i=0; i< wf.GetDim(0); i++)  if(wf(i)) delete wf(i);
    for(int i=0; i< sample.GetDim(0); i++) if(sample(i)) delete sample(i);
  }

  void savecheckpoint(string & filename, Sample_point *);
  void restorecheckpoint(string & filename, System * sys,
			 Wavefunction_data * wfdata,Pseudopotential * pseudo);

  doublevar propagate_walker(int walker);
  Rmc_corr_history add_point(int walker); //calculate energies and other values for a new point
  doublevar get_green_weight(Rmc_corr_history & a,
                             Rmc_corr_history & b,
                             int sys, doublevar & branching);
    
  
  
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
  int nhist; //!< amount of history to go back 
  int pc_gf; //!< whether or not to use the Pierleoni-Ceperley correlated sampling
  int dmc_gf; //!< use the green's function suggested by Umrigar (suggested not to use, actually!)
  
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

  Array1 <Rmc_corr_point> pts;
  Array2 <Local_density_accumulator *> local_dens;
  vector <vector <string> > denswords;
  Array2 <Average_generator *> avggen;
  vector <vector < string> > avgwords;
  
  Space_warper warper;
  Properties_manager myprop;
  Properties_manager centerprop;
  Properties_gather mygather;
  
};




#endif //RMC_CORR_METHOD_H_INCLUDED
//------------------------------------------------------------------------
