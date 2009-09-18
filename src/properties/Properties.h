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
#ifndef PROPERTIES_H_INCLUDED
#define PROPERTIES_H_INCLUDED
#include "Qmc_std.h"
#include "System.h"
#include "Wavefunction.h"
#include "Guiding_function.h"
#include "Properties_point.h"
#include "Properties_block.h"
#include "Properties_average.h"
#include "Properties_gather.h"
#include <list>


//--------------------------------------------------


/*!


 */

class Properties_manager {
 public:


  Properties_manager();
  ~Properties_manager();
  void read(vector <string> & words, 
            vector <string> & systxt, 
            vector < string> & wftxt);

  void setSize(int nwf_, int nblocks, int nsteps, int maxwalkers, System *, 
               Wavefunction_data *,int naux_,int n_aux_cvg=1);
  void setFirstStep(int s) {
    start_avg_step=s;
    int nsteps=trace.GetDim(0);

    if(nsteps-start_avg_step-2 < autocorr_depth) 
      autocorr_depth=nsteps-start_avg_step-2;
  }

  void setLog(string & file, string & label) {
    log_file=file;
    log_label=label;
  }

  void getLog(string & file, string & label) {
    file=log_file;
    label=log_label;
  }

  void insertPoint(int step, 
                   int walker, 
                   const Properties_point & pt);

  void endBlock();

  void endBlock_per_step();

  void initializeLog(Array1 <Average_generator*> & avg_gen);

  void printBlockSummary(ostream & os);
  
  void printSummary(ostream & os) { 
    error("update call to printSummary");
  }
  void printSummary(ostream & os, Array1 <Average_generator*> & avg_gen);

  void getFinal(Properties_final_average & final) {
    final=final_avg;
  }

  void getLastBlock(Properties_block & block) {
    assert(current_block > 0);
    block=block_avg(current_block-1);
  }
 

  


 private:
  void autocorrelation(Array2 <doublevar> &,Array2 <doublevar> & ,  int);
  Array2 < Properties_point > trace;
  //control
  int start_avg_step;

  int max_autocorr_depth; 

  string log_file;
  string log_label;


  //Limits
  int nwf; //!< the number of wave functions we can have in the primary walk
  int maxchildren; //!< maximum number of children a walker can have
  int autocorr_depth;
  int num_aux_converge; //!<convergence points in the auxillary walk

  int naux;  //!< number of auxillary energies to collect(forces)
  Array1 <doublevar> aux_size; //!< size of the distortion

  //Averaging 
  int current_block;
  Array1 < Properties_block > block_avg;
  Properties_block local_sum;

  Properties_final_average final_avg;



};



//####################################################################################
//---------------------------------------------------------------------
//The Local_density_accumulators are for *local* averages that may be very large 
//and shouldn't be in the normal averaging substructure.  One can still evaluate error bars
//by printing out the average every block and then estimating from there, although 
//in some cases, even that may be prohibitively expensive.

class Local_density_accumulator { 
 public:
    virtual void init(vector <string> & , System *, string & runid)=0;
    virtual void accumulate(Sample_point *, doublevar weight)=0;
    /*!
      Write out the density.  Note that when running in parallel, this
      must be called by all processes!
    */
    virtual void write()=0;
    virtual ~Local_density_accumulator() { }
};

class Nonlocal_density_accumulator {
 public:
    virtual void init(vector <string> & , System *, string & runid)=0;
    virtual void accumulate(Sample_point *, doublevar weight,
			    Wavefunction_data *, Wavefunction *)=0;
    virtual void write(string & log_label)=0;
    virtual ~Nonlocal_density_accumulator() { }
};


class One_particle_density:public Local_density_accumulator { 
 public:
  void init(vector <string> & , System *, string & runid);
  void accumulate(Sample_point *, doublevar weight);

  /*!
    Write out the density.  Note that when running in parallel, this
    must be called by all processes!
   */
  void write(); 

 protected:
  Array3 <doublevar> bin; //!< count of hits
  Array1 <doublevar> min_; //!< min position
  Array1 <doublevar> max; //!< max position
  Array1 <int> npoints;
  doublevar resolution;
  doublevar nsample; //double because we can have weights
  int nup;
  string outputfile;
  Array2 <doublevar> atominfo;//!< (ion#, [charge, x,y,z])
  doublevar norm; //!< renormalization for the density output

  int start_electron; //!< which electrons to average over(for spin density)
  int end_electron; 
};




void allocate(vector<string> & words, System * sys, string & runid,
	 Local_density_accumulator  *& denspt);
void allocate(vector<string> & words, System * sys, string & runid,
	 Nonlocal_density_accumulator  *& nldenspt);

//---------------------------------------------------------------------


#endif //PROPERTIES_H_INCLUDED
//----------------------------------------------------------------------
