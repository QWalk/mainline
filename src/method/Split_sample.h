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


#ifndef SPLIT_SAMPLE_H_INCLUDED
#define SPLIT_SAMPLE_H_INCLUDED

#include "Qmc_std.h"
#include "Wavefunction.h"
#include "System.h"
#include "Sample_point.h"
#include "Guiding_function.h"
struct Point {
  Array1 <doublevar> drift; //!< total drift(deterministic move)
  Array1 <doublevar> pos;
  doublevar sign; //!< sign from the sample_point
  Wf_return lap;
  Array1 <doublevar> gauss; //!< gaussian part that generated this move
  Array1 <doublevar> translation; //!< total translation of the move

  Point() {
    drift.Resize(3);
    pos.Resize(3);
    gauss.Resize(3);
    translation.Resize(3);
    translation=0;
    gauss=0;
    pos=0; drift=0;
  }
};

enum drift_type {drift_cutoff, drift_cyrus};

struct Dynamics_info {
  doublevar green_forward;
  doublevar green_backward;
  Array1 <doublevar> orig_pos;
  Array1 <doublevar> gauss; //!< gaussian random numbers used for the move
  Array1 <doublevar> new_pos;
  Array1 <doublevar> diffuse_start; //!< position of start of diffusion
  Array1 <doublevar> diffuse_end; //!< position of end of diffusion
  doublevar acceptance;
  int accepted;
  //Array2 <doublevar> newlap;
  doublevar diffusion_rate;
  
  doublevar symm_gf; //!< ratio versus symmetrized green's function
  doublevar resample_gf; //!< green's function that we want to resample
  Dynamics_info() { 
     symm_gf=1;
  }
};


class Dynamics_generator {
 public:
  virtual void read(vector <string> & words)=0;
  virtual int sample(int e,
                     Sample_point * sample, 
                     Wavefunction * wf, 
                     Wavefunction_data * wfdata, 
                     Guiding_function * guidewf,
                     Dynamics_info & info,
                     doublevar & efftimestep
                     )=0;


  //returns the acceptance ratio                                  
  virtual doublevar greenFunction(Sample_point * sample, Wavefunction * wf,
                     Wavefunction_data * wfdata, Guiding_function * guidewf,
                             int e,
                             Array1 <doublevar> & newpos, doublevar timestep,
                             Dynamics_info & info, Dynamics_info & oldinfo);

  virtual void enforceNodes(int i) {
    assert(i==0 || i==1);
    restrict_nodes=i;
  };

  virtual int showinfo(string & indent, ostream & os)=0;
  virtual void showStats(ostream & os)=0;
  virtual  void resetStats()=0;


  virtual ~Dynamics_generator() {}
 protected:
  int restrict_nodes;
};


/*!
\todo
Make this average the acceptance ratio over all processors

 */
class Split_sampler:public Dynamics_generator {
 public:

  Split_sampler() {
    divide_=1.0;
    recursion_depth_=1;
    restrict_nodes=0;
    dtype=drift_cyrus;
    //drift_2pt=0;

    acceptances.Resize(recursion_depth_);
    tries.Resize(recursion_depth_);
    acceptances=0;
    tries=0;
  }

  void read(vector <string> & words);

  void setDriftType(drift_type dtype_) {
    dtype=dtype_;
  }
  
  int showinfo(string & indent, ostream & os) {
    os << indent << "recursion depth " << recursion_depth_ << endl;
    os << indent << "timestep divider " << divide_ << endl;
    if(restrict_nodes) os << indent << "restricting node crossings" << endl;
    os << indent << "drift type " << dtype << endl;
    return 1;
  }

  void setDivider(doublevar divide) {
    divide_=divide;
  }

  void setRecursionDepth(int depth) {
    recursion_depth_=depth;
    acceptances.Resize(recursion_depth_);
    tries.Resize(recursion_depth_);
    acceptances=0;
    tries=0;
  }

  int getRecursionDepth() {
    return recursion_depth_;
  }

  /*!
    Returns the step that was accepted, or 0 if 
    there was a rejection.
  */
  int sample(int e,
             Sample_point * sample, 
             Wavefunction * wf, 
             Wavefunction_data * wfdata, 
             Guiding_function * guidewf,
             Dynamics_info & info,
             doublevar & efftimestep
             );
  int split_driver(int e,
                   Sample_point * sample,
                   Wavefunction * wf, 
                   Wavefunction_data * wfdata,
                   Guiding_function * guidewf,
                   int depth,
                   Dynamics_info & info,
                   doublevar & efftimestep);

  //virtual doublevar greenFunction(Sample_point * sample, Wavefunction * wf,
  //                   Wavefunction_data * wfdata, Guiding_function * guidewf,
  //                           int e,
  //                           Array1 <doublevar> & newpos, doublevar timestep,
  //                           Dynamics_info & info);

  doublevar get_acceptance(Guiding_function * guidingwf, int x, int y);
  
  void showStats(ostream & os);
  void resetStats();
 private:
 
  doublevar transition_prob(int point1, int point2,
                            doublevar timestep, 
                            drift_type dtype);
  
  drift_type dtype;
  Array1 <Point> trace;
  Array1 <doublevar> timesteps;
  int recursion_depth_;
  doublevar divide_;
  
  Array1 <doublevar> acceptances;
  Array1 <long int > tries;

  string indent; //for debugging..
  Storage_container wfStore;
};


doublevar transition_prob(Point &  point1, Point & point2,
                          doublevar timestep, 
                          drift_type dtype);


/*!
sampler from Umrigar, Nightingale, and Runge. 
JCP 99 2865 (1993)
 */
class UNR_sampler:public Dynamics_generator {
 public:

  UNR_sampler() {
    restrict_nodes=0;
    acceptance=0;
    tries=0;
  }

  void read(vector <string> & words) {}

  void setDriftType(drift_type dtype_) {
  }
  
  int showinfo(string & indent, ostream & os) {
    os << "UNR sampler " << endl;
    if(restrict_nodes) os << indent << "restricting node crossings" << endl;
    return 1;
  }

  void getDriftEtc(Point & pt, Sample_point * sample,
		   doublevar tstep, int e,
		   doublevar & exponent, //!< exponent at which to do an exponential move
		   doublevar & probability,//!< probability to do an exponential move
		   Array1 <doublevar> & rnuc //!<position of closest ion
		   );
  int sample(int e,
             Sample_point * sample, 
             Wavefunction * wf, 
             Wavefunction_data * wfdata, 
             Guiding_function * guidewf,
             Dynamics_info & info,
             doublevar & efftimestep
             );

  void showStats(ostream & os);
  void resetStats();
 private:
  doublevar acceptance;
  long int tries;
};

//----------------------------------------------------------------------

class SRK_dmc:public Dynamics_generator {
  public:
    SRK_dmc() {
      resetStats();
      resample=0;
      resample_tol=.001;
      dtype=drift_cyrus;
    }
    int sample(int e,
               Sample_point * sample, 
               Wavefunction * wf, 
               Wavefunction_data * wfdata, 
               Guiding_function * guidewf,
               Dynamics_info & info,
               doublevar & efftimestep
              );
    virtual void read(vector <string> & words);
    virtual doublevar greenFunction(Sample_point * sample, Wavefunction * wf,
                             Wavefunction_data * wfdata, Guiding_function * guidewf,
                             int e, Array1 <doublevar> & newpos, doublevar timestep,
                             Dynamics_info & info, Dynamics_info & oldinfo);

    int showinfo(string & indent, ostream & os) {
      os << indent << "Stochastic Runge-Kutta algorithm (PRA 45 600) " << endl;
      return 1;
    }
  void setDriftType(drift_type dtype_) {
    dtype=dtype_;
  }
  
    void showStats(ostream & os) {
      os << "acceptance " << acceptances/tries <<  " average retries " << doublevar(retries)/doublevar(tries) << endl;
      os << "nbottom " << nbottom << endl;
    }
    void resetStats() {
      acceptances=0; tries=0;
      retries=0; nbottom=0;
    }
  private:
    int rk_step(int e, Sample_point * sample, 
        Wavefunction * wf,Wavefunction_data * wfdata, 
        Guiding_function * guidingwf, Dynamics_info & info,
        doublevar & timestep, Array1 <Point> & trace);
    //Array1 <Point> trace;
    doublevar acceptances;
    long int tries;
    long int retries;
    long int nbottom; //!< number of times we bottom out with the retries
    int resample; //!< whether to try to improve the approximation to the symmetric green's function
    doublevar resample_tol; //!< tolerance for the ratio of green_symm/green_forward

    drift_type dtype;
    bool diagnostics;//!< Whether to print out a ton of diagonostic information on the moves
    bool third_move; //!< Whether to do a third predictor-corrector move
    ofstream diagnostics_print;

 
};



int allocate(vector <string> & words, Dynamics_generator *& sam);
void limDrift(Array1 <doublevar> & drift, doublevar tau, drift_type dtype);


#endif //SPLIT_SAMPLE_H_INCLUDED

//----------------------------------------------------------------------


