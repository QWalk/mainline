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

#ifndef PSEUDOPOTENTIAL_H_INCLUDED
#define PSEUDOPOTENTIAL_H_INCLUDED

#include "Array.h"
#include "Qmc_std.h"
#include "Wavefunction.h"
#include "Force_fitter.h"
class Program_options;
class System;
class Wavefunction;
class Sample_point;
class Wavefunction_data;
class Wavefunction_storage;
class Basis_function;


//---------
class Pseudo_buffer { 
public:
  void push_value(double d) { 
    data.push_back(d);
  }
  double next_value() {
    assert(current_point < data.size());
    return data[current_point++];
  }
  void clear() { data.clear(); current_point=0; } 
  void start_again() { current_point=0; } 
private:
  vector <double> data;
  int current_point;
};

//-----------

struct Tmove { 
public:
  Array1 <doublevar> pos;
  int e;
  doublevar vxx;
};

//-----------

/*!

*/
class Pseudopotential
{
public:
  Pseudopotential():maxaip(85)
  {deterministic=0;}
  


  void setDeterministic(int i) {
    assert(i==0 || i==1);
    deterministic=i;
  }
  void read(vector < vector <string> > & pseudotext,
            System * sys);
  int showinfo(ostream & os);

  /*!
    \brief
    Standard calculation of the pseudopotential.

    Uses random evaluation of the pseudopotential automatically.
   */
  void calcNonloc(Wavefunction_data * wfdata, System *,
                  Sample_point * sample, Wavefunction * wf,
                  Array1 <doublevar> & totalv);


  void calcNonlocTmove(Wavefunction_data * wfdata, System *,
                       Sample_point * sample,
                       Wavefunction * wf,
                       Array1 <doublevar> & totalv,  //total p.e. from the psp
                       vector <Tmove> & tmoves  //variables for T-moves of Casula
                       );
  /*!
    \brief
    Provide your own random numbers for random evaluation of psp

    If all of them are set to one, everything is evaluated.  If it's set
    to zero, nothing is evaluated.
  */
  void calcNonlocWithTest(Wavefunction_data *, System *,Sample_point *, Wavefunction *,
                          const Array1 <doublevar> & accept_var,
                          Array1 <doublevar> & totalv);

  /*!
    \brief
    Uses file created by initializeStatic to calculated the nonlocal energy

    Uses a cutoff radius for evaluation, completely nonrandom
   */
  void calcNonlocWithFile(Wavefunction_data *, System *,Sample_point *, Wavefunction *,
                          Array1 <doublevar> &, Pseudo_buffer & input);


  
  void calcNonlocParmDeriv(Wavefunction_data * wfdata, System *,
                                            Sample_point * sample,
                                            Wavefunction * wf,
                                            const Array1 <doublevar> & accept_var,
                                            Array1 <doublevar> & totalv, Array1 <doublevar> & parm_deriv);
    
  /*!
    The worker function; the rest just provide simple defaults when functions don't need everything
   
    */
  void calcNonlocWithAllvariables(Wavefunction_data * wfdata, System *,
                                  Sample_point * sample,
                                  Wavefunction * wf,
                                  const Array1 <doublevar> & accept_var, //random variables for stochastic evaluation
                                  Array1 <doublevar> & totalv,  //total p.e. from the psp
                                  bool do_tmoves,vector <Tmove> & tmoves,  //variables for T-moves of Casula
                                  bool parm_derivatives, Array1 <doublevar> & parm_deriv //derivatives wrt wf parameters
                                  );
  /*!
    \brief
    Initialize a file for use with calcNonlocWithFile
   */
  int initializeStatic(Wavefunction_data *, Sample_point *,
                       Wavefunction *, Pseudo_buffer & output);



  
  /*!
    \brief
    number of values between zero and one for random evaluation
    of the pseudopotential
   */
  int nTest();
  
  /*!
    \brief
    Randomly rotate the axes to unbias the calculation
  */
  void randomize();
  
  /*!
    \brief
    Rotate the quadrature to the axes specified in x,y,z
  */
  void rotateQuadrature(Array1 <doublevar> & x,
                        Array1 <doublevar> & y,
                        Array1 <doublevar> & z);

  ~Pseudopotential();


private:
  int deterministic;

  void getRadial(int at, int spin, Sample_point * sample,
                                Array1 <doublevar> & r, Array1 <doublevar> & v_l);
  void getRadial(int at, int spin, Sample_point * sample,
                                  Array1 <doublevar> & r, 
                                  Array2 <doublevar> & v_l);
  const int maxaip; //!< Maximum number of atomic integration points
  int nelectrons;
  Array1 <int> aip;
  Array3 <doublevar> integralpt;
  Array3 <doublevar> integralpt_orig;
  Array2 <doublevar> integralweight;
  Array1 <doublevar> cutoff;
  Array1 <int> numL;
  vector <string> atomnames;
  Array1 <bool> addzeff; //!< whether or not to add Z_eff/r to the local function

  Storage_container wfStore;

  Array2 <Basis_function *> radial_basis;

};

/*!
Computes the integration points for a quadrature for spherical harmonics
*/
void gesqua(int & nq, Array1 <doublevar> & xq,Array1 <doublevar> & yq,
            Array1 <doublevar> & zq, Array1 <doublevar> & wq);


#endif //PSEUDOPOTENTIAL_H_INCLUDED
//------------------------------------------------------------------------
