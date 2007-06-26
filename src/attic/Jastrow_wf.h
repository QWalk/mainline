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
//--------------------------------------------------------------------------
// include/Jastrow_wf.h
//
//
#ifndef JASTROW_WF_H_INCLUDED

#define JASTROW_WF_H_INCLUDED

#include "Qmc_std.h"
#include "Array.h"
#include "Wavefunction.h"
#include "Jastrow_wf_data.h"
#include "Array45.h"

class System;

class Jastrow_wf_storage : public Wavefunction_storage
{
private:
  friend class Jastrow_wf;

  //!The natural logarithm of the value(U)
  doublevar value_temp;

  doublevar valPartialSum_temp;

  Array1 <doublevar> elecPartialSum_temp;
  doublevar ionPartialSum_temp;
  Array1 < Array1 <doublevar> > a_kval_temp;
};


/*!
 
Jastrow factor = \f$e^U\f$, where
\f[ U=\sum_{i<j, I}u(r_{iI}, r_{jI}, r_{ij}) \f]
and
 
\f[
u(r_{iI}, r_{jI}, r_{ij}) = \frac{-c}{\gamma} e^{-\gamma r_{ij}}
+ \sum_{k<l,m} c_{klm}v_{kl}(r_{iI}, r_{jI})b_m(r_{ij})
\f]
 
\f[
v_{kl}(r_{iI}, r_{jI})=a_k(r_{iI})a_l(r_{jI})+a_k(r_{jI})a_l(r_{iI})
\f]
 
\f[
\nabla_i u_{ijI} =
\frac{\vec{r_{ij}}}{r_{ij}}c e^{-\gamma r_{ij}}
+ \sum c_{klm} v_{kl} (\nabla_i b_m) + (\nabla_i v_{kl} )b_m
\f]
 
\f[
\nabla_i v_{kl} = (\nabla_i a_k(r_{iI}))a_l(r_{jI})+a_k(r_{jI})(\nabla_i a_l(r_{iI}))
\f]
 
\f[
\nabla_j u_{ijI} =
\frac{\vec{-r_{ij}}}{r_{ij}}c e^{-\gamma r_{ij}}
+ \sum c_{klm} \left[ v_{kl} (-\nabla_i b_m) + (\nabla_i v_{kl} )b_m \right]
\f]
 
\f[
\nabla_i^2 u_{ijI} =
\left( \frac{2}{r_{ij}} - \gamma \right) c e^{-\gamma r_{ij}}
+ \sum c_{klm} ( v_{kl}(\nabla_i^2b_m)
+ 2 (\nabla_i v_{kl}) \cdot (\nabla_i b_m)
+ (\nabla_i^2 v_{kl}) b_m )
\f]
 
So,
\f[
\nabla_n U = \sum_{I} \left[ \sum_{i=0}^{n-1} \nabla_n u_{inI}
+ \sum_{i=n+1}^{jmax} \nabla_n u_{njI} \right]
\f]
 
and
\f[
\nabla_n^2 U
= \sum_{I} \left[ \sum_{j=n+1}^{jmax} \nabla_n^2 u_{njI}
+ \sum_{j=0}^{n-1} \nabla_n^2 u_{jnI} \right]
\f]
 
 
Optimizable parameters are the \f$c_{klm}\f$'s and  the
\f$ \gamma \f$'s.

 
*/
class Jastrow_wf : public  Wavefunction
{

public:

  virtual void notify(change_type , int );

  virtual void updateVal(Wavefunction_data *, Sample_point *);
  virtual void updateLap(Wavefunction_data *, Sample_point *);
  virtual void updateForceBias(Wavefunction_data *, Sample_point *);

  virtual void getVal(Wavefunction_data *, int, Array2 <doublevar> &,
                      int startpos=0);
  virtual void getLap(Wavefunction_data *, int, Array2 <doublevar> &,
                      int startpos=0 );
  virtual void getForceBias(Wavefunction_data *, int, Array2 <doublevar> &,
                            int startpos=0);

  virtual void getDensity(Wavefunction_data *,int,  Array2 <doublevar> &);


  virtual void storeParmIndVal(Wavefunction_data *, Sample_point *,
                               int, Array1 <doublevar> & );
  virtual void getParmDepVal(Wavefunction_data *,  Sample_point *,
                             int,
                             Array1 <doublevar> &,
                             Array2 <doublevar> &);

  virtual int getParmDeriv(Wavefunction_data *, 
			   Sample_point *,
			   Parm_deriv_return & );


  virtual void saveUpdate(Sample_point *, int e, Wavefunction_storage *);
  virtual void restoreUpdate(Sample_point *, int e, Wavefunction_storage *);

  void generateStorage(Wavefunction_storage * & wfstore);
  void init(Wavefunction_data *);

private:
  int nelectrons, nions;

  doublevar value;
  //!The natural logarithm of the value(U)

  Array2 <doublevar> derivatives;
  //!the derivatives of the Jastrow factor, with indices(electron#, [null, dx, dy, dz, lap])

  Array2 <doublevar> elecPartialSum;
  Array1 <doublevar> ionPartialSum;

  Array1 <doublevar> valPartialSum;
  //!The contribution from each electron to the value.  Note that sum(partialSum) != value, since there are two-body interactions.

  Array1 <int> electronIsStaleVal;
  int updateEverythingLap;
  int updateEverythingVal;
  int sampleAttached;
  int dataAttached;

  int staticSample;
  //!< Whether we have a unmoving sample.

  int updateStatic;
  //!<If we need to update the saved variables used in static calcs

  int parmChanged;//!<Whether a parameter has changed

  Array2 < Array1 <doublevar> > a_kval;

  //saved variables for a static calculation.
  Array3 <doublevar> deriv_static;
  Array2 <doublevar> partialValues_static;
  Array1 <doublevar> values_static;

  void calcVal(Jastrow_wf_data * , Sample_point *);
  void updateVal(Jastrow_wf_data *, Sample_point *,int);
  void calcLap(Jastrow_wf_data *, Sample_point *);
};

#endif //JASTROW_WF_H_INCLUDED
//--------------------------------------------------------------------------
