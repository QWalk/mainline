/*
 
Copyright (C) 2007 Lucas K. Wagner, 2008 Jindrich Kolorenc

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

/*
JK: Present form of this BCS wave function was (re)created as
a simplification and subsequent improvements of an abandoned
piece of code written originally by Lucas.
*/

#ifndef BCS_WF_H_INCLUDED

#define BCS_WF_H_INCLUDED

#include "Qmc_std.h"
#include "Array.h"
#include "Wavefunction.h"
#include "Jastrow2_wf.h"
class Wavefunction_data;
class BCS_wf_data;
class System;
class BCS_jastrow_cofactor;

class BCS_wf_storage : public Wavefunction_storage
{
public:
  BCS_wf_storage() { jast_store=NULL; } 
  virtual ~BCS_wf_storage()
  {if(jast_store) delete jast_store; }
  
private:
  friend class BCS_wf;
  Wavefunction_storage * jast_store;
  Array1 <Array2 <doublevar> > inverse;
  //Array2 <doublevar> moVal;
  Array2 <doublevar> derivatives;
  //Added by Matous
  Array2 <doublevar> derivatives_2;
  Array1 <doublevar> detVal;

};


/*!
\brief
BCS wave function. At present it supports only non-polarized systems, i.e.
n_up=n_down, in which case the wave function form is a determinant constructed
out of a two-particle orbital, \f$\Psi=det_(\phi_{ij})\f$ with
\f$\phi_{ij}=\phi(r_i,r_j)\f$.

The two-particle orbital is constructed in the same way as is constructed
the two-body part of the Jastrow factor. (?BUG: It is possible to input
also other parts of the Jastrow factor, but those will be silently ignored.)

In this form, the BCS wave function is only good for homogeneous systems
(HEG). See Pfaffian for more general situations.

Partial support for multi-determinantal form is implemented. TODO: input
and storage of multiple pair orbitals.
*/

class BCS_wf : public  Wavefunction
{

public:

  BCS_wf()
  {jastcof=NULL;}

  ~BCS_wf(); 

  virtual int nfunc() {
    return 1;
  }


  virtual void notify(change_type , int );

  virtual void updateVal(Wavefunction_data *, Sample_point *);
  virtual void updateLap(Wavefunction_data *, Sample_point *);

  virtual void getVal(Wavefunction_data *, int, Wf_return &);
  virtual void getLap(Wavefunction_data *, int, Wf_return &);

  virtual void getDensity(Wavefunction_data *,int, Array2 <doublevar> &);

  virtual void saveUpdate(Sample_point *, int e, Wavefunction_storage *);
  virtual void restoreUpdate(Sample_point *, int e, Wavefunction_storage *);

  //Added by Matous
  virtual void saveUpdate(Sample_point *, int e1, int e2, Wavefunction_storage *);
  virtual void restoreUpdate(Sample_point *, int e1, int e2, Wavefunction_storage *);

  virtual void storeParmIndVal(Wavefunction_data *, Sample_point *,
                               int, Array1 <doublevar> & );
  virtual void getParmDepVal(Wavefunction_data *,
                             Sample_point *,
                             int,
                             Array1 <doublevar> &,
                             Wf_return &);

  virtual int getParmDeriv(Wavefunction_data *, 
			   Sample_point *,
			   Parm_deriv_return & );

  virtual void getSymmetricVal(Wavefunction_data *, 
			       int, 
			       Wf_return &);

  virtual void plot1DInternals(Array1 <doublevar> &,
			       vector <Array1 <doublevar> > &,
			       vector <string> &,
			       string );

  void generateStorage(Wavefunction_storage * & wfstore);


  void init(Wavefunction_data *);

  //--
private:


  void calcLap(Sample_point *);
  void updateLap(int e, Sample_point *);
  Array1 <int> electronIsStaleVal;
  Array1 <int> electronIsStaleLap;
  int updateEverythingVal;
  int updateEverythingLap;

  int sampleAttached;
  int dataAttached;
  BCS_wf_data * parent;

  //Saved variables for electron updates
  //Array3 <doublevar>  moVal;//(electron,mo,[val grad lap])
  //Array2 <doublevar> updatedMoVal;//(mo,[val grad lap])
  Array1 < Array2 <doublevar> > inverse;
  //!<inverse of the value part of the mo_values array transposed
  //(det)(elec_up,elec_down)

  Array1 <doublevar> detVal;

  Array3 <doublevar> twobody; 
  //!< twobody correlation from the jastrow factor

  Array3 <doublevar> derivatives;
  //(tot#electrons, ndown, gradlap)
  //Keeps the derivatives of the pairing wave functions:
  //D_ij = grad_i phi(r_i,r_j) for i < nup and
  //D_ij = grad_j phi(r_i,r_j) for i >= nup

  //int staticSample;

  //int nmo;      //!<Number of molecular orbitals
  int ndet;       //!<Number of determinants
  int ndim;       //!<Number of (spacial) dimensions each electron has
  Array1 <int> nelectrons; //!< 2 elements, for spin up and down
  Array1 <int> spin;       //!< lookup table for the spin of a given electron

  Jastrow2_wf jast;
  BCS_jastrow_cofactor * jastcof;
  
};

#endif //BCS_WF_H_INCLUDED
//--------------------------------------------------------------------------
