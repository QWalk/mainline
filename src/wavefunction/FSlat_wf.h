/*
 
Based on Slat_wf.h Copyright (C) 2007 Lucas K. Wagner

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

#ifndef FSLAT_WF_H_INCLUDED

#define FSLAT_WF_H_INCLUDED

#include "Qmc_std.h"
#include "Array.h"
#include "Wavefunction.h"
class Wavefunction_data;
class Slat_wf_data;
class System;


class FSlat_wf_storage : public Wavefunction_storage
{
public:
  virtual ~FSlat_wf_storage()
  {}
private:
  friend class FSlat_wf;

  //dimensions are [value gradient lap, MO]
  Array2 <doublevar>  moVal_temp;

  // Added by Matous
  Array2 <doublevar>  moVal_temp_2;
 
  Array3 < Array2 <doublevar> > inverse_temp;
  Array3 <doublevar> detVal_temp;

};


/*!
A "fast" slater wavefunction; \f$\Psi=\sum_i det_i(\Phi_1\Phi_2...)\f$
where the \f$\Phi\f$'s are one-particle molecular orbitals.
Also supports multiple states as a vector of wavefunction values.  Just
specify multiple STATE keywords. "Fast" implementation avoids storage
of all inverse determinants via iterative update procedure.
*/
class FSlat_wf : public  Wavefunction
{

public:

  FSlat_wf()
  {}


  virtual int nfunc() {
    return nfunc_;
  }


  virtual void notify(change_type , int );

  virtual void updateVal(Wavefunction_data *, Sample_point *);
  virtual void updateLap(Wavefunction_data *, Sample_point *);

  virtual void getVal(Wavefunction_data *, int, Wf_return &);
  virtual void getLap(Wavefunction_data *, int, Wf_return &);

  virtual void getDensity(Wavefunction_data *,int, Array2 <doublevar> &);

  virtual void saveUpdate(Sample_point *, int e, Wavefunction_storage *);
  virtual void restoreUpdate(Sample_point *, int e, Wavefunction_storage *);

  // Added by Matous
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


  void generateStorage(Wavefunction_storage * & wfstore);


  void init(Wavefunction_data *);

  //--
private:

  void save_for_static();

  void updateInverse(Slat_wf_data *, int e);
  int updateValNoInverse(Slat_wf_data *, int e); 
  //!< update the value, but not the inverse.  Returns 0 if the determinant is zero and updates aren't possible
  
  void calcVal(Slat_wf_data *, Sample_point *);
  void updateVal(Slat_wf_data *, Sample_point *, int);
  void calcLap(Slat_wf_data *, Sample_point *);
  void updateLap(Slat_wf_data *, Sample_point *, int);

  Array1 <int> electronIsStaleVal;
  Array1 <int> electronIsStaleLap;
  int updateEverythingVal;
  int updateEverythingLap;

  int sampleAttached;
  int dataAttached;
  Slat_wf_data * parent;

  //lazy updates of the determinant(which saves a lot of time in pseudopotentials, etc)
  int inverseStale;
  int lastValUpdate;
  Array3<doublevar> lastDetVal;
  
  //Saved variables for electron updates
  Array3 <doublevar>  moVal;

  Array2 <doublevar> updatedMoVal;

  Array3 < Array2 <doublevar> > inverse;
  //!<inverse of the value part of the mo_values array transposed

  Array3 <doublevar> detVal;

  //Variables for a static(electrons not moving) calculation
  int staticSample;
  Array3 <doublevar> saved_laplacian;
  //!<Saved laplacian for a static calculation (electron, function, [val grad lap])

  Array1<doublevar> invgammastore;
  Array1< Array1<doublevar> > ustore, veestore;

  int nmo;        //!<Number of molecular orbitals
  int ndet;       //!<Number of determinants
  int ndim;       //!<Number of (spacial) dimensions each electron has
  int nfunc_;      //!<Number of functions this class represents.
  Array1 <int> nelectrons; //!< 2 elements, for spin up and down
  Array1 <int> spin;       //!< lookup table for the spin of a given electron

};

#endif //FSLAT_WF_H_INCLUDED
//--------------------------------------------------------------------------
