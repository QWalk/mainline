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

#ifndef CSLAT_WF_H_INCLUDED

#define CSLAT_WF_H_INCLUDED

#include "Qmc_std.h"
#include "Array.h"
#include "Wavefunction.h"
class Wavefunction_data;
class Slat_wf_data;
class System;


class Cslat_wf_storage : public Wavefunction_storage
{
public:
  virtual ~Cslat_wf_storage()
  {}
private:
  friend class Cslat_wf;

  //dimensions are [value gradient lap, MO]
  Array2 <dcomplex>  moVal_temp;
  Array3 < Array2 <dcomplex> > inverse_temp;
  Array3 <dcomplex> detVal_temp;

};


/*!
A slater wavefunction; \f$\Psi=\sum_i det_i(\Phi_1\Phi_2...)\f$
where the \f$\Phi\f$'s are one-particle molecular orbitals.
Also supports multiple states as a vector of wavefunction values.  Just
specify multiple STATE keywords.
*/
class Cslat_wf : public  Wavefunction
{

public:

  Cslat_wf()
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


  virtual void storeParmIndVal(Wavefunction_data *, Sample_point *,
                               int, Array1 <doublevar> & );
  virtual void getParmDepVal(Wavefunction_data *,
                             Sample_point *,
                             int,
                             Array1 <doublevar> &,
                             Wf_return &);


  void generateStorage(Wavefunction_storage * & wfstore);


  void init(Wavefunction_data *);

  //--
private:

  void save_for_static();

  void calcVal(Sample_point *);
  void updateVal(Sample_point *, int);
  void calcLap(Sample_point *);
  void updateLap(Sample_point *, int);

  Array1 <int> electronIsStaleVal;
  Array1 <int> electronIsStaleLap;
  int updateEverythingVal;
  int updateEverythingLap;

  int sampleAttached;
  int dataAttached;
  Slat_wf_data * parent;

  //Saved variables for electron updates
  Array3 <dcomplex>  moVal;

  Array2 <dcomplex> updatedMoVal;

  Array3 < Array2 <dcomplex> > inverse;
  //!<inverse of the value part of the mo_values array transposed

  Array3 <dcomplex> detVal;

  //Variables for a static(electrons not moving) calculation
  int staticSample;
  Array3 <dcomplex> saved_laplacian;
  //!<Saved laplacian for a static calculation (electron, function, [val grad lap])


  int nmo;        //!<Number of molecular orbitals
  int ndet;       //!<Number of determinants
  int ndim;       //!<Number of (spacial) dimensions each electron has
  int nfunc_;      //!<Number of functions this class represents.
  Array1 <int> nelectrons; //!< 2 elements, for spin up and down
  Array1 <int> spin;       //!< lookup table for the spin of a given electron

};

#endif //CSLAT_WF_H_INCLUDED
//--------------------------------------------------------------------------
