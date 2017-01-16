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
// include/Ext_slat.h
//
//
#ifndef EXT_SLAT_H_INCLUDED

#define EXT_SLAT_H_INCLUDED

#include "Qmc_std.h"
#include "Array.h"
#include "Wavefunction.h"
class Wavefunction_data;
class Ext_slat_data;
class System;


class Ext_slat_storage : public Wavefunction_storage
{
public:
  virtual ~Ext_slat_storage()
  {}
private:
  friend class Ext_slat;

  //dimensions are [value gradient lap, MO]
  Array2 <doublevar>  moVal_temp;
  Array2 <doublevar>  inverse_temp;
  doublevar detVal;

};


/*!
A slater wavefunction; \f$\Psi=\sum_i det_i(\Phi_1\Phi_2...)\f$
where the \f$\Phi\f$'s are one-particle molecular orbitals.
Also supports multiple states as a vector of wavefunction values.  Just
specify multiple STATE keywords.
*/
class Ext_slat : public  Wavefunction
{

public:

  Ext_slat()
  {}


  virtual int nfunc()
  {
    return 1;
  }


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
  virtual void getDensity(Wavefunction_data *,int, Array2 <doublevar> &);

  virtual void saveUpdate(Sample_point *, int e, Wavefunction_storage *);
  virtual void restoreUpdate(Sample_point *, int e, Wavefunction_storage *);


  virtual void storeParmIndVal(Wavefunction_data *, Sample_point *,
                               int, Array1 <doublevar> & );
  virtual void getParmDepVal(Wavefunction_data *,
                             Sample_point *,
                             int,
                             Array1 <doublevar> &,
                             Array2 <doublevar> &);


  void generateStorage(Wavefunction_storage * & wfstore);


  void init(Wavefunction_data *);

  //--
private:


  void calcVal(Ext_slat_data *, Sample_point *);
  void updateVal(Ext_slat_data *, Sample_point *, int);
  void calcLap(Ext_slat_data *, Sample_point *);
  void updateLap(Ext_slat_data *, Sample_point *, int);

  Array1 <doublevar> electronIsStaleVal;
  Array1 <doublevar> electronIsStaleLap;
  int updateEverythingVal;
  int updateEverythingLap;

  int sampleAttached;
  int dataAttached;
  Ext_slat_data * parent;

  int staticSample;

  //Saved variables for electron updates
  Array3 <doublevar>  moVal;

  Array2 <doublevar> updatedMoVal;

  Array2 <doublevar> inverse;
  //!<inverse of the value part of the mo_values array transposed

  doublevar detVal;
  int nelectrons;

  
};

#endif //EXT_SLAT_H_INCLUDED
//--------------------------------------------------------------------------
