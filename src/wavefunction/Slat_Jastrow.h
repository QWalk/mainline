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

#ifndef SLAT_JASTROW_H_INCLUDED
#define SLAT_JASTROW_H_INCLUDED

#include "Wavefunction.h"
#include "Slat_Jastrow_data.h"

/*!
 
*/
class Slat_Jastrow_storage : public Wavefunction_storage
{
  Slat_Jastrow_storage()
  {
    //cout << "creating Slat_jastrow" << endl;
    slat_store=NULL;
    jast_store=NULL;
  }
  ~Slat_Jastrow_storage()
  {
    //cout << "destroying Slat_jastrow" << endl;
    if(slat_store)
      delete slat_store;
    if(jast_store)
      delete jast_store;
  }
private:
  friend class Slat_Jastrow;
  Wavefunction_storage * slat_store;
  Wavefunction_storage * jast_store;

};

/*!
While this class is named Slat_Jastrow, it can be used with any
wavefunction that's made from a multiplication of two others. There
shouldn't have to be any changes to do it.
*/
class Slat_Jastrow : public Wavefunction
{
public:

  Slat_Jastrow() : slater_wf(NULL), jastrow_wf(NULL)
  {}
  ~Slat_Jastrow()
  {
    deallocate(slater_wf);
    deallocate(jastrow_wf);
  }

  void init(Wavefunction_data *);

  int nfunc()
  {
    return nfunc_;
  }

  virtual void notify(change_type , int );

  virtual void updateVal(Wavefunction_data *, Sample_point *);
  virtual void updateLap(Wavefunction_data *, Sample_point *);
  virtual void updateForceBias(Wavefunction_data *, Sample_point *);


  virtual void getVal(Wavefunction_data *, int, Wf_return &);
  virtual void getLap(Wavefunction_data *, int, Wf_return &);
  virtual void getForceBias(Wavefunction_data *, int, Wf_return &);

  //virtual void getDensity(Wavefunction_data *,int,  Array2 <doublevar> &);
  virtual void evalTestPos(Array1 <doublevar> & pos, Sample_point * sample,Array1 <Wf_return> & wf);
  

  virtual void generateStorage(Wavefunction_storage * & wfstore);
  virtual void saveUpdate(Sample_point *, int e, Wavefunction_storage *);
  virtual void restoreUpdate(Sample_point *, int e, Wavefunction_storage *);

  //Added by Matous
  virtual void saveUpdate(Sample_point *, int e1, int e2, Wavefunction_storage *);
  virtual void restoreUpdate(Sample_point *, int e1, int e2, Wavefunction_storage *);

  
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

private:
  friend class Slat_Jastrow_data;
  Wavefunction * slater_wf;
  Wavefunction * jastrow_wf;
  int nfunc_;
};

#endif //SLAT_JASTROW_H_INCLUDED
//--------------------------------------------------------------------------
