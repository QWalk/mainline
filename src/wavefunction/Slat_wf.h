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

#ifndef SLAT_WF_H_INCLUDED

#define SLAT_WF_H_INCLUDED

#include "Qmc_std.h"
#include "Array.h"
#include "Wavefunction.h"
#include "MatrixAlgebra.h"
#include "MO_matrix.h"
#include "clark_updates.h"
class Wavefunction_data;
class Slat_wf_data;
class System;


//----------------------------------------------------------------------

template <class T> class Slat_wf_storage : public Wavefunction_storage
{
public:
  virtual ~Slat_wf_storage()
  {}
private:
  friend class Slat_wf<T>;

  //dimensions are [value gradient lap, MO]
  Array2 <T>  moVal_temp;

  // Added by Matous
  Array2 <T>  moVal_temp_2;
 
  Array3 < Array2 <T> > inverse_temp;
  Array3 <log_value<T> > detVal_temp;

};


/*!
A slater wavefunction; \f$\Psi=\sum_i det_i(\Phi_1\Phi_2...)\f$
where the \f$\Phi\f$'s are one-particle molecular orbitals.
Also supports multiple states as a vector of wavefunction values.  Just
specify multiple STATE keywords.
*/
template <class T> class Slat_wf : public  Wavefunction
{

public:

  Slat_wf()
  {}


  virtual int nfunc() {
    return nfunc_;
  }


  virtual void notify(change_type , int );

  virtual void updateVal(Wavefunction_data *, Sample_point *);
  virtual void updateLap(Wavefunction_data *, Sample_point *);

  virtual void getVal(Wavefunction_data *, int, Wf_return &);
  virtual void getLap(Wavefunction_data *, int, Wf_return &);
  virtual void evalTestPos(Array1 <doublevar> & pos, Sample_point *, Array1 <Wf_return> & wf);

  virtual void saveUpdate(Sample_point *, int e, Wavefunction_storage *);
  virtual void restoreUpdate(Sample_point *, int e, Wavefunction_storage *);

  // Added by Matous
  virtual void saveUpdate(Sample_point *, int e1, int e2, Wavefunction_storage *);
  virtual void restoreUpdate(Sample_point *, int e1, int e2, Wavefunction_storage *);
  

  virtual int getParmDeriv(Wavefunction_data *, 
			   Sample_point *,
			   Parm_deriv_return & );

  virtual void getSymmetricVal(Wavefunction_data *, 
			       int, 
			       Wf_return &);


  void generateStorage(Wavefunction_storage * & wfstore);


  void init(Wavefunction_data *,Templated_MO_matrix<T> * molecorb);

  //--
private:


  void updateInverse(Slat_wf_data *, int e);
  int updateValNoInverse(Slat_wf_data *, int e); 
  //!< update the value, but not the inverse.  Returns 0 if the determinant is zero and updates aren't possible
  
  void calcVal(Slat_wf_data *, Sample_point *);
  void updateVal(Slat_wf_data *, Sample_point *, int);
  void calcLap(Slat_wf_data *, Sample_point *);
  void updateLap(Slat_wf_data *, Sample_point *, int);
  void getDetLap(int e, Array3<log_value <T> > & vals );
  

  Array1 <int> electronIsStaleVal;
  Array1 <int> electronIsStaleLap;
  int updateEverythingVal;
  int updateEverythingLap;

  int sampleAttached;
  int dataAttached;
  Slat_wf_data * parent;
  Templated_MO_matrix<T> * molecorb;
  //lazy updates of the determinant(which saves a lot of time in pseudopotentials, etc)
  int inverseStale;
  int lastValUpdate;
  Array3<log_value<T> > lastDetVal;
  
  //Saved variables for electron updates
  Array3 <T>  moVal;

  Array2 <T> updatedMoVal;

  Array3 < Array2 <T> > inverse;
  //!<inverse of the value part of the mo_values array transposed

  Array3 <log_value<T> > detVal; //function #, determinant #, spin


  int nmo;        //!<Number of molecular orbitals
  int ndet;       //!<Number of determinants
  int ndim;       //!<Number of (spacial) dimensions each electron has
  int nfunc_;      //!<Number of functions this class represents.
  Array1 <int> nelectrons; //!< 2 elements, for spin up and down
  Array1 <int> spin;       //!< lookup table for the spin of a given electron
  Array1 <int> rede;       //!< number of the electron within its spin channel
  Array1 <int> opspin;


  Array2 <T> work1, work2; //!< Work matrices

};


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

#include "Qmc_std.h"
#include "MatrixAlgebra.h"
#include "Sample_point.h"
#include "Slat_wf_data.h"

//----------------------------------------------------------------------


template <class T> 
inline void Slat_wf<T>::generateStorage(Wavefunction_storage * & wfstore)
{
  wfstore=new Slat_wf_storage<T>;
  Slat_wf_storage<T> * store;
  recast(wfstore, store);
  store->moVal_temp.Resize (5,   nmo);
  
  // Added by Matous
  store->moVal_temp_2.Resize (5,   nmo);

  store->detVal_temp.Resize(nfunc_, ndet, 2);
  store->inverse_temp.Resize(nfunc_, ndet, 2);
  for(int i=0; i< nfunc_; i++)
  {
    for(int det=0; det < ndet; det++)
    {
      for(int s=0; s<2; s++)
      {
        store->inverse_temp(i,det,s).Resize(nelectrons(s), nelectrons(s));
        store->inverse_temp(i,det,s)=0;
        //store->detVal_temp(i,det,s)=1;
      }
    }
  }
}


//----------------------------------------------------------------------

template <class T> inline void Slat_wf<T>::init(Wavefunction_data * wfdata,
    Templated_MO_matrix<T> * orb)
{

  molecorb=orb;
  Slat_wf_data * dataptr;
  recast(wfdata, dataptr);
  recast(wfdata, parent);
  nfunc_=dataptr->nfunc;
  nmo=dataptr->nmo;
  ndet=dataptr->ndet;
  nelectrons.Resize(2);
  nelectrons=dataptr->nelectrons;

  int tote=nelectrons(0)+nelectrons(1);
  ndim=3;

  spin.Resize(tote);
  rede.Resize(tote);
  opspin.Resize(tote);
  for(int e=0; e < nelectrons(0); e++) {
    spin(e)=0;
    rede(e)=e;
    opspin(e)=1;
  }
  for(int e=nelectrons(0); e< nelectrons(0)+nelectrons(1); e++) {
    rede(e)=e-nelectrons(0);
    spin(e)=1;
    opspin(e)=0;
  }
  //Properties and intermediate calculation storage.
  moVal.Resize(5,   tote, nmo);
  updatedMoVal.Resize(nmo,5);

  detVal.Resize (nfunc_, ndet, 2);
  inverse.Resize(nfunc_, ndet, 2);

  for(int i=0; i< nfunc_; i++) {
    for(int det=0; det < ndet; det++) {
      for(int s=0; s<2; s++) {
        inverse(i,det,s).Resize(nelectrons(s), nelectrons(s));
        inverse(i,det,s)=0;
        for(int e=0; e< nelectrons(s); e++) {
          inverse(i,det,s)(e,e)=1;
          inverse(i,det,s)(e,e)=1;
        }

        detVal(i,det,s)=T(1.0);
      }
    }
  }


  electronIsStaleVal.Resize(tote);
  electronIsStaleLap.Resize(tote);

  electronIsStaleVal=0;
  electronIsStaleLap=0;
  updateEverythingVal=1;
  updateEverythingLap=1;
  sampleAttached=0;
  dataAttached=0;
  
  inverseStale=0;
  lastValUpdate=0;
}

//----------------------------------------------------------------------

/*!
 */
template<class T> inline void Slat_wf<T>::notify(change_type change, int num)
{
  switch(change) {
    case electron_move:
      electronIsStaleVal(num)=1;
      electronIsStaleLap(num)=1;
      break;
    case all_electrons_move:
      updateEverythingVal=1;
      updateEverythingLap=1;
      break;
    case wf_parm_change:
    case all_wf_parms_change:
      if(parent->optimize_mo  ) {
        updateEverythingVal=1;
        updateEverythingLap=1;
      }
      break;
    case sample_attach:
      sampleAttached=1;
      updateEverythingVal=1;
      updateEverythingLap=1;
      break;
    case data_attach:
      dataAttached=1;
      updateEverythingVal=1;
      updateEverythingLap=1;
      break;
    default:
      updateEverythingVal=1;
      updateEverythingLap=1;
  }
}


//----------------------------------------------------------------------

template<class T>inline void Slat_wf<T>::saveUpdate(Sample_point * sample, int e,
                                                    Wavefunction_storage * wfstore) {
  
  Slat_wf_storage<T> * store;
  recast(wfstore, store);
  
  //presumably, if we care enough to save the update, we care enough
  //to have the inverse up to date
  if(inverseStale) {
    detVal=lastDetVal;
    updateInverse(parent, lastValUpdate);
    inverseStale=0;
  }
  int s=spin(e);
  
  int ndet_save=ndet;
  if(parent->use_clark_updates) ndet_save=1;
  for(int f=0; f< nfunc_; f++) {
    for(int det=0; det<ndet_save; det++) {
      store->inverse_temp(f,det,s)=inverse(f,det,s);
    }
    for(int det=0; det < ndet; det++) {
      store->detVal_temp(f,det,s)=detVal(f,det,s);
    }
  }
  
  
  int norb=moVal.GetDim(2);
  for(int d=0; d< 5; d++) {
    for(int i=0; i< norb; i++) {
      store->moVal_temp(d,i)=moVal(d,e,i);
    }
  }
  

}

//----------------------------------------------------------------------

template<class T>inline void Slat_wf<T>::restoreUpdate(Sample_point * sample, int e,
                            Wavefunction_storage * wfstore) {

  Slat_wf_storage<T> * store;
  recast(wfstore, store);
  int s=spin(e);
  inverseStale=0;
  
  for(int j=0; j<5; j++) {
    for(int i=0; i<moVal.GetDim(2); i++) {
      moVal(j,e,i)=store->moVal_temp(j,i);
    }
  }
  int ndet_save=ndet;
  if(parent->use_clark_updates) ndet_save=1;
  
  for(int f=0; f< nfunc_; f++) {
    for(int det=0; det < ndet_save; det++) {
      inverse(f,det,s)=store->inverse_temp(f,det,s);
    }
    for(int det=0; det < ndet; det++) {
      detVal(f,det,s)=store->detVal_temp(f,det,s);
    }
  }
  //It seems to be faster to update the inverse than to save it and
  //recover it.  However, it complicates the implementation too much.
  //For now, we'll disable it.
  //updateInverse(parent,e);
  
  electronIsStaleVal(e)=0;
  electronIsStaleLap(e)=0;

}

//----------------------------------------------------------------------

// Added by Matous
template <class T>inline void Slat_wf<T>::saveUpdate(Sample_point * sample, int e1, int e2,
                         Wavefunction_storage * wfstore) {

  Slat_wf_storage<T> * store;
  recast(wfstore, store);
  
  //presumably, if we care enough to save the update, we care enough
  //to have the inverse up to date
  if(inverseStale) {
    detVal=lastDetVal;
    updateInverse(parent, lastValUpdate);
    inverseStale=0;
  }
  
  int s1=spin(e1), s2=spin(e2);
  
  for(int f=0; f< nfunc_; f++) {
    for(int det=0; det<ndet; det++) {
      if ( s1 == s2 ) {
        store->inverse_temp(f,det,s1)=inverse(f,det,s1);
        store->detVal_temp(f,det,s1)=detVal(f,det,s1);
      }
      else {
        store->inverse_temp(f,det,s1)=inverse(f,det,s1);
        store->inverse_temp(f,det,s2)=inverse(f,det,s2);
        store->detVal_temp(f,det,s1)=detVal(f,det,s1);
        store->detVal_temp(f,det,s2)=detVal(f,det,s2);
      }
    }
  }
  
  
  for(int d=0; d< 5; d++) {
    for(int i=0; i< moVal.GetDim(2); i++) {
      store->moVal_temp(d,i)=moVal(d,e1,i);
      store->moVal_temp_2(d,i)=moVal(d,e2,i);
    }
  }
  

}

//----------------------------------------------------------------------

// Added by Matous
template<class T> inline void Slat_wf<T>::restoreUpdate(Sample_point * sample, int e1, int e2,
                            Wavefunction_storage * wfstore)
{

  Slat_wf_storage<T> * store;
  recast(wfstore, store);
  
  int s1=spin(e1), s2=spin(e2);
  inverseStale=0;
  
  for(int j=0; j<5; j++) {
    for(int i=0; i<moVal.GetDim(2); i++) {
      moVal(j,e1,i)=store->moVal_temp(j,i);
      moVal(j,e2,i)=store->moVal_temp_2(j,i);
    }
  }
  for(int f=0; f< nfunc_; f++) {
    for(int det=0; det < ndet; det++) {
      if ( s1 == s2 ) {
		      inverse(f,det,s1)=store->inverse_temp(f,det,s1);
		      detVal(f,det,s1)=store->detVal_temp(f,det,s1);
      }
      else {
		      inverse(f,det,s1)=store->inverse_temp(f,det,s1);
		      inverse(f,det,s2)=store->inverse_temp(f,det,s2);
		      detVal(f,det,s1)=store->detVal_temp(f,det,s1);
		      detVal(f,det,s2)=store->detVal_temp(f,det,s2);
      }
    }
  }
  
  electronIsStaleVal(e1)=0;
  electronIsStaleLap(e1)=0;
  electronIsStaleVal(e2)=0;
  electronIsStaleLap(e2)=0;

}

//----------------------------------------------------------------------

template <class T> inline void Slat_wf<T>::updateVal(Wavefunction_data * wfdata,
                        Sample_point * sample)
{

  assert(sampleAttached);
  assert(dataAttached);

  if(updateEverythingVal==1) {
    calcVal(parent, sample);
    updateEverythingVal=0;
    electronIsStaleVal=0;
  }
  else {
    for(int e=0; e< nelectrons(0)+nelectrons(1); e++) {
      if(electronIsStaleVal(e)) {
        updateVal(parent, sample, e);
        electronIsStaleVal(e)=0;
      }
    }
  }

}

//----------------------------------------------------------------------

template <class T> inline void Slat_wf<T>::updateLap( Wavefunction_data * wfdata,
                        Sample_point * sample)
{
  assert(sampleAttached);
  assert(dataAttached);
  
  if(updateEverythingLap==1) {
    calcLap(parent, sample);
    updateEverythingVal=0;
    updateEverythingLap=0;
    electronIsStaleLap=0;
    electronIsStaleVal=0;
  }
  else {
    for(int e=0; e< nelectrons(0)+nelectrons(1); e++) {
      if(electronIsStaleLap(e)) {
        updateLap(parent, sample, e);
        electronIsStaleLap(e)=0;
        electronIsStaleVal(e)=0;
      }
    }
    
  }

}




//-----------------------------------------------------------------------


template <> inline int Slat_wf<dcomplex>::getParmDeriv(Wavefunction_data *  wfdata, 
			  Sample_point * sample ,
			  Parm_deriv_return & derivatives){
  error("parmderiv not supported for complex orbitals yet");
  return 0;
}

template <> inline int Slat_wf<doublevar>::getParmDeriv(Wavefunction_data *  wfdata, 
			  Sample_point * sample ,
			  Parm_deriv_return & derivatives){
// error("parmderiv not supported yet!"); 
  
  if(inverseStale) { 
    detVal=lastDetVal;
    updateInverse(parent, lastValUpdate);
    inverseStale=0;
  }
  
  int nparms_start=derivatives.nparms_start;
  int nparms_end=derivatives.nparms_end;
  int nparms=nparms_end-nparms_start;
  int tote=nelectrons(0)+nelectrons(1);
  
  derivatives.gradient.Resize(nparms);
  derivatives.hessian.Resize(nparms, nparms);
  derivatives.gradderiv.Resize(nparms,tote,4);
  derivatives.val_gradient.Resize(tote,4);
  if(parent->optimize_mo) {
    error("don't support optimizing mo yet");
  }
  else if(parent->optimize_det) {
    log_value<doublevar> detsum=0;
    Array1 <log_value<doublevar> > detvals(ndet);
    Array3 <log_value <doublevar> > detgrads(ndet,tote,5);
    Array2 <log_value <doublevar> > totgrads(tote,5);
    for(int det=0; det < ndet; det++) {
      log_value<doublevar> thisdet=detVal(0,det,0)*detVal(0,det,1);
      detvals(det)=parent->detwt(det)*thisdet;
    }

    Array3 <log_value<doublevar> > tmp_detgrads;
    for(int e=0; e< tote; e++) { 
      getDetLap(e,tmp_detgrads);
      for(int det=0; det < ndet; det++) { 
        for(int d=1; d< 5; d++) {
          detgrads(det,e,d)=tmp_detgrads(0,det,d);
        }
      }
    }

    
    detsum=sum(detvals);
    detsum.logval*=-1;
    derivatives.gradient=0.0;
    derivatives.gradderiv=0.0;
    //---------------  Assign the value gradient
    for(int e=0; e< tote; e++) {
      for(int d=1; d< 5; d++) { 

        Array1 <log_value<doublevar> > tmpgrad(ndet);
        for(int det=0; det < ndet; det++) 
          tmpgrad(det)=parent->detwt(det)*detgrads(det,e,d);
        totgrads(e,d)=sum(tmpgrad);
        totgrads(e,d)*=detsum;
      }
    }

    for(int e=0; e< tote; e++) {
      for(int d=1; d< 5; d++) {
        derivatives.val_gradient(e,d-1)=totgrads(e,d).val();
      }
    }

    //---------------
    int counter=0;
    derivatives.gradderiv=0.0;
    for(int csf=0; csf < parent->ncsf; csf++) { 
      if(csf > nparms_start && csf <= nparms_end ){
        for(int j=1;j<parent->CSF(csf).GetDim(0);j++){
          int det=counter;
          doublevar coeff=parent->CSF(csf)(j);
          int index=csf-1-nparms_start;
          log_value<doublevar> thisdet=detVal(0,det,0)*detVal(0,det,1);
          derivatives.gradient(index)+=
            coeff*thisdet.val();
          for(int e=0; e< tote; e++) { 
            for(int d=1; d< 5; d++) { 
              derivatives.gradderiv(index,e,d-1)+=coeff
                *(detgrads(det,e,d).val()-totgrads(e,d).val()*thisdet.val())*detsum.val();
              //cout << "coeff " << coeff << " detgrad " << detgrads(det,e,d).val()
              //  << " totgrad " << totgrads(e,d).val() << " thisdet " << thisdet.val()
              //  << " detsum " << detsum.val() << endl;
              //cout << "deriv " << derivatives.gradderiv(index,e,d-1) << endl;
            }
          }
          counter++;
        }
      }
      else{
        counter+=parent->CSF(csf).GetDim(0)-1;
      }
    }
    for(int csf=0; csf< nparms; csf++) {
      derivatives.gradient(csf)*=detsum.val(); 
    }
    derivatives.hessian=0;
    //for(int csf=0; csf < parent->ncsf; csf++) { 
    //  for(int e=0; e< tote; e++) { 
    //    for(int d=0; d< 4; d++) { 
    //      cout << "deriv " << derivatives.gradderiv(csf,e,d) << endl;
    //    }
    //  }
    //}
    return 1;
  }
  else { 
    derivatives.gradient=0;
    derivatives.hessian=0;
    return 1;
  }
  
  return 0;
}


//------------------------------------------------------------------------


template <class T> inline void Slat_wf<T>::calcVal(Slat_wf_data * dataptr, Sample_point * sample)
{
  //Hmm, I don't completely understand why, but something is not 
  //completely clean, so we can't just cycle through and updateVal
  //This is actually probably the best way to do it anyway, since it 
  //should in theory be faster, and it gives us a clean start
  calcLap(dataptr, sample);

}

//------------------------------------------------------------------------
inline doublevar real_qw(doublevar & a) { return a; } 
inline doublevar real_qw(dcomplex & a) { return real(a); } 

template <class T>inline void Slat_wf<T>::updateInverse(Slat_wf_data * dataptr, int e) { 
  int maxmatsize=max(nelectrons(0),nelectrons(1));
  Array1 <T> modet(maxmatsize);
  int s=spin(e);
  int ndet_update=ndet;
  if(parent->use_clark_updates) ndet_update=1;
  for(int f=0; f< nfunc_; f++)  {
    for(int det=0; det< ndet_update; det++)  {
      //fill the molecular orbitals for this
      //determinant
      if(real_qw(detVal(f,det,s).logval) < -1e200) { 
        Array2 <T> allmos(nelectrons(s), nelectrons(s));
        for(int e=0; e< nelectrons(s); e++) {
          int curre=s*nelectrons(0)+e;
          for(int i=0; i< nelectrons(s); i++) {
            allmos(e,i)=moVal(0,curre, dataptr->occupation(f,det,s)(i));
          }
        }

#ifdef SUPERDEBUG
        cout << "Slat_wf::updateInverse: near-zero determinant " 
          << " f " << f << " det " << det << " old det " << detVal(f,det,s).logval
          << endl;
#endif


        detVal(f,det,s)=
          TransposeInverseMatrix(allmos,inverse(f,det,s), nelectrons(s));
#ifdef SUPERDEBUG
        cout << "Slat_wf::updateInverse: near-zero determinant " 
          << " f " << f << " det " << det << " new det " << detVal(f,det,s).logval
          << endl;
#endif

        

      }
      else { 
        for(int i = 0; i < nelectrons(s); i++) {
          modet(i)=moVal(0,e,dataptr->occupation(f,det,s)(i));
        }
        T ratio=1./InverseUpdateColumn(inverse(f,det,s),
            modet, rede(e),
            nelectrons(s));

        detVal(f,det, s)=ratio*detVal(f,det, s);
      }
    }
    if(parent->use_clark_updates) { 

      Array2 <T> M(nelectrons(s),updatedMoVal.GetDim(0));
      for(int i=0; i< nelectrons(s); i++){ 
        int elec=i+s*nelectrons(0);
        for(int j=0; j< updatedMoVal.GetDim(0); j++) { 
          M(i,j)=moVal(0,elec,j);
        }
      }
      Array1 <T> ratios;
      parent->excitations.clark_updates(inverse(0,0,s),M,s,ratios);
      for(int d=0; d< ndet; d++) { 
        detVal(0,d,s)=ratios(d)*detVal(0,0,s); 
      }

    }

  }
  
}

//------------------------------------------------------------------------

template <class T> inline int Slat_wf<T>::updateValNoInverse(Slat_wf_data * dataptr, int e) { 
  int maxmatsize=max(nelectrons(0),nelectrons(1));
  Array1 <T> modet(maxmatsize);
  int s=spin(e);
  for(int f=0; f< nfunc_; f++)  {
    for(int det=0; det< ndet; det++)  {
      //fill the molecular orbitals for this
      //determinant
      if(real_qw(detVal(f,det,s).logval) < -1e200) return 0;
    }
  }
  
  
  for(int f=0; f< nfunc_; f++)  {
    for(int det=0; det< ndet; det++)  {
      //fill the molecular orbitals for this
      //determinant
      for(int i = 0; i < nelectrons(s); i++) {
        modet(i)=moVal(0,e,dataptr->occupation(f,det,s)(i));
      }
      
      
      T ratio=1./InverseGetNewRatio(inverse(f,det,s),
                                            modet, rede(e),
                                            nelectrons(s));
#ifdef SUPERDEBUG
      T tmpratio=InverseGetNewRatio(inverse(f,det,s),
                                            modet, rede(e),
                                            nelectrons(s));

      cout << "Slat_wf::updateValNoInverse: " << "ratio " << ratio 
        << " inv ratio " << tmpratio << " old detVal " << detVal(f,det,s).logval << endl;
#endif

      detVal(f,det, s)=ratio*detVal(f,det, s);
      
    }
  }
  return 1;
}

//------------------------------------------------------------------------
/*!

*/
template <class T> inline void Slat_wf<T>::updateVal( Slat_wf_data * dataptr, Sample_point * sample,int e) {

  if(inverseStale && lastValUpdate!=e) { 
    inverseStale=0;
    detVal=lastDetVal;
    updateInverse(dataptr, lastValUpdate);
  }
  if(inverseStale && lastValUpdate==e) { 
    inverseStale=0;
    detVal=lastDetVal;
  }

  assert(dataptr != NULL);
  sample->updateEIDist();
  int s=spin(e);

  //update all the mo's that we will be using.
  molecorb->updateVal(sample, e, s,
                              updatedMoVal);

  for(int i=0; i< updatedMoVal.GetDim(0); i++)
    moVal(0,e,i)=updatedMoVal(i,0);


  inverseStale=1;
  lastValUpdate=e;
  lastDetVal=detVal;
  if(!parent->use_clark_updates) { 
    if(!updateValNoInverse(dataptr, e)) { 
      inverseStale=0;
      updateInverse(dataptr,e);
    }
  }
  else { 
    updateInverse(dataptr,e);
    inverseStale=0;
  } 
//  for(int d=0; d< ndet; d++) { 
//    cout << "orig " << detVal(0,d,s).val()  
//       << " update " << ratios(d)*detVal(0,0,s).val()
//       << " ratio " << detVal(0,d,s).val()/detVal(0,0,s).val() 
//       << " computed ratio " << ratios(d) << endl;
//  }
  //----------------------------------
}

//------------------------------------------------------------------------


template <class T>inline void Slat_wf<T>::getVal(Wavefunction_data * wfdata, int e,
                     Wf_return & val)
{
  //Array1 <doublevar> si(nfunc_, 0.0);
  Array2 <log_value<T> > vals(nfunc_,1,T(0.0));

  assert(val.amp.GetDim(0) >=nfunc_);
  assert(val.amp.GetDim(1) >= 1);
 
  
  Slat_wf_data * dataptr;
  recast(wfdata, dataptr);
  
  int s=spin(e);
  int opp=opspin(e);
  
  for(int f=0; f< nfunc_; f++) {
    Array1 <log_value<T> > detvals(ndet);
    for(int det=0; det < ndet; det++) {
      detvals(det) = dataptr->detwt(det)*detVal(f,det,s)*detVal(f,det,opp);
    }
    log_value<T> totval=sum(detvals);
    //vals(f,0)=totval.logval;
    //si(f)=totval.sign;
    vals(f,0)=totval;
  }

  val.setVals(vals);

}

//----------------------------------------------------------
template <class T> inline void Slat_wf<T>::getSymmetricVal(Wavefunction_data * wfdata,
		     int e, Wf_return & val){
  val.phase(0, 0)=0;
  val.amp(0, 0)=0;
  val.cvals(0,0)=0;
} 

//----------------------------------------------------------------------



template <class T> inline void Slat_wf<T>::calcLap(Slat_wf_data * dataptr, Sample_point * sample)
{
  //cout << "calcLap " << endl;
  inverseStale=0;
  for(int e=0; e< nelectrons(0)+nelectrons(1); e++)  {
    int s=spin(e);
    sample->updateEIDist();

    //update all the mo's that we will be using, using the lists made in
    //Slat_wf_data(one for each spin).
    //cout << "mo_updatelap " << endl;
    molecorb->updateLap(sample, e, s,
                                updatedMoVal);
    //cout << "done " << endl;
    for(int d=0; d< 5; d++)  {
      for(int i=0; i< updatedMoVal.GetDim(0); i++) {
        moVal(d,e,i)=updatedMoVal(i,d);
      }
    }
  }

  int maxmatsize=max(nelectrons(0),nelectrons(1));
  Array2 <T> modet(maxmatsize, maxmatsize);
  //ofstream matout("matrix_out", ios::app);
  //matout.precision(15);
  //matout << "initial_matrix " << nelectrons(0) << " rows are electrons, columns are orbital values " << endl;
  //cout << "here " << endl;
  for(int f=0; f< nfunc_; f++)   {
    for(int det=0; det < ndet; det++ ) {
      for(int s=0; s< 2; s++ ) {


        for(int e=0; e< nelectrons(s); e++) {
          int curre=s*nelectrons(0)+e;
          for(int i=0; i< nelectrons(s); i++) {
            modet(e,i)=moVal(0,curre, dataptr->occupation(f,det,s)(i));
            //if(f==0 && det==0 && s==0) { 
            //  matout << modet(e,i) << " ";
            //}
          }
          //if(f==0 && det==0 && s==0) matout << endl;
        }
        
        if(nelectrons(s) > 0) { 
          detVal(f,det,s)=
          TransposeInverseMatrix(modet,inverse(f,det,s), nelectrons(s));
        }
        else detVal(f,det,s)=T(1.0);
#ifdef SUPERDEBUG
        cout << "Slat_wf::calcLap: f " << f<< " det " << det 
          << " detVal " << detVal(f,det,s).logval << endl;
#endif
        //if(f==0 && det==0 && s==0) matout << "determinant " << detVal(f,det,s) 
         //   << " should be " << modet(0,0)*modet(1,1)-modet(0,1)*modet(1,0) << endl;
      }
    }
  }
  //cout << "done " << endl;
}

//------------------------------------------------------------------------


template <class T> void Slat_wf<T>::getDetLap(int e, Array3<log_value <T> > &  vals ) { 
  vals.Resize(nfunc_,ndet,5);
  int s=spin(e);
  int opp=opspin(e);

  //Prepare the matrices we need for the inverse.
  //There is likely a way to do this via updates, but we'll
  //leave it for now since it doesn't seem to cost too much.
  Array2 <T> & lapvec=work2;
  Array1 <T> tmplapvec(nelectrons(s));
  int n=moVal.GetDim(2);
  if(parent->use_clark_updates) { 
    lapvec.Resize(nelectrons(s),n);
    for(int e1=0; e1< nelectrons(s); e1++) {
      int shift=e1+s*nelectrons(0);
      for(int j=0; j< n; j++) {  
        lapvec(e1,j)=moVal(0,shift,j);
      }
    }
  }


  for(int f=0; f< nfunc_; f++) {
    Array1 <log_value <T> > detvals(ndet);
    for(int det=0; det < ndet; det++) {
      detvals(det) = detVal(f,det,s)*detVal(f,det,opp);

      vals(f,det,0)=detvals(det);
    }
    
    Array1 <log_value <T> > detgrads(ndet);
    for(int i=1; i< 5; i++) {
      if(!parent->use_clark_updates) {   //Sherman-Morrison updates
        for(int det=0; det < ndet; det++) {
          T temp=0;
          for(int j=0; j<nelectrons(s); j++) {
            temp+=moVal(i , e, parent->occupation(f,det,s)(j) )
              *inverse(f,det,s)(rede(e), j);
          }
          detgrads(det)=temp; 
          detgrads(det)*=detVal(f,det,s);
        }
      } //-------
      else {  //clark updates

        Array2 <T> &  tmpinverse=work1;
        tmpinverse=inverse(f,0,s);

        for(int j=0; j< n; j++) 
          lapvec(rede(e),j)=moVal(i,e,j);

        for(int j=0; j< nelectrons(s); j++) 
          tmplapvec(j)=lapvec(rede(e),parent->occupation(f,0,s)(j));

        T baseratio=1.0/InverseUpdateColumn(tmpinverse,tmplapvec,
            rede(e),nelectrons(s));
        Array1 <T> ratios;
        parent->excitations.clark_updates(tmpinverse,lapvec,s,ratios);
        detgrads(0)=baseratio*detVal(f,0,s);
        for(int d=1; d< ndet; d++) { 
          //detgrads(d)=baseratio*ratios(d)*detVal(f,0,s);
          detgrads(d)=ratios(d)*detgrads(0);
        }
      } //------Done clark updates
      for(int d=0; d< ndet; d++) {
        detgrads(d)*=detVal(f,d,opp);
      }

      //--------------------------------
      for(int d=0; d< ndet; d++) {
        vals(f,d,i)=detgrads(d);
      }
    }

  }
}

/*!
*/

template <class T> void Slat_wf<T>::getLap(Wavefunction_data * wfdata,
                     int e, Wf_return & lap)
{

  //Array1 <doublevar> si(nfunc_, 0.0);
  //Array2 <doublevar> vals(nfunc_,5,0.0);
  Array2 <log_value <T> > vals(nfunc_,5);
  
  Array3 <log_value<T> > detvals;
  getDetLap(e,detvals);
  Array1 <log_value<T> > tempsum(ndet);
  for(int f=0; f< nfunc_; f++) {
    for(int i=0; i< 5; i++) {
      for(int d=0;d < ndet; d++) {
        tempsum(d)=parent->detwt(d)*detvals(f,d,i);
      }
      vals(f,i)=sum(tempsum);
    }
    log_value<T> inv=vals(f,0);
    inv.logval*=-1;
    for(int i=1; i< 5; i++) vals(f,i)*=inv;
  }
  
  
  lap.setVals(vals);
  
}

//-------------------------------------------------------------------------

/*!
*/
template <class T> inline void Slat_wf<T>::updateLap(Slat_wf_data * dataptr,
                        Sample_point * sample,
                        int e ) {
  
  if(inverseStale && lastValUpdate!=e) { 
    inverseStale=0;
    detVal=lastDetVal;
    updateInverse(dataptr, lastValUpdate);
  }
  if(inverseStale && lastValUpdate==e) { 
    detVal=lastDetVal;
    inverseStale=0;
  }
  assert(dataptr != NULL);

  int s=spin(e);
  sample->updateEIDist();


  //update all the mo's that we will be using.
  molecorb->updateLap(sample,e,s,updatedMoVal);

  for(int d=0; d< 5; d++)
    for(int i=0; i< updatedMoVal.GetDim(0); i++)
      moVal(d,e,i)=updatedMoVal(i,d);
  
  updateInverse(dataptr,e);
}

//-------------------------------------------------------------------------

template <class T> inline void Slat_wf<T>::evalTestPos(Array1 <doublevar> & pos, 
    Sample_point * sample, Array1 <Wf_return> & wf) {
  
  if(inverseStale) { 
    inverseStale=0;
    detVal=lastDetVal;
    updateInverse(parent, lastValUpdate);
  }

  int nspin=2;
  Array1<Array2 <T> > movals(nspin);
  Array1 <doublevar> oldpos(ndim);
  Array1 <T> modet(nmo);
  
  sample->getElectronPos(0,oldpos);
  for(int s=0; s< nspin; s++) {
    movals(s).Resize(nmo,1);
    sample->setElectronPosNoNotify(0,pos);
    sample->updateEIDist();
    molecorb->updateVal(sample,0,s,movals(s));
  }
  sample->setElectronPosNoNotify(0,oldpos);

  int tote=sample->electronSize();
  wf.Resize(tote);

  for(int e=0; e< tote; e++) { 
    wf(e).Resize(nfunc_,1);
    Array2 <log_value<T> > vals(nfunc_,1,T(0.0));
     
    int s=spin(e);
    int opps= s==0?1:0;
    Array1 <log_value <T> > new_detVals(ndet);
    int f=0;
    if(!parent->use_clark_updates) {  //Sherman-morrison updates
      for(int det=0; det< ndet; det++)  {
        for(int i = 0; i < nelectrons(s); i++) {
          modet(i)=movals(s)(parent->occupation(f,det,s)(i),0);
        }
        T ratio=1./InverseGetNewRatio(inverse(f,det,s),
            modet, rede(e),
            nelectrons(s));
        new_detVals(det)=parent->detwt(det)*detVal(f,det,s);
        new_detVals(det)*=ratio;
        new_detVals(det)*=detVal(f,det,opps);
      }
    }
    else { //Clark updates 
      Array2 <T> & motmp=work1;
      Array2 <T> & invtmp=work2;
      invtmp=inverse(f,0,s);
      int n=moVal.GetDim(2);
      motmp.Resize(nelectrons(s),n);
      for(int e1=0; e1< nelectrons(s); e1++) {
        int shift=e1+s*nelectrons(0);
        for(int j=0; j< n; j++) {  
          motmp(e1,j)=moVal(0,shift,j);
        }
      }
      for(int j=0; j< n; j++) motmp(rede(e),j)=movals(s)(j,0);
      Array1 <T> tmpvec(nelectrons(s));
      for(int j=0; j< nelectrons(s); j++) 
        tmpvec(j)=movals(s)(parent->occupation(f,0,s)(j),0);
      T baseratio=1.0/InverseUpdateColumn(invtmp,tmpvec,rede(e),nelectrons(s));
      Array1 <T> ratios;
      parent->excitations.clark_updates(invtmp,motmp,s,ratios);
      for(int d=0; d< ndet; d++) {
        new_detVals(d)=parent->detwt(d)*detVal(f,0,s);
        new_detVals(d)*=baseratio*ratios(d);
        new_detVals(d)*=detVal(f,d,opps);
      }
      
    }
    log_value<T> totval=sum(new_detVals);
    vals(f,0)=totval;
    wf(e).setVals(vals);

  }

  

}


//----------------------------------------------------------------------
#endif //SLAT_WF_H_INCLUDED
//--------------------------------------------------------------------------
