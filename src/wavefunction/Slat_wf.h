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
class Wavefunction_data;
class Slat_wf_data;
class System;


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


  void init(Wavefunction_data *,Templated_MO_matrix<T> * molecorb);

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

  //Variables for a static(electrons not moving) calculation
  int staticSample;
  Array3 <T> saved_laplacian;
  //!<Saved laplacian for a static calculation (electron, function, [val grad lap])


  int nmo;        //!<Number of molecular orbitals
  int ndet;       //!<Number of determinants
  int ndim;       //!<Number of (spacial) dimensions each electron has
  int nfunc_;      //!<Number of functions this class represents.
  Array1 <int> nelectrons; //!< 2 elements, for spin up and down
  Array1 <int> spin;       //!< lookup table for the spin of a given electron

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
  spin=dataptr->spin;
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
  staticSample=0;
  
  inverseStale=0;
  lastValUpdate=0;
}

//----------------------------------------------------------------------

/*!
Behavior under staticSample:  if sample_static is set, then
we ignore all electron moves in the main algorithm, and the only way
to update based on a move is by the getParmDepVal function.  The
regular update functions will not work in this case.
*/
template<class T> inline void Slat_wf<T>::notify(change_type change, int num)
{
  switch(change)
  {
  case electron_move:
    if(staticSample==0) {
      electronIsStaleVal(num)=1;
      electronIsStaleLap(num)=1;
    }
    break;
  case all_electrons_move:
    if(staticSample==0) {
      updateEverythingVal=1;
      updateEverythingLap=1;
    }
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
  case sample_static:
    if(!parent->optimize_mo) {
      save_for_static();
      staticSample=1;
    }
    break;
  case sample_dynamic:
    staticSample=0;
    updateEverythingVal=1;
    updateEverythingLap=1;
    break;
  default:
    updateEverythingVal=1;
    updateEverythingLap=1;
  }
}

//----------------------------------------------------------------------

template<class T> inline void Slat_wf<T>::save_for_static() {
  assert(staticSample==0);

  if(!parent->optimize_mo && !parent->optimize_det) {
    
    int totelectrons=nelectrons(0)+nelectrons(1);
    saved_laplacian.Resize(totelectrons, nfunc_, 6);
    Wf_return templap(nfunc_, 5);
    
    for(int e=0; e< totelectrons; e++) {
      getLap(parent, e, templap);
      for(int f=0; f< nfunc_; f++) {
      for(int i=0; i< 5; i++) saved_laplacian(e,f,i)=templap.amp(f,i);
      //saved_laplacian(e,f,5)=templap.phase(f,0);
      saved_laplacian(e,f,5)=templap.sign(f);
      }
    }
  }

}

//----------------------------------------------------------------------

template<class T>inline void Slat_wf<T>::saveUpdate(Sample_point * sample, int e,
                         Wavefunction_storage * wfstore)
{

  if(staticSample==0) {
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

    for(int f=0; f< nfunc_; f++) {
      for(int det=0; det<ndet; det++) {
        store->inverse_temp(f,det,s)=inverse(f,det,s);
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


}

//----------------------------------------------------------------------

template<class T>inline void Slat_wf<T>::restoreUpdate(Sample_point * sample, int e,
                            Wavefunction_storage * wfstore)
{

  if(staticSample==0) {
    Slat_wf_storage<T> * store;
    recast(wfstore, store);
    int s=spin(e);
    inverseStale=0;
    
    for(int j=0; j<5; j++) {
      for(int i=0; i<moVal.GetDim(2); i++) {
        moVal(j,e,i)=store->moVal_temp(j,i);
      }
    }
    for(int f=0; f< nfunc_; f++) {
      for(int det=0; det < ndet; det++) {
        inverse(f,det,s)=store->inverse_temp(f,det,s);
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

}

//----------------------------------------------------------------------

// Added by Matous
template <class T>inline void Slat_wf<T>::saveUpdate(Sample_point * sample, int e1, int e2,
                         Wavefunction_storage * wfstore)
{

  if(staticSample==0) {
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


}

//----------------------------------------------------------------------

// Added by Matous
template<class T> inline void Slat_wf<T>::restoreUpdate(Sample_point * sample, int e1, int e2,
                            Wavefunction_storage * wfstore)
{

  if(staticSample==0) {
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

}

//----------------------------------------------------------------------

template <class T> inline void Slat_wf<T>::updateVal(Wavefunction_data * wfdata,
                        Sample_point * sample)
{

  assert(sampleAttached);
  assert(dataAttached);

  Slat_wf_data * slatdata;
  recast(wfdata, slatdata);

  if(staticSample==0 || parent->optimize_mo ) {
    if(updateEverythingVal==1) {
      calcVal(slatdata, sample);
      updateEverythingVal=0;
      electronIsStaleVal=0;
    }
    else {
      for(int e=0; e< nelectrons(0)+nelectrons(1); e++) {
        if(electronIsStaleVal(e)) {
          updateVal(slatdata, sample, e);
          electronIsStaleVal(e)=0;
        }
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

  Slat_wf_data * slatdata;
  recast(wfdata, slatdata);
  if(staticSample==0 || parent->optimize_mo ) {
    if(updateEverythingLap==1) {
      calcLap(slatdata, sample);
      updateEverythingVal=0;
      updateEverythingLap=0;
      electronIsStaleLap=0;
      electronIsStaleVal=0;
    }
    else {
      for(int e=0; e< nelectrons(0)+nelectrons(1); e++) {
        if(electronIsStaleLap(e)) {
          assert(!staticSample);
          updateLap(slatdata, sample, e);
          electronIsStaleLap(e)=0;
          electronIsStaleVal(e)=0;
        }
      }

    }
  }

}

//----------------------------------------------------------------------


template <class T> inline void Slat_wf<T>::storeParmIndVal(Wavefunction_data * wfdata, Sample_point * sample,
                              int e, Array1 <doublevar> & vals )
{

  if(parent->optimize_mo) {

  }
  else { error("parmindval not implemented yet!"); } 
  /*
  else if(parent->optimize_det) {
    assert(vals.GetDim(0) >= 2*parent->detwt.GetDim(0));
    updateVal(wfdata, sample);
    int count=0;
    for(int det=0; det < ndet; det++) {
      vals(count++)=detVal(0,det,0).val();
      vals(count++)=detVal(0,det,1).val();
    }
  }
  else {
    assert(vals.GetDim(0) >=2);
    Wf_return newval(nfunc_,1);
    updateVal(wfdata, sample);
    getVal(wfdata, e, newval);
    vals(0)=newval.amp(0,0);
    vals(1)=newval.phase(0,0);
  }
  */
}

//----------------------------------------------------------------------

template<class T> inline void Slat_wf<T>::getParmDepVal(Wavefunction_data * wfdata,
                            Sample_point * sample,
                            int e,
                            Array1 <doublevar> & oldval,
                            Wf_return & newval)
{

  if(parent->optimize_mo) {
    updateVal(wfdata, sample);
    getVal(wfdata, e, newval);
  }
  else if(parent->optimize_det) {
    assert(oldval.GetDim(0) >=2*ndet);
    if(nfunc_ != 1) error("Don't support several functions and ndet yet");
    doublevar tempval=0;
    int count=0;
    for(int det=0; det < ndet; det++) {
      tempval+=parent->detwt(det)*oldval(count)*oldval(count+1);
      count+=2;
    }
    newval.phase(0,0)=.5*pi*(1-sign(tempval));// pi if the function is negative

    if(fabs(tempval) > 0)
      newval.amp(0,0)=log(fabs(tempval));
    else
      newval.amp(0,0)=-1e3;
  }
  else { 
    assert(oldval.GetDim(0) >=2);
    assert(newval.amp.GetDim(1) >= 1);
    assert(newval.amp.GetDim(0) >= nfunc_);
    int counter=0;
    newval.amp(0,0)=oldval(counter++);
    newval.phase(0,0)=oldval(counter++);
    
  }
}


//-----------------------------------------------------------------------


template <> inline int Slat_wf<dcomplex>::getParmDeriv(Wavefunction_data *  wfdata, 
			  Sample_point * sample ,
			  Parm_deriv_return & derivatives){
  error("parmderiv not supported for complex orbitals yet");
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
  
  int nparms_full=parent->nparms();
  int nparms_start=derivatives.nparms_start;
  int nparms_end=derivatives.nparms_end;
  int nparms=nparms_end-nparms_start;
  
  derivatives.gradient.Resize(nparms);
  derivatives.hessian.Resize(nparms, nparms);
  derivatives.lapderiv.Resize(nparms);
 /* 
  if(parent->optimize_mo) {
    //get values of determinats with 1 and 2 rows differentiated
    Array3 < Array2 <T> > detGrad(nfunc_,ndet,2);
    Array2 < Array1 <doublevar> > BasisVal;
    Array3 < Array2 <doublevar> > doubletrace(nfunc_,ndet,2);
    
    int totmo=parent->molecorb->getNmo();
    int totcoeff=parent->molecorb->nMoCoeff();
    int nmocoeff=totcoeff/totmo;
    assert(nparms_full==parent->orbitals_for_optimize_mo.GetSize()*nmocoeff);
    int etotal=nelectrons(0)+nelectrons(1);

    if(nparms%nmocoeff)
      error("The number of MAXNPARMS_AT_ONCE != n*nmocoeff/orbital ");

    if(ndet>1){
      cout<<"WARNING!!! The analytic derivatives of orbital coeficients are calculated without normalization condition."<<endl;
      cout<<"This might be an issue for overal normalization if you have more than one determinant in your w.f."<<endl;
      cout<<"Please contact M.B. about details/fix"<<endl;
    }

    int totmo_start=nparms_start/nmocoeff;
    int totmo_end=nparms_end/nmocoeff;
    int totmo_diff=totmo_end-totmo_start;


    Array3 < Array1 <int> > which_orbitals(nfunc_,ndet,2), which_occupation(nfunc_,ndet,2);
    Array3 <int> dfuncdet(nfunc_,ndet,2);
    for(int f=0; f< nfunc_; f++)  
      for(int det=0; det< ndet; det++)
        for(int s=0; s< 2; s++) {
          dfuncdet(f,det,s)=0;
          vector <int> orbitals_tmp;
          vector <int> occupation_tmp;
          for(int i=totmo_start; i<totmo_end; i++) //because only these will be used
            for(int ii = 0; ii < nelectrons(s); ii++) 
              if (parent->occupation(f,det,s)(ii)==parent->orbitals_for_optimize_mo(i)){
                dfuncdet(f,det,s)=1;
                orbitals_tmp.push_back(i-totmo_start); //so this one starts from 0
                occupation_tmp.push_back(ii);
              }
          which_orbitals(f,det,s).Resize(orbitals_tmp.size());
          which_occupation(f,det,s).Resize(orbitals_tmp.size());
          for(int k=0;k<which_orbitals(f,det,s).GetSize();k++){
            which_orbitals(f,det,s)(k)=orbitals_tmp[k];
            which_occupation(f,det,s)(k)=occupation_tmp[k];
            //cout <<det<<"  "<<s<<"  "<<k<<" which_orbitals(f,det,s)(k) "<<which_orbitals(f,det,s)(k)
            //	 <<" which_occupation(f,det,s)(k) " <<which_occupation(f,det,s)(k)<<endl;
          }
        }
     
    //cout <<"setting up temp arrays"<<endl;
    for(int f=0; f< nfunc_; f++)  
      for(int det=0; det< ndet; det++)  
        for(int s=0;s<2;s++){
          detGrad(f,det,s).Resize(totmo_diff,nmocoeff);
          detGrad(f,det,s)=0.0;
        }
         
    
    //cout <<"getting BasisVal for all electrons"<<endl;
    BasisVal.Resize(2,etotal);
    for(int e=0; e< etotal; e++) {
      int s=parent->spin(e);
      int ee=parent->rede(e);
      parent->molecorb->getBasisVal(sample, e, BasisVal(s,ee));
      assert(BasisVal(s)(ee).GetSize()==nmocoeff);
    }

    //cout <<"getting gradients (single trace)"<<endl;
    for(int f=0; f< nfunc_; f++)  
      for(int det=0; det< ndet; det++) 
        for(int s=0; s< 2; s++)
          if(dfuncdet(f,det,s))
            for(int orb=0;orb<which_orbitals(f,det,s).GetSize();orb++){
              //cout <<" orb "<<orb<<endl;
              for(int j=0;j<nmocoeff;j++)
                for(int ee = 0; ee < nelectrons(s); ee++) {
                  detGrad(f,det,s)(which_orbitals(f,det,s)(orb),j)+=
                    inverse(f,det,s)(ee,which_occupation(f,det,s)(orb))*BasisVal(s,ee)(j);
                }   
            }
    
    //get the wf's val
    T sum=0;
    for(int det=0; det < ndet; det++) {
      sum+=parent->detwt(det)*(detVal(0,det,0)*detVal(0,det,1)).val();
    }
    
    //calculate the actual gradient and hessian
    //cout <<"calculating the actual gradient and hessian"<<endl;
    derivatives.gradient=0.0;
    for(int orb=0;orb<totmo_diff;orb++)
      for(int coef=0;coef<nmocoeff;coef++){
        int i=orb*nmocoeff+coef;
        derivatives.gradient(i)=0.0;
        for(int f=0; f< nfunc_; f++)  
          for(int det=0; det< ndet; det++){ 
            T temp=0.0;
            for(int s=0; s< 2; s++)
              //if(dfuncdet(f,det,s))
              temp+=detGrad(f,det,s)(orb,coef);
            derivatives.gradient(i)+=(parent->detwt(det)*detVal(0,det,0)*detVal(0,det,1)*temp).val();
          }
        derivatives.gradient(i)/=sum;
      }
    
    derivatives.hessian=0.0;
    if(derivatives.need_hessian){
      for(int orbi=0;orbi<totmo_diff;orbi++)
        for(int coefi=0;coefi<nmocoeff;coefi++){
          int i=orbi*nmocoeff+coefi;
          for(int orbj=0;orbj<totmo_diff;orbj++)
            for(int coefj=0;coefj<nmocoeff;coefj++){
              int j=orbj*nmocoeff+coefj;
              derivatives.hessian(i,j)=0.0;
              for(int f=0; f< nfunc_; f++)  
                for(int det=0; det< ndet; det++){
                  doublevar temp=0;
                  for(int s=0; s< 2; s++)
                    if(dfuncdet(f,det,s) && orbi!=orbj){
                      temp+=detGrad(f,det,s)(orbi,coefi)*detGrad(f,det,s)(orbj,coefj)
                        -detGrad(f,det,s)(orbi,coefj)*detGrad(f,det,s)(orbj,coefi);
                    }
                  temp+=detGrad(0,det,0)(orbi,coefi)*detGrad(0,det,1)(orbj,coefj)
                    +detGrad(0,det,0)(orbj,coefj)*detGrad(0,det,1)(orbi,coefi);
                  derivatives.hessian(i,j)+=(parent->detwt(det)*detVal(0,det,0)*detVal(0,det,1)).val()*temp;
                }
              //cout <<"derivatives.hessian(i,j) "<<derivatives.hessian(i,j)<<endl;
              derivatives.hessian(i,j)/=sum;
              derivatives.hessian(j,i)=derivatives.hessian(i,j);
            }
        }
      return 1;
    }
  }
*/
  if(parent->optimize_mo) { 
    error("don't support optimizing mo yet");
  }
  else if(parent->optimize_det) {
    log_value<doublevar> detsum=0;
    Array1 <log_value<doublevar> > detvals(ndet);
    for(int det=0; det < ndet; det++) {
      detvals(det)=parent->detwt(det)*detVal(0,det,0)*detVal(0,det,1);
    }
    detsum=sum(detvals);
    detsum.logval*=-1;
    derivatives.gradient=0.0;
    int counter=0;
    for(int csf=0; csf < parent->ncsf; csf++) { 
      if(csf > nparms_start && csf <= nparms_end ){
        for(int j=1;j<parent->CSF(csf).GetDim(0);j++){
          derivatives.gradient(csf-1-nparms_start)+=
            parent->CSF(csf)(j)*(detVal(0,counter,0)*detVal(0,counter,1)).val();
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
    return 1;

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
inline doublevar real(doublevar & a) { return a; } 

template <class T>inline void Slat_wf<T>::updateInverse(Slat_wf_data * dataptr, int e) { 
  int maxmatsize=max(nelectrons(0),nelectrons(1));
  Array1 <T> modet(maxmatsize);
  int s=spin(e);
  for(int f=0; f< nfunc_; f++)  {
    for(int det=0; det< ndet; det++)  {
      //fill the molecular orbitals for this
      //determinant
      if(real(detVal(f,det,s).logval) < -1e200) { 
        Array2 <T> allmos(nelectrons(s), nelectrons(s));
        for(int e=0; e< nelectrons(s); e++) {
          int curre=s*nelectrons(0)+e;
          for(int i=0; i< nelectrons(s); i++) {
            allmos(e,i)=moVal(0,curre, dataptr->occupation(f,det,s)(i));
          }
        }


        detVal(f,det,s)=
          TransposeInverseMatrix(allmos,inverse(f,det,s), nelectrons(s));

      }
      else { 
        for(int i = 0; i < nelectrons(s); i++) {
          modet(i)=moVal(0,e,dataptr->occupation(f,det,s)(i));
        }
        T ratio=1./InverseUpdateColumn(inverse(f,det,s),
            modet, dataptr->rede(e),
            nelectrons(s));

        detVal(f,det, s)=ratio*detVal(f,det, s);
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
      if(abs(detVal(f,det,s).logval) < -1e200) return 0;
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
                                            modet, dataptr->rede(e),
                                            nelectrons(s));

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
  int s=dataptr->spin(e);

  //update all the mo's that we will be using.
  molecorb->updateVal(sample, e, s,
                              updatedMoVal);

  for(int i=0; i< updatedMoVal.GetDim(0); i++)
    moVal(0,e,i)=updatedMoVal(i,0);


  inverseStale=1;
  lastValUpdate=e;
  lastDetVal=detVal;
  if(!updateValNoInverse(dataptr, e)) { 
    inverseStale=0;
    updateInverse(dataptr,e);
  }


}

//------------------------------------------------------------------------


template <class T>inline void Slat_wf<T>::getVal(Wavefunction_data * wfdata, int e,
                     Wf_return & val)
{
  //Array1 <doublevar> si(nfunc_, 0.0);
  Array2 <log_value<T> > vals(nfunc_,1,T(0.0));

  assert(val.amp.GetDim(0) >=nfunc_);
  assert(val.amp.GetDim(1) >= 1);
 

  if(staticSample==1 && parent->optimize_mo==0 && parent->optimize_det==0) {
    //lap.phase=0;
    for(int f=0; f< nfunc_; f++) {
      for(int i=0; i< 1; i++) 
        vals(f,i)=saved_laplacian(e,f,i);

      //si(f)=saved_laplacian(e,f,5);
    }
  }
  else {
    Slat_wf_data * dataptr;
    recast(wfdata, dataptr);

    int s=dataptr->spin(e);
    int opp=dataptr->opspin(e);

    for(int f=0; f< nfunc_; f++) {
      Array1 <log_value<T> > detvals(ndet);
      for(int det=0; det < ndet; det++) {
        detvals(det) = T(dataptr->detwt(det))*detVal(f,det,s)*detVal(f,det,opp);
      }
      log_value<T> totval=sum(detvals);
      //vals(f,0)=totval.logval;
      //si(f)=totval.sign;
      vals(f,0)=totval;
    }
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

template <class T> inline void Slat_wf<T>::getDensity(Wavefunction_data * wfdata, int e,
                         Array2 <doublevar> & dens)
{

  assert(dens.GetDim(0) >= nfunc_);
  Slat_wf_data * dataptr;
  recast(wfdata, dataptr);

  int s=dataptr->spin(e);

  if(ndet > 1)
  {
    error("Haven't done density for several determinants yet.");
  }
  dens=0;
  int det=0;
  for(int f=0; f< nfunc_; f++)
  {
    for(int j=0; j< nelectrons(s); j++)
    {
      dens(f,0)+=abs(moVal(0 , e, dataptr->occupation(f,det,s)(j) )
                 *moVal(0,e,dataptr->occupation(f,det,s)(j)));
    }
  }

}

//----------------------------------------------------------------------------


template <class T> inline void Slat_wf<T>::calcLap(Slat_wf_data * dataptr, Sample_point * sample)
{
  //cout << "calcLap " << endl;
  inverseStale=0;
  for(int e=0; e< nelectrons(0)+nelectrons(1); e++)  {
    int s=dataptr->spin(e);
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
        //if(f==0 && det==0 && s==0) matout << "determinant " << detVal(f,det,s) 
         //   << " should be " << modet(0,0)*modet(1,1)-modet(0,1)*modet(1,0) << endl;
      }
    }
  }
  //cout << "done " << endl;
}

//------------------------------------------------------------------------


/*!
*/

template <class T> void Slat_wf<T>::getLap(Wavefunction_data * wfdata,
                     int e, Wf_return & lap)
{

  //Array1 <doublevar> si(nfunc_, 0.0);
  //Array2 <doublevar> vals(nfunc_,5,0.0);
  Array2 <log_value <T> > vals(nfunc_,5);

  if(staticSample==1 && parent->optimize_mo==0 && parent->optimize_det==0) {
    //lap.phase=0;
    for(int f=0; f< nfunc_; f++) {
      for(int i=0; i< 5; i++) 
        vals(f,i)=saved_laplacian(e,f,i);

      //si(f)=saved_laplacian(e,f,5);
    }
  }
  else {
    Slat_wf_data * dataptr;
    recast(wfdata, dataptr);

    int s=dataptr->spin(e);
    int opp=dataptr->opspin(e);

    for(int f=0; f< nfunc_; f++) {
      Array1 <log_value <T> > detvals(ndet);
      for(int det=0; det < ndet; det++) {
        detvals(det) = T(dataptr->detwt(det))*detVal(f,det,s)*detVal(f,det,opp);
      }
      log_value<T> totval=sum(detvals);
      //vals(f,0)=totval.logval;
      //si(f)=totval.sign;
      vals(f,0)=totval;

      log_value<T> invtotval=totval;
      invtotval.logval*=-1;
      Array1 <log_value <T> > detgrads(ndet);
      for(int i=1; i< 5; i++) {
        //lap.amp(f,i)=0;
        for(int det=0; det < ndet; det++) {
          T temp=0;
          for(int j=0; j<nelectrons(s); j++) {
            temp+=moVal(i , e, dataptr->occupation(f,det,s)(j) )
                 *inverse(f,det,s)(dataptr->rede(e), j);
          }
          detgrads(det)=temp*dataptr->detwt(det);//*detVal(f,det,s)*detVal(f,det,opp)*invtotval;
          detgrads(det)*=detVal(f,det,s);
          detgrads(det)*=detVal(f,det,opp);
          detgrads(det)*=invtotval;
        }

        //log_real_value totgrad=sum(detgrads);

        //vals(f,i)=totgrad.val();
        vals(f,i)=sum(detgrads);
      }
      
    }
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

  int s=dataptr->spin(e);
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
    for(int det=0; det< ndet; det++)  {
      for(int i = 0; i < nelectrons(s); i++) {
        modet(i)=movals(s)(parent->occupation(f,det,s)(i),0);
      }
      T ratio=1./InverseGetNewRatio(inverse(f,det,s),
          modet, parent->rede(e),
          nelectrons(s));
      new_detVals(det)=T(parent->detwt(det))*ratio*detVal(f,det, s)*detVal(f,det,opps);
    }
    log_value<T> totval=sum(new_detVals);
    vals(f,0)=totval;
    wf(e).setVals(vals);

  }

  

}


#endif //SLAT_WF_H_INCLUDED
//--------------------------------------------------------------------------
