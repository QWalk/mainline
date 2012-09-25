/*
 
Based of Slat_wf_calc.cpp Copyright (C) 2007 Lucas K. Wagner
Modified for "Fast" updates of multideterminant wavefunctions
using iterative updates using only ground state determinant by
Paul Kent / CNMS / ORNL 2009.

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

// Debugging defines to enable iterative updates for specific methods only

#define RECDET
#ifdef RECDET
#define RECDETUI
#define RECDETUVNI
#define RECDETCL
#define RECDETGL
#endif

// Debugging:
//#define RECDETDBG
//#define RECDETUI
//#define RECDETUVNI
//#define RECDETCL
//#define RECDETGL

#include "Qmc_std.h"
#include "FSlat_wf.h"
#include "MatrixAlgebra.h"
#include "Sample_point.h"
#include "Slat_wf_data.h"

//----------------------------------------------------------------------


void FSlat_wf::generateStorage(Wavefunction_storage * & wfstore)
{
  wfstore=new FSlat_wf_storage;
  FSlat_wf_storage * store;
  recast(wfstore, store);
  store->moVal_temp.Resize (5,   nmo);
  
  // Added by Matous
  store->moVal_temp_2.Resize (5,   nmo);

  store->detVal_temp.Resize(nfunc_, ndet, 2);
#ifndef RECDET
  store->inverse_temp.Resize(nfunc_, ndet, 2); 
#else
  store->inverse_temp.Resize(nfunc_, 1, 2); 
#endif

#ifndef RECDET
  for(int i=0; i< nfunc_; i++)
  {
    for(int det=0; det < ndet; det++)
    {
      for(int s=0; s<2; s++)
      {
        store->inverse_temp(i,det,s).Resize(nelectrons(s), nelectrons(s));
        store->inverse_temp(i,det,s)=0;
        store->detVal_temp(i,det,s)=1;
      }
    }
  }
#else
  for(int i=0; i< nfunc_; i++)
  {
    for(int s=0; s<2; s++)
      {
	for(int det=0; det < ndet; det++)
	{
	  store->detVal_temp(i,det,s)=1;
	}
	store->inverse_temp(i,0,s).Resize(nelectrons(s), nelectrons(s));
	store->inverse_temp(i,0,s)=0;
      }
  }
#endif

}


//----------------------------------------------------------------------

void FSlat_wf::init(Wavefunction_data * wfdata)
{

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
#ifdef RECDET
  inverse.Resize(nfunc_, ndet, 2);
#else
  inverse.Resize(nfunc_, 1, 2);
#endif

  for(int i=0; i< nfunc_; i++) {
#ifdef RECDET    
    int tmpndet=1;
#else
    int tmpndet=ndet
#endif
    for(int det=0; det < tmpndet; det++) {
      for(int s=0; s<2; s++) {
        inverse(i,det,s).Resize(nelectrons(s), nelectrons(s));
        inverse(i,det,s)=0;
        for(int e=0; e< nelectrons(s); e++) {
          inverse(i,det,s)(e,e)=1;
          inverse(i,det,s)(e,e)=1;
        }
      }
    }
  }

  for(int i=0; i< nfunc_; i++) {
    for(int det=0; det < ndet; det++) {
      for(int s=0; s<2; s++) {
        detVal(i,det,s)=1;
      }
    }
  }

  invgammastore.Resize(dataptr->max_occupation_changes);
  ustore.Resize(dataptr->max_occupation_changes);
  veestore.Resize(dataptr->max_occupation_changes);
  int maxe;
  if (nelectrons(0)>nelectrons(1)) { maxe=nelectrons(0); } else { maxe=nelectrons(1); }
  for (int k=0; k<dataptr->max_occupation_changes; k++) {
    ustore(k).Resize(maxe);
    veestore(k).Resize(maxe);
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
void FSlat_wf::notify(change_type change, int num)
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

void FSlat_wf::save_for_static() {
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

void FSlat_wf::saveUpdate(Sample_point * sample, int e,
                         Wavefunction_storage * wfstore)
{

  if(staticSample==0) {
    FSlat_wf_storage * store;
    recast(wfstore, store);

    //presumably, if we care enough to save the update, we care enough
    //to have the inverse up to date
    if(inverseStale) { 
      detVal=lastDetVal;
      updateInverse(parent, lastValUpdate);
      inverseStale=0;
    }
    int s=spin(e);

    for(int f=0; f< nfunc_; f++)
    {
#ifdef RECDET
      int tmpndet=1;
#else
      int tmpndet=ndet;
#endif
      for(int det=0; det<tmpndet; det++) {
        store->inverse_temp(f,det,s)=inverse(f,det,s);
      }
      for(int det=0; det<ndet; det++) {
        store->detVal_temp(f,det,s)=detVal(f,det,s);
      }
    }


    for(int d=0; d< 5; d++)
    {
      for(int i=0; i< moVal.GetDim(2); i++)
      {
        store->moVal_temp(d,i)=moVal(d,e,i);
      }
    }
  }

}

//----------------------------------------------------------------------

void FSlat_wf::restoreUpdate(Sample_point * sample, int e,
                            Wavefunction_storage * wfstore)
{

  if(staticSample==0) {
    FSlat_wf_storage * store;
    recast(wfstore, store);
    int s=spin(e);
    inverseStale=0;
    
    for(int j=0; j<5; j++)
    {
      for(int i=0; i<moVal.GetDim(2); i++)
      {
        moVal(j,e,i)=store->moVal_temp(j,i);
      }
    }

    for(int f=0; f< nfunc_; f++)
    {
#ifdef RECDET
      int tmpndet=1;
#else
      int tmpndet=ndet;
#endif
      for(int det=0; det < tmpndet; det++)
      {
        inverse(f,det,s)=store->inverse_temp(f,det,s);
      }
      for(int det=0; det < ndet; det++)
      {
        detVal(f,det,s)=store->detVal_temp(f,det,s);
      }
    }

    electronIsStaleVal(e)=0;
    electronIsStaleLap(e)=0;
  }

}

//----------------------------------------------------------------------

// Added by Matous
void FSlat_wf::saveUpdate(Sample_point * sample, int e1, int e2,
                         Wavefunction_storage * wfstore)
{

  if(staticSample==0) {
    FSlat_wf_storage * store;
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
#ifdef RECDET
      int tmpndet=1;
#else
      int tmpndet=ndet;
#endif
      for(int det=0; det<tmpndet; det++) {
        if ( s1 == s2 ) {
		store->inverse_temp(f,det,s1)=inverse(f,det,s1);
	}
	else {
		store->inverse_temp(f,det,s1)=inverse(f,det,s1);
		store->inverse_temp(f,det,s2)=inverse(f,det,s2);
	} 
      }
      for(int det=0; det<ndet; det++) {
        if ( s1 == s2 ) {
		store->detVal_temp(f,det,s1)=detVal(f,det,s1);
	}
	else {
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
void FSlat_wf::restoreUpdate(Sample_point * sample, int e1, int e2,
                            Wavefunction_storage * wfstore)
{

  if(staticSample==0) {
    FSlat_wf_storage * store;
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
#ifdef RECDET
      int tmpndet=1;
#else
      int tmpndet=ndet;
#endif
      for(int det=0; det < tmpndet; det++)
      {
         if ( s1 == s2 )
         {
	   inverse(f,det,s1)=store->inverse_temp(f,det,s1);
         }
         else
         {
	   inverse(f,det,s1)=store->inverse_temp(f,det,s1);
	   inverse(f,det,s2)=store->inverse_temp(f,det,s2);
	 }
      }
      for(int det=0; det < ndet; det++)
      {
         if ( s1 == s2 )
         {
	    detVal(f,det,s1)=store->detVal_temp(f,det,s1);
	 }
	 else
	 {
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

void FSlat_wf::updateVal(Wavefunction_data * wfdata,
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

void FSlat_wf::updateLap( Wavefunction_data * wfdata,
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


void FSlat_wf::storeParmIndVal(Wavefunction_data * wfdata, Sample_point * sample,
                              int e, Array1 <doublevar> & vals )
{

  if(parent->optimize_mo) {

  }
  else if(parent->optimize_det) {
    assert(vals.GetDim(0) >= 2*parent->detwt.GetDim(0));
    updateVal(wfdata, sample);
    int count=0;
    for(int det=0; det < ndet; det++) {
      vals(count++)=detVal(0,det,0);
      vals(count++)=detVal(0,det,1);
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
}

//----------------------------------------------------------------------

void FSlat_wf::getParmDepVal(Wavefunction_data * wfdata,
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
      tempval+=parent->detwt(det).val()*oldval(count)*oldval(count+1);
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

int FSlat_wf::getParmDeriv(Wavefunction_data *  wfdata, 
			  Sample_point * sample ,
			  Parm_deriv_return & derivatives){
  error("Need to update Fslat_wf parmderiv"); 
 /* 
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
  
  if(parent->optimize_mo) {
    //get values of determinats with 1 and 2 rows differentiated
    Array3 < Array2 <doublevar> > detGrad(nfunc_,ndet,2);
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
		  error("LAZY PROGRAMMER in FSlat_wf_calc.cpp:getParmDeriv  not implemented for fast updates");
		  detGrad(f,det,s)(which_orbitals(f,det,s)(orb),j)+=inverse(f,det,s)(ee,which_occupation(f,det,s)(orb))*BasisVal(s,ee)(j);
		}   
	    }
    
    //get the wf's val
    doublevar sum=0;
    for(int det=0; det < ndet; det++) {
      sum+=parent->detwt(det)*detVal(0,det,0)*detVal(0,det,1);
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
	    doublevar temp=0;
	    for(int s=0; s< 2; s++)
	      //if(dfuncdet(f,det,s))
	      temp+=detGrad(f,det,s)(orb,coef);
	    derivatives.gradient(i)+=parent->detwt(det)*detVal(0,det,0)*detVal(0,det,1)*temp;
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
		  derivatives.hessian(i,j)+=parent->detwt(det)*detVal(0,det,0)*detVal(0,det,1)*temp;
		}
	      //cout <<"derivatives.hessian(i,j) "<<derivatives.hessian(i,j)<<endl;
	      derivatives.hessian(i,j)/=sum;
	      derivatives.hessian(j,i)=derivatives.hessian(i,j);
	    }
	}
      return 1;
    }
  }
  else if(parent->optimize_det) {
    doublevar sum=0;
    for(int det=0; det < ndet; det++) {
      sum+=parent->detwt(det)*detVal(0,det,0)*detVal(0,det,1);
    }
    if(parent->use_csf){
      if(parent->all_weights)
        assert(nparms<parent->ncsf+1);
      else
        assert(nparms<parent->ncsf);
      int counter=0;
      derivatives.gradient=0.0;
      for(int csf=0; csf< parent->ncsf; csf++) {
        if(parent->all_weights){
          if(csf >= nparms_start && csf <= nparms_end ){
            for(int j=1;j<parent->CSF(csf).GetDim(0);j++){
              derivatives.gradient(csf-nparms_start)+=parent->CSF(csf)(j)*detVal(0,counter,0)*detVal(0,counter,1);
              counter++;
            }
          }
          else{
            counter+=parent->CSF(csf).GetDim(0)-1;
          }
        }
        else{
          if(csf > nparms_start && csf <= nparms_end ){
            for(int j=1;j<parent->CSF(csf).GetDim(0);j++){
              derivatives.gradient(csf-1-nparms_start)+=parent->CSF(csf)(j)*detVal(0,counter,0)*detVal(0,counter,1);
              counter++;
            }
          }
          else{
            counter+=parent->CSF(csf).GetDim(0)-1;
          }
        }
      }
      for(int csf=0; csf< nparms; csf++) {
        derivatives.gradient(csf)/=sum; 
      }
      derivatives.hessian=0;
      return 1;
    }
    else{
      for(int det=0; det < ndet; det++) {
        if(parent->all_weights){
          if(det >= nparms_start && det <= nparms_end )
            derivatives.gradient(det-nparms_start)=detVal(0,det,0)*detVal(0,det,1);
        }
        else{
          if(det > nparms_start && det <= nparms_end )
            derivatives.gradient(det-1-nparms_start)=detVal(0,det,0)*detVal(0,det,1);
        }
      }

      for(int det=0; det < nparms; det++) {
        derivatives.gradient(det)/=sum; 
      }
      derivatives.hessian=0;
      return 1;
    }
  }
  else { 
    derivatives.gradient=0;
    derivatives.hessian=0;
    return 1;
  }
  return 0;
  */
}


//------------------------------------------------------------------------


void FSlat_wf::calcVal(Slat_wf_data * dataptr, Sample_point * sample)
{
  //Hmm, I don't completely understand why, but something is not 
  //completely clean, so we can't just cycle through and updateVal
  //This is actually probably the best way to do it anyway, since it 
  //should in theory be faster, and it gives us a clean start
  calcLap(dataptr, sample);

}

//------------------------------------------------------------------------


void FSlat_wf::updateInverse(Slat_wf_data * dataptr, int e) { 
  int maxmatsize=max(nelectrons(0),nelectrons(1));
  Array1 <doublevar> modet(maxmatsize);
  int s=spin(e);
#ifdef RECDETUI
  for(int f=0; f< nfunc_; f++)  { 

    // If initial determinant is zero attempt rebuild
    if ( fabs(fabs(detVal(f,0,s))) == 0 ) {
        Array2 <doublevar> allmos(nelectrons(s), nelectrons(s));
        for(int e=0; e< nelectrons(s); e++) {
          int curre=s*nelectrons(0)+e;
          for(int i=0; i< nelectrons(s); i++) {
            allmos(e,i)=moVal(0,curre, dataptr->occupation(f,0,s)(i));
          }
        }
        detVal(f,0,s)=
          TransposeInverseMatrix(allmos,inverse(f,0,s), nelectrons(s)).val();
    }


    if (fabs(detVal(f,0,s))>0) {
      // First compute ratio for zeroth determinant, updating the inverse
      for(int i = 0; i < nelectrons(s); i++) {
	modet(i)=moVal(0,e,dataptr->occupation(f,0,s)(i));
      }
      doublevar zero_ratio=1./InverseUpdateColumn(inverse(f,0,s),
					     modet, dataptr->rede(e),
					     nelectrons(s));
      for(int idx=1; idx< ndet; idx++)  {
	int recdet=dataptr->evaluation_order(f,idx,s);
	int nchanges=dataptr->occupation_nchanges(f,recdet,s);
	doublevar newval;
	if (nchanges==0) 
	{
           newval=zero_ratio*detVal(f,recdet,s);
	}
	else
	{
	  doublevar calcq,calcqtot;
	  calcqtot=1.0;
	  int ichstart=dataptr->occupation_first_diff_change_from_last_det(f,recdet,s);
	  for (int ichange=0; ichange<ichstart; ichange++)
	  {
	    calcqtot*=invgammastore(ichange);
	  }
	  for (int ichange=ichstart; ichange<nchanges; ichange++)
	  {
	    // Factor the ichangeth change
	    int erow=dataptr->occupation_changes(f,recdet,s)(ichange);
	    for (int j=0; j<nelectrons(s); j++)
	    {
	      ustore(ichange)(j)=inverse(f,0,s)(j,erow);
	    }
	    for (int jchange=0;jchange<ichange;jchange++) 
	    {
	      doublevar prod=dot(ustore(ichange),veestore(jchange))/invgammastore(jchange);
	      for (int j=0; j<nelectrons(s); j++)
	      {
		ustore(ichange)(j)=ustore(ichange)(j)-prod*ustore(jchange)(j);
	      }
	    }
		      
	    for(int j=0; j<nelectrons(s); j++)
	    {
	      int jr=j+s*nelectrons(0);
	      veestore(ichange)(j)=moVal(0,jr,dataptr->occupation(f,recdet,s)(erow))-moVal(0,jr,dataptr->occupation(f,0,s)(erow));
	    }
	    calcq=1.0+dot(veestore(ichange),ustore(ichange));
	    calcqtot*=calcq;
	    invgammastore(ichange)=calcq;
	  } // end for ichange
	  newval=calcqtot*zero_ratio*detVal(f,0,s);
	}

#ifdef RECDETDBG
	      Array2 <doublevar> allmos(nelectrons(s), nelectrons(s));
	      for(int e=0; e< nelectrons(s); e++) {
		int curre=s*nelectrons(0)+e;
		for(int i=0; i< nelectrons(s); i++) {
		  allmos(e,i)=moVal(0,curre, dataptr->occupation(f,recdet,s)(i));
		}
	      }
	      doublevar tdetVal=TransposeInverseMatrix(allmos,inverse(f,recdet,s), nelectrons(s));
	      if (abs(tdetVal-curratio*detVal(f,recdet,s))>1e-5) {
		cout << "FSlat_wf_calc.cpp:UpdateInverse RECDETDBG accuracy problem" << f << " " << recdet << " " << s << endl;
		cout << "Reference="<< tdetVal << " actual="<< detVal(f,recdet,s) << endl;
		error ("FSlat_wf_calc.cpp:UpdateInverse RECDETDBG accuracy problem");
	      }

#endif
	if (recdet>0)
	{
	  detVal(f,recdet,s)=newval;
	} 
      }
      detVal(f,0,s)=zero_ratio*detVal(f,0,s);
    }
    else
    {
      error("FSlat_wf_calc.cpp:UpdateInverse not implemented for ZERO valued ground state (1st) determinant");
    }
  }
#else

  for(int f=0; f< nfunc_; f++)  {
    for(int det=0; det< ndet; det++)  {
      //fill the molecular orbitals for this
      //determinant
      if(fabs(detVal(f,det,s)) > 0) { 
        for(int i = 0; i < nelectrons(s); i++) {
          modet(i)=moVal(0,e,dataptr->occupation(f,det,s)(i));
        }
        
        
        doublevar ratio=1./InverseUpdateColumn(inverse(f,det,s),
                                     modet, dataptr->rede(e),
                                     nelectrons(s));
        detVal(f,det, s)=ratio*detVal(f,det, s);
      }
      else { 
        
        Array2 <doublevar> allmos(nelectrons(s), nelectrons(s));
        for(int e=0; e< nelectrons(s); e++) {
          int curre=s*nelectrons(0)+e;
          for(int i=0; i< nelectrons(s); i++) {
            allmos(e,i)=moVal(0,curre, dataptr->occupation(f,det,s)(i));
          }
        }
        
        
        detVal(f,det,s)=
          TransposeInverseMatrix(allmos,inverse(f,det,s), nelectrons(s));
      }
    }
  }
#endif  
}

//------------------------------------------------------------------------

int FSlat_wf::updateValNoInverse(Slat_wf_data * dataptr, int e) { 
  int maxmatsize=max(nelectrons(0),nelectrons(1));
  Array1 <doublevar> modet(maxmatsize);
  int s=spin(e);
  for(int f=0; f< nfunc_; f++)  {
    for(int det=0; det< ndet; det++)  {
      //fill the molecular orbitals for this
      //determinant
      if(!(fabs(detVal(f,det,s)) > 0)) return 0;
    }
  }
        
  
  for(int f=0; f< nfunc_; f++)
  {
#ifdef RECDETUVNI
      // First compute ratio for zeroth determinant, updating the inverse
#ifdef RECDETDBG
    for(int det=0; det< ndet; det++)  {
      Array2 <doublevar> tinverse=inverse(f,det,s);

      //fill the molecular orbitals for this
      //determinant
      for(int i = 0; i < nelectrons(s); i++) {
        modet(i)=moVal(0,e,dataptr->occupation(f,det,s)(i));
      }
      
      
      doublevar ratio=1./InverseGetNewRatio(tinverse,
                                            modet, dataptr->rede(e),
                                            nelectrons(s));
      cout << "UVNI REF " << s << " " << det << " IVAL=" << detVal(f,det,s) << " VAL/IVAL=" << ratio << " VAL=" << ratio*detVal(f,det,s) << endl;
    }
#endif
    for(int i = 0; i < nelectrons(s); i++) {
	modet(i)=moVal(0,e,dataptr->occupation(f,0,s)(i));
    }
    Array2 <doublevar> zeroinverse=inverse(f,0,s);
    doublevar zero_ratio=1./InverseUpdateColumn(zeroinverse,
					     modet, dataptr->rede(e),
					     nelectrons(s));
    for(int idx=1; idx< ndet; idx++)
    {
      int recdet=dataptr->evaluation_order(f,idx,s);
      int nchanges=dataptr->occupation_nchanges(f,recdet,s);
      doublevar newval;
      if (nchanges==0) 
      {
	newval=zero_ratio*detVal(f,recdet,s);
      }
      else
      {
	doublevar calcq,calcqtot;
	calcqtot=1.0;
	int ichstart=dataptr->occupation_first_diff_change_from_last_det(f,recdet,s);
	for (int ichange=0; ichange<ichstart; ichange++)
	{
	  calcqtot*=invgammastore(ichange);
	}
	for (int ichange=ichstart; ichange<nchanges; ichange++)
	{
	  // Factor the ichangeth change
	  int erow=dataptr->occupation_changes(f,recdet,s)(ichange);
	  for (int j=0; j<nelectrons(s); j++)
	  {
	    ustore(ichange)(j)=zeroinverse(j,erow);
	  }
	  for (int jchange=0;jchange<ichange;jchange++) 
	  {
	    doublevar prod=dot(ustore(ichange),veestore(jchange))/invgammastore(jchange);
	    for (int j=0; j<nelectrons(s); j++)
	    {
	      ustore(ichange)(j)=ustore(ichange)(j)-prod*ustore(jchange)(j);
	    }
	  }
	  
	  for(int j=0; j<nelectrons(s); j++)
	  {
	    int jr=j+s*nelectrons(0);
	    veestore(ichange)(j)=moVal(0,jr,dataptr->occupation(f,recdet,s)(erow))-moVal(0,jr,dataptr->occupation(f,0,s)(erow));
	  }
	  calcq=1.0+dot(veestore(ichange),ustore(ichange));
	  calcqtot*=calcq;
	  invgammastore(ichange)=calcq;
	} // end for ichange
	newval=calcqtot*zero_ratio*detVal(f,0,s);
      }

      if (recdet>0) {
	detVal(f,recdet, s)=newval;
      }
      
    }
    detVal(f,0, s)=zero_ratio*detVal(f,0, s);	    
#else
    for(int det=0; det< ndet; det++)  {

      //fill the molecular orbitals for this
      //determinant
      for(int i = 0; i < nelectrons(s); i++) {
        modet(i)=moVal(0,e,dataptr->occupation(f,det,s)(i));
      }
      
      
      doublevar ratio=1./InverseGetNewRatio(inverse(f,det,s),
                                            modet, dataptr->rede(e),
                                            nelectrons(s));
      detVal(f,det, s)=ratio*detVal(f,det, s);

    }
#endif
  }
  return 1;
}

//------------------------------------------------------------------------
/*!

*/
void FSlat_wf::updateVal( Slat_wf_data * dataptr, Sample_point * sample,int e) {

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
  dataptr->molecorb->updateVal(sample, e, s,
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


void FSlat_wf::getVal(Wavefunction_data * wfdata, int e,
                     Wf_return & val)
{
  Array1 <doublevar> si(nfunc_, 0.0);
  Array2 <doublevar> vals(nfunc_,1,0.0);

  assert(val.amp.GetDim(0) >=nfunc_);
  assert(val.amp.GetDim(1) >= 1);
 

  if(staticSample==1 && parent->optimize_mo==0 && parent->optimize_det==0) {
    //lap.phase=0;
    for(int f=0; f< nfunc_; f++) {
      for(int i=0; i< 1; i++) 
        vals(f,i)=saved_laplacian(e,f,i);

      si(f)=saved_laplacian(e,f,5);
    }
  }
  else {
    Slat_wf_data * dataptr;
    recast(wfdata, dataptr);

    int s=dataptr->spin(e);
    int opp=dataptr->opspin(e);

    for(int f=0; f< nfunc_; f++) {
      doublevar funcval=0;

      for(int det=0; det < ndet; det++) {
        funcval += dataptr->detwt(det).val()*detVal(f,det,s)*detVal(f,det,opp);
      }
      if(fabs(funcval) > 0) 
        vals(f,0)=log(fabs(funcval));
      else vals(f,0)=-1e3;
      si(f)=sign(funcval);
    }
  }

  val.setVals(vals, si);

}

//----------------------------------------------------------
void FSlat_wf::getSymmetricVal(Wavefunction_data * wfdata,
		     int e, Wf_return & val){
  val.phase(0, 0)=0;
  val.amp(0, 0)=0;
  val.cvals(0,0)=0;
} 

//----------------------------------------------------------------------

void FSlat_wf::getDensity(Wavefunction_data * wfdata, int e,
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
      dens(f,0)+=moVal(0 , e, dataptr->occupation(f,det,s)(j) )
                 *moVal(0,e,dataptr->occupation(f,det,s)(j));
    }
  }

}

//----------------------------------------------------------------------------


void FSlat_wf::calcLap(Slat_wf_data * dataptr, Sample_point * sample)
{

  inverseStale=0;
  for(int e=0; e< nelectrons(0)+nelectrons(1); e++)  {
    int s=dataptr->spin(e);
    sample->updateEIDist();

    //update all the mo's that we will be using, using the lists made in
    //Slat_wf_data(one for each spin).
    dataptr->molecorb->updateLap(sample, e, s,
                                updatedMoVal);
    for(int d=0; d< 5; d++)  {
      for(int i=0; i< updatedMoVal.GetDim(0); i++) {
        moVal(d,e,i)=updatedMoVal(i,d);
      }
    }
  }

  int maxmatsize=max(nelectrons(0),nelectrons(1));
  Array2 <doublevar> modet(maxmatsize, maxmatsize);

#ifndef RECDETCL
  for(int f=0; f< nfunc_; f++)   {
    for(int det=0; det < ndet; det++ ) {
      for(int s=0; s< 2; s++ ) {


        for(int e=0; e< nelectrons(s); e++) {
          int curre=s*nelectrons(0)+e;
          for(int i=0; i< nelectrons(s); i++) {
            modet(e,i)=moVal(0,curre, dataptr->occupation(f,det,s)(i));
          }
        }
        
        detVal(f,det,s)=
          TransposeInverseMatrix(modet,inverse(f,det,s), nelectrons(s));

      }
    }
  }
#else
  for(int f=0; f< nfunc_; f++)   {
    int det=0;
      for(int s=0; s< 2; s++ ) {
	// Fill the zeroth Slater determinant, setup inverse and value

	for(int e=0; e< nelectrons(s); e++) {
	  int curre=s*nelectrons(0)+e;
	  for(int i=0; i< nelectrons(s); i++) {
	    modet(e,i)=moVal(0,curre, dataptr->occupation(f,det,s)(i));
	  }
	}
	detVal(f,det,s)=
	  TransposeInverseMatrix(modet,inverse(f,det,s), nelectrons(s)).val();
	    
	for(int idx=1; idx< ndet; idx++)
        {
	  int recdet=dataptr->evaluation_order(f,idx,s);
	  int nchanges=dataptr->occupation_nchanges(f,recdet,s);
	  doublevar curratio;
	  if (nchanges==0) 
	  {
	    curratio=1.0;
	  }
	  else
	  {
	    doublevar calcq,calcqtot;
	    calcqtot=1.0;
	    int ichstart=dataptr->occupation_first_diff_change_from_last_det(f,recdet,s);
	    for (int ichange=0; ichange<ichstart; ichange++)
	    {
	      calcqtot*=invgammastore(ichange);
	    }

	    for (int ichange=ichstart; ichange<nchanges; ichange++)
	    {
		// Factor the ichangeth change
	      int erow=dataptr->occupation_changes(f,recdet,s)(ichange);
	      for (int j=0; j<nelectrons(s); j++)
	      {
		ustore(ichange)(j)=inverse(f,0,s)(j,erow);
	      }
	      for (int jchange=0;jchange<ichange;jchange++) 
	      {
		doublevar prod=dot(ustore(ichange),veestore(jchange))/invgammastore(jchange);
		for (int j=0; j<nelectrons(s); j++)
		{
		  ustore(ichange)(j)=ustore(ichange)(j)-prod*ustore(jchange)(j);
		}
	      }
		      
	      for(int j=0; j<nelectrons(s); j++)
	      {
		int jr=j+s*nelectrons(0);
		veestore(ichange)(j)=moVal(0,jr,dataptr->occupation(f,recdet,s)(erow))-moVal(0,jr,dataptr->occupation(f,0,s)(erow));
	      }
	      
	      calcq=1.0+dot(veestore(ichange),ustore(ichange));
	      calcqtot*=calcq;
	      invgammastore(ichange)=calcq;
	    } // end for ichange
	    curratio=calcqtot;
	  }
	  detVal(f,recdet,s)=curratio*detVal(f,0,s);
	}
      }
  }

#endif  

}

//------------------------------------------------------------------------


/*!
\todo
Add a lazy evaluation to this, so it only evaluates the
derivatives if something has really changed.

\bug
When the system is sufficiently large, the wavefunction
is either too small when we begin sampling, or too large
when we actually equilibrate, so we get under/overflows.
Solution:
Put the logarithms further up the chain, since in some
cases, we overflow on the value

*/

void FSlat_wf::getLap(Wavefunction_data * wfdata,
                     int e, Wf_return & lap)
{

  Array1 <doublevar> si(nfunc_, 0.0);
  Array2 <doublevar> vals(nfunc_,5,0.0);




  if(staticSample==1 && parent->optimize_mo==0 && parent->optimize_det==0) {
    //lap.phase=0;
    for(int f=0; f< nfunc_; f++) {
      for(int i=0; i< 5; i++) 
        vals(f,i)=saved_laplacian(e,f,i);

      si(f)=saved_laplacian(e,f,5);
    }
  }
  else {
    Slat_wf_data * dataptr;
    recast(wfdata, dataptr);

    int s=dataptr->spin(e);
    int opp=dataptr->opspin(e);

#ifndef RECDETGL
    for(int f=0; f< nfunc_; f++) {
      doublevar funcval=0;

      for(int det=0; det < ndet; det++) {
        funcval += dataptr->detwt(det)*detVal(f,det,s)*detVal(f,det,opp);
      }


      if(fabs(funcval) > 0) 
        vals(f,0)=log(fabs(funcval));
      else vals(f,0)=-1e3;
      si(f)=sign(funcval);


      for(int i=1; i< 5; i++) {
        lap.amp(f,i)=0;
        for(int det=0; det < ndet; det++) {
          doublevar temp=0;
          for(int j=0; j<nelectrons(s); j++) {
            temp+=moVal(i , e, dataptr->occupation(f,det,s)(j) )
                 *inverse(f,det,s)(dataptr->rede(e), j);
          }

          //Prevent catastrophe with a singular matrix.
          //Shouldn't happen much.
          if(detVal(f,det,s)*detVal(f,det,opp)==0)
            temp=0;
          
          if(ndet >1) 
            temp*=dataptr->detwt(det)*detVal(f,det, s)*detVal(f,det, opp);
          vals(f,i)+=temp;
          
          //vals(f,i)+=dataptr->detwt(det)*temp
          //          *detVal(f,det, s)*detVal(f,det, opp);
        }
       
        if(funcval==0)
          vals(f,i)=0;  //Prevent division by zero
        else if(ndet > 1) 
          vals(f,i)/=funcval;

      }
      
    }
#else
    for(int f=0; f< nfunc_; f++) {
      doublevar funcval=0;

      for(int det=0; det < ndet; det++) {
        funcval += dataptr->detwt(det).val()*detVal(f,det,s)*detVal(f,det,opp);
      }


      if(fabs(funcval) > 0) 
        vals(f,0)=log(fabs(funcval));
      else vals(f,0)=-1e3;
      si(f)=sign(funcval);

      for(int i=1; i< 5; i++) {
        lap.amp(f,i)=0;

	Array1 <doublevar> modet(nelectrons(s));	
	Array2 <doublevar> tinverse=inverse(f,0,s);

	for(int ic = 0; ic < nelectrons(s); ic++) {
	  modet(ic)=moVal(i,e,dataptr->occupation(f,0,s)(ic));
	}

	doublevar zero_ratio=1./InverseUpdateColumn(tinverse,
						    modet, dataptr->rede(e),
						    nelectrons(s));

	for(int idx=0; idx< ndet; idx++)
	{
	  int recdet=dataptr->evaluation_order(f,idx,s);
	  int nchanges=dataptr->occupation_nchanges(f,recdet,s);
	  doublevar curratio;
	  if (nchanges==0) 
	  {
	    curratio=zero_ratio;
	  }
	  else
	  {
	    doublevar calcq,calcqtot;
	    calcqtot=1.0;
	    int ichstart=dataptr->occupation_first_diff_change_from_last_det(f,recdet,s);
	      
	    for (int ichange=0; ichange<ichstart; ichange++)
	    {
	      calcqtot*=invgammastore(ichange);
	    }
	      
	    for (int ichange=ichstart; ichange<nchanges; ichange++)
	    {
		  // Factor the ichangeth change
	      int erow=dataptr->occupation_changes(f,recdet,s)(ichange);
	      for (int j=0; j<nelectrons(s); j++)
	      {
		ustore(ichange)(j)=tinverse(j,erow);
	      }
	      for (int jchange=0;jchange<ichange;jchange++) 
	      {
		doublevar prod=dot(ustore(ichange),veestore(jchange))/invgammastore(jchange);
		for (int j=0; j<nelectrons(s); j++)
		{
		  ustore(ichange)(j)=ustore(ichange)(j)-prod*ustore(jchange)(j);
		}
	      }
		  
	      for(int j=0; j<nelectrons(s); j++)
	      {
		int jr=j+s*nelectrons(0);
		veestore(ichange)(j)=moVal(0,jr,dataptr->occupation(f,recdet,s)(erow))-moVal(0,jr,dataptr->occupation(f,0,s)(erow));
	      }
	      int er=e-s*nelectrons(0);
	      veestore(ichange)(er)=moVal(i , e, dataptr->occupation(f,recdet,s)(erow) )-moVal(i , e, dataptr->occupation(f,0,s)(erow) );
	      
	      calcq=1.0+dot(veestore(ichange),ustore(ichange));
	      calcqtot*=calcq;
	      invgammastore(ichange)=calcq;
	    } // end for ichange
	    curratio=calcqtot*zero_ratio*detVal(f,0,s)/detVal(f,recdet,s);
	  }

	  doublevar temp=curratio;
          //Prevent catastrophe with a singular matrix.
          //Shouldn't happen much.
          if(detVal(f,recdet,s)*detVal(f,recdet,opp)==0)
            temp=0;
          
          if(ndet >1) 
            temp*=dataptr->detwt(recdet).val()*detVal(f,recdet, s)*detVal(f,recdet, opp);
          vals(f,i)+=temp;
        }

        if(funcval==0)
          vals(f,i)=0;  //Prevent division by zero
        else if(ndet > 1) 
          vals(f,i)/=funcval;
      }
      
    }
#endif
  }
  lap.setVals(vals, si);
  
}

//-------------------------------------------------------------------------

/*!
*/
void FSlat_wf::updateLap(Slat_wf_data * dataptr,
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
  dataptr->molecorb->updateLap(sample,e,s,updatedMoVal);

  for(int d=0; d< 5; d++)
    for(int i=0; i< updatedMoVal.GetDim(0); i++)
      moVal(d,e,i)=updatedMoVal(i,d);
  
  updateInverse(dataptr,e);
}

//-------------------------------------------------------------------------
