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
#include "Slat_wf.h"
#include "MatrixAlgebra.h"
#include "Sample_point.h"
#include "Slat_wf_data.h"

//----------------------------------------------------------------------


void Slat_wf::generateStorage(Wavefunction_storage * & wfstore)
{
  wfstore=new Slat_wf_storage;
  Slat_wf_storage * store;
  recast(wfstore, store);
  store->moVal_temp.Resize (5,   nmo);
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
        store->detVal_temp(i,det,s)=1;
      }
    }
  }
}


//----------------------------------------------------------------------

void Slat_wf::init(Wavefunction_data * wfdata)
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

        detVal(i,det,s)=1;
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
}

//----------------------------------------------------------------------

/*!
Behavior under staticSample:  if sample_static is set, then
we ignore all electron moves in the main algorithm, and the only way
to update based on a move is by the getParmDepVal function.  The
regular update functions will not work in this case.
*/
void Slat_wf::notify(change_type change, int num)
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

void Slat_wf::save_for_static() {
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

void Slat_wf::saveUpdate(Sample_point * sample, int e,
                         Wavefunction_storage * wfstore)
{

  if(staticSample==0) {
    Slat_wf_storage * store;
    recast(wfstore, store);

    int s=spin(e);

    for(int f=0; f< nfunc_; f++) {
      for(int det=0; det<ndet; det++) {
        store->inverse_temp(f,det,s)=inverse(f,det,s);
        store->detVal_temp(f,det,s)=detVal(f,det,s);
      }
    }


    for(int d=0; d< 5; d++) {
      for(int i=0; i< moVal.GetDim(2); i++) {
        store->moVal_temp(d,i)=moVal(d,e,i);
      }
    }
  }


}

//----------------------------------------------------------------------

void Slat_wf::restoreUpdate(Sample_point * sample, int e,
                            Wavefunction_storage * wfstore)
{

  if(staticSample==0) {
    Slat_wf_storage * store;
    recast(wfstore, store);
    int s=spin(e);

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

    electronIsStaleVal(e)=0;
    electronIsStaleLap(e)=0;
  }

}

//----------------------------------------------------------------------

void Slat_wf::updateVal(Wavefunction_data * wfdata,
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

void Slat_wf::updateLap( Wavefunction_data * wfdata,
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


void Slat_wf::storeParmIndVal(Wavefunction_data * wfdata, Sample_point * sample,
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

void Slat_wf::getParmDepVal(Wavefunction_data * wfdata,
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

int Slat_wf::getParmDeriv(Wavefunction_data *  wfdata, 
			  Sample_point * sample ,
			  Parm_deriv_return & derivatives){
  
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
      assert(nparms<parent->ncsf);
      int counter=0;
      derivatives.gradient=0.0;
      for(int csf=0; csf< parent->ncsf; csf++) {
        if(csf > nparms_start && csf <= nparms_end ){
          for(int j=1;j<parent->CSF(csf).GetDim(0);j++){
            //cout <<" csf "<<csf<<" j "<<j<<" counter "<<counter<<endl;
            derivatives.gradient(csf-1-nparms_start)+=parent->CSF(csf)(j)*detVal(0,counter,0)*detVal(0,counter,1);
            counter++;
          }
        }
        else{
          counter+=parent->CSF(csf).GetDim(0)-1;
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
        if(det > nparms_start && det <= nparms_end )
          derivatives.gradient(det-1-nparms_start)=detVal(0,det,0)*detVal(0,det,1);
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
}


//------------------------------------------------------------------------


void Slat_wf::calcVal(Slat_wf_data * dataptr, Sample_point * sample)
{
  //Hmm, I don't completely understand why, but something is not 
  //completely clean, so we can't just cycle through and updateVal
  //This is actually probably the best way to do it anyway, since it 
  //should in theory be faster, and it gives us a clean start
  calcLap(dataptr, sample);

}

//------------------------------------------------------------------------
/*!

*/
void Slat_wf::updateVal( Slat_wf_data * dataptr, Sample_point * sample,int e)
{

  /* Note on the updates(from LKW):  we don't actually need to update the inverse(which is N^2) for the next value;
     we can get away with saving the MO values and updating the inverse only when another electron is 
     updated.  This could save some time particularly for the pseudopotential evaluation.  For the moment, 
     I haven't seen the inverse update take a significant amount of time even up to around 1000 electrons, 
     so I'll leave the simpler implementation.
   */
  assert(dataptr != NULL);
  sample->updateEIDist();
  int s=dataptr->spin(e);
  /*
  for(int f=0; f< nfunc_; f++) {
    for(int det=0; det < ndet; det++) {
      if(fabs(detVal(f, det, s))< 1e-15) {
	
        cout << "updateVal::WARNING: determinant zero: " 
	     << detVal(f,det,s)
	     << " func " << f << " spin " << s << "  det " << det 
	     << endl;
        calcLap(dataptr, sample);
        return;
      }
    }
  }
  */

  doublevar ratio;

  int maxmatsize=max(nelectrons(0),nelectrons(1));
  Array1 <doublevar> modet(maxmatsize);

  //update all the mo's that we will be using.
  dataptr->molecorb->updateVal(sample, e, s,
                              updatedMoVal);

  for(int i=0; i< updatedMoVal.GetDim(0); i++)
    moVal(0,e,i)=updatedMoVal(i,0);



  //If the determinant is zero, we can't do the updates,
  //so instead we have to take the total determinant.
  //This will cost a bit, but it should be relatively rare.
  for(int f=0; f< nfunc_; f++)  {
    for(int det=0; det< ndet; det++)  {
      //fill the molecular orbitals for this
      //determinant
      if(fabs(detVal(f,det,s)) > 0) { 
	for(int i = 0; i < nelectrons(s); i++) {
	  modet(i)=updatedMoVal(dataptr->occupation(f,det,s)(i),0 );
	}
	
	
	ratio=1./InverseUpdateColumn(inverse(f,det,s),
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
	//cout << "determinant was zero: new one " << detVal(f,det,s) << endl;
        //calcLap(dataptr, sample);
	//cout << "test " << detVal(f,det,s) << endl;
      }
    }
  }


  //for(int i=0; i< updatedMoVal.GetDim(0); i++)
  //  moVal(0,e,i)=updatedMoVal(i,0);

}

//------------------------------------------------------------------------


void Slat_wf::getVal(Wavefunction_data * wfdata, int e,
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
        funcval += dataptr->detwt(det)*detVal(f,det,s)*detVal(f,det,opp);
      }
      if(fabs(funcval) > 0) 
        vals(f,0)=log(fabs(funcval));
      else vals(f,0)=-1e3;
      si(f)=sign(funcval);
    }
  }

  val.setVals(vals, si);

}

//----------------------------------------------------------------------

void Slat_wf::getDensity(Wavefunction_data * wfdata, int e,
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


void Slat_wf::calcLap(Slat_wf_data * dataptr, Sample_point * sample)
{

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

void Slat_wf::getLap(Wavefunction_data * wfdata,
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
  }

  lap.setVals(vals, si);

}

//-------------------------------------------------------------------------

/*!
*/
void Slat_wf::updateLap(Slat_wf_data * dataptr,
                        Sample_point * sample,
                        int e ) {
  
  assert(dataptr != NULL);

  int s=dataptr->spin(e);
  sample->updateEIDist();
  doublevar ratio;

  int maxmatsize=max(nelectrons(0),nelectrons(1));
  static Array1 <doublevar> modet(maxmatsize);

  //check to make sure the determinant isn't zero
  for(int f=0; f< nfunc_; f++)
  {
    for(int det=0; det < ndet; det++)
    {
      if(!(fabs(detVal(f, det, s)))> 0)
      {
        cout << "updateLap::WARNING: determinant zero!"<< detVal(f,det,s) << endl;
        calcLap(dataptr, sample);
        return;
      }
    }
  }


  //update all the mo's that we will be using.
  dataptr->molecorb->updateLap(sample, e,
                              s,
                              updatedMoVal);

  for(int f=0; f<nfunc_; f++)
  {
    for(int det=0; det< ndet; det++)
    {

      //fill the molecular orbitals for this
      //determinant
      for(int i = 0; i < nelectrons(s); i++)
        modet(i)=updatedMoVal(dataptr->occupation(f,det,s)(i),0);

      doublevar tmpratio=InverseUpdateColumn(inverse(f,det,s),
                                             modet,dataptr->rede(e),
                                             nelectrons(s));
      if(tmpratio==0)
        ratio=0;
      else ratio=1./tmpratio;

      detVal(f,det, s)=ratio*detVal(f,det, s);
    }
  }


  for(int d=0; d< 5; d++)
    for(int i=0; i< updatedMoVal.GetDim(0); i++)
      moVal(d,e,i)=updatedMoVal(i,d);

}

//-------------------------------------------------------------------------
