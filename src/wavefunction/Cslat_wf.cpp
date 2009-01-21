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
#include "Cslat_wf.h"
#include "MatrixAlgebra.h"
#include "Sample_point.h"
#include "Slat_wf_data.h"

//----------------------------------------------------------------------


void Cslat_wf::generateStorage(Wavefunction_storage * & wfstore)
{
  wfstore=new Cslat_wf_storage;
  Cslat_wf_storage * store;
  recast(wfstore, store);
  store->moVal_temp.Resize (5,   nmo);
  store->detVal_temp.Resize(nfunc_, ndet, 2);
  store->inverse_temp.Resize(nfunc_, ndet, 2);
  for(int i=0; i< nfunc_; i++) {
    for(int det=0; det < ndet; det++) {
      for(int s=0; s<2; s++) {
        store->inverse_temp(i,det,s).Resize(nelectrons(s), nelectrons(s));
        store->inverse_temp(i,det,s)=0;
        store->detVal_temp(i,det,s)=1;
      }
    }
  }
}


//----------------------------------------------------------------------

void Cslat_wf::init(Wavefunction_data * wfdata)
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
void Cslat_wf::notify(change_type change, int num)
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
    init(parent);
    break;
  default:
    updateEverythingVal=1;
    updateEverythingLap=1;
  }
}

//----------------------------------------------------------------------

void Cslat_wf::save_for_static() {
  assert(staticSample==0);

  if(!parent->optimize_mo && !parent->optimize_det) {
    
    int totelectrons=nelectrons(0)+nelectrons(1);
    saved_laplacian.Resize(totelectrons, nfunc_, 6);
    Wf_return templap(nfunc_, 5);
    
    for(int e=0; e< totelectrons; e++) {
      getLap(parent, e, templap);
      for(int f=0; f< nfunc_; f++) {
        saved_laplacian(e,f,0)=templap.amp(f,0);
        saved_laplacian(e,f,5)=templap.phase(f,0);
        
      for(int i=1; i< 5; i++) {
        saved_laplacian(e,f,i)=templap.cvals(f,i);
      }
      }
    }
  }

}

//----------------------------------------------------------------------

void Cslat_wf::saveUpdate(Sample_point * sample, int e,
                         Wavefunction_storage * wfstore) {
  if(staticSample==0) {
    Cslat_wf_storage * store;
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

void Cslat_wf::restoreUpdate(Sample_point * sample, int e,
                            Wavefunction_storage * wfstore)
{

  if(staticSample==0) {
    Cslat_wf_storage * store;
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

void Cslat_wf::updateVal(Wavefunction_data * wfdata,
                        Sample_point * sample)
{
  assert(sampleAttached);
  assert(dataAttached);

  //Slat_wf_data * slatdata;
  //recast(wfdata, slatdata);

  if(staticSample==0 || parent->optimize_mo ) {
    if(updateEverythingVal==1) {
      calcVal(sample);
      updateEverythingVal=0;
      electronIsStaleVal=0;
    }
    else {
      for(int e=0; e< nelectrons(0)+nelectrons(1); e++) {
        if(electronIsStaleVal(e)) {
          updateVal(sample, e);
          electronIsStaleVal(e)=0;
        }
      }
    }
  }
}

//----------------------------------------------------------------------

void Cslat_wf::updateLap( Wavefunction_data * wfdata,
                        Sample_point * sample)
{
  assert(sampleAttached);
  assert(dataAttached);

  

  //Slat_wf_data * slatdata;
  //recast(wfdata, slatdata);
  if(staticSample==0 || parent->optimize_mo ) {
    if(updateEverythingLap==1) {
      calcLap(sample);
      updateEverythingVal=0;
      updateEverythingLap=0;
      electronIsStaleLap=0;
      electronIsStaleVal=0;
    }
    else {
      for(int e=0; e< nelectrons(0)+nelectrons(1); e++) {
        if(electronIsStaleLap(e)) {
          assert(!staticSample);
          updateLap(sample, e);
          electronIsStaleLap(e)=0;
          electronIsStaleVal(e)=0;
        }
      }

    }
  }
}

//----------------------------------------------------------------------


void Cslat_wf::storeParmIndVal(Wavefunction_data * wfdata, 
                               Sample_point * sample,
                              int e, Array1 <doublevar> & vals ) {
  assert(vals.GetDim(0) >=2);
  Wf_return newval(nfunc_,1);
  updateVal(wfdata, sample);
  getVal(wfdata, e, newval);
  vals(0)=newval.amp(0,0);
  vals(1)=newval.phase(0,0);

}

//----------------------------------------------------------------------

void Cslat_wf::getParmDepVal(Wavefunction_data * wfdata,
                            Sample_point * sample,
                            int e,
                            Array1 <doublevar> & oldval,
                            Wf_return & newval) {

  assert(oldval.GetDim(0) >=2);
  assert(newval.amp.GetDim(1) >= 1);
  assert(newval.amp.GetDim(0) >= nfunc_);
  int counter=0;
  newval.amp(0,0)=oldval(counter++);
  newval.phase(0,0)=oldval(counter++);
}


//------------------------------------------------------------------------


void Cslat_wf::calcVal(Sample_point * sample)
{
  //Hmm, I don't completely understand why, but something is not 
  //completely clean, so we can't just cycle through and updateVal
  //This is actually probably the best way to do it anyway, since it 
  //should in theory be faster, and it gives us a clean start
  calcLap(sample);

}

//------------------------------------------------------------------------
/*!

*/
void Cslat_wf::updateVal(Sample_point * sample,int e)
{

  sample->updateEIDist();
  int s=parent->spin(e);

  for(int f=0; f< nfunc_; f++) {
    for(int det=0; det < ndet; det++) {
      if(cabs(detVal(f, det, s))==0) {
        cout << "updateVal::WARNING: determinant zero!" << endl;
        calcLap(sample);
        return;
      }
    }
  }

  dcomplex ratio;

  int maxmatsize=max(nelectrons(0),nelectrons(1));
  Array1 <dcomplex> modet(maxmatsize);

  //update all the mo's that we will be using.
  parent->cmolecorb->updateVal(sample, e, s,
                              updatedMoVal);


  for(int f=0; f< nfunc_; f++)  {
    for(int det=0; det< ndet; det++)  {
      //fill the molecular orbitals for this
      //determinant
      for(int i = 0; i < nelectrons(s); i++) {
        modet(i)=updatedMoVal(parent->occupation(f,det,s)(i),0 );
      }
      ratio=1./InverseUpdateColumn(inverse(f,det,s),
                                   modet, parent->rede(e),
                                   nelectrons(s));
      detVal(f,det, s)=ratio*detVal(f,det, s);
    }
  }


  for(int i=0; i< updatedMoVal.GetDim(0); i++)
    moVal(0,e,i)=updatedMoVal(i,0);

}

//------------------------------------------------------------------------


void Cslat_wf::getVal(Wavefunction_data * wfdata, int e,
                     Wf_return & val)
{

  assert(val.amp.GetDim(0) >= nfunc_);

  Array2 <dcomplex> vals(nfunc_, 1);
  Array1 <doublevar> phase(nfunc_);


  if(staticSample==1 && parent->optimize_mo==0 && parent->optimize_det==0) {
    for(int f=0; f< nfunc_; f++) {
      vals(f,0)=saved_laplacian(e,f,0);
      phase(f)=saved_laplacian(e,f,5).real();
    }
  }
  else {
    int s=parent->spin(e);
    int opp=parent->opspin(e);
    for(int f=0; f< nfunc_; f++) {

      Array1 <dcomplex> funcval(5);
      funcval=dcomplex(0.0, 0.0);
      for(int det=0; det < ndet; det++) {
        funcval(0) += parent->detwt(det)*detVal(f,det,s)*detVal(f,det,opp);
      }

      doublevar amp=cabs(funcval(0));


      vals(f,0)=dcomplex(log(amp),0.0);
      if(fabs(funcval(0).imag()) > 1e-8) {
        phase(f)=atan(funcval(0).imag()/funcval(0).real());
	// JK: atan alone has a period of only pi but phase has 2pi period,
	// following should fix this
	if ( funcval(0).real() < 0.0 ) phase(f)+=pi;
      }
      else {
        phase(f)=-.5*pi*(1-sign(funcval(0).real()));
      }
    }
  }

  val.setVals(vals,phase);
}

//----------------------------------------------------------------------

void Cslat_wf::getDensity(Wavefunction_data * wfdata, int e,
                         Array2 <doublevar> & dens)
{

  assert(dens.GetDim(0) >= nfunc_);
  //Slat_wf_data * dataptr;
  //recast(wfdata, dataptr);

  int s=parent->spin(e);

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
      doublevar tmp=cabs(moVal(0,e,parent->occupation(f,det,s)(j)));
      dens(f,0)+=tmp*tmp;
    }
  }

}

//----------------------------------------------------------------------------


void Cslat_wf::calcLap(Sample_point * sample)
{

  for(int e=0; e< nelectrons(0)+nelectrons(1); e++)  {
    int s=parent->spin(e);
    sample->updateEIDist();

    //update all the mo's that we will be using, using the lists made in
    //Slat_wf_data(one for each spin).
    parent->cmolecorb->updateLap(sample, e, s,
                                updatedMoVal);
    for(int d=0; d< 5; d++)  {
      for(int i=0; i< updatedMoVal.GetDim(0); i++) {
        moVal(d,e,i)=updatedMoVal(i,d);
      }
    }
  }

  int maxmatsize=max(nelectrons(0),nelectrons(1));
  Array2 <dcomplex> modet(maxmatsize, maxmatsize);

  for(int f=0; f< nfunc_; f++)   {
    for(int det=0; det < ndet; det++ ) {
      for(int s=0; s< 2; s++ ) {


        for(int e=0; e< nelectrons(s); e++) {
          int curre=s*nelectrons(0)+e;
          for(int i=0; i< nelectrons(s); i++) {
            modet(e,i)=moVal(0,curre, parent->occupation(f,det,s)(i));
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

void Cslat_wf::getLap(Wavefunction_data * wfdata,
                     int e, Wf_return & lap)
{

  assert(lap.amp.GetDim(1) >= 5);
  assert(lap.amp.GetDim(0) >= nfunc_);

  Array2 <dcomplex> vals(nfunc_, 5,0.0);
  Array1 <doublevar> phase(nfunc_,0.0);
  int zero_func=0;

  if(staticSample==1 && parent->optimize_mo==0 && parent->optimize_det==0) {
    for(int f=0; f< nfunc_; f++) {
      for(int i=0; i< 5; i++)
        vals(f,i)=saved_laplacian(e,f,i);
      phase(f)=saved_laplacian(e,f,5).real();
    }
  }
  else {
    int s=parent->spin(e);
    int opp=parent->opspin(e);

    for(int f=0; f< nfunc_; f++)
    {

      Array1 <dcomplex> funcval(5);
      funcval=dcomplex(0.0, 0.0);

      for(int det=0; det < ndet; det++) {
        funcval(0) += parent->detwt(det)*detVal(f,det,s)*detVal(f,det,opp);
      }

      for(int i=1; i< 5; i++) {
        for(int det=0; det < ndet; det++) {
          dcomplex temp=0;
          for(int j=0; j<nelectrons(s); j++) {
            temp+=moVal(i , e, parent->occupation(f,det,s)(j) )
                 *inverse(f,det,s)(parent->rede(e), j);
          }

          //Prevent catastrophe with a singular matrix.
          //Shouldn't happen much.
          if(detVal(f,det,s)*detVal(f,det,opp)==dcomplex(0.0, 0.0))
          temp=0;

          funcval(i)+=parent->detwt(det)*temp
                         *detVal(f,det, s)*detVal(f,det, opp);
        }
      }

      doublevar amp=cabs(funcval(0));
      //cout << "amp " << amp << endl;
      if(fabs(amp) > 1e-16) {
        vals(f,0)=dcomplex(log(amp),0.0);
        if(fabs(funcval(0).imag()) > 1e-8) {
          phase(f)=atan(funcval(0).imag()/funcval(0).real());
	  // JK: atan alone has a period of only pi but phase has 2pi period,
	  // following should fix this
	  if ( funcval(0).real() < 0.0 ) phase(f)+=pi;
        }
        else {
          phase(f)=-.5*pi*(1-sign(funcval(0).real()));
        }
        for(int i=1; i< 5; i++) {
          vals(f,i)=funcval(i)/funcval(0);

        }
      }
      else {
        zero_func=1;
      }
         
    }
  }

  if(zero_func) {
        lap.cvals=0;
        lap.amp=0;
        lap.phase=0;
        for(int w=0; w< nfunc_; w++) 
          lap.amp(w,0)=-1e3;
          
  }
  else 
    lap.setVals(vals,phase);


}

//-------------------------------------------------------------------------

/*!
*/
void Cslat_wf::updateLap(Sample_point * sample,
                        int e ) {
  
  assert(parent != NULL);

  int s=parent->spin(e);
  sample->updateEIDist();
  dcomplex ratio;

  int maxmatsize=max(nelectrons(0),nelectrons(1));
  static Array1 <dcomplex> modet(maxmatsize);

  //check to make sure the determinant isn't zero
  for(int f=0; f< nfunc_; f++) {
    for(int det=0; det < ndet; det++) {
      if(detVal(f, det, s)==dcomplex(0.,0.)) {
        cout << "updateLap::WARNING: determinant zero!" << endl;
        Array1 <doublevar> pos(3);
        sample->getElectronPos(e,pos);
        cout << "e= " << e<<" pos " <<  pos(0) << endl;
        calcLap(sample);
        return;
      }
    }
  }


  //update all the mo's that we will be using.
  parent->cmolecorb->updateLap(sample, e,
                              s,
                              updatedMoVal);

  for(int f=0; f<nfunc_; f++)
  {
    for(int det=0; det< ndet; det++)
    {

      //fill the molecular orbitals for this
      //determinant
      for(int i = 0; i < nelectrons(s); i++) {
        modet(i)=updatedMoVal(parent->occupation(f,det,s)(i),0);
      }

      dcomplex tmpratio=InverseUpdateColumn(inverse(f,det,s),
                                             modet,parent->rede(e),
                                             nelectrons(s));
      if(tmpratio==dcomplex(0.,0.))
        ratio=0;
      else ratio=1./tmpratio;

      detVal(f,det, s)=ratio*detVal(f,det, s);
    }
  }


  for(int d=0; d< 5; d++)
    for(int i=0; i< updatedMoVal.GetDim(0); i++)
      moVal(d,e,i)=updatedMoVal(i,d);


  //cout << "MO determinant " << endl;
  //for(int i=0; i< 2; i++) {
  //  for(int j=0; j< 2; j++) {
  //    cout << moVal(0,i,j) <<"       " ;
  //  }
  //  cout << endl;
  //}

  //cout << "determinant value " << detVal(0,0,0) << endl;

}

//-------------------------------------------------------------------------
