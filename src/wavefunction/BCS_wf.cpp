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

#include "Qmc_std.h"
#include "BCS_wf.h"
#include "MatrixAlgebra.h"
#include "Sample_point.h"
#include "BCS_wf_data.h"

//----------------------------------------------------------------------

BCS_wf::~BCS_wf() {
  if(jastcof) delete jastcof;
}

void BCS_wf::generateStorage(Wavefunction_storage * & wfstore)
{
  wfstore=new BCS_wf_storage;
  BCS_wf_storage * store;
  recast(wfstore, store);
  jast.generateStorage(store->jast_store);
  store->inverse=inverse;
  store->detVal=detVal;
  //store->moVal.Resize(moVal.GetDim(1),moVal.GetDim(2));
  store->derivatives.Resize(derivatives.GetDim(0),derivatives.GetDim(2));
  store->derivatives_2.Resize(derivatives.GetDim(0),derivatives.GetDim(2));
}


//----------------------------------------------------------------------

void BCS_wf::init(Wavefunction_data * wfdata)
{
  
  BCS_wf_data * dataptr;
  recast(wfdata, dataptr);
  recast(wfdata, parent);
  
  //nmo=dataptr->molecorb->getNmo();
  ndet=1;
  nelectrons.Resize(2);
  nelectrons=dataptr->nelectrons;
  
  int tote=nelectrons(0)+nelectrons(1);
  ndim=3;
  
  if(nelectrons(0) != nelectrons(1)) {
    error("For BCS wave function, the number of spin up electrons must"
          " be equal to the number of spin down electrons.");
  }
  
  spin.Resize(tote);
  spin=dataptr->spin;
  
  //moVal.Resize(tote, nmo,5);
  //updatedMoVal.Resize(nmo,5);
  
  detVal.Resize (ndet);
  inverse.Resize(ndet);
  
  for(int det=0; det < ndet; det++) {
    inverse(det).Resize(nelectrons(0), nelectrons(0));
    inverse(det)=0;
    for(int e=0; e< nelectrons(0); e++) {
      inverse(det)(e,e)=1;
      inverse(det)(e,e)=1;
    }
    
    detVal(det)=1;
  }
  
  
  electronIsStaleVal.Resize(tote);
  electronIsStaleLap.Resize(tote);
  
  electronIsStaleVal=0;
  electronIsStaleLap=0;
  updateEverythingVal=1;
  updateEverythingLap=1;
  sampleAttached=0;
  dataAttached=0;
  
  jast.init(&parent->jastdata);
  
  derivatives.Resize(2*nelectrons(0),nelectrons(1),4);
  
  assert(jastcof==NULL);
  jastcof=new BCS_jastrow_u;
  
}

//----------------------------------------------------------------------

/*!
 Behavior under staticSample:  if sample_static is set, then
 we ignore all electron moves in the main algorithm, and the only way
 to update based on a move is by the getParmDepVal function.  The
 regular update functions will not work in this case.
 */
void BCS_wf::notify(change_type change, int num)
{
  
  jast.notify(change,num);
  switch(change)
  {
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
      updateEverythingVal=1;
      updateEverythingLap=1;
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


void BCS_wf::saveUpdate(Sample_point * sample, int e,
                        Wavefunction_storage * wfstore)
{
  BCS_wf_storage * store;
  recast(wfstore, store);
  jast.saveUpdate(sample,e,store->jast_store);
  store->inverse=inverse;
  store->detVal=detVal;
  
  //
  // unpaired electrons go into extra one-body orbitals, this
  // functionality will be resurrected later
  //
  //int nmo=moVal.GetDim(1);
  //int nd=moVal.GetDim(2);
  //for(int m=0; m< nmo; m++) {
  //  for(int d=0; d< nd; d++) {
  //    store->moVal(m,d)=moVal(e,m,d);
  //  }
  //}
  
  if(spin(e)==0) {
    int nj=derivatives.GetDim(1);
    assert(nj==nelectrons(0));
    int nd=derivatives.GetDim(2);
    for(int j=0; j< nj; j++) {
      for(int d=0; d< nd; d++) {
        store->derivatives(j,d)=derivatives(e,j,d);
      }
    }
    int ep=e+nelectrons(0);
    for(int j=0; j< nj; j++) {
      for(int d=0; d< nd; d++) {
        store->derivatives(j+nelectrons(0),d)=derivatives(ep,j,d);
      }
    }
  }
  else {
    int ni=derivatives.GetDim(0);
    assert(ni==2*nelectrons(0));
    int nd=derivatives.GetDim(2);
    int rede=e-nelectrons(0);
    for(int i=0; i< ni; i++) {
      for(int d=0; d< nd; d++) {
        store->derivatives(i,d)=derivatives(i,rede,d);
      }
    }
  }
  
}

//----------------------------------------------------------------------

void BCS_wf::restoreUpdate(Sample_point * sample, int e,
                           Wavefunction_storage * wfstore)
{
  
  BCS_wf_storage * store;
  recast(wfstore, store);
  
  jast.restoreUpdate(sample,e,store->jast_store);
  detVal=store->detVal;
  inverse=store->inverse;
  
  //
  // unpaired electrons go into extra one-body orbitals, this
  // functionality will be resurrected later
  //
  //int nmo=moVal.GetDim(1);
  //int nd=moVal.GetDim(2);
  //for(int m=0; m< nmo; m++) {
  //  for(int d=0; d< nd; d++) {
  //    moVal(e,m,d)=store->moVal(m,d);
  //  }
  //}
  
  if(spin(e)==0) {
    int nj=derivatives.GetDim(1);
    assert(nj==nelectrons(0));
    int nd=derivatives.GetDim(2);
    for(int j=0; j< nj; j++) {
      for(int d=0; d< nd; d++) {
        derivatives(e,j,d)=store->derivatives(j,d);
      }
    }
    int ep=e+nelectrons(0);
    for(int j=0; j< nj; j++) {
      for(int d=0; d< nd; d++) {
        derivatives(ep,j,d)=store->derivatives(j+nelectrons(0),d);
      }
    }
  }
  else {
    int ni=derivatives.GetDim(0);
    assert(ni==2*nelectrons(0));
    int nd=derivatives.GetDim(2);
    int rede=e-nelectrons(0);
    for(int i=0; i< ni; i++) {
      for(int d=0; d< nd; d++) {
        derivatives(i,rede,d)=store->derivatives(i,d);
      }
    }
  }
  
  
  electronIsStaleLap=0;
  electronIsStaleVal=0;
  updateEverythingLap=0;
  updateEverythingVal=0;
  
  //calcLap(sample);
}

//----------------------------------------------------------------------
//Added by Matous

void BCS_wf::saveUpdate(Sample_point * sample, int e1, int e2,
                        Wavefunction_storage * wfstore)
{
  BCS_wf_storage * store;
  recast(wfstore, store);
  jast.saveUpdate(sample,e1,e2,store->jast_store);
  store->inverse=inverse;
  store->detVal=detVal;
  
  //
  // unpaired electrons go into extra one-body orbitals, this
  // functionality will be resurrected later
  //
  //int nmo=moVal.GetDim(1);
  //int nd=moVal.GetDim(2);
  //for(int m=0; m< nmo; m++) {
  //  for(int d=0; d< nd; d++) {
  //    store->moVal(m,d)=moVal(e,m,d);
  //  }
  //}
  
  // We assume the two electrons have opposite spins
  assert( spin(e1) != spin(e2) );
  
  int e_up, e_down;
  
  e_up = spin(e1) ? e2 : e1;
  e_down = spin(e1) ? e1 : e2;
  
  
  {
    int nj=derivatives.GetDim(1);
    assert(nj==nelectrons(0));
    int nd=derivatives.GetDim(2);
    for(int j=0; j< nj; j++) {
      for(int d=0; d< nd; d++) {
        store->derivatives(j,d)=derivatives(e_up,j,d);
      }
    }
    int ep=e_up+nelectrons(0);
    for(int j=0; j< nj; j++) {
      for(int d=0; d< nd; d++) {
        store->derivatives(j+nelectrons(0),d)=derivatives(ep,j,d);
      }
    }
  }
  {
    int ni=derivatives.GetDim(0);
    assert(ni==2*nelectrons(0));
    int nd=derivatives.GetDim(2);
    int rede=e_down-nelectrons(0);
    for(int i=0; i< ni; i++) {
      for(int d=0; d< nd; d++) {
        store->derivatives_2(i,d)=derivatives(i,rede,d);
      }
    }
  }
  
}

//----------------------------------------------------------------------

void BCS_wf::restoreUpdate(Sample_point * sample, int e1, int e2,
                           Wavefunction_storage * wfstore)
{
  
  BCS_wf_storage * store;
  recast(wfstore, store);
  
  jast.restoreUpdate(sample,e1,e2,store->jast_store);
  
  detVal=store->detVal;
  inverse=store->inverse;
  
  //
  // unpaired electrons go into extra one-body orbitals, this
  // functionality will be resurrected later
  //
  //int nmo=moVal.GetDim(1);
  //int nd=moVal.GetDim(2);
  //for(int m=0; m< nmo; m++) {
  //  for(int d=0; d< nd; d++) {
  //    moVal(e,m,d)=store->moVal(m,d);
  //  }
  //}
  
  // We assume the two electrons have opposite spins
  assert( spin(e1) != spin(e2) );
  
  int e_up, e_down;
  
  e_up = spin(e1) ? e2 : e1;
  e_down = spin(e1) ? e1 : e2;
  
  {
    int ni=derivatives.GetDim(0);
    assert(ni==2*nelectrons(0));
    int nd=derivatives.GetDim(2);
    int rede=e_down-nelectrons(0);
    for(int i=0; i< ni; i++) {
      for(int d=0; d< nd; d++) {
        derivatives(i,rede,d)=store->derivatives_2(i,d);
      }
    }
  }
  {
    int nj=derivatives.GetDim(1);
    assert(nj==nelectrons(0));
    int nd=derivatives.GetDim(2);
    for(int j=0; j< nj; j++) {
      for(int d=0; d< nd; d++) {
        derivatives(e_up,j,d)=store->derivatives(j,d);
      }
    }
    int ep=e_up+nelectrons(0);
    for(int j=0; j< nj; j++) {
      for(int d=0; d< nd; d++) {
        derivatives(ep,j,d)=store->derivatives(j+nelectrons(0),d);
      }
    }
  }
  
  
  electronIsStaleLap=0;
  electronIsStaleVal=0;
  updateEverythingLap=0;
  updateEverythingVal=0;
  
  
  //calcLap(sample);
}

//----------------------------------------------------------------------

void BCS_wf::updateVal(Wavefunction_data * wfdata,
                       Sample_point * sample)
{
  updateLap(wfdata,sample);
  return;
  
}

//----------------------------------------------------------------------

void BCS_wf::updateLap( Wavefunction_data * wfdata,
                       Sample_point * sample)
{
  assert(sampleAttached);
  assert(dataAttached);
  
  
  
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
        if ( abs(detVal(0))>0 ) {
          updateLap(e,sample);
        } else {
          calcLap(sample);
        }
        electronIsStaleLap(e)=0;
        electronIsStaleVal(e)=0;
      }
    }
    
  }
  
  
}

//-----------------------------------------------------------------------


void BCS_wf::getVal(Wavefunction_data * wfdata, int e,
                    Wf_return & val)
{
  assert(val.amp.GetDim(0) >=1);
  assert(val.amp.GetDim(1) >= 1);
  
  doublevar funcval=0;
  
  for(int det=0; det < ndet; det++) {
    funcval += detVal(det);
  }
  if(fabs(funcval) > 0)
    val.amp(0,0)=log(fabs(funcval));
  else val.amp(0,0)=-1e3;
  double si=sign(funcval);
  
  if(si>0)
    val.phase(0,0)=0;
  else val.phase(0,0)=pi;
  
}

//-------------------------------------------------------------------------
void BCS_wf::getSymmetricVal(Wavefunction_data * wfdata,
                             int e, Wf_return & val){
  
  Array1 <doublevar> si(1, 0.0);
  Array2 <doublevar> vals(1,1,0.0);
  val.setVals(vals, si);
}


//----------------------------------------------------------------------


//Note that the matrix is almost symmetric--
//we can use that to halve the computational time
void BCS_wf::calcLap(Sample_point * sample)
{
  
  int tote=nelectrons(0)+nelectrons(1);
  Array3 <doublevar> temp_der;
  Array2 <doublevar> temp_lap;
  for(int e=0; e< tote; e++) {
    sample->updateEIDist();
    //
    //parent->molecorb->updateLap(sample,e,0,updatedMoVal);
    //for(int i=0; i< updatedMoVal.GetDim(0); i++) {
    //  for(int d=0; d< 5; d++) {
    //    moVal(e,i,d)=updatedMoVal(i,d);
    //  }
    //}
  }
  
  Array2 <doublevar> onebody;
  jast.updateLap(&parent->jastdata,sample);
  jast.get_twobody(twobody);
  jast.get_onebody_save(onebody);
  
  int maxmatsize=max(nelectrons(0),nelectrons(1));
  Array2 <doublevar> modet(maxmatsize, maxmatsize);
  modet=0;
  for(int det=0; det < ndet; det++ ) {
    for(int i=0; i< nelectrons(0); i++) {
      for(int j=0; j< nelectrons(1); j++) {
        int jshift=nelectrons(0)+j;
        
        Array2 <doublevar> ugrad;
        jastcof->valGradLap(twobody,onebody,i,jshift,ugrad);
        modet(i,j)=parent->magnification_factor*ugrad(0,0);
        
        for(int d=0; d< 4; d++) {
          derivatives(i,j,d)=ugrad(0,d+1);
          derivatives(i+nelectrons(0),j,d)=ugrad(1,d+1);
        }
        
      }
      //
      // unpaired electrons go into extra one-body orbitals, this
      // functionality will be resurrected later
      //
      //for(int j=nelectrons(1); j< nelectrons(0); j++) {
      //  modet(i,j)=moVal(i,parent->occupation(det,0)(j),0);
      //}
    }
    
    
    detVal(det)=
    TransposeInverseMatrix(modet,inverse(det), nelectrons(0)).val();
    
  }
  
}

//------------------------------------------------------------------------

void BCS_wf::updateLap(int e, Sample_point * sample ) {
  
  int s=spin(e);
  sample->updateEIDist();
  //
  //parent->molecorb->updateLap(sample,e,0,updatedMoVal);
  //for(int i=0; i< updatedMoVal.GetDim(0); i++) {
  //  for(int d=0; d< 5; d++) {
  //    moVal(e,i,d)=updatedMoVal(i,d);
  //  }
  //}
  
  Array2 <doublevar> onebody;
  jast.updateLap(&parent->jastdata,sample);
  jast.get_twobody(twobody);
  jast.get_onebody_save(onebody);
  
  
  for(int det=0; det < ndet; det++ ) {
    Array1 <doublevar> detUpdate(nelectrons(0));
    
    if(s==0) { //------------spin up electron
      for(int j=0; j< nelectrons(1); j++) {
        int jshift=nelectrons(0)+j;
        Array2 <doublevar> ugrad;
        jastcof->valGradLap(twobody,onebody,e,jshift,ugrad);
        detUpdate(j)=parent->magnification_factor*ugrad(0,0);
        
        for(int d=0; d< 4; d++) {
          derivatives(e,j,d)=ugrad(0,d+1);
          derivatives(e+nelectrons(0),j,d)=ugrad(1,d+1);
        }
      }
      //
      // unpaired electrons go into extra one-body orbitals, this
      // functionality will be resurrected later
      //
      //for(int j=nelectrons(1); j< nelectrons(0); j++) {
      //  detUpdate(j)=moVal(e,parent->occupation(det,0)(j),0);
      //}
      detVal(det)/=InverseUpdateColumn(inverse(det),detUpdate,e,nelectrons(0));
      
    }
    else { //--------------------spin down electrons
      int eshift=e;
      int edown=e-nelectrons(0);
      
      for(int i=0; i< nelectrons(0); i++) {
        Array2 <doublevar> ugrad;
        jastcof->valGradLap(twobody,onebody,i,eshift,ugrad);
        detUpdate(i)=parent->magnification_factor*ugrad(0,0);
        
        for(int d=0; d< 4; d++) {
          derivatives(i,edown,d)=ugrad(0,d+1);
          derivatives(i+nelectrons(0),edown,d)=ugrad(1,d+1);
        }
      }
      detVal(det)/=InverseUpdateRow(inverse(det),detUpdate,
                                    edown,nelectrons(0));
      
    }
    
  }
  
}

//----------------------------------------------------------------------

/*!
 */

void BCS_wf::getLap(Wavefunction_data * wfdata,
                    int e, Wf_return & lap)
{
  
  getVal(wfdata,e,lap);
  
  for(int d=1; d< 5; d++) {
    if(spin(e)==0) {
      for(int det=0; det < ndet; det++) {
        doublevar temp=0;
        for(int j=0; j< nelectrons(1); j++) {
          temp+=derivatives(e,j,d-1)
          *inverse(det)(e,j);
        }
        //
        // unpaired electrons go into extra one-body orbitals, this
        // functionality will be resurrected later
        //
        //for(int j=nelectrons(1); j < nelectrons(0); j++) {
        //  //error("need to do spin-polarized derivatives");
        // temp+=moVal(e,parent->occupation(det,0)(j),d)
        //   *inverse(det)(e,j);
        //}
        if(ndet > 1) error("update BCS::getLap() for ndet > 1");
        lap.amp(0,d)=parent->magnification_factor*temp;
      }
    }
    else {
      for(int det=0; det < ndet; det++) {
        doublevar temp=0;
        int red=parent->rede(e);
        for(int i=0; i< nelectrons(0); i++) {
          temp+=derivatives(i+nelectrons(0),red,d-1)
          *inverse(det)(i,red);
        }
        lap.amp(0,d)=parent->magnification_factor*temp;
      }
      
    }
  }
  
  for(int d=1; d< 5; d++) {
    lap.phase(0,d)=0.0;
  }
  
}

//-------------------------------------------------------------------------

int BCS_wf::getParmDeriv(Wavefunction_data *wfdata , Sample_point * sample,
                         Parm_deriv_return & parm_derivatives ) {
  
  // analytic expressions below are valid only if jastrow can provide
  // analytic derivatives (not sure if it is enough to test this)
  if ( !parent->jastdata.supports(parameter_derivatives) ) return 0;
  
  //cout << "BCS_wf: doing analytic derivatives with respect to parameters."
  //     << endl ;
  
  int nparms_full=parent->nparms();
  int nparms_start=0;
  int nparms_end=nparms_full;
  int nparms=nparms_end-nparms_start;
  
  parm_derivatives.gradient.Resize(nparms);
  parm_derivatives.hessian.Resize(nparms, nparms);
  
  // gradients of the two-body orbital
  Array3 <doublevar> twobody_parm_deriv;
  jast.updateLap(&parent->jastdata,sample);
  jast.get_twobody_ParmDeriv(sample,twobody_parm_deriv);
  
  // gradient
  parm_derivatives.gradient=0.0;
  for(int det=0; det < ndet; det++) {
    for (int alpha=0; alpha<nparms; alpha++) {
      for(int i=0; i< nelectrons(0); i++) {
        for (int j=0; j< nelectrons(1); j++) {
          parm_derivatives.gradient(alpha)+=
          twobody_parm_deriv(i,nelectrons(0)+j,alpha)*inverse(det)(i,j)*
          parent->magnification_factor;
        }
      }
    }
  }
  
  // hessian
  parm_derivatives.hessian=0.0;
  for(int det=0; det < ndet; det++) {
    for (int alpha=0; alpha<nparms; alpha++) {
      for (int beta=0; beta<nparms; beta++) {
        parm_derivatives.hessian(alpha,beta)+=
        parm_derivatives.gradient(alpha)*parm_derivatives.gradient(beta);
        for(int i=0; i< nelectrons(0); i++) {
          for (int j=0; j< nelectrons(1); j++) {
            for(int k=0; k< nelectrons(0); k++) {
              for (int l=0; l< nelectrons(1); l++) {
                parm_derivatives.hessian(alpha,beta)-=
                twobody_parm_deriv(i,nelectrons(0)+j,alpha)*inverse(det)(k,j)*
                twobody_parm_deriv(k,nelectrons(0)+l,beta)*inverse(det)(i,l)*
                parent->magnification_factor*parent->magnification_factor;
              }
            }
          }
        }
      }
    }
  }
  
  return 1;
  
}


//-------------------------------------------------------------------------

void BCS_wf::plot1DInternals(Array1 <doublevar> & xdata,
                             vector <Array1 <doublevar> > & data,
                             vector <string> & desc,
                             string desc0) {
  
  jast.plot1DInternals(xdata,data,desc,"BCS, ");
  
}

//-------------------------------------------------------------------------
