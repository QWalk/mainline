//--------------------------------------------------------------------------
// src/Backflow_wf_calc.cpp
//
//
#include "Qmc_std.h"
#include "Backflow_wf.h"
#include "MatrixAlgebra.h"
#include "Sample_point.h"
#include "Backflow_wf_data.h"

//----------------------------------------------------------------------


void Backflow_wf::generateStorage(Wavefunction_storage * & wfstore)
{
  wfstore=new Backflow_wf_storage;
  Backflow_wf_storage * store;
  recast(wfstore, store);
  store->gradlap=gradlap;
  store->detVal=detVal;
}


//----------------------------------------------------------------------

void Backflow_wf::init(Wavefunction_data * wfdata)
{

  Backflow_wf_data * dataptr;
  recast(wfdata, dataptr);
  recast(wfdata, parent);

  nmo=dataptr->bfwrapper.nmo();
  ndet=dataptr->dkeeper.ndet;
  nelectrons.Resize(2);
  nelectrons=dataptr->dkeeper.nelectrons;

  int tote=nelectrons(0)+nelectrons(1);
  ndim=3;

  spin.Resize(tote);
  spin=dataptr->dkeeper.spin;

  //Properties and intermediate calculation storage.
  moVal.Resize(tote, nmo,10);
  updatedMoVal.Resize(nmo,10);

  detVal.Resize (ndet, 2);
  inverse.Resize(ndet, 2);

  for(int det=0; det < ndet; det++) {
    for(int s=0; s<2; s++) {
      inverse(det,s).Resize(nelectrons(s), nelectrons(s));
      inverse(det,s)=0;
      for(int e=0; e< nelectrons(s); e++) {
	inverse(det,s)(e,e)=1;
	inverse(det,s)(e,e)=1;
      }
      
      detVal(det,s)=1;
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


  coor_grad.Resize(tote,tote,ndim,ndim);
  coor_lap.Resize(tote,tote,ndim);

  gradlap.Resize(tote,5);

  jast.init(&parent->bfwrapper.jdata);
  jast.keep_ion_dep();
}

//----------------------------------------------------------------------

/*!
Behavior under staticSample:  if sample_static is set, then
we ignore all electron moves in the main algorithm, and the only way
to update based on a move is by the getParmDepVal function.  The
regular update functions will not work in this case.
*/
void Backflow_wf::notify(change_type change, int num)
{
  
  jast.notify(change,num);
  switch(change)
  {
  case electron_move:
    //if(staticSample==0) {
      electronIsStaleVal(num)=1;
      electronIsStaleLap(num)=1;
      //}
    break;
  case all_electrons_move:
    //if(staticSample==0) {
      updateEverythingVal=1;
      updateEverythingLap=1;
      //}
    break;
  case wf_parm_change:  
  case all_wf_parms_change:
    //if(parent->optimize_backflow  ) {
      updateEverythingVal=1;
      updateEverythingLap=1;
      //}
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
    //if(!parent->optimize_backflow) {
      //save_for_static();
      //staticSample=1;
      //}
    break;
  case sample_dynamic:
    staticSample=0;
    //init(parent);
    updateEverythingVal=1;
    updateEverythingLap=1;
    break;
  default:
    updateEverythingVal=1;
    updateEverythingLap=1;
  }

}

//----------------------------------------------------------------------


void Backflow_wf::saveUpdate(Sample_point * sample, int e,
                         Wavefunction_storage * wfstore)
{
  Backflow_wf_storage * store;
  recast(wfstore, store);
  store->gradlap=gradlap;
  store->detVal=detVal;
	
}

//----------------------------------------------------------------------

void Backflow_wf::restoreUpdate(Sample_point * sample, int e,
                            Wavefunction_storage * wfstore)
{
  Backflow_wf_storage * store;
  recast(wfstore, store);
  gradlap=store->gradlap;
  detVal=store->detVal;
	  
  electronIsStaleLap=0;
  electronIsStaleVal=0;
  updateEverythingLap=0;
  updateEverythingVal=0;
}

//----------------------------------------------------------------------

void Backflow_wf::updateVal(Wavefunction_data * wfdata,
                        Sample_point * sample)
{

  assert(sampleAttached);
  assert(dataAttached);

  Backflow_wf_data * slatdata;
  recast(wfdata, slatdata);

  if(updateEverythingVal){ 
    calcVal(sample);
    updateEverythingVal=0;
    electronIsStaleVal=0;
  }
  else { 
    for(int e=0; e< nelectrons(0)+nelectrons(1); e++) {
      if(electronIsStaleVal(e)) {
	calcVal(sample);
	//updateVal(e,sample);
	electronIsStaleVal=0;
      }
    }
  }

}

//----------------------------------------------------------------------

void Backflow_wf::updateLap( Wavefunction_data * wfdata,
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
	calcLap(sample);
	electronIsStaleLap(e)=0;
	electronIsStaleVal(e)=0;
      }
    }
    
  }


}

//----------------------------------------------------------------------


void Backflow_wf::storeParmIndVal(Wavefunction_data * wfdata, Sample_point * sample,
                              int e, Array1 <doublevar> & vals )
{
  /*
  if(parent->optimize_backflow) {

  }
  else if(parent->optimize_det) {
    assert(vals.GetDim(0) >= 2*parent->detwt.GetDim(0));
    updateVal(wfdata, sample);
    int count=0;
    for(int det=0; det < ndet; det++) {
      vals(count++)=detVal(det,0);
      vals(count++)=detVal(det,1);
    }
  }
  else {
    assert(vals.GetDim(0) >=2);
    Wf_return newval(1,1);
    updateVal(wfdata, sample);
    getVal(wfdata, e, newval);
    vals(0)=newval.amp(0,0);
    vals(1)=newval.phase(0,0);
  }
  */
}

//----------------------------------------------------------------------

void Backflow_wf::getParmDepVal(Wavefunction_data * wfdata,
                            Sample_point * sample,
                            int e,
                            Array1 <doublevar> & oldval,
                            Wf_return & newval)
{
  updateVal(wfdata,sample);
  getVal(wfdata,e,newval);

  /*
  if(parent->optimize_backflow) {
    updateVal(wfdata, sample);
    getVal(wfdata, e, newval);
  }
  else if(parent->dkeeper.optimize_det) {
    assert(oldval.GetDim(0) >=2*ndet);
    doublevar tempval=0;
    int count=0;
    for(int det=0; det < ndet; det++) {
      tempval+=parent->dkeeper.detwt(det)*oldval(count)*oldval(count+1);
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
  */
}


//-----------------------------------------------------------------------

int Backflow_wf::getParmDeriv(Wavefunction_data *  wfdata, 
			  Sample_point * sample ,
			  Parm_deriv_return & derivatives){
  return 0; 
}


//------------------------------------------------------------------------

void Backflow_wf::updateVal(int e, Sample_point * sample) { 
  int nlist; Array1 <int> list;
  parent->bfwrapper.getNeighbors(sample,jast,e,list,nlist);
  cout << "nlist " << nlist << endl;
  calcVal(sample);
}

void Backflow_wf::calcVal(Sample_point * sample)
{
  int tote=nelectrons(0)+nelectrons(1);

  for(int e=0; e< tote; e++) { 
    int s=spin(e);
    sample->updateEIDist();
    parent->bfwrapper.updateVal(sample,jast, e,s,updatedMoVal);
    for(int i=0; i< updatedMoVal.GetDim(0); i++) {
      for(int d=0; d< 1; d++) { 
        moVal(e,i,d)=updatedMoVal(i,d);
	//if(s==1) cout << "moval " << i << " d " << d 
	//	      << "   " << updatedMoVal(i,d) << endl;
      }
    }

  } 
     

  int maxmatsize=max(nelectrons(0),nelectrons(1));
  Array2 <doublevar> modet(maxmatsize, maxmatsize);
  //cout << "doing determinant " << endl;
  for(int det=0; det < ndet; det++ ) {
    for(int s=0; s< 2; s++ ) {
      
      
      for(int e=0; e< nelectrons(s); e++) {
	int curre=s*nelectrons(0)+e;
	for(int i=0; i< nelectrons(s); i++) {
	  modet(e,i)=moVal(curre, parent->occupation(det,s)(i),0);
	}
      }
      log_real_value tmp=TransposeInverseMatrix(modet,inverse(det,s), nelectrons(s));
      detVal(det,s)=tmp.val();
    }
  }
}

//------------------------------------------------------------------------
/*!

*/

//------------------------------------------------------------------------


void Backflow_wf::getVal(Wavefunction_data * wfdata, int e,
                     Wf_return & val)
{

  Array1 <doublevar> si(1, 0.0);
  Array2 <doublevar> vals(1,1,0.0);

  assert(val.amp.GetDim(0) >=1);
  assert(val.amp.GetDim(1) >= 1);


  doublevar funcval=0;
  
  for(int det=0; det < ndet; det++) {
    //cout << "getVal: detVal "<<detVal(det,0)<<"  "<<detVal(det,1)<<endl;
    funcval += parent->dkeeper.detwt(det)*detVal(det,0)*detVal(det,1);
  }
  if(fabs(funcval) > 0)
    vals(0,0)=log(fabs(funcval));
  else vals(0,0)=-1e3;
  si(0)=sign(funcval);

  if(si(0)>0) 
    val.phase(0,0)=0;
  else val.phase(0,0)=pi;
  val.amp(0,0)=vals(0,0);

}

//-----------------------------------------------------------------------
void Backflow_wf::getSymmetricVal(Wavefunction_data * wfdata,
		     int e, Wf_return & val){

  Array1 <doublevar> si(1, 0.0);
  Array2 <doublevar> vals(1,1,0.0);
  val.setVals(vals, si);
} 


//----------------------------------------------------------------------

void Backflow_wf::getDensity(Wavefunction_data * wfdata, int e,
                         Array2 <doublevar> & dens)
{
  error("No density for Backflow now..");

}

//----------------------------------------------------------------------------

//Watch out that some of the indices are reversed from Kwon--
//this is to make the updates a bit more convienent.
void Backflow_wf::calcLap(Sample_point * sample)
{


  

  int tote=nelectrons(0)+nelectrons(1);
  Array3 <doublevar> temp_der;
  Array2 <doublevar> temp_lap;
  for(int e=0; e< tote; e++) { 
    int s=spin(e);
    sample->updateEIDist();
    parent->bfwrapper.updateLap(sample,jast,e,s,updatedMoVal,
				temp_der, temp_lap);
    for(int i=0; i< updatedMoVal.GetDim(0); i++) {
      for(int d=0; d< 10; d++) { 
        moVal(e,i,d)=updatedMoVal(i,d);
	//if(s==1) cout << "moval " << i << " d " << d 
	//	      << "   " << updatedMoVal(i,d) << endl;
      }
    }

    for(int i=0; i< tote; i++) { 
      for(int a=0; a< 3; a++) {
	for(int b=0; b < 3; b++) { 
	  coor_grad(e,i,a,b)=temp_der(i,a,b);
	}
	coor_lap(e,i,a)=temp_lap(i,a);
      }
    }
  } 
     

  int maxmatsize=max(nelectrons(0),nelectrons(1));
  Array2 <doublevar> modet(maxmatsize, maxmatsize);
  //cout << "doing determinant " << endl;
  for(int det=0; det < ndet; det++ ) {
    for(int s=0; s< 2; s++ ) {
      
      
      for(int e=0; e< nelectrons(s); e++) {
	int curre=s*nelectrons(0)+e;
	for(int i=0; i< nelectrons(s); i++) {
	  modet(e,i)=moVal(curre, parent->occupation(det,s)(i),0);
	}
      }
      log_real_value tmp=TransposeInverseMatrix(modet,inverse(det,s), nelectrons(s));
      detVal(det,s)=tmp.val();
    }
  }


  gradlap=0;
  Array1 <doublevar> gradmod2(tote);
  Array2 <Array3 <doublevar> > F(ndet,2); //(det,spin)
  doublevar funcval=0;
  gradmod2=0;
  Array2 <doublevar> grad(tote,ndim);
  for(int det=0;det<ndet;det++){
    grad=0;
    //cout << "detVal "<<detVal(det,0)<<"  "<<detVal(det,1)<<endl;
    funcval+=parent->dkeeper.detwt(det)*detVal(det,0)*detVal(det,1);
    for(int s=0; s< 2; s++) { 
      F(det,s).Resize(nelectrons(s),nelectrons(s),ndim);
      F(det,s)=0;
    }
    if(detVal(det,0)*detVal(det,1)!=0){
      for(int s=0; s< 2; s++) { 
	Array3 <doublevar> & fref(F(det,s));
	for(int i=0; i< nelectrons(s); i++) { 
	  for(int j=0; j< nelectrons(s); j++) { 
	    int jp=j+s*nelectrons(0);
	    for(int a=0;a < ndim; a++) { 
	      for(int k=0; k< nelectrons(s);k++) 
		fref(i,j,a)+=inverse(det,s)(i,k)*
		moVal(jp,parent->occupation(det,s)(k),a+1);
	    }
	  }
	}
      }
      for(int i=0; i< tote; i++) { 
	int ir=i-spin(i)*nelectrons(0);
	for(int j=0; j< tote; j++) { 
	  int s=spin(j);
	  Array3 <doublevar> & fref(F(det,s));
	  int jr=j-s*nelectrons(0);
	  for(int a=0; a< ndim; a++) { 
	    for(int b=0; b< ndim; b++) { 
	      grad(i,a)+=fref(jr,jr,b)*coor_grad(j,i,a,b);
	    }
	  }
	}
	for(int a=0; a< ndim; a++) {
	  doublevar weight=1.0;
	  if(ndet >1) 
	    weight=parent->dkeeper.detwt(det)*detVal(det,0)*detVal(det,1);
	  gradmod2(i)+=grad(i,a)*grad(i,a)*weight;
	  gradlap(i,a)+=grad(i,a)*weight;
	}
      }
    }
  }//det
  
  for(int i=0; i< tote; i++) { 
    for(int a=0; a< ndim; a++) { 
      if(funcval==0)
	gradlap(i,a)=0;
      else if(ndet>1)
	gradlap(i,a)*=1/funcval;
    }
    if(funcval==0)
      gradmod2(i)=0;
    else if(ndet>1)
      gradmod2(i)*=1/funcval;
  }


  
  //Hessian lookup for MO's
  Array2 <int> hess_lookup(3,3);
  hess_lookup(0,0)=4;
  hess_lookup(1,1)=5;
  hess_lookup(2,2)=6;
  hess_lookup(0,1)=hess_lookup(1,0)=7;
  hess_lookup(0,2)=hess_lookup(2,0)=8;
  hess_lookup(1,2)=hess_lookup(2,1)=9;


  //Laplacian
  for(int i=0; i< tote; i++) { //coordinate whose lap we're taking
    for(int det=0;det<ndet;det++){
      doublevar lap=0;
      if(detVal(det,0)*detVal(det,1)!=0){
	for(int j=0; j< tote; j++) { 
	  int s=spin(j);
	  Array3 <doublevar> & fref(F(det,s));
	  int jr=j-s*nelectrons(0);
	  
	  for(int a=0; a< ndim; a++) { 
	    lap+=fref(jr,jr,a)*coor_lap(j,i,a);
	  }
	}
	//cout << "i= " << i << endl;
	//doublevar lap1=lap;
	//cout << "lap1 " << lap << endl;
	
	for(int s=0; s< 2; s++) { 
	  Array3 <doublevar> & fref(F(det,s));
	  for(int jr=0; jr < nelectrons(s); jr++) { 
	    int j=jr+s*nelectrons(0);
	    for(int kr=0; kr < nelectrons(s); kr++) { 
	      int k=kr+s*nelectrons(0);
	      for(int a=0; a< ndim; a++) { 
		for(int b=0; b< ndim; b++) { 
		  for(int g=0; g< ndim; g++) { 
		    lap-=coor_grad(j,i,a,b)*coor_grad(k,i,a,g)*fref(kr,jr,b)
		      *fref(jr,kr,g);
		    //if(fabs(coor_grad(i,j,a,b)) > 1e-16 
		    //   && fabs(coor_grad(i,k,a,g)) > 1e-16 ) 
		    //cout << "j " << j << " k " << k << " a " << a 
		    //     << " b " << b << " g " << g << " cijab " << coor_grad(i,j,a,b)
		    //     << " cikag " << coor_grad(i,k,a,g) 
		    //     << " fkjb " << fref(kr,jr,b) << " fjkg " << fref(jr,kr,g)
		    //     << endl;
		  }
		}
	      }
	    }
	  }
	}
	
	//doublevar lap2=lap-lap1;
	//cout << "lap2 " << lap2 << endl;
	
	int si=spin(i);
	for(int s=0; s< 2; s++) { 
	  Array2 <doublevar> & invref(inverse(det,s));
	  for(int jr=0; jr < nelectrons(s); jr++) { 
	    int j=jr+s*nelectrons(0);
	    for(int mr=0; mr < nelectrons(s); mr++) { 
	      int m=mr+s*nelectrons(0);
	      int occ=parent->occupation(det,s)(mr);
	      for(int a=0; a< ndim; a++) { 
		for(int b=0; b< ndim; b++) { 
		  for(int g=0; g< ndim; g++) { 
		    lap+=coor_grad(j,i,a,b)*coor_grad(j,i,a,g)
		      *invref(jr,mr)*moVal(j,occ,hess_lookup(b,g));
		    
		    //if(fabs(coor_grad(i,j,a,b)) > 1e-16 
		    //   && fabs(coor_grad(i,j,a,g)) > 1e-16 ) 
		    //  cout << "j " << j << " m " << m << " a " << a 
		    //       << " b " << b << " g " << g 
		    //       << " cijab " << coor_grad(i,j,a,b)
		    //       << " cijag " << coor_grad(i,j,a,g) 
		    //       << " invjm " << invref(jr,mr) 
		    //       << " hess "  << hess_lookup(b,g)
		    //       << " occ "   << occ 
		    //       << " moval " << moVal(jr,occ,hess_lookup(b,g)) 
		    //       << endl;
		  }
		}
	      }
	    }
	  }
	}
	//doublevar lap3=lap-lap1-lap2;
	//cout << "lap3 " << lap3 << endl;
      }
      if(ndet>1)
	lap*=parent->dkeeper.detwt(det)*detVal(det,0)*detVal(det,1);
      gradlap(i,3)+=lap;
    }//det
    
    if(funcval==0)
      gradlap(i,3)=0;
    else if(ndet>1)
      gradlap(i,3)*=1/funcval;
    
    gradlap(i,3)+=gradmod2(i);
    if(isnan(gradlap(i,3))) { 
      cout << "NAN! " << "gradmod2 " << gradmod2(i) <<endl;
      // << " lap1 " << lap1 << " lap2 " << lap2 
      //   << " lap3 " << lap3 << endl;
    }
  }//i-th electron

  //cout << "done " << endl;
  

}

//------------------------------------------------------------------------


/*!
*/

void Backflow_wf::getLap(Wavefunction_data * wfdata,
                     int e, Wf_return & lap)
{
  getVal(wfdata,e,lap);
  for(int d=0; d< ndim+1; d++) {
    lap.amp(0,d+1)=gradlap(e,d);
    lap.phase(0,d+1)=0;
  }

}

//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
