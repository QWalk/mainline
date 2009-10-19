//--------------------------------------------------------------------------
// src/Backflow_pf_wf.cpp
//
//
#include "Qmc_std.h"
#include "Backflow_pf_wf.h"
#include "MatrixAlgebra.h"
#include "Sample_point.h"
#include "Backflow_pf_wf_data.h"

//----------------------------------------------------------------------


void Backflow_pf_wf::generateStorage(Wavefunction_storage * & wfstore)
{
  wfstore=new  Backflow_pf_wf_storage;
  Backflow_pf_wf_storage * store;
  recast(wfstore, store);
  store->gradlap=gradlap;
  store->pfaffVal=pfaffVal;
}




//----------------------------------------------------------------------

void Backflow_pf_wf::init(Wavefunction_data * wfdata)
{

  Backflow_pf_wf_data * dataptr;
  recast(wfdata, dataptr);
  recast(wfdata, parent);

  nmo=dataptr->bfwrapper.nmo();
  npf=dataptr->pfkeeper.npf;
  coef_eps=dataptr->pfkeeper.coef_eps;
  //nsfunc=dataptr->pfkeeper.nsfunc;
  nelectrons.Resize(2);
  nelectrons=dataptr->pfkeeper.nelectrons;

  int tote=nelectrons(0)+nelectrons(1);
  ndim=3;

  int upuppairs=dataptr->pfkeeper.npairs(0);
  int downdownpairs=dataptr->pfkeeper.npairs(1);
  int nopairs=dataptr->pfkeeper.npairs(2);
  
  npairs=upuppairs+downdownpairs+nopairs;
    
  //Properties and intermediate calculation storage.
  moVal.Resize(10,tote, dataptr->pfkeeper.totoccupation(0).GetSize());
  updatedMoVal.Resize(dataptr->pfkeeper.totoccupation(0).GetSize(),10);
  
  inverse.Resize(npf);
  pfaffVal.Resize(npf);
  mopfaff_tot.Resize(npf);
  for (int pf=0;pf<npf;pf++){
    inverse(pf).Resize(npairs, npairs);
    inverse(pf)=0;
    
    for(int e=0; e< npairs; e++){
        inverse(pf)(e,e)=1;
    }
    pfaffVal(pf)=1;
  
    mopfaff_tot(pf).Resize(npairs,npairs);
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
void Backflow_pf_wf::notify(change_type change, int num)
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


void Backflow_pf_wf::saveUpdate(Sample_point * sample, int e,
                         Wavefunction_storage * wfstore)
{
  Backflow_pf_wf_storage * store;
  recast(wfstore, store);
  store->gradlap=gradlap;
  store->pfaffVal=pfaffVal;
	
}

//----------------------------------------------------------------------

void Backflow_pf_wf::restoreUpdate(Sample_point * sample, int e,
                            Wavefunction_storage * wfstore)
{
  Backflow_pf_wf_storage * store;
  recast(wfstore, store);
  gradlap=store->gradlap;
  pfaffVal=store->pfaffVal;
	  
  electronIsStaleLap=0;
  electronIsStaleVal=0;
  updateEverythingLap=0;
  updateEverythingVal=0;
}

//----------------------------------------------------------------------

void Backflow_pf_wf::updateVal(Wavefunction_data * wfdata,
                        Sample_point * sample)
{

  assert(sampleAttached);
  assert(dataAttached);

  Backflow_pf_wf_data * pfaffdata;
  recast(wfdata, pfaffdata);

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
  
  assert(sampleAttached);
  assert(dataAttached);
}

//----------------------------------------------------------------------

void Backflow_pf_wf::updateLap( Wavefunction_data * wfdata,
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


void Backflow_pf_wf::storeParmIndVal(Wavefunction_data * wfdata, Sample_point * sample,
                              int e, Array1 <doublevar> & vals )
{
}

//----------------------------------------------------------------------

void Backflow_pf_wf::getParmDepVal(Wavefunction_data * wfdata,
                            Sample_point * sample,
                            int e,
                            Array1 <doublevar> & oldval,
                            Wf_return & newval)
{
  updateVal(wfdata,sample);
  getVal(wfdata,e,newval);
}


//-----------------------------------------------------------------------

int Backflow_pf_wf::getParmDeriv(Wavefunction_data *  wfdata, 
			  Sample_point * sample ,
			  Parm_deriv_return & derivatives){
  return 0; 
}


//------------------------------------------------------------------------

void Backflow_pf_wf::updateVal(int e, Sample_point * sample) { 
  int nlist; Array1 <int> list;
  parent->bfwrapper.getNeighbors(sample,jast,e,list,nlist);
  cout << "nlist " << nlist << endl;
  calcVal(sample);
}

void Backflow_pf_wf::calcVal(Sample_point * sample)
{
  int tote=nelectrons(0)+nelectrons(1);

  for(int e=0; e< tote; e++) { 
    //int s=spin(e);
    sample->updateEIDist();
    parent->bfwrapper.updateVal(sample,jast, e, 0, updatedMoVal);
    for(int i=0; i< updatedMoVal.GetDim(0); i++) {
      for(int d=0; d< 1; d++) { 
        moVal(d,e,i)=updatedMoVal(i,d);
      }
    }

  } 

  //cout << "doing pfaffian " << endl;
  for (int pf=0;pf<npf;pf++){
    FillPfaffianMatrix( mopfaff_tot(pf),
                        moVal, 
                        parent->pfkeeper.occupation_pos,
                        parent->pfkeeper.npairs,
                        parent->pfkeeper.order_in_pfaffian(pf),
                        parent->pfkeeper.tripletorbuu, 
                        parent->pfkeeper.tripletorbdd,
                        parent->pfkeeper.singletorb,
                        parent->pfkeeper.unpairedorb,
                        parent->pfkeeper.normalization,
                        coef_eps
                        );
    pfaffVal(pf) = PfaffianInverseMatrix(mopfaff_tot(pf), inverse(pf));
    //cout << "Pfaffian value:             "<< pfaffVal(pf) << endl;
  }
}

//------------------------------------------------------------------------

void Backflow_pf_wf::getVal(Wavefunction_data * wfdata, int e,
                     Wf_return & val)
{

  Array1 <doublevar> si(1, 0.0);
  Array2 <doublevar> vals(1,1,0.0);

  assert(val.amp.GetDim(0) >=1);
  assert(val.amp.GetDim(1) >= 1);


  doublevar funcval=0;
  for (int pf=0;pf<npf;pf++)
    funcval += parent->pfkeeper.pfwt(pf)*pfaffVal(pf);

  if(fabs(funcval) > 0)
    vals(0,0)=log(fabs(funcval));
  else 
    vals(0,0)=-1e3;

  si(0)=sign(funcval);

  val.setVals(vals, si);
}

//-----------------------------------------------------------------------
void Backflow_pf_wf::getSymmetricVal(Wavefunction_data * wfdata,
		     int e, Wf_return & val){

  Array1 <doublevar> si(1, 0.0);
  Array2 <doublevar> vals(1,1,0.0);
  val.setVals(vals, si);
} 


//----------------------------------------------------------------------

void Backflow_pf_wf::getDensity(Wavefunction_data * wfdata, int e,
                         Array2 <doublevar> & dens)
{
  error("No density for Backflow now..");

}

//----------------------------------------------------------------------------

//Watch out that some of the indices are reversed from Kwon--
//this is to make the updates a bit more convienent.
void Backflow_pf_wf::calcLap(Sample_point * sample)
{

  //cout <<"Start Backflow_pf_wf::calcLap"<<endl;
  int tote=nelectrons(0)+nelectrons(1);
  Array3 <doublevar> temp_der;
  Array2 <doublevar> temp_lap;
  Array3 <doublevar> moVal_gradonly;
  moVal_gradonly.Resize(3,tote,updatedMoVal.GetDim(0)); 
  for(int e=0; e< tote; e++) { 
    //int s=spin(e);
    sample->updateEIDist();
    parent->bfwrapper.updateLap(sample,jast,e,0,updatedMoVal,
				temp_der, temp_lap);
    
    for(int i=0; i< updatedMoVal.GetDim(0); i++) {
      for(int d=0; d< 10; d++) { 
        moVal(d,e,i)=updatedMoVal(i,d);
	if(d>0 && d<4)
	  moVal_gradonly(d-1,e,i)=moVal(d,e,i);
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
     
  //cout << "doing pfaffian " << endl;
  for (int pf=0;pf<npf;pf++){
    FillPfaffianMatrix( mopfaff_tot(pf),
                        moVal, 
                        parent->pfkeeper.occupation_pos,
                        parent->pfkeeper.npairs,
                        parent->pfkeeper.order_in_pfaffian(pf),
                        parent->pfkeeper.tripletorbuu, 
                        parent->pfkeeper.tripletorbdd,
                        parent->pfkeeper.singletorb,
                        parent->pfkeeper.unpairedorb,
                        parent->pfkeeper.normalization,
                        coef_eps
                        );
    pfaffVal(pf) = PfaffianInverseMatrix(mopfaff_tot(pf), inverse(pf));
    //cout << "Pfaffian value:             "<< pfaffVal(pf) << endl;
  }
  //cout << "done pfaffian" << endl;


  gradlap=0;
  Array1 <doublevar> gradmod2(tote);
  Array1 <Array3 <doublevar> > F(npf); 
  Array1 <Array4 <doublevar> > H(npf);
  doublevar funcval=0;
  gradmod2=0;
  Array2 <doublevar> grad(tote,ndim);
  Array1 < Array1 < Array1 < Array1 <doublevar> >  > > mopfaff_row;
  Array1 < Array1 < Array1 < Array2 <doublevar> >  > > mopfaff_row_hess;
  mopfaff_row.Resize(npf);
  mopfaff_row_hess.Resize(npf);

  //cout <<"Getting all pairing orbital derivatives"<<endl;
  for(int pf=0;pf<npf;pf++){
    mopfaff_row(pf).Resize(tote);
    mopfaff_row_hess(pf).Resize(tote);
    for(int i=0; i<tote; i++) {
      UpdatePfaffianRowLap(mopfaff_row(pf)(i), 
			   i,  
			   moVal, 
			   parent->pfkeeper.occupation_pos,
			   parent->pfkeeper.npairs,
			   parent->pfkeeper.order_in_pfaffian(pf),
			   parent->pfkeeper.tripletorbuu, 
			   parent->pfkeeper.tripletorbdd,
			   parent->pfkeeper.singletorb,
			   parent->pfkeeper.unpairedorb,
			   parent->pfkeeper.normalization,
                           coef_eps
			   );
      
      UpdatePfaffianRowHess(mopfaff_row_hess(pf)(i), 
			   i,  
			   moVal_gradonly, 
			   parent->pfkeeper.occupation_pos,
			   parent->pfkeeper.npairs,
			   parent->pfkeeper.order_in_pfaffian(pf),
			   parent->pfkeeper.tripletorbuu, 
			   parent->pfkeeper.tripletorbdd,
			   parent->pfkeeper.singletorb,
			   parent->pfkeeper.unpairedorb,
                           parent->pfkeeper.normalization,
                           coef_eps
			   );	

      //  for(int k=0;k<mopfaff_row(pf)(i).GetSize();k++){
      //     cout<<"i and k "<<i<<" "<<k<<endl;
      //     for(int a=0;a<mopfaff_row(pf)(i)(k).GetDim(0);a++){
      //  for(int b=0;b<mopfaff_row_hess(pf)(i)(k).GetDim(0);b++)
      //       cout << mopfaff_row(pf)(i)(k)(a) <<"  ";
      //       cout <<endl;
      //     }
      //}
    }
  }
  //cout <<"Done getting all pairing orbital derivatives"<<endl;

  for(int pf=0;pf<npf;pf++){
    grad=0;
    //cout << "pfaffVal "<<pfaffVal(det,0)<<"  "<<pfaffVal(det,1)<<endl;
    funcval+=parent->pfkeeper.pfwt(pf)*pfaffVal(pf);
    F(pf).Resize(tote,tote,ndim);
    F(pf)=0;
    H(pf).Resize(tote,tote,ndim,ndim);
    H(pf)=0;

    if(pfaffVal(pf)!=0){
      Array3 <doublevar> & fref(F(pf));
      Array4 <doublevar> & href(H(pf));

      for(int i=0; i<tote; i++)
	for(int j=0; j< tote; j++) {
	  for(int k=0; k< npairs; k++){
	    for(int a=0;a < ndim; a++) { 
	      fref(i,j,a)+=mopfaff_row(pf)(j)(k)(a+1)*inverse(pf)(k,i);
	      for(int l=0; l< npairs; l++){
		for(int b=0;b < ndim; b++) { 
		  href(i,j,a,b)+=mopfaff_row(pf)(i)(k)(a+1)*mopfaff_row(pf)(j)(l)(b+1)*inverse(pf)(k,l);
		}
	      }
	    }
	  }
	}

      //for(int i=0; i<tote; i++)
      //	for(int j=0; j< tote; j++) {
      //  cout <<i<<" "<<j<<endl;
      //  for(int a=0;a < ndim; a++)
      //    for(int b=0;b < ndim; b++) 
      //      cout <<href(i,j,a,b)<<endl;
      //}
      
      for(int i=0; i< tote; i++) { 
	for(int j=0; j< tote; j++) { 
	  for(int a=0; a< ndim; a++) { 
	    for(int b=0; b< ndim; b++) { 
	      grad(i,a)+=fref(j,j,b)*coor_grad(j,i,a,b);
	    }
	  }
	}
	for(int a=0; a< ndim; a++) {
	  doublevar weight=1.0;
	  if(npf >1) 
	    weight=parent->pfkeeper.pfwt(pf)*pfaffVal(pf);
	  gradmod2(i)+=grad(i,a)*grad(i,a)*weight;
	  gradlap(i,a)+=grad(i,a)*weight;
	}
      }
    }  
  }//pf
  
  for(int i=0; i< tote; i++) { 
    for(int a=0; a< ndim; a++) { 
      if(funcval==0)
	gradlap(i,a)=0;
      else if(npf>1)
	gradlap(i,a)*=1/funcval;
    }
    if(funcval==0)
      gradmod2(i)=0;
    else if(npf>1)
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
    for(int pf=0;pf<npf;pf++){
      Array3 <doublevar> & fref(F(pf));
      Array4 <doublevar> & href(H(pf));
      Array2 <doublevar> & invref(inverse(pf)); 
      
      doublevar lap=0;
      if(pfaffVal(pf)!=0){
	for(int j=0; j< tote; j++) { 
	  for(int a=0; a< ndim; a++) { 
	    lap+=fref(j,j,a)*coor_lap(j,i,a);
	  }
	}
	//cout << "i= " << i << endl;
	//doublevar lap1=lap;
	//cout << "lap1 " << lap << endl;
	
	for(int j=0; j < tote; j++) { 
	  for(int k=0; k < npairs; k++) {
	    for(int a=0; a< ndim; a++) { 
	      for(int b=0; b< ndim; b++) { 
		for(int g=0; g< ndim; g++) { 
		  if(k<tote){
		    lap-=coor_grad(j,i,a,b)*coor_grad(k,i,a,g)*
		      (fref(k,j,b)*fref(j,k,g)-
		       (href(j,k,b,g)+mopfaff_row_hess(pf)(j)(k)(b,g))*invref(k,j));
		  }
		  lap+=coor_grad(j,i,a,b)*coor_grad(j,i,a,g)
		    *invref(k,j)*mopfaff_row(pf)(j)(k)(hess_lookup(b,g));
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
      if(npf>1)
	lap*=parent->pfkeeper.pfwt(pf)*pfaffVal(pf);
      gradlap(i,3)+=lap;
    }//pf

    //cout <<"gradlap(i,3) "<<gradlap(i,3)<<endl;
    
    if(funcval==0)
      gradlap(i,3)=0;
    else if(npf>1)
      gradlap(i,3)*=1/funcval;
    
    gradlap(i,3)+=gradmod2(i);
    //cout <<"After grad correction: gradlap(i,3) "<<gradlap(i,3)<<endl;
    if(isnan(gradlap(i,3))) { 
      cout << "NAN! " << "gradmod2 " << gradmod2(i) <<endl;
      //<< " lap1 " << lap1 << " lap2 " << lap2 
      //<< " lap3 " << lap3 << endl;
    }
  }//i-th electron

  //cout << "done " << endl;
  

}

//------------------------------------------------------------------------


/*!
*/

void Backflow_pf_wf::getLap(Wavefunction_data * wfdata,
                     int e, Wf_return & lap)
{
  getVal(wfdata,e,lap);
  for(int d=0; d< ndim+1; d++) {
    lap.amp(0,d+1)=gradlap(e,d);
    //cout <<"lap.amp(0,d+1) "<<lap.amp(0,d+1)<<endl;
    lap.phase(0,d+1)=0;
  }
}

//-------------------------------------------------------------------------

