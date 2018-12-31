//--------------------------------------------------------------------------
// src/Pfaff_wf_calc.cpp
//
//
#include "Qmc_std.h"
#include "Pfaff_wf.h"
#include "MatrixAlgebra.h"
#include "Sample_point.h"
#include "Pfaff_wf_data.h"





//----------------------------------------------------------------------

void Pfaff_wf::init(Wavefunction_data * wfdata)
{

  recast(wfdata, parent);
  nmo=parent->nmo;
  npf=parent->npf;
  nsfunc=parent->nsfunc;
  coef_eps=parent->coef_eps;
  nelectrons.Resize(2);
  nelectrons=parent->nelectrons;
  int upuppairs=parent->npairs(0);
  int downdownpairs=parent->npairs(1);
  int nopairs=parent->npairs(2);
  
  npairs=upuppairs+downdownpairs+nopairs;
  int tote=nelectrons(0)+nelectrons(1);
  ndim=3;

  
  //Properties and intermediate calculation storage.

  moVal.Resize(5, tote, parent->totoccupation(0).GetSize());
  updatedMoVal.Resize(parent->totoccupation(0).GetSize(),5);
  
  inverse.Resize(npf);
  pfaffVal.Resize(npf);
  mopfaff_tot.Resize(npf);
  for (int pf=0;pf<npf;pf++){
    inverse(pf).Resize(upuppairs+downdownpairs+nopairs, upuppairs+downdownpairs+nopairs);
    inverse(pf)=0;
    
    for(int e=0; e< upuppairs+downdownpairs+nopairs; e++){
        inverse(pf)(e,e)=1;
    }
    pfaffVal(pf)=1;
  
    mopfaff_tot(pf).Resize(upuppairs+downdownpairs+nopairs, 
                       upuppairs+downdownpairs+nopairs);
  }

  electronIsStaleVal.Resize(tote);
  electronIsStaleLap.Resize(tote);

  electronIsStaleVal=0;
  electronIsStaleLap=0;
  updateEverythingVal=1;
  updateEverythingLap=1;
  updateMoLap=1;
  sampleAttached=0;
  dataAttached=0;
  firstime=true;
  


}




//----------------------------------------------------------------------

/*!
 */
void Pfaff_wf::notify(change_type change, int num)
{
  //  cout<<"notify\n";
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
      updateMoLap=0;
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

void Pfaff_wf::updateVal(Sample_point * sample)
{
  // cout << "updateVal-main part \n";
  
  assert(sampleAttached);
  assert(dataAttached);

  if(parent->optimize_pf ) {
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

  //cout << "This is end of updateVal-main part"<<endl;;

}

//----------------------------------------------------------------------


void Pfaff_wf::updateLap(Sample_point * sample)
{
  //Array1 <doublevar> elecpos(3);
  assert(sampleAttached);
  assert(dataAttached);

  //  cout << "updateLap-main part \n";

  if(updateEverythingLap==1){
    //cout <<"calcLap" <<endl;
    calcLap(sample);
    updateEverythingVal=0;
    updateEverythingLap=0;
    electronIsStaleLap=0;
    electronIsStaleVal=0;
  }
  else {
    for(int e=0; e< nelectrons(0)+nelectrons(1); e++){
      // sample->getElectronPos(e, elecpos);
      //cout <<  "elecpos:  "<< elecpos(0)<<"  "<< elecpos(1)<<"  "<< elecpos(2)<<endl;
      
      
      //  cout << "Checking electron: "<<e<<endl;
      if(electronIsStaleLap(e)){
        //cout <<"updateLap1" <<endl;
        updateLap(sample, e);
        electronIsStaleLap(e)=0;
        electronIsStaleVal(e)=0;
      }
    }
    
  }
 // cout <<"pfaffVal inside UpdateLap" <<pfaffVal<<endl;
 //cout << "This is end of updateLap-main part"<<endl;;
}

//----------------------------------------------------------------------

doublevar TripletPairFunction(int i, int j, 
                              const Array3 <doublevar> & moVal, 
                              const Array1 <int> & occupation_pos,
                              const Array1 <doublevar> & tripletorb,
                              const doublevar & norm,
			      doublevar coef_eps
                              )
{

  if (norm<coef_eps )
    return 0.0;
  else {
    doublevar temp;
    if (i==j) { return 0.0; }
    else {
      temp=0.0;
      int count=0;
      int d=0;
      int dim=occupation_pos.GetSize();
      //cout <<"dim "<<dim<<endl;
      for (int k=0;k<dim;k++)
	for (int l=k+1;l<dim;l++){
	  if (fabs(tripletorb(count))>coef_eps){
	    temp+=tripletorb(count)*
	          (moVal(d,i,occupation_pos(k))*moVal(d,j,occupation_pos(l))
		   -moVal(d,i,occupation_pos(l))*moVal(d,j,occupation_pos(k)));
	  }
	  count++;
	}
      return 0.70710678118654752440*temp/norm; 
    }
  }
}

//----------------------------------------------------------------------


doublevar SingletPairFunction(int i, int j, 
                              const Array3 <doublevar> & moVal,
                              const Array1 <int> & occupation_pos,
                              const Array1 <doublevar> & singletorb,
                              const doublevar & norm,
			      doublevar coef_eps
                              )
{
  if (norm<coef_eps )
    return 0.0;
  else {
    //cout << "singletpairfunc\n";
    doublevar temp=0.0;
    int count=0;
    int d=0;
    int dim=occupation_pos.GetSize();
    for (int k=0;k<dim;k++)
      for (int l=k;l<dim;l++){
	if (fabs(singletorb(count))>coef_eps){
	  if (k!=l)
	    temp+=0.70710678118654752440*
	      singletorb(count)*
	      (moVal(d,i,occupation_pos(k))*moVal(d,j,occupation_pos(l))
	       +moVal(d,i,occupation_pos(l))*moVal(d,j,occupation_pos(k)));
	  else
	    temp+=singletorb(count)*moVal(d,i,occupation_pos(k))*moVal(d,j,occupation_pos(l));
	}
	count++;
      }
    return temp/norm; 
  }
}

doublevar UnpairedOrbFunction(int i, int l,  
                              const Array3 <doublevar> & moVal, 
                              const Array1 <int> & occupation_pos,
                              const Array1 <doublevar> & unpairedorb,
                              const doublevar & norm,
			      doublevar coef_eps
                              )
{ 
  if (norm<coef_eps )
    return 0.0;
  else {
    int lenght=occupation_pos.GetSize();
    doublevar temp=0.0;
    for (int k=0;k<lenght;k++){
      //cout << "Index: "<<k+l*lenght<< "k "<< k<< "l "<<l<<"lenght " <<lenght;
      if (fabs(unpairedorb(k+l*lenght))>coef_eps)
	temp+=unpairedorb(k+l*lenght)*(moVal(0,i,occupation_pos(k)));
    }
    //cout << endl;
    return temp/norm; 
  }
}





void TripletPairFunctionLap(int i, int j, 
                            const Array3 <doublevar> & moVal, 
                            const Array1 <int> & occupation_pos,
                            const Array1 <doublevar> & tripletorb,
                            const doublevar & norm,
                            Array1 <doublevar> & func_lap,
			    doublevar coef_eps
                            )
  //note: differentiation is with respect to i-th electron only!!!
{
  int dmax=moVal.GetDim(0);
  Array1 <doublevar> temp(dmax);
  if (norm<coef_eps )
    func_lap=0.0;
  else {
    func_lap(0)=0.0;
    if (i==j) {
      func_lap=0;
    }
    else {
      int count=0;
      temp=0.0;
      int dim=occupation_pos.GetSize();
      for (int k=0;k<dim;k++)
	for (int l=k+1;l<dim;l++){    
	  if (abs(tripletorb(count))>coef_eps)
	    for (int d=1;d<dmax;d++){
	      temp(d)+=tripletorb(count)*
		(moVal(d,i,occupation_pos(k))*moVal(0,j,occupation_pos(l))-
                 moVal(d,i,occupation_pos(l))*moVal(0,j,occupation_pos(k)));
	      //action of operator only on i-th electron
	    }
	  count++;
	}
      for (int d=1;d<dmax;d++)
	func_lap(d)=0.70710678118654752440*temp(d)/norm; 
    }
  }
}

void SingletPairFunctionLap(int i, int j, 
                            const Array3 <doublevar> & moVal, 
                            const Array1 <int> & occupation_pos,
                            const Array1 <doublevar> & singletorb,
                            const doublevar & norm,
                            Array1 <doublevar> & func_lap,
			    doublevar coef_eps
                            )

  //note: differentiation is with respect to i-th electron only!!!
{ 
  if (norm<coef_eps )
    func_lap=0.0;
  else {
    int dmax=moVal.GetDim(0);
    // cout <<"SingletpairfuncLap\n";
    Array1 <doublevar> temp(dmax);
    func_lap(0)=0.0;
    temp=0.0;
    int count=0;
    int dim=occupation_pos.GetSize();
    for (int k=0;k<dim;k++)
      for (int l=k;l<dim;l++){
	if (abs(singletorb(count))>coef_eps){
	  if (k!=l)
	    for (int d=1;d<dmax;d++)
	      temp(d)+=0.70710678118654752440*
		singletorb(count)*
		(moVal(d,i,occupation_pos(k))*moVal(0,j,occupation_pos(l))+
		 moVal(d,i,occupation_pos(l))*moVal(0,j,occupation_pos(k)));
	  //action of operator only on i-th electron
	  else 
	    for (int d=1;d<dmax;d++)
	      temp(d)+=singletorb(count)*moVal(d,i,occupation_pos(k))*moVal(0,j,occupation_pos(l));
	}
	count++;
      }
    
    for (int d=1;d<dmax;d++)
      func_lap(d)=temp(d)/norm;
  }

}


void UnpairedOrbFunctionLap(int i, int l,  
                            const Array3 <doublevar> &  moVal, 
                            const Array1 <int> & occupation_pos,
                            const Array1 <doublevar> & unpairedorb,
                            const doublevar & norm,
                            Array1 <doublevar> & func_lap,
			    doublevar coef_eps
                            )
  //note: differentiation is with respect to i-th electron only!!!
  //second index (l) is excited  orbital number
{ 
  if (norm<coef_eps )
    func_lap=0.0;
  else {
    int dmax=moVal.GetDim(0);
    Array1 <doublevar> temp(dmax);
    int lenght=occupation_pos.GetSize();
    func_lap(0)=0.0;
    temp=0.0;
    for (int k=0;k<lenght;k++){
      if (abs(unpairedorb(k+l*lenght))>coef_eps)
        //cout << "Index1: "<<k+l*lenght<<endl;
	for (int d=1;d<dmax;d++)
	  temp(d)+=unpairedorb(k+l*lenght)*moVal(d,i,occupation_pos(k));
    }
    
    for (int d=1;d<dmax;d++)
      func_lap(d)=temp(d)/norm; 
  }
}


void TripletPairFunctionHess(int i, int j, 
			     const Array3 <doublevar> & moVal, 
			     const Array1 <int> & occupation_pos,
			     const Array1 <doublevar> & tripletorb,
			     const doublevar & norm,
			     Array2 <doublevar> & func_lap,
			     doublevar coef_eps
                            )
{
  int dmax=moVal.GetDim(0);
  Array2 <doublevar> temp(dmax,dmax);
  if (norm<coef_eps )
    func_lap=0.0;
  else {
    func_lap(0)=0.0;
    if (i==j) {
      func_lap=0;
    }
    else {
      int count=0;
      temp=0.0;
      int dim=occupation_pos.GetSize();
      for (int k=0;k<dim;k++)
	for (int l=k+1;l<dim;l++){    
	  if (abs(tripletorb(count))>coef_eps)
	    for (int a=0;a<dmax;a++)
	      for (int b=0;b<dmax;b++){
	      temp(a,b)+=tripletorb(count)*
		(moVal(a,i,occupation_pos(k))*moVal(b,j,occupation_pos(l))-
                 moVal(a,i,occupation_pos(l))*moVal(b,j,occupation_pos(k)));
	      //action of operator only on i-th electron
	    }
	  count++;
	}
      for (int a=0;a<dmax;a++)
	for (int b=0;b<dmax;b++)
	  func_lap(a,b)=0.70710678118654752440*temp(a,b)/norm; 
    }
  }
}

void SingletPairFunctionHess(int i, int j, 
			     const Array3 <doublevar> & moVal, 
			     const Array1 <int> & occupation_pos,
			     const Array1 <doublevar> & singletorb,
			     const doublevar & norm,
			     Array2 <doublevar> & func_lap,
			     doublevar coef_eps
                            )

{ 
  if ( norm<coef_eps )
    func_lap=0.0;
  else {
    int dmax=moVal.GetDim(0);
    // cout <<"SingletpairfuncLap\n";
    Array2 <doublevar> temp(dmax,dmax);
    func_lap(0)=0.0;
    temp=0.0;
    int count=0;
    int dim=occupation_pos.GetSize();
    for (int k=0;k<dim;k++)
      for (int l=k;l<dim;l++){
	if (abs(singletorb(count))>coef_eps){
	  if (k!=l)
	    for (int a=0;a<dmax;a++)
	      for (int b=0;b<dmax;b++)
		temp(a,b)+=0.70710678118654752440*singletorb(count)*
		  (moVal(a,i,occupation_pos(k))*moVal(b,j,occupation_pos(l))+
		   moVal(a,i,occupation_pos(l))*moVal(b,j,occupation_pos(k)));
	  //action of operator only on i-th electron
	  else 
	    for (int a=0;a<dmax;a++)
	      for (int b=0;b<dmax;b++)
		temp(a,b)+=singletorb(count)*moVal(a,i,occupation_pos(k))*moVal(b,j,occupation_pos(l));
	}
	count++;
      }
    for (int a=0;a<dmax;a++)
      for (int b=0;b<dmax;b++)
	func_lap(a,b)=temp(a,b)/norm;
  }
}


void UnpairedOrbFunctionHess(int i, int l,  
			     const Array3 <doublevar> &  moVal, 
			     const Array1 <int> & occupation_pos,
			     const Array1 <doublevar> & unpairedorb,
			     const doublevar & norm,
			     Array2 <doublevar> & func_lap,
			     doublevar coef_eps
                            )
{ 
  func_lap=0.0;
}


//----------------------------------------------------------------------

doublevar TripletPairFunctionGrad(int i, int j, 
				  const Array3 <doublevar> & moVal, 
				  const Array1 <int> & occupation_pos,
				  const Array1 <doublevar> & tripletorb,
				  const doublevar & norm,
				  int coef,
				  doublevar coef_eps
				  )
{
  if ( norm<coef_eps )
    return 0.0;
  else {
    doublevar temp;
    if (i==j) { return 0.0; }
    else {
      doublevar norm2=norm*norm;
      int count=0;
      int d=0;
      int dim=occupation_pos.GetSize();
      temp=0.0;
      for (int k=0;k<dim;k++)
	for (int l=k+1;l<dim;l++){
          if (abs(tripletorb(count))>coef_eps){
            temp-=(tripletorb(coef)/norm2)*tripletorb(count)*
              (moVal(d,i,occupation_pos(k))*moVal(d,j,occupation_pos(l))
               -moVal(d,i,occupation_pos(l))*moVal(d,j,occupation_pos(k)));
          }
	  if(count==coef)
            temp+=moVal(d,i,occupation_pos(k))*moVal(d,j,occupation_pos(l))
              -moVal(d,i,occupation_pos(l))*moVal(d,j,occupation_pos(k));
	  count++;
	}
      return 0.70710678118654752440*temp/norm; 
    }
  }
}

//----------------------------------------------------------------------


doublevar SingletPairFunctionGrad(int i, int j, 
				  const Array3 <doublevar> & moVal,
				  const Array1 <int> & occupation_pos,
				  const Array1 <doublevar> & singletorb,
				  const doublevar & norm,
				  int coef,
				  doublevar coef_eps
				  )
{
  if (norm<coef_eps )
    return 0.0;
  else {
    doublevar temp=0.0;
    int count=0;
    int d=0;
    int dim=occupation_pos.GetSize();
    doublevar norm2=norm*norm;
    for (int k=0;k<dim;k++)
      for (int l=k;l<dim;l++){
        if (abs(singletorb(count))>coef_eps){
          if (k!=l)
	    temp-=0.70710678118654752440*
	      (singletorb(coef)/norm2)*singletorb(count)*
	      (moVal(d,i,occupation_pos(k))*moVal(d,j,occupation_pos(l))
	       +moVal(d,i,occupation_pos(l))*moVal(d,j,occupation_pos(k)));
	  else
	    temp-=(singletorb(coef)/norm2)*singletorb(count)*moVal(d,i,occupation_pos(k))*moVal(d,j,occupation_pos(l));
        }
        if(count==coef){
          if (k!=l)
            temp+=0.70710678118654752440*
              (moVal(d,i,occupation_pos(k))*moVal(d,j,occupation_pos(l))
               +moVal(d,i,occupation_pos(l))*moVal(d,j,occupation_pos(k)));
          else{
            temp+=moVal(d,i,occupation_pos(k))*moVal(d,j,occupation_pos(l));
          }
        }
        count++;
      }
    return temp/norm; 
  }
}
//----------------------------------------------------------------------

doublevar UnpairedOrbFunctionGrad(int i, int l,  
				  const Array3 <doublevar> & moVal, 
				  const Array1 <int> & occupation_pos,
				  const Array1 <doublevar> & unpairedorb,
				  const doublevar & norm,
				  int coef,
				  doublevar coef_eps
				  )
{ 
  if ( norm<coef_eps )
    return 0.0;
  else {
    int lenght=occupation_pos.GetSize();
    doublevar norm2=norm*norm;
    doublevar temp=0.0;
    for (int k=0;k<lenght;k++){
      if (abs(unpairedorb(k+l*lenght))>coef_eps)
        temp-=(unpairedorb(coef)/norm2)*unpairedorb(k+l*lenght)*(moVal(0,i,occupation_pos(k)));
      if(coef==k+l*lenght)
        temp+=(moVal(0,i,occupation_pos(k)));
    }
    //cout << endl;
    return temp/norm; 
  }
}

//----------------------------------------------------------------------


doublevar TripletPairFunctionHess(int i, int j, 
				  const Array3 <doublevar> & moVal, 
				  const Array1 <int> & occupation_pos,
				  const Array1 <doublevar> & tripletorb,
				  const doublevar & norm,
				  int coef1,
				  int coef2,
				  doublevar coef_eps
				  )
{
  if ( norm<coef_eps )
    return 0.0;
  else {
    doublevar temp;
    doublevar norm2=norm*norm;
    doublevar norm4=norm2*norm2;
    if (i==j) { return 0.0; }
    else {
      int count=0;
      int d=0;
      int dim=occupation_pos.GetSize();
      temp=0.0;
      for (int k=0;k<dim;k++)
	for (int l=k+1;l<dim;l++){
          if (abs(tripletorb(count))>coef_eps){
            if(coef1==coef2)
              temp+=(3.0*tripletorb(coef1)*tripletorb(coef2)/norm4-1)*tripletorb(count)*
                (moVal(d,i,occupation_pos(k))*moVal(d,j,occupation_pos(l))
                 -moVal(d,i,occupation_pos(l))*moVal(d,j,occupation_pos(k)));
            else
              temp+=(3.0*tripletorb(coef1)*tripletorb(coef2)/norm4)*tripletorb(count)*
                (moVal(d,i,occupation_pos(k))*moVal(d,j,occupation_pos(l))
                 -moVal(d,i,occupation_pos(l))*moVal(d,j,occupation_pos(k)));
          }
	  if(count==coef1)
	   temp-=(tripletorb(coef2)/norm2)*(moVal(d,i,occupation_pos(k))*moVal(d,j,occupation_pos(l))
	                                 -moVal(d,i,occupation_pos(l))*moVal(d,j,occupation_pos(k)));
	  if(count==coef2)
	   temp-=(tripletorb(coef1)/norm2)*(moVal(d,i,occupation_pos(k))*moVal(d,j,occupation_pos(l))
	                                 -moVal(d,i,occupation_pos(l))*moVal(d,j,occupation_pos(k)));
	  
	  count++;
	}
      return 0.70710678118654752440*temp/norm; 
    }
  }
}

//----------------------------------------------------------------------


doublevar SingletPairFunctionHess(int i, int j, 
				  const Array3 <doublevar> & moVal,
				  const Array1 <int> & occupation_pos,
				  const Array1 <doublevar> & singletorb,
				  const doublevar & norm,
				  int coef1,
				  int coef2,
				  doublevar coef_eps
				  )
{
  if ( norm<coef_eps )
    return 0.0;
  else {
    doublevar temp=0.0;
    doublevar norm2=norm*norm;
    doublevar norm4=norm2*norm2;
    int count=0;
    int d=0;
    int dim=occupation_pos.GetSize();
    for (int k=0;k<dim;k++)
      for (int l=k;l<dim;l++){
        if (abs(singletorb(count))>coef_eps){
	  if (k!=l){
	    if(coef1==coef2)
	      temp+=0.70710678118654752440*
		(3.0*singletorb(coef1)*singletorb(coef2)/norm4-1)*singletorb(count)*
		(moVal(d,i,occupation_pos(k))*moVal(d,j,occupation_pos(l))
		 +moVal(d,i,occupation_pos(l))*moVal(d,j,occupation_pos(k)));
	    else
	      temp+=0.70710678118654752440*
		(3.0*singletorb(coef1)*singletorb(coef2)/norm4)*singletorb(count)*
		(moVal(d,i,occupation_pos(k))*moVal(d,j,occupation_pos(l))
		 +moVal(d,i,occupation_pos(l))*moVal(d,j,occupation_pos(k)));
	  }
	  else{
	    if(coef1==coef2)
	      temp+=(3.0*singletorb(coef1)*singletorb(coef2)/norm4-1)
		*singletorb(count)*moVal(d,i,occupation_pos(k))*moVal(d,j,occupation_pos(l));
	    else
	      temp+=(3.0*singletorb(coef1)*singletorb(coef2)/norm4)
		*singletorb(count)*moVal(d,i,occupation_pos(k))*moVal(d,j,occupation_pos(l));
	  }
        }
        if(count==coef1){
          if (k!=l)
            temp-=0.70710678118654752440*(singletorb(coef2)/norm2)*
              (moVal(d,i,occupation_pos(k))*moVal(d,j,occupation_pos(l))
               +moVal(d,i,occupation_pos(l))*moVal(d,j,occupation_pos(k)));
          else
            temp-=(singletorb(coef2)/norm2)*moVal(d,i,occupation_pos(k))*moVal(d,j,occupation_pos(l));
        }
        if(count==coef2){
          if (k!=l)
            temp-=0.70710678118654752440*(singletorb(coef1)/norm2)*
              (moVal(d,i,occupation_pos(k))*moVal(d,j,occupation_pos(l))
               +moVal(d,i,occupation_pos(l))*moVal(d,j,occupation_pos(k)));
          else
            temp-=(singletorb(coef1)/norm2)*moVal(d,i,occupation_pos(k))*moVal(d,j,occupation_pos(l));
        }
	count++;
      }
    return temp/norm; 
  }
}
//----------------------------------------------------------------------

doublevar UnpairedOrbFunctionHess(int i, int l,  
				  const Array3 <doublevar> & moVal, 
				  const Array1 <int> & occupation_pos,
				  const Array1 <doublevar> & unpairedorb,
				  const doublevar & norm,
				  int coef1,
				  int coef2,
				  doublevar coef_eps
				  )
{ 
  
  if ( norm<coef_eps )
    return 0.0;
  else {
    
    int lenght=occupation_pos.GetSize();
    doublevar temp=0.0;
    doublevar norm2=norm*norm;
    doublevar norm4=norm2*norm2;
    for (int k=0;k<lenght;k++){
      if (abs(unpairedorb(k+l*lenght))>coef_eps){
        if(coef1==coef2){
          temp+=(3.0*unpairedorb(coef1)*unpairedorb(coef2)/norm4-1)*unpairedorb(k+l*lenght)*(moVal(0,i,occupation_pos(k)));
        }
        else{
          temp+=(3.0*unpairedorb(coef1)*unpairedorb(coef2)/norm4)*unpairedorb(k+l*lenght)*(moVal(0,i,occupation_pos(k)));
        }
      }
      if(coef1==k+l*lenght)
	temp-=(unpairedorb(coef2)/norm2)*(moVal(0,i,occupation_pos(k)));
      if(coef2==k+l*lenght)
	temp-=(unpairedorb(coef1)/norm2)*(moVal(0,i,occupation_pos(k)));
     
    }
    //cout << endl;
    return temp/norm; 
  }
  
}

//----------------------------------------------------------------------




void UpdatePfaffianRowVal(Array1 <doublevar> & mopfaff_row, 
                          int & e,  
                          const Array3 <doublevar> &  moVal, 
                          const Array1 < Array1 <int> > & occupation_pos,
                          const Array1 <int> & npairs, 
                          const Array2 <int> & order_in_pfaffian,
                          const Array1 < Array1 <doublevar> > & tripletorbuu,
			  const Array1 < Array1 <doublevar> > & tripletorbdd,
                          const Array1 < Array1 <doublevar> > & singletorb,
                          const Array1 < Array1 <doublevar> > & unpairedorb,
                          const Array1 < Array1 <doublevar> > & normalization,
                          doublevar coef_eps
                          )
{
  //Updates only actual value!!!
    int pos=0;
  if (e<npairs(0)){
    //upup start in row
    for(int j=0;j<npairs(0)+npairs(1)+npairs(2);j++){
      pos=order_in_pfaffian(e,j);
      //cout <<"j= "<<j<<endl;
      if (j<npairs(0)){
	mopfaff_row(j)=TripletPairFunction(e, j, moVal, occupation_pos(pos),
                                           tripletorbuu(pos), normalization(pos)(0),coef_eps); 
      }
      else if (j<npairs(0)+npairs(1)){
	mopfaff_row(j)=SingletPairFunction(e, j, moVal, 
                                           occupation_pos(pos), 
                                           singletorb(pos), 
                                           normalization(pos)(2),coef_eps);
      }
      else {
	mopfaff_row(j)=UnpairedOrbFunction(e, j- npairs(0)-npairs(1), 
                                           moVal, occupation_pos(pos), unpairedorb(pos), 
                                           normalization(pos)(3+j- npairs(0)-npairs(1)),coef_eps);
	  //mopfaff_row(j)=moVal(0,e,j-npairs(0));
      }
    }
  }
  else { // if (e<upuppairs+downdownpairs) 
    //-up-down start in row 
    for(int j=0;j<npairs(0)+npairs(1)+npairs(2);j++){
      pos=order_in_pfaffian(e,j);
      if (j<npairs(0)){
	mopfaff_row(j)=-SingletPairFunction(e, j, moVal, 
                                            occupation_pos(pos), 
                                            singletorb(pos), 
                                            normalization(pos)(2),coef_eps);
      }
      else if (j<npairs(0)+npairs(1)){
	mopfaff_row(j)=TripletPairFunction(e, j, moVal, occupation_pos(pos), 
                                           tripletorbdd(pos), 
					   normalization(pos)(1),coef_eps);
      }
      else {
	mopfaff_row(j)=UnpairedOrbFunction(e, j- npairs(0)-npairs(1), 
                                           moVal,occupation_pos(pos),
                                           unpairedorb(pos), 
                                           normalization(pos)(3+j- npairs(0)-npairs(1)),coef_eps);
        //mopfaff_row(j)=moVal(0,e,j-npairs(0));
      } 
     
    }     
  }
}


void UpdatePfaffianRowLap( Array1 < Array1 <doublevar> > & mopfaff_row, 
                           int & e,  
                           const Array3 <doublevar> &  moVal, 
                           const Array1 < Array1 <int> > & occupation_pos,
                           const Array1 <int> & npairs, 
                           const Array2 <int> & order_in_pfaffian,
                           const Array1 < Array1 <doublevar> > & tripletorbuu,
			   const Array1 < Array1 <doublevar> > & tripletorbdd,
                           const Array1 < Array1 <doublevar> > & singletorb,
                           const Array1 < Array1 <doublevar> > & unpairedorb,
                           const Array1 < Array1 <doublevar> > & normalization,
                           doublevar coef_eps
                           )
{
  //Updates only gradients and laplacian, not actual value!!!
  int dmax=moVal.GetDim(0);
  //cout <<"dmax "<<dmax<<endl;
  mopfaff_row.Resize(npairs(0)+npairs(1)+npairs(2));
  Array1 <doublevar> temp(dmax);
  //cout <<"electron number: "<<e<<endl;
  int pos=0;
  if (e<npairs(0)){
    //upup start in row
    for(int j=0;j<npairs(0)+npairs(1)+npairs(2);j++){
      mopfaff_row(j).Resize(dmax);
      pos=order_in_pfaffian(e,j);
      if (j<npairs(0)){
        TripletPairFunctionLap(e, j, moVal,occupation_pos(pos), 
                               tripletorbuu(pos), normalization(pos)(0) , mopfaff_row(j),coef_eps);
      }
      else if (j<npairs(0)+npairs(1)){
        SingletPairFunctionLap(e, j, moVal,occupation_pos(pos), 
                               singletorb(pos), 
                               normalization(pos)(2),
                               mopfaff_row(j),coef_eps);
      }
      else {
        UnpairedOrbFunctionLap(e, j- npairs(0)-npairs(1), moVal, occupation_pos(pos),
                               unpairedorb(pos), normalization(pos)(3+j- npairs(0)-npairs(1)), mopfaff_row(j),coef_eps);
	//for (int d=1;d<5;d++)
	// mopfaff_row(j)(d)=moVal(d,e,j-npairs(0));
      }
    }
  }
  else {
    //-up-down start in row 
    for(int j=0;j<npairs(0)+npairs(1)+npairs(2);j++){
      mopfaff_row(j).Resize(dmax);
      pos=order_in_pfaffian(e,j);
      if (j<npairs(0)){
        SingletPairFunctionLap(e, j, moVal, occupation_pos(pos),
                               singletorb(pos), 
                               normalization(pos)(2),
                               temp,coef_eps);
        for(int d=0;d<dmax;d++)
          mopfaff_row(j)(d)=-temp(d);
      }
      else if (j<npairs(0)+npairs(1)){
        TripletPairFunctionLap(e, j, moVal, occupation_pos(pos),
                               tripletorbdd(pos), normalization(pos)(1), mopfaff_row(j),coef_eps);
      }
      else {
	//mopfaff_row(j)=0.0;
        UnpairedOrbFunctionLap(e, j- npairs(0)-npairs(1), moVal, occupation_pos(pos), 
                               unpairedorb(pos), normalization(pos)(3+j- npairs(0)-npairs(1)), mopfaff_row(j),coef_eps);
	//for (int d=1;d<5;d++)
	//mopfaff_row(j)(d)=moVal(d,e,j-npairs(0));
      } 
    }     
  }
  
}

//----------------------------------------------------------------------

void UpdatePfaffianRowHess( Array1 < Array2 <doublevar> > & mopfaff_row, 
			    int & e,  
			    const Array3 <doublevar> &  moVal, 
			    const Array1 < Array1 <int> > & occupation_pos,
			    const Array1 <int> & npairs, 
			    const Array2 <int> & order_in_pfaffian,
			    const Array1 < Array1 <doublevar> > & tripletorbuu,
			    const Array1 < Array1 <doublevar> > & tripletorbdd,
			    const Array1 < Array1 <doublevar> > & singletorb,
			    const Array1 < Array1 <doublevar> > & unpairedorb,
			    const Array1 < Array1 <doublevar> > & normalization,
                            doublevar  coef_eps
			    )
{
  
  //cout << "UpdatePfaffianRowHess "<<endl;
    //int upuppairs=npairs(0);
  //int downdownpairs=npairs(1);
  //int nopairs=npairs(2);
  int dmax=moVal.GetDim(0);
  //cout <<"dmax "<<dmax<<endl;
  mopfaff_row.Resize(npairs(0)+npairs(1)+npairs(2));

  Array2 <doublevar> temp(dmax,dmax);
  //cout <<"electron number: "<<e<<endl;
  int pos=0;
  if (e<npairs(0)){
    //upup start in row
    for(int j=0;j<npairs(0)+npairs(1)+npairs(2);j++){
      mopfaff_row(j).Resize(dmax,dmax);
      pos=order_in_pfaffian(e,j);
      if (j<npairs(0)){
        TripletPairFunctionHess(e, j, moVal,occupation_pos(pos), 
				tripletorbuu(pos), 
				normalization(pos)(0), 
				mopfaff_row(j),
                                coef_eps);
      }
      else if (j<npairs(0)+npairs(1)){
        SingletPairFunctionHess(e, j, moVal,occupation_pos(pos), 
                                singletorb(pos), 
                                normalization(pos)(2),
                                mopfaff_row(j),
                                coef_eps);
      }
      else {
        UnpairedOrbFunctionHess(e, j- npairs(0)-npairs(1), moVal, occupation_pos(pos),
				unpairedorb(pos), 
				normalization(pos)(3+j- npairs(0)-npairs(1)), 
				mopfaff_row(j),coef_eps);
	//for (int d=1;d<5;d++)
	// mopfaff_row(j)(d)=moVal(d,e,j-npairs(0));
      }
    }
  }
  else {
    //-up-down start in row 
    for(int j=0;j<npairs(0)+npairs(1)+npairs(2);j++){
      mopfaff_row(j).Resize(dmax,dmax);
      pos=order_in_pfaffian(e,j);
      if (j<npairs(0)){
        SingletPairFunctionHess(e, j, moVal, occupation_pos(pos),
                                singletorb(pos), 
                                normalization(pos)(2),
                                temp,coef_eps);
        for(int a=0;a<dmax;a++)
	  for(int b=0;b<dmax;b++)
	    mopfaff_row(j)(a,b)=-temp(a,b);
      }
      else if (j<npairs(0)+npairs(1)){
        TripletPairFunctionHess(e, j, moVal, occupation_pos(pos),
				tripletorbdd(pos), 
				normalization(pos)(1), 
				mopfaff_row(j),coef_eps);
      }
      else {
	//mopfaff_row(j)=0.0;
        UnpairedOrbFunctionHess(e, j- npairs(0)-npairs(1), moVal, occupation_pos(pos), 
				unpairedorb(pos), 
				normalization(pos)(3+j- npairs(0)-npairs(1)), 
				mopfaff_row(j),coef_eps);
	//for (int d=1;d<5;d++)
	//mopfaff_row(j)(d)=moVal(d,e,j-npairs(0));
      } 
    }     
  }
  
}


//----------------------------------------------------------------------
void FillPfaffianMatrix( Array2 <doublevar> & mopfaff_tot,  
                         const Array3 <doublevar> &  moVal,
                         const Array1 < Array1 <int> > & occupation_pos,
                         const Array1 <int> & npairs, 
                         const Array2 <int> & order_in_pfaffian,
                         const Array1 < Array1 <doublevar> > & tripletorbuu,
			 const Array1 < Array1 <doublevar> > & tripletorbdd,
                         const Array1 < Array1 <doublevar> > & singletorb,
                         const Array1 < Array1 <doublevar> > & unpairedorb,
                         const Array1 < Array1 <doublevar> > & normalization,
                         doublevar  coef_eps
                         )
{
  mopfaff_tot=0.0;
  int pos=0;
  for(int i = 0; i < npairs(0)+npairs(1); i++){
    for(int j =i+1; j < npairs(0)+npairs(1)+ npairs(2); j++){
      pos=order_in_pfaffian(i,j);
      //cout <<"pos "<<pos<<endl;
      if (j<npairs(0)){
	mopfaff_tot(i,j)=TripletPairFunction(i, j, moVal, occupation_pos(pos), 
                                             tripletorbuu(pos), 
                                             normalization(pos)(0),coef_eps); 
      }
      else if (j<npairs(0)+npairs(1)){
        if (i<npairs(0)){
	  mopfaff_tot(i,j)=SingletPairFunction(i, j, moVal,
                                               occupation_pos(pos),
                                               singletorb(pos), 
                                               normalization(pos)(2),coef_eps);
	}
        else { // filling row(i) with i> npairs(0)
	  mopfaff_tot(i,j)=TripletPairFunction(i, j, moVal,
                                               occupation_pos(pos), 
                                               tripletorbdd(pos), 
                                               normalization(pos)(1),coef_eps);
	}
      }
      else {
	//if (i<npairs(0))
	mopfaff_tot(i,j)=UnpairedOrbFunction(i, j- npairs(0)-npairs(1), 
                                             moVal, 
                                             occupation_pos(pos), 
                                             unpairedorb(pos), 
                                             normalization(pos)(3+j- npairs(0)-npairs(1)),coef_eps);
        // mopfaff_tot(i,j)=moVal(0,i,j-npairs(0));
         //start with first unpared orbital, ie when second index=N
      }
                              
      mopfaff_tot(j,i)=-mopfaff_tot(i,j); 
    }
  }

}


//----------------------------------------------------------------------
void FillPfaffianMatrixGrad( Array2 <doublevar> & mopfaff_tot,  
			     const Array3 <doublevar> &  moVal,
			     const Array1 < Array1 <int> > & occupation_pos,
			     const Array1 <int> & npairs, 
			     const Array2 <int> & order_in_pfaffian,
			     const int which_pfaffian_func,
			     const int which_part_of_func,
			     const int which_index_of_part_of_func,
			     const Array1 < Array1 <doublevar> > & tripletorbuu,
			     const Array1 < Array1 <doublevar> > & tripletorbdd,
			     const Array1 < Array1 <doublevar> > & singletorb,
			     const Array1 < Array1 <doublevar> > & unpairedorb,
			     const Array1 < Array1 <doublevar> > & normalization,
                             doublevar  coef_eps
			     )
{

  mopfaff_tot=0.0;
  int pos=0;
  for(int i = 0; i < npairs(0)+npairs(1); i++){
    for(int j =i+1; j < npairs(0)+npairs(1)+ npairs(2); j++){
      pos=order_in_pfaffian(i,j);
      if(pos==which_pfaffian_func){
	//cout <<"pos "<<pos<<endl;
	if (j<npairs(0)){
	  if(which_part_of_func==0)
	    mopfaff_tot(i,j)=TripletPairFunctionGrad(i, j, moVal, occupation_pos(pos), 
						     tripletorbuu(pos), 
						     normalization(pos)(0),
						     which_index_of_part_of_func,coef_eps); 
	}
	else if (j<npairs(0)+npairs(1)){
	  if (i<npairs(0)){
	    if(which_part_of_func==2)
	      mopfaff_tot(i,j)=SingletPairFunctionGrad(i, j, moVal,
						       occupation_pos(pos),
						       singletorb(pos), 
						       normalization(pos)(2),
						       which_index_of_part_of_func,coef_eps);
	  }
	  else{ // filling row(i) with i> npairs(0)
	    if(which_part_of_func==1)
	      mopfaff_tot(i,j)=TripletPairFunctionGrad(i, j, moVal,
						       occupation_pos(pos), 
						       tripletorbdd(pos), 
						       normalization(pos)(1),
						       which_index_of_part_of_func,coef_eps);
	  }    
	}
	else {
	  //if (i<npairs(0))
	  if(which_part_of_func==3)
	    mopfaff_tot(i,j)=UnpairedOrbFunctionGrad(i, j- npairs(0)-npairs(1), 
						     moVal, 
						     occupation_pos(pos), 
						     unpairedorb(pos), 
						     normalization(pos)(3+j- npairs(0)-npairs(1)),
						     which_index_of_part_of_func,coef_eps);
	}
      }
      mopfaff_tot(j,i)=-mopfaff_tot(i,j); 
    }
  }
  //cout <<"end: FillPfaffianMatrix"<<endl;
}

//----------------------------------------------------------------------
void FillPfaffianMatrixHess( Array2 <doublevar> & mopfaff_tot,  
			     const Array3 <doublevar> &  moVal,
			     const Array1 < Array1 <int> > & occupation_pos,
			     const Array1 <int> & npairs, 
			     const Array2 <int> & order_in_pfaffian,
			     const int which_pfaffian_func,
			     const int which_part_of_func,
			     const int which_index_of_part_of_func1,
			     const int which_index_of_part_of_func2,
			     const Array1 < Array1 <doublevar> > & tripletorbuu,
			     const Array1 < Array1 <doublevar> > & tripletorbdd,
			     const Array1 < Array1 <doublevar> > & singletorb,
			     const Array1 < Array1 <doublevar> > & unpairedorb,
			     const Array1 < Array1 <doublevar> > & normalization,
                             doublevar coef_eps
			     )
{

  mopfaff_tot=0.0;
  int pos=0;
  for(int i = 0; i < npairs(0)+npairs(1); i++){
    for(int j =i+1; j < npairs(0)+npairs(1)+ npairs(2); j++){
      pos=order_in_pfaffian(i,j);
      //cout <<"pos "<<pos<<endl;
       if(pos==which_pfaffian_func){
	 if (j<npairs(0)){
	   if(which_part_of_func==0)
	     mopfaff_tot(i,j)=TripletPairFunctionHess(i, j, moVal, occupation_pos(pos), 
						      tripletorbuu(pos), 
						      normalization(pos)(0),
						      which_index_of_part_of_func1,
						      which_index_of_part_of_func2,
                                                      coef_eps
						      ); 
	 }
	 else if (j<npairs(0)+npairs(1)){
	   if (i<npairs(0)){
	     if(which_part_of_func==2)
	       mopfaff_tot(i,j)=SingletPairFunctionHess(i, j, moVal,
							occupation_pos(pos),
							singletorb(pos), 
							normalization(pos)(2),
							which_index_of_part_of_func1,
							which_index_of_part_of_func2,
                                                        coef_eps
							);
	   }
	   else{ // filling row(i) with i> npairs(0)
	     if(which_part_of_func==1)
	       mopfaff_tot(i,j)=TripletPairFunctionHess(i, j, moVal,
							occupation_pos(pos), 
							tripletorbdd(pos), 
							normalization(pos)(1),
							which_index_of_part_of_func1,
							which_index_of_part_of_func2,
                                                        coef_eps
							);
	   }    
	 }
	 else {
	   //if (i<npairs(0))
	   if(which_part_of_func==3)
	     mopfaff_tot(i,j)=UnpairedOrbFunctionHess(i, j- npairs(0)-npairs(1), 
						      moVal, 
						      occupation_pos(pos), 
						      unpairedorb(pos), 
						      normalization(pos)(3+j- npairs(0)-npairs(1)),
						      which_index_of_part_of_func1,
						      which_index_of_part_of_func2,
                                                      coef_eps
						      );
	 }
       }
       mopfaff_tot(j,i)=-mopfaff_tot(i,j); 
    }
  }
  //cout <<"end: FillPfaffianMatrix"<<endl;
}



//-----------------------------------------------------------------------
int Pfaff_wf::getParmDeriv(Wavefunction_data *  wfdata, 
			   Sample_point * sample ,
			   Parm_deriv_return & derivatives){

  int nparms_full=parent->nparms();
  int nparms_start=0;
  int nparms_end=nparms_full;
  int nparms=nparms_end-nparms_start;
  if(nparms_end > nparms_full)
    error("nparms_end > nparms_full in Pfaff_wf::getParmDeriv");

  derivatives.gradient.Resize(nparms);
  derivatives.hessian.Resize(nparms, nparms);
 
  if(parent->optimize_pf){

    Array1 <int> which_pfaffian_func(nparms);
    Array1 <int> which_part_of_func(nparms);
    Array1 <int> which_index_of_part_of_func(nparms);
    Array1 < Array1 <doublevar> > pfGrad(nparms);
    Array1 < Array2 <doublevar> > doubletrace(npf);
    Array1 < Array1 <int> > deltacheck(nparms);
    Array1 < Array1 < Array2 <doublevar> > > mopfaff_tot(nparms);
    
    //cout <<"find the pf, k, l indexes for each parameter"<<endl; 
    int counter=0;
    int counter2=0;
    for(int pffunc=0; pffunc< nsfunc; pffunc++){
      for(int k=0;k<parent->optimize_total(pffunc).GetSize();k++){
	for(int l=0;l<parent->optimize_total(pffunc)(k).GetSize();l++)
	  if(parent->optimize_total(pffunc)(k)(l)){
	    if(counter >=nparms_start && counter <nparms_end){
	      which_pfaffian_func(counter2)=pffunc;
	      which_part_of_func(counter2)=k;
	      which_index_of_part_of_func(counter2)=l;
	      counter2++;
	    }
	    counter++;
	  }
      }
    }
    
    assert(nparms==counter2);
    assert(nparms_full==counter);

    //cout << "get pfGrad"<<endl;
    for(int parms=0;parms<nparms;parms++){
      pfGrad(parms).Resize(npf);
      int pffunc=which_pfaffian_func(parms);
      int k=which_part_of_func(parms);
      int l= which_index_of_part_of_func(parms);
      int size=inverse(0).GetDim(0);
      mopfaff_tot(parms).Resize(npf);
      deltacheck(parms).Resize(npf);
      for (int pf=0;pf<npf;pf++){
	deltacheck(parms)(pf)=0;
	pfGrad(parms)(pf)=0;
	mopfaff_tot(parms)(pf).Resize(size,size);
	mopfaff_tot(parms)(pf)=0;
	int yes=0;
	for(int ii=0;ii<parent->order_in_pfaffian(pf).GetDim(0);ii++)
	  for(int jj=ii+1;jj<parent->order_in_pfaffian(pf).GetDim(1);jj++)
	    if(pffunc==parent->order_in_pfaffian(pf)(ii,jj)){
	      yes++;
	      break;
	    }
	if(yes){
	  deltacheck(parms)(pf)=1;
	  FillPfaffianMatrixGrad( mopfaff_tot(parms)(pf),
				  moVal, 
				  parent->occupation_pos,
				  parent->npairs,
				  parent->order_in_pfaffian(pf),
				  pffunc,
				  k,
				  l,
				  parent->tripletorbuu, 
				  parent->tripletorbdd,
				  parent->singletorb,
				  parent->unpairedorb,
				  parent->normalization,
                                  coef_eps
				  );
	  
	  //cout<<" parms "<<parms<<" pf "<<pf<<endl;
	  //for(int i=0;i<size;i++){
	  //for(int j=0;j<size;j++)
	  //  cout <<mopfaff_tot(parms)(pf)(i,j)<<"  ";
	  // cout <<endl;
	  //}
	  
	
	  //multiply by inverse
	  doublevar tmp=0.0;
	  for(int i=0;i<size;i++)
	    for(int j=i+1;j<size;j++){
	      tmp+=2.0*inverse(pf)(i,j)*mopfaff_tot(parms)(pf)(j,i);
	    }
	  pfGrad(parms)(pf)=tmp;
	}
      }
    }

    //cout <<"Get pfHess"<<endl;
    Array2 < Array1 <doublevar> > pfHess(nparms, nparms);
    if(derivatives.need_hessian){
      Array2 <doublevar> mopfaff_tot_dx_dy;
      for(int i=0;i<nparms;i++){
	int pffunc1=which_pfaffian_func(i);
	int k1=which_part_of_func(i);
	int l1= which_index_of_part_of_func(i);
	for(int j=i;j<nparms;j++){
	  int pffunc2=which_pfaffian_func(j);
	  int k2=which_part_of_func(j);
	  int l2= which_index_of_part_of_func(j);
	  if(pffunc1==pffunc2 && k1==k2){
	    int k=k1;
	    int pffunc=pffunc1;
	    int goahead=1;

	    for (int pf=0;pf<npf;pf++){
	      pfHess(i,j).Resize(npf);
	      pfHess(i,j)(pf)=0;
	      int size=inverse(0).GetDim(0);
	      if(goahead){
		mopfaff_tot_dx_dy.Resize(size,size);
		mopfaff_tot_dx_dy=0.0;
		int yes=0;
		for(int ii=0;ii<parent->order_in_pfaffian(pf).GetDim(0);ii++)
		  for(int jj=ii+1;jj<parent->order_in_pfaffian(pf).GetDim(1);jj++)
		    if(pffunc==parent->order_in_pfaffian(pf)(ii,jj)){
		      yes++;
		      break;
		    }
		if(yes){
		  FillPfaffianMatrixHess( mopfaff_tot_dx_dy,
					  moVal, 
					  parent->occupation_pos,
					  parent->npairs,
					  parent->order_in_pfaffian(pf),
					  pffunc,
					  k,
					  l1,
					  l2,
					  parent->tripletorbuu, 
					  parent->tripletorbdd,
					  parent->singletorb,
					  parent->unpairedorb,
					  parent->normalization,
                                          coef_eps
					  );
		  doublevar tmp=0.0;
		  for(int m=0;m<size;m++){
		    for(int n=m+1;n<size;n++){
		      tmp+=2.0*inverse(pf)(m,n)*mopfaff_tot_dx_dy(n,m);
		    }
		  }
		  pfHess(i,j)(pf)=tmp;
		}
	      }
	    }
	  }
	  else{
	    for (int pf=0;pf<npf;pf++){
	      pfHess(i,j).Resize(npf);
	      pfHess(i,j)(pf)=0;
	    }
	  }
	}
      }
      
    
      //cout <<"get 2 body terms"<<endl;
      doublevar temp1,temp2;
      
      for (int pf=0;pf<npf;pf++){
	doubletrace(pf).Resize(nparms,nparms);
	doubletrace(pf)=0;
	int size=inverse(pf).GetDim(0);
	for(int i=0;i<nparms;i++)
	  if(deltacheck(i)(pf))
	    for(int j=i;j<nparms;j++)
	      if(deltacheck(j)(pf))
		for(int ii=0;ii<size;ii++)
		  for(int jj=0;jj<size;jj++)
                    for(int kk=0;kk<size;kk++){
                      temp1=inverse(pf)(ii,jj)*mopfaff_tot(i)(pf)(jj,kk);
		      for(int ll=0;ll<size;ll++){
			temp2=temp1*inverse(pf)(kk,ll)*mopfaff_tot(j)(pf)(ll,ii);
			doubletrace(pf)(i,j)+=temp2;
		      }
		    }
      }
    }
        
    //get full value
    doublevar sum=0;
    for (int pf=0;pf<npf;pf++){
      sum+=parent->pfwt(pf)*pfaffVal(pf);
    }
        
    //cout <<"calculating the actual gradient and hessian"<<endl;
    for(int i=0;i<nparms+1;i++){
      for(int j=i;j<nparms+1;j++){
	if(i==0 && j<nparms ){
	  derivatives.gradient(j)=0.0;
	  for(int pf=0; pf < npf; pf++) {
	     if(deltacheck(j)(pf))
	       derivatives.gradient(j)+=parent->pfwt(pf)*0.5*pfaffVal(pf)*pfGrad(j)(pf);
	  }
	  derivatives.gradient(j)/=sum;
	}
	else if (i>0){
	  derivatives.hessian(i-1,j-1)=0.0;
	  if(derivatives.need_hessian){
	    for(int pf=0; pf < npf; pf++) {
	      derivatives.hessian(i-1,j-1)+=parent->pfwt(pf)*pfaffVal(pf)*
		(0.25*pfGrad(i-1)(pf)*pfGrad(j-1)(pf)
		  +0.5*pfHess(i-1,j-1)(pf)
		 -0.5*doubletrace(pf)(i-1,j-1));
	      //      cout << i-1 <<" "<<j-1<<"  "<<0.25*pfGrad(i-1)(pf)*pfGrad(j-1)(pf)<<"  "
	      //   << 0.5*pfHess(i-1,j-1)(pf) <<"  "<<-0.5*doubletrace(pf)(i-1,j-1)<<endl;
	    }
	  }
	  derivatives.hessian(i-1,j-1)/=sum;
	  derivatives.hessian(j-1,i-1)=derivatives.hessian(i-1,j-1);
	}
      }//end of j 
    }//end of i
    //cout << "end of getParmDeriv"<<endl;
    
    return 1;
  }
  else if (parent->optimize_pfwt ){
    doublevar sum=0;
    for (int pf=0;pf<npf;pf++){
      sum+=parent->pfwt(pf)*pfaffVal(pf);
      if(pf > nparms_start && pf <= nparms_end )
	derivatives.gradient(pf-1-nparms_start)=pfaffVal(pf);
    }
    for(int pf=0; pf < nparms; pf++) {      
      derivatives.gradient(pf)/=sum; 
    }
    derivatives.hessian=0;
    return 1;
  }
  else{
    derivatives.gradient=0;
    derivatives.hessian=0;
    return 1;
  }
  
  return 0;
}

//------------------------------------------------------------------------

void Pfaff_wf::calcVal(Sample_point * sample) {
  calcLap(sample);
}


void Pfaff_wf::updateVal(Sample_point * sample, int e)
{
  int upuppairs=parent->npairs(0);
  int downdownpairs=parent->npairs(1);
  int nopairs=parent->npairs(2);
  int totelectrons=nelectrons(0)+nelectrons(1);
  Array3 <doublevar> moVal_temmp(5, totelectrons, updatedMoVal.GetDim(0));
  doublevar ratio;
  

  //cout << "updateVal" << endl;
  assert(parent != NULL);
  for (int pf=0;pf<npf;pf++){
    if(pfaffVal(pf)==0){
      cout << "updateVal::WARNING: pfaffian  zero!" << endl;
      calcLap(sample);
      return;
    }
  }
  
  if (updateMoLap==1){
    sample->updateEIDist();
    updatedMoVal=0;
    //update all the mo's that we will be using.
    parent->molecorb->updateVal(sample, e,
				 0,
				 updatedMoVal);
    
    //  for(int i=0; i< updatedMoVal.GetDim(0); i++) {
    //  cout << "updatedMoVal " << updatedMoVal(i,0) << endl;
    //}
    
    int size=updatedMoVal.GetDim(0);
    for(int i=0; i< size; i++) {
      moVal(0,e,i)=updatedMoVal(i,0);
    }
  }

   
  Array1 <doublevar> mopfaff_row(upuppairs+downdownpairs+nopairs);
  Array1 <doublevar> mopfaff_column(upuppairs+downdownpairs+nopairs);

  

 
  for (int pf=0;pf<npf;pf++){
    //cout <<" Update e-th row& column for pfaffian"<<endl;
    UpdatePfaffianRowVal(mopfaff_row, 
                         e,  
                         moVal,
                         parent->occupation_pos,
                         parent->npairs, 
                         parent->order_in_pfaffian(pf),
                         parent->tripletorbuu, 
                         parent->tripletorbdd,
                         parent->singletorb,
                         parent->unpairedorb,
                         parent->normalization,
                         coef_eps
                         );
   
    //cout <<" After Update e-th row& column for pfaffian"<<endl;

    ratio= UpdateInversePfaffianMatrix(inverse(pf), mopfaff_row, mopfaff_column, e);
  
    // cout << "Pfaffval before " << pfaffVal << "ratio " << ratio << endl;
    //update detVal
    pfaffVal(pf)=ratio*pfaffVal(pf);
    //cout << "Pfaffian value in UpdateVal is "<<pfaffVal(pf)<<endl;
  }

  
  //cout << "End updateVal " << endl;

}

//------------------------------------------------------------------------

void Pfaff_wf::getVal(Wf_return & val)
{
  //cout << "getVal"<<endl;
  Array1 <doublevar> si(1, 0.0);
  Array2 <doublevar> vals(1,1,0.0);

  assert(val.amp.GetDim(0) >=1);
  assert(val.amp.GetDim(1) >= 1);
  assert(val.amp.GetDim(1) >= 1);
  
  doublevar tempval=0;
  
  for (int pf=0;pf<npf;pf++){
    tempval+=parent->pfwt(pf)*pfaffVal(pf);
  }
  si(0)=sign(tempval);
  
  if(fabs(tempval) > 0)
    vals(0,0)=log(fabs(tempval));
  else
    vals(0,0)=-1e3;
  val.setVals(vals, si);
  
  //cout << "End getVal"<<endl;
}

//-----------------------------------------------------------------------
void Pfaff_wf::getSymmetricVal(Wavefunction_data * wfdata,
		     int e, Wf_return & val){

  Array1 <doublevar> si(1, 0.0);
  Array2 <doublevar> vals(1,1,0.0);
  val.setVals(vals, si);
} 


//----------------------------------------------------------------------


void Pfaff_wf::calcLap(Sample_point * sample)
{
  Array1 <doublevar> elecpos(3);
  //cout << "calcLap" << endl;
  // int upuppairs=parent->npairs(0);
  // int downdownpairs=parent->npairs(1);
  // int nopairs=parent->npairs(2);
  //int ntote_pairs=parent->ntote_pairs;
  int tote=parent->tote;

  //cout << "calcLap " << endl;
  for(int e=0; e< tote; e++)
  {
    //cout << " e " << e<< endl;
   
    if (updateMoLap==1){
      sample->updateEIDist();
      //update all the mo's that we will be using, using the lists made in
      //Pfaff_wf_data(one for each spin).

      parent->molecorb->updateLap(sample, e,
                                   0,
                                   updatedMoVal);
      

      //   for(int i=0; i< updatedMoVal.GetDim(0); i++) {
      // cout << "updatedMoVal " << updatedMoVal(i,0) << endl;
      //}
      
      assert(parent->totoccupation(0).GetSize()==updatedMoVal.GetDim(0));
    
      int size=updatedMoVal.GetDim(0);
      for(int d=0; d< 5; d++){
        for(int i=0; i< size; i++){
          moVal(d,e,i)=updatedMoVal(i,d);
        }
      }
      // sample->getElectronPos(e, elecpos);
    }
  }

  //cout << "moVal(0,0,2) "<<moVal(0,0,2);
  //cout <<  moVal(0,0,0)*moVal(0,0,0)+ moVal(0,0,1)*moVal(0,0,1)<<endl;

  //cout <<"calculate the whole pfaffian matrix"<<endl;
  

  for (int pf=0;pf<npf;pf++){
  //for(int e=0; e< tote; e++){
  //  for(int i=0;i<updatedMoVal.GetDim(0);i++){
       //cout <<parent->occupation_pos(pf)(i)<<" , "<<endl;
  //    cout <<moVal(0,e,i)<<" , "<<endl;
       
  //}
  //   cout <<endl;
  //}
    //cout <<"parent->occupation_pos(pf).GetSize() "<<parent->occupation_pos(pf).GetSize()<<endl;
    
    FillPfaffianMatrix( mopfaff_tot(pf),
                        moVal, 
                        parent->occupation_pos,
                        parent->npairs,
                        parent->order_in_pfaffian(pf),
                        parent->tripletorbuu, 
                        parent->tripletorbdd,
                        parent->singletorb,
                        parent->unpairedorb,
                        parent->normalization,
                        coef_eps
                        );
    // cout << endl;
    
    /*
    cout << "Pfaffian number "<<pf+1<<" total_pf_number "<<npf<<endl;
    
    cout.setf(ios::scientific| ios:: showpos);
    for(int i = 0; i < mopfaff_tot(pf).GetDim(0); i++){
      for(int j =0 ; j <  mopfaff_tot(pf).GetDim(1); j++){
        cout << mopfaff_tot(pf)(i,j) <<"  ";  
      }
      cout << endl;
    }
    cout.unsetf(ios::scientific| ios:: showpos);
    cout <<endl;
    */
    // cout << "Sqrt of Pfaffian matrix is: "
    //	 <<sqrt(Determinant(mopfaff_tot,mopfaff_tot.GetDim(0)))<<endl;
    pfaffVal(pf) = PfaffianInverseMatrix(mopfaff_tot(pf), inverse(pf));
    //cout << "Pfaffian value:             "<< pfaffVal(pf) << endl;

    
    
  }
  //cout << "----------"<<endl;
  
  //cout << parent->check_pfwt_sign<< endl;
  if(parent->check_pfwt_sign && firstime ){
    for (int pf=0;pf<npf;pf++)
      cout << parent->pfwt(pf)*pfaffVal(pf)<<endl;
    firstime=false;
  }
       
}
 
//------------------------------------------------------------------------


void Pfaff_wf::getLap(int e, Wf_return & lap) {

  //cout << "getLap"<<endl;
  Array1 <doublevar> si(1, 0.0);
  Array2 <doublevar> vals(1,5,0.0);
  
  
  int upuppairs=parent->npairs(0);
  int downdownpairs=parent->npairs(1);
  int nopairs=parent->npairs(2);
  //int ntote_pairs=parent->ntote_pairs;
  
  int shiftf=0;
  
  Array1 < Array1 < Array1 <doublevar> >  > mopfaff_row;
  mopfaff_row.Resize(npf);
  
  
  doublevar funcval=0;
  for (int pf=0;pf<npf;pf++){
    funcval+=parent->pfwt(pf)*pfaffVal(pf);
  }
  
  si(shiftf)=sign(funcval);
  
  if(fabs(funcval) > 0)
    vals(shiftf,0)=log(fabs(funcval));
  else vals(shiftf,0)=-1e3;
  
  //   cout << "lap(0,0)" << lap(0,0) << endl;
  //cout << "funcval " << funcval << endl;
  
  for (int pf=0;pf<npf;pf++){
    //mopfaff_row(pf).Resize(upuppairs+downdownpairs+nopairs);
    //for (int i=0;i<upuppairs+downdownpairs+nopairs;i++)
    // mopfaff_row(pf)(i).Resize(5); //now done inside UpdatePfaffianRowLap
    UpdatePfaffianRowLap(mopfaff_row(pf),
                         e,
                         moVal,
                         parent->occupation_pos,
                         parent->npairs,
                         parent->order_in_pfaffian(pf),
                         parent->tripletorbuu,
                         parent->tripletorbdd,
                         parent->singletorb,
                         parent->unpairedorb,
                         parent->normalization,
                         coef_eps
                         );
  }
  // cout.setf(ios::scientific| ios:: showpos);
  for(int d=1; d< 5; d++){
    vals(shiftf,d)=0.0;
    for (int pf=0;pf<npf;pf++){
      doublevar temp=0.0;
      for(int j=0; j<upuppairs+downdownpairs+nopairs; j++){
        temp+=mopfaff_row(pf)(j)(d)*inverse(pf)(j,e);
        //updated row*inverse matrix;
      }
      if(parent->pfwt(pf)==0)
        temp=0.0;
      if(npf >1)
        temp*=parent->pfwt(pf)*pfaffVal(pf);
      vals(shiftf,d)+=temp;
    }
    if(funcval==0)
      vals(shiftf,d)=0;
    else if (npf >1)
      vals(shiftf,d)/=funcval;
  }
  
  lap.setVals(vals, si);
  /*
  cout.unsetf(ios::scientific| ios:: showpos);
  cout << "lap  ";
  for(int d=0; d< 5; d++) {
   cout << lap.amp(0, d) << "   ";
  }
  cout <<endl;
  */
  //cout << "End of getLap"<<endl;
}

//-------------------------------------------------------------------------

/*!
*/
void Pfaff_wf::updateLap(
  Sample_point * sample,
  int e)
{
  
  //cout << "updateLap: Single row-column update"<<endl;
  assert(parent != NULL);

  int upuppairs=parent->npairs(0);
  int downdownpairs=parent->npairs(1);
  int nopairs=parent->npairs(2);
  //  int ntote_pairs=parent->ntote_pairs;
  Array1 <doublevar> elecpos(3);

  
  doublevar ratio;
  for (int pf=0;pf<npf;pf++){
    if(pfaffVal(pf)==0){
      cout << "updateLap::WARNING: Pfaffian zero!" << endl;
      calcLap(sample);
      return;
    }
  }
 
  if (updateMoLap==1){
    sample->updateEIDist();
    //cout << "mo update\n";
    //update all the mo's that we will be using.
    parent->molecorb->updateLap(sample, e,
                                 0,
                                 updatedMoVal);
  

  // sample->getElectronPos(e, elecpos);
  //cout <<  "elecpos:  "<< elecpos(0)<<"  "<< elecpos(1)<<"  "<< elecpos(2)<<endl;
    
  //  for(int i=0; i< updatedMoVal.GetDim(0); i++) {
  // cout << "updatedMoVal " << updatedMoVal(i,0) << endl;
  //}

    int size=updatedMoVal.GetDim(0);
    for(int d=0; d< 5; d++){
      for(int i=0; i< size; i++){
        moVal(d,e,i)=updatedMoVal(i,d);
      }
    }
  
  }

  //cout << "moVal(0,0,2) "<<moVal(0,0,2);

  Array1 <doublevar> mopfaff_row(upuppairs+downdownpairs+nopairs);
  Array1 <doublevar> mopfaff_column(upuppairs+downdownpairs+nopairs);
 
  for (int pf=0;pf<npf;pf++){
    UpdatePfaffianRowVal(mopfaff_row, 
                         e,  
                         moVal, 
                         parent->occupation_pos,
                         parent->npairs,
                         parent->order_in_pfaffian(pf),
                         parent->tripletorbuu, 
                         parent->tripletorbdd,
                         parent->singletorb,
                         parent->unpairedorb,
                         parent->normalization,
                         coef_eps
                         );
  
    //ratio pf(new)/pf(old);
    ratio= UpdateInversePfaffianMatrix(inverse(pf), mopfaff_row, mopfaff_column, e);
  
    //cout << "Pfaffval before# " << pfaffVal << "  ratio#  " << ratio << endl;

    pfaffVal(pf)*=ratio;
  }
  //cout << "New value#  "<<pfaffVal<<endl;
}
//-------------------------------------------------------------------------

