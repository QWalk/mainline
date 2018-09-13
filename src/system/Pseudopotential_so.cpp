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


#include "Pseudopotential_so.h"
#include "qmc_io.h"
#include "Wavefunction.h"
#include "Sample_point.h"
#include <iomanip>
#include "ulec.h"
#include "Wavefunction_data.h"
#include "Basis_function.h"
#include "System.h"
using namespace std;
/*!
Some legendre polynomials and spherical harmonics..
*/

doublevar legendre_so(doublevar x, int n)
{
  switch(n)
  {
  case 0:
    return 1;
  case 1:
    return x;
  case 2:
    return .5*(3*x*x-1);
  case 3:
    return .5*(5*x*x*x - 3*x);
  case 4:
    return 0.125*(35*x*x*x*x -30*x*x +3);
  case 5:
    return 0.125*(63*x*x*x*x*x - 70*x*x*x + 15*x);
  default:
    error("Do not have legendre polynomial of order ", n);
    return 0; //shouldn't get here, but gets rid of a compiler message
  }
}

doublevar legendre_der_so(doublevar x, int n) {
  switch(n) {
  case 0: 
    return 0;
  case 1:
    return 1;
  case 2:
    return 3*x;
  case 3:
    return .5*(15*x*x-3);
  case 4:
    return 35*x*x*x-15*x;
  case 5:
    return .125*(315*x*x*x*x-210*x*x+15);
  default:
    error("Do not have legendre polynomial of order ", n);
    return 0;
  }
}

dcomplex Ylm(int l, int ml, Array1 <doublevar> r){
  assert(r.GetDim(0)==5);
  switch(l)
  {
  case 0:
    return 1;

  case 1:
    switch(ml){
      case -1:
        return sqrt(3.00/2.00)*dcomplex(r(2)/r(0),-r(3)/r(0));
      case 0:
        return sqrt(3.00)*dcomplex(r(4)/r(0), 0);
      case 1:
        return -sqrt(3.00/2.00)*dcomplex(r(2)/r(0),r(3)/r(0));
      default:
 //       error("Do not have spherical harmonics of order", ml);
        return 0;
    }

  case 2:
    switch(ml){
      case -2:
        return sqrt(15.00/8.00)*dcomplex( (r(2)*r(2)-r(3)*r(3))/r(1), -2*r(2)*r(3)/r(1) );
      case -1:
        return sqrt(15.00/2.00)*dcomplex(r(2)*r(4)/r(1),-r(3)*r(4)/r(1));
      case 0:
        return sqrt(5.00/4.00)*dcomplex((2*r(4)*r(4)-r(2)*r(2)-r(3)*r(3))/r(1),0);
      case 1:
        return -sqrt(15.00/2.00)*dcomplex(r(2)*r(4)/r(1),r(3)*r(4)/r(1));
      case 2:
        return sqrt(15.00/8.00)*dcomplex( (r(2)*r(2)-r(3)*r(3))/r(1), 2*r(2)*r(3)/r(1) );
      default:
   //     error("Do not have spherical harmonics of order", ml);
        return 0;
    }

  case 3:
    switch(ml){
      case -3:
        return sqrt(35.00/16.00)*dcomplex( r(2)*(r(2)*r(2)-3*r(3)*r(3))/(r(0)*r(1)),
                                           -r(3)*(3*r(2)*r(2)-r(3)*r(3))/(r(0)*r(1)) );
      case -2:
        return sqrt(105.00/8.00)*dcomplex( r(4)*(r(2)*r(2)-r(3)*r(3))/(r(0)*r(1)),
                                           -2*r(2)*r(3)*r(4)/(r(0)*r(1))  );
      case -1:
        return sqrt(21.00/16.00)*dcomplex(r(2)*(4*r(4)*r(4)-r(2)*r(2)-r(3)*r(3))/(r(0)*r(1)),
                                          -r(3)*(4*r(4)*r(4)-r(2)*r(2)-r(3)*r(3))/(r(0)*r(1)));
      case 0:
        return sqrt(7.00/4.00)*dcomplex( r(4)*(2*r(4)*r(4)-3*r(2)*r(2)-3*r(3)*r(3))/(r(0)*r(1)), 0 );
      case 1:
        return -sqrt(21.00/16.00)*dcomplex(r(2)*(4*r(4)*r(4)-r(2)*r(2)-r(3)*r(3))/(r(0)*r(1)),
                                          r(3)*(4*r(4)*r(4)-r(2)*r(2)-r(3)*r(3))/(r(0)*r(1)));
      case 2:
        return sqrt(105.00/8.00)*dcomplex( r(4)*(r(2)*r(2)-r(3)*r(3))/(r(0)*r(1)),
                                           2*r(2)*r(3)*r(4)/(r(0)*r(1))  );
      case 3:
        return -sqrt(35.00/16.00)*dcomplex( r(2)*(r(2)*r(2)-3*r(3)*r(3))/(r(0)*r(1)),
                                           r(3)*(3*r(2)*r(2)-r(3)*r(3))/(r(0)*r(1)) );
      default:
     //   error("Do not have spherical harmonics of order", ml);
        return 0;
    }

  case 4:
    switch(ml){
      case -4:
        return 3.00*sqrt(70.00)/16.00*pow(dcomplex(r(2),-r(3)),4)/pow(r(1),2);
      case -3:
        return 3.00*sqrt(35.00)/4.00*pow(dcomplex(r(2),-r(3)),3)*r(4)/pow(r(1),2);
      case -2:
        return 3.00*sqrt(10.00)/8.00*pow(dcomplex(r(2),-r(3)),2)*(7*r(4)*r(4)-r(1))/pow(r(1),2);
      case -1:
        return 3.00*sqrt(5.00)/4.00*dcomplex(r(2),-r(3))*r(4)*(7*r(4)*r(4)-3*r(1))/pow(r(1),2);
      case 0:
        return 3.00/8.00*dcomplex( (35*r(4)*r(4)*r(4)*r(4)-30*r(1)*r(4)*r(4)+3*r(1)*r(1))/(r(1)*r(1)), 0 );
      case 1:
        return -3.00*sqrt(5.00)/4.00*dcomplex(r(2),r(3))*r(4)*(7*r(4)*r(4)-3*r(1))/pow(r(1),2);
      case 2:
        return 3.00*sqrt(10.00)/8.00*pow(dcomplex(r(2),r(3)),2)*(7*r(4)*r(4)-r(1))/pow(r(1),2);
      case 3:
        return -3.00*sqrt(35.00)/4.00*pow(dcomplex(r(2),r(3)),3)*r(4)/pow(r(1),2);
      case 4:
        return 3.00*sqrt(70.00)/16.00*pow(dcomplex(r(2),r(3)),4)/pow(r(1),2);
      default:
 //       error("Do not have spherical harmonics of order", ml);
        return 0;
    }

  default:
   // error("Do not have spherical harmonics of order", l);
    return 0;
  }
}

dcomplex kron_del(int a, int b){
  if (a==b){ return dcomplex(1.,0.);}
  else { return dcomplex(0.,0.);}
}

dcomplex Ylm_c(int l, int m, Array1 <doublevar> r){
  return conj(Ylm(l,m,r));
}

dcomplex integral_ls_x(int l, int m, int m_p){
  return 0.5*sqrt((l-m_p)*(l+m_p+1))*kron_del(m,m_p+1)+0.5*sqrt((l+m_p)*(l-m_p+1))*kron_del(m,m_p-1);
}

dcomplex integral_ls_y(int l, int m, int m_p){
  return dcomplex(0.,-0.5)*sqrt((l-m_p)*(l+m_p+1))*kron_del(m,m_p+1)
         +dcomplex(0.,0.5)*sqrt((l+m_p)*(l-m_p+1))*kron_del(m,m_p-1);
}

dcomplex integral_ls_z(int l, int m, int m_p){
  return dcomplex(m_p,0.0)*kron_del(m,m_p);
}



dcomplex Yl1(Array1 <doublevar> r, int l){
  assert(r.GetDim(0)==5); 
  switch(l)
  {
  case 0:
    return 0;
  case 1:
    return -sqrt(3.00)*dcomplex(r(2)/r(0),r(3)/r(0));
  case 2:
    return -sqrt(15.00)*dcomplex(r(2)*r(4)/(r(1)),r(3)*r(4)/(r(1)));
  case 3:
    return -sqrt(21.00/16.00)*dcomplex(r(2)*(4*r(4)*r(4)-r(2)*r(2)-r(3)*r(3))/(r(0)*r(1)), 
                                       r(3)*(4*r(4)*r(4)-r(2)*r(2)-r(3)*r(3))/(r(0)*r(1)));
  case 4:
    return -sqrt(45.00/16.00)*dcomplex(r(2)*r(4)*(7*r(4)*r(4)-3*r(1))/(r(1)*r(1)), 
                                       r(3)*r(4)*(7*r(4)*r(4)-3*r(1))/(r(1)*r(1)));
  default:
    error("Do not have spherical harmonics of order", l);
    return 0;
  }
}

dcomplex Yl_1(Array1 <doublevar> r, int l){
  assert(r.GetDim(0)==5); 
  switch(l)
  {
  case 0:
    return 0;
  case 1:
    return sqrt(3.00)*dcomplex(r(2)/r(0),-r(3)/r(0));
  case 2:
    return sqrt(15.00)*dcomplex(r(2)*r(4)/(r(1)),-r(3)*r(4)/(r(1)));
  case 3:
    return sqrt(21.00/16.00)*dcomplex(r(2)*(4*r(4)*r(4)-r(2)*r(2)-r(3)*r(3))/(r(0)*r(1)), 
                                      -r(3)*(4*r(4)*r(4)-r(2)*r(2)-r(3)*r(3))/(r(0)*r(1)));
  case 4:
    return sqrt(45.00/16.00)*dcomplex(r(2)*r(4)*(7*r(4)*r(4)-3*r(1))/(r(1)*r(1)), 
                                      -r(3)*r(4)*(7*r(4)*r(4)-3*r(1))/(r(1)*r(1)));
  default:
    error("Do not have spherical harmonics of order", l);
    return 0;
  }
}
//----------------------------------------------------------------------


Pseudopotential_so::~Pseudopotential_so()  {
  for(int i=0; i< radial_basis.GetDim(0); i++)
    for(int j=0; j < radial_basis.GetDim(1); j++)  
      for(int k=0; k < radial_basis.GetDim(2); k++)  
        if(radial_basis(i,j,k)) delete radial_basis(i,j,k);
}


//----------------------------------------------------------------------


/*!
Evaluates the radial part of the pseudopotential for a given
atom and l-value.

 */

void Pseudopotential_so::getRadial(int at, int spin,Sample_point * sample,
                                Array1 <doublevar> & r, Array2 <doublevar> & v_jl) {
  assert(radial_basis(at,spin,0) != NULL);
  //cout << "v_jl size=" << v_jl.GetDim(1) << endl;
  Array1 <doublevar> tempv_jl(v_jl.GetDim(1));
  for(int j=0;j<=1;j++){
    tempv_jl=v_jl(j);
    radial_basis(at,spin,j)->calcVal(r, tempv_jl);
    for(int l=0;l<tempv_jl.GetDim(0);l++){ 
      v_jl(j,l)=tempv_jl(l);
   //   cout << j <<" "<< l <<" " << tempv_jl(l) << " " << v_jl(j,l) << endl;
    }
    if(addzeff(at)) {   
    //for(int l=0; l < numL(at); l++) {
      int l=numL(at)-1;           // the actual l=0
      doublevar cutoff_rad=radial_basis(at,spin,j)->cutoff(l);
      if(r(0) < cutoff_rad) {
        v_jl(j,l) += sample->getIonCharge(at)/r(0);
   //     cout << v_jl(j,l) << endl;
      }
    }
  }
}


void Pseudopotential_so::getRadial(int at, int spin, Sample_point * sample,
                                Array1 <doublevar> & r, 
                                Array3 <doublevar> & v_jl){
  assert(radial_basis(at,spin,0) != NULL);
  Array2 <doublevar> tempv_jl(v_jl.GetDim(1),v_jl.GetDim(2));
  for(int j=0;j<=1;j++){
    tempv_jl=v_jl(j);
    radial_basis(at,spin,j)->calcLap(r, tempv_jl);
    for(int l0=0;l0<tempv_jl.GetDim(0);l0++)
      for(int l1=0;l1<tempv_jl.GetDim(1);l1++)
        v_jl(j,l0,l1)=tempv_jl(l0,l1);
    
    if(addzeff(at)) {
      int l=numL(at)-1;
      doublevar cutoff_rad=radial_basis(at,spin,j)->cutoff(l);
      if(r(0) < cutoff_rad) {
        v_jl(j,l,0) += sample->getIonCharge(at)/r(0);
        for(int d=0; d< 3; d++) {
          v_jl(j,l,d+1) -= sample->getIonCharge(at)*r(d+2)/(r(0)*r(1));
        }
      }  
    }
  }
}



//----------------------------------------------------------------------


int Pseudopotential_so::nTest() {
  int tot=0;
  int natoms=numL.GetDim(0);
  for(int at=0; at < natoms; at++) {
    if(numL(at) != 0) {
      tot+=nelectrons;
    }
  }
  return tot;
}


//----------------------------------------------------------------------

void Pseudopotential_so::calcNonloc(Wavefunction_data * wfdata, System * sys,
                                 Sample_point * sample, Wavefunction * wf,
                                 Array2 <doublevar> & totalv) {
  int tot=nTest();
  Array1 <doublevar> test(tot);
  for(int i=0; i< tot; i++) {
    test(i)=rng.ulec();
  }
  calcNonlocWithTest(wfdata, sys, sample, wf, test, totalv);
}

//----------------------------------------------------------------------


void Pseudopotential_so::calcNonlocTmove(Wavefunction_data * wfdata, System * sys,
                     Sample_point * sample,
                     Wavefunction * wf,
                     Array2 <doublevar> & totalv,  //total p.e. from the psp
                     vector <Tmove> & tmoves  //variables for T-moves of Casula
                     ) { 
  int tot=nTest();
  Array1 <doublevar> test(tot);
  for(int i=0; i< tot; i++) {
    test(i)=rng.ulec();
  }
  Array1 <doublevar> parm_deriv;
  Array3 <doublevar> totalv_alle(nelectrons, wf->nfunc(),3);  // the last dim: x,y,z
  totalv = 0.0; 
  calcNonlocWithAllvariables(wfdata,sys,sample, wf, test,totalv_alle, Tmoves::negative_tmove, tmoves,false, parm_deriv);
  for (int e=0; e<nelectrons; e++) 
    for (int w = 0; w < wf->nfunc(); w++)
      for (int idim = 0; idim <=2; idim++)
        totalv(w,idim) += totalv_alle(e, w, idim);
}

//----------------------------------------------------------------------

void Pseudopotential_so::calcNonlocSeparated(Wavefunction_data * wfdata, System * sys,
					  Sample_point * sample,
					  Wavefunction * wf,
					  Array3 <doublevar> & totalv
					  ) 
{ 
  int tot=nTest();
  Array1 <doublevar> test(tot);
  for(int i=0; i< tot; i++) {
    test(i)=rng.ulec();
  }
  Array1 <doublevar> parm_deriv;
  //  Array1 <doublevar> totalv(wf->nfunc(),0.0);
  vector<Tmove> tmoves; 
  calcNonlocWithAllvariables(wfdata,sys,sample, wf, test, totalv, Tmoves::no_tmove, tmoves,false, parm_deriv);
}

//----------------------------------------------------------------------

void Pseudopotential_so::calcNonlocWithTest(Wavefunction_data *wfdata , System * sys, 
                                         Sample_point * sample, Wavefunction *wf ,
                                         const Array1 <doublevar> & accept_var,
                                         Array2 <doublevar> & totalv) { 
  vector<Tmove>  tmoves;
  Array1 <doublevar> parm_deriv;
  Array3 <doublevar> totalv_alle(nelectrons, wf->nfunc(),3);
  totalv = 0.0;
  //cout << "before calling withallvariables, see pspso:393" << endl; 
  calcNonlocWithAllvariables(wfdata,sys, sample, wf, accept_var,totalv_alle, Tmoves::no_tmove, tmoves,false, parm_deriv);
  for (int w = 0; w < wf->nfunc(); w++) 
    for (int e=0; e<nelectrons; e++) 
      for (int idim=0; idim<=2; idim++)
        totalv(w,idim) += totalv_alle(e, w,idim);
}

void Pseudopotential_so::calcNonlocParmDeriv(Wavefunction_data * wfdata, System * sys,
                                          Sample_point * sample,
                                          Wavefunction * wf,
                                          const Array1 <doublevar> & accept_var,
                                          Array2 <doublevar> & totalv, Array1 <doublevar> & parm_deriv) { 
  vector<Tmove>  tmoves;
  assert(totalv.GetDim(0)>=wf->nfunc());
  Array3 <doublevar> totalv_alle(nelectrons,wf->nfunc(),3);
  calcNonlocWithAllvariables(wfdata,sys,sample, wf, accept_var,totalv_alle, Tmoves::no_tmove, tmoves,true, parm_deriv);
  totalv=0.0;
  for(int e=0; e< nelectrons; e++) 
    for(int w=0; w< wf->nfunc(); w++) 
      for (int idim=0; idim<=2; idim++)
        totalv(w,idim)+=totalv_alle(e,w,idim);
  
}



/*
not used in Psp_so
*/
void Pseudopotential_so::calcPseudoSeparated(Wavefunction_data * wfdata,
					  System * sys,
					  Sample_point * sample,
					  Wavefunction * wf,
					  const Array1 <doublevar> & accept_var,
					  Array2 <doublevar> & totalv)//, 
{
  int natoms=sample->ionSize();
  int nwf=wf->nfunc();
  assert(accept_var.GetDim(0) >= nTest());
  //assert(totalv.GetDim(0) >= nwf);
  assert(nelectrons == sample->electronSize());

  Array1 <doublevar> ionpos(3), oldpos(3), newpos(3);
  Array1 <doublevar> newdist(5), olddist(5);
  Wf_return val(nwf,2);
  
  //totalv=0;
  //  doublevar accum_local=0;
  //  doublevar accum_nonlocal=0;
  
  wf->updateVal(wfdata, sample);
  wfStore.initialize(sample, wf);
  
  int accept_counter=0;
  //deriv.Resize(natoms, 3);
  //deriv=0;
  Array1 <doublevar> nonlocal(nwf);
  Array2 <doublevar> integralpts_real(nwf, maxaip);
  Array2 <doublevar> integralpts_imag(nwf, maxaip);
  Array1 <doublevar> rDotR(maxaip);
  for(int at = 0; at< natoms; at++){
    if(numL(at) != 0) {
      Array2 <doublevar> v_jl(2,numL(at));
      sample->getIonPos(at, ionpos);
      for (int e=0; e<sample->electronSize(); e++) {
	sample->getElectronPos(e, oldpos);
	
	//note: this updateEIDist might become inefficient..
	//depends on how often we're rejecting/how much it costs
	//to do it.  If needed, we should add an interface to 
	//Sample_point
	sample->updateEIDist();
	sample->getEIDist(e, at, olddist);
	nonlocal=0.0; 
	
	int spin=1;
	if(e < sys->nelectrons(0)) spin=0;
          getRadial(at,spin, sample, olddist, v_jl);
	//----------------------------------------
	//Start integral
	
	int accept;
	if(deterministic) {
	  accept= olddist(0) < cutoff(at);
	}
	else {
	  doublevar strength=0;
	  const doublevar calculate_threshold=10;
	  
	  for(int l=0; l<numL(at)-1; l++)
            for(int j=0;j<=1;j++)
	      strength+=calculate_threshold*(2*l+1)*fabs(v_jl(j,l));
	  
	  strength=min((doublevar) 1.0, strength);
	  
	  doublevar rand=accept_var(accept_counter++);

	  if (strength > 0.0) {
	    for(int l=0; l<numL(at)-1; l++)
              for(int j=0;j<=1;j++)
	        v_jl(j,l)/=strength;
	  }
	  accept=strength>rand;
	}
	
	//bool localonly = true;
	if(accept)  {
	  wfStore.saveUpdate(sample, wf, e);
	  Wf_return  oldWfVal(nwf,2);
	  wf->getVal(wfdata, e,oldWfVal);
	  
	  for(int i=0; i< aip(at); i++) {
	    sample->setElectronPos(e, oldpos);
	    doublevar base_sign=sample->overallSign();
	    doublevar base_phase=sample->overallPhase();
         
	    for(int d=0; d < 3; d++) 
	      newpos(d)=integralpt(at,i,d)*olddist(0)-olddist(d+2);
	    sample->translateElectron(e, newpos);
	    sample->updateEIDist();
	    sample->getEIDist(e,at,newdist);
	    rDotR(i)=0;
	    for(int d=0; d < 3; d++)
	      rDotR(i)+=newdist(d+2)*olddist(d+2);
	    doublevar new_sign=sample->overallSign();
	    doublevar new_phase=sample->overallPhase();
	    
	    rDotR(i)/=(newdist(0)*olddist(0));  //divide by the magnitudes
            

//            doublevar theta_i=acos(olddist(4)/olddist(0)); 
 //           doublevar phi_i=atan(olddist(3)/olddist(2));
//            Array1 <doublevar> newdistrot;   //rotate newdist to frame where olddist is along the z axis
//            newdistrot(0)=newdist(0); newdistrot(1)=newdist(1);
//            newdistrot(2)=-sin(phi_i)*newdist(2)+cos(phi_i)*newdist(3);
 //           newdistrot(3)=-cos(theta_i)*cos(phi_i)*newdist(2) 
 //                         -cos(theta_i)*sin(phi_i)*newdist(3)
 //                         +sin(theta_i)*newdist(4); 
 //           newdistrot(4)=sin(theta_i)*cos(phi_i)*newdist(2) 
 //                         +sin(theta_i)*sin(phi_i)*newdist(3)
 //                         +cos(theta_i)*newdist(4); 
	    
            wf->updateVal(wfdata, sample);
	    wf->getVal(wfdata, e, val); 
	    //----
	    for(int w=0; w< nwf; w++) {
	      integralpts_real(w,i)=exp(val.amp(w, 0) - oldWfVal.amp(w, 0))*integralweight(at, i);
	      integralpts_imag(w,i)=integralpts_real(w,i);
	      if (val.is_complex==1) {
		integralpts_real(w, i)*=cos(val.phase(w, 0)+new_phase - oldWfVal.phase(w, 0) - base_phase); 
		integralpts_imag(w, i)*=sin(val.phase(w, 0)+new_phase - oldWfVal.phase(w, 0) - base_phase); 
	      } else {
		integralpts_real(w, i)*=val.sign(w)*oldWfVal.sign(w)*base_sign*new_sign;
                integralpts_imag(w, i)=0.0; 
	      }
	    }
	    
	    for(int w=0; w< nwf; w++)  {	   
              doublevar tempsum_real=0;
              doublevar tempsum_imag=0;
	      for(int l=0; l< numL(at)-1; l++) {
                for(int j=0;j<=1;j++){
		  doublevar j_actual=l-j+1.0/2.0;
		  for(int m=-l;m<=l;m++){
		    dcomplex Y_o=Ylm(l,m,olddist);
		    dcomplex Y_n=Ylm_c(l,m,newdist);
		    int temp_sign=pow(-1,j)*pow(-1,spin);
		    tempsum_real+=temp_sign/(2*l+1.0)*v_jl(j,l)*m*(Y_o.real()*Y_n.real()-Y_o.imag()*Y_n.imag());
		    tempsum_imag+=temp_sign/(2*l+1.0)*v_jl(j,l)*m*(Y_o.imag()*Y_n.real()+Y_o.real()*Y_n.imag());
		    
                  }//m	


                }//j
	      }//l
	      doublevar vxx=tempsum_real*integralpts_real(w,i)-tempsum_imag*integralpts_imag(w,i);

	      //if (vxx>=0.0)
	      nonlocal(w) += vxx; 
	    }
	    sample->setElectronPos(e, oldpos);
	  } 
	  
	  //--------------------
	  wfStore.restoreUpdate(sample, wf, e);
	}

      //----------------------------------------------
      //now do the local part
	Array1 <doublevar> vLocal(2);
        vLocal=0.0;
	int localL=numL(at)-1; //The l-value of the local part is
	//the last part.
	
        for(int j=0;j<=1;j++)
	  vLocal(j)=v_jl(j,localL);
	//accum_local+=vLocal(0);
	//accum_nonlocal+=nonlocal(0);
	
	for(int w=0; w< nwf; w++) {
 //         for(int j=0;j<=1;j++){
          totalv(e, w)+=vLocal(0)+nonlocal(w);
 //      }
        }
      } //electron loop
    }  //if atom has any psp's
    
  }  //atom loop
  
}

void Pseudopotential_so::calcNonlocWithAllvariables(Wavefunction_data * wfdata,
                                                 System * sys,
                                                 Sample_point * sample,
                                                 Wavefunction * wf,
                                                 const Array1 <doublevar> & accept_var,
                                                 Array3 <doublevar> & totalv, //TODO check dimension of this
                                                 Tmoves::tmove_type do_tmoves,vector <Tmove> & tmoves,
                                                 bool parm_derivatives, Array1 <doublevar> & parm_deriv
                                          )
{
  //Note: I left the derivative stuff commented out.
  //I don't know if it's even really correct, so beware.
  
  int natoms=sample->ionSize();
  int nwf=wf->nfunc();

  assert(accept_var.GetDim(0) >= nTest());
  assert(totalv.GetDim(0) >= nelectrons);
  assert(totalv.GetDim(1) >= nwf);
  assert(totalv.GetDim(2) >= 3);
  assert(nelectrons == sample->electronSize());

  Array1 <doublevar> ionpos(3), oldpos(3), newpos(3);
  Array1 <doublevar> newdist(5), olddist(5);
  Wf_return val(nwf,2);

  totalv=0;
  doublevar accum_local=0;
  doublevar accum_nonlocal=0;

  wf->updateVal(wfdata, sample);
  wfStore.initialize(sample, wf);
  Parm_deriv_return base_deriv;
  if(parm_derivatives) { 
    parm_deriv.Resize(wfdata->nparms());
    parm_deriv=0;
    base_deriv.need_hessian=0;
    wf->getParmDeriv(wfdata, sample, base_deriv);
  }
  int accept_counter=0;
  //deriv.Resize(natoms, 3);
  //deriv=0;
  Array2 <doublevar>  nonlocal(nwf,3);
  Array2 <doublevar> integralpts_real(nwf, maxaip);
  Array2 <doublevar> integralpts_imag(nwf, maxaip);
  Array1 <doublevar> rDotR(maxaip);
          
  
  for(int at=0; at< natoms; at++){
    if(numL(at) != 0) {
      Array2 <doublevar> v_jl(2,numL(at));
      
      sample->getIonPos(at, ionpos);

      for(int e=0; e < sample->electronSize(); e++)  {
        sample->getElectronPos(e, oldpos);
        
        //note: this updateEIDist might become inefficient..
        //depends on how often we're rejecting/how much it costs
        //to do it.  If needed, we should add an interface to 
        //Sample_point
        sample->updateEIDist();
        sample->getEIDist(e,at, olddist);
        nonlocal=0;

        int spin=1;
        if(e < sys->nelectrons(0)) spin=0;   // for spin up electrons
        getRadial(at,spin, sample, olddist, v_jl);
         
        //----------------------------------------
        //Start integral

        int accept;
        if(deterministic) {
          accept= olddist(0) < cutoff(at);
        }
        else {
          doublevar strength=0;
          const doublevar calculate_threshold=10;
  
          for(int l=0; l<numL(at)-1; l++)
            for(int j=0;j<=1;j++)
              strength+=calculate_threshold*(l-j+1)*fabs(v_jl(j,l));
   
          strength=min((doublevar) 1.0, strength);
  
          doublevar rand=accept_var(accept_counter++);
        //
          if(fabs(strength)> 0)
            for(int l=0; l<numL(at)-1; l++)
              for(int j=0;j<=1;j++)
                v_jl(j,l)/=strength;
          accept=strength>rand;
        }
      

        if(accept)  {
          wfStore.saveUpdate(sample, wf, e);
          Wf_return  oldWfVal(nwf,2);
          wf->getVal(wfdata, e,oldWfVal);  //ignore the index e here.

          for(int i=0; i< aip(at); i++) {
            sample->setElectronPos(e, oldpos);
            doublevar base_sign=sample->overallSign();
            doublevar base_phase=sample->overallPhase();

            //Make sure to move the electron relative to the nearest neighbor
            //in a periodic calculation(so subtract the distance rather than
            //adding to the ionic position).  This actually only matters 
            //when we're doing non-zero k-points.
            for(int d=0; d < 3; d++) 
              newpos(d)=integralpt(at,i,d)*olddist(0)-olddist(d+2);
            
            sample->translateElectron(e, newpos); 
            sample->updateEIDist(); 
            sample->getEIDist(e,at,newdist);


            rDotR(i)=0;
            for(int d=0; d < 3; d++)
              rDotR(i)+=newdist(d+2)*olddist(d+2);
            doublevar new_sign=sample->overallSign();
            doublevar new_phase=sample->overallPhase();
            
            rDotR(i)/=(newdist(0)*olddist(0));  //divide by the magnitudes
            

 
            wf->updateVal(wfdata, sample);
            wf->getVal(wfdata, e, val); 
            //----
            for(int w=0; w< nwf; w++) {
              integralpts_real(w,i)=exp(val.amp(w,0)-oldWfVal.amp(w,0))
                *integralweight(at, i);
              integralpts_imag(w,i)=integralpts_real(w,i);
	      //cout << exp(val.amp(w,0)) << endl;
              if ( val.is_complex==1 ) {
                integralpts_real(w,i)*=cos(val.phase(w,0)+new_phase
                    -oldWfVal.phase(w,0)-base_phase);
                integralpts_imag(w,i)*=sin(val.phase(w,0)+new_phase
                    -oldWfVal.phase(w,0)-base_phase);
                //cout << val.phase(w,0)+new_phase-oldWfVal.phase(w,0)-base_phase << endl;
	      } else {
                integralpts_real(w,i)*=val.sign(w)*oldWfVal.sign(w)
                  *base_sign*new_sign;
                integralpts_imag(w,i)=0.00;
              }
            }
            
            for(int w=0; w< nwf; w++)  {
              doublevar tempsum_real_x=0;
              doublevar tempsum_real_y=0;
              doublevar tempsum_real_z=0;
          //    doublevar tempsum_real_AREP=0;
              doublevar tempsum_imag_x=0;
              doublevar tempsum_imag_y=0;
              doublevar tempsum_imag_z=0;
          //    doublevar tempsum_imag_AREP=0;
              for(int l=0; l< numL(at)-1; l++) {
                for(int j=0;j<=1;j++){
                  doublevar j_actual=l-j+1.0/2.0;

       	
		// new method l \dot s
		    //for(doublevar mj=-j_actual;mj<=j_actual;mj++){
                    //    int m;
		    //	m=round(mj+spin-1.0/2.0); 
		    for(int m=-l;m<=l;m++){
			dcomplex Y_o=Ylm(l,m,olddist);
                        for(int m_p=-l;m_p<=l;m_p++){
			    dcomplex Y_n=Ylm_c(l,m_p,newdist);
		            int temp_sign=pow(-1,j)*pow(-1,spin); //TODO: check this term
		            //int temp_sign=pow(-1,j);  // didn't take the spin direction into account.
                            dcomplex int_ls_x = integral_ls_x(l,m,m_p);
                            dcomplex int_ls_y = integral_ls_y(l,m,m_p);
                            dcomplex int_ls_z = integral_ls_z(l,m,m_p);
                         //   cout << "l=" << l << " m=" << m << " m_p=" << m_p << " ";
			 //   cout << "ls_x " << int_ls_x << " ";
			 //   cout << "ls_y " << int_ls_y << " ";
			 //   cout << "ls_z " << int_ls_z << " ";
 			    
  			    tempsum_real_x+=temp_sign/(2*l+1.0)*v_jl(j,l)*(int_ls_x.real()
                                            *(Y_o.real()*Y_n.real()-Y_o.imag()*Y_n.imag())
					    -int_ls_x.imag()
					    *(Y_o.imag()*Y_n.real()+Y_o.real()*Y_n.imag()));
			    tempsum_imag_x+=temp_sign/(2*l+1.0)*v_jl(j,l)*(int_ls_x.real()
					    *(Y_o.imag()*Y_n.real()+Y_o.real()*Y_n.imag())
					    +int_ls_x.imag()
					    *(Y_o.real()*Y_n.real()-Y_o.imag()*Y_n.imag()));
  			    tempsum_real_y+=temp_sign/(2*l+1.0)*v_jl(j,l)*(int_ls_y.real()
                                            *(Y_o.real()*Y_n.real()-Y_o.imag()*Y_n.imag())
					    -int_ls_y.imag()
					    *(Y_o.imag()*Y_n.real()+Y_o.real()*Y_n.imag()));
			    tempsum_imag_y+=temp_sign/(2*l+1.0)*v_jl(j,l)*(int_ls_y.real()
					    *(Y_o.imag()*Y_n.real()+Y_o.real()*Y_n.imag())
					    +int_ls_y.imag()
					    *(Y_o.real()*Y_n.real()-Y_o.imag()*Y_n.imag()));
  			    tempsum_real_z+=temp_sign/(2*l+1.0)*v_jl(j,l)*(int_ls_z.real()
                                            *(Y_o.real()*Y_n.real()-Y_o.imag()*Y_n.imag())
					    -int_ls_z.imag()
					    *(Y_o.imag()*Y_n.real()+Y_o.real()*Y_n.imag()));
			    tempsum_imag_z+=temp_sign/(2*l+1.0)*v_jl(j,l)*(int_ls_z.real()
					    *(Y_o.imag()*Y_n.real()+Y_o.real()*Y_n.imag())
					    +int_ls_z.imag()
					    *(Y_o.real()*Y_n.real()-Y_o.imag()*Y_n.imag()));

                        } //m_p
		    }  //m
	 	
		 
	/*    
		    for(doublevar mj=-j_actual;mj<=j_actual;mj++){
 			int m=round(mj+spin-1.0/2.0);
			dcomplex Y_o=Ylm(l,m,olddist);
			dcomplex Y_n=Ylm_c(l,m,newdist);
		        doublevar temp_coeff;
			if (((spin==0)&&(j==0))||((spin==1)&&(j==1))){
			  temp_coeff=l+1.0/2.0+mj;
			}	
			else {temp_coeff=l+1.0/2.0-mj;}
			tempsum_real+=temp_coeff/(2*l+1.0)*v_jl(j,l)*(Y_o.real()*Y_n.real()
				      -Y_o.imag()*Y_n.imag());
			tempsum_imag+=temp_coeff/(2*l+1.0)*v_jl(j,l)*(Y_o.imag()*Y_n.real()
				      +Y_o.real()*Y_n.imag());  	 
		    }  //m
	 


		// AREP
             
                    for (int ml=-l;ml<=l;ml++){
                      dcomplex Y_o_AREP=Ylm(l,ml,olddist);
                      dcomplex Y_n_AREP=Ylm_c(l,ml,newdist);
   		      tempsum_real_AREP+=(l+1-j)/(2*l+1.0)*v_jl(j,l)*(Y_o_AREP.real()*Y_n_AREP.real()-Y_o_AREP.imag()*Y_n_AREP.imag());
                      tempsum_imag_AREP+=(l+1-j)/(2*l+1.0)*v_jl(j,l)*(Y_o_AREP.real()*Y_n_AREP.imag()+Y_o_AREP.imag()*Y_n_AREP.real());
                    }  //ml
          */   
		    //tempsum_real+=(l)*v_jl(j,l)*legendre(rDotR(i), l);  
                } //j
              }  //l
      
              doublevar vxx_x=tempsum_real_x*integralpts_real(w,i)-tempsum_imag_x*integralpts_imag(w,i);
              doublevar vxx_y=tempsum_real_y*integralpts_real(w,i)-tempsum_imag_y*integralpts_imag(w,i);
              doublevar vxx_z=tempsum_real_z*integralpts_real(w,i)-tempsum_imag_z*integralpts_imag(w,i);
	//		   -(tempsum_real_AREP*integralpts_real(w,i)-tempsum_imag_AREP*integralpts_imag(w,i));
            //  cout << "vxx=" << vxx << endl;
              doublevar vxx_imag_x=tempsum_imag_x*integralpts_real(w,i)+tempsum_real_x*integralpts_imag(w,i);
              doublevar vxx_imag_y=tempsum_imag_y*integralpts_real(w,i)+tempsum_real_y*integralpts_imag(w,i);
              doublevar vxx_imag_z=tempsum_imag_z*integralpts_real(w,i)+tempsum_real_z*integralpts_imag(w,i);
	//	             +tempsum_real_AREP*integralpts_imag(w,i)+tempsum_imag_AREP*integralpts_real(w,i);

              if(do_tmoves==Tmoves::no_tmove || w!=0 || (do_tmoves==Tmoves::negative_tmove && vxx_x >= 0.0 && vxx_y >=0.0 && vxx_z >= 0.0) ) { 
                nonlocal(w,0)+=vxx_x;
                nonlocal(w,1)+=vxx_y;
                nonlocal(w,2)+=vxx_z;
                //cout << "no tmoves" << endl;
              }
              else { 
                Tmove nwtmove; nwtmove.pos.Resize(3);
                //sample->getElectronPos(e,nwtmove.pos);
                nwtmove.pos=newpos;
                nwtmove.e=e;
                nwtmove.vxx=vxx_z;  //TODO check this
                tmoves.push_back(nwtmove);
                //cout << "do tmoves" << endl;
              }
              
              //-----------parameter derivatives
              if(parm_derivatives) { 
                Parm_deriv_return deriv;
                deriv.need_hessian=0;
                wf->getParmDeriv(wfdata, sample, deriv);
                int np=wfdata->nparms();
                for(int p=0; p < np; p++) { 
                  parm_deriv(p)+=(deriv.gradient(p)-base_deriv.gradient(p))*vxx_x; //TODO: this is probably wrong
                }
              }
              //------
            }
            sample->setElectronPos(e, oldpos);
          } //aip

          //--------------------



          wfStore.restoreUpdate(sample, wf, e);
        }

        //----------------------------------------------
        //now do the local part
        Array1 <doublevar> vLocal(2);
        vLocal=0.0;
        int localL=numL(at)-1; 
                               

        for(int j=0;j<=1;j++)//{
          vLocal(j)=v_jl(j,localL);
          //accum_local+=vLocal(j);}
       
        //accum_nonlocal+=nonlocal(w);

        for(int w=0; w< nwf; w++) {
          totalv(e,w,0)+=nonlocal(w,0);   
          totalv(e,w,1)+=nonlocal(w,1);   
          totalv(e,w,2)+=nonlocal(w,2);   
        }



      }  //electron loop

    }  //if atom has any psp's

  }  //atom loop

  //cout << "psp: local part " << accum_local
  // << "  nonlocal part " << accum_nonlocal << endl;

}

//------------------------------------------------------------------------


void Pseudopotential_so::randomize() {
  Array1 <doublevar> x(3), y(3), z(3);
  generate_random_rotation(x,y,z);
  rotateQuadrature(x,y,z);
}


//----------------------------------------------------------------------
void Pseudopotential_so::rotateQuadrature(Array1 <doublevar> & x,
                                       Array1 <doublevar> & y,
                                       Array1 <doublevar> & z)
{


  int natoms=aip.GetDim(0);

  //cout << "x1, x2, x3" << x1 << "   " << x2 << "   " << x3 << endl;
  //cout << "y1, y2, y3" << y1 << "   " << y2 << "   " << y3 << endl;
  //cout << "z1, z2, z3" << z1 << "   " << z2 << "   " << z3 << endl;

  for(int at=0; at<natoms; at++)
  {
    //cout << "-----------atom   " << at << endl;
    for(int i=0; i< aip(at); i++)
    {
      //cout << "quadrature points before:\n"
      //   << integralpt(at, i, 0) << "   "
      //   <<integralpt(at, i, 1) << "    "
      //   <<integralpt(at, i, 2) << "    \n";
      integralpt(at,i,0)=integralpt_orig(at,i,0)*x(0)
                         +integralpt_orig(at,i,1)*y(0)
                         +integralpt_orig(at,i,2)*z(0);
      integralpt(at,i,1)=integralpt_orig(at,i,0)*x(1)
                         +integralpt_orig(at,i,1)*y(1)
                         +integralpt_orig(at,i,2)*z(1);
      integralpt(at,i,2)=integralpt_orig(at,i,0)*x(2)
                         +integralpt_orig(at,i,1)*y(2)
                         +integralpt_orig(at,i,2)*z(2);
      //cout << "quadrature points after:\n"
      //   << integralpt(at, i, 0) << "   "
      //	   <<integralpt(at, i, 1) << "    "
      //   <<integralpt(at, i, 2) << "    \n";
    }
  }
}



//------------------------------------------------------------------------

#include "System.h"

void Pseudopotential_so::read(vector <vector <string> > & pseudotext, System * sys){

  //cout << "trigers reading pspso" << endl;
  sys->getAtomicLabels(atomnames);
  //cout << atomnames[0] << endl;
  //cout << "after getting atomnames" << endl;
  int natoms=atomnames.size();
  nelectrons=sys->nelectrons(0)+sys->nelectrons(1); //up+down
  //cout << "after geting nelec" << endl;
  aip.Resize(natoms);
  aip=6;
  //cout << "after reading aip" << endl;
  integralpt.Resize(natoms,maxaip,3); 
  integralpt_orig.Resize(natoms,maxaip,3);
  integralweight.Resize(natoms,maxaip);
  addzeff.Resize(natoms);
  addzeff=false;
  Array1 <int> atom_has_psp(natoms);
  atom_has_psp=0;
  int maxL=0;

  //cout << "begin reading pspso" << endl;

  radial_basis.Resize(natoms,2,2);
  radial_basis=NULL;

  //cout << "psptxtsize=" << pseudotext.size() << endl;

  for(unsigned int i=0; i<pseudotext.size(); i++){ //i=j
    for (int at=0; at<natoms; at++) {
      unsigned int pos=0;
      if( pseudotext[i][0] == atomnames[at] ){
        atom_has_psp(at)=1;
        vector <string> basistmp;
        if(!readsection(pseudotext[i], pos=0, basistmp, "BASIS"))
            error("Need Basis section in pseudopotential for ", pseudotext[i][0]);
        allocate(basistmp, radial_basis(at,0,i%2)); 
        int basistmpSize=basistmp.size();
      //  for (int ibasistmp=0;ibasistmp<basistmpSize;ibasistmp++){
       //   cout << "basistmp[" << ibasistmp <<"]=" << basistmp[ibasistmp] << endl;}

        allocate(basistmp, radial_basis(at,1,i%2));
      
        assert(radial_basis(at,0,0)->nfunc()==radial_basis(at,1,0)->nfunc());
        int nlval=radial_basis(at,0,0)->nfunc();
        if(maxL < nlval) maxL=nlval;
        pos=0;
        readvalue(pseudotext[i],pos,aip(at),"AIP");
        if(haskeyword(pseudotext[i],pos=0,"ADD_ZEFF")){addzeff(at)=true;}
      }
    }
  }

  for(int at=0;at<natoms;at++) {
    int aiptemp=aip(at);
    //cout << "aiptemp=" << aiptemp;
    //cout << "maxaip=" << maxaip;
    Array1 <doublevar> xpt(maxaip);
    Array1 <doublevar> ypt(maxaip);
    Array1 <doublevar> zpt(maxaip);
    Array1 <doublevar> weight(maxaip);

    gesqua(aiptemp, xpt, ypt, zpt, weight);

    aip(at)=aiptemp;
    for(int i=0;i<aip(at);i++){
      integralpt(at,i,0)=integralpt_orig(at,i,0)=xpt(i);
      integralpt(at,i,1)=integralpt_orig(at,i,1)=ypt(i);
      integralpt(at,i,2)=integralpt_orig(at,i,2)=zpt(i);
      integralweight(at,i)=weight(i);
    }
  }

  numL.Resize(natoms);
  numL=0;
  for(int at=0;at<natoms;at++){
    if(radial_basis(at,0,0)!=NULL)
      numL(at)=radial_basis(at,0,0)->nfunc();
  }

  Sample_point * tempsample=NULL;
  sys->generateSample(tempsample);
  cutoff.Resize(natoms);
  const doublevar cutoff_threshold=1e-5;
  const doublevar cutoff_max=20.0;
  cutoff=0.0;
  const doublevar cutoff_interval=.05;
  for(int at=0; at< natoms; at++)
  {
    if(numL(at)>0) {
      Array1 <doublevar> cutoffL(numL(at));
      Array2 <doublevar> v_jl(2,numL(at)); //spin up
      Array2 <doublevar> v_jl2(2,numL(at)); //spin down
      Array1 <doublevar> tempr(5);//r,r*r,x,y,z
      Array1 <int> foundcutoff(numL(at));
      foundcutoff=0;
      cutoffL=cutoff_max;
      for(doublevar r=cutoff_max;r>0;r-=cutoff_interval){
        tempr(0)=r; tempr(1)=r*r; tempr(4)=r;
        getRadial(at,0,tempsample,tempr,v_jl);
        getRadial(at,1,tempsample,tempr,v_jl);

        for(int l=0; l< numL(at)-1; l++)
        {
          if( (fabs(v_jl(0,l))>cutoff_threshold || fabs(v_jl2(0,l))>cutoff_threshold) && !foundcutoff(l) ){
            cutoffL(l)=r-cutoff_interval;
            foundcutoff(l)=1;
          }
        }
      }

      for(int l=1; l<numL(at)-1;l++){
        cutoff(at)=max(cutoff(at),cutoffL(l));
      }
    }
  }


  //cout << "end of reading psp_so" << endl;


}

//----------------------------------------------------------------------

int Pseudopotential_so::showinfo(ostream & os)
{
  os << "Pseudopotential_so " << endl;
  string indent="  ";
  vector <string> uniquenames;
  int natoms=aip.GetDim(0);
  for(int at=0; at < natoms; at++) {
    int unique=1;
    for(unsigned int i=0; i< uniquenames.size(); i++) {
      if(uniquenames[i]==atomnames[at]) {
        unique=0;
        break;
      }
    }
    if(unique) {
      uniquenames.push_back(atomnames[at]);
      
      os << "atom " << atomnames[at] << endl;
      if(numL(at)==0) { os << "No pseudopotential" << endl; }
      else {
        os << "Integration points " << aip(at) << endl;
//os << setw(10) << "x" 
//   << setw(10) << "y" 
//   << setw(10) << "z"
//   << setw(10) << "weight" << endl;
//	for(int i=0; i< aip(at); i++)
//	  {
//	    os << setw(10) << integralpt(at,i,0)
//	       << setw(10) << integralpt(at,i,1)
//	       << setw(10) << integralpt(at,i,2)
//	       << setw(10)  << integralweight(at, i) << endl;
//	  }	
        os << "Cutoff for static calculation "<<cutoff(at)<<endl;
        os << "Pseudopotential for spin up,j=l+1/2: \n";
        radial_basis(at,0,0)->showinfo(indent, os);
        os << "Pseudopotential for spin down,j=l+1/2: \n";
        radial_basis(at,1,0)->showinfo(indent, os);
        os << "Pseudopotential for spin up,j=l-1/2: \n";
        radial_basis(at,0,1)->showinfo(indent, os);
        os << "Pseudopotential for spin down,j=l-1/2: \n";
        radial_basis(at,1,1)->showinfo(indent, os);
      }
    }
  }

  return 1;
}


//------------------------------------------------------------------------
