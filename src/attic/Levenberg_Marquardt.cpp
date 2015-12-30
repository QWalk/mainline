/*
 
Copyright (C) 2007 Michal Bajdich

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
#include "Optimize_method2.h"
#include "Array.h"
#include "MatrixAlgebra.h"

void Optimize_method2::wf_printout(int iter, doublevar value, doublevar energy, doublevar variance, int min_nconfig, 
				   doublevar mu, ostream & output){
  if (output){
  string indent="";
  ofstream wfoutput(wfoutputfile.c_str());
  wfoutput.precision(15);
  wfdata->writeinput(indent,wfoutput);
  wfoutput.close();
  
  if (use_extended_output){
    char strbuff[40];
    sprintf(strbuff, "%d", iter);
    string wfoutputfile2=wfoutputfile+"_"+strbuff;
    ofstream wfoutput2(wfoutputfile2.c_str());
    wfoutput2.precision(15);
    wfdata->writeinput(indent,wfoutput2);
    wfoutput2.close();
  }
  output << "iteration # " << iter;
  switch(min_function)
    {
    case min_variance:
      output <<" dispersion= " << sqrt(fabs(value)) << " energy= " << energy 
             <<" nconfig "<< min_nconfig << "  damping "<<mu<< endl;
      break;
    case min_energy:
      output <<" energy= " << value <<" dispersion= "<< sqrt(variance) 
             <<" nconfig "<<min_nconfig << " damping "<<mu <<endl;
      break;
    case min_mixed:
      output <<" mixing*energy+ (1-mixing)*variance= " << value 
             <<" energy= " << energy <<" dispersion= "<< sqrt(variance)
             <<" nconfig "<<min_nconfig << "  damping "<<mu <<endl;
      break;  
    default:
      error("Optimize_method::variance() : min_function has a very strange value");
    }
  }
}

void make_posite_definite2(Array2 <doublevar> & function_hessian, ostream & output){
  int dim=function_hessian.GetDim(0);
  Array1 <doublevar> eigenvals(dim);
  Array2 <doublevar> eigenvecs(dim,dim);
  EigenSystemSolverRealSymmetricMatrix(function_hessian, eigenvals, eigenvecs);
  if(output){
    cout << "EigenValues of Hessian are: "<<endl;
    for(int i=0; i<dim; ++i){
      cout << i+1<<": "<<eigenvals(i)<<endl;
    }
  }
  //if(eigenvals(dim-1)<0){
  // for(int i=0; i<dim; ++i)
  //  function_hessian(i,i)-=eigenvals(dim-1);
  //}
}


int Optimize_method2::LEVMAR_DER(Array1 <double> & parms, int nparms_start, int nparms_end,  
                                 Array1 <double> & delta,
                                 int & min_nconfig, int iter_min, 
                                 const int itmax, int & iter, ostream & output)
{
  doublevar EPSILON=1e-12;
  doublevar ONE_THIRD=0.3333333334; 
  doublevar LM_INIT_MU=1e-01;
  doublevar LM_STOP_THRESH=1e-6;
  doublevar LM_STOP_THRESH2=1e-10;
  
  int m=parms.GetSize();
  int m_part=nparms_end-nparms_start;

  Array1 <doublevar> grad(m_part);    //gradient 
  Array2 <doublevar> hessian(m_part,m_part), hessian_inverse(m_part,m_part); //hessian and its inverse
  Array1 <doublevar> Dp(m_part);       //change of parms 
  Array1 <doublevar> diag_hess(m_part);  //diagonal of hessian 
  Array1 <doublevar> pDp(m_part);       //partial new parms 
  Array1 <doublevar> pDp_full(m);       //full new parms 
 

  doublevar mu=0.0;  //damping constant 
  doublevar tmp; // mainly used in matrix & vector multiplications 
  doublevar grad_inf;// 
  doublevar p_L2, Dp_L2, dF, dL;
  doublevar tau, eps1, eps2, eps2_sq;
  int nu=2, nu2, stop;

  doublevar value, energy, variance;
     
  tau=LM_INIT_MU;
  eps1=LM_STOP_THRESH;
  eps2=LM_STOP_THRESH2;
  eps2_sq=eps2*eps2;
  stop=0;
  iter=0;
  int startfull=0;
  int counter1=0;
  doublevar variance_old;

  func_val(m, parms, value, energy, variance, min_nconfig, output);
  variance_old=variance;

  if(use_weights)
    eref=energy;

  wf_printout(iter, value, energy, variance, min_nconfig, mu, output);
    
  for(int k=0; k< itmax && !stop; ++k){
    if(analytic_wf_ders)
      func_hessian_analytical(parms, nparms_start, nparms_end, hessian, grad, energy, delta, min_nconfig, output);
    else
      func_hessian(parms, nparms_start, nparms_end, hessian, grad, energy, delta, min_nconfig, output);
      
    make_posite_definite2(hessian, output);

    p_L2=grad_inf=0.0; 
    for(int i=0; i<m_part; ++i){
      if(grad_inf < (tmp=abs(grad(i)))) 
	grad_inf=tmp;
      diag_hess(i)=hessian(i,i); 
      // save diagonal entries so that augmentation can be later canceled 
      p_L2+=parms(i+nparms_start)*parms(i+nparms_start);
    }
    //p_L2=sqrt(p_L2);

    if(!(k%1)){
      if (output) {
	cout << "Current estimate of parms: "<<endl;;
	for(int i=0; i<m; ++i){
	  cout << parms(i)<<"  ";
	  if(i%6==5)
	    cout << endl;
	}
	cout <<endl;
	cout << "size of the largest component of gradient: "<<grad_inf<<"  value: "<<value<<endl;;
      }
    }
    
    // check for convergence 
    if (output)
      cout << "check for convergence"<<endl;
    if((grad_inf <= eps1)){
      Dp_L2=0.0;  //no increment for p in this case 
      stop=1;
      cout <<"exit"<<endl;
      break;
    }

    // compute initial damping factor 
    if(k==0){
      tmp=0.0;
      //double tmp2=0.0;
      for(int i=0; i<m_part; ++i){
        if(diag_hess(i)>tmp) tmp=diag_hess(i);  //find max diagonal element 
	// if(diag_jacTjac[i]<tmp2) tmp2=diag_jacTjac[i];  find min diagonal element 
      }
      mu=tau*tmp;//-tmp2;
      //mu=0.2;
    }

    // determine increment using adaptive damping 
    if (output)
      cout << "determine increment using adaptive damping"<<endl;
    while(1){
      if (output)
	cout <<"Damping factor mu: "<<mu<<endl;
      for(int i=0; i<m_part; ++i)
        hessian(i,i)+=mu;

      
      // find new parametres by Newton method using inverse of hessian
      InvertMatrix(hessian, hessian_inverse, m_part);
      Dp_L2=0.0; 
      for (int i=0;i<m_part;i++){
	tmp=0;
	for (int j=0;j<m_part;j++)
	  tmp-=hessian_inverse(i,j)*grad(j);
	Dp(i)=tmp;
	pDp[i]=parms[i+nparms_start]+tmp;
	Dp_L2+=tmp*tmp;
      }
      
      for (int i=0;i<m;i++){
	if(i>= nparms_start && i<nparms_end)
	  pDp_full[i]=pDp[i-nparms_start];
	else
	  pDp_full[i]=parms[i];
      }
      
      

      if(Dp_L2<=eps2_sq*p_L2){ //  relative change in p is small, stop 
      //if(Dp_L2<=eps2*(p_L2 + eps2)){  // relative change in p is small, stop 
	if (output)
	  cout <<"relative change in p is small, stop"<<endl;
	stop=2;
	break;
      }

      if(Dp_L2>=(p_L2+eps2)/(EPSILON*EPSILON)){ // almost singular 
	if (output)
	  cout <<"almost singular, stop"<<endl;
	stop=4;
	break;
      }
      // evaluate function at p + Dp 
      doublevar new_value,new_energy,new_variance;
      
      func_val(m, pDp_full, new_value,new_energy,new_variance,min_nconfig, output);
      
      if (output)
	cout << "new value is " <<new_value<<endl;
      
      dL=0.0;
      for(int i=0; i<m_part; ++i)
	dL+=Dp(i)*(mu*Dp(i)-grad(i));
      dF=value-new_value;
      if(output)
	cout << "dL "<<dL<<" and dF "<<dF<<endl;
      if(dL>0.0 && dF>0.0 ){ // reduction in error, increment is accepted 
      	 if (output)
	   cout << "reduction in error, increment is accepted"<<endl;  

	 value=new_value;
         energy=new_energy;
         if(use_weights)
           eref=energy;
         variance=new_variance;
         for(int i=0 ; i<m; ++i)  //update p's estimate 
           parms[i]=pDp_full[i];
	 
	 wf_printout(iter+1, value, energy, variance, min_nconfig, mu, output);
         if(sqrt(variance) > 1.1*sqrt(variance_old))
           output <<"  WARNING: dispersion is growing too fast! Using too few/old configurations ???"<<endl;
         variance_old=variance;
           

	 tmp=(4.0*dF/dL-1.0);
         tmp=1.0-tmp*tmp*tmp;
         mu=mu*( (tmp>=ONE_THIRD)? tmp : ONE_THIRD );
         nu=2;
	 break;
      }
      else 
	if (  startfull<0 && abs(dF)< abs(value*0.001) ){
	  if (output)
	   cout << "no reduction in error, but increment is accepted, keeping the same damping"<<endl;
	  value=new_value;
          energy=new_energy;
          if(use_weights)
            eref=energy;
          variance=new_variance;
	  for(int i=0 ; i<m; ++i)  //update p's estimate 
	    parms[i]=pDp_full[i];
	  wf_printout(iter+1, value, energy, variance, min_nconfig, mu, output);
          if(sqrt(variance) > 1.1*sqrt(variance_old))
            output <<"  WARNING: dispersion is growing too fast! Using too few/old configurations ???"<<endl;
          variance_old=variance;
	  break;
	  }
     
      

      // the error did not reduce; in any case, the increment must be rejected
      
      mu*=nu;
      nu2=nu<<1; // 2*nu;
      if(nu2<=nu){ // nu has wrapped around (overflown). 
	//Thanks to Frank Jordan for spotting this case 
	cout << "nu has wrapped around (overflown)"<<endl;
        stop=5;
        break;
      }
      nu=nu2;

      for(int i=0; i<m_part; ++i) // restore diagonal hessian entries 
        hessian(i,i)=diag_hess(i);
    } // inner loop of conditional 
    // if (!stop)
    // if (output)
    //cout <<"Damping parameter mu=           "<<mu<<endl;
    //cout <<"iter_min           "<<iter_min<<endl;
    

    //min_nconfig= should change here
    if (counter1+1 == iter_min && min_nconfig < nconfig){
      if (min_nconfig*multiply < nconfig){
	min_nconfig*=multiply;
      }
      else 
	min_nconfig=nconfig;
      counter1=0;
    }
    else {
      counter1++;
    }

    if ( min_nconfig==nconfig )
      startfull++;

    iter++;
   }// loop over iterations

  //cout << iter<<" "<<itmax<<endl;
  if(iter>=itmax) {
    stop=3;
    if (output)
    cout << "Too many iterations in Levenberg-Marquardt routine"<<endl;
  }
  
  if (stop!=4 || stop!=5)
    return 1;
  else 
    return -1;
}


