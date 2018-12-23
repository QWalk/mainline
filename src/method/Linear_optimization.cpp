#include "ulec.h"
#include "Wavefunction_data.h"
#include "Sample_point.h"
#include "Guiding_function.h"
#include "System.h"
#include "Pseudopotential.h"
#include "Linear_optimization.h"
#include "qmc_io.h"
#include "Program_options.h"
#include "MatrixAlgebra.h"
#include "Vmc_method.h"
#include "Properties_average.h"
#include "Properties_block.h"
#include <string>
#include <algorithm>

void Linear_optimization_method::read(vector <string> words,
            unsigned int & pos, Program_options & options_) { 

  options=options_;

  int total_nstep,total_fit;
  if(!readvalue(words,pos=0,total_nstep,"TOTAL_NSTEP"))
    total_nstep=16384;
  if(!readvalue(words,pos=0,total_fit,"TOTAL_FIT"))
    total_fit=1024;

  if(!readvalue(words,pos=0, iterations, "ITERATIONS")) 
    iterations=30;
  if(! readvalue(words, pos=0, wfoutputfile, "WFOUTPUT") )
      wfoutputfile=options.runid+".wfout";
  if(!readvalue(words,pos=0, vmc_nstep,"VMC_NSTEP"))
    vmc_nstep=max(total_nstep/mpi_info.nprocs,10);
  if(!readvalue(words,pos=0,nconfig_eval,"FIT_NCONFIG")) 
    nconfig_eval=max(total_fit/mpi_info.nprocs,1);
  if(!readvalue(words,pos=0,en_convergence,"EN_CONVERGENCE"))
    en_convergence=0.001;
  if(!readvalue(words,pos=0,sig_H_threshold,"SIG_H_THRESHOLD"))
    sig_H_threshold=0.5;
  if(!readvalue(words,pos=0,minimum_psi0,"MINIMUM_PSI0"))
    minimum_psi0=0.8;
  if(!readvalue(words, pos=0, max_vmc_nstep, "MAX_VMC_NSTEP"))
    max_vmc_nstep=2*vmc_nstep;
  if(!readvalue(words, pos=0, max_nconfig_eval, "MAX_FIT_NCONFIG"))
    max_nconfig_eval=4*nconfig_eval;
  do_uncorrelated_evaluation=false;
  if(haskeyword(words, pos=0, "UNCORRELATED_EVALUATION"))
    do_uncorrelated_evaluation=true;
  pseudopotential_derivatives=haskeyword(words,pos=0,"PSEUDOPOTENTIAL_DERIVATIVES");
  if(!readvalue(words, pos=0, max_nconfig_eval, "MAX_ZERO_ITERATIONS"))
    max_zero_iterations=2;
  if(!readvalue(words,pos=0,nodal_cutoff, "NODAL_CUTOFF"))
    nodal_cutoff=0;

  if(!readvalue(words, pos=0, max_vmc_nstep, "SVD_TOLERANCE"))
    svd_tolerance=-1;
  
  allocate(options.systemtext[0],  sys);
  sys->generatePseudo(options.pseudotext, pseudo);
  wfdata=NULL;
  allocate(options.twftext[0], sys, wfdata);

  if(wfdata->nparms() <= 0 ) 
    error("There appear to be no parameters to optimize!");

  
}


//----------------------------------------------------------------------

int Linear_optimization_method::showinfo(ostream & os) { 
  os << "     System " << endl;
  sys->showinfo(os);
  os << endl << endl;
  os << "     Wavefunction " << endl;
  wfdata->showinfo(os);
  os << endl << endl;
  pseudo->showinfo(os);
  os << endl << endl;
  os << "-----------------------------" << endl;
  os << "Linear wave function optimization:  " << endl;
  os << "Number of processors: " << mpi_info.nprocs << endl;
  os << "Number of MC steps : " << vmc_nstep*mpi_info.nprocs << endl;
  os << "Number of correlated sampling walkers : " << nconfig_eval*mpi_info.nprocs << endl;
  os << "Wave function output to file  : " << wfoutputfile << endl;
  if(do_uncorrelated_evaluation) 
    os << "Uncorrelated evaluation of energy differences using above number of walkers" << endl;
  if(pseudopotential_derivatives) 
    os << "Computing derivatives using the pseudopotential" << endl;
  os << "nparms " << wfdata->nparms() << endl;
  os << "----------------------------" << endl;
  return 1;
}

//----------------------------------------------------------------------

void Linear_optimization_method::run(Program_options & options, ostream & output) { 
  int nparms=wfdata->nparms();
  Array1 <doublevar> x(nparms+1); 
  Array2 <doublevar> S(nparms+1,nparms+1),Sinv(nparms+1,nparms+1);
  Array2 <doublevar> H(nparms+1,nparms+1);
  Array1 <doublevar> alpha(nparms); 
 
  if(!wfdata->supports(parameter_derivatives))
    error("Wavefunction needs to supports analytic parameter derivatives");
  
  wfdata->getVarParms(alpha); 
  
  Array1 <doublevar> en,tmp_en;
  Array2 <doublevar> energy_step(iterations,2);
  Array2 <doublevar> alpha_step(iterations,nparms);
  
  wfdata->getVarParms(alpha);
  int nzero_iterations=0;


  for(int it=0; it< iterations; it++) {
    //cout<< "wf derivative" <<endl;
    Array1 <doublevar> olden=en;
    
    wavefunction_derivative(H,S,en);
    output << "####################\n";
    output << "step " << it << ": current energy " << setprecision(10) << en(0) 
           << " +/- " << setprecision(10) << en(1) << endl;
    output.flush();

    /* Debug output
    if(mpi_info.node==0) {
    stringstream outfname; 
    outfname << "linear_matrices" << it;
    ofstream matrixout(outfname.str().c_str());
    int n=H.GetDim(0);
    for(int i=0; i< n; i++) { 
      for(int j=0; j < n; j++) { 
        matrixout << H(i,j) << " ";
      }
      matrixout << endl;
    }
    for(int i=0; i< n; i++) { 
      for(int j=0; j < n; j++) { 
        matrixout << S(i,j) << " ";
      }
      matrixout << endl;
    }
    matrixout.close();
    }
    */
    

    //if(it>0)
    //  output << "  energy change " << en(0)-olden(0) 
    //    << " +/- " << sqrt(en(1)*en(1)+olden(1)*olden(1)) << endl;
    
    doublevar endiff= line_minimization(S,Sinv,H,alpha,output);
    if(endiff >= 0 && vmc_nstep < max_vmc_nstep) { 
      vmc_nstep*=4;
      output << "Did not find a downhill move; increasing total vmc steps to "
        << vmc_nstep*mpi_info.nprocs << endl;
    }
    else if(endiff >= 0) {
      nzero_iterations++;
      output << "Iterations without a downhill move:" << nzero_iterations << endl;
    }
    wfdata->setVarParms(alpha);
    wfdata->lockInParms(alpha);
    wfdata->renormalize();
   
    if(mpi_info.node==0) { 
      string indentation="";
      ofstream wfoutput(wfoutputfile.c_str());
      wfoutput.precision(15);
      wfdata->writeinput(indentation,wfoutput);
      wfoutput.close();
    }
   

    for(int i=0; i< nparms; i++) {
      alpha_step(it,i)=alpha(i);
    }
    if(nzero_iterations >= max_zero_iterations) {
      output << "Reached " << max_zero_iterations
             << " without a downhill move and at the maximum number of samples.\n"
             << "Halting." << endl;
      break;
    }
  }


#ifdef USE_MPI
  MPI_Barrier(MPI_Comm_grp);
#endif
  

  ifstream wfinput(wfoutputfile.c_str());
  options.twftext[0].clear();
  parsefile(wfinput, options.twftext[0]);
  wfinput.close();

}

//----------------------------------------------------------------------

//return the proportion of the move that was Psi_0
doublevar find_directions(Array2 <doublevar> & S, 
    Array2 <doublevar> & Hin, 
    Array1 <doublevar> & delta_alpha,doublevar stabilization,
    Array1 <bool> & linear) { 
  int n=S.GetDim(0);
  assert(S.GetDim(1)==n);
  assert(Hin.GetDim(0)==n);
  assert(Hin.GetDim(1)==n);
  Array2 <doublevar> H=Hin;
  for(int i=1; i< n; i++) H(i,i)+=stabilization;

  Array1 <dcomplex> W(n);
  Array2 <doublevar> VL(n,n), VR(n,n);
  DGGEV(H,S,W,VL,VR);


  int min_index=0;
  doublevar min_eigenval=1e8; //W(0).real();
  doublevar max_p0=0.0;
  
  for(int j=0; j< n; j++) { 
    doublevar dpnorm=0.0;
    for(int i=0; i< n; i++) dpnorm+=VR(j,i)*VR(j,i);
    dpnorm=sqrt(dpnorm);
    for(int i=0; i< n; i++) VR(j,i)/=dpnorm;
  }

  single_write(cout,"eigenvals ");
  for(int i=0; i< n; i++) {
    single_write(cout,W(i)," ");
  }
  single_write(cout,"\n");

  for(int i=0; i< n; i++) { 
    
    if( fabs(VR(i,0)) > max_p0) {
    //if(W(i).real() < min_eigenval) { 
      min_index=i;
      max_p0=fabs(VR(i,0));
      min_eigenval=W(i).real();
    }
  }

  for(int i=0; i<n; i++) { 
    if(fabs(VR(i,0)) > 0.1)  { 
      single_write(cout,"i ",i," stab ",stabilization); 
      single_write(cout," p0 ",fabs(VR(i,0)));
      single_write(cout," eigenval ",W(i),"\n");
    }
  }
  single_write(cout,"chose ",min_index, " eigenval ", min_eigenval);
  single_write(cout,"\n");

  Array1 <doublevar> dp(n);
  for(int i=0; i < n; i++) { 
    dp(i)=VR(min_index,i);
  }
  
  doublevar dpnorm=0.0;
  for(int i=0; i< n; i++) dpnorm+=dp(i)*dp(i);
  dpnorm=sqrt(dpnorm);
  for(int i=0; i< n; i++) dp(i)/=dpnorm;

  doublevar dp0=dp(0);
  for(int i=1; i< n; i++) dp(i)/=dp(0);
  dp(0)=1.0;
  delta_alpha.Resize(n-1);
  
  for(int i=0; i< n-1; i++) delta_alpha(i)=dp(i+1);

  
  return fabs(dp0);
  
}

//--------------------------------------------------------------------
//

bool operator<( const pair<double,double> & l, const pair<double,double>& r) { return l.first < r.first; }

doublevar Linear_optimization_method::fit_stabil(Array1 <doublevar> & stabil_in, 
                                                 Array2 <doublevar> & energies_in,
                                                 int nfit) {

  vector<pair <double,double> > enstabil;
  for(int i=0; i< energies_in.GetDim(0); i++) { 
    enstabil.push_back(pair<double,double>(energies_in(i,0),stabil_in(i)));
  }
  sort(enstabil.begin(),enstabil.end());

  Array2 <doublevar> energies(nfit,2);
  Array1 <doublevar> stabil(nfit);

  for(int i=0; i< nfit; i++) { 
    energies(i,0)=enstabil[i].first;
    stabil(i)=enstabil[i].second;
    single_write(cout , "fitting to " , energies(i,0));
    single_write(cout," " , stabil(i),"\n");
  }
  
  int nstabil=energies.GetDim(0);
  Array2 <doublevar> descriptors(nstabil,3);
  for(int i=0; i< nstabil; i++) { 
    descriptors(i,0)=1.;
    descriptors(i,1)=stabil[i];
    descriptors(i,2)=stabil[i]*stabil[i];
  }

  Array1 <doublevar> dty(3,0.0);
  for(int j=0 ;j< 3; j++) { 
    for(int i=0; i< nstabil; i++) {
      dty(j)+=descriptors(i,j)*energies(i,0);
    }
  }

  //Now we need D^TD
  Array2 <doublevar> tmp(3,3,0.0),inv(3,3,0.0);
  for(int j1=0; j1 < 3; j1++) { 
    for(int j2=0; j2 < 3; j2++) { 
      for(int i=0; i< nstabil; i++) { 
        tmp(j1,j2)+=descriptors(i,j1)*descriptors(i,j2);
        
      }
    }
  }
  //Solving D^T D p = D^T y 
  Array1 <doublevar> optparms(3,0.0);
  Array1 <int> index(3);doublevar d;
  inv=tmp;
  if(ludcmp(inv,3,index,d)) { 
    optparms=dty;
    lubksb(inv,3,index,optparms);
  }
  else { optparms(2)=-50; } 

  //Find the predicted values
  Array1 <doublevar> predict(nstabil,0.0);
  doublevar err2=0,err=0;
  doublevar val2=0,val=0;
  for(int i=0; i< nstabil; i++) { 
    for(int j=0; j< 3; j++) { 
      predict(i)+=descriptors(i,j)*optparms(j);
    }
    double delta=predict(i)-energies(i,0);
    err2+=(delta*delta)/nstabil;
    err+=delta/nstabil;
    val2+=(energies(i,0)*energies(i,0))/nstabil;
    val+=energies(i,0)/nstabil;
  }
  doublevar r2= 1.0 - (err2-err*err)/(val2-val*val);
  single_write(cout,"fit r2 ",r2);
  single_write(cout," curvature ",optparms(2),"\n");


  doublevar min_stabil= - optparms(1)/(2*optparms(2));  
  // If we didn't find a minimum or if the fit is very bad, 
  // just take the lowest energy wave function
  if(optparms(2) < 0 or r2 < 0.5) { 
    int min_stabil_i=0;
    for(int i=0; i< nstabil; i++) { 
      if(energies(i,0) < energies(min_stabil_i,0)) { 
        min_stabil_i=i;
      }
    }
    min_stabil=stabil(min_stabil_i);
  }

  //Limit our stabilization to the max/min in the range.
  doublevar stabilmin=1e8,stabilmax=0;
  for(int i=0; i< nstabil; i++) { 
    if(stabilmin > stabil(i)) stabilmin=stabil(i);
    if(stabilmax < stabil(i)) stabilmax=stabil(i);
  }
  if(min_stabil < stabilmin) min_stabil=stabilmin;
  if(min_stabil > stabilmax) min_stabil=stabilmax;
  return min_stabil; 
  
  
}


//----------------------------------------------------------------

void dimension_reduction(Array2 <doublevar> & H, //input
                         Array2 <doublevar> & S, //input
                         Array2 <doublevar> & Pinv, //output
                         Array2 <doublevar> & Htilde, //reduced dimension H
                         Array2 <doublevar> & Stilde,  //reduced dimension S

                         doublevar cutoff=1e-5
                         ) {
  int n=H.GetDim(0);
  assert(n==H.GetDim(1));
  assert(n==S.GetDim(0));
  assert(n==S.GetDim(1));
  Array1 <doublevar> sigma;
  Array2 <doublevar> U;
  DGESVD(S,sigma,U,Pinv);
  for(int i=0; i< n; i++) single_write(cout, sigma(i)," ");
  single_write(cout,"\n");

  //Find the 1 eigenvalue and put it first
  int ibase=-1;
  for(int i=0; i< n; i++) { 
    if(fabs(sigma(i)-1) < 1e-8) { 
      if(ibase!=-1) { 
        single_write(cout,"WARNING: found more than one unit eigenvalue in S\n");
      }
      ibase=i;
    }
  }
  single_write(cout,"ibase ",ibase,"\n");
  
  //This is not particularly efficient, but I don't think it matters too much.
  Array2 <doublevar> permute(n,n,0.0);
  for(int i=0; i< n; i++) {
    permute(i,i)=1.0;
  }
  permute(0,0)=0.0; 
  permute(ibase,ibase)=0.0;
  permute(0,ibase)=1.0;
  permute(ibase,0)=1.0;
  Array2 <doublevar> t=Pinv;
  MultiplyMatrices(t,permute,Pinv,n);
  t=U;
  MultiplyMatrices(permute,t,U,n);
  doublevar tt=sigma(0);
  sigma(0)=sigma(ibase);
  sigma(ibase)=tt;

  int ntilde=n;
  for(int i=0; i< n; i++) {
    if(sigma(i) < cutoff) {
      ntilde=i;
      break;
    }
  }

  Array2 <doublevar> tmp(n,n),Ht(n,n);
  tmp=0.0; Ht=0.0;
  MultiplyMatrices(H,U,tmp,n);
  MultiplyMatrices(Pinv,tmp,Ht,n);

  Array2 <doublevar> Stest(n,n,0.0);
  MultiplyMatrices(S,Pinv,tmp,n);
  MultiplyMatrices(U,tmp,Stest,n);

  
  Htilde.Resize(ntilde,ntilde);
  Stilde.Resize(ntilde,ntilde);
  Stilde=0;
  for(int i=0; i< ntilde; i++) { 
    for(int j=0; j < ntilde; j++) { 
      Htilde(i,j)=Ht(i,j);
    }
    Stilde(i,i)=sigma(i);
  }

  single_write(cout,"Stilde\n");
  for(int i=0; i< ntilde; i++) {
    for(int j=0; j< ntilde; j++) { 
      single_write(cout,Stilde(i,j), " ");
    }
    single_write(cout,"\n");
  }

  single_write(cout,"Htilde\n");
  for(int i=0; i< ntilde; i++) {
    for(int j=0; j< ntilde; j++) { 
      single_write(cout,Htilde(i,j), " ");
    }
    single_write(cout,"\n");
  }

  single_write(cout,"Stest\n");
  for(int i=0; i< ntilde; i++) {
    for(int j=0; j< ntilde; j++) { 
      single_write(cout,Stest(i,j), " ");
    }
    single_write(cout,"\n");
  }
  
  
}

//----------------------------------------------------------------

void recover_deltap(Array2 <doublevar> & Pinv,
                    Array1 <doublevar> & dptilde, //this should be the output from find_directions (nparm)
                    Array1 <doublevar> & dpout) { 
  int ntilde=dptilde.GetDim(0);
  int n=Pinv.GetDim(0);
  Array1 <doublevar> tmp(n,0.0);
  tmp(0)=1.0;
  for(int i=1; i < ntilde; i++) tmp(i)=dptilde(i-1);
  Array1 <doublevar> dp(n,0.0);
  for(int i=0; i< n; i++) {
    for(int j=0; j< n; j++) { 
      dp(i)+=Pinv(j,i)*tmp(j);
    }
  }

  dpout.Resize(n-1);
  for(int i=1; i<n; i++) dpout(i-1)=dp(i)/dp(0);
  
  single_write(cout,"dp ");
  for(int i=0; i< n-1; i++) single_write(cout,dpout(i)," ");
  single_write(cout,"\n");
}

//----------------------------------------------------------------
doublevar Linear_optimization_method::line_minimization(Array2 <doublevar> & S, 
    Array2 <doublevar> & Sinv, Array2 <doublevar> & H, Array1 <doublevar> & alpha, ostream & os) { 


  int n=S.GetDim(0);

  Array2 <doublevar> Pinv,Htilde,Stilde;
  if(svd_tolerance > 0) 
    dimension_reduction(H,S,Pinv,Htilde,Stilde,svd_tolerance);
  else { 
    Pinv.Resize(n,n);
    Htilde=H;
    Stilde=S;
    Pinv=0.0;
    for(int i=0; i< n; i++) Pinv=1.0;
  }
  
  Array1 <bool> linear;
  wfdata->linearParms(linear);
  Array1 <doublevar> alpha_tmp;
  alpha_tmp=alpha;
  doublevar stabilmax=100.0*fabs(H(0,0));
  vector <doublevar> acc_stabils;

  int nstabil_test=15;
  doublevar stabilbase=fabs(H(0,0))*1e-6;
  for(int i=0; i< nstabil_test; i++) {
    doublevar stabil=stabilbase*pow(5,i);
    doublevar guesspsi=find_directions(Stilde,Htilde,alpha_tmp,stabil,linear);
    if(guesspsi > minimum_psi0)
      acc_stabils.push_back(stabil);
  }

  if(acc_stabils.size() < 5) { 
    os << "Did not find enough acceptable stabilizations. Exiting." << endl;
    return 0.0;
  }
  
  int nstabil=acc_stabils.size()+1;
  Array1 <Array1 <doublevar> > alphas(nstabil);
  Array1 <doublevar> prop_psi(nstabil,1.0);
  alphas(nstabil-1)=alpha;
  for(int i=0; i < nstabil-1; i+=1) { 
    Array1 <doublevar> dptilde(Htilde.GetDim(0));
    doublevar prop_psi0=find_directions(Stilde,Htilde,dptilde,acc_stabils[i],linear);
    if(svd_tolerance > 0)
      recover_deltap(Pinv,dptilde,alphas(i));
    else
      alphas(i)=dptilde;
    prop_psi(i)=prop_psi0;
  }
  
  for(int i=0; i< nstabil-1; i++) { 
    for(int j=0; j< alpha.GetDim(0); j++) { 
      alphas(i)(j)+=alpha(j);
    }
  }
 
  Array2 <doublevar> energies_corr2(nstabil,2);
  if(do_uncorrelated_evaluation)
    uncorrelated_evaluation(alphas,energies_corr2);
  else 
    correlated_evaluation(alphas,nstabil-1,energies_corr2);

  for(int n=0; n< nstabil; n++) {
    single_write(os,"prop_psi0 ",prop_psi(n));
    if(n<nstabil-1) single_write(os," stabilization ",acc_stabils[n]);
    single_write(os," energy ", energies_corr2(n,0));
    single_write(os," +/- ",energies_corr2(n,1));
    os << endl;
  }

  ///----------
  Array1 <doublevar> stable_fit(nstabil-1,0.0);
  Array2 <doublevar> energies_fit(nstabil-1,2,0.0);
  for(int n=0; n< nstabil-1; n++) {
    stable_fit(n)=acc_stabils[n];
    energies_fit(n,0)=energies_corr2(n,0);
    energies_fit(n,1)=energies_corr2(n,1);
  }

  doublevar min_stabil=fit_stabil(stable_fit,energies_fit);
  os << "Minimum stabilization " << min_stabil << endl;
  Array1 <doublevar> alphanew(alpha.GetDim(0));
  Array1 <doublevar> dptilde(Htilde.GetDim(0));
  doublevar prop_psi0=find_directions(Stilde,Htilde,dptilde,min_stabil,linear);
  if(svd_tolerance > 0) 
    recover_deltap(Pinv,dptilde,alphanew);
  else
    alphanew=dptilde;
  
  for(int i=0; i< alpha.GetDim(0); i++) alphanew(i)+=alpha(i);
  
  //for(int i=0; i< alpha.GetDim(0); i++) os << alphanew(i) << "  ";
  //os << endl;
  //os.flush();
  
  Array1 <Array1<doublevar> > alphasfinal(2);
  alphasfinal(0)=alpha;
  alphasfinal(1)=alphanew;
  correlated_evaluation(alphasfinal,0,energies_corr2);
  if(energies_corr2(1,0) < 0)
    alpha=alphanew;
  os << "Estimated change " << energies_corr2(1,0) << " +/- " << energies_corr2(1,1) << endl;
//alpha=alphas(min_alpha);
  return energies_corr2(1,0);
}
//----------------------------------------------------------------------
#include "Generate_sample.h"

void Linear_optimization_method::uncorrelated_evaluation(Array1 <Array1 <doublevar> > & alphas,Array2 <doublevar> & energies) {
  Sample_point * sample=NULL;
  Wavefunction * wf=NULL;
  sys->generateSample(sample);
  wfdata->generateWavefunction(wf);
  sample->attachObserver(wf);
  Array1 <Config_save_point> config_pos(nconfig_eval);
  Primary guide;
  
  int nwfs=alphas.GetDim(0);
  energies.Resize(nwfs,2);  
  Array1<doublevar> alpha_save;
  wfdata->getVarParms(alpha_save);
  Properties_point pt;
  Array1 <doublevar> kinetic(1),ecp(1,0.0),pseudo_test(pseudo->nTest());
  Array2 <doublevar> all_energies(nwfs,nconfig_eval);
  doublevar local;

  for(int w=0; w< nwfs; w++) { 
    wfdata->setVarParms(alphas(w));
    generate_sample(sample,wf,wfdata,&guide,nconfig_eval,config_pos,20,2);

    for(int config=0; config < nconfig_eval; config++) { 
      config_pos(config).restorePos(sample);
      wf->updateLap(wfdata,sample);
      local=sys->calcLoc(sample);
      sys->calcKinetic(wfdata,sample,wf,kinetic);
      pseudo->calcNonloc(wfdata,sys,sample,wf,ecp);
      all_energies(w,config)=local+kinetic(0)+ecp(0);
    }
  }
  energies=0.0;
  for(int w=0; w< nwfs; w++) { 
    for(int config=0; config < nconfig_eval; config++) { 
      energies(w,0)+=all_energies(w,config)/nconfig_eval;
    }
    for(int config=0; config < nconfig_eval; config++) { 
      doublevar d=all_energies(w,config)-energies(w,0);
      energies(w,1)+=d*d/nconfig_eval;
    }
    energies(w,1)=sqrt(energies(w,1)/nconfig_eval);
  }
  wfdata->setVarParms(alpha_save);
  wfdata->clearObserver();
  delete wf;
  delete sample;  
}

//----------------------------------------------------------------------

void Linear_optimization_method::correlated_evaluation(Array1 <Array1 <doublevar> > & alphas,int ref_alpha,Array2 <doublevar> & energies) {

  Sample_point * sample=NULL;
  Wavefunction * wf=NULL;
  sys->generateSample(sample);
  wfdata->generateWavefunction(wf);
  sample->attachObserver(wf);
  Array1 <Config_save_point> config_pos(nconfig_eval);
  Primary guide;
  int nwfs=alphas.GetDim(0);
  energies.Resize(nwfs,2);  
  
  Array1<doublevar> alpha_save;
  wfdata->getVarParms(alpha_save);
  

  wfdata->setVarParms(alphas(ref_alpha));

  generate_sample(sample,wf,wfdata,&guide,nconfig_eval,config_pos,10,5);
    

  Array2 <doublevar> all_energies(nwfs,nconfig_eval),all_weights(nwfs,nconfig_eval);
  Array2 <Wf_return> wf_vals(nwfs,nconfig_eval);
  Properties_point pt;
  Array1 <doublevar> kinetic(1,0.0),ecp(1,0.0),pseudo_test(pseudo->nTest(),0.0);
  doublevar local;
  for(int config=0; config < nconfig_eval; config++) { 
     config_pos(config).restorePos(sample);
     for(int i=0; i < pseudo_test.GetDim(0); i++) 
       pseudo_test(i)=rng.ulec();
     for(int w=0; w< nwfs; w++) {
       wfdata->setVarParms(alphas(w));
       wf->updateLap(wfdata,sample);
       wf_vals(w,config).Resize(wf->nfunc(),2);
       wf->getVal(wfdata,0,wf_vals(w,config));
       local=sys->calcLoc(sample);
       sys->calcKinetic(wfdata,sample,wf,kinetic);
       pseudo->calcNonlocWithTest(wfdata,sys,sample,wf,pseudo_test,ecp);
       all_energies(w,config)=local+kinetic(0)+ecp(0);
     }
  }

  for(int w=0; w< nwfs; w++) {
    for(int config=0; config < nconfig_eval; config++){ 
      doublevar weight=exp(2*(wf_vals(w,config).amp(0,0)-wf_vals(ref_alpha,config).amp(0,0)));
      //if(fabs(weight) > 3.0) {
      //  cout << "WARNING: cut off weight of " << weight << endl;        
      //  weight=3.0;
      //}
      all_weights(w,config)=weight;
    }
  }

  int totnconfig_eval=nconfig_eval*mpi_info.nprocs;
#ifdef USE_MPI
  Array2<doublevar> allproc_energies(nwfs,totnconfig_eval);
  Array2<doublevar> allproc_weights(nwfs,totnconfig_eval);
  for(int w=0; w< nwfs; w++) {
    MPI_Allgather(all_energies.v+w*nconfig_eval,nconfig_eval,
                  MPI_DOUBLE, allproc_energies.v+w*totnconfig_eval,
                  nconfig_eval,MPI_DOUBLE,MPI_Comm_grp);
    MPI_Allgather(all_weights.v+w*nconfig_eval,nconfig_eval,
                  MPI_DOUBLE, allproc_weights.v+w*totnconfig_eval,
                  nconfig_eval,MPI_DOUBLE,MPI_Comm_grp);
  }
#else
  Array2 <doublevar> allproc_energies=all_energies;
  Array2 <doublevar> allproc_weights=all_weights;
#endif 

  int npartitions=20;
  Array2 <doublevar> diffpartition(nwfs,npartitions);
  for(int w=0; w< nwfs; w++) { 
    for(int p=0; p < npartitions; p++) { 
      doublevar weightsum=0.0;
      doublevar ensum=0.0;
      for(int config=p; config<totnconfig_eval; config+=npartitions) { 
         ensum+=allproc_energies(w,config)*allproc_weights(w,config);
         weightsum+=allproc_weights(w,config);
      }
      diffpartition(w,p)=ensum/weightsum;
      //single_write(cout, "partitions " , w , " ");
      //single_write(cout, p , " " , diffpartition(w,p),"\n");
    }
  }

  Array1<doublevar> ref_wf(npartitions);
  for(int p=0; p < npartitions; p++) 
    ref_wf(p)=diffpartition(ref_alpha,p);
  for(int w=0; w< nwfs; w++) { 
    for(int p=0; p < npartitions; p++) { 
      diffpartition(w,p)-=ref_wf(p);
      //cout << "partitions " << w << " "<< p << " " << diffpartition(w,p) 
      //  << endl;
    }
  }

  
  energies=0.0;
  for(int w=0; w< nwfs; w++) { 
    for(int p=0; p < npartitions; p++) { 
      energies(w,0)+=diffpartition(w,p)/npartitions;
    }
    for(int p=0; p < npartitions; p++) { 
      doublevar d=diffpartition(w,p)-energies(w,0);
      energies(w,1)+=d*d/npartitions;
    }
    energies(w,1)=sqrt(energies(w,1)/npartitions);
  }
  

  wfdata->setVarParms(alpha_save);

  wfdata->clearObserver();
  delete wf;
  delete sample;
}

//----------------------------------------------------------------------



//----------------------------------------------------------------------

void Linear_optimization_method::wavefunction_derivative(
    Array2 <doublevar> & H,Array2<doublevar> & S, Array1 <doublevar> & en) { 
  int n=wfdata->nparms();

  Properties_final_average final;
  

  string vmc_section="VMC nstep ";
  append_number(vmc_section,vmc_nstep);
  vmc_section+="  nblock 20 average { WF_PARMDERIV "; 
  if(pseudopotential_derivatives) vmc_section+="EVALUATE_PSEUDOPOTENTIAL";
  vmc_section+=" nodal_cutoff ";
  vmc_section+=to_string(nodal_cutoff);
  vmc_section+="} ";
  vector <string> words;
  string sep=" ";
  split(vmc_section,sep,words);
  unsigned int pos=0;
  Vmc_method vmc;
  vmc.read(words,pos,options);
  Properties_manager prop;
  string name=options.runid+".vmcout";
  ofstream vmcout;
  if(mpi_info.node==0) vmcout.open(name.c_str());
  vmc.runWithVariables(prop,sys,wfdata,pseudo,vmcout);
  prop.getFinal(final);


  Average_return &  deriv_avg=final.avgavg(0,0);
  Average_return & deriv_err=final.avgerr(0,0);

  
  for(int i=0; i< n; i++) { 
    single_write(cout, "energy derivative ",deriv_avg.vals(2*n+i)," +/- ");
    single_write(cout, deriv_err.vals(2*n+i), "\n");
  }
  
  H.Resize(n+1,n+1);
  S.Resize(n+1,n+1);
  H=0.; S=0.;
  en.Resize(2);
  en(0)=final.avg(Properties_types::total_energy,0);
  en(1)=sqrt(final.err(Properties_types::total_energy,0));
  S(0,0)=1;
  for(int i=0; i< n; i++) { 
    for(int j=0; j< n; j++) { 
      S(i+1,j+1)=deriv_avg.vals(3*n+i*n+j)-deriv_avg.vals(n+i)*deriv_avg.vals(n+j);
    }
  }

  H(0,0)=en(0);
  for(int i=0; i < n; i++) { 
    H(i+1,0)=(deriv_avg.vals(i)-en(0)*deriv_avg.vals(n+i));
    H(0,i+1)=(H(i+1,0)+deriv_avg.vals(2*n+i));
  }
  
  for(int i=0; i< n; i++) { 
    for(int j=0; j< n; j++) { 
      H(i+1,j+1)=deriv_avg.vals(3*n+n*n+i*n+j)
        -deriv_avg.vals(n+i)*deriv_avg.vals(j)
        -deriv_avg.vals(n+j)*deriv_avg.vals(i)
        +deriv_avg.vals(n+i)*deriv_avg.vals(n+j)*en(0)
        +deriv_avg.vals(3*n+2*n*n+i*n+j)
        -deriv_avg.vals(2*n+j)*deriv_avg.vals(n+i);
     
    }
  }
}


//----------------------------------------------------------------------

void Linear_optimization_method::wavefunction_energy(
    Array1 <doublevar> & energies) {
  string vmc_section="VMC nconfig 1 nstep ";
  append_number(vmc_section,vmc_nstep);
  vmc_section+=" timestep 1.0 nblock 32  ";
  vector <string> words;
  string sep=" ";
  split(vmc_section,sep,words);
  unsigned int pos=0;
  Vmc_method vmc;
  vmc.read(words,pos,options);
  Properties_manager prop;
  string name=options.runid+"vmcout";
  ofstream vmcout;
  if(mpi_info.node==0) vmcout.open(name.c_str());
  vmc.runWithVariables(prop,sys,wfdata,pseudo,vmcout);
  Properties_final_average final;
  prop.getFinal(final);
  energies.Resize(2);
  energies(0)=final.avg(Properties_types::total_energy,0);
  energies(1)=sqrt(final.err(Properties_types::total_energy,0));
  
}
//----------------------------------------------------------------------

