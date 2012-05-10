
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

void Linear_optimization_method::read(vector <string> words,
            unsigned int & pos, Program_options & options_) { 

  options=options_;
  if(!readvalue(words,pos=0, iterations, "ITERATIONS")) 
    iterations=30;
  if(!readvalue(words,pos=0, vmc_nstep,"VMC_NSTEP"))
    vmc_nstep=100;
  if(! readvalue(words, pos=0, wfoutputfile, "WFOUTPUT") )
      wfoutputfile=options.runid+".wfout";
  if(!readvalue(words,pos=0,nconfig_eval,"FIT_NCONFIG")) 
    nconfig_eval=200;
  
  if(!readvalue(words, pos=0, tau, "TAU") )
    tau=0.1;
  
  allocate(options.systemtext[0],  sys);
  sys->generatePseudo(options.pseudotext, pseudo);
  wfdata=NULL;
  allocate(options.twftext[0], sys, wfdata);

  if(wfdata->nparms() <= 0 ) 
    error("There appear to be no parameters to optimize!");

  
}

//----------------------------------------------------------------------

int Linear_optimization_method::showinfo(ostream & os) { 
  os << "Tau  " << tau << endl;
  return 1;
}

//----------------------------------------------------------------------

void Linear_optimization_method::run(Program_options & options, ostream & output) { 
  int nparms=wfdata->nparms();
  Array1 <doublevar> x(nparms+1); 
  Array2 <doublevar> S(nparms+1,nparms+1),Sinv(nparms+1,nparms+1);
  Array2 <doublevar> H(nparms+1,nparms+1);
  Array1 <doublevar> alpha(nparms); 
 
  wfdata->getVarParms(alpha); 
  Array1 <doublevar> en,tmp_en;
  Array2 <doublevar> energy_step(iterations,2);
  Array2 <doublevar> alpha_step(iterations,nparms);
  
  wfdata->getVarParms(alpha);
  //cout << "alpha0 ";
  //for(int i=0; i<nparms; i++) 
  //  cout << alpha(i) << " ";
  //cout << endl;
  for(int it=0; it< iterations; it++) {
    //cout << "wf derivative " << endl;
    wavefunction_derivative(H,S,en);

    output << "energy " << en(0) << " +/- " << en(1) << endl;
    line_minimization(S,Sinv,H,alpha);
    //cout << "alpha ";
    //for(int i=0; i< nparms; i++) cout << alpha(i) << " ";
    //cout << endl;
    //cout << "setvarparms " << alpha.GetDim(0) <<  endl;
    wfdata->setVarParms(alpha);
    wfdata->renormalize();
    string indentation="";
    ofstream wfoutput(wfoutputfile.c_str());
    wfoutput.precision(15);
    wfdata->writeinput(indentation,wfoutput);
    wfoutput.close();
    
    for(int i=0; i< nparms; i++) {
      alpha_step(it,i)=alpha(i);
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

void Linear_optimization_method::line_minimization(Array2 <doublevar> & S, 
    Array2 <doublevar> & Sinv, Array2 <doublevar> & H, Array1 <doublevar> & alpha) { 
  int nparms=wfdata->nparms();
  Sinv.Resize(nparms+1,nparms+1);
  doublevar stabilization=0.0;
  for(int i=1; i< nparms+1; i++) H(i,i)+=stabilization;
  
  cout << "here " << endl;
  for(int i=0; i< nparms+1; i++) { 
    cout << "S ";
    for(int j=0; j< nparms+1; j++) { 
      cout <<  S(i,j) << " ";
    }
    cout << endl;
  }
  for(int i=0; i< nparms+1; i++) { 
    cout << "H ";
    for(int j=0; j< nparms+1; j++) { 
      cout <<  H(i,j) << " ";
    }
    cout << endl;
  }
  
  Sinv=0.0;
  InvertMatrix(S,Sinv,nparms+1);
  
  for(int i=0; i< nparms+1; i++) { 
    cout << "Sinv ";
    for(int j=0; j< nparms+1; j++) { 
      cout <<  Sinv(i,j) << " ";
    }
    cout << endl;
  }
  
  
  Array2 <doublevar> prodmatrix(nparms+1,nparms+1);
  prodmatrix=0.0;
  for(int i=0; i< nparms+1; i++) {
    for(int j=0; j< nparms+1; j++) { 
      for(int k=0; k< nparms+1; k++) { 
        prodmatrix(i,k)+=Sinv(i,j)*H(j,k);
      }
    }
  }

  for(int i=0; i< nparms+1; i++) { 
    cout << "P ";
    for(int j=0; j< nparms+1; j++) { 
      cout << setprecision(2) << prodmatrix(i,j) << " ";
    }
    cout << endl;
  }
  
  Array1 <dcomplex> W(nparms+1);
  Array2 <doublevar> VL(nparms+1,nparms+1), VR(nparms+1,nparms+1);
  GeneralizedEigenSystemSolverRealGeneralMatrices(prodmatrix,W,VL,VR);
  
  int min_index=0;
  int min_eigenval=W(0).real();
  
  for(int i=0; i< nparms+1; i++) { 
    cout << "eigenvalue " << i << " " << W(i) << endl;
    if(W(i).real() < min_eigenval) { 
      min_index=i;
      min_eigenval=W(i).real();
    }
  }

  for(int i=0; i< nparms+1; i++ ) {
    cout << "VR ";
    for(int j=0; j< nparms+1; j++) { 
      cout << VR(i,j) << " ";
    }
    cout << endl;
  }


  for(int i=0; i< nparms+1; i++ ) {
    cout << "VL ";
    for(int j=0; j< nparms+1; j++) { 
      cout << VL(i,j) << " ";
    }
    cout << endl;
  }
  

  Array1 <doublevar> dp(nparms+1);
  cout << "initial eigenvector ";
  for(int i=0; i < nparms+1; i++) { 
    dp(i)=VL(min_index,i);
    cout << dp(i) << " ";
  }
  cout << endl;

  for(int i=1; i< nparms+1; i++) dp(i)/=dp(0);
  dp(0)=1.0;

  cout << "Parameter variation ";
  for(int i=0; i< nparms+1; i++) 
    cout << dp(i) << " ";
  cout << endl;

  doublevar xi=0.5;
  Array1 <doublevar> norm(nparms+1);
  /*
  norm=0.;
  doublevar denominator_sum=0.0;
  for(int j=1; j < nparms+1; j++) { 
    for(int k=1; k< nparms+1; k++) { 
      denominator_sum+=dp(j)*dp(k)*S(j,k);
    }
  }
  
  for(int i=1; i< nparms+1; i++) { 
    double numerator_sum=0.0;
    double denominator_sum=0.0;
    for(int j=1; j < nparms+1; j++) { 
      numerator_sum+=dp(j)*S(i,j);
    }
    norm(i)=-(1-xi)*numerator_sum/(1-xi+xi*sqrt(1+denominator_sum));
  }
  */
  doublevar D=1.0;
  for(int j=1; j< nparms+1; j++) { 
    D+=2*S(0,j)*dp(j);
    for(int k=0; k< nparms+1; k++) {
      D+=S(j,k)*dp(j)*dp(k);
    }
  }
  for(int i=1; i< nparms+1;  i++) { 
    doublevar num_sum=0.0;
    doublevar denom_sum=0.0;
    for(int j=1; j < nparms+1; j++) { 
      num_sum+=S(i,j)*dp(j);
      denom_sum+=S(0,j)*dp(j);
    }
    norm(i)=-(xi*D*S(0,i)+(1-xi)*(S(0,i)+num_sum))
      /(xi*D+(1-xi)*(1+denom_sum));
  }




  doublevar renorm_dp=0.;
  for(int i=1; i< nparms+1; i++) { 
    renorm_dp+=norm(i)*dp(i);
  }
  cout << "renorm_dp " << renorm_dp << endl;

  for(int i=0; i< nparms; i++) { 
    alpha(i)+=dp(i+1)/(1-renorm_dp);
    cout << "new alpha " << alpha(i) << endl;
  }


/*
  Array2 <doublevar> energies_corr2(ntau,2);
  correlated_evaluation(alphas,0,energies_corr2);
  doublevar mixing=0.1;
  doublevar min_en=energies_corr(0,0)+energies_corr(0,1)*energies_corr(0,1)*mixing;
  doublevar min_n=0;
  for(int n=0; n< ntau; n++) { 
    single_write(cout,tau_prefactor[n]*tau," " ,energies_corr2(n,0)," ");
    single_write(cout,energies_corr2(n,1),"\n");
    doublevar opt_val=energies_corr2(n,0);
    if(opt_val < min_en && energies_corr2(n,1) < 1.5*energies_corr2(0,1) ) { 
      min_en=opt_val;
      min_n=n;
    }
  }
  alpha=alphas(min_n);
  cout << "alpha ";
  for(int i=0; i< nparms;i++) {
    cout << alpha(i) << " ";
  }
  cout << endl;
*/
  
}
//----------------------------------------------------------------------

#include "Guiding_function.h"
#include "Split_sample.h"
void Linear_optimization_method::correlated_evaluation(Array1 <Array1 <doublevar> > & alphas,int ref_alpha,Array2 <doublevar> & energies) {
  Sample_point * sample=NULL;
  Wavefunction * wf=NULL;
  sys->generateSample(sample);
  wfdata->generateWavefunction(wf);
  sample->attachObserver(wf);
  int nstep=10;
  doublevar timestep=0.3;
  Array1 <Config_save_point> config_pos(nconfig_eval);
  Split_sampler sampler;
  sampler.setRecursionDepth(2);
  Dynamics_info dinfo;
  Primary guide;
  sample->randomGuess();
  int nelectrons=sample->electronSize();
  wfdata->setVarParms(alphas(ref_alpha));
  for(int config=0; config < nconfig_eval; config++) { 
    for(int step=0; step < nstep; step++) { 
      for(int e=0; e < nelectrons; e++) { 
        sampler.sample(e,sample,wf,wfdata,&guide,dinfo,timestep);
      }
    }
    config_pos(config).savePos(sample);
  }

  Properties_gather mygather;
  int nwfs=alphas.GetDim(0);
  Array2 <doublevar> all_energies(nwfs,nconfig_eval);
  Array2 <Wf_return> wf_vals(nwfs,nconfig_eval);
  Properties_point pt;
  for(int config=0; config < nconfig_eval; config++) { 
     config_pos(config).restorePos(sample);
     for(int w=0; w< nwfs; w++) {
       wfdata->setVarParms(alphas(w));
       wf->updateVal(wfdata,sample);
       wf_vals(w,config).Resize(wf->nfunc(),2);
       wf->getVal(wfdata,0,wf_vals(w,config));
       mygather.gatherData(pt,pseudo,sys,wfdata,wf,sample,&guide);
       all_energies(w,config)=pt.energy(0);
     }
  }
  Array1 <doublevar> avg_energies(nwfs),avg_weight(nwfs),avg_var(nwfs);
  avg_energies=0.0;
  avg_weight=0.0;
  avg_var=0.0;
  for(int w=0; w< nwfs; w++) {
    doublevar avg_en_unweight=0;
    for(int config=0; config < nconfig_eval; config++)  { 
      doublevar weight=exp(2*(wf_vals(w,config).amp(0,0)-wf_vals(ref_alpha,config).amp(0,0)));
      avg_energies(w)+=weight*all_energies(w,config)/nconfig_eval;
      avg_en_unweight+=all_energies(w,config)/nconfig_eval;
      avg_weight(w)+=weight/nconfig_eval;
    }
    for(int config=0; config < nconfig_eval; config++) { 
      avg_var(w)+=(all_energies(w,config)-avg_en_unweight)*(all_energies(w,config)-avg_en_unweight);
    }
    avg_var(w)=sqrt(avg_var(w))/nconfig_eval;
  }
 
  energies.Resize(nwfs,2);
  for(int w=0; w< nwfs; w++) { 
    energies(w,0)=avg_energies(w)/avg_weight(w);//+0.1*avg_var(w);
    energies(w,1)=avg_var(w);
    //cout << w << " " << avg_energies(w)/avg_weight(w) <<
    //   "  " << avg_energies(w) << "  " << avg_weight(w) << endl;
    
  }
  //cout << "done " << endl;
  wfdata->clearObserver();
  delete wf;
  delete sample;
  
}

//----------------------------------------------------------------------



void Linear_optimization_method::wavefunction_derivative(
    Array2 <doublevar> & H,Array2<doublevar> & S, Array1 <doublevar> & en) { 
  string vmc_section="VMC nconfig 1 nstep ";
  append_number(vmc_section,vmc_nstep);
  vmc_section+=" timestep 1.0 nblock 20 average { WF_PARMDERIV } ";
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
  Average_return &  deriv_avg=final.avgavg(0,0);
  Average_return & deriv_err=final.avgerr(0,0);

  bool nonzero_element=false;
  for(int i=0; i< deriv_avg.vals.GetDim(0); i++) { 
    //cout << "avg deriv " << deriv_avg.vals(i) << " " << deriv_err.vals(i) << endl;
    //holding the significance to some number of sigmas.
    //if( fabs(deriv_avg.vals(i))/deriv_err.vals(i) < 3)  
    //  deriv_avg.vals(i)=0.0;
    //else 
    //  nonzero_element=true;
    nonzero_element=true;
  }
  if(!nonzero_element) { 
    cout << "WARNING: set all elements to zero because they are not significant."
      << " Increasing vmc_nstep to " << vmc_nstep*4<< endl;
    vmc_nstep*=4;
    wavefunction_derivative(H, S,en);
    return;
  }

  int n=wfdata->nparms();
  for(int i=0; i< n; i++) { 
    cout << "energy derivative " << deriv_avg.vals(2*n+i) << endl;
  }
  H.Resize(n+1,n+1);
  S.Resize(n+1,n+1);
  H=0.; S=0.;
  en.Resize(2);
  en(0)=final.avg(Properties_types::total_energy,0);
  en(1)=sqrt(final.err(Properties_types::total_energy,0));
  S(0,0)=1;
  //S(0,0)=deriv_avg.vals(3*n+3*n*n);
  doublevar s_renorm=sqrt(deriv_avg.vals(3*n+3*n*n));
  for(int i=0; i < n; i++) { 
    //S(0,i+1)=0.0;
    //S(i+1,0)=0.0;
    S(0,i+1)=deriv_avg.vals(n+i);
    S(i+1,0)=deriv_avg.vals(n+i);
  }
  for(int i=0; i< n; i++) { 
    for(int j=0; j< n; j++) { 
      //S(i+1,j+1)=(deriv_avg.vals(3*n+i*n+j)-deriv_avg.vals(n+i)*deriv_avg.vals(n+j))/s_renorm;
      S(i+1,j+1)=deriv_avg.vals(3*n+i*n+j);
    }
  }
  H(0,0)=en(0);
  for(int i=0; i < n; i++) { 
    //H(i+1,0)=deriv_avg.vals(i)-en(0)*deriv_avg.vals(n+i);
    //H(0,i+1)=H(i+1,0)+deriv_avg.vals(2*n+i);
    H(i+1,0)=deriv_avg.vals(i);
    H(0,i+1)=deriv_avg.vals(2*n+i);
  }
  
  for(int i=0; i< n; i++) { 
    for(int j=0; j< n; j++) { 
      H(i+1,j+1)=deriv_avg.vals(3*n+2*n*n+i*n+j);
     // H(i+1,j+1)=deriv_avg.vals(3*n+n*n+i*n+j)
     //   -deriv_avg.vals(n+i)*deriv_avg.vals(j)
     //   -deriv_avg.vals(n+j)*deriv_avg.vals(i)
     //   +deriv_avg.vals(n+i)*deriv_avg.vals(n+j)*en(0)
     //   +deriv_avg.vals(3*n+2*n*n+i*n+j)
     //   -deriv_avg.vals(2*n+j)*deriv_avg.vals(n+i);
      
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

