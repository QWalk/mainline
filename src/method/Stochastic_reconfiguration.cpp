
#include "Wavefunction_data.h"
#include "Sample_point.h"
#include "Guiding_function.h"
#include "System.h"
#include "Pseudopotential.h"
#include "Stochastic_reconfiguration.h"
#include "qmc_io.h"
#include "Program_options.h"
#include "MatrixAlgebra.h"
#include "Vmc_method.h"
#include "Properties_average.h"
#include "Properties_block.h"

void Stochastic_reconfiguration_method::read(vector <string> words,
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
    tau=0.5;
  
  allocate(options.systemtext[0],  sys);
  sys->generatePseudo(options.pseudotext, pseudo);
  wfdata=NULL;
  allocate(options.twftext[0], sys, wfdata);

  if(wfdata->nparms() <= 0 ) 
    error("There appear to be no parameters to optimize!");

  
}

//----------------------------------------------------------------------

int Stochastic_reconfiguration_method::showinfo(ostream & os) { 
  os << "Tau  " << tau << endl;
  return 1;
}

//----------------------------------------------------------------------

void Stochastic_reconfiguration_method::run(Program_options & options, ostream & output) { 
  int nparms=wfdata->nparms();
  Array1 <doublevar> x(nparms+1); 
  Array2 <doublevar> S(nparms+1,nparms+1),Sinv(nparms+1,nparms+1),tmp_S(nparms+1,nparms+1);
  Array1 <doublevar> alpha(nparms); 
  Array1 <doublevar> energies(nparms+1),tmp_energies(nparms+1);
  wfdata->getVarParms(alpha); 
  Array1 <doublevar> en,tmp_en;
  Array2 <doublevar> energy_step(iterations,2);
  Array2 <doublevar> alpha_step(iterations,nparms);
  
  wfdata->getVarParms(alpha);
  //cout << "alpha0 ";
  //for(int i=0; i<nparms; i++) 
  //  cout << alpha(i) << " ";
  //cout << endl;
  cout.precision(1);
  int nit_completed=0;
  for(int it=0; it< iterations; it++) {
    cout << "_________start iteration_______" << endl;
    //cout << "wf derivative " << endl;
    wavefunction_derivative(energies,S,en);

    output << "energy " << en(0) << " +/- " << en(1) << endl;
    line_minimization(S,Sinv,energies,alpha);
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

  //cout << "alphas " << endl;
  //
  /*
  Array1 <Array1 <doublevar> > alphas(iterations);
  for(int it=0; it< iterations; it++) { 
    alphas(it).Resize(nparms);
    for(int p=0; p< nparms; p++) {
      alphas(it)(p)=alpha_step(it,p);
    }
  }
  Array2 <doublevar> energies_iteration(iterations,2);
  correlated_evaluation(alphas,iterations/2,energies_iteration);
  */

#ifdef USE_MPI
  MPI_Barrier(MPI_Comm_grp);
#endif
  

  ifstream wfinput(wfoutputfile.c_str());
  options.twftext[0].clear();
  parsefile(wfinput, options.twftext[0]);
  wfinput.close();

}

//----------------------------------------------------------------------

void Stochastic_reconfiguration_method::line_minimization(Array2 <doublevar> & S, 
    Array2 <doublevar> & Sinv, Array1 <doublevar> & energies, Array1 <doublevar> & alpha) { 
  int nparms=S.GetDim(0)-1;
  Array1 <doublevar> save_alpha=alpha;

  int ntau=5;

  InvertMatrix(S,Sinv,nparms+1);
  double  tau_prefactor[5]={0.0,0.25,0.5,0.75,1.0};
  Array1 <doublevar> x(nparms+1);
  Array1 <Array1 <doublevar> > alphas(ntau);

  for(int n=0; n < ntau; n++) { 
    doublevar tmp_tau=tau_prefactor[n]*tau;
    x=0.0;
    for(int i=0; i <nparms+1; i++) { 
      for(int j=0; j< nparms+1; j++) { 
        x(i)+=Sinv(i,j)*(S(0,j)-tmp_tau*energies(j));
      }
    }
    //cout << "x ";
    //for(int i=0; i< nparms+1; i++) cout << x(i) << " ";
    //cout << endl;
    for(int i=0; i< nparms; i++) { 
      alpha(i)=save_alpha(i)+x(i+1)/x(0);
    }
    //wfdata->setVarParms(alpha);
    //wavefunction_energy(en);
   // cout << tmp_tau << "  " << en(0) << "  " << en(1) << endl;
    alphas(n)=alpha;
  }
 
  Array2 <doublevar> energies_corr(ntau,2);
  correlated_evaluation(alphas,3,energies_corr);

  //Now do a least-squares fit to a quadratic model
  Array2 <doublevar> tmatrix(3,3),tmatrixinv(3,3);
  Array1 <doublevar> coeff(3),data(3);
  tmatrix=0.0; data=0.0;
  tmatrix(0,0)=ntau;
  for(int n=0; n< ntau; n++) { 
    doublevar t=tau_prefactor[n]*tau;
    doublevar t2=t*t, t3=t*t*t,t4=t*t*t*t;
    tmatrix(0,1)+=t;
    tmatrix(1,0)+=t;
    tmatrix(1,1)+=t2;
    tmatrix(0,2)+=t2;
    tmatrix(2,0)+=t2;
    tmatrix(1,2)+=t3;
    tmatrix(2,1)+=t3;
    tmatrix(2,2)+=t4;
    doublevar d=energies_corr(n,0);
    data(0)+=d;
    data(1)+=d*t;
    data(2)+=d*t*t;
  }
  InvertMatrix(tmatrix,tmatrixinv,3);
  coeff=0;
  for(int i=0; i< 3; i++) {
    for(int j=0; j< 3; j++) { 
      coeff(i)+=tmatrixinv(i,j)*data(j);
    }
  }
  for(int n=0; n< ntau; n++) { 
    doublevar t=tau_prefactor[n]*tau;
    doublevar f=coeff(0)+coeff(1)*t+coeff(2)*t*t;
    cout << t << " " << f << " " << energies_corr(n,0) << endl;
  }
  doublevar min_tau=-coeff(1)/(2*coeff(2));
  if(coeff(2) < 0) min_tau=0;
  min_tau=min(min_tau,tau_prefactor[ntau-1]*tau);
  min_tau=max(min_tau,0.0);
  cout << "min_tau " << min_tau << endl;
 
  for(int i=0; i <nparms+1; i++) { 
    for(int j=0; j< nparms+1; j++) { 
      x(i)+=Sinv(i,j)*(S(0,j)-min_tau*energies(j));
    }
  }
  cout << "save_alpha ";
  for(int i=0; i< nparms+1; i++) cout << save_alpha(i) << " ";
  cout << endl;
  
  cout << "x ";
  for(int i=0; i< nparms+1; i++) cout << x(i) << " ";
  cout << endl;
  if(min_tau > 1e-4) { 
    for(int i=0; i< nparms; i++) { 
      alpha(i)=save_alpha(i)+x(i+1)/x(0);
    }
  }
  else alpha=save_alpha;
  
  //cout << "done line minimization" << endl;
  
}
//----------------------------------------------------------------------

#include "Guiding_function.h"
#include "Split_sample.h"
void Stochastic_reconfiguration_method::correlated_evaluation(Array1 <Array1 <doublevar> > & alphas,int ref_alpha,Array2 <doublevar> & energies) {
  Sample_point * sample=NULL;
  Wavefunction * wf=NULL;
  sys->generateSample(sample);
  wfdata->generateWavefunction(wf);
  sample->attachObserver(wf);
  int nstep=5;
  doublevar timestep=0.7;
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
    energies(w,0)=avg_energies(w)/avg_weight(w)+0.1*avg_var(w);
    //cout << w << " " << avg_energies(w)/avg_weight(w) <<
    //   "  " << avg_energies(w) << "  " << avg_weight(w) << endl;
    
  }
  //cout << "done " << endl;
  wfdata->clearObserver();
  delete wf;
  delete sample;
  
}

//----------------------------------------------------------------------

void Stochastic_reconfiguration_method::output_average_wf(Array2 <doublevar> & alpha_step, 
    Array2 <doublevar> & energy_step, int nit_completed) { 
  int nparms=alpha_step.GetDim(0);
  assert(alpha_step.GetDim(1) >= nit_completed);
  int navg=0;
  Array1 <doublevar> avg_alpha(nparms);
  avg_alpha=0;
  Array1 <doublevar> ref_en(2);
  ref_en(0)=energy_step(nit_completed-1,0);
  ref_en(1)=energy_step(nit_completed-1,1);

  for(int  it=0; it < nit_completed; it++) { 
    if(energy_step(it,0)-ref_en(0) < 2*ref_en(1)) { 
      navg++;
      for(int i=0; i< nparms; i++) 
        avg_alpha(i)+=alpha_step(it,i);
    }
  }

  for(int i=0; i< nparms; i++) avg_alpha(i)/=navg;
  wfdata->setVarParms(avg_alpha);

  wfdata->renormalize();
  string indentation="";
  ofstream wfoutput(wfoutputfile.c_str());
  wfoutput.precision(15);
  wfdata->writeinput(indentation,wfoutput);
  wfoutput.close();

}

//----------------------------------------------------------------------


void Stochastic_reconfiguration_method::wavefunction_derivative(
    Array1 <doublevar> & energies,Array2<doublevar> & S, Array1 <doublevar> & en) { 
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
  int n=wfdata->nparms();
  energies.Resize(n+1);
  for(int i=0; i< n; i++) {
    energies(i+1)=deriv_avg.vals(i);
  }
  energies(0)=final.avg(Properties_types::total_energy,0);
  en.Resize(2);
  en(0)=energies(0);
  en(1)=sqrt(final.err(Properties_types::total_energy,0));
  cout << "energy " << en(0) <<  " +/- " << en(1) << endl;
  S.Resize(n+1,n+1);
  S(0,0)=1;
  for(int i=0; i< n; i++) { 
    S(0,i+1)=deriv_avg.vals(n+i);
    S(i+1,0)=deriv_avg.vals(n+i);
  }
  for(int i=0; i< n; i++) { 
    for(int j=0; j< n; j++) { 
      S(i+1,j+1)=deriv_avg.vals(2*n+i*n+j);
    }
  }
  
}


//----------------------------------------------------------------------

void Stochastic_reconfiguration_method::wavefunction_energy(
    Array1 <doublevar> & energies) {
  string vmc_section="VMC nconfig 1 nstep ";
  append_number(vmc_section,vmc_nstep);
  vmc_section+=" timestep 1.0 nblock 20  ";
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

