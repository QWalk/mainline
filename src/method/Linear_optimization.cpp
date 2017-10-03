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
    minimum_psi0=0.95;
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
    output << "step " << it << ": current energy " << setprecision(10) << en(0) 
           << " +/- " << setprecision(10) << en(1);
    output.flush();
    //if(it>0)
    //  output << "  energy change " << en(0)-olden(0) 
    //    << " +/- " << sqrt(en(1)*en(1)+olden(1)*olden(1)) << endl;
    
    doublevar endiff= line_minimization(S,Sinv,H,alpha);

    output << " energy change  " << endiff << endl;
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
doublevar find_directions(Array2 <doublevar> & S, Array2 <doublevar> & Sinv, 
    Array2 <doublevar> & Hin, 
    Array1 <doublevar> & delta_alpha,doublevar stabilization,
    Array1 <bool> & linear) { 
  int n=S.GetDim(0);
  assert(S.GetDim(1)==n);
  assert(Hin.GetDim(0)==n);
  assert(Hin.GetDim(1)==n);
  Array2 <doublevar> H=Hin;
  //Sinv.Resize(n,n);
  for(int i=1; i< n; i++) H(i,i)+=stabilization;
  //InvertMatrix(S,Sinv,n);

  //Array2 <doublevar> prodmatrix(n,n);

  //MultiplyMatrices(Sinv,H,prodmatrix,n);
  Array1 <dcomplex> W(n);
  Array2 <doublevar> VL(n,n), VR(n,n);
  //GeneralizedEigenSystemSolverRealGeneralMatrices(prodmatrix,W,VL,VR);
  DGGEV(H,S,W,VL,VR);
  
  int min_index=0;
  doublevar min_eigenval=W(0).real();
  doublevar max_p0=0.0;
  for(int i=0; i< n; i++) { 
    single_write(cout,"eigenvalue ",i," ");
    single_write(cout,W(i),"\n");
    
    //if(W(i).real() < min_eigenval) { 
    //  min_index=i;
    //  min_eigenval=W(i).real();
    //}
    // Here we select the eigenvector with the highest percentage of the
    // original wave function. It's possible this makes the algorithm go a little
    // slower than necessary, but it should improve the stability, especially with 
    // near-degenerate eigenvectors.
    if(fabs(VR(i,0)) > max_p0) { 
      min_index=i;
      max_p0=fabs(VR(i,0));
      min_eigenval=W(i).real();
    }
  }

  Array1 <doublevar> dp(n);
  single_write(cout,"eigenvector ");
  for(int i=0; i < n; i++) { 
    dp(i)=VR(min_index,i);
    single_write(cout,dp(i)," ");
  }
  single_write(cout,"\n");
  single_write(cout,"eigenvalue ",min_eigenval,"\n");
  
  doublevar dpnorm=0.0;
  for(int i=0; i< n; i++) dpnorm+=dp(i)*dp(i);
  dpnorm=sqrt(dpnorm);
  for(int i=0; i< n; i++) dp(i)/=dpnorm;

  doublevar dp0=dp(0);
  for(int i=1; i< n; i++) dp(i)/=dp(0);
  dp(0)=1.0;

  doublevar xi=0.5;
  Array1 <doublevar> norm(n);
  doublevar D=1.0;
  for(int j=1; j< n; j++) { 
    D+=2*S(0,j)*dp(j);
    for(int k=0; k< n; k++) {
      D+=S(j,k)*dp(j)*dp(k);
    }
  }
  D=sqrt(D);
  norm=0.0;
  
  for(int i=1; i< n;  i++) { 
    doublevar num_sum=0.0;
    doublevar denom_sum=0.0;
    for(int j=1; j < n; j++) { 
      if(!linear(j-1)) {
        num_sum+=S(i,j)*dp(j);
        denom_sum+=S(0,j)*dp(j);
      }
    }
    if(!linear(i-1)) 
      norm(i)=-(xi*D*S(0,i)+(1-xi)*(S(0,i)+num_sum))
        /(xi*D+(1-xi)*(1+denom_sum));
    else norm(i)=0.0;
  }
  



  doublevar renorm_dp=0.;
  for(int i=1; i< n; i++) { 
    renorm_dp+=norm(i)*dp(i);
  }


  delta_alpha.Resize(n-1);
  for(int i=0; i< n-1; i++) { 
    if(!linear(i))
      delta_alpha(i)=dp(i+1)/(1-renorm_dp);
    else delta_alpha(i)=dp(i+1);
  }

  
  return fabs(dp0);
  
}

//--------------------------------------------------------------------

doublevar Linear_optimization_method::line_minimization(Array2 <doublevar> & S, 
    Array2 <doublevar> & Sinv, Array2 <doublevar> & H, Array1 <doublevar> & alpha) { 
  Array1 <bool> linear;
  wfdata->linearParms(linear);
  Array1 <doublevar> alpha_tmp;
  alpha_tmp=alpha;
  doublevar stabilmax=10.0*fabs(H(0,0));
  vector <doublevar> acc_stabils;
  //cout << "checking the following alphas"<< endl;

  
  doublevar psi_0_min=find_directions(S,Sinv,H,alpha_tmp,0.0,linear);
  doublevar psi_0_max=find_directions(S,Sinv,H,alpha_tmp,stabilmax,linear);
  if(psi_0_min > psi_0_max) { 
    if(psi_0_min > minimum_psi0) { 
      single_write(cout,"Found psi_0_min > 0.95, so going forward with just that.");
      acc_stabils.push_back(0.0);
    }
    else 
      error("In LINEAR, encountered psi_0_min > psi_0_max. Perhaps you should increase TOTAL_NSTEP?",psi_0_min,psi_0_max);
  }
  psi_0_min=max(psi_0_min,minimum_psi0);

  doublevar step=(psi_0_max-psi_0_min)/5;
  doublevar tol=step/10;
  //cout << "max " << psi_0_max << " min " << psi_0_min << endl;
  for(doublevar psi0=psi_0_min+step; psi0 < psi_0_max; psi0+=step) { 
    doublevar bracket_above=stabilmax;
    doublevar bracket_below=0.0;
    doublevar guesstab,guesspsi;
    int count=0;
    do { 
      guesstab=(bracket_above-bracket_below)/2.0+bracket_below;
      guesspsi=find_directions(S,Sinv,H,alpha_tmp,guesstab,linear);
      if(guesspsi > psi0) bracket_above=guesstab;
      else bracket_below=guesstab;
      //cout << "psi0 " << psi0 << " guesstab " << guesstab << " guesspsi " << guesspsi 
      //     << " brackets " << bracket_below << " " << bracket_above << endl;
      if(++count > 100) { 
        single_write(cout,"In Linear_optimization::line_minimization: took more than 100 iterations to find sample. Giving up.");
        break;
      }
    } while(fabs(guesspsi-psi0) > tol);
    acc_stabils.push_back(guesstab);
  }


  int nstabil=acc_stabils.size()+1;
  Array1 <Array1 <doublevar> > alphas(nstabil);
  Array1 <doublevar> prop_psi(nstabil,1.0);
  alphas(0)=alpha;
  for(int i=1; i < nstabil; i+=1) { 
    doublevar prop_psi0=find_directions(S,Sinv,H,alphas(i),acc_stabils[i-1],linear);  
    cout << "i " << i << " " << acc_stabils[i-1] << " prop_psi0 " << prop_psi0 << endl;
    prop_psi(i)=prop_psi0;
    
  }

  for(int i=1; i< nstabil; i++) { 
    for(int j=0; j< alpha.GetDim(0); j++) { 
      alphas(i)(j)+=alpha(j);
    }
  }
  
  Array2 <doublevar> energies_corr2(nstabil,2);
  bool significant_stabil=false;
  while(!significant_stabil) { 
    if(do_uncorrelated_evaluation)
      uncorrelated_evaluation(alphas,energies_corr2);
    else 
      correlated_evaluation(alphas,0,energies_corr2);
    
    for(int n=0; n< nstabil; n++) { 
      single_write(cout,"alpha ");
      for(int i=0; i< alpha.GetDim(0); i++) { 
        single_write(cout,alphas(n)(i), " ");
      }
      single_write(cout,"\n");
    }
    for(int n=0; n< nstabil; n++) {
      single_write(cout,"prop_psi0 ",prop_psi(n));
      single_write(cout," energy ", energies_corr2(n,0));
      single_write(cout," +/- ",energies_corr2(n,1),"\n");
    }

    for(int n=1; n< nstabil; n++) {
      doublevar diff=energies_corr2(n,0)-energies_corr2(0,0);
      if( fabs(diff)/energies_corr2(n,1) > 3.0) {
        significant_stabil=true;
        single_write(cout, "Significant change in energy ",diff);
         single_write(cout," +/- ", energies_corr2(n,1),"\n");
      }
    }

    if(nconfig_eval >= max_nconfig_eval and !significant_stabil){
      single_write(cout, "Hit energy convergence\n");
      alpha=alphas(0);
      return 0;
    }


    if(!significant_stabil) {
      nconfig_eval*=4;
      single_write(cout,"fit not significant, increasing resolution: ",nconfig_eval,"\n");
    }
    
  }
  doublevar min_en=energies_corr2(0,0);
  int min_alpha=0;
  for(int n=0; n< nstabil; n++) { 
//single_write(cout,stabilization[n]," ",energies_corr2(n,0),"\n");
    if(energies_corr2(n,0) < min_en
        && fabs(energies_corr2(n,0))/energies_corr2(n,1) > 2.0 ) { 
      min_en=energies_corr2(n,0);
      min_alpha=n;
    }
  }
  alpha=alphas(min_alpha);
  return energies_corr2(min_alpha,0)-energies_corr2(0,0);
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
  


  generate_sample(sample,wf,wfdata,&guide,nconfig_eval,config_pos,10,5);
    

  Array2 <doublevar> all_energies(nwfs,nconfig_eval),all_weights(nwfs,nconfig_eval);
  Array2 <Wf_return> wf_vals(nwfs,nconfig_eval);
  Properties_point pt;
  Array1 <doublevar> kinetic(1),ecp(1,0.0),pseudo_test(pseudo->nTest());
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

  Array1 <doublevar> avg_energies(nwfs,0.0),avg_weight(nwfs,0.0);
  Array1 <doublevar> diff_var(nwfs,0.0),weight_var(nwfs,0.0);
  Array1 <doublevar> avg_en_unweight(nwfs,0.0);
  Array2 <doublevar> diff_en(nwfs,nconfig_eval,0.0);
  for(int w=0; w< nwfs; w++) {
    for(int config=0; config < nconfig_eval; config++){ 
      doublevar weight=exp(2*(wf_vals(w,config).amp(0,0)-wf_vals(ref_alpha,config).amp(0,0)));
      
      avg_energies(w)+=weight*all_energies(w,config)/nconfig_eval;
      avg_en_unweight(w)+=all_energies(w,config)/nconfig_eval;
      avg_weight(w)+=weight/nconfig_eval;
      diff_en(w,config)=all_energies(w,config)-all_energies(0,config);
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
      //cout << "partitions " << w << " "<< p << " " << diffpartition(w,p) 
      //  << endl;      
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
  


  /*
  
  for(int w=0; w< nwfs; w++) { 
    avg_energies(w)=parallel_sum(avg_energies(w))/mpi_info.nprocs;
    avg_weight(w)=parallel_sum(avg_weight(w))/mpi_info.nprocs;
  }

  for(int w=0; w< nwfs; w++) { 
    diff_var(w)=0.0;
    doublevar diff=avg_energies(w)/avg_weight(w)-avg_energies(0)/avg_weight(0);
//    cout << w << " diff " << diff << endl;
    for(int config=0; config < nconfig_eval; config++) { 
      diff_var(w)+=(diff_en(w,config)-diff)
                  *(diff_en(w,config)-diff);
      weight_var(w)+=(all_weights(w,config)-avg_weight(w))
                     *(all_weights(w,config)-avg_weight(w));
    }
    int npts=nconfig_eval*mpi_info.nprocs;    
    diff_var(w)=sqrt(parallel_sum(diff_var(w))/(npts*(npts)));
    weight_var(w)=sqrt(parallel_sum(weight_var(w))/(npts*(npts)));
  }

  for(int w=0; w< nwfs; w++) { 
    energies(w,0)=avg_energies(w)/avg_weight(w);//+0.1*avg_var(w);
    doublevar diff=energies(w,0)-energies(0,0);
    if(w!=0)
      energies(w,1)=sqrt(diff_var(w)*diff_var(w)*avg_weight(w)*avg_weight(w)+diff*diff*weight_var(w)*weight_var(w));
    if(mpi_info.node==0) { 
      cout << "endiff " << diff
        <<  " estimated error " << energies(w,1) 
        << " average weight " << avg_weight(w) 
        << " weight variance " << weight_var(w) << endl;
    }
  }
  */
  //cout << "done " << endl;
  wfdata->clearObserver();
  delete wf;
  delete sample;
}

//----------------------------------------------------------------------


bool  Linear_optimization_method::deriv_is_significant(Average_return & avg,
    Average_return & err,
    int n) { 
  doublevar thresh=3; //Number of sigmas above which we consider this significant
  
  int nsig_S0=0;
  for(int i=0; i< n; i++) {
    if(fabs(avg.vals(n+i)/err.vals(n+i)) > thresh) nsig_S0++;
  }

  int nsig_S=0;
  for(int i=0; i< n; i++) { 
    for(int j=0; j< n; j++) { 
      if(fabs(avg.vals(3*n+i*n+j)/err.vals(3*n+i*n+j)) > thresh) nsig_S++;
    }
  }

  int nsig_H0=0;
  int nsig_enderiv=0;
  for(int i=0; i< n; i++) {
    if(fabs(avg.vals(i)/err.vals(i)) > thresh) nsig_H0++;
    if(fabs(avg.vals(2*n+i)/err.vals(2*n+i)) > thresh) nsig_enderiv++;
  }
  int nsig_H=0;
  for(int i=0; i< n; i++) { 
    for(int j=0; j< n; j++) { 
      if(fabs(avg.vals(3*n+2*n*n+i*n+j)/err.vals(3*n+2*n*n+i*n+j)) > thresh) 
        nsig_H++;
    }
  }
  single_write(cout, "proportions significant\n");
  single_write(cout, "S0 ", double(nsig_S0)/double(n), "\n");
  single_write(cout, "S ", double(nsig_S0)/double(n*n), "\n");
  single_write(cout, "H0 ", double(nsig_H0)/double(n), "\n");
  single_write(cout, "enderiv ", double(nsig_enderiv)/double(n), "\n");
  single_write(cout, "H ", double(nsig_H0)/double(n*n), "\n");
  //if(double(nsig_H0)/double(n) > sig_H_threshold) return true;
  //else return false;
  return true;

}

//----------------------------------------------------------------------

void Linear_optimization_method::wavefunction_derivative(
    Array2 <doublevar> & H,Array2<doublevar> & S, Array1 <doublevar> & en) { 
  int n=wfdata->nparms();
  Properties_final_average final;
  
  bool significant=false;
  while(!significant) { 

    string vmc_section="VMC nstep ";
    append_number(vmc_section,vmc_nstep);
    vmc_section+="  nblock 20 average { WF_PARMDERIV "; 
    if(pseudopotential_derivatives) vmc_section+="EVALUATE_PSEUDOPOTENTIAL";
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

    significant=deriv_is_significant(final.avgavg(0,0), final.avgerr(0,0),n);
    if(!significant && vmc_nstep <=max_vmc_nstep) {
      single_write(cout, " didn't find significant derivatives: increasing vmc_nstep.\n");
      vmc_nstep*=4;
    }
    else if(!significant) { 
      single_write(cout, "WARNING: did not find significant derivatives and vmc_nstep > max_vmc_nstep.  Continuing.\n");
    }
  }

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
  //S(0,0)=deriv_avg.vals(3*n+3*n*n);
  for(int i=0; i < n; i++) { 
    S(0,i+1)=0.0;
    S(i+1,0)=0.0;
    //S(0,i+1)=deriv_avg.vals(n+i);
    //S(i+1,0)=deriv_avg.vals(n+i);
  }
  for(int i=0; i< n; i++) { 
    for(int j=0; j< n; j++) { 
      S(i+1,j+1)=deriv_avg.vals(3*n+i*n+j)-deriv_avg.vals(n+i)*deriv_avg.vals(n+j);
      //S(i+1,j+1)=deriv_avg.vals(3*n+i*n+j);
    }
  }
//  for(int i=0; i< n+1; i++) { 
//    cout <<"S-- ";
//    for(int j=0; j< n+1; j++) { 
//      cout << S(i,j) << " ";
//    }
//    cout << endl;
//  }
  H(0,0)=en(0);
  for(int i=0; i < n; i++) { 
    H(i+1,0)=(deriv_avg.vals(i)-en(0)*deriv_avg.vals(n+i));
    H(0,i+1)=(H(i+1,0)+deriv_avg.vals(2*n+i));
    //H(i+1,0)=deriv_avg.vals(i);
    //H(0,i+1)=deriv_avg.vals(2*n+i);
  }
  
  for(int i=0; i< n; i++) { 
    for(int j=0; j< n; j++) { 
      //int indx=3*n+2*n*n+i*n+j;

     // H(i+1,j+1)=deriv_avg.vals(3*n+2*n*n+i*n+j);
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

