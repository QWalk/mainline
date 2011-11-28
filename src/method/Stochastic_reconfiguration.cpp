
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
  cout << "alpha0 ";
  for(int i=0; i<nparms; i++) 
    cout << alpha(i) << " ";
  cout << endl;
  cout.precision(1);
  wavefunction_derivative(energies,S,en);
  doublevar tau_threshold=1e-3;
  int nit_completed=0;
  for(int it=0; it< iterations; it++) { 
    if(mpi_info.node==0) { 
    cout << "S " << endl;
    for(int i=0; i< nparms+1; i++) { 
      for(int j=0; j< nparms+1; j++) { 
        cout << S(i,j)  << " ";
      }
    }
    cout << endl;
    }


    Array1 <doublevar> save_alpha=alpha;
    x=0;
    doublevar tmp_tau=tau;
    bool energy_lowered=false;
    while(!energy_lowered) { 
      InvertMatrix(S,Sinv,nparms+1);
      
      for(int i=0; i <nparms+1; i++) { 
        for(int j=0; j< nparms+1; j++) { 
          x(i)+=Sinv(i,j)*(S(0,j)-tmp_tau*energies(j));
        }
      }

      cout << "alpha ";
      for(int i=0; i< nparms; i++) { 
        alpha(i)=save_alpha(i)+x(i+1)/x(0);
        cout << alpha(i) << " ";
      }
      cout << endl;
      wfdata->setVarParms(alpha);
      wavefunction_derivative(tmp_energies,tmp_S,tmp_en);
      cout << mpi_info.node << ":tmp_en " << tmp_en(0) << " en " << en(0) << "+/- " << en(1) << endl;
      if(mpi_info.node==0 and tmp_en(0) < -500) { 
        cout << "#######################################" << endl;
        cout << "S " << endl;
        for(int i=0; i< nparms+1; i++) { 
          for(int j=0; j< nparms+1; j++) { 
            cout << tmp_S(i,j)  << " ";
          }
        }
        cout << endl;
        cout << "alpha ";
        for(int i=0; i< nparms; i++) { 
          cout << alpha(i) << " ";
        }
        cout << endl;
        cout << "energies ";
        for(int i=0; i< nparms+1; i++) { 
          cout << tmp_energies(i) << " ";
        }
        cout << endl;
        exit(122);
      }

      if((tmp_en(0) < en(0) or fabs(tmp_en(0) - en(0)) < en(1)*2) and tmp_en(1) < 2*en(1)) {
        energy_lowered=true;
        energies=tmp_energies;
        en=tmp_en;
        S=tmp_S;
      }
      else { 
        wfdata->setVarParms(save_alpha);
        wavefunction_derivative(energies,S,en);
        tmp_tau*=0.5;
        tau=tmp_tau;
        output << "lowering tau to " << tau << endl;
        cout << "tmp_tau " << tmp_tau << endl;
      }
      if(tmp_tau < tau_threshold) { 
        output << "Tau is now below " << tau_threshold << ". Stopping" << endl;
        break;
      }
    }
    output << "energy " << energies(0) << " tau " << tau << endl;
    if(tau < tau_threshold) break;
    for(int i=0; i< 2; i++) { 
      energy_step(it,i)=en(i);
    }
    for(int i=0; i< nparms; i++) { 
      alpha_step(it,i)=alpha(i);
    }
    nit_completed++;
  }


  int navg=0;
  Array1 <doublevar> avg_alpha(nparms);
  avg_alpha=0;
  for(int  it=0; it < nit_completed; it++) { 
    if(energy_step(it,0)-en(0) < 2*en(1)) { 
      navg++;
      for(int i=0; i< nparms; i++) 
        avg_alpha(i)+=alpha_step(it,i);
    }
  }

  for(int i=0; i< nparms; i++) avg_alpha(i)/=navg;
  wfdata->setVarParms(avg_alpha);
  wavefunction_derivative(tmp_energies,tmp_S,tmp_en);
  output << "energy of average position " << tmp_en(0) << " +/- " << tmp_en(1) << endl;

  if(output) { 
    wfdata->renormalize();
    string indentation="";
    ofstream wfoutput(wfoutputfile.c_str());
    wfoutput.precision(15);
    wfdata->writeinput(indentation,wfoutput);
    wfoutput.close();
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

