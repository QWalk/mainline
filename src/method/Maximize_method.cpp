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


#include "Maximize_method.h"
#include "qmc_io.h"
#include "ulec.h"
#include <iomanip>
#include "Program_options.h"
#include "System.h"
#include "Split_sample.h"
#include "Properties.h"
#include "Generate_sample.h"

//----------------------------------------------------------------------

void Maximize_method::read(vector <string> words,
            unsigned int & pos,
            Program_options & options) { 
  allocate(options.systemtext[0],  sys);
  sys->generatePseudo(options.pseudotext, pseudo);

  debug_write(cout, "wfdata allocate\n");
  wfdata=NULL;
  if(options.twftext.size() < 1) error("Need TRIALFUNC section for OPTIMIZE");
  allocate(options.twftext[0], sys, wfdata);

  if(!readvalue(words,pos=0,nconfig,"NCONFIG"))
    nconfig=100;

  sample_only=false;
  if(haskeyword(words,pos=0,"SAMPLE"))
    sample_only=true;

  json_file = options.runid+".json";
  cout << "filename " << json_file << endl;
#ifdef USE_MPI
  int extra_configs = nconfig % mpi_info.nprocs;
  if(mpi_info.node < extra_configs) {
    nconfigs_per_node = (nconfig/mpi_info.nprocs) + 1;
  }
  else {
    nconfigs_per_node = (nconfig/mpi_info.nprocs);
  }
#else
  nconfigs_per_node = nconfig;
#endif
}

//----------------------------------------------------------------------
void Maximize_method::run(Program_options & options, ostream & output) { 
  Wavefunction * wf=NULL;
  Sample_point * sample=NULL;
  sys->generateSample(sample);
  wfdata->generateWavefunction(wf);
  sample->attachObserver(wf);
  Properties_gather mygather;
  Primary guidewf;

  if(sample_only) {
    run_sample(options, output);
    delete sample;
    delete wf;
    return;
  }
#ifdef USE_MPI
  if(mpi_info.nprocs<2) error("MAXIMIZE must be run with at least 2 processes to be run in parallel.");
  if(mpi_info.node==0) {
    master(wf, sample, output);
  }
  else {
    worker(wf, sample);
  }
#else
  run_old(options, output);
#endif //USE_MPI

  delete sample;
  delete wf;
}

int Maximize_method::gen_max_conf(Wavefunction * wf, Sample_point * sample,
    Config_save_point & config_pos, Maximize_config & maximize_config) {
  Primary guidewf;
  Properties_gather mygather;
  Properties_point pt;
  sys->generateSample(sample);
  wfdata->generateWavefunction(wf);
  sample->attachObserver(wf);
  
  Wf_return lap(1,5);
  int nelectrons=sample->electronSize();
  Array1 <doublevar> epos(3);
  Array2 <doublevar> tempconfig(nelectrons,3);
  Array1 <doublevar> tempgrad(3*nelectrons);
  Array2 <doublevar> temphessian(3*nelectrons,3*nelectrons);
  Array2 <doublevar> inverse_hessian(3*nelectrons,3*nelectrons);
  pseudo->setDeterministic(1); 

  config_pos.restorePos(sample);
  
  // Get initial info
  for(int e=0; e< nelectrons; e++) {
    sample->getElectronPos(e,epos);
    for(int d=0; d< 3; d++) {
      tempconfig(e,d) = epos(d);
    }
  }
  mygather.gatherData(pt, pseudo, sys, wfdata, wf, sample, &guidewf);
  
  maximize_config.psi_init = pt.wf_val.amp(0,0);
  maximize_config.energy_init = pt.energy(0);
  maximize_config.config_init = tempconfig;

  // Maximize Sample
  maximize(sample,wf,temphessian);
  
  // Get maximized info
  for(int e=0; e< nelectrons; e++) {
    sample->getElectronPos(e,epos);
    for(int d=0; d< 3; d++) {
      tempconfig(e,d) = epos(d);
    }
  }
  mygather.gatherData(pt, pseudo, sys, wfdata, wf, sample, &guidewf);
  
  // find gradient
  int count=0;
  doublevar psi_error=0;
  wf->updateLap(wfdata,sample);
  for(int e=0; e< nelectrons; e++) { 
    wf->getLap(wfdata,e,lap);
    for(int d=0; d< 3; d++) {
      doublevar grad=-lap.amp(0,d+1);
      tempgrad[count++]=grad;
    }
  }
  // estimate error in psi by 0.5 * g.T * H_inv * g
  InvertMatrix(temphessian, inverse_hessian, 3*nelectrons);
  for(int j=0; j<3*nelectrons; j++) {
    for(int k=0; k<3*nelectrons; k++) {
      psi_error += 0.5*tempgrad(j)*inverse_hessian(j,k)*tempgrad(k);
    }
  }

  maximize_config.nelectrons = nelectrons;
  maximize_config.psi = lap.amp(0,0);
  maximize_config.energy = pt.energy(0);
  maximize_config.config = tempconfig;
  maximize_config.hessian= temphessian;
  maximize_config.error = psi_error;
  
}

int Maximize_method::worker(Wavefunction * wf, Sample_point * sample) {
#ifdef USE_MPI
  Config_save_point tmpconfig;
  tmpconfig.mpiReceive(0);
  Maximize_config maximize_config;

  while(true) {
    gen_max_conf(wf, sample, tmpconfig, maximize_config);
    int done=1;
    MPI_Send(done,0);
    maximize_config.mpiSend(0);
    MPI_Recv(done,0);
    if(done==0) break;
    tmpconfig.mpiReceive(0);
  }
  cout << mpi_info.node <<  " : done " << endl;
#endif //USE_MPI
}

int Maximize_method::master(Wavefunction * wf, Sample_point * sample, ostream & output) {
#ifdef USE_MPI
  Array1 <Config_save_point> config_pos(nconfig);
  Properties_point pt;
  pt.setSize(1);
  Maximize_config maximize_config;
  MPI_Status status;

  Primary guidewf;
  generate_sample(sample,wf,wfdata,&guidewf,nconfig,config_pos);
  int nelectrons=sample->electronSize();
  maximize_config.set_nelectrons(nelectrons);
  int configcounter=0;
  int totcount=0;
  
  cout << "Writing to " << json_file << endl;
  ofstream os(json_file.c_str());
  stringstream tempstream;
  tempstream.precision(15);
  os << '[';
  //Get everyone started with data
  for(int r=1; r < mpi_info.nprocs; r++) {
    if(nconfig<r) {break;}
    config_pos(configcounter++).mpiSend(r);
  }

  while( configcounter < config_pos.GetDim(0)) {
    //Is anyone done?
    //When done, receive completed maximize_config and send out new point
    int done;
    MPI_Recv(&done,1,MPI_INT,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_Comm_grp,&status);
    done=1;
    maximize_config.mpiReceive(status.MPI_SOURCE);
    MPI_Send(done,status.MPI_SOURCE);
    config_pos(configcounter++).mpiSend(status.MPI_SOURCE);
    //write maximize_config
    tempstream.clear();
    tempstream.str("");
    maximize_config.write(tempstream,false);
    os << tempstream.str();
    os << ',';
    totcount++;
    if(totcount%(1)==0) cout << "Completed " << totcount << " samples " << endl;
  }

  //Loop through all the nodes and collect their last maximize_configs
  for(int r=1; r < mpi_info.nprocs; r++) {
    int done;
    MPI_Recv(done,r);
    done=0;
    maximize_config.mpiReceive(r);
    MPI_Send(done,r);
    maximize_config.write(os,false);
    if(r!=mpi_info.nprocs-1) {
      os << ',';
    }
    totcount++;
    if(totcount%(1)==0) cout << "Completed " << totcount << " samples " << endl;
  }
  os << ']';
  os.close();
#endif //USE_MPI
}

void Maximize_method::run_sample(Program_options & options, ostream & output) { 
  Wavefunction * wf=NULL;
  Sample_point * sample=NULL;
  sys->generateSample(sample);
  wfdata->generateWavefunction(wf);
  sample->attachObserver(wf);
  
  Wf_return lap(1,5);

  Array1 <Config_save_point> config_pos(nconfigs_per_node);
  Array1 <Maximize_config> maximize_config(nconfigs_per_node);
  Array1 <doublevar> epos(3);

  Primary guidewf;
  generate_sample(sample,wf,wfdata,&guidewf,nconfig,config_pos);
  int nelectrons=sample->electronSize();
  Properties_gather mygather;
  
  Array2 <doublevar> tempconfig(nelectrons,3);
  Array1 <doublevar> tempgrad(3*nelectrons);
  Array2 <doublevar> temphessian(3*nelectrons,3*nelectrons);
  pseudo->setDeterministic(1); 
  
  for(int i=0; i < nconfigs_per_node; i++) { 
    config_pos(i).restorePos(sample);
    stringstream tableout;
    
    // Get initial info
    for(int e=0; e< nelectrons; e++) {
      sample->getElectronPos(e,epos);
      for(int d=0; d< 3; d++) {
        tempconfig(e,d) = epos(d);
      }
    }

    Properties_point pt0;
    mygather.gatherData(pt0, pseudo, sys, wfdata, wf, 
                            sample, &guidewf);
    
    maximize_config(i).psi_init = pt0.wf_val.amp(0,0);
    maximize_config(i).energy_init = pt0.energy(0);
    maximize_config(i).config_init = tempconfig;
    
    maximize_config(i).nelectrons = nelectrons;
    maximize_config(i).psi = pt0.wf_val.amp(0,0);
    maximize_config(i).energy = pt0.energy(0);
    maximize_config(i).config = tempconfig;
    maximize_config(i).hessian= temphessian;
    maximize_config(i).error = 0;
     
    config_pos(i).write(cout);
    cout << "node" << mpi_info.node << " " << "sample " << i << " of " << nconfigs_per_node << " finished" << endl; 
  }
  
  write_configurations_maximize_json(json_file, maximize_config);
  delete wf;
  delete sample;
}
//----------------------------------------------------------------------


void Maximize_method::run_old(Program_options & options, ostream & output) { 
  Wavefunction * wf=NULL;
  Sample_point * sample=NULL;
  sys->generateSample(sample);
  wfdata->generateWavefunction(wf);
  sample->attachObserver(wf);
  
  Wf_return lap(1,5);

  Array1 <Config_save_point> config_pos(nconfigs_per_node);
  Array1 <Maximize_config> maximize_config(nconfigs_per_node);
  Array1 <doublevar> epos(3);

  Primary guidewf;
  generate_sample(sample,wf,wfdata,&guidewf,nconfig,config_pos);
  int nelectrons=sample->electronSize();
  Properties_gather mygather;
  
  Array2 <doublevar> tempconfig(nelectrons,3);
  Array1 <doublevar> tempgrad(3*nelectrons);
  Array2 <doublevar> temphessian(3*nelectrons,3*nelectrons);
  Array2 <doublevar> inverse_hessian(3*nelectrons,3*nelectrons);
  pseudo->setDeterministic(1); 
  
  for(int i=0; i < nconfigs_per_node; i++) { 
    config_pos(i).restorePos(sample);
    stringstream tableout;
    
    // Get initial info
    for(int e=0; e< nelectrons; e++) {
      sample->getElectronPos(e,epos);
      for(int d=0; d< 3; d++) {
        tempconfig(e,d) = epos(d);
      }
    }

    Properties_point pt0;
    mygather.gatherData(pt0, pseudo, sys, wfdata, wf, 
                            sample, &guidewf);
    
    maximize_config(i).psi_init = pt0.wf_val.amp(0,0);
    maximize_config(i).energy_init = pt0.energy(0);
    maximize_config(i).config_init = tempconfig;
    
    // Maximize Sample
    maximize(sample,wf,temphessian);
    
    // Get maximized info
    for(int e=0; e< nelectrons; e++) {
      sample->getElectronPos(e,epos);
      for(int d=0; d< 3; d++) {
        tempconfig(e,d) = epos(d);
      }
    }

    Properties_point pt;
    mygather.gatherData(pt, pseudo, sys, wfdata, wf, 
                            sample, &guidewf);
    
    // find gradient
    int count=0;
    doublevar psi_error=0;
    wf->updateLap(wfdata,sample);
    for(int e=0; e< nelectrons; e++) { 
      wf->getLap(wfdata,e,lap);
      for(int d=0; d< 3; d++) {
        doublevar grad=-lap.amp(0,d+1);
        tempgrad[count++]=grad;
      }
    }
    // estimate error in psi by 0.5 * g.T * H_inv * g
    InvertMatrix(temphessian, inverse_hessian, 3*nelectrons);
    for(int j=0; j<3*nelectrons; j++) {
      for(int k=0; k<3*nelectrons; k++) {
        psi_error += 0.5*tempgrad(j)*inverse_hessian(j,k)*tempgrad(k);
      }
    }

    maximize_config(i).nelectrons = nelectrons;
    maximize_config(i).psi = lap.amp(0,0);
    maximize_config(i).energy = pt.energy(0);
    maximize_config(i).config = tempconfig;
    maximize_config(i).hessian= temphessian;
    maximize_config(i).error = psi_error;
     
    config_pos(i).write(cout);
    cout << "node" << mpi_info.node << " " << "sample " << i << " of " << nconfigs_per_node << " finished" << endl; 
  }
  
  write_configurations_maximize_json(json_file, maximize_config);
  
  // Hessian steps
  // write Hessians calculated with different step sizes for the first configuration
  // int nsteps=10;
  // Array1 <Hessian_step> hessian_steps(nsteps);
  // hessian_vary_step(sample,wf,config_pos(0),hessian_steps);
  // write_hessian_vary_step_json(outfilename, hessian_steps,maximize_config(0).config,maximize_config(0).psi);
  ////

  delete wf;
  delete sample;
}
//----------------------------------------------------------------------

#include "macopt.h"
class Maximize_fn:public Macopt {
public:
  Maximize_fn(int _n,int _verbose, double _tol,
              int _itmax,int _rich):Macopt(_n,_verbose,_tol,_itmax,_rich) { 
              } 
  ~Maximize_fn(){ };


  //-----------------------------------------
  void update_positions(double * _p) { 
    int count=1; 
    Array1 <doublevar> epos(3);
    int nelec=sample->electronSize();
    for(int e=0; e< nelec; e++) { 
      for(int d=0; d< 3; d++) { 
        epos(d)=_p[count++];
      }
      sample->setElectronPos(e,epos);
    }
    
  }
  //-----------------------------------------

  double func(double * _p) { 
    update_positions(_p);
    wf->updateVal(wfdata,sample);
    Wf_return wfret(1,2);
    wf->getVal(wfdata,0,wfret);
    return -wfret.amp(0,0);
  }
  //-----------------------------------------
  
  double dfunc(double * _p, double * _g) { 
    update_positions(_p);
    wf->updateLap(wfdata,sample);
    int count=0;
    int nelec=sample->electronSize();
    Wf_return lap(1,5);
    wf->getVal(wfdata,0,lap);
    //doublevar cutoff=0.2;
    _g[count++]=-lap.amp(0,0);
    for(int e=0; e< nelec; e++) { 
      wf->getLap(wfdata,e,lap);
      for(int d=0; d< 3; d++) {
        doublevar grad=-lap.amp(0,d+1);
        //if(abs(grad) > cutoff) {
        //  cout << "big grad " << grad << endl;
        //  grad=sign(grad)*cutoff;
        //}
        _g[count++]=grad;
      }
    }
    return _g[0];
  }
  //-----------------------------------------
  //move along the gradient for time tstep and evaluate the function
  doublevar eval_tstep(double * x, doublevar tstep,Array1 <doublevar> & grad,
                       int n) {
    Array1 <doublevar> xnew(n+1);
    for(int i=1; i<=n; i++)
      xnew(i)=x[i]-tstep*grad(i);
    return func(xnew.v);
  }

  //-----------------------------------------
  // check if the magnitude of the gradient is within tol
  bool grad_under_tol(Array1 <doublevar> & grad, int n, doublevar tol) {
    doublevar total = 0;
    for(int i=1; i<=n; i++)
      total += grad(i)*grad(i);
    return total < tol*tol;
  }

  doublevar grad_abs(Array1 <doublevar> & grad, int n) {
    doublevar total = 0;
    for(int i=1; i<=n; i++)
      total += grad(i)*grad(i);
    return sqrt(total);
  }
  
  //-----------------------------------------
  void macoptII(double * x,int n) {
    Array1 <doublevar> grad(n+1);
    Array1 <doublevar> prev_grad(n+1);
    Array1 <doublevar> step_dir(n+1);
    //Array1 <doublevar> xnew(n+1);
    stringstream os;
    os.precision(15);
    int max_it = 100;
    int max_big_it = n*5;
    int beta_counter=0;
    doublevar tol = 1e-9;
    doublevar outer_tol = 1e-4; // make outer loop tolerance looser than inner
    doublevar beta;
    os << "node" << mpi_info.node << " max_it=" << max_it << " tol=" << tol <<" max_big_it=" << max_big_it << " outer_tol=" << outer_tol << endl;
    for(int big_it=0; big_it < max_big_it; big_it++) {
      //find the direction of the gradient at x
      //*prev_grad.v = *grad.v;
      for(int i=1; i<=n;i++) {
        prev_grad(i) = grad(i);
      }
      dfunc(x,grad.v);
      doublevar gradlen = grad_abs(grad,n);
      doublevar num=0;
      doublevar den=0;
      //find conjugate gradient direction
      if(beta_counter%n==0) {
        os << "node" << mpi_info.node << " big_it=" << big_it << " reset beta_counter=" << beta_counter << endl;
        beta = 0.0;
      } else {
        for(int i=1; i<=n; i++) {
          num += grad(i)*(grad(i)-prev_grad(i));
          den += prev_grad(i)*prev_grad(i);
        }
        beta = num/den;
        if(beta<0) {
          os << "node" << mpi_info.node << " big_it=" << big_it << " reset beta=" << beta << endl;
          beta = 0.0;
          beta_counter = 0;
        }
      }
      beta_counter++;
      for(int i=1; i<=n; i++) {
        step_dir(i) = grad(i) + beta*step_dir(i);
      }
      
      //bracket the minimum in this direction (in 1D units of tstep) for bisection method
      if (gradlen < outer_tol){ os << "node" << mpi_info.node << " func " << grad(0) << " gradlen " << gradlen << endl; break;}
      doublevar fbase=grad(0);
      doublevar bracket_tstep=0.0,last_func=fbase;
      doublevar f_tol = 1e-14;
      for(doublevar tstep=1e-4; tstep < 1e3; tstep*=2.0) {
        doublevar f=eval_tstep(x,tstep,step_dir,n);
        if(f > fbase+f_tol or f > last_func+f_tol) {
          bracket_tstep=tstep;
          last_func=f;
          break;
        }
        else last_func=f;
      }
      os << "node" << mpi_info.node << " big_it=" << big_it << " bracket_tstep=" << bracket_tstep << " gradlen=" << gradlen << " fbase=" << fbase << " f=" << last_func << " f-fbase=" << last_func-fbase << endl;
      if (bracket_tstep==0){os << "node" << mpi_info.node << " gradient too small, bracket_tstep=0, exiting loop" << endl; break;}
      
      //bisection method works best using the golden ratio
      doublevar resphi=2.-(1.+sqrt(5.))/2.;
      doublevar a=0, b=resphi*bracket_tstep, c=bracket_tstep;
      doublevar af=fbase,
      bf=eval_tstep(x,b,step_dir,n),
      cf=eval_tstep(x,c,step_dir,n);
      doublevar unit;
      os << "node" << mpi_info.node << " big_it=" << big_it << " " << "first step  a,b,c " 
      << a << " " << b << "  " << c << " funcs " << af << " " << bf << " " << cf << endl;
      //bisection method iteration, (a, b, c)
      for(int it=0; it < max_it; it++) {
        doublevar d,df;
        //choose point d on the bigger interval (a,b) or (b,c)
        if( (c-b) > (b-a))
          d=b+resphi*(c-b);
        else
          d=b-resphi*(b-a);
        df=eval_tstep(x,d,step_dir,n);
        //check if the function is smaller at d than at b to choose next bracket
        if(df < bf) {
          if( (c-b) > (b-a) ) {
            // How to check tolerance: Want to make sure to check relative change in psi. 
            // grad psi / psi = grad log psi ~ (log psi_a-log psi_b)/(a-b) ~ log(psi_a/psi_b)/(a-b) ~ (psi_a-psi_b)/psi_b/(a-b) (since log(1+x) ~ x).
            //(af-bf) = log psi_a - log psi_b = log(psi_a/psi_b) ~ (psi_a-psi_b)/psi_b < tol
            //gradlen is abs(grad) at start, timesteps a, b, c, d are in units of grad, need to multiply to get actual distance a-b.
            unit = gradlen*(b-a);
            if(((af-bf)/unit<tol and (bf-af)/unit<tol) or unit<1e-15) { os << "node" << mpi_info.node << " big_it=" << big_it << " break step " << it << " (af-bf)/(b-a)=" << af/unit-bf/unit << " unit=" << unit << endl; break; }
            a=b;
            af=bf;
            b=d;
            bf=df;
          }
          else {
            unit = gradlen*(c-b);
            if(((cf-bf)/unit<tol and (bf-cf)/unit<tol) or unit<1e-15) { os << "node" << mpi_info.node << " big_it=" << big_it << " break step " << it << " (cf-bf)/(c-b)=" << cf/unit-bf/unit << " unit=" << unit << endl; break; }
            c=b;
            cf=bf;
            b=d;
            bf=df;
          }
        }
        //if function is bigger at d than at b, make d the new bracket boundary
        else {
          if( (c-b) > (b-a) ) {
            unit = gradlen*(c-d);
            if(((cf-df)/unit<tol and (df-cf)/unit<tol) or unit<1e-15) { os << "node" << mpi_info.node << " big_it=" << big_it << " break step " << it << " (df-cf)/(d-c)=" << df/unit-cf/unit << " unit=" << unit << endl; break; }
            c=d;
            cf=df;
          }
          else {
            unit = gradlen*(d-a);
            if(((af-df)/unit<tol and (df-af)/unit<tol) or unit<1e-15) { os << "node" << mpi_info.node << " big_it=" << big_it << " break step " << it << " (af-df)/(d-a)=" << af/unit-df/unit << " unit=" << unit << endl; break; }
            a=d;
            af=df;
          }
        }
        if(it==max_it-1) {
          os << "node" << mpi_info.node << " big_it=" << big_it << " " << "Warning: inner loop did not reach tolerance" << endl;
        }
      }
      os << "node" << mpi_info.node << " big_it=" << big_it << " " << "last step b-a,c-b " << b-a << " " << c-b << " func diffs " << (af-bf)/((b-a)*gradlen) << " " << (cf-bf)/((c-b)*gradlen) << " " << endl;
      //finished bisection search, minimum at b; compute x for tstep b
      doublevar best_tstep=b;
      for(int i=1; i<=n; i++)
        x[i]=x[i]-best_tstep*step_dir[i];
      os << "node" << mpi_info.node << " big_it=" << big_it << " took step of size " << best_tstep*grad_abs(step_dir,n) << endl;
      
      if(big_it==max_big_it-1) {
        os << "node" << mpi_info.node << " " << "Warning: outer loop did not reach tolerance." << endl;
      }
    }
    newton_iteration(x, n, 0, os);
    cout << os.str(); 
  }
  
  //-----------------------------------------
  void newton_iteration(double * x, int n, int nit, stringstream& os) {
    Array2 <doublevar> hessian(n,n);
    Array2 <doublevar> inverse_hessian(n,n);
    Array1 <doublevar> grad(n+1);
    Array1 <doublevar> delta_x(n+1);
    for (int it=0; it<=nit; it++) {
      doublevar psi_error = 0;
      dfunc(x,grad.v);
      calc_hessian(x, hessian, n, 1e-6);
      InvertMatrix(hessian, inverse_hessian, n);
      doublevar grad_squared = 0, delta_x_squared = 0;
      for (int i=1;i<=n;i++){
          delta_x(i)=0;
        for (int j=1;j<=n;j++) {
          delta_x(i) += inverse_hessian(i-1,j-1)*grad(j); // Estimate x error: dx = H_inv * g
          psi_error += 0.5*grad(i)*inverse_hessian(i-1,j-1)*grad(j); // Estimate psi error: dpsi = g.T * x / 2 = g.T * H_inv * g / 2
        }
        grad_squared += grad(i)*grad(i);
        delta_x_squared += delta_x(i)*delta_x(i);
      }
      os << "node" << mpi_info.node << " " << "newton" << it << ", |grad| = " << sqrt(grad_squared);
      os << ", psi error = " << psi_error << ", x error = " << sqrt(delta_x_squared) << endl;
      if(it<nit){ // don't update x on the last step, so that the ouput represents the final estimated error (and not one that has been corrected)
        for(int i=1; i<=n; i++)
          x[i]=x[i]-delta_x(i);
      }
    }
  }

  //-----------------------------------------
  void iteration_print(double f, double gg, double tol,  int itn) {
    cout << "It " << endl;
  }
  
  //-----------------------------------------
  void calc_hessian(double * x, Array2 <doublevar> & hessian, int n, double h=1e-6) {
    Array1 <doublevar> grad_plus(n+1);
    Array1 <doublevar> grad_minus(n+1);
    Array1 <doublevar> xnew(n+1);
    if(h<=0) { h = 1e-3; } // default step size 
    
    for(int i=1; i<=n; i++) 
      xnew(i)=x[i];
    
    for(int i=1; i<=n; i++) {
      // forward step
      xnew(i) = x[i] + h;
      dfunc(xnew.v,grad_plus.v);
      // backward step
      xnew(i) = x[i] - h;
      dfunc(xnew.v,grad_minus.v);
      // reset position
      xnew(i) = x[i];
      // compute derivative
      for(int j=1; j<=n; j++) {
        hessian(i-1,j-1) = (grad_plus(j) - grad_minus(j)) / (2*h);
      }   
    }
  }
  //-----------------------------------------

  Wavefunction * wf;
  Wavefunction_data * wfdata;
  Sample_point * sample;
  System * system;
};

//----------------------------------------------------------------------

void Maximize_method::maximize(Sample_point * sample,Wavefunction * wf,Array2 <doublevar> & hessian) { 
  int nelectrons=sample->electronSize();
  Maximize_fn maximizer(nelectrons*3,1,1e-12,1000,1);
  maximizer.wf=wf;
  maximizer.wfdata=wfdata;
  maximizer.sample=sample;
  maximizer.system=sys;

  int count=1;
  Array1 <doublevar> allpos(nelectrons*3+1);
  Array1 <doublevar> epos(3);
  for(int e=0; e< nelectrons; e++) {
    sample->getElectronPos(e,epos);
    for(int d=0; d< 3; d++) { 
      allpos(count++)=epos(d);
    }
  }

  maximizer.maccheckgrad(allpos.v,nelectrons*3,0.001,nelectrons*3);
  maximizer.macoptII(allpos.v,nelectrons*3);  
  maximizer.calc_hessian(allpos.v,hessian,nelectrons*3);
}

//----------------------------------------------------------------------

void Maximize_method::hessian_vary_step(Sample_point * sample,Wavefunction * wf,Config_save_point & pt,Array1 <Hessian_step> & hessian_steps) { 
  int nelectrons=sample->electronSize();
  int nsteps=hessian_steps.GetDim(0);
  pt.restorePos(sample);
  Maximize_fn maximizer(nelectrons*3,1,1e-12,1000,1);
  maximizer.wf=wf;
  maximizer.wfdata=wfdata;
  maximizer.sample=sample;
  maximizer.system=sys;

  int count=1;
  Array2 <doublevar> temphessian(nelectrons*3,nelectrons*3);
  Array1 <doublevar> allpos(nelectrons*3+1);
  Array1 <doublevar> epos(3);
  for(int e=0; e< nelectrons; e++) {
    sample->getElectronPos(e,epos);
    for(int d=0; d< 3; d++) { 
      allpos(count++)=epos(d);
    }
  }
  // compute hessians using different step sizes
  for(int i=0;i<nsteps;i++) {
    hessian_steps(i).step = pow(10.0,-i-1);
    maximizer.calc_hessian(allpos.v,temphessian,nelectrons*3,hessian_steps(i).step);
    hessian_steps(i).hessian = temphessian;
  }
}

//----------------------------------------------------------------------
  
int Maximize_method::showinfo(ostream & os) { 
  return 1;
}

//----------------------------------------------------------------------
void write_hessian_vary_step_json(string & outfilename, Array1 <Hessian_step> hessian_steps,Array2 <doublevar> config, doublevar psi) {
  int nsteps = hessian_steps.GetDim(0);
  int nelectrons = config.GetDim(0);
  // write to file
  ofstream os(outfilename.c_str());
  os.precision(15);
  os << "[{";
  os << "\"psi\": " << psi << ", ";
  os << "\"config\": " << "[";
  for(int e=0; e<nelectrons; e++) {
    if(e!=0) { os << ", "; }
    os << "[";
    for(int d=0; d<3; d++) {
      if(d!=0) { os << ", "; }
      os << config(e,d); 
    }
    os << "]";
  }
  os << "], ";
  os << "\"steps\": " << "[";
  for(int i=0;i<nsteps;i++) {
    if(i!=0) { os << ", "; }
    os << hessian_steps(i).step;
  }
  os << "], ";
  os << "\"hessians\": " << "[";
  for(int i=0;i<nsteps;i++) {
    if(i!=0) { os << ", "; }
    os << "[";
    for(int j=0;j<3*nelectrons;j++) {
      if(j!=0) { os << ", "; }
      os << "[";
      for(int k=0;k<3*nelectrons;k++) {
        if(k!=0) { os << ", "; }
          os << hessian_steps(i).hessian(j,k);
      }
      os << "]";
    }
    os << "]";
  }
  os << "]";
  os << "}]";
  os.close();
}

//----------------------------------------------------------------------
void write_configurations_maximize_json(string & filename, Array1 <Maximize_config> configs) { 
  int nconfigs=configs.GetDim(0);
  time_t starttime;
  time(&starttime);
  string tmpfilename=filename; //+".qw_tomove";
  string backfilename=filename+".backup";
  if(mpi_info.node==0) { rename(tmpfilename.c_str(),backfilename.c_str()); }

#ifdef USE_MPI
  stringstream os;
  os.precision(15);
  for(int i=0; i< nconfigs; i++) {
    if(i!=0) { os << ", "; }
    configs(i).write(os,false);
  }

  string walkstr=os.str();
  int nthis_string=walkstr.size();
  if(mpi_info.node==0) {
    ofstream os(tmpfilename.c_str());
    os << "[";
    os << walkstr;
    MPI_Status status;
    for(int i=1; i< mpi_info.nprocs; i++) {
      MPI_Recv(nthis_string,i);
      char * buf=new char[nthis_string+1];
      MPI_Recv(buf,nthis_string,MPI_CHAR, i, 0, MPI_Comm_grp, & status);
      buf[nthis_string]='\0';
      os << ", " << buf;
      delete [] buf;
    }
    os << "]";
    os.close();
  }
  else {
    MPI_Send(nthis_string,0);
    //we know that MPI_Send obeys const-ness, but the interfaces are not clean 
    // and so...casting!
    MPI_Send((char *) walkstr.c_str(),nthis_string, MPI_CHAR, 0,0,MPI_Comm_grp);
  }
#else
  ofstream os(tmpfilename.c_str());
  os.precision(15);
  os << "[";
  for(int i=0; i< nconfigs; i++) {
    if(i!=0) { os << ", "; }
    configs(i).write(os,false);
  }
  os << "]";
  os.close();
#endif
  time_t endtime;
  time(&endtime);
  if(mpi_info.node==1) {
    debug_write(cout, "Write took ", difftime(endtime, starttime), " seconds\n");
  }
}


void write_configurations_maximize_yaml(string & filename, Array1 <Maximize_config> configs) { 
  int nconfigs=configs.GetDim(0);
  time_t starttime;
  time(&starttime);
  string tmpfilename=filename; //+".qw_tomove";
  string backfilename=filename+".backup";
  if(mpi_info.node==0) { rename(tmpfilename.c_str(),backfilename.c_str()); }

#ifdef USE_MPI
  stringstream os;
  os.precision(15);
  for(int i=0; i< nconfigs; i++) {
     configs(i).write(os,false);
  }

  string walkstr=os.str();
  int nthis_string=walkstr.size();
  if(mpi_info.node==0) {
    ofstream os(tmpfilename.c_str());
    os << "results:" << endl;
    os << walkstr;
    MPI_Status status;
    for(int i=1; i< mpi_info.nprocs; i++) {
      MPI_Recv(nthis_string,i);
      char * buf=new char[nthis_string+1];
      MPI_Recv(buf,nthis_string,MPI_CHAR, i, 0, MPI_Comm_grp, & status);
      buf[nthis_string]='\0';
      os << buf;
      delete [] buf;
    }
  }
  else {
    MPI_Send(nthis_string,0);
    //we know that MPI_Send obeys const-ness, but the interfaces are not clean 
    // and so...casting!
    MPI_Send((char *) walkstr.c_str(),nthis_string, MPI_CHAR, 0,0,MPI_Comm_grp);
  }
#else
  ofstream os(tmpfilename.c_str());
  os.precision(15);
  os << "results:" << endl;
  for(int i=0; i< nconfigs; i++) {
    configs(i).write(os,false);
  }
  os.close();
#endif
  time_t endtime;
  time(&endtime);
  if(mpi_info.node==1) {
    debug_write(cout, "Write took ", difftime(endtime, starttime), " seconds\n");
  }
}

//----------------------------------------------------------------------
void write_configurations_maximize(string & filename, 
                              Array1 <Maximize_config> configs) { 
  int nconfigs=configs.GetDim(0);
  time_t starttime;
  time(&starttime);
  string tmpfilename=filename; //+".qw_tomove";
  string backfilename=filename+".backup";
  if(mpi_info.node==0) { rename(tmpfilename.c_str(),backfilename.c_str()); }

#ifdef USE_MPI
  stringstream os;
  os.precision(15);
  for(int i=0; i< nconfigs; i++) {
     //os << " walker { \n";
     configs(i).write(os,false);
     //os << " } \n";
   }

  string walkstr=os.str();
  int nthis_string=walkstr.size();
  if(mpi_info.node==0) {
    ofstream os(tmpfilename.c_str());
    for(int e=0; e< configs(0).nelectrons; e++) {
      for(int d=0; d< 3; d++) {
        os << "e" << e << "_"<< d << " ";
      }
    }
    os << "elocal" << " ";
    os << "psi" << endl;
    
    os << walkstr;
    MPI_Status status;
    for(int i=1; i< mpi_info.nprocs; i++) {
      MPI_Recv(nthis_string,i);
      char * buf=new char[nthis_string+1];
      MPI_Recv(buf,nthis_string,MPI_CHAR, i, 0, MPI_Comm_grp, & status);
      buf[nthis_string]='\0';
      os << buf;
      delete [] buf;
    }
  }
  else {
    MPI_Send(nthis_string,0);
    //we know that MPI_Send obeys const-ness, but the interfaces are not clean 
    // and so...casting!
    MPI_Send((char *) walkstr.c_str(),nthis_string, MPI_CHAR, 0,0,MPI_Comm_grp);
  }
#else
    ofstream os(tmpfilename.c_str());
    os.precision(15);
    for(int i=0; i< nconfigs; i++) {
      //os << " walker { \n";
      configs(i).write(os,false);
      //os << " } \n";
    }
    os.close();
#endif
    time_t endtime;
    time(&endtime);
    if(mpi_info.node==1) {
      debug_write(cout, "Write took ", difftime(endtime, starttime), " seconds\n");
    }
}

//---------------------------------------------------------------------
//YAML format
//void Maximize_config::write(ostream & os) {
//  os << "- " << "psi: " << psi << endl;
//  os << "  " << "energy: " << energy << endl;
//  os << "  " << "config: " << endl;
//  for(int e=0; e<nelectrons; e++) {
//    os << "  - - " << config(e,0) << endl;
//    for(int d=1; d<3; d++) {
//      os << "    - " << config(e,d) << endl;
//    }
//  }
//  os << "  " << "hessian: " << endl;
//  for(int c=0; c<3*nelectrons; c++) {
//    os << "  - - " << hessian(c,0) << endl;
//    for(int d=1; d<3*nelectrons; d++) {
//      os << "    - " << hessian(c,d) << endl;
//    }
//  }
//}

// JSON format
void Maximize_config::write(ostream & os,bool write_hessian) {
  os << "{";
  os << "\"psi\": " << psi;
  os << "," << "\"psi_init\": " << psi_init;
  os << "," << "\"error\": " << error;
  os << "," << "\"energy\": " << energy;
  os << "," << "\"energy_init\": " << energy_init;
  os << "," << "\"config\": " << "[";
  for(int e=0; e<nelectrons; e++) {
    if(e!=0) { os << ", "; }
    os << "[";
    for(int d=0; d<3; d++) {
      if(d!=0) { os << ", "; }
      os << config(e,d); 
    }
    os << "]";
  }
  os << "]";
  os << "," << "\"config_init\": " << "[";
  for(int e=0; e<nelectrons; e++) {
    if(e!=0) { os << ", "; }
    os << "[";
    for(int d=0; d<3; d++) {
      if(d!=0) { os << ", "; }
      os << config_init(e,d); 
    }
    os << "]";
  }
  os << "]";
  if(write_hessian){
    os << "," << "\"hessian\": " << "[";
    for(int c=0; c<3*nelectrons; c++) {
      if(c!=0) { os << ", "; }
      os << "[";
      for(int d=0; d<3*nelectrons; d++) {
        if(d!=0) { os << ", "; }
        os << hessian(c,d); 
      }
      os << "]";
    }
    os << "]";
  }
  os << "}";
}

void Maximize_config::read(istream & is) {}
void Maximize_config::mpiSend(int node) {
#ifdef USE_MPI
  MPI_Send(&psi, 1, MPI_DOUBLE, node, 0, MPI_Comm_grp);
  MPI_Send(&psi_init, 1, MPI_DOUBLE, node, 0, MPI_Comm_grp);
  MPI_Send(&error, 1, MPI_DOUBLE, node, 0, MPI_Comm_grp);
  MPI_Send(&energy, 1, MPI_DOUBLE, node, 0, MPI_Comm_grp);
  MPI_Send(&energy_init, 1, MPI_DOUBLE, node, 0, MPI_Comm_grp);
  MPI_Send(config.v, config.GetSize(), MPI_DOUBLE, node, 0, MPI_Comm_grp);
  MPI_Send(config_init.v, config_init.GetSize(), MPI_DOUBLE, node, 0, MPI_Comm_grp);
  MPI_Send(hessian.v, hessian.GetSize(), MPI_DOUBLE, node, 0, MPI_Comm_grp);
#else
  //error("Maximize_config::mpiSend: not using MPI, this is most likely a bug");
#endif
}
void Maximize_config::mpiReceive(int node) {
#ifdef USE_MPI
  MPI_Status status;
  MPI_Recv(&psi, 1, MPI_DOUBLE, node, 0, MPI_Comm_grp, &status);
  MPI_Recv(&psi_init, 1, MPI_DOUBLE, node, 0, MPI_Comm_grp, &status);
  MPI_Recv(&error, 1, MPI_DOUBLE, node, 0, MPI_Comm_grp, &status);
  MPI_Recv(&energy, 1, MPI_DOUBLE, node, 0, MPI_Comm_grp, &status);
  MPI_Recv(&energy_init, 1, MPI_DOUBLE, node, 0, MPI_Comm_grp, &status);
  MPI_Recv(config.v, config.GetSize(), MPI_DOUBLE, node, 0, MPI_Comm_grp, &status);
  MPI_Recv(config_init.v, config_init.GetSize(), MPI_DOUBLE, node, 0, MPI_Comm_grp, &status);
  MPI_Recv(hessian.v, hessian.GetSize(), MPI_DOUBLE, node, 0, MPI_Comm_grp, &status);
#else
  //error("Maximize_config::mpiReceive: not using MPI, this is most likely a bug");
#endif
}
void Maximize_config::set_nelectrons(int nelec) {
  nelectrons = nelec;
  config.Resize(nelectrons,3);
  config_init.Resize(nelectrons,3);
  hessian.Resize(nelectrons*3,nelectrons*3);
}
