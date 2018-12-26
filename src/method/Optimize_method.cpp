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
#include "Optimize_method.h"
#include "qmc_io.h"
#include "Program_options.h"
#include "ulec.h"
#include "System.h"
#include "Generate_sample.h"
void Optimize_method::read(vector <string> words,
                           unsigned int & pos,
                           Program_options & options)
{

  int totconfig=2048;
  readvalue(words,pos=0,totconfig,"TOT_CONFIG");

  if(!readvalue(words, pos=0, nconfig, "NCONFIG"))
    nconfig=max(totconfig/mpi_info.nprocs,1);
  if(!readvalue(words,pos=0, eref, "EREF"))
    guess_eref=1;
  else guess_eref=0;
  

  //Optional options

  if(!readvalue(words,pos=0, iterations, "ITERATIONS")) 
    iterations=30;
  

  pos=0;
  if( ! readvalue(words, pos=0, pseudostore, "PSEUDOTEMP") )
    pseudostore=options.runid+".pseudo";

  if(! readvalue(words, pos=0, wfoutputfile, "WFOUTPUT") )
      wfoutputfile=options.runid+".wfout";
  
  if(haskeyword(words, pos=0, "USE_WEIGHTS") )
  {
    use_weights=1;
  }
  else { use_weights=0; }
  if(haskeyword(words, pos=0, "EXTENDED_WFOUTPUT") )
  {
    use_extended_output=1;
  }
  else { use_extended_output=0; }
  if(!readvalue(words,pos=0, mixing, "MIXING"))
  {
    mixing=0.95;
  }

  string functiontype_str;
  pos=0;
  if(readvalue(words, pos, functiontype_str, "MINFUNCTION")) {
    if(caseless_eq(functiontype_str, "VARIANCE"))
      min_function=min_variance;
    else if(caseless_eq(functiontype_str,"ABSOLUTE"))
      min_function=min_abs;
    else if(caseless_eq(functiontype_str,"LORENTZ"))
      min_function=min_lorentz;
    else if(caseless_eq(functiontype_str,"ENERGY")) 
      min_function=min_energy;
    else if(caseless_eq(functiontype_str,"MIXED")) 
      min_function=min_mixed;
    else
      error("I don't know ", functiontype_str, " for MINFUNCTION.");
  }
  else {
    min_function=min_variance;
  }

  readvalue(words, pos=0, readconfig, "READCONFIG");
  update_psp=0;
  if(haskeyword(words, pos=0, "UPDATE_PSP")) update_psp=1;

  //--Set up variables
  allocate(options.systemtext[0],  sys);
  sys->generatePseudo(options.pseudotext, pseudo);
  sys->generateSample(sample);

  debug_write(cout, "wfdata allocate\n");
  wfdata=NULL;
  if(options.twftext.size() < 1) error("Need TRIALFUNC section for OPTIMIZE");
  allocate(options.twftext[0], sys, wfdata);

  if(wfdata->nparms() <= 0 ) 
    error("There appear to be no parameters to optimize!");

}

int Optimize_method::showinfo(ostream & os)
{
  os << "     System " << endl;
  sys->showinfo(os);
  os << endl << endl;
  os << "     Wavefunction " << endl;
  wfdata->showinfo(os);
  os << endl << endl;
  pseudo->showinfo(os);
  os << endl << endl;
  os << "-----------------------------" << endl;
  os << "Wave function optimization:  " << endl;
  
  os << "Number of configurations : " << nconfig*mpi_info.nprocs << endl;
  os << "Maximum iterations: "  << iterations
     << endl;
  os << "nparms " << wfdata->nparms() << endl;
  os << "Minimization function : ";
  switch(min_function) {
  case min_variance:
    os << "Variance \n";
    break;
  case min_abs:
    os << "Absolute value\n";
    break;
  case min_lorentz:
    os << "Lorentz\n";
    break;
  case min_energy:
    os << "Energy\n";
    break;
  case min_mixed:
    os << "Mixed energy and variance\n";
    break;
  default:
    os << "Unknown--add to showinfo()\n";
  }
 
  if(use_weights) 
    os << "Reweighting using correlated sampling " << endl;
  os << "----------------------------" << endl;
  
  
  return 1;
}

/*!
*/
void Optimize_method::run(Program_options & options, ostream & output)
{
  int nparms=wfdata->nparms(); //Number of variables
  
  Array1 <double> x(nparms+1); //Parameters
  double f=0;  //Best value of the function
 
  Array1 <doublevar> temp_parms(nparms);
  wfdata->getVarParms(temp_parms);
  lastparms.Resize(nparms);
  lastparms=temp_parms;
  cout.precision(16);

  

  //cout << "pseudopotential " << endl;

   
  psp_buff.clear();

  //Here we save the original values (in case we do weighted correlated sampling)
  // and save the values for pseudopotential evaluation, in case the wf doesn't have
  //analytic derivatives

  psp_test.Resize(nconfig);
  orig_vals.Resize(nconfig);
  local_energy.Resize(nconfig);
  wfdata->generateWavefunction(wf);
  sample->attachObserver(wf);

  if(readconfig=="") {
    Primary guidewf;
    generate_sample(sample,wf,wfdata,&guidewf,nconfig,config_pos);
  }
  else read_configurations(readconfig,config_pos);
  int configsread=config_pos.GetDim(0);


  if(configsread < nconfig)
  {
    nconfig=configsread;
    cout << "processor " << mpi_info.node << " : "
    << "WARNING: Didn't find enough configurations in the "
    << "file.  Running optimization with only " << nconfig
    << " sample points." << endl;
  }




  for(int walker=0; walker < nconfig; walker++)  {
    config_pos(walker).restorePos(sample);
    nfunctions=wf->nfunc();
    orig_vals(walker).Resize(nfunctions, 2);
    wf->updateLap(wfdata, sample);
    wf->getVal(wfdata, 0, orig_vals(walker));
    psp_test(walker).Resize(pseudo->nTest());
    for(int i=0; i< pseudo->nTest(); i++) { 
      psp_test(walker)(i)=rng.ulec();
    }
    local_energy(walker)=sys->calcLoc(sample);
    if(!update_psp) {
      Array1 <doublevar> nonloc(wf->nfunc());
      pseudo->calcNonlocWithTest(wfdata,sys, sample, wf,psp_test(walker),nonloc );    
      local_energy(walker)+=nonloc(0);
    }
    
  }
  nfunctions=wf->nfunc();
    
  for(int i=0; i< nparms; i++) {
    x(i+1)=temp_parms(i);
  }
  
  if(guess_eref) { 
    eref=variance(nparms,x,f,0);
    output << "Running with reference energy " << eref << endl;
  }

  //cout << nparms << "   " << sqrt(f) << " " << dfn << "  " << eps << "   "
  //    << maxfn << "   " << iexit << endl;

  int verbose=1;
  double tolerance=0.0001;
  int rich=0;
  Optimize_fn optimizer(nparms,verbose,tolerance,iterations,rich);
  optimizer.opt_method=this;
  optimizer.output= &output;
//optimizer.maccheckgrad(x.v,nparms, .0001, nparms);
  optimizer.macoptII(x.v,nparms);  
  variance(nparms, x,f);

  //remove(pseudostore.c_str());
  if(output) {
    output << "Optimization finished.  ";
    output << "New wave function is in " << wfoutputfile << endl;
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

  //Put the wave function into options.trialfunc so further methods can use it.
  ifstream wfinput(wfoutputfile.c_str());
  options.twftext[0].clear();
  parsefile(wfinput, options.twftext[0]);
  wfinput.close();
}


//----------------------------------------------------------------------

/*!
\bug
When optimizing with weights, we need to renormalize the function every
time; otherwise the weights blow up.

This also returns the average energy..
 */
doublevar Optimize_method::variance(int n, Array1 <double> & parms, double & val, int check)
{

  int nparms=wfdata->nparms();
  Array1 <doublevar> temp_parms(nparms);

  for(int i=0; i< nparms; i++)  {
    temp_parms(i)=parms(i+1);
  }

  Wf_return  wfval(nfunctions,2);

  wfdata->setVarParms(temp_parms);

  Array1 <doublevar> variance(nfunctions);
  variance=0;
  Array1 <doublevar> weightsum(nfunctions);
  weightsum=0;
  Array1 <doublevar> kinetic(nfunctions), nonloc(nfunctions);
  //psp_buff.start_again();

  doublevar avgen=0;
  doublevar avgkin=0;
  doublevar avgpot=0;
  doublevar avgnon=0;
  //cout << "loop " << endl;
  for(int walker=0; walker< nconfig; walker++) {
    config_pos(walker).restorePos(sample);
    wf->updateLap(wfdata, sample);
    sys->calcKinetic(wfdata, sample, wf, kinetic);
    doublevar coulpot=local_energy(walker);
    nonloc=0;
    if(update_psp) 
      pseudo->calcNonlocWithTest(wfdata,sys, sample, wf,psp_test(walker),nonloc );    
    wf->getVal(wfdata, 0, wfval);
    for(int w=0; w< nfunctions; w++) {
      doublevar energy=kinetic(w) + coulpot+ nonloc(w);
      avgkin+=kinetic(w);
      avgpot+=coulpot;
      avgnon+=nonloc(w);
      avgen+=energy;
      doublevar reweight=1.0;
      if(use_weights) {
        reweight=exp(2*(wfval.amp(w,0)-orig_vals(walker).amp(w,0)));
      }
      doublevar weight=reweight*guide_wf.getWeight(wfval,wfval, w);
      switch(min_function) {
      case min_variance:
        variance(w)+=(energy-eref)*(energy-eref)*weight;
        break;
      case min_abs:
        variance(w)+=fabs(energy-eref)*weight;
        break;
      case min_lorentz:
        variance(w)+=log(1+(energy-eref)*(energy-eref)/2)*weight;
        break;
      case min_energy:
        variance(w)+=energy*weight;
        break;
      case min_mixed:
        variance(w)+=(mixing*energy+(1.0-mixing)*(energy-eref)*(energy-eref))*weight;
        break;
      default:
        error("Optimize_method::variance() : min_function has a very strange value");
      }
      //cout << "variance " << variance(w) << endl;
      weightsum(w)+=weight;
    }
  }


  doublevar weight_tot=parallel_sum(weightsum(0));
  doublevar var_tot=parallel_sum(variance(0));
  doublevar en_tot=parallel_sum(avgen);
  
  val=var_tot/weight_tot;
  
  
  if(check) { 
    if(fabs(weightsum(0)) < 1e-14)
      error("sum of weights is ", weightsum(0), " this is way too small");
    if(val > 1e14) { 
      cout << "en " << en_tot/weight_tot << " kin " << avgkin/weight_tot 
        << " pot " << avgpot/weight_tot << " nonloc " << avgnon/weight_tot << endl;
      
      error("variance is too large: ", val );
    }
  }
  //cout << "avg en " << en_tot/weight_tot  << "  " << << endl;
  
  return en_tot/weight_tot;  
}


//------------------------------------------------------------------------

doublevar Optimize_method::derivatives(int n, Array1 <double> & parms, Array1 <double> & deriv, 
                                       double & val, int check)
{
  
  int nparms=wfdata->nparms();
  Array1 <doublevar> temp_parms(nparms);
  for(int i=0; i< nparms; i++)  {
    temp_parms(i)=parms(i+1);
  }
  
  Wf_return  wfval(nfunctions,2);
  
  wfdata->setVarParms(temp_parms);
  
  Array1 <doublevar> variance(nfunctions);  variance=0;
  Array1 <doublevar> weightsum(nfunctions);  weightsum=0;
  Array1 <doublevar> kinetic(nfunctions), nonloc(nfunctions);
  deriv.Resize(nparms); deriv=0;
  psp_buff.start_again();
  
  doublevar avgen=0;
  doublevar avgkin=0;
  doublevar avgpot=0;
  doublevar avgnon=0;
  //pseudo->setDeterministic(1);
  //cout << "loop " << endl;
  last_ens.Resize(nconfig);
  for(int walker=0; walker< nconfig; walker++) {
    config_pos(walker).restorePos(sample);
    wf->updateLap(wfdata, sample);
    Array1 <doublevar> kin_deriv(nparms);
    sys->calcKinetic(wfdata, sample, wf, kinetic);
    
    Array1 <doublevar> nonloc_deriv(nparms);
    nonloc_deriv=0;
    nonloc=0.0;
    if(update_psp) { 
      if(wfdata->supports(parameter_derivatives)) { 
        pseudo->calcNonlocParmDeriv(wfdata, sys,sample, wf,
            psp_test(walker),nonloc,nonloc_deriv );    
      }
      else { 
        pseudo->calcNonlocWithTest(wfdata, sys, sample, wf,psp_test(walker),nonloc );    
      }
    }
    
    //Take the derivative of the kinetic energy by finite difference
    Array1 <doublevar> kin_tmp(nfunctions);
    doublevar del=0.0001;
    for(int p=0; p< nparms; p++) { 
      Array1<doublevar> tparms=temp_parms;
      tparms(p)+=del;
      wfdata->setVarParms(tparms);
      wf->updateLap(wfdata,sample);
      sys->calcKinetic(wfdata, sample,wf,kin_tmp);
      if(!wfdata->supports(parameter_derivatives) && update_psp) { 
        //cout << "calculating pseudo non-analytically! " << endl;
        Array1<doublevar> tnon(nfunctions); 
        pseudo->calcNonlocWithTest(wfdata, sys, sample, wf,psp_test(walker),tnon );    
        nonloc_deriv(p)=(tnon(0)-nonloc(0))/del;
      }
      //double deriv_non_chk=(tnon(0)-nonloc(0))/del;
      //cout << "p " << p << " analytic " << nonloc_deriv(p) << " numeric " << deriv_non_chk << endl;
      kin_deriv(p)=(kin_tmp(0)-kinetic(0))/del;
      wfdata->setVarParms(temp_parms);
    }
    //------end kinetic energy derivative
    
    doublevar coulpot=local_energy(walker);
    wf->getVal(wfdata, 0, wfval);
    assert(nfunctions==1); //note that we implicitly assume this throughout
    for(int w=0; w< nfunctions; w++)
    {
      doublevar energy=kinetic(w) + coulpot+ nonloc(w);
      last_ens(walker)=energy;

      avgkin+=kinetic(w);
      avgpot+=coulpot;
      avgnon+=nonloc(w);
      avgen+=energy;
      doublevar reweight=1.0;
      if(use_weights) {
        reweight=exp(2*(wfval.amp(w,0)-orig_vals(walker).amp(w,0)));
      }
      doublevar weight=reweight*guide_wf.getWeight(wfval,wfval, w);
      switch(min_function) {
        case min_variance:
          variance(w)+=(energy-eref)*(energy-eref)*weight;
          for(int p=0; p < nparms; p++) { 
            deriv(p)+=2.0*(energy-eref)*(kin_deriv(p)+nonloc_deriv(p))*weight;
          }
          break;
        case min_abs:
          variance(w)+=fabs(energy-eref)*weight;
          for(int p=0; p < nparms; p++) { 
            if(energy > eref) deriv(p)+=(kin_deriv(p)+nonloc_deriv(p))*weight;
            else deriv(p)-=(kin_deriv(p)+nonloc_deriv(p))*weight;
          }
          break;
        case min_lorentz:
          variance(w)+=log(1+(energy-eref)*(energy-eref)/2)*weight;
          for(int p=0; p< nparms; p++) { 
            deriv(p)+=weight*2.0*(energy-eref)*(kin_deriv(p)+nonloc_deriv(p))/(1+(energy-eref)*(energy-eref)/2.0);
          }
          break;
        case min_energy:
          variance(w)+=energy*weight;
          for(int p=0; p< nparms; p++) { 
            deriv(p)+=weight*(kin_deriv(p)+nonloc_deriv(p));
          }
          break;
        case min_mixed:
          variance(w)+=(mixing*energy+(1.0-mixing)*(energy-eref)*(energy-eref))*weight;
          for(int p=0; p < nparms; p++) { 
            for(int p=0; p < nparms; p++) { 
              deriv(p)+=(mixing*(kin_deriv(p)+nonloc_deriv(p))+(1.0-mixing)*2.0*(energy-eref)*(kin_deriv(p)+nonloc_deriv(p)))*weight;
            }            
          }
          break;
        default:
          error("Optimize_method::variance() : min_function has a very strange value");
      }
      //cout << "variance " << variance(w) << endl;
      weightsum(w)+=weight;
    }
  }
  
  
  doublevar weight_tot=parallel_sum(weightsum(0));
  doublevar var_tot=parallel_sum(variance(0));
  for(int p=0; p < nparms; p++) { 
    deriv(p)=parallel_sum(deriv(p))/weight_tot;
  }
  doublevar en_tot=parallel_sum(avgen);
  
  val=var_tot/weight_tot;
  
  if(check) { 
    if(fabs(weightsum(0)) < 1e-14)
      error("sum of weights is ", weightsum(0), " this is way too small");
    if(val > 1e14) {
      cout << "en " << en_tot/weight_tot << " kin " << avgkin/weight_tot 
        << " pot " << avgpot/weight_tot << " nonloc " << avgnon/weight_tot << endl;
      
      
      error("variance is too large: ", val );
    }
  }
  
  return en_tot/weight_tot;  
}

//------------------------------------------------------------------------


void Optimize_method::iteration_print(double f, double gg, double tol,  int itn, ostream & output) {
  if(output)
  {
    string indent="";
    ofstream wfoutput(wfoutputfile.c_str());
    wfoutput.precision(15);
    wfdata->writeinput(indent,wfoutput);
    wfoutput.close();
    
    if (use_extended_output){
      char strbuff[40];
      sprintf(strbuff, "%d", itn);
      string wfoutputfile2=wfoutputfile+"_"+strbuff;
      ofstream wfoutput2(wfoutputfile2.c_str());
      wfoutput2.precision(15);
      wfdata->writeinput(indent,wfoutput2);
      wfoutput2.close();
    }
    
    output << "iteration # " << itn;
    if(min_function==min_energy) 
      output << " energy= " << f;
    else if (min_function==min_mixed)
      output << " mix of energy and variance= " << f;
    else 
      output << "  dispersion= " << sqrt(fabs(f)) ;
    output << endl << "  average gradient: " << gg << " tolerance " << tol;
    output << endl;
    if(global_options::rappture) { 
      ofstream rapout("rappture.log");
      rapout << "=RAPPTURE-PROGRESS=>" << int(100.0*doublevar(itn+1)/doublevar(iterations) )
             << "  Optimizing wave function (percentage is worst-case)" << endl;
      cout << "=RAPPTURE-PROGRESS=>" << int(100.0*doublevar(itn+1)/doublevar(iterations) )
             << "  Optimizing wave function (percentage is worst-case)" << endl;
      rapout.close();
    }
    /*  This is useful for studying how the optimization goes..not sure if it should
        be a regular feature or not, though.
    int nc=last_ens.GetDim(0);
    double avg=0;
    for(int i=0; i< nc; i++) { 
      cout << "iteration_energy_disp: " << itn <<  "  " << last_ens(i) << endl;
      avg+=last_ens(i)/nc;
    }
    cout << "iteration_energy_avg: " << itn << " " << avg << endl;
    */
  }
  
}

//------------------------------------------------------------------------
//########################################################################

double Optimize_fn::func(double * _p) { 
  Array1 <doublevar> temp_parms(param_n+1);
  for(int i=0; i< param_n+1; i++) 
    temp_parms[i]=_p[i];
  double val;
  opt_method->variance(param_n,temp_parms,val);
  //cout << "val " << val << endl;
  return val;
}

double Optimize_fn::dfunc(double * _p, double * _g) {
/*
  int m=param_n+1;
  double base=func(_p);
  double step=1e-4;
  for(int i=1; i< m; i++) {
    _p[i]+=step;
    double nwfunc=func(_p);
    _p[i]-=step;
    _g[i]=(nwfunc-base)/step;
//    cout << _g[i] << "   ";
  }
  _g[0]=base;
 */
  double base;
  Array1 <double> deriv(param_n);
  Array1 <doublevar> temp_parms(param_n+1);
  for(int i=0; i< param_n+1; i++) 
    temp_parms[i]=_p[i];
  opt_method->derivatives(param_n, temp_parms,deriv,base,1);
  
  _g[0]=base;
  for(int i=0;i< param_n; i++) {
    _g[i+1]=deriv[i];
  }
  
  return base;
//  cout << endl;
}
