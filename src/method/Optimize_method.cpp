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
#include "System.h"
void Optimize_method::read(vector <string> words,
                           unsigned int & pos,
                           Program_options & options)
{

  pos=0;
  if(!readvalue(words, pos=0, nconfig, "NCONFIG"))
  {
    error("Need NCONFIG in METHOD section");
  }
  pos=0;
 
  pos=0;
  if(!readvalue(words,pos, eref, "EREF"))
    guess_eref=1;
  else guess_eref=0;
  

  //Optional options

  if(!readvalue(words,pos=0, iterations, "ITERATIONS")) 
    iterations=30;
  

  pos=0;
  if( ! readvalue(words, pos=0, pseudostore, "PSEUDOTEMP") )
    pseudostore=options.runid+".pseudo";
  canonical_filename(pseudostore);

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
  if(readvalue(words, pos, functiontype_str, "MINFUNCTION"))
  {
    if(functiontype_str== "VARIANCE")
    {
      min_function=min_variance;
    }
    else if(functiontype_str=="ABSOLUTE")
    {
      min_function=min_abs;
    }
    else if(functiontype_str=="LORENTZ")
    {
      min_function=min_lorentz;
    }
    else if(functiontype_str=="ENERGY") {
      min_function=min_energy;
    }
    else if(functiontype_str=="MIXED") {
      min_function=min_mixed;
    }
    else
    {
      error("I don't know ", functiontype_str, " for MINFUNCTION.");
    }
  }
  else
  {
    min_function=min_variance;
  }

  string readconfig;
  if(!readvalue(words, pos=0, readconfig, "READCONFIG"))
    error("READCONFIG required for OPTIMIZE method!");
  canonical_filename(readconfig);

  string oldreadconfig;
  if(readvalue(words, pos=0, oldreadconfig, "OLDREADCONFIG"))
    canonical_filename(oldreadconfig);
  

  //--Set up variables
  sysprop=NULL;
  allocate(options.systemtext[0],  sysprop);
  sysprop->notify(sample_static,0);
  sysprop->generatePseudo(options.pseudotext, pseudo);

  electrons.Resize(nconfig);
  electrons=NULL;
  wf.Resize(nconfig);
  wf=NULL;


  for(int i=0; i< nconfig; i++)
  {
    electrons(i)=NULL;
    sysprop->generateSample(electrons(i));
  }

  
  int configsread=0;
  if(readconfig !="") {
    ifstream checkfile(readconfig.c_str());
    if(!checkfile) 
      error("Couldn't open ", readconfig);
    
    long int is1, is2;
    string dummy;
    checkfile >> dummy;
    if(dummy != "RANDNUM") error("Expected RANDNUM in checkfile");
    checkfile >> is1 >> is2;
    
    
    while(checkfile >>dummy && configsread < nconfig) {
      if(read_config(dummy, checkfile, electrons(configsread)))
        configsread++;
    }
    checkfile.close();    
  }

  if(readconfig =="") {
    ifstream configin(oldreadconfig.c_str());
    if(!configin)
      error("Couldn't open ", oldreadconfig);
    configsread=read_array(electrons, configin);
    configin.close();
  }

  if(configsread < nconfig)
  {
    nconfig=configsread;
    cout << "processor " << mpi_info.node << " : "
    << "WARNING: Didn't find enough configurations in the "
    << "file.  Running optimization with only " << nconfig
    << " sample points." << endl;
  }


  debug_write(cout, "wfdata allocate\n");
  wfdata=NULL;
  allocate(options.twftext[0], sysprop, wfdata);

  if(wfdata->nparms() <= 0 ) 
    error("There appear to be no parameters to optimize!");

}

int Optimize_method::showinfo(ostream & os)
{
  os << "     System " << endl;
  sysprop->showinfo(os);
  os << endl << endl;
  os << "     Wavefunction " << endl;
  wfdata->showinfo(os);
  os << endl << endl;
  pseudo->showinfo(os);
  os << endl << endl;
  os << "-----------------------------" << endl;
  os << "Wave function optimization: version $Date: 2006/12/07 23:31:45 $ " << endl;
  
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
  
  //This was stuff for va10a, but we should probably remove it.
  //Array1 <double> g(nparms+1);  //Storage for the derivative
  //Array1 <double> h(nparms*(nparms+1)/2+1); //Storage for the Hessian
  //Array1 <double> a(3*nparms+1); //working space for va10a
  //double dfn=-0.5; //estimate of likely reduction in f(x)(order of mag)
  //dfn >0 : likely reduction
  //==0 : minimum value of F has been set in f
  // < 0 : a multiple |dfn| of the mod of init function
  //       will be taken as estimate of likely reduction
  //Array1 <double> xm(nparms+1); //indication of the magnitude of x
  //double eps=0.00001;  //Accuracy required in x(i)=eps*xm(i)
  //int mode=1; //type of storage in H.  makes va10a make the estimate
  //int maxfn=iterations; //Maximum number of function calls
  //int iprint=1; //Printing occurs every iprint iterations
  //int iexit=-1;//reason for exiting:
  // 0: estimate of hessian isn't positive definite
  // 1: Normal exit
  // 2: error because of rounding or eps is too small
  //    or truncation error is too big
  // 3: funct has been called maxfn times

  //if(min_function==min_energy || min_function==min_mixed ) {
  //  dfn=-.01;
  //  f=1.1*eref;  //I don't think this is important(it gets overwritten later)
  //}
    

  Array1 <doublevar> temp_parms(nparms);
  wfdata->getVarParms(temp_parms);
  lastparms.Resize(nparms);
  lastparms=temp_parms;
  cout.precision(16);

  

  //cout << "pseudopotential " << endl;

  FILE * pseudoout;
  pseudoout=fopen(pseudostore.c_str(), "w");
  if(!pseudoout) {
    error("couldn't open pseudopotential temporary file ", pseudostore,
          " for writing.");
  }
  

  //Here we save the original values (in case we do weighted correlated sampling)
  // and save the values for pseudopotential evaluation.
  orig_vals.Resize(nconfig);
  for(int walker=0; walker < nconfig; walker++)  {
    wfdata->generateWavefunction(wf(walker));
    nfunctions=wf(walker)->nfunc();
    orig_vals(walker).Resize(nfunctions, 2);

    electrons(walker)->attachObserver(wf(walker));
    wf(walker)->updateLap(wfdata, electrons(walker));

    wf(walker)->getVal(wfdata, 0, orig_vals(walker));

    pseudo->initializeStatic(wfdata, electrons(walker), wf(walker), pseudoout);
    wf(walker)->notify(sample_static,0);
  }

  fclose(pseudoout);
  nfunctions=wf(0)->nfunc();
  
  local_energy.Resize(nconfig);
  for(int i=0; i< nconfig; i++) 
    local_energy(i)=sysprop->calcLoc(electrons(i));
  
   //More va10a stuff
  for(int i=0; i< nparms; i++) {
    x(i+1)=temp_parms(i);
  //  g(i+1)=0;
  //  xm(i+1)=1.0;
  }
  
  //cout << "starting minimization: "<<endl;
  
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
  optimizer.macoptII(x.v,nparms);
  
  
  //va10a(nparms,x,f,g,dfn,xm,eps,mode,maxfn,iprint,iexit, output);

  variance(nparms, x,f);

  remove(pseudostore.c_str());
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
  MPI_Barrier(MPI_COMM_WORLD);
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

  //single_write(cout ,"parms\n");
  for(int i=0; i< nparms; i++)
  {
    temp_parms(i)=parms(i+1);
    
    //single_write(cout,i+1, "   " ,temp_parms(i),"\n");
    //if(temp_parms(i) != lastparms(i)) {
    // cout << "parm " << i << " changed " << lastparms(i)
    //     << " to " << temp_parms(i) << endl;
    //  lastparms(i)=temp_parms(i);
    //}

  }

  Wf_return  wfval(nfunctions,2);

  wfdata->setVarParms(temp_parms);

  Array1 <doublevar> variance(nfunctions);
  variance=0;
  Array1 <doublevar> weightsum(nfunctions);
  weightsum=0;
  Array1 <doublevar> kinetic(nfunctions), nonloc(nfunctions);


  //Pseudopotential file
  FILE * pseudoin;
  pseudoin=fopen(pseudostore.c_str(), "r");

  doublevar avgen=0;
  doublevar avgkin=0;
  doublevar avgpot=0;
  doublevar avgnon=0;
  //cout << "loop " << endl;
  for(int walker=0; walker< nconfig; walker++)
  {
    //cout << "updateLap " << endl;
    wf(walker)->updateLap(wfdata, electrons(walker));
    //cout << "calcKinetic " << endl;
    sysprop->calcKinetic(wfdata, electrons(walker), wf(walker), kinetic);
    //cout << "coulpot " << endl;
    doublevar coulpot=local_energy(walker);
    //cout << "pseudopotential " << endl;
    pseudo->calcNonlocWithFile(wfdata, electrons(walker), wf(walker),
                               nonloc, pseudoin);
    //cout << "getVal " << endl;
    wf(walker)->getVal(wfdata, 0, wfval);
    //cout << "done calc " << endl;
    for(int w=0; w< nfunctions; w++)
    {

      doublevar energy=kinetic(w) + coulpot+ nonloc(w);
      //single_write(cout ,"walker " ,walker ,"   value " ,wfval(w,0));
      //single_write(cout, "value   " ,wfval.amp(w,0) ,"\n");
      //single_write(cout, "kinetic ", kinetic(0), "\n");
      //single_write(cout, "coulomb ", coulpot, "\n");
      //single_write(cout, "nonloc ", nonloc(0), "\n");
      //single_write(cout, "energy " , energy, "\n");
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

  //cout << "Average energy " << avgen/nconfig 
  //     << " kinetic " << avgkin/nconfig
  //     << " nonloc  " << avgnon/nconfig
  //     << " potential " << avgpot/nconfig
  //     << endl;

  fclose(pseudoin);



  doublevar weight_tot=parallel_sum(weightsum(0));
  doublevar var_tot=parallel_sum(variance(0));
  doublevar en_tot=parallel_sum(avgen);
  
  val=var_tot/weight_tot;
  
  if(check) { 
    if(fabs(weightsum(0)) < 1e-14)
      error("sum of weights is ", weightsum(0), " this is way too small");
    if(val > 1e14)
      error("variance is too large: ", val );
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
  return base;
//  cout << endl;
}
