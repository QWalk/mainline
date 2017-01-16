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
#include "qmc_io.h"
#include "Program_options.h"
#include "ulec.h"
#include "System.h"
#include "Generate_sample.h"

void Optimize_method2::read(vector <string> words,
                           unsigned int & pos,
                           Program_options & options)
{
 
  if(!readvalue(words, pos=0, nconfig, "NCONFIG"))
    nconfig=max(2048/mpi_info.nprocs,1);
  
  if(!readvalue(words,pos=0, iterations, "ITERATIONS"))
    iterations=30;
    //Optional options

  if(readvalue(words,pos=0, eref, "EREF")){
    eref_exists=1;
  }
  else{
    eref_exists=0;
    eref=0;
  }
 
  if(!readvalue(words, pos=0, wfoutputfile, "WFOUTPUT") )
    wfoutputfile=options.runid+".wfout";
  else  use_weights=0; 
  if(haskeyword(words, pos=0, "EXTENDED_WFOUTPUT") )
    use_extended_output=1;
  else  use_extended_output=0; 
  
  if(haskeyword(words, pos=0, "DETERMINISTIC_PP") )
    dynamic_pp=0;
  else  dynamic_pp=1; 
  
  if(haskeyword(words, pos=0, "DYNAMIC_WF") )
     dynamic_wf=1;
  else  dynamic_wf=0; 


  if(!readvalue(words,pos=0, mixing, "MIXING"))
    mixing=0.95;
  if(haskeyword(words, pos=0, "USE_WEIGHTS") )
    use_weights=1;
  else  use_weights=0; 

  if(!readvalue(words,pos=0, iter_min_read, "NFIXED_ITERATIONS"))
    iter_min_read=4;
  if(!readvalue(words,pos=0, multiply, "MULTIPLICATIVE_FACTOR"))
    multiply=2;
  
  string functiontype_str;
  if(readvalue(words, pos=0, functiontype_str, "MINFUNCTION")) {
    if(caseless_eq(functiontype_str, "VARIANCE")) {
      min_function=min_variance;
      //cout <<"Are you sure you want to use weights for variance optimization?"<<endl;
    }
    else if(caseless_eq(functiontype_str,"ENERGY")) {
      if(!use_weights) { 
        single_write(cout, "Turning on USE_WEIGHTS for energy optimization \n");
        use_weights=1;
      }
      min_function=min_energy;
    }
    else if(caseless_eq(functiontype_str,"MIXED")) {
      if(!use_weights) { 
	single_write(cout, "Turning on USE_WEIGHTS for mixed optimization \n");
	use_weights=1;
      }
      min_function=min_mixed;
    }
    else {
      error("I don't know ", functiontype_str, " for MINFUNCTION.");
    }
  }
  else {
    min_function=min_variance;
  }

  string readconfig;
  if(!readvalue(words, pos=0, readconfig, "READCONFIG"))
    readconfig="";


  //--Set up variables
  sysprop=NULL;
  allocate(options.systemtext[0],  sysprop);
  debug_write(cout, "wfdata allocate\n");
  wfdata=NULL;
  allocate(options.twftext[0], sysprop, wfdata);

  analytic_wf_ders=0;
  if(wfdata->supports(parameter_derivatives))
    analytic_wf_ders=1;

  
  
  sysprop->generatePseudo(options.pseudotext, pseudo);

  //alocate wf object for each config
  //might be expensive to store all the inverse matrices
  //if dynamic_wf=1, only one wf is stored
  if(dynamic_wf){
    sample.Resize(1);
    sample(0)=NULL;
    sysprop->generateSample(sample(0));
    wf.Resize(1); 
    wf(0)=NULL;
    wfdata->generateWavefunction(wf(0));
    sample(0)->attachObserver(wf(0));
    //sysprop->notify(sample_dynamic,0);
  }
  else{
    sample.Resize(nconfig);
    wf.Resize(nconfig);
    for(int walker=0; walker<nconfig; walker++){
      sample(walker)=NULL;
      sysprop->generateSample(sample(walker));
      wf(walker)=NULL;
      wfdata->generateWavefunction(wf(walker));
      sample(walker)->attachObserver(wf(walker));
    }
    //sysprop->notify(sample_static,0);
  }

  if(!readvalue(words,pos=0, maxnparmsatonce, "MAXNPARMS_AT_ONCE")) {
    maxnparmsatonce=wfdata->nparms();
  }

  readcheck(readconfig);
  
  if(!readvalue(words,pos=0, min_nconfig_read, "START_NCONFIG")) {
    min_nconfig_read=nconfig;
  }
  else 
    if( min_nconfig_read > nconfig ){
      cout << "Number of starting configurations is higher than all configurations!"<<endl;
      cout << "setting it >nconfig "<<endl;
      min_nconfig_read=nconfig;
    }
  //set the printout level
  debug_out=0;
}

//-------------------------------------------------------------------

void Optimize_method2::readcheck(string & readconfig) {
  //This can be modified to get DMC weights; we just need to make a 
  //new structure that contains config_pos and the weights. 
  //If someone actually uses this, they may reimplement it.
  dmc_weight.Resize(nconfig);
  dmc_weight=1.0;
  
  if(readconfig ==""){
    Primary guidewf;
    generate_sample(sample(0),wf(0),wfdata,&guidewf,nconfig,config_pos);

    //error("No file name given for READCONFIG ", readconfig);
  }
  else { 
    read_configurations(readconfig, config_pos);
  }
  if(config_pos.GetDim(0) < nconfig) { 
    error("Not enough configurations in ", readconfig);
  }
  
  if(!dynamic_wf) { 
    for(int i=0; i< nconfig; i++) { 
      config_pos(i).restorePos(sample(i));
    }
  }
  else { config_pos(0).restorePos(sample(0)); } 
}

//-------------------------------------------------------------------

int Optimize_method2::showinfo(ostream & os)
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
  os << "Wavefunction optimization " << endl;
  os << "Iterations : "  << iterations << endl;
  os << "Configurations per processor : " <<  nconfig   << endl;
  os << "Number of processors : "        <<  mpi_info.nprocs << endl;
  os << "Total configurations : " <<          nconfig*mpi_info.nprocs << endl;
  os << "Number of starting configurations : " << min_nconfig_read << endl;
  os << "Number of fixed iterations : "  << iter_min_read << endl;
  os << "Maximum number of parameters done at once : " <<  maxnparmsatonce<<endl;
  os << "Multiplicative factor : "  <<  multiply<< endl;
  os << "nparms " << wfdata->nparms() << endl;
  os << "Wavefunction parameter derivatives : ";
  if(analytic_wf_ders)
    os << "analytic" << endl;
  else
    os << "finite difference" << endl;
  os << "Wavefunction storage mode : ";
  if(dynamic_wf)
    os << "one wf for all walkers" << endl;
  else
    os << "one wf for each walker" << endl;
  os << "Nonlocal PP calculation mode : ";
  if(dynamic_pp)
    os << "dynamic" << endl;
  else
    os << "static" << endl;
  if(eref_exists && !use_weights)
    os << "Reference energy : "           << eref << endl;
  else
    os << "Reference energy : energy in each iteration" << endl;
  os << "Minimization function : ";
  switch(min_function) {
  case min_variance:
    os << "Variance \n";
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

void Optimize_method2::run(Program_options & options, ostream & output)
{
  int nparms=wfdata->nparms(); //Number of variables
  Array1 <double> delta(nparms);
  Array1 <double> grad(nparms);  //Storage for the derivative
  //best delta is 0.0001 if parameters are det weights 
  delta=0.0001;
  //called also a measurement vector
  grad=0.0; //should be zero for when looking for grad=0;

  doublevar value, energy, variance;  //Best value of the function
  Array1 <doublevar> temp_parms(nparms);
  
  wfdata->getVarParms(temp_parms);
  lastparms.Resize(nparms);
  lastparms=temp_parms;
  cout.precision(16);
  
  if(nparms<= 0 ) error("There appear to be no parameters to optimize!");

  orig_vals.Resize(nconfig);
  local_energy.Resize(nconfig);
  if(dynamic_pp){
    psp_test.Resize(nconfig);
  }
  doublevar sum_tmp=0.0;
  for(int walker=0; walker < nconfig; walker++)  {
    if(walker==0)
      nfunctions=wf(0)->nfunc();
    orig_vals(walker).Resize(nfunctions, 2);
    if(dynamic_wf){
      config_pos(walker).restorePos(sample(0));
      wf(0)->updateLap(wfdata, sample(0));
      wf(0)->getVal(wfdata, 0, orig_vals(walker));
      local_energy(walker)=sysprop->calcLoc(sample(0));
    }
    else{
      wf(walker)->updateLap(wfdata, sample(walker));
      wf(walker)->getVal(wfdata, 0, orig_vals(walker));
      local_energy(walker)=sysprop->calcLoc(sample(walker));
    }
    sum_tmp+=orig_vals(walker).amp(0,0);
    if(dynamic_pp){
      //generate random numbers for PP integration
      psp_test(walker).Resize(pseudo->nTest());
      for(int i=0; i< pseudo->nTest(); i++) { 
        psp_test(walker)(i)=rng.ulec();
      }
    } 
    else{ //not dynamic_pp
      if(dynamic_wf){
        pseudo->initializeStatic(wfdata, sample(0), wf(0), psp_buff);
        //wf(0)->notify(sample_static,0);
      }
      else{
        pseudo->initializeStatic(wfdata, sample(walker), wf(walker), psp_buff);
        wf(walker)->notify(sample_static,0);
      }
    }
  }//walker

  ln_norm_orig_vals=parallel_sum(sum_tmp/nconfig)/parallel_sum(1);
  if(use_weights || !eref_exists){
    // cout << "1st function calculation to get estimation of ref energy";
    func_val(nparms, temp_parms, value, energy, variance, nconfig, output);
    eref=energy;
  }

  //starting minimization routine
  int iter=0;
  int iter_min=iter_min_read;
  int min_nconfig= min_nconfig_read;
  int divideto=int(nparms/maxnparmsatonce);
  if (output){
     cout << "Starting min procedure"<<endl;
     cout << "divinding to "<< divideto<<" parts"<<endl; 
  }
  int nparms_start=0;
  int nparms_end=maxnparmsatonce;
  //int nparms_now=0;
  while(nparms_start < nparms){
    if(nparms_end > nparms)
      nparms_end=nparms;
    if(output)
      cout << "nparms_start "<<nparms_start<<" and  nparms_end "<< nparms_end<<endl;
    //perform the optimization
    if (LEVMAR_DER(temp_parms, nparms_start, nparms_end, delta,  min_nconfig, iter_min, iterations, iter, output)!=1)
      error ("Died in Lev-Mar min routine");
    if (output)
      cout << "Iterations needed: "<<iter-1<<endl;
    
    nparms_start+=maxnparmsatonce;
    nparms_end+=maxnparmsatonce;
  }

  //final evaluation before exit
  func_val(nparms, temp_parms, value, energy, variance, nconfig, output);

  if(output)
  {
    output << "Optimization finished.  ";
    output << "New wave function is in " << wfoutputfile << endl;
    output << "Function's final value " << value<< " , energy: "
           << energy << " and variance: "<<variance << endl;
    wfdata->renormalize();
    //wfdata->showinfo(output);
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

void Optimize_method2::func_val(int n, const Array1 <double> & parms, double & val, double & energy_mean, double & variance,
                                int & min_nconfig,ostream & output)
{

  int nparms=wfdata->nparms();
  Array1 <doublevar>   temp_parms(nparms);

  // setting up parameters
  for(int i=0; i< nparms; i++){
    temp_parms(i)=parms(i);
  }
  
  Array1 < Wf_return > wfval(min_nconfig);
  wfdata->setVarParms(temp_parms);
    
  doublevar sum_tmp=0.0;
  for(int walker=0; walker< min_nconfig; walker++){
    wfval(walker).Resize(nfunctions,2);
    if(dynamic_wf){
      config_pos(walker).restorePos(sample(0)); 
      wf(0)->updateVal(wfdata, sample(0));
      wf(0)->getVal(wfdata, 0, wfval(walker));
    }
    else{
      wf(walker)->updateVal(wfdata, sample(walker));
      wf(walker)->getVal(wfdata, 0, wfval(walker));
    }
    sum_tmp+=wfval(walker).amp(0,0);
  }
  
  doublevar ln_norm_new_vals=parallel_sum(sum_tmp/min_nconfig)/parallel_sum(1);
  
  //if(output)
  //cout << "Func_val: Aver. Sum of ln(psi)/min_nconfig= "<<ln_norm_new_vals<<endl;
  
  doublevar variance_sum=0;
  doublevar energy_sum=0;
  doublevar weightsum=0;
  Array1 <doublevar> kinetic(nfunctions), nonloc(nfunctions);
  if(!dynamic_pp){
    psp_buff.start_again();
  }
  doublevar weight_max=1;
  for(int walker=0; walker< min_nconfig; walker++){
    if(dynamic_wf){
      config_pos(walker).restorePos(sample(0)); 
      wf(0)->updateLap(wfdata, sample(0));
      sysprop->calcKinetic(wfdata, sample(0), wf(0), kinetic);
      //wf(0)->getVal(wfdata, 0, wfval);
    }
    else{
      wf(walker)->updateLap(wfdata, sample(walker));
      sysprop->calcKinetic(wfdata, sample(walker), wf(walker), kinetic);
      //wf(walker)->getVal(wfdata, 0, wfval);
    }
 
    doublevar coulpot=local_energy(walker);
    if(dynamic_pp){
      if(dynamic_wf)
	pseudo->calcNonlocWithTest(wfdata, sysprop,sample(0), wf(0), psp_test(walker), nonloc );    
      else{
	pseudo->calcNonlocWithTest(wfdata, sysprop,sample(walker), wf(walker), psp_test(walker), nonloc );  
      }
    }
    else{
      if(dynamic_wf)
	pseudo->calcNonlocWithFile(wfdata,sysprop, sample(0), wf(0),
				   nonloc, psp_buff);
      else
	pseudo->calcNonlocWithFile(wfdata, sysprop, sample(walker), wf(walker),
				   nonloc, psp_buff);
    }

    if(debug_out)
      cout << "kinetic(0) "<<kinetic(0)<<" coulpot "<<coulpot<<" nonloc(0) "<<nonloc(0)<<endl;

    doublevar energy=kinetic(0) + coulpot+ nonloc(0);
    doublevar reweight=1.0;
    //get the reweighting factor
    if(use_weights) {
      reweight=exp(2*(wfval(walker).amp(0,0)-orig_vals(walker).amp(0,0)-ln_norm_new_vals+ln_norm_orig_vals));
      //reweight=2.0*reweight/(1+reweight);
    }
    doublevar weight=reweight;
    variance_sum+=(energy-eref)*(energy-eref)*weight;
    energy_sum+=energy*weight;
    weightsum+=weight;
    if(use_weights){
      //if(weight <0 )
      //  error("I have a negative weight")
      if(weight > weight_max)
        weight_max=weight;
    }
  }//walker
   
  if(use_weights){
    //doublevar average_weight=weightsum/min_nconfig;
    //cout << "Average weight= "<<average_weight<<" and maximum weight "<<weight_max<<endl;
    if(fabs(weightsum) < 1e-14)
      error("sum of weights is ", weightsum, " this is way too small");
  }
  
  doublevar weight_sum;
  weight_sum=parallel_sum(weightsum);
  variance=parallel_sum(variance_sum)/weight_sum;
  if(variance > 1e14)
    error("variance is too large: ", variance );
  energy_mean=parallel_sum(energy_sum)/weight_sum;
    switch(min_function)
    {
    case min_variance:
      val=variance;
       
      break;
    case min_energy:
      val=energy_mean;
      break;
    case min_mixed:
      val=mixing*energy_mean+(1-mixing)*variance;
      break;  
    default:
      error("Optimize_method::variance() : min_function has a very strange value");
    }
}

//-------------------------------------------------------------------------------
/*!
\bug
When optimizing with weights, we need to renormalize the function every
time; otherwise the weights blow up.
 */

void Optimize_method2::energy_grad(Array1 <double> & parms, int nparms_start, int nparms_end, Array1 <double> & grad, doublevar & energy_mean,
				   Array1 <double> & delta, int & min_nconfig, ostream & output)
{
  int nparms_full=wfdata->nparms();
  int nparms=nparms_end-nparms_start;
  Array1 <doublevar> temp_parms(nparms_full);
  Wf_return  wfval(nfunctions,2);
  Array1 <doublevar> kinetic(nfunctions), nonloc(nfunctions);
   Array1 <doublevar> weightsum(nparms);
  Array2 <doublevar> psi0(nconfig,2);
  Array1 <doublevar> e_local0(nconfig);
  doublevar e_local, psi;
  grad=0;
  weightsum=0;
  // setting up parameters
  temp_parms=parms;
 
  for (int i=0;i<=nparms;i++){
    if (i>0){
      temp_parms(i-1+nparms_start)+=delta(i-1+nparms_start);
      }
    wfdata->setVarParms(temp_parms);
    doublevar sum_tmp=0.0;
    for(int walker=0; walker< min_nconfig; walker++){
      if(dynamic_wf){
	config_pos(walker).restorePos(sample(0)); 
	wf(0)->updateVal(wfdata, sample(0));
	wf(0)->getVal(wfdata, 0, wfval);
      }
      else{
	wf(walker)->updateVal(wfdata, sample(walker));
	wf(walker)->getVal(wfdata, 0, wfval);
      }
      sum_tmp+=wfval.amp(0,0);
    }
    doublevar ln_norm_new_vals=parallel_sum(sum_tmp/min_nconfig)/parallel_sum(1);
    
    //if(output && i==0)
    //  cout << "(energy grad)Aver. Sum of ln(psi)/min_nconfig= "<<ln_norm_new_vals<<endl;
    if(!dynamic_pp)
      psp_buff.start_again();
    for(int walker=0; walker< min_nconfig; walker++){
      if(dynamic_wf){
	config_pos(walker).restorePos(sample(0)); 
	wf(0)->updateLap(wfdata, sample(0));
	sysprop->calcKinetic(wfdata, sample(0), wf(0), kinetic);
	wf(0)->getVal(wfdata, 0, wfval);
      }
      else{
	wf(walker)->updateLap(wfdata, sample(walker));
	sysprop->calcKinetic(wfdata, sample(walker), wf(walker), kinetic);	
	wf(walker)->getVal(wfdata, 0, wfval);
      }
      doublevar coulpot=local_energy(walker);

      if(dynamic_pp){
	if(dynamic_wf){
	  pseudo->calcNonlocWithTest(wfdata, sysprop, sample(0), wf(0), psp_test(walker),nonloc );
	}
	else{
	  pseudo->calcNonlocWithTest(wfdata, sysprop, sample(walker), wf(walker), psp_test(walker),nonloc );
	}
      }
      else{
	if(dynamic_wf){
	  pseudo->calcNonlocWithFile(wfdata, sysprop, sample(0), wf(0),
				     nonloc, psp_buff);
	}
	else{  
	  pseudo->calcNonlocWithFile(wfdata, sysprop,sample(walker), wf(walker),
				     nonloc, psp_buff);
	}
      }
      
      
      if (i==0) {
        e_local0(walker)=kinetic(0) + coulpot+ nonloc(0);
        //cout <<"&& "<<kinetic(0)<<" "<<coulpot<<" "<<nonloc(0)<<endl;
        psi0(walker,0)= wfval.sign(0);
        psi0(walker,1)= wfval.amp(0,0);
      }
      else {
        psi= wfval.sign(0)*psi0(walker,0)*exp(wfval.amp(0,0)-psi0(walker,1));
        //cout <<"psi/psi0(walker) "<<psi<<endl; 
        e_local=kinetic(0) + coulpot+ nonloc(0);
        doublevar reweight=1.0;
        if(use_weights) {
          reweight=exp(2*(wfval.amp(0,0)-orig_vals(walker).amp(0,0)-ln_norm_new_vals+ln_norm_orig_vals));
        }
        doublevar weight=reweight;
        grad(i-1)+=2.0*((psi-1)/delta(i-1+nparms_start))*(e_local0(walker)-eref)*weight;
        weightsum(i-1)+=weight;
      }
    }
    if (i > 0)
      temp_parms(i-1+nparms_start)-=delta(i-1+nparms_start);
  }//end of i loop 

  parallel_sum(weightsum);
  parallel_sum(grad);
  for (int i=0;i<nparms;i++){
    grad(i)/=weightsum(i);
    if(fabs(weightsum(i)) < 1e-14)
      error("sum of weights is ", weightsum(i), " this is way too small");
    if(abs(grad(i)) > 1e14)
      error("gradient is too large: ", grad(i) );
  }
}
//------------------------------------------------------------------------

void Optimize_method2::func_hessian(Array1 <double> & parms, int nparms_start, int nparms_end, Array2 <double> & hessian,
				    Array1 <doublevar>  & grad_var_tmp, doublevar & energy_mean, 
				    Array1 <double> & delta, int & min_nconfig, ostream & output)
{
  int nparms_full=wfdata->nparms();
  int nparms=nparms_end-nparms_start;
  
  Array1 <doublevar> temp_parms(nparms_full);
  Array1 <doublevar> grad_eng_tmp(nparms);

  if (min_function==min_variance || min_function==min_mixed){
    energy_grad(parms, nparms_start, nparms_end, grad_eng_tmp, energy_mean, delta, min_nconfig, output);
  }

  // setting up parameters
  temp_parms=parms;

  Wf_return  wfval(nfunctions,2);
  Array1 <doublevar> kinetic(nfunctions), nonloc(nfunctions); 
  Array2 <doublevar> hess1_eng(nparms,nparms); 
  Array2 <doublevar> hess2_eng(nparms,nparms); 
  Array2 <doublevar> hess3_1eng(nparms,nparms);
  Array2 <doublevar> hess3_2eng(nparms,nparms);
  Array2 <doublevar> hess(nparms,nparms);
  Array1 <doublevar> grad_wf(nparms);
  Array1 <doublevar> grad_e_local(nparms);
  Array1 <doublevar> grad_var(nparms);
  Array1 <doublevar> grad_eng(nparms);
  Array1 <doublevar> grad_wf_tmp(nparms);
  Array1 <doublevar> grad_e_local_tmp(nparms);
  Array1 < Array2 <doublevar> > Psi(nparms+1); 
  Array1 < Array1 <doublevar> > e_local(nparms+1);
  Array2 <doublevar> Weightsum(nparms,nparms);
  Array1 <doublevar> weightsum(nparms);
  hess1_eng=hess2_eng=hess3_1eng=hess3_2eng=hess=0.0;
  grad_wf=grad_e_local=grad_var=grad_eng=0;
  weightsum=0;
  Weightsum=0;
  
  doublevar reweight,psi,weight;;
  Array1 <doublevar> ddelta(nparms);
  Array1 < Array1 <doublevar> > e_local_gradient(nparms);
  Array1 < Array1 <doublevar> > wf_gradient(nparms);
  for(int i=0;i<=nparms;i++){
    if (i<nparms){
      ddelta(i)=delta(i+nparms_start);
      e_local_gradient(i).Resize(min_nconfig);
      wf_gradient(i).Resize(min_nconfig);
    }
    //init of Psi and e_local array
    Psi(i).Resize(nconfig,2);
    e_local(i).Resize(nconfig);
    Psi(i)=0.0;
    e_local(i)=0.0;
  }

  //starting calculation the derivatives using finite difference
  for (int i=0;i<=nparms;i++){
    if (i>0){
      temp_parms(i-1+nparms_start)+=ddelta(i-1);
    }
    for (int j=i;j<=nparms;j++){
      if (j>0){
        temp_parms(j-1+nparms_start)+=ddelta(j-1);
      }
      //set new parameters
      wfdata->setVarParms(temp_parms);
      doublevar sum_tmp=0.0;
      for(int walker=0; walker< min_nconfig; walker++){
	if(dynamic_wf){
	  config_pos(0).restorePos(sample(0)); 
	  wf(0)->updateVal(wfdata, sample(0));
	  wf(0)->getVal(wfdata, 0, wfval);
	}
	else{
	  wf(walker)->updateVal(wfdata, sample(walker));
	  wf(walker)->getVal(wfdata, 0, wfval);
	}
        sum_tmp+=wfval.amp(0,0);
      }
      doublevar ln_norm_new_vals=parallel_sum(sum_tmp/min_nconfig)/parallel_sum(1);
      //cout << "(hess) Sum of ln(psi)/min_nconfig= "<<sum_tmp/min_nconfig<<endl;
      //if(output)
      //	cout << "(hess) Aver. Sum of ln(psi)/min_nconfig= "<<ln_norm_new_vals<<endl;

      if(!dynamic_pp)
	psp_buff.start_again();
      for(int walker=0; walker< min_nconfig; walker++){
	if(dynamic_wf)
	  config_pos(walker).restorePos(sample(0)); 
	if(i==0){
	  if(dynamic_wf){
	    wf(0)->updateLap(wfdata, sample(0));
	    sysprop->calcKinetic(wfdata, sample(0), wf(0), kinetic);
	    wf(0)->getVal(wfdata, 0, wfval);
	  }
	  else{
	    wf(walker)->updateLap(wfdata, sample(walker));
	    sysprop->calcKinetic(wfdata, sample(walker), wf(walker), kinetic);
	    wf(walker)->getVal(wfdata, 0, wfval);
	  }
	  doublevar coulpot=local_energy(walker);
	  if(dynamic_pp){
	    if(dynamic_wf)
	      pseudo->calcNonlocWithTest(wfdata,sysprop, sample(0), wf(0), psp_test(walker),nonloc );
	    else
	      pseudo->calcNonlocWithTest(wfdata,sysprop, sample(walker), wf(walker), psp_test(walker),nonloc );
	  }
	  else{
	    if(dynamic_wf){
	      pseudo->calcNonlocWithFile(wfdata,sysprop, sample(0), wf(0),
					 nonloc, psp_buff);
	    }
	    else{
	      pseudo->calcNonlocWithFile(wfdata,sysprop, sample(walker), wf(walker),
					 nonloc, psp_buff);
	    }
	  }
	  e_local(j)(walker)=kinetic(0) + coulpot+ nonloc(0);

          reweight=1.0;
          if(use_weights) {
            reweight=exp(2*(wfval.amp(0,0)-orig_vals(walker).amp(0,0)-ln_norm_new_vals+ln_norm_orig_vals));
          }
          weight=reweight;
         
          
          //if j=0 we store Psi0(walker))
          if(j==0){
            Psi(0)(walker,0)=wfval.sign(0);
            Psi(0)(walker,1)=wfval.amp(0,0);
          }
          //if j>0 we calc. averages of single index quantities
          else {
            //cout <<"j>0 we calc. averages of single index quantities"<<endl;
            Psi(j)(walker,0)=wfval.sign(0)*Psi(0)(walker,0);
            Psi(j)(walker,1)=wfval.sign(0)*Psi(0)(walker,0)*exp(wfval.amp(0,0)-Psi(0)(walker,1));
            grad_wf(j-1)+=(  (Psi(j)(walker,1)-1)/ddelta(j-1)  )*weight;
            e_local_gradient(j-1)(walker)=((e_local(j)(walker)-e_local(0)(walker))/ddelta(j-1))*weight;
            grad_e_local(j-1)+=e_local_gradient(j-1)(walker);
            wf_gradient(j-1)(walker)=(Psi(j)(walker,1)-1)/ddelta(j-1);
            switch(min_function)
            {
              case min_energy:
                grad_eng(j-1)+=2.0*wf_gradient(j-1)(walker)*(e_local(0)(walker)-energy_mean)*weight;
                break;
              case min_variance:
                grad_var(j-1)+=2.0*(e_local_gradient(j-1)(walker)-grad_eng_tmp(j-1))*(e_local(0)(walker)-eref)*weight; 
                break;  
              case min_mixed:
                grad_eng(j-1)+=2.0*wf_gradient(j-1)(walker)*(e_local(0)(walker)-energy_mean)*weight;
                grad_var(j-1)+=2.0*(e_local_gradient(j-1)(walker)-grad_eng_tmp(j-1))*(e_local(0)(walker)-eref)*weight;
                break;	
              default:
                error("Optimize_method2::variance() : min_function has a very strange value");
              }
            
            weightsum(j-1)+=weight;
          }
          
        }
        // i>0 j>=i calculation of double index quantities
        else {
          if(i==1 && walker==0 ){
	    //doing averages over processor for single index quantities
            if(j==1){
              if (output)
                cout << "Gradients of: minimized quantity, vmc energy, wavefunction and local energy"<<endl;
	      switch(min_function){
	        case min_energy:
		  parallel_sum(grad_eng);
		  break;
		case min_variance:  
		  parallel_sum(grad_var);
		  break;
		case min_mixed:
		  parallel_sum(grad_eng);
		  parallel_sum(grad_var);
		  break;
		}
	      parallel_sum(weightsum);
	      parallel_sum(grad_wf);
	      parallel_sum(grad_e_local);
	    }//end of j==1
            
	    switch(min_function){
	      case min_energy:
		grad_var_tmp(j-1)=grad_eng(j-1)/weightsum(j-1);
		grad_eng_tmp(j-1)=grad_var_tmp(j-1);
		break;
	      case min_variance:  
		grad_var_tmp(j-1)=grad_var(j-1)/weightsum(j-1);
		break;
	      case min_mixed:
		grad_var_tmp(j-1)=mixing*grad_eng(j-1)/weightsum(j-1)
		  +(1-mixing)*grad_var(j-1)/weightsum(j-1);
	        break;
	    }
	    grad_wf_tmp(j-1)=grad_wf(j-1)/weightsum(j-1);
	    grad_e_local_tmp(j-1)=grad_e_local(j-1)/weightsum(j-1);
	    if (output)
	      cout <<j<<":  "<<  grad_var_tmp(j-1)<<"  "<<grad_eng_tmp(j-1)<<"  "<<grad_wf_tmp(j-1)<<"  "<<grad_e_local_tmp(j-1)<<endl;
            
	  }//i==1 && walker==0 
	  //collecting  double index quantities
	  if(dynamic_wf){
	    wf(0)->updateVal(wfdata, sample(0));
	    wf(0)->getVal(wfdata, 0, wfval);
	  }
	  else{
	    wf(walker)->updateVal(wfdata, sample(walker));
	    wf(walker)->getVal(wfdata, 0, wfval);
	  }
          reweight=1.0;
          if(use_weights) {
            reweight=exp(2*(wfval.amp(0,0)-orig_vals(walker).amp(0,0)-ln_norm_new_vals+ln_norm_orig_vals));
          }
          weight=reweight;
          doublevar hessian_wf;
	  switch(min_function)
          {
            case min_energy:
              psi=wfval.sign(0)*Psi(0)(walker,0)*exp(wfval.amp(0,0)-Psi(0)(walker,1));
              hessian_wf=(psi-Psi(i)(walker,1)-Psi(j)(walker,1)+1)/(ddelta(i-1)*ddelta(j-1));
	      hess1_eng(i-1,j-1)+=(hessian_wf- wf_gradient(i-1)(walker)* wf_gradient(j-1)(walker)) *(e_local(0)(walker)-energy_mean)*weight;
              hess2_eng(i-1,j-1)+=(wf_gradient(i-1)(walker)-grad_wf_tmp(i-1))*(wf_gradient(j-1)(walker)-grad_wf_tmp(j-1))*
                (e_local(0)(walker)-energy_mean)*weight;
              hess3_1eng(i-1,j-1)+=wf_gradient(i-1)(walker)*e_local_gradient(j-1)(walker)*weight;
              hess3_2eng(i-1,j-1)+=wf_gradient(j-1)(walker)*e_local_gradient(i-1)(walker)*weight;
              break;
            case min_variance: 
              hess(i-1,j-1)+=((e_local_gradient(i-1)(walker)-grad_eng_tmp(i-1))*
                              (e_local_gradient(j-1)(walker)-grad_eng_tmp(j-1)))*weight;
              break;
            case min_mixed:
              psi=wfval.sign(0)*Psi(0)(walker,0)*exp(wfval.amp(0,0)-Psi(0)(walker,1));
              hessian_wf=(psi-Psi(i)(walker,1)-Psi(j)(walker,1)+1)/(ddelta(i-1)*ddelta(j-1));
              hess1_eng(i-1,j-1)+=(hessian_wf- wf_gradient(i-1)(walker)* wf_gradient(j-1)(walker)) *(e_local(0)(walker)-energy_mean)*weight;
              hess2_eng(i-1,j-1)+=(wf_gradient(i-1)(walker)-grad_wf_tmp(i-1))*(wf_gradient(j-1)(walker)-grad_wf_tmp(j-1))*
                (e_local(0)(walker)-energy_mean)*weight;
              hess3_1eng(i-1,j-1)+=wf_gradient(i-1)(walker)*e_local_gradient(j-1)(walker)*weight;
              hess3_2eng(i-1,j-1)+=wf_gradient(j-1)(walker)*e_local_gradient(i-1)(walker)*weight;
              hess(i-1,j-1)+=((e_local_gradient(i-1)(walker)-grad_eng_tmp(i-1))*
                              (e_local_gradient(j-1)(walker)-grad_eng_tmp(j-1)))*weight;
              break;
          }
          Weightsum(i-1,j-1)+=weight;
          //if(walker==0) cout <<"#";
        }// i>0 j>=i
      }//end of loop over walkers
      
           
      if (j>0){
        // cout <<"["<<  temp_parms(j-1)-parms(j-1) << "] ";
        temp_parms(j-1+nparms_start)-=ddelta(j-1); 
      }
     
    }//end of loop over j
    //cout <<endl;
    if(i>0)
      temp_parms(i-1+nparms_start)-=ddelta(i-1); 
  }//end of loop over i
  //cout <<endl;

  // doing averages over processor for double index quantities
  // + calculation of final hessian.
  switch(min_function){
  case min_energy:
    parallel_sum(hess1_eng);
    parallel_sum(hess2_eng);
    parallel_sum(hess3_1eng);
    parallel_sum(hess3_2eng);
    break;
  case min_variance:
    parallel_sum(hess);
    break;
  case min_mixed:
    parallel_sum(hess1_eng);
    parallel_sum(hess2_eng);
    parallel_sum(hess3_1eng);
    parallel_sum(hess3_2eng);
    parallel_sum(hess);
    break;
  default:
    error("Optimize_method2::hessian error!");
  }
  parallel_sum(Weightsum);
  doublevar h_tmp, h_tmp2;

  for (int i=0;i<nparms;i++){
    for (int j=i;j<nparms;j++){
      switch(min_function)
        {
        case min_energy:
          h_tmp=2.0*(hess1_eng(i,j)+2.0*hess2_eng(i,j))/Weightsum(i,j);
          hessian(i,j)=h_tmp+hess3_1eng(i,j)/Weightsum(i,j) -grad_wf_tmp(i)*grad_e_local_tmp(j)
                            +hess3_2eng(i,j)/Weightsum(i,j) -grad_wf_tmp(j)*grad_e_local_tmp(i);
          hessian(j,i)=hessian(i,j);
          break;
        case min_variance:
          hessian(i,j)=2.0*hess(i,j)/Weightsum(i,j);
          hessian(j,i)=hessian(i,j);
          break;
	case min_mixed:
	  h_tmp=2.0*(hess1_eng(i,j)+2.0*hess2_eng(i,j))/Weightsum(i,j)
	    +hess3_1eng(i,j)/Weightsum(i,j) -grad_wf_tmp(i)*grad_e_local_tmp(j)
	    +hess3_2eng(i,j)/Weightsum(i,j) -grad_wf_tmp(j)*grad_e_local_tmp(i);
	  h_tmp2=2.0*hess(i,j)/Weightsum(i,j);
	  hessian(i,j)=mixing*h_tmp+(1-mixing)*h_tmp2;
          hessian(j,i)=hessian(i,j);  
	  break;
        default:
          error("Optimize_method2::hessian error!");
        }
      // if (output)
      //  cout << hessian(i,j)<<"  ";
    }
    //if (output)
    // cout <<endl;
  }
  

  //cout <<"node "<< mpi_info.node  << " End: Optimize_method2::func_hessian"<<endl;
}

//------------------------------------------------------------------------

void Optimize_method2::energy_grad_analytical(Array1 <double> & parms, int nparms_start, int nparms_end, 
					      Array1 <double> & grad,  Array1 <double> & grad_wf, doublevar & energy_mean,
					      Array1 <double> & delta, int & min_nconfig, ostream & output)
{
  //cout << "Start: Optimize_method2::func_grad"<<endl;
  int nparms_full=wfdata->nparms();
  int nparms=nparms_end-nparms_start;

  Array1 <doublevar> temp_parms(nparms_full);
  Array1 < Wf_return > wfval(min_nconfig);
  Parm_deriv_return derivatives;
  derivatives.nparms_start=nparms_start;
  derivatives.nparms_end=nparms_end;

  Array1 <doublevar> kinetic(nfunctions), nonloc(nfunctions);
  Array1 <doublevar> grad_var(nparms);
  Array1 <doublevar> grad_wf_tmp(nparms);
  doublevar weightsum;
  grad_var=0;
  grad_wf_tmp=0;
  weightsum=0;
  
   // setting up parameters
  temp_parms=parms;
  wfdata->setVarParms(temp_parms);
  doublevar sum_tmp=0.0;
  weight.Resize(min_nconfig);
  E_local.Resize(min_nconfig);

  for(int walker=0; walker< min_nconfig; walker++){
    wfval(walker).Resize(nfunctions,2);
    if(dynamic_wf){
      config_pos(walker).restorePos(sample(0)); 
      wf(0)->updateVal(wfdata, sample(0));
      wf(0)->getVal(wfdata, 0, wfval(walker));
    }
    else{
      wf(walker)->updateVal(wfdata, sample(walker));
      wf(walker)->getVal(wfdata, 0, wfval(walker));
    }
    sum_tmp+=wfval(walker).amp(0,0);
  }
  doublevar ln_norm_new_values=parallel_sum(sum_tmp/min_nconfig)/parallel_sum(1);
  //cout << "(energy grad)Sum of ln(psi)/min_nconfig= "<<sum_tmp/min_nconfig<<endl;
  if(output)
    cout << "(energy grad)Aver. Sum of ln(psi)/min_nconfig= "<<ln_norm_new_values<<endl;
  if(!dynamic_pp){
    psp_buff.start_again();
  }
  for(int walker=0; walker< min_nconfig; walker++){
    if(dynamic_wf){
      config_pos(walker).restorePos(sample(0));
      wf(0)->updateLap(wfdata, sample(0));
      sysprop->calcKinetic(wfdata, sample(0), wf(0), kinetic);
      wf(0)->getParmDeriv(wfdata, sample(0), derivatives);
    }
    else {
      wf(walker)->updateLap(wfdata, sample(walker));
      sysprop->calcKinetic(wfdata, sample(walker), wf(walker), kinetic);
      wf(walker)->getParmDeriv(wfdata, sample(walker), derivatives);
    }
    doublevar coulpot=local_energy(walker);
    if(dynamic_pp){
      if(dynamic_wf)
	pseudo->calcNonlocWithTest(wfdata,sysprop, sample(0), wf(0), psp_test(walker),nonloc );    
      else{
	pseudo->calcNonlocWithTest(wfdata, sysprop, sample(walker), wf(walker), psp_test(walker),nonloc );
      }
    }
    else { //not dynamic_pp
      if(dynamic_wf)
	pseudo->calcNonlocWithFile(wfdata,sysprop, sample(0), wf(0),nonloc, psp_buff);
      else
	pseudo->calcNonlocWithFile(wfdata, sysprop,sample(walker), wf(walker),nonloc, psp_buff);
    }
    
    E_local(walker)=kinetic(0) + coulpot+ nonloc(0);
    doublevar reweight=1.0;
    if(use_weights) {
      reweight=exp(2*(wfval(walker).amp(0,0)-orig_vals(walker).amp(0,0)-ln_norm_new_values+ln_norm_orig_vals));
    }
    weight(walker)=reweight;
    for(int i=0; i< nparms; i++){
      //be carefull how you get grad for partial number of parms. 
      grad_var(i)+=2.0*derivatives.gradient(i)*(E_local(walker)-eref)*weight(walker);
      grad_wf_tmp(i)+=derivatives.gradient(i)*weight(walker);
    }
    weightsum+=weight(walker);
  }//walker

  //if (output)
  //  cout <<"Calculated Energy gradient and WF's gradient: "<<endl;

  doublevar sum_of_weights=parallel_sum(weightsum);
  parallel_sum(grad_var);
  parallel_sum(grad_wf_tmp);
  for (int i=0;i<nparms;i++){
    grad(i)=grad_var(i)/sum_of_weights;
    grad_wf(i)=grad_wf_tmp(i)/sum_of_weights;
    if(fabs(weightsum) < 1e-14)
      error("sum of weights is ", weightsum, " this is way too small");
    if(abs(grad(i)) > 1e14)
      error("gradient is too large: ", grad(i) );
  }
  //cout << "End: Optimize_method2::func_grad"<<endl;

}
//------------------------------------------------------------------------

void Optimize_method2::func_hessian_analytical(Array1 <double> & parms, int nparms_start, int nparms_end, Array2 <double> & hessian,
                                    Array1 <doublevar>  & grad_var, doublevar & energy_mean, 
                                    Array1 <double> & delta, int & min_nconfig, ostream & output)
{
  //cout <<"node "<< mpi_info.node << " Start: Optimize_method2::func_hessian"<<endl;
  int nparms_full=wfdata->nparms();
  int nparms=nparms_end-nparms_start;
  
  Array1 <doublevar> temp_parms(nparms_full);
  Array1 <doublevar> grad_eng(nparms);
  Array1 <doublevar> grad_wf(nparms);

  //get energy and wf's gradient  
  energy_grad_analytical(parms, nparms_start, nparms_end, grad_eng, grad_wf, energy_mean, delta, min_nconfig, output);
  
  // setting up parameters
  temp_parms=parms;
  Wf_return  wfval(nfunctions,2);
  Parm_deriv_return wfders;
  wfders.nparms_start=nparms_start;
  wfders.nparms_end=nparms_end;
  
  
  if(min_function!=min_variance){
     wfders.need_hessian=1;
  }

  Array1 <doublevar> kinetic(nfunctions), nonloc(nfunctions); 
  Array2 <doublevar> hess1_eng(nparms,nparms); 
  Array2 <doublevar> hess2_eng(nparms,nparms); 
  Array2 <doublevar> hess3_1eng(nparms,nparms);
  Array2 <doublevar> hess3_2eng(nparms,nparms);
  Array2 <doublevar> hess1_var(nparms,nparms);
  Array1 <doublevar> grad_e_local(nparms);
  Array1 <doublevar> grad_var_tmp(nparms);
  Array1 <doublevar> grad_e_local_tmp(nparms);

  
  doublevar weightsum;
 
  //setting everything to zero
  hess1_eng=hess2_eng=hess3_1eng=hess3_2eng=hess1_var=0.0;
  grad_e_local_tmp=grad_var_tmp=0;
  weightsum=0;
  
  Array1 <doublevar> ddelta(nparms);
  doublevar sum_of_weights=0;
    
  //set variables
  Array1 <Array1 <doublevar> > e_local_gradient(nparms);
  

  for(int i=0;i<=nparms;i++){
    if (i<nparms){
      ddelta(i)=delta(i+nparms_start);
      e_local_gradient(i).Resize(min_nconfig);
    }
  }

  
  Array1 <doublevar> kinetic_old(min_nconfig);
  Array1 <doublevar> nonloc_old(min_nconfig);
  Array1 < Array1 <doublevar> > nonloc_deriv(min_nconfig);
  

  //get the local energy components for ach walker
  if(!dynamic_pp){
    psp_buff.start_again();
  }
  for(int walker=0; walker< min_nconfig; walker++){
    if(dynamic_wf){
      config_pos(walker).restorePos(sample(0));
      wf(0)->updateLap(wfdata, sample(0));
      sysprop->calcKinetic(wfdata, sample(0), wf(0), kinetic);
    }
    else{
      wf(walker)->updateLap(wfdata, sample(walker));
      sysprop->calcKinetic(wfdata, sample(walker), wf(walker), kinetic);
    }
    kinetic_old(walker)=kinetic(0);
    doublevar coulpot=local_energy(walker);

    if(dynamic_pp){
      nonloc_deriv(walker).Resize(nparms);
      if(dynamic_wf)
	pseudo->calcNonlocParmDeriv(wfdata, sysprop, sample(0), wf(0),
				    psp_test(walker), nonloc, nonloc_deriv(walker));
      else{
	pseudo->calcNonlocParmDeriv(wfdata, sysprop, sample(walker), wf(walker),
				    psp_test(walker), nonloc, nonloc_deriv(walker));
      }
    }
    else{ //not dynamic_pp
      if(dynamic_wf){
	pseudo->calcNonlocWithFile(wfdata, sysprop, sample(0), wf(0),nonloc, psp_buff);
      }
      else{
	pseudo->calcNonlocWithFile(wfdata,sysprop, sample(walker), wf(walker),nonloc, psp_buff);
      }
      nonloc_old(walker)=nonloc(0);
    }
    E_local(walker)=kinetic(0)+coulpot+nonloc(0);
  }//walker
  

  //getting derivatives of kinetic energy by finite difference
  for (int i=0;i<nparms;i++){
    temp_parms(i+nparms_start)+=ddelta(i);
    wfdata->setVarParms(temp_parms);
    if(!dynamic_pp){
      psp_buff.start_again();
    }
    for(int walker=0; walker< min_nconfig; walker++){
      
      if(dynamic_pp){
	if(dynamic_wf){
	  config_pos(walker).restorePos(sample(0)); 
	  wf(0)->updateLap(wfdata, sample(0));
	  sysprop->calcKinetic(wfdata, sample(0), wf(0), kinetic);
	}
	else{
	  wf(walker)->updateLap(wfdata, sample(walker));
	  sysprop->calcKinetic(wfdata, sample(walker), wf(walker), kinetic);
	}
	e_local_gradient(i)(walker)=nonloc_deriv(walker)(i)+(kinetic(0)-kinetic_old(walker))/ddelta(i);
      }
      else{//not dynamic_pp
	if(dynamic_wf){
	  config_pos(walker).restorePos(sample(0)); 
	  wf(0)->updateLap(wfdata, sample(0));
	  sysprop->calcKinetic(wfdata, sample(0), wf(0), kinetic);
	  pseudo->calcNonlocWithFile(wfdata,sysprop, sample(0), wf(0),nonloc, psp_buff);
	}
	else{
	  wf(walker)->updateLap(wfdata, sample(walker));
	  sysprop->calcKinetic(wfdata, sample(walker), wf(walker), kinetic);
	  pseudo->calcNonlocWithFile(wfdata, sysprop, sample(walker), wf(walker),nonloc, psp_buff);
	}
	e_local_gradient(i)(walker)=(kinetic(0)+nonloc(0)-kinetic_old(walker)-nonloc_old(walker))/ddelta(i);
      }
    }//walker
    temp_parms(i+nparms_start)-=ddelta(i);
  }
 
  grad_var_tmp=0;
  grad_e_local_tmp=0;
  
  for(int walker=0; walker< min_nconfig; walker++){
    for (int i=0;i<nparms;i++){
      switch(min_function)
	{
	case min_energy:	
	  grad_e_local_tmp(i)+=e_local_gradient(i)(walker)*weight(walker);
	  break;
	case min_variance:
	  grad_var_tmp(i)+=2.0*(e_local_gradient(i)(walker)-grad_eng(i))*
	                       (E_local(walker)-energy_mean)*weight(walker);
	  grad_e_local_tmp(i)+=e_local_gradient(i)(walker)*weight(walker);
	  break;
	case min_mixed:
	  grad_var_tmp(i)+=2.0*(e_local_gradient(i)(walker)-grad_eng(i))*(E_local(walker)-energy_mean)*weight(walker);
	  grad_e_local_tmp(i)+=e_local_gradient(i)(walker)*weight(walker);
	  break;
	}
    }
    weightsum+=weight(walker); 
  }

  sum_of_weights=parallel_sum(weightsum);
  parallel_sum(grad_var_tmp);
  parallel_sum(grad_e_local_tmp);

  if (output)
    cout << "Gradients of: minimized quantity, vmc energy, wavefunction and local energy"<<endl;
  for (int i=0;i<nparms;i++){
    switch(min_function)
      {
      case min_energy:
	grad_var(i)=grad_eng(i);
	break;
      case min_variance:  
	grad_var(i)=grad_var_tmp(i)/sum_of_weights;
	break;
      case min_mixed:
	grad_var(i)=mixing*grad_eng(i)
	  +(1-mixing)*grad_var_tmp(i)/sum_of_weights;
	break;
      }
    grad_e_local(i)=grad_e_local_tmp(i)/sum_of_weights;
    if (output)
      cout <<i+1<<":  "<<  grad_var(i)<<"   "<<grad_eng(i)<<"  "<<grad_wf(i)<<"  "<<grad_e_local(i)<<endl;
  }

  wfdata->setVarParms(temp_parms);
  for(int walker=0; walker< min_nconfig; walker++){
    if(dynamic_wf){
      config_pos(walker).restorePos(sample(0)); 
      wf(0)->updateVal(wfdata, sample(0));
      wf(0)->getParmDeriv(wfdata, sample(0),  wfders);
    }
    else{
      wf(walker)->updateVal(wfdata, sample(walker));
      wf(walker)->getParmDeriv(wfdata, sample(walker),  wfders);
    }
    for (int i=0;i<nparms;i++){
      for (int j=i;j<nparms;j++){
	switch(min_function)
	  {
	  case min_energy:
	    hess1_eng(i,j)+=(wfders.hessian(i,j)
			     -wfders.gradient(i)*wfders.gradient(j))*
	      (E_local(walker)-energy_mean)*weight(walker);
	    hess2_eng(i,j)+=(wfders.gradient(i)-grad_wf(i))*
	      (wfders.gradient(j)-grad_wf(j))*(E_local(walker)-energy_mean)*weight(walker);
	    hess3_1eng(i,j)+=wfders.gradient(i)*e_local_gradient(j)(walker)*weight(walker);
	    hess3_2eng(i,j)+=wfders.gradient(j)*e_local_gradient(i)(walker)*weight(walker);	
	    break;
	  case min_variance:   
	    hess1_var(i,j)+=(e_local_gradient(i)(walker)-grad_eng(i))*
		                  (e_local_gradient(j)(walker)-grad_eng(j))*weight(walker);
	    break;
	  case min_mixed:
	    hess1_eng(i,j)+=(wfders.hessian(i,j)
			     -wfders.gradient(i)*wfders.gradient(j))*
	      (E_local(walker)-energy_mean)*weight(walker);
	    hess2_eng(i,j)+=(wfders.gradient(i)-grad_wf(i))*
	      (wfders.gradient(j)-grad_wf(j))*(E_local(walker)-energy_mean)*weight(walker);
	    hess3_1eng(i,j)+=wfders.gradient(i)*e_local_gradient(j)(walker)*weight(walker);
	    hess3_2eng(i,j)+=wfders.gradient(j)*e_local_gradient(i)(walker)*weight(walker);
	    hess1_var(i,j)+=(e_local_gradient(i)(walker)-grad_eng(i))*
		                  (e_local_gradient(j)(walker)-grad_eng(j))*weight(walker);
	    break;
	  }
      }//end of loop over j
    }//end of loop over i
  }//end of loop over walkers
  


  // doing averages over processor for double index quantities
  // + calculation of final hessian.
  //if (output)
  // cout << "Our Hessian matrix"<<endl;
  switch(min_function){
    case min_energy:
      parallel_sum(hess1_eng); 
      parallel_sum(hess2_eng);
      parallel_sum(hess3_1eng);
      parallel_sum(hess3_2eng);
      break;
    case min_variance:
      parallel_sum(hess1_var);
      break;
    case min_mixed:
      parallel_sum(hess1_eng); 
      parallel_sum(hess2_eng);
      parallel_sum(hess3_1eng);
      parallel_sum(hess3_2eng);
      parallel_sum(hess1_var);
      break;
    default:
      error("Optimize_method2::hessian error!");
    }

  doublevar h_tmp, h_tmp2;
  for (int i=0;i<nparms;i++){
    for (int j=i;j<nparms;j++){
      switch(min_function)
        {
        case min_energy:
          h_tmp=2.0*(hess1_eng(i,j)+2.0*hess2_eng(i,j))/sum_of_weights;
          hessian(i,j)=h_tmp+hess3_1eng(i,j)/sum_of_weights-grad_wf(i)*grad_e_local(j)
                            +hess3_2eng(i,j)/sum_of_weights-grad_wf(j)*grad_e_local(i);
          hessian(j,i)=hessian(i,j);
          break;
        case min_variance:
          hessian(i,j)=2.0*hess1_var(i,j)/sum_of_weights;
          hessian(j,i)=hessian(i,j);
          break;
	case min_mixed:
          h_tmp=2.0*(hess1_eng(i,j)+2.0*hess2_eng(i,j))/sum_of_weights+
               +hess3_1eng(i,j)/sum_of_weights-grad_wf(i)*grad_e_local(j)
               +hess3_2eng(i,j)/sum_of_weights-grad_wf(j)*grad_e_local(i);
	  h_tmp2=2.0*hess1_var(i,j)/sum_of_weights;
	  hessian(i,j)=mixing*h_tmp+(1-mixing)*h_tmp2;
          hessian(j,i)=hessian(i,j);  
	  break;
        default:
          error("Optimize_method2::hessian error!");
        }
      // if (output)
      //cout << hessian(i,j)<<"  ";
    }
    //if (output)
    //cout <<endl;
  }
  
  //cout <<"node "<< mpi_info.node  << " End: Optimize_method2::func_hessian"<<endl;
}


