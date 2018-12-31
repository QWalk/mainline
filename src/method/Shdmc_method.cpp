/*
Copyright (C) 2009 Michal Bajdich and Fernando A. Reboredo

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

#include "Shdmc_method.h"
#include "qmc_io.h"
#include "Program_options.h"
#include "System.h"
#include "Array.h"
#include "ulec.h"
#include "MatrixAlgebra.h"

void Shdmc_method::read(vector <string> words,
                           unsigned int & pos,
                           Program_options & options)
{

  pos=0;
  if(!readvalue(words, pos=0, nconfig, "NCONFIG"))
  {
    error("Need NCONFIG in METHOD section");
  }
  pos=0;
  if(!readvalue(words,pos, iterations, "ITERATIONS"))
  {
    error("Need ITERATIONS in METHOD section");
  }

  //Optional options
  if(! readvalue(words, pos=0, wfoutputfile, "WFOUTPUT") )
  {
    wfoutputfile=options.runid+".wfout";
  }
  else { use_weights=0; }
  
  if(haskeyword(words, pos=0, "NO_EXTENDED_WFOUTPUT") )
  {
    use_extended_output=0;
  }
  else { use_extended_output=1; }
  
  if(!readvalue(words, pos=0, block_increase_factor, "BLOCK_MULTIPLY") )
  {
    block_increase_factor=1.5;
  }
  
  if(!readvalue(words, pos=0, scale_step, "SCALE_DLAMBDA") )
  {
    scale_step=1.8;
  }
  

  if(!(block_increase_factor>0))
    error("BLOCK_MULTIPLY should be larger than 0");


  if(!readvalue(words, pos=0, nblocks_max, "MAX_NBLOCK") )
  {
    nblocks_max=100;
  }


  if(!readvalue(words,pos=0, mixing, "MIXING"))
  {
    mixing=0.95;
  }
  if(haskeyword(words, pos=0, "USE_WEIGHTS") )
  {
    use_weights=1;
  }
  else { use_weights=0; } 
  
 
  string functiontype_str;
  pos=0;
  if(readvalue(words, pos, functiontype_str, "MINFUNCTION"))
  {
    if(functiontype_str== "VARIANCE")
    {
      min_function=min_variance;
      //cout <<"Are you sure you want to use weights for variance optimization?"<<endl;
    }
    else if(functiontype_str=="ENERGY") {
     
      if(!use_weights) { 
	single_write(cout, "Turning on USE_WEIGHTS for energy optimization\n");
	use_weights=1;
      }
      min_function=min_energy;
    }
    else if(functiontype_str=="MIXED") {
      if(!use_weights) { 
	single_write(cout, "Turning on USE_WEIGHTS for mixed optimization\n");
	use_weights=1;
      }
      min_function=min_mixed;
    }
    else if(functiontype_str=="WEIGHT_VARIANCE")
    {
      min_function=min_weight_variance;
    }
    else
    {
      error("I don't know ", functiontype_str, " for MINFUNCTION.");
    }
  }
  else
  {
    min_function=min_mixed;
  }

  if(!readvalue(words, pos=0, readconfig, "READCONFIG"))
    error("Need READCONFIG file in method section");
  //readconfig_non_cannonical=readconfig;
  //canonical_filename(readconfig);
  
  if(!readvalue(words, pos=0, storeconfig, "STORECONFIG"))
    storeconfig=readconfig;
  //storeconfig_non_cannonical=storeconfig;
  //canonical_filename(storeconfig);

  //read input stuff for MC run inside optimization
  pos=0;
  vector <string> mc_words_tmp;
  while(readsection(words, pos, mc_words_tmp, "MC"))
    mc_words.push_back(mc_words_tmp);
  
  if(!mc_words.size())
    error("Need MC section inside SHDMC method");

  /*
    old way of doing things
  if (!readsection(words, pos=0, mc_words, "MC")){
      mc_words.push_back("VMC");
      mc_words.push_back("NBLOCK");
      mc_words.push_back("10");
      mc_words.push_back("NSTEP");
      mc_words.push_back("10");
      mc_words.push_back("TIMESTEP");
      mc_words.push_back("1.0");
      single_write(cout,"Using default version of MC"); 
      //error("Need MC section inside NEWTON_OPT method");
  }
  */

  //--Set up variables
  sysprop=NULL;
  allocate(options.systemtext[0],  sysprop);
  sysprop->generatePseudo(options.pseudotext, pseudo);
  sample=NULL;
  wf=NULL;
  sysprop->generateSample(sample);
  wfdata=NULL;
  debug_write(cout, "wfdata allocate\n");
  allocate(options.twftext[0], sysprop, wfdata);
  wfdata->generateWavefunction(wf);
  sample->attachObserver(wf);
  //--Read configurations
  readcheck(readconfig);
}
/*
void Shdmc_method::readcheck(string & readconfig) {
  int configsread=0;
  int nwread=0;
  config_pos.Resize(nconfig);
  dmc_weight.Resize(nconfig);
  dmc_weight=1.0;
  if(readconfig !="") {
    ifstream checkfile(readconfig.c_str());
    if(!checkfile) 
      error("Couldn't open ", readconfig);
    
    long int is1, is2;
    string dummy;
    checkfile >> dummy;
    if(dummy != "RANDNUM") error("Expected RANDNUM in checkfile");
    checkfile >> is1 >> is2;

    while(checkfile >>dummy && configsread < nconfig && nwread < nconfig) {
      if(read_config(dummy, checkfile, sample)) {
        config_pos(configsread++).savePos(sample);
      }
      if(dummy=="DMC") {
	checkfile >> dummy;
	if(dummy != "{") error("Need a { after DMC");
	checkfile >> dummy >> dmc_weight(nwread);
	if(dummy != "DMCWEIGHT") {
	  error("expected DMCWEIGHT, got ", dummy);
	}
	int nwf_temp;
	checkfile >> dummy >> nwf_temp;
	if(dummy != "VALEN") {
	  error("expected VALEN, got ", dummy);
	}
	//if(nwf_temp != nwf) {
	// error("Wrong number of wavefunctions in the checkpoint file");
	nwread++;
      }
    }
    checkfile.close();    
  }

  if(readconfig ==""){
    error("No file name given for READCONFIG ", readconfig);
  }

  if(configsread < nconfig)
  {
    nconfig=configsread;
    cout << "processor " << mpi_info.node << " : "
    << "WARNING: Didn't find enough configurations in the "
    << "file.  Running optimization with only " << nconfig
    << " sample points." << endl;
  }
}
*/

void Shdmc_method::readcheck(string & filename) {
  config_pos.Resize(0);
  if(filename!="") { 
    read_configurations(filename, config_pos);
  }
  if(config_pos.GetDim(0) < nconfig) { 
    Array1 <Config_save_point> tmpconfig=config_pos;
    config_pos.Resize(nconfig);
    for(int i=0; i< tmpconfig.GetDim(0); i++) { config_pos(i)=tmpconfig(i);} 
    for(int i=tmpconfig.GetDim(0); i< nconfig; i++) {
      sample->randomGuess();
      config_pos(i).savePos(sample);
    }
  } 
}


void write_wf(Array1 <doublevar> & parms, string & wfoutputfile, Wavefunction_data * wfdata){
  if( mpi_info.node==0 ){
    string indent="";
    wfdata->setVarParms(parms);
    ofstream wfoutput(wfoutputfile.c_str());
    wfoutput.precision(15);
    wfdata->writeinput(indent,wfoutput);
    wfoutput.close();
  }
}

void Shdmc_method::wf_printout(Array1 <doublevar> & parms, int iter, 
			       doublevar & value, 
			       doublevar & energy, 
			       doublevar & energy_err, 
			       doublevar & variance, 
			       doublevar & weight_variance, 
			       ostream & output){
  wfdata->setVarParms(parms);
  int field=output.precision()+2;
  if( mpi_info.node==0 ){
    if (use_extended_output){
      output <<flush;
      string indent="";
      char strbuff[40];
      sprintf(strbuff, "%d", iter);
      string wfoutputfile2=wfoutputfile+"_"+strbuff;
      ofstream wfoutput2(wfoutputfile2.c_str());
      wfoutput2.precision(15);
      wfdata->writeinput(indent,wfoutput2);
      wfoutput2.close();
    }
    output << "iteration # " << setw(3)<<iter;
    switch(min_function)
      {
      case min_variance:
	output << setw(field)<<" dispersion= " <<setw(field)<<setiosflags(ios::fixed | ios::showpoint)<<sqrt(fabs(value)) 
	       << " energy= " <<setw(field)<<setiosflags(ios::fixed | ios::showpoint)<< energy
	       <<" +/- "<< setw(field)<< setiosflags(ios::fixed | ios::showpoint) <<energy_err
	       <<" weight dispersion= " <<setw(field)<<setiosflags(ios::fixed | ios::showpoint)<<sqrt(fabs(weight_variance))
	       << endl;
	break;
      case min_energy:
	output <<" energy= " <<setw(field)<<setiosflags(ios::fixed | ios::showpoint)<< value <<" +/- "<<setw(field)<<setiosflags(ios::fixed | ios::showpoint)
	       << energy_err <<" dispersion= "<< setw(field)<< sqrt(variance) 
	       <<" weight dispersion= " <<setw(field)<<setiosflags(ios::fixed | ios::showpoint)<<sqrt(fabs(weight_variance))
	       <<endl;
	break;
      case min_mixed:
	output <<" mixing*energy+ (1-mixing)*variance= " <<setw(field)<<setiosflags(ios::fixed | ios::showpoint)<< value 
	       <<" energy= " <<setw(field)<<setiosflags(ios::fixed | ios::showpoint)<< energy 
	       <<" +/- "<<setw(field)<<setiosflags(ios::fixed | ios::showpoint)<< energy_err
	       <<" dispersion= "<<setw(field)<<setiosflags(ios::fixed | ios::showpoint)<< sqrt(variance)
	       <<" weight dispersion= " <<setw(field)<<setiosflags(ios::fixed | ios::showpoint)<<sqrt(fabs(weight_variance))
	       <<endl;
	break;  
      case min_weight_variance:
	output <<setw(field)<<" weight dispersion= " <<setw(field)<<setiosflags(ios::fixed | ios::showpoint)<<sqrt(fabs(value)) 
			  << " energy= " <<setw(field)<<setiosflags(ios::fixed | ios::showpoint)<< energy
			  <<" +/- "<< setw(field)<< setiosflags(ios::fixed | ios::showpoint) <<energy_err
	                  <<" dispersion= "<<setw(field)<<setiosflags(ios::fixed | ios::showpoint)<< sqrt(variance)
			  << endl;
	break;
      default:
	error("Optimize_method::variance() : min_function has a very strange value");
      }
  }
}


void Shdmc_method::get_averages(Array1 <doublevar> & parms,
				int & iter,
				ostream & output, 
				Program_options & options,
				doublevar & function,
				doublevar & energy,
				doublevar & energy_err,
				doublevar & variance,
				doublevar & weight_variance,
				Array1 <doublevar> & new_parms
				){
  char strbuff[40];
  sprintf(strbuff,  "%d", nconfig);
  int ndim=parms.GetSize();

  for (int mc_size=0; mc_size<mc_words.size(); mc_size++){
    wfdata->setVarParms(parms);
    string log_label=mc_words[mc_size][0];
    char strbuff2[40];
    sprintf(strbuff2,  "%d", iter);
    log_label+=strbuff2;
    vector <string>  run_words;
    run_words.push_back(mc_words[mc_size][0]);
    run_words.push_back("NCONFIG");
    run_words.push_back(strbuff);
    run_words.push_back("READCONFIG");
    if(iter<1)
      run_words.push_back(readconfig);
    else
      run_words.push_back(storeconfig);
    run_words.push_back("STORECONFIG");
    run_words.push_back(storeconfig);

    
    for(int s=1;s<mc_words[mc_size].size();s++)
      run_words.push_back(mc_words[mc_size][s]);
    if(mpi_info.node==0 ){
      for(int s=0;s<run_words.size();s++)
	cout << run_words[s]<<" ";
      cout <<endl;
    }
  
    avg_method=NULL;
    allocate(run_words, options, avg_method);
    string logfile=options.runid+".log";
    myprop.setLog(logfile, log_label);
    ofstream logout;
    if(mpi_info.node==0 ) {
      logout.open(logfile.c_str(), ios::app);
      logout << "#-------------------------------------------------\n";
       logout << "#"<<mc_words[mc_size][0]<<" run:" <<endl;
    }
    avg_method->runWithVariables(myprop, 
				 sysprop,
				 wfdata,
				 pseudo,
				 output);
    if(mpi_info.node==0 ) {
      logout.close();
    }


    Properties_final_average  finavg;
    myprop.getFinal(finavg);
    energy=finavg.avg(Properties_types::total_energy,0);
    variance=finavg.avgvar(Properties_types::total_energy,0);
    energy_err=sqrt(finavg.err(Properties_types::total_energy,0));
    weight_variance=finavg.avgvar(Properties_types::weight,0);
    
    int navg_vals=finavg.avgavg.GetDim(0);
    for(int i=0; i< navg_vals; i++) { 
      //output << "      "<<  "average_generator { " << finavg2.avgavg(i).type << " ";
      if(finavg.avgavg(0,i).type=="linear_der"){
	int ndimm=finavg.avgavg(0,i).vals.GetDim(0)-1;
	assert(ndim==ndimm);
	for(int j=0; j< ndim; j++) { 
	  doublevar value=finavg.avgavg(0,i).vals(j)/finavg.avgavg(0,i).vals(ndim);
	  doublevar error=finavg.avgerr(0,i).vals(j)/finavg.avgavg(0,i).vals(ndim);
	  //output << value << " +/- "<<error<<endl;
	    if(abs(value)>2*error)
	      new_parms(j)=value;
	    else
	      new_parms(j)=0;
	}
	  //output << " } " << endl;
      }//end if 
      else if(finavg.avgavg(0,i).type=="linear_delta_der"){
	int ndimm=finavg.avgavg(0,i).vals.GetDim(0)-1;
	assert(ndim==ndimm);
	for(int j=0; j< ndim; j++) { 
	  doublevar value=parms(j)+scale_step*finavg.avgavg(0,i).vals(j)/finavg.avgavg(0,i).vals(ndim);
	  if(fabs(parms(j))>0)
	    delta_parms(j)=scale_step*finavg.avgavg(0,i).vals(j)/finavg.avgavg(0,i).vals(ndim); 
	  else
	    delta_parms(j)=0.0; 
	  doublevar error=finavg.avgerr(0,i).vals(j)/finavg.avgavg(0,i).vals(ndim);
	  //output << value << " +/- "<<error<<endl;
	  if(abs(value)>2*error){
	    new_parms(j)=value;
	  }
	  else
	    new_parms(j)=0;
	}
      }//end if 
    }//i
    wfdata->attachObserver(wf);
  }//mc_size
  
  switch(min_function)
    {
    case min_energy:
      function=energy;
      break;
    case min_variance:
      function=variance;
      break;
    case min_mixed:
      function=mixing*energy+(1-mixing)*variance;
      break;
    case min_weight_variance:
      function=weight_variance;
      break;
    }
  
}


int Shdmc_method::showinfo(ostream & os)
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
  os << "SHDMC wave function optimization" << endl;
  os << "Configurations per processor: " <<  nconfig   << endl;
  os << "Number of processors: "        <<  mpi_info.nprocs << endl;
  os << "Total configurations: " <<          nconfig*mpi_info.nprocs << endl;
  os << "Iterations : "  << iterations << endl;
  os << "Number of parameters " << wfdata->nparms() << endl;
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
  case min_weight_variance:
    os << "Weight variance\n";
    break;
  default:
    os << "Unknown--add to showinfo()\n";
  }
  if(use_weights) 
    os << "Reweighting using correlated sampling " << endl;
  os << "----------------------------" << endl;
  return 1;
}


void Shdmc_method::run(Program_options & options, ostream & output)
{
  output.precision(10);
  nparms=wfdata->nparms(); //Number of variables

  if(nparms<= 0 ) error("There appear to be no parameters to optimize!");
  delta_parms.Resize(nparms);
  delta_parms=0.0; 
  old_delta_parms.Resize(nparms);
  old_delta_parms=0.0;

  Array1 <doublevar> parms(nparms);
  Array1 <doublevar> new_parms(nparms);
  new_parms=0.0;
  Array1 <doublevar> local_energy(nconfig);

  wfdata->getVarParms(parms);
  nfunctions=wf->nfunc();

  Array1 <doublevar> bestparms(nparms);
  doublevar energy_mean; 
  doublevar energy_mean_err=0;
  doublevar variance_mean=0;
  doublevar weight_variance_mean=0;
  doublevar function_mean=0;


  string vmcout=options.runid+"_mc.o";
  ofstream vmcoutput;
  if( mpi_info.node==0 )
    vmcoutput.open(vmcout.c_str());
  doublevar best_function, best_energy, best_energy_err, best_variance, best_weight_variance;
  best_function=best_energy=best_energy_err=best_variance=best_weight_variance=0;
  int iter_best=0;

  //start iterations
  int iter=0;
  int reached_max_blocks=0;
  if( mpi_info.node==0 )
    cout <<"start iterations"<<endl;

  while(iter<iterations && !reached_max_blocks){
    if( mpi_info.node==0 ){
      cout <<endl;
      cout <<"############### iteration "<<iter<<" #############################"<<endl<<endl;;
    }

    /*
    if(iter==0){
      //do things for the first time
      get_averages(parms, iter, vmcoutput, options, function_mean, energy_mean, energy_mean_err, variance_mean, weight_variance_mean, new_parms);

      //print the wavefunction
      if( mpi_info.node==0 )
	wf_printout(parms, iter, function_mean, energy_mean, energy_mean_err, variance_mean, weight_variance_mean, output);

      readcheck(storeconfig);
      parms=new_parms;
      bestparms=parms;
      best_function=function_mean+2.0*energy_mean_err;
      best_energy=energy_mean;
      best_variance=variance_mean;
      best_weight_variance=weight_variance_mean;
      iter++;
    }
    */

    old_delta_parms=delta_parms;
    get_averages(parms, iter, vmcoutput, options, function_mean, energy_mean, energy_mean_err, variance_mean, weight_variance_mean, new_parms);
    /*
    if( mpi_info.node==0 ){
      for(int k=0; k<nparms; k++){
	cout<<k+1<<": old= "<< old_delta_parms(k)<<" new= "<<delta_parms(k)<<endl;
      }
    }
    */
    doublevar scalar_product_of_delta_lamda=dot(old_delta_parms,delta_parms);
    doublevar norm_old=dot(old_delta_parms,old_delta_parms);
    doublevar norm_new=dot(delta_parms,delta_parms);
    scalar_product_of_delta_lamda/=sqrt(norm_old*norm_new);
    if( mpi_info.node==0 ){
      //cout <<" dlabmba_old.dlabmba_new "<<scalar_product_of_delta_lamda<<endl;
      wf_printout(parms, iter, function_mean, energy_mean, energy_mean_err, variance_mean, weight_variance_mean, output);
    }
    
    if( mpi_info.node==0 ){
      cout << "iteration  function_mean      energy_mean     variance_mean   weight_variance_mean"<<endl;
      cout <<"%%  "<<iter<<"  "<<function_mean << "  "<<energy_mean<<" +/- "<<energy_mean_err<<"   "<<variance_mean<<"  "<<weight_variance_mean<<endl;
    }
        
    //change number of blocks in the last MC
    int nblock_pos=1;
    int last_mc_words=mc_words.size()-1;
    while(mc_words[last_mc_words][nblock_pos++]!="NBLOCK"){}
    //cout <<"NBLOCK "<<mc_words[last_mc_words][nblock_pos]<<endl;
    int nblocks=atoi(mc_words[last_mc_words][nblock_pos].c_str());
    if(nblocks==nblocks_max){
      output<<"maximum number of blocks reached with no further improvement, exiting"<<endl;
      reached_max_blocks=1;
      //break;
    }
    if(nblocks < nblocks_max ){
      if(scalar_product_of_delta_lamda>=0)
	nblocks+=1;
      else
	nblocks=int(block_increase_factor*nblocks);
      if(nblocks>nblocks_max)
	nblocks=nblocks_max;
      char strbuff3[40];
      sprintf(strbuff3,  "%d", nblocks);
      mc_words[last_mc_words][nblock_pos]=strbuff3;
    }
    
    //if((function_mean+2.0*energy_mean_err) < best_function){
    if((function_mean) < best_function || iter==0){
      if( mpi_info.node==0 )
	cout << "Found better parameters"<<endl;
      
      output<< "Found better parameters"<<endl;
      best_function=function_mean;//+2.0*energy_mean_err;
      best_energy=energy_mean;
      best_energy_err=energy_mean_err;
      best_variance=variance_mean;
      best_weight_variance=weight_variance_mean;
      bestparms=parms;
      iter_best=iter;
      write_wf(bestparms, wfoutputfile, wfdata);
    }
    //overwrite the old parameters with new
    parms=new_parms;
    readcheck(storeconfig);
    iter++; //go for another iteration
  }
 
  if( mpi_info.node==0 )
    vmcoutput.close();
  write_wf(bestparms, wfoutputfile, wfdata);
  if( mpi_info.node==0 )
  {
    output << "Optimization finished.  ";
    output << "New wave function from iteration "<<iter_best<<" is in " << wfoutputfile << endl;
    output << "Final objective function: " << best_function
	   <<" and energy: "<<  best_energy <<" +/-" <<best_energy_err
	   <<" and dispersion: "<<sqrt(best_variance)
	   <<" and weight dispersion: "<<sqrt(best_weight_variance)
	   <<endl;
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
