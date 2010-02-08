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

#include "Newton_opt_method.h"
#include "qmc_io.h"
#include "Program_options.h"
#include "System.h"
#include "Array.h"
#include "ulec.h"
#include "MatrixAlgebra.h"
void Newton_opt_method::read(vector <string> words,
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
  pos=0;
  if(readvalue(words,pos, eref, "EREF")){
    eref_exists=1;
  }
  else{
    eref_exists=0;
    eref=0;
  }
  
  //Optional options
  pos=0;
  if( ! readvalue(words, pos, hessianstore, "HESSIANTEMP") )
  {
    hessianstore=options.runid+".hessian";
  }
  canonical_filename(hessianstore);

  if(! readvalue(words, pos=0, wfoutputfile, "WFOUTPUT") )
  {
    wfoutputfile=options.runid+".wfout";
  }
  else { use_weights=0; }
  if(haskeyword(words, pos=0, "EXTENDED_WFOUTPUT") )
  {
    use_extended_output=1;
  }
  else { use_extended_output=0; }
  
  if(haskeyword(words, pos=0, "SHDMC") )
    {
      calculate_projectors=1;
    }
  else {  calculate_projectors=0; }
  
  if(!readvalue(words, pos=0, vmc_nblocks, "VMC_NBLOCK") )
  {
    vmc_nblocks=1;
  }

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
  
  if(haskeyword(words, pos=0, "USE_HESSIAN_PLUS_TYPE") )
  {
    plus_version_of_hessian_of_energy=1;
  }
  else { plus_version_of_hessian_of_energy=0; } 

  if(haskeyword(words, pos=0, "STORE_HESSIAN_TO_FILE") )
  {
    store_hessian_to_file=1;
  }
  else { store_hessian_to_file=0; } 
  
  

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
  readconfig_non_cannonical=readconfig;
  canonical_filename(readconfig);
  
  if(!readvalue(words, pos=0, storeconfig, "STORECONFIG"))
    storeconfig=readconfig_non_cannonical;
  storeconfig_non_cannonical=storeconfig;
  canonical_filename(storeconfig);

  //read input stuff for MC run inside optimization
  if(haskeyword(words, pos=0, "CORRELATED_SAMPLING") )
  {
    use_correlated_sampling=1;
  }
  else { use_correlated_sampling=0; } 

  pos=0;
  vector <string> mc_words_tmp;
  while(readsection(words, pos, mc_words_tmp, "MC"))
    mc_words.push_back(mc_words_tmp);
  
  if(!mc_words.size())
    error("Need MC section inside NEWTON_OPT method");
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

void Newton_opt_method::readcheck(string & readconfig) {
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



void copy_of_2D_sym_array_to_1D_array(Array2 <doublevar> & A, Array1 <doublevar> & B){
  int size=A.GetDim(0);
  assert(A.GetDim(1)==size);
  assert(B.GetSize()>=int(size*(size+1)/2));
  int counter=0;
  for(int i=0;i<size;i++){
    for(int j=i;j<size;j++){
      B(counter++)=A(i,j);
    }
  }
}

void copy_of_1D_array_to_2D_sym_array(Array1 <doublevar> & B, Array2 <doublevar> & A){
  int size=A.GetDim(0);
  assert(A.GetDim(1)==size);
  assert(B.GetSize()>=int(size*(size+1)/2));
  int counter=0;
  for(int i=0;i<size;i++){
    for(int j=i;j<size;j++){
      A(i,j)=B(counter++);
      A(j,i)=A(i,j);
    }
  }
}



void Newton_opt_method::generate_all_quanities(Array1 <doublevar> & parms,
					       Array1 <doublevar> & delta, 
					       Array1 <doublevar> & local_energy,
					       Array1 < Array1 <doublevar> > & local_energy_gradient,
					       Array1 < Array1 <doublevar> > & wf_gradient,
					       Array1 < Array1 <doublevar> > & wf_hessian
					       ){
  Array1 <doublevar> kinetic(nfunctions), nonloc(nfunctions);
  Array1 <doublevar> Psi(nparms);
  doublevar energy;
  doublevar psi;
  Wf_return wfval(nfunctions,2);
  Parm_deriv_return derivatives;
  derivatives.nparms_start=0;
  derivatives.nparms_end=nparms;
  if(min_function!=min_variance){
     derivatives.need_hessian=1;
  }
  
  Pseudo_buffer psp_buff;
  doublevar ln_wf_value, sign_wf_value;

  FILE * hessianout;
  hessianout=NULL;
  if(store_hessian_to_file){
    hessianout=fopen(hessianstore.c_str(), "w");
    if(!hessianout) { error("couldn't open hessian temporary file ", hessianstore," for writing.");}
  }
    
  for(int walker=0; walker < nconfig; walker++)  {
    wfdata->setVarParms(parms);
    config_pos(walker).restorePos(sample); 
    wf->updateLap(wfdata, sample);
    wf->getVal(wfdata, 0, wfval);
    ln_wf_value=wfval.amp(0,0);
    sign_wf_value=wfval.sign(0);
    calcpot(walker)=sysprop->calcLoc(sample);
    psp_buff.clear();
    pseudo->initializeStatic(wfdata, sample, wf, psp_buff);
    wf->notify(sample_static,0);
    /* remove comments to make use of nonloc_deriv
    Array1 <doublevar> psp_test(pseudo->nTest());
    for(int i=0; i< pseudo->nTest(); i++) { 
      psp_test(i)=rng.ulec();
      }
    */

    if(wfdata->supports(parameter_derivatives)){
      wf->getParmDeriv(wfdata, sample, derivatives);
      wf_gradient(walker)=derivatives.gradient;
      if(store_hessian_to_file){
	Array1 <doublevar> wf_hessian_store(int(nparms*(nparms+1)/2));
	copy_of_2D_sym_array_to_1D_array(derivatives.hessian, wf_hessian_store);
	for(int i=0; i< wf_hessian_store.GetDim(0); i++){
	  fwrite(&wf_hessian_store(i), sizeof(doublevar), 1, hessianout);
	}
      }
      else{
	copy_of_2D_sym_array_to_1D_array(derivatives.hessian,wf_hessian(walker));
      }
      /* remove comments to make use of nonloc_deriv
      Array1 <doublevar> nonloc_deriv(nparms);
      pseudo->calcNonlocParmDeriv(wfdata, sample, wf,psp_test, nonloc, nonloc_deriv);
      */
      doublevar old_kinetic;
      for (int i=0;i<=nparms;i++){
        if(i>0){
          parms(i-1)+=delta(i-1);
          wfdata->setVarParms(parms);
        }
	
        wf->updateLap(wfdata, sample);
        sysprop->calcKinetic(wfdata, sample, wf, kinetic);
        psp_buff.start_again();
        pseudo->calcNonlocWithFile(wfdata, sample, wf,nonloc, psp_buff);
	energy=kinetic(0)+calcpot(walker)+nonloc(0);
	if(i==0){
	  /* remove comments to make use of nonloc_deriv
	     energy=kinetic(0)+calcpot(walker)+nonloc(0);
	  */
	  local_energy(walker)=energy;
	  old_kinetic=kinetic(0);
	}
        else 
	  local_energy_gradient(walker)(i-1)=(energy-local_energy(walker))/delta(i-1);
	/* remove comments to make use of nonloc_deriv
	   local_energy_gradient(walker)(i-1)=nonloc_deriv(i-1)+(kinetic(0)-old_kinetic)/delta(i-1);
	*/
        if(i>0)
          parms(i-1)-=delta(i-1);
      }
    }
    else{
      Array2 <doublevar> hessian(nparms,nparms);
      //hessian=0;
      for (int i=0;i<=nparms;i++){
        if (i>0)
          parms(i-1)+=delta(i-1);
        for (int j=i;j<=nparms;j++){
          if (j>0)
            parms(j-1)+=delta(j-1);
          wfdata->setVarParms(parms);
          
          if(i==0){
            wf->updateLap(wfdata, sample);
            wf->getVal(wfdata, 0, wfval);
            
            sysprop->calcKinetic(wfdata, sample, wf, kinetic);
            psp_buff.start_again();
            pseudo->calcNonlocWithFile(wfdata, sample, wf,nonloc, psp_buff);
            energy=kinetic(0)+calcpot(walker)+nonloc(0);
            if(j==0){
              local_energy(walker)=energy;
            }
            else{
              Psi(j-1)=wfval.sign(0)*sign_wf_value*exp(wfval.amp(0,0)-ln_wf_value);
              local_energy_gradient(walker)(j-1)=(energy-local_energy(walker))/delta(j-1);
              wf_gradient(walker)(j-1)=(Psi(j-1)-1)/delta(j-1);
            }
          }// i>0 j>=i calculation of double index quantities
          else {
            wf->updateVal(wfdata, sample);
            wf->getVal(wfdata, 0, wfval);
            psi=wfval.sign(0)*sign_wf_value*exp(wfval.amp(0,0)-ln_wf_value);
            hessian(i-1,j-1)=(psi-Psi(i-1)-Psi(j-1)+1)/(delta(i-1)*delta(j-1));
            //cout<< wf_hessian(walker)(i-1,j-1)<< "  ";
          }
          if (j>0)
            parms(j-1)-=delta(j-1);
        }//j loop
         //cout <<endl;
        if (i>0)
          parms(i-1)-=delta(i-1);
      }//i loop
      
      if(store_hessian_to_file){
	Array1 <doublevar> wf_hessian_store(int(nparms*(nparms+1)/2));
	copy_of_2D_sym_array_to_1D_array(hessian,wf_hessian_store);
	for(int i=0; i< wf_hessian_store.GetDim(0); i++){
	  fwrite(&wf_hessian_store(i), sizeof(doublevar), 1, hessianout);
	}
      }
      else{
	copy_of_2D_sym_array_to_1D_array(hessian,wf_hessian(walker));
      }
    }//end of else (wfdata->supports(parameter_derivatives))
    wf->notify(sample_dynamic,0);
  }//walker loop
  if(store_hessian_to_file){
    fclose(hessianout);
  }
}



void Newton_opt_method::calculate_first_averages(Array1 <doublevar> & parms,
						 Array1 <doublevar> & local_energy, 
						 Array1 < Array1 <doublevar> > local_energy_gradient,
						 Array1 < Array1 <doublevar> > wf_gradient,
						 doublevar & energy_mean,
						 doublevar & function_mean_est,
						 doublevar & energy_mean_est, 
						 doublevar & variance_mean_est,
						 Array1 <doublevar> & wf_gradient_mean,
						 Array1 <doublevar> & local_energy_gradient_mean,
						 Array1 <doublevar> & energy_gradient_mean
						 ){
  Array1 <doublevar> wf_grad_sum(nparms);
  Array1 <doublevar> local_energy_grad_sum(nparms);
  Array1 <doublevar> energy_grad_sum(nparms);
  wf_grad_sum=0;
  local_energy_grad_sum=0;
  energy_grad_sum=0;
  doublevar energy_sum=0;
  doublevar variance_sum=0;
  doublevar weightsum=0;
  doublevar energy,weight;
  for(int walker=0; walker < nconfig; walker++){
    energy=local_energy(walker);
    weight=dmc_weight(walker);
    energy_sum+=energy*weight;
    variance_sum+=(energy-eref)*(energy-eref)*weight;
    for(int i=0;i<nparms;i++){
      energy_grad_sum(i)+=2.0*wf_gradient(walker)(i)*(energy-energy_mean)*weight;
      wf_grad_sum(i)+=wf_gradient(walker)(i)*weight;
      local_energy_grad_sum(i)+=local_energy_gradient(walker)(i)*weight;
    }
    weightsum+=weight;
  }
  doublevar weight_sum=parallel_sum(weightsum);
  energy_mean_est=parallel_sum(energy_sum)/weight_sum;
  variance_mean_est=parallel_sum(variance_sum)/weight_sum;
  switch(min_function)
    {
    case min_variance:
      function_mean_est=variance_mean_est;
      break;
    case min_energy:
      function_mean_est=energy_mean_est;
      break;
    case min_mixed:
      function_mean_est=mixing*energy_mean_est+(1-mixing)*variance_mean_est;
      break; 
    case min_weight_variance:
      error("Newton_opt_method::min_function=WEIGHT_VARIANCE not implemented");
    default:
      error("Newton_opt_method::min_function has a very strange value");
    }

  parallel_sum(energy_grad_sum);
  parallel_sum(wf_grad_sum);
  parallel_sum(local_energy_grad_sum);
  for(int i=0;i<nparms;i++){
    energy_gradient_mean(i)=energy_grad_sum(i)/weight_sum;
    wf_gradient_mean(i)=wf_grad_sum(i)/weight_sum;
    local_energy_gradient_mean(i)=local_energy_grad_sum(i)/weight_sum;
  }
}

void Newton_opt_method::calculate_second_averages(Array1 <doublevar> & parms, 
						  Array1 <doublevar> & local_energy, 
						  Array1 < Array1 <doublevar> > local_energy_gradient,
						  Array1 < Array1 <doublevar> > wf_gradient,
						  Array1 < Array1 <doublevar> > wf_hessian,
						  doublevar & energy_mean,
						  Array1 <doublevar> & wf_gradient_mean,
						  Array1 <doublevar> & local_energy_gradient_mean,
						  Array1 <doublevar> & energy_gradient_mean,
						  Array1 <doublevar> & variance_grad_mean,
						  Array2 <doublevar> & hess1_eng_mean,
						  Array2 <doublevar> & hess2_eng_mean,
						  Array2 <doublevar> & hess3_1eng_mean,
						  Array2 <doublevar> & hess3_2eng_mean,
						  Array2 <doublevar> & hess1_var_mean
						  ){
  
  Array1 <doublevar> variance_grad(nparms);
  Array2 <doublevar> hess1_eng(nparms,nparms);
  Array2 <doublevar> hess2_eng(nparms,nparms); 
  Array2 <doublevar> hess3_1eng(nparms,nparms);
  Array2 <doublevar> hess3_2eng(nparms,nparms);
  Array2 <doublevar> hess1_var(nparms,nparms);
  Array2 <doublevar> hessian_tmp(nparms,nparms);
  variance_grad=0;
  hess1_eng=hess2_eng=hess3_1eng=hess3_2eng=hess1_var=0;
  doublevar weightsum=0;
  doublevar energy,weight;
  
  FILE * hessianin;
  hessianin=NULL;
  if(store_hessian_to_file){
    hessianin=fopen(hessianstore.c_str(), "r");
    if(!hessianin) { error("couldn't open hessian temporary file ", hessianstore," for reading.");}
  }

  for(int walker=0; walker < nconfig; walker++){
    energy=local_energy(walker);
    weight=dmc_weight(walker);
    if(store_hessian_to_file){
      Array1 <doublevar> wf_hessian_store(int(nparms*(nparms+1)/2));
      for(int j=0; j< wf_hessian_store.GetDim(0); j++) {
	fread(&wf_hessian_store(j), sizeof(doublevar), 1, hessianin);
      }
      copy_of_1D_array_to_2D_sym_array(wf_hessian_store,hessian_tmp);
    }
    else{
      copy_of_1D_array_to_2D_sym_array(wf_hessian(walker),hessian_tmp);
    }
    for(int i=0;i<nparms;i++){
      for (int j=i;j<nparms;j++){
	if(i==0){
	  variance_grad(j)+=2.0*(local_energy_gradient(walker)(j)-energy_gradient_mean(j))*(energy-eref)*weight;
	}
	if(!plus_version_of_hessian_of_energy){
	  hess1_eng(i,j)+=(hessian_tmp(i,j)-wf_gradient(walker)(i)*wf_gradient(walker)(j))*
	                  (energy-energy_mean)*weight;
	  hess2_eng(i,j)+=(wf_gradient(walker)(i)-wf_gradient_mean(i))*
	                  (wf_gradient(walker)(j)-wf_gradient_mean(j))*
	                  (energy-energy_mean)*weight;
	}
	else{
	  hess1_eng(i,j)+=(hessian_tmp(i,j)+wf_gradient(walker)(i)*wf_gradient(walker)(j))*
	                  (energy-energy_mean)*weight;
	  //hess2_eng(i,j)+=wf_gradient_mean(i)*energy_gradient_mean(j)+
	  //              wf_gradient_mean(j)*energy_gradient_mean(i);
	}
	hess3_1eng(i,j)+=wf_gradient(walker)(i)*local_energy_gradient(walker)(j)*weight;
	hess3_2eng(i,j)+=wf_gradient(walker)(j)*local_energy_gradient(walker)(i)*weight;
	hess1_var(i,j)+=(local_energy_gradient(walker)(i)-energy_gradient_mean(i))*
	  (local_energy_gradient(walker)(j)-energy_gradient_mean(j))*weight;
	 
       }
    }
    weightsum+=weight;
  }//walker loop
  if(store_hessian_to_file){
    fclose(hessianin);
  }
  doublevar weight_sum=parallel_sum(weightsum);
  parallel_sum(variance_grad);
  parallel_sum(hess1_eng);
  parallel_sum(hess2_eng);
  parallel_sum(hess3_1eng);
  parallel_sum(hess3_2eng);
  parallel_sum(hess1_var);
  for(int i=0;i<nparms;i++){
    for (int j=i;j<nparms;j++){
      if(i==0){
	 variance_grad_mean(j)=variance_grad(j)/weight_sum;
      }
      hess1_eng_mean(i,j)=hess1_eng(i,j)/weight_sum;
      if(!plus_version_of_hessian_of_energy)
	hess2_eng_mean(i,j)=hess2_eng(i,j)/weight_sum;
      else
	hess2_eng_mean(i,j)=wf_gradient_mean(i)*energy_gradient_mean(j)+
	                    wf_gradient_mean(j)*energy_gradient_mean(i);
      hess3_1eng_mean(i,j)=hess3_1eng(i,j)/weight_sum;
      hess3_2eng_mean(i,j)=hess3_2eng(i,j)/weight_sum;
      hess1_var_mean(i,j)=hess1_var(i,j)/weight_sum;
    }
  }
  
}




void Newton_opt_method::build_gradient(Array1 <doublevar> & energy_gradient_mean,
				       Array1 <doublevar> & variance_grad_mean,
				       Array1 <doublevar> & grad_var
				       ){
  for (int i=0;i<nparms;i++){
    switch(min_function)
      {
      case min_energy:
	grad_var(i)=energy_gradient_mean(i);
	break;
      case min_variance:
	grad_var(i)=variance_grad_mean(i);
	break;
      case min_mixed:
	grad_var(i)=mixing*energy_gradient_mean(i)
	  +(1-mixing)*variance_grad_mean(i);
	break;
      }
  }
}

void Newton_opt_method::build_hessian(Array1 <doublevar> & wf_gradient_mean,
				      Array1 <doublevar> & local_energy_gradient_mean,
				      Array2 <doublevar> & hess1_eng_mean,
				      Array2 <doublevar> & hess2_eng_mean,
				      Array2 <doublevar> & hess3_1eng_mean,
				      Array2 <doublevar> & hess3_2eng_mean,
				      Array2 <doublevar> & hess1_var_mean,
				      Array2 <doublevar> & hessian
				      ){
  doublevar h_tmp, h_tmp2, h_tmp3;
  //if( mpi_info.node==0 ){
  // cout << "Our hessian matrix"<<endl;
  //}
  for (int i=0;i<nparms;i++){
    for (int j=i;j<nparms;j++){
      switch(min_function)
        {
        case min_energy:
	  if(!plus_version_of_hessian_of_energy){
	    h_tmp=2.0*(hess1_eng_mean(i,j)+2.0*hess2_eng_mean(i,j));
	  }
	  else{
	    h_tmp=2.0*(hess1_eng_mean(i,j)-hess2_eng_mean(i,j));
	  }
          hessian(i,j)=h_tmp+hess3_1eng_mean(i,j)-wf_gradient_mean(i)*local_energy_gradient_mean(j)+
	                     hess3_2eng_mean(i,j)-wf_gradient_mean(j)*local_energy_gradient_mean(i);
          hessian(j,i)=hessian(i,j);
          break;
        case min_variance:
          hessian(i,j)=2.0*hess1_var_mean(i,j);
          hessian(j,i)=hessian(i,j);
          break;
	case min_mixed:
	  if(!plus_version_of_hessian_of_energy){
	    h_tmp=2.0*(hess1_eng_mean(i,j)+2.0*hess2_eng_mean(i,j));
	  }
	  else{
	    h_tmp=2.0*(hess1_eng_mean(i,j)-hess2_eng_mean(i,j));
	    
	  }
	  h_tmp2=h_tmp+hess3_1eng_mean(i,j)-wf_gradient_mean(i)*local_energy_gradient_mean(j)+
	               hess3_2eng_mean(i,j)-wf_gradient_mean(j)*local_energy_gradient_mean(i);
	  h_tmp3=2.0*hess1_var_mean(i,j);
	  hessian(i,j)=mixing*h_tmp2+(1-mixing)*h_tmp3;
          hessian(j,i)=hessian(i,j);  
	  break;
        default:
          error("Optimize_method2::hessian error!");
        }
      //  if( mpi_info.node==0 ){
      //cout <<  hessian(i,j)<<"  ";
      //}
    }
    //if( mpi_info.node==0 ){
    //cout << endl;
    //}
  }
  
}


void make_posite_definite(Array2 <doublevar> & function_hessian){
  int dim=function_hessian.GetDim(0);
  Array1 <doublevar> eigenvals(dim);
  Array2 <doublevar> eigenvecs(dim,dim);
  EigenSystemSolverRealSymmetricMatrix(function_hessian, eigenvals, eigenvecs);
  if( mpi_info.node==0 ){
    cout << "Eigen-values of hessian are: "<<endl;
    for(int i=0; i<dim; ++i){
      cout << i+1<<": "<<eigenvals(i)<<endl;
    }
  }
 
  if(eigenvals(dim-1)<0){
    for(int i=0; i<dim; ++i)
      function_hessian(i,i)-=eigenvals(dim-1);
  }
  
}


void newton_step(Array1 <doublevar> & function_gradient, Array2 <doublevar> & function_hessian, 
		 doublevar & damping, Array1 <doublevar> & parms, Array1 <doublevar> & newparms){
  doublevar tmp;
  int dim=function_hessian.GetDim(0);
  Array2 <doublevar> hessian_inverse(dim,dim);
  Array2 <doublevar> hessian(dim,dim);
  hessian=function_hessian;
  for(int i=0; i<dim; ++i)
    hessian(i,i)+=damping;
  
  // find new parametres by Newton method using inverse of hessian
  
  InvertMatrix(hessian, hessian_inverse, dim);
  for (int i=0;i<dim;i++){
    tmp=0;
    for (int j=0;j<dim;j++)
      tmp-=hessian_inverse(i,j)*function_gradient(j);
    newparms(i)=parms(i)+tmp;
  }
}

void write_wf(Array1 <doublevar> & parms, string & wfoutputfile, Wavefunction_data * wfdata){
  if( mpi_info.node==0 ){
    string indent="";
    wfdata->setVarParms(parms);
    wfdata->renormalize();
    ofstream wfoutput(wfoutputfile.c_str());
    wfoutput.precision(15);
    wfdata->writeinput(indent,wfoutput);
    wfoutput.close();
  }
}

void Newton_opt_method::Get_correlated_energies_and_variances(Array1 < Array1 <doublevar> > & parms, 
							      int & iter,
							      Array1 <doublevar> & y, 
							      Array1 <doublevar> & energies,
							      Array1 <doublevar> & variancies,
							      ostream & output, 
							      Program_options & options
							      ){

  int size=y.GetSize();
  char strbuff[40];
  int ii=0;
  vector < vector <string> > aux_wftext(size-1);
  for(int i=0;i<size;i++){
    if(i!=1){
      sprintf(strbuff, "%d", ii);
      string wfoutputfile_tmp=wfoutputfile+"_aux"+strbuff;
      write_wf(parms(i), wfoutputfile_tmp, wfdata);
#ifdef USE_MPI
      MPI_Barrier(MPI_Comm_grp);
#endif
      ifstream wfinput(wfoutputfile_tmp.c_str());
      parsefile(wfinput, aux_wftext[ii]);
      wfinput.close();
      ii++;
    }
  }
  string wfoutputfile_tmp=wfoutputfile+"_ref";
  write_wf(parms(1), wfoutputfile_tmp, wfdata);
#ifdef USE_MPI
      MPI_Barrier(MPI_Comm_grp);
#endif

  //do correlated run with parms(0) and get energy diffrences: 
  avg_method=NULL;
  //modify mc_words before alocation READCONFIG=STORECONFIG, and PROPERTIES section
  vector <string> total_mc_words;
  char strbuff2[40];
  sprintf(strbuff2, "%d", nconfig);
  total_mc_words=mc_words[0]; 
  total_mc_words.insert(total_mc_words.begin(),storeconfig_non_cannonical);
  total_mc_words.insert(total_mc_words.begin(),"STORECONFIG");
  total_mc_words.insert(total_mc_words.begin(),storeconfig_non_cannonical);
  total_mc_words.insert(total_mc_words.begin(),"READCONFIG");
  total_mc_words.insert(total_mc_words.begin(),strbuff2);
  total_mc_words.insert(total_mc_words.begin(),"NCONFIG");
  total_mc_words.insert(total_mc_words.begin(),mc_words[0][0]);
  total_mc_words.push_back("PROPERTIES") ;
  total_mc_words.push_back("{");
  for(unsigned int i=0;i<aux_wftext.size();i++){
    total_mc_words.push_back("AUX_SYS");
    total_mc_words.push_back("{");
    for(unsigned int j=0;j<options.systemtext[0].size();j++)
      total_mc_words.push_back(options.systemtext[0][j]);
    total_mc_words.push_back("}");
    total_mc_words.push_back("AUX_WF");
    total_mc_words.push_back("{");
    for(unsigned int j=0;j<aux_wftext[i].size();j++)
      total_mc_words.push_back(aux_wftext[i][j]);
    total_mc_words.push_back("}");
  }
  total_mc_words.push_back("}");

  //for(unsigned int i=0;i<total_mc_words.size();i++)
  //  cout <<total_mc_words[i]<<" ";

  allocate(total_mc_words, options, avg_method);
  //avg_method->showinfo(output);
  string logfile=options.runid+".log";
  char strbuff3[40];
  sprintf(strbuff3, "%d", iter);
  string log_label="mc_with2_aux_func_";
  log_label+=strbuff3;
  myprop.setLog(logfile, log_label);
  ofstream logout;
  if(mpi_info.node==0 ) {
    logout.open(logfile.c_str(), ios::app);
    logout << "#-------------------------------------------------\n";
    logout << "#"<<mc_words[0][0]<<"run: iteration " <<iter<<endl;;
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
  Array1 <doublevar> rescaled_eng(size);

  wfdata->attachObserver(wf);

  energies(1)=finavg.avg(Properties_types::total_energy,0);
  variancies(1)=finavg.avgvar(Properties_types::total_energy,0);
  rescaled_eng(1)=energies(1);
  if( mpi_info.node==0 )
    cout <<"reference energy "<<energies(1)<<" +/- "<<sqrt(finavg.err(Properties_types::total_energy,0))<<" and sigma "<<sqrt(variancies(1))<<endl;
  int naux=finavg.aux_energy.GetDim(0);
  int n_cvg=finavg.aux_energy.GetDim(1);
  if(n_cvg>1)
    error("Not supported yet");
  doublevar diffrence;
  doublevar diffrence_err;
  for(int i=0; i< naux; i++) {
    for(int w=0; w< n_cvg; w++) {
      if(i==0)
	ii=i;
      else
	ii=i+1;
      energies(ii)=finavg.aux_energy(i,w);
      variancies(ii)=finavg.aux_energyvar(i,w);
      diffrence=finavg.aux_diff(i,w)/finavg.aux_size(i);
      diffrence_err=sqrt(finavg.aux_differr(i,w)/finavg.aux_size(i));
      if( mpi_info.node==0 )
	cout <<"energy_aux"<<i<<"-"<<w<<" "<<setw(15)<<setiosflags(ios::fixed | ios::showpoint)<<finavg.aux_energy(i,w)
	     <<" var_aux"<<i<<"-"<<w<<" "<<setw(15)<<setiosflags(ios::fixed | ios::showpoint)
	     <<finavg.aux_energyvar(i,w)<<" diff "<<setw(15)<<setiosflags(ios::fixed | ios::showpoint)<<diffrence<<" +/- "
	     <<setw(15)<<setiosflags(ios::fixed | ios::showpoint)<<diffrence_err<<endl;
      if(diffrence_err > abs(diffrence)){
	if( mpi_info.node==0 ){
	  cout <<" WARNING: Newton_opt_method::Get_correlated_energies_and_variances: could not get accurate energy diffrences from correlated run"<<endl; 
	  cout<<" using  more walkers could help "<<endl;
	}
      }
      rescaled_eng(ii)=energies(1)+diffrence;
    }
  }
  
  switch(min_function)
    {
    case min_energy:
      for(int i=0;i<size;i++)
	y(i)=rescaled_eng(i);
      break;
    case min_variance:
      for(int i=0;i<size;i++)
	y(i)=variancies(i);
      break;
    case min_mixed:
      for(int i=0;i<size;i++)
	y(i)=mixing*rescaled_eng(i)+(1-mixing)*variancies(i);
      break;
    default:
      error("Optimize_method2::hessian error!");
    }
}

void Newton_opt_method::adjust_distribution(Array1 <doublevar> & parms,
					    int & iter,
					    string log_label,
					    ostream & output, 
					    Program_options & options,
					    doublevar & function,
					    doublevar & energy,
					    doublevar & energy_err,
					    doublevar & variance,
					    doublevar & weight_variance
					    ){
  char strbuff[40];
  sprintf(strbuff,  "%d", nconfig);
  for (int mc_size=0; mc_size<mc_words.size(); mc_size++){
    if(calculate_projectors){
      int ndim=parms.GetSize();
      doublevar norm=0.0;
      for(int j=0; j< ndim; j++) {
	norm+=parms(j)*parms(j);
      }
      doublevar sqrtnorm=1.0/sqrt(norm);
      for(int j=0; j< ndim; j++) {
	parms(j)*=sqrtnorm;
      }
    }

    wfdata->setVarParms(parms);
    wfdata->renormalize();
    log_label=mc_words[mc_size][0];
    char strbuff2[40];
    sprintf(strbuff2,  "%d", iter);
    log_label+=strbuff2;
    vector <string>  run_words;
    run_words.push_back(mc_words[mc_size][0]);
    run_words.push_back("NCONFIG");
    run_words.push_back(strbuff);
    run_words.push_back("READCONFIG");
    if(iter<1)
      run_words.push_back(readconfig_non_cannonical);
    else
      run_words.push_back(storeconfig_non_cannonical);
    run_words.push_back("STORECONFIG");
    run_words.push_back(storeconfig_non_cannonical);


    /*
      to double the size of the block
    int nblock_pos=1;
    while(mc_words[mc_size][nblock_pos++]!="NBLOCK"){
    }
    //cout <<"NBLOCK "<<mc_words[mc_size][nblock_pos]<<endl;
    int nblocks=atoi(mc_words[mc_size][nblock_pos].c_str());
    if(mc_size==mc_words.size()-1 && nblocks < 101 ){
    nblocks*=2;
    char strbuff3[40];
    sprintf(strbuff3,  "%d", nblocks);
    mc_words[mc_size][nblock_pos]=strbuff3;
    }


    */
    
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
    
    if(calculate_projectors){
      int navg_vals=finavg.avgavg.GetDim(0);
      for(int i=0; i< navg_vals; i++) { 
	//output << "      "<<  "average_generator { " << finavg2.avgavg(i).type << " ";
	if(finavg.avgavg(i).type=="linear_der"){
	  int ndim=finavg.avgavg(i).vals.GetDim(0)-1;
	  for(int j=0; j< ndim; j++) { 
	    doublevar value=finavg.avgavg(i).vals(j)/finavg.avgavg(i).vals(ndim);
	    doublevar error=finavg.avgerr(i).vals(j)/finavg.avgavg(i).vals(ndim);
	    //output << value << " +/- "<<error<<endl;
	    if(abs(value)>2*error)
	      parms(j)=value;
	    else
	      parms(j)=0;
	  }
	  //output << " } " << endl;
	}//end if 
	else if(finavg.avgavg(i).type=="linear_delta_der"){
	  int ndim=finavg.avgavg(i).vals.GetDim(0)-1;
	  doublevar norm=0.0;
          doublevar multiply=1.8;
	  for(int j=0; j< ndim; j++) { 
	    doublevar value=parms(j)+multiply*finavg.avgavg(i).vals(j)/finavg.avgavg(i).vals(ndim);
	    if(fabs(parms(j))>0)
	      delta_parms(j)=multiply*finavg.avgavg(i).vals(j)/finavg.avgavg(i).vals(ndim); 
	    else
	      delta_parms(j)=0.0; 
	    doublevar error=finavg.avgerr(i).vals(j)/finavg.avgavg(i).vals(ndim);
	    //output << value << " +/- "<<error<<endl;
	    if(abs(value)>2*error){
	      parms(j)=value;
	      norm+=parms(j)*parms(j);
	    }
	    else
	      parms(j)=0;
	  }
	  doublevar sqrtnorm=1.0/sqrt(norm);
	  for(int j=0; j< ndim; j++) {
	    parms(j)*=sqrtnorm;
	    delta_parms(j)*=sqrtnorm;
	  }
	  //output << " } " << endl;
	}//end if 
      }//i
    }
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

doublevar parabola(doublevar x, Array1 <doublevar> & c){
  return c(2)*x*x+c(1)*x+c(0);
}

void fit_to_polynomial(Array1 <doublevar> & y,
		       Array1 <doublevar> & x,
		       Array1 <doublevar> & cof){
  int k,j,i;
  doublevar phi,ff,b;

  int n=x.GetSize();
  Array1 <doublevar> s(n);
  for (i=0;i<n;i++) s[i]=cof[i]=0.0;
  s[n-1]= -x[0];
  for (i=1;i<n;i++) {
    for (j=n-1-i;j<n-1;j++)
      s[j] -= x[i]*s[j+1];
    s[n-1] -= x[i];
  }
  for (j=0;j<n;j++) {
    phi=n;
    for (k=n-1;k>0;k--)
      phi=k*s[k]+x[j]*phi;
    ff=y[j]/phi;
    b=1.0;
    for (k=n-1;k>=0;k--) {
      cof[k] += b*ff;
      b=s[k]+x[j]*b;
    }
  }
}

doublevar parabola_min(Array1 <doublevar> & y,
		       Array1 <doublevar> & x,
		       doublevar & minimalvalue
		       ){
  Array1 <doublevar> c(y.GetSize());
  doublevar min=0;

  if(y(1)<y(0)){ 
    if(y(1)<y(2)){
      fit_to_polynomial(y,x,c);
      //if( mpi_info.node==0 )
      //	cout <<"fitting function"<<endl;
      //for(int i=0;i<y.GetSize();i++)
      //cout << x(i) <<" "<<y(i)<<endl;
      //if( mpi_info.node==0 )
      //cout<< c(2)<<"*x*x+"<<c(1)<<"*x+"<<c(0)<<endl;
      if(c(2)>0){
	min=-c(1)/(2*c(2));
	minimalvalue=parabola(min,c);
      }
      else
	error("error in parabola_min routine");
    }
    else{
      min=x(2);
      minimalvalue=y(2);
    }
  }
  else{
    if(y(2)<y(0)){
      min=x(2);
      minimalvalue=y(2);
    }
    else{
      min=x(0);
      minimalvalue=y(0);
    }
  }
  //if( mpi_info.node==0 )
    //cout <<"Optimal dumping: "<<exp(min)<<"  ";
  //if( mpi_info.node==0 )
  // cout <<"estimated minimal value "<< minimalvalue<<endl;
  return min; 
}

doublevar parabola_value(Array1 <doublevar> & y,
			 Array1 <doublevar> & x,
			 doublevar & pos
			 ){
  Array1 <doublevar> c(y.GetSize());
  fit_to_polynomial(y,x,c);
  //if( mpi_info.node==0 )
  //  cout <<"energy/variance fit"<<endl;
  //for(int i=0;i<y.GetSize();i++)
  // cout << x(i) <<" "<<y(i)<<endl;
  //if( mpi_info.node==0 )
  // cout<< c(2)<<"*x*x+"<<c(1)<<"*x+"<<c(0)<<endl;
  return parabola(pos,c);
  
}





void Newton_opt_method::wf_printout(Array1 <doublevar> & parms, int iter, 
				    doublevar & value, 
				    doublevar & energy, 
				    doublevar & energy_err, 
				    doublevar & variance, 
				    doublevar & weight_variance, 
				    int min_nconfig, 
				    doublevar mu, ostream & output){
  wfdata->setVarParms(parms);
  wfdata->renormalize();
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
	       <<" nconfig " << min_nconfig << "  damping "<<mu<< endl;
	break;
      case min_energy:
	output <<" energy= " <<setw(field)<<setiosflags(ios::fixed | ios::showpoint)<< value <<" +/- "<<setw(field)<<setiosflags(ios::fixed | ios::showpoint)
	       << energy_err <<" dispersion= "<< setw(field)<< sqrt(variance) 
	       <<" weight dispersion= " <<setw(field)<<setiosflags(ios::fixed | ios::showpoint)<<sqrt(fabs(weight_variance))
	       <<" nconfig " << min_nconfig << " damping "<< setw(field)<<mu <<endl;
	break;
      case min_mixed:
	output <<" mixing*energy+ (1-mixing)*variance= " <<setw(field)<<setiosflags(ios::fixed | ios::showpoint)<< value 
	       <<" energy= " <<setw(field)<<setiosflags(ios::fixed | ios::showpoint)<< energy 
	       <<" +/- "<<setw(field)<<setiosflags(ios::fixed | ios::showpoint)<< energy_err
	       <<" dispersion= "<<setw(field)<<setiosflags(ios::fixed | ios::showpoint)<< sqrt(variance)
	       <<" weight dispersion= " <<setw(field)<<setiosflags(ios::fixed | ios::showpoint)<<sqrt(fabs(weight_variance))
	       <<" nconfig " << min_nconfig << "  damping "<< setw(field)<<mu <<endl;
	break;  
      case min_weight_variance:
	output <<setw(field)<<" weight dispersion= " <<setw(field)<<setiosflags(ios::fixed | ios::showpoint)<<sqrt(fabs(value)) 
			  << " energy= " <<setw(field)<<setiosflags(ios::fixed | ios::showpoint)<< energy
			  <<" +/- "<< setw(field)<< setiosflags(ios::fixed | ios::showpoint) <<energy_err
	                  <<" dispersion= "<<setw(field)<<setiosflags(ios::fixed | ios::showpoint)<< sqrt(variance)
			  <<" nconfig " << min_nconfig << "  damping "<<mu<< endl;
	break;
      default:
	error("Optimize_method::variance() : min_function has a very strange value");
      }
  }
}


int Newton_opt_method::Levenberg_marquad(Array1 <doublevar> & gradient, 
					 Array2 <doublevar> & hessian, 
					 int & iter, 
					 Array1 <doublevar> & parms, 
					 doublevar & function_mean,
					 doublevar & energy_mean,
					 doublevar & energy_mean_err,
					 doublevar & variance_mean,
					 int & nu, 
					 doublevar & damping,
					 ostream & vmcoutput,
					 ostream & output,
					 Program_options & options
					 ){
  doublevar tau=0.01;
  doublevar dL;
  int nu2;
  doublevar ONE_THIRD=0.3333333334;
  doublevar LM_STOP_THRESH=1e-6;
  doublevar LM_STOP_THRESH2=1e-10;
  doublevar EPSILON=1e-12;
  doublevar eps1=LM_STOP_THRESH;
  doublevar eps2=LM_STOP_THRESH2;
  doublevar tmp;
  int dim=parms.GetSize();
  Array1 <doublevar> parms_new(dim);
  doublevar dF=0;
  doublevar function_mean_new, energy_mean_new, energy_mean_err_new, variance_mean_new;
  int flag;

  if(iter==1 || use_correlated_sampling){
    nu=2;
    tmp=0.0;
    for(int i=0; i<dim; ++i){
      if(hessian(i,i)>tmp) tmp=hessian(i,i);  //find max diagonal element 
      // if(function_hessian(i,i)<tmp2) tmp2=function_hessian(i,i);  find min diagonal element 
    }
    damping=tau*tmp;
  }
  //continue with given damping:
  while(1){
  if( mpi_info.node==0 )
    cout << "Used reference damping: "<<damping<<endl;
  if(use_correlated_sampling){
    doublevar multyplier=10.0;
    Array1 < Array1 <doublevar> > tmpparms(3);
    Array1 <doublevar> mu(3);// mulog(3);
    Array1 <doublevar> y(3), energies(3), variances(3);
    mu(0)=damping/multyplier;
    mu(0)=damping;
    //mulog(0)=log(mu(0));
    for(int i=0;i<tmpparms.GetSize();i++){
      if(i>0){
	mu(i)=multyplier*mu(i-1);
	multyplier*=multyplier;
	//mulog(i)=log(multyplier)+mulog(i-1);
      }
      //mu(i)=exp(mulog(i));
      //mulog(i)=log(mu(i));
      tmpparms(i).Resize(nparms);
      newton_step(gradient, hessian, mu(i), parms, tmpparms(i));
    }
    Get_correlated_energies_and_variances(tmpparms, iter, y, energies, variances, vmcoutput, options);
    for(int i=0;i<tmpparms.GetSize();i++){
      if( mpi_info.node==0 )
	cout << i <<"  "<<mu(i)<<"   "<< y(i) << endl;
    }
    doublevar newvalue;
    damping=parabola_min(y,mu,newvalue);
  }
  else {
     if( mpi_info.node==0 ){
       cout <<"correlated_sampling not used; keeping reference dumping"<<endl;
     }
    //damping=damping;
  }

  //found new optimal dumping get values for this step 
  if( mpi_info.node==0 )
    cout << "Damping with lowest function value: "<<damping<<endl;
  newton_step(gradient, hessian, damping, parms, parms_new);
  char strbuff[40];
  sprintf(strbuff, "%d", iter);
  string log_label="vmc_";
  log_label+=strbuff;
  doublevar weight_variance;
  adjust_distribution(parms_new, iter, log_label, vmcoutput, options, function_mean_new, energy_mean_new, energy_mean_err_new, variance_mean_new, weight_variance);
  if( mpi_info.node==0 ){
    cout << "iteration  function_mean      energy_mean     variance_mean "<<endl;
    cout <<"%%  "<<iter<<"  "<<function_mean_new << "  "<<energy_mean_new<<" +/- "<<energy_mean_err_new<<"   "<<variance_mean_new<<"  "<<endl;
  }
  
  if(use_correlated_sampling){
    wf_printout(parms_new, iter, function_mean_new, energy_mean_new, energy_mean_err_new, variance_mean_new, weight_variance, nconfig, damping, output);
    function_mean=function_mean_new;
    energy_mean=energy_mean_new;
    energy_mean_err=energy_mean_err_new;
    variance_mean=variance_mean_new;
    parms=parms_new;
    flag=0;
    break;
  }
  

  dF=function_mean-function_mean_new;
  //do convergence checks
  Array1 <doublevar>  Dp(dim);
  doublevar Dp_L2=0;
  doublevar grad_inf=0;
  doublevar p_L2=0;
  for(int i=0; i<dim; ++i){
    if(grad_inf < (tmp=abs(gradient(i)))) 
      grad_inf=tmp;
    Dp(i)=parms_new(i)-parms(i);
    Dp_L2+=Dp(i)*Dp(i);
    p_L2+=parms_new(i)*parms_new(i);
  }
  
  if((grad_inf <= eps1)){
    if(mpi_info.node==0 )
      cout <<"gradient is almost zero, you are at minimum, stop"<<endl;
    flag=1;
    break;
  }
  
  if(Dp_L2<=eps2*(p_L2+eps2)){  // relative change in p is small, stop 
      if(mpi_info.node==0 )
	cout <<"relative change in p is small, stop"<<endl;
      flag=1;
      break;
  }
  
  if(Dp_L2>=(p_L2+eps2)/(EPSILON*EPSILON)){ // almost singular 
    if(mpi_info.node==0 )
      cout <<"almost singular, stop"<<endl;
    flag=1;
    break;
    }
  
  dL=0.0;
  for(int i=0; i<dim; ++i)
    dL+=Dp(i)*(damping*Dp(i)-gradient(i));
  if( mpi_info.node==0 )
    cout << "dL "<<dL<<" and dF "<<dF<<endl;
  //check for acceptance
  if(dL>0.0 && dF>0.0 ){
    if(mpi_info.node==0 )
      cout << "reduction in error, step is accepted"<<endl;  
    wf_printout(parms_new, iter, function_mean_new, energy_mean_new, energy_mean_err_new, variance_mean_new, weight_variance, nconfig, damping, output);
    function_mean=function_mean_new;
    energy_mean=energy_mean_new;
    energy_mean_err=energy_mean_err_new;
    variance_mean=variance_mean_new;
    parms=parms_new;
    //adjust dumping
    tmp=(4.0*dF/dL-1.0);
    tmp=1.0-tmp*tmp*tmp;
    damping=damping*( (tmp>=ONE_THIRD)? tmp : ONE_THIRD );
    nu=2;
    flag=0;
    break;
  }
  else if(dL> -5.0*abs(energy_mean_err_new) && dF>-5.0*abs(energy_mean_err_new)){
    if(mpi_info.node==0 )
	cout << "No reduction in error, but step is accepted"<<endl;  
    wf_printout(parms_new, iter, function_mean_new, energy_mean_new, energy_mean_err_new, variance_mean_new, weight_variance, nconfig, damping, output);
    function_mean=function_mean_new;
    energy_mean=energy_mean_new;
    energy_mean_err=energy_mean_err_new;
    variance_mean=variance_mean_new;
    parms=parms_new;
    nu=2;
    flag=0;
    break;
  }

  if(mpi_info.node==0 )
    cout << "No reduction in error, step is not accepted, retrying with larger dumping"<<endl;  
  damping*=nu;
  nu2=nu<<1; // 2*nu;
  if(nu2<=nu){ // nu has wrapped around (overflown). 
    // Thanks to Frank Jordan for spotting this case 
    if( mpi_info.node==0 )
      cout << "nu has wrapped around (overflown)"<<endl;
    flag=1;
    break;
  }
  nu=nu2;    
  }
  return flag;
}





int Newton_opt_method::showinfo(ostream & os)
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
  os << "Wave function optimization " << endl;
  if(calculate_projectors)
    os << "using the SH-DMC projectors for parameters"<<endl;
  os << "Configurations per processor: " <<  nconfig   << endl;
  os << "Number of processors: "        <<  mpi_info.nprocs << endl;
  os << "Total configurations: " <<          nconfig*mpi_info.nprocs << endl;
  os << "Iterations : "  << iterations << endl;
  os << "Number of parameters " << wfdata->nparms() << endl;
  os << "Reference energy : "           << eref << endl;
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
  if(use_correlated_sampling)
    os << "correlated sampling used for damping" << endl; 
  else
    os << "self-adapting method for damping" << endl; 
    
  os << "----------------------------" << endl;
  
  
  return 1;
}

void Newton_opt_method::run(Program_options & options, ostream & output)
{
  output.precision(10);
  nparms=wfdata->nparms(); //Number of variables
  if(nparms<= 0 ) error("There appear to be no parameters to optimize!");
  Array1 <doublevar> delta(nparms);
  delta=0.0001; //decent value fo finite differences
  delta_parms.Resize(nparms);
  delta_parms=0.0; 
  old_delta_parms.Resize(nparms);
  old_delta_parms=0.0;
  Array1 <doublevar> parms(nparms);
  Array1 <doublevar> local_energy(nconfig);
  Array1 < Array1 <doublevar> > local_energy_gradient(nconfig);
  Array1 < Array1 <doublevar> > wf_gradient;
  Array1 < Array1 <doublevar> > wf_hessian;
  calcpot.Resize(nconfig);
  wf_gradient.Resize(nconfig);
  wf_hessian.Resize(nconfig);
  
  for(int walker=0;walker<nconfig;walker++){
    local_energy_gradient(walker).Resize(nparms);
    wf_gradient(walker).Resize(nparms);
    if(!store_hessian_to_file)
      wf_hessian(walker).Resize(int(nparms*(nparms+1)/2));
  }
  
  wfdata->getVarParms(parms);
  nfunctions=wf->nfunc();
  
  
  Array1 <doublevar> wf_gradient_mean(nparms);
  Array1 <doublevar> local_energy_gradient_mean(nparms);
  Array1 <doublevar> energy_gradient_mean(nparms);
  Array1 <doublevar> variance_gradient_mean(nparms);
  Array2 <doublevar> hess1_eng_mean(nparms,nparms),
                     hess2_eng_mean(nparms,nparms),
                     hess3_1eng_mean(nparms,nparms),
                     hess3_2eng_mean(nparms,nparms),
                     hess1_var_mean(nparms,nparms);

  Array2 <doublevar> function_hessian(nparms,nparms);
  Array1 <doublevar> function_gradient(nparms);

  doublevar damping;
  int nu=2;
  Array1 <doublevar> bestparms(nparms);
  doublevar bestdamping=0;
  doublevar energy_mean; 
  doublevar energy_mean_err=0;
  //doublevar energy_mean_fit=0;
  doublevar energy_mean_est=0;
  //doublevar variance_mean_fit=0;
  doublevar variance_mean=0;
  doublevar variance_mean_est=0;
  doublevar weight_variance_mean=0;
  doublevar function_mean=0;
  //doublevar function_mean_fit=0; 
  doublevar function_mean_est=0;
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
    
    if(iter==0){
      string log_label="vmc_-1";
      adjust_distribution(parms, iter, log_label, vmcoutput, options, function_mean, energy_mean, energy_mean_err, variance_mean, weight_variance_mean);
      readcheck(storeconfig);
      if(!eref_exists)
	eref=energy_mean;
      wf_printout(parms, iter, function_mean, energy_mean, energy_mean_err, variance_mean, weight_variance_mean, nconfig, bestdamping, output);
      bestparms=parms;
      best_function=function_mean+2.0*energy_mean_err;
      best_energy=energy_mean;
      best_variance=variance_mean;
      best_weight_variance=weight_variance_mean;
      iter++;
    }
    
    if(!calculate_projectors){
      if( mpi_info.node==0 )
	cout << "generate_all_quanities" <<endl;
      generate_all_quanities(parms,
                             delta, 
                             local_energy,
                             local_energy_gradient,
                             wf_gradient,
                             wf_hessian
                             );
      
      if(use_weights || !eref_exists)
	eref=energy_mean;
      
      if( mpi_info.node==0 )
	cout << "calculate_first_averages" <<endl;
      calculate_first_averages(parms,
			       local_energy, 
			       local_energy_gradient,
			       wf_gradient,
			       energy_mean,
			       function_mean_est,
			       energy_mean_est,
			       variance_mean_est,
			       wf_gradient_mean,
			       local_energy_gradient_mean,
			       energy_gradient_mean
			       );
      if( mpi_info.node==0 )
	cout << "calculate_second_averages" <<endl;
      calculate_second_averages(parms, 
				local_energy, 
				local_energy_gradient,
				wf_gradient,
				wf_hessian,
				energy_mean,
				wf_gradient_mean,
				local_energy_gradient_mean,
				energy_gradient_mean,
				variance_gradient_mean,
				hess1_eng_mean,
				hess2_eng_mean,
				hess3_1eng_mean,
				hess3_2eng_mean,
				hess1_var_mean
				);
      
      if( mpi_info.node==0 )
	cout << "build gradient" <<endl;
      build_gradient(energy_gradient_mean,
		     variance_gradient_mean,
		     function_gradient);
      
      if(mpi_info.node==0 ){
	cout << "Gradients of: minimized quantity, vmc energy, wavefunction and local energy"<<endl;
	for (int i=0;i<nparms;i++){
	  cout <<i+1<<":  "<<  function_gradient(i)<<"   "<<energy_gradient_mean(i)<<"  "<<wf_gradient_mean(i)<<"  "<<local_energy_gradient_mean(i)<<endl;
	}
      }
      if( mpi_info.node==0 )
	cout << "build hessian" <<endl;
      build_hessian(wf_gradient_mean,
		    local_energy_gradient_mean,
		    hess1_eng_mean,
		    hess2_eng_mean,
		    hess3_1eng_mean,
		    hess3_2eng_mean,
		    hess1_var_mean,
		    function_hessian);
      
      make_posite_definite(function_hessian);
      
      if( mpi_info.node==0 )
	cout << "determine reference damping by adaptive method" <<endl;
      
      if(Levenberg_marquad(function_gradient, 
			   function_hessian, 
			   iter, 
			   parms, 
			   function_mean,
			   energy_mean,
			   energy_mean_err,
			   variance_mean,
			   nu, 
			   damping,
			   vmcoutput,
			   output,
			   options)){
	break;
      }
      
    }//calculate only if projectors in not set
    else{ //calculate only if projectors is set
      char strbuff[40];
      sprintf(strbuff, "%d", iter);
      string log_label="vmc_";
      log_label+=strbuff;
      old_delta_parms=delta_parms;
      adjust_distribution(parms, iter, log_label, vmcoutput, options, function_mean, energy_mean, energy_mean_err, variance_mean, weight_variance_mean);
      if( mpi_info.node==0 ){
	for(int k=0; k<nparms; k++){
	  cout<<k+1<<": old= "<< old_delta_parms(k)<<" new= "<<delta_parms(k)<<endl;
	}
      }
      doublevar scalar_product_of_delta_lamda=dot(old_delta_parms,delta_parms);
      doublevar norm_old=dot(old_delta_parms,old_delta_parms);
      doublevar norm_new=dot(delta_parms,delta_parms);
      scalar_product_of_delta_lamda/=sqrt(norm_old*norm_new);
      if( mpi_info.node==0 ){
	cout <<" dlabmba_old.dlabmba_new "<<scalar_product_of_delta_lamda<<endl;
	wf_printout(parms, iter, function_mean, energy_mean, energy_mean_err, variance_mean, weight_variance_mean, nconfig, damping, output);
      }
      
      if( mpi_info.node==0 ){
	cout << "iteration  function_mean      energy_mean     variance_mean   weight_variance_mean"<<endl;
	cout <<"%%  "<<iter<<"  "<<function_mean << "  "<<energy_mean<<" +/- "<<energy_mean_err<<"   "<<variance_mean<<"  "<<weight_variance_mean<<endl;
      }

      
      //change number of blocks in the last MC
      //if(function_mean > best_function){
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
	    nblocks=int(1.5*nblocks);
	  //  nblocks*=2;
	  if(nblocks>nblocks_max)
	    nblocks=nblocks_max;
	  char strbuff3[40];
	  sprintf(strbuff3,  "%d", nblocks);
	  mc_words[last_mc_words][nblock_pos]=strbuff3;
	}
	
	//}
    }


    //if((function_mean+2.0*energy_mean_err) < best_function){
    if((function_mean) < best_function){
      if( mpi_info.node==0 )
	cout << "Found better parameters"<<endl;
      
      output<< "Found better parameters"<<endl;
      best_function=function_mean;//+2.0*energy_mean_err;
      best_energy=energy_mean;
      best_energy_err=energy_mean_err;
      best_variance=variance_mean;
      best_weight_variance=weight_variance_mean;
      bestdamping=damping;
      bestparms=parms;
      iter_best=iter;
      write_wf(bestparms, wfoutputfile, wfdata);
    }
    readcheck(storeconfig);
    iter++; //go for another iteration
  } 
  if( mpi_info.node==0 )
    vmcoutput.close();
  remove(hessianstore.c_str());
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


