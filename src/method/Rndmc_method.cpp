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



#include "Rndmc_method.h"
#include "qmc_io.h"
#include "ulec.h"
#include "Program_options.h"
#include "average.h"


void Rndmc_method::read(vector <string> words,
                      unsigned int & pos,
                      Program_options & options)
{
  ndim=3;

  have_read_options=1;

  vector <string> offsetwords;

  //required options

  
  if(!readvalue(words, pos=0, nblock, "NBLOCK"))
    error("Need NBLOCK in METHOD section");
  
  if(!readvalue(words, pos=0, nconfig, "NCONFIG"))
    error("Need NCONFIG in METHOD section");

  if(!readvalue(words, pos=0, nstep, "NSTEP"))
    error("Need NSTEP in METHOD section");

  if(!readvalue(words, pos=0, timestep, "TIMESTEP"))
    error("Need TIMESTEP in METHOD section");

  if(readvalue(words, pos=0, readconfig, "READCONFIG"))
    canonical_filename(readconfig);
  else
    error("Must give READCONFIG for DMC");

  //optional options

  if(!readvalue(words, pos=0, eref, "EREF"))
    eref=0.0;

  if(haskeyword(words, pos=0, "CDMC")) do_cdmc=1;
  else do_cdmc=0;
  if(haskeyword(words, pos=0, "TMOVES")) tmoves=1; 
  else tmoves=0;

  if(!readvalue(words, pos=0, nhist, "CORR_HIST")) 
    nhist=-1;

  if(readvalue(words, pos=0, storeconfig, "STORECONFIG"))
    canonical_filename(storeconfig);

  if(!readvalue(words, pos=0, log_label, "LABEL"))
    log_label="dmc";

  if(!readvalue(words, pos=0, start_feedback, "START_FEEDBACK"))
    start_feedback=1;

  if(readvalue(words, pos=0, feedback_interval, "FEEDBACK_INTERVAL")) {
    if(feedback_interval < 1) 
      error("FEEDBACK_INTERVAL must be greater than or equal to 1");
  }
  else feedback_interval=5;

  if(!readvalue(words, pos=0, feedback, "FEEDBACK"))
    feedback=1.0;

  if(!readvalue(words, pos=0, branch_start_cutoff, "BRANCH_START_CUTOFF")) 
    branch_start_cutoff=10;
  

  branch_stop_cutoff=branch_start_cutoff*1.5;
  
  
  vector <string> proptxt;
  if(readsection(words, pos=0, proptxt, "PROPERTIES")) 
    //    myprop.read(proptxt, options.systemtext[0], options.twftext[0]);
    mygather.read(proptxt);

  vector<string> tmp_dens;
  pos=0;
  while(readsection(words, pos, tmp_dens, "DENSITY")) {
    dens_words.push_back(tmp_dens);
  }

  pos=0;
  while(readsection(words, pos, tmp_dens, "NONLOCAL_DENSITY")) {
    nldens_words.push_back(tmp_dens);
  }

  pos=0;
  while(readsection(words, pos, tmp_dens, "AVERAGE")) {
    avg_words.push_back(tmp_dens);
  }
  
  vector <string> dynamics_words;
  if(!readsection(words, pos=0, dynamics_words, "DYNAMICS") ) 
    dynamics_words.push_back("SPLIT");

  low_io=0;
  if(haskeyword(words, pos=0,"LOW_IO")) low_io=1;


  //MB: reading in the alpha relative value
  if(!readvalue(words, pos=0, alpha, "ALPHA"))
    error("supply value of ALPHA for node release ");

  //MB: needs to be possitive
  if(alpha<0 && alpha >1 )
    error("need positive value of ALPHA for node release");

  allocate(dynamics_words, dyngen);

  //MB: by default enforceNodes(0) so if 
  //MB: alpha =0 we want the fixed node
  if(alpha==0.0){
    dyngen->enforceNodes(1);
  }
}

//----------------------------------------------------------------------

int Rndmc_method::generateVariables(Program_options & options) {

  if(!have_read_options) 
    error("need to call Dmc_method::read before generateVariables");
  if(have_allocated_variables) 
    error("already allocated variables in Dmc_method");

  have_allocated_variables=1;
  allocate(options.systemtext[0], mysys);
  mysys->generatePseudo(options.pseudotext, mypseudo);
  allocate(options.twftext[0], mysys, mywfdata);
  
  densplt.Resize(dens_words.size());
  for(int i=0; i< densplt.GetDim(0); i++) {
    allocate(dens_words[i], mysys, options.runid,densplt(i));
  }
  nldensplt.Resize(nldens_words.size());
  for(int i=0; i< nldensplt.GetDim(0); i++) {
    allocate(nldens_words[i], mysys, options.runid,nldensplt(i));
  }
  
  return 1;
}

//----------------------------------------------------------------------



int Rndmc_method::allocateIntermediateVariables(System * sys,
                                              Wavefunction_data * wfdata) {
  
  if(wf) delete wf;
  wf=NULL;
  if(sample) delete sample;
  sample=NULL;
  nelectrons=sys->nelectrons(0)+sys->nelectrons(1);
  wfdata->generateWavefunction(wf);
  sys->generateSample(sample);
  sample->attachObserver(wf);
  nwf=wf->nfunc();

  if(wf->nfunc() >1) 
    error("DMC doesn't support more than one guiding wave function");

  //MB: new node release guiding function defined in Guiding_wavefunction.h
  guidingwf= new Primary_noderelease;
  doublevar tmp_alpha=0.0;
  guidingwf->set_alpha_for_noderelease(tmp_alpha);
  //MB: allocating the rndm points
  pts.Resize(nconfig);
  pts_unbiased.Resize(nconfig, nstep);
  for(int i=0; i < nconfig; i++) {
    pts(i).age.Resize(nelectrons);
    pts(i).age=0;
    for(int s=0; s < nstep; s++){ 
      pts_unbiased(i,s).age.Resize(nelectrons);
      pts_unbiased(i,s).age=0;
    }
  }
  
  average_var.Resize(avg_words.size());
  average_var=NULL;
  for(int i=0; i< average_var.GetDim(0); i++) { 
    allocate(avg_words[i], sys, wfdata, average_var(i));
  }
  
  return 1;
}

//----------------------------------------------------------------------

int Rndmc_method::showinfo(ostream & os)
{

  if(have_allocated_variables) {
    mysys->showinfo(os);
    mypseudo->showinfo(os);
    mywfdata->showinfo(os);
  }
  os << "###########################################################\n";
  os << "Release Node Diffusion Monte Carlo:\n";
  os << "Number of processors " <<           mpi_info.nprocs << endl;
  os << "Blocks: " <<                        nblock    << endl;
  os << "Steps per block: " <<               nstep     << endl;
  os << "Timestep: " <<                      timestep  << endl;
  if(tmoves) 
    os << "T-moves turned on" << endl;
  os << "Alpha: "<< alpha << endl;
  string indent="  ";

  dyngen->showinfo(indent, os);

  os << "###########################################################" << endl;
  return 1;
}

//----------------------------------------------------------------------

void Rndmc_method::find_cutoffs() {
  doublevar eaverage=0;
  
  doublevar totweight=0;
  for(int i=0; i< nconfig; i++) {
      //eaverage+=(prop.trace(step,i).energy(w)+offset(w))
      //  *guidingwf->getOperatorWeight(prop.trace(step,i).wf_val,w);
    eaverage+=pts(i).prop.energy(0)*pts(i).weight;
    totweight+=pts(i).weight;
    //cout << "en " << pts(i).prop.energy(0) << endl;
  }
  //cout << mpi_info.node << "par sum " << nconfig << endl; 
  totweight=parallel_sum(totweight);
  //cout << mpi_info.node << "eaverage " << endl;
  eaverage=parallel_sum(eaverage)/totweight;

  eref=eaverage;
  single_write(cout, " setting eref= ", eref, "\n");
  int totconf=parallel_sum(nconfig);
  doublevar eaverage2=0;
  for(int i=0; i < nconfig; i++){
    doublevar effenergy=pts[i].prop.energy(0);
    eaverage2+=(effenergy-eaverage)
               *(effenergy-eaverage);
  }
  
  eaverage2=parallel_sum(eaverage2)/totconf;

  //The variance of the starting distribution.
  doublevar sigmac=sqrt(eaverage2);

  branchcut_start=branch_start_cutoff*sigmac; //start of cutoff region for branching
  branchcut_stop=branch_stop_cutoff*sigmac; //end of cutoff region  

  single_write(cout, " start branch cut at ", branchcut_start, "\n");
  single_write(cout, " stop branch cut at ", branchcut_stop, "\n");


}

//----------------------------------------------------------------------

void Rndmc_method::run(Program_options & options, ostream & output) {
  if(!have_allocated_variables) 
    error("Must generate variables to use Dmc_method::run");
  string logfile=options.runid+".log";
  string logfile2=options.runid+"_unbiased.log";
  string logfile3=options.runid+"_abs_weights.log";
  
  //MB: adding 2 new log files
  //MB: one with unbiased weights according to lubos
  //MB: one with absolute weights to calculate the efficiency 
  if(mpi_info.node==0 ) {
    ofstream logout(logfile.c_str(), ios::app);
    logout << "#-------------------------------------------------\n";
    logout << "#RN-DMC run: biased weights, timestep " << timestep 
           << endl;
    logout << "#-------------------------------------------------\n\n\n";
    logout.close();
    ofstream logout2(logfile2.c_str(), ios::app);
    logout2 << "#-------------------------------------------------\n";
    logout2 << "#RN-DMC run: unbiased weights, timestep " << timestep 
           << endl;
    logout2 << "#-------------------------------------------------\n\n\n";
    logout2.close();
    
    ofstream logout3(logfile3.c_str(), ios::app);
    logout3 << "#-------------------------------------------------\n";
    logout3 << "#RN-DMC run: absolute weights, timestep " << timestep 
           << endl;
    logout3 << "#-------------------------------------------------\n\n\n";
    logout3.close();
  }

  //MB: setting up the property managers
  myprop.setLog(logfile, log_label);
  myprop_unbiased.setLog(logfile2, log_label);
  myprop_absolute.setLog(logfile3, log_label);
  runWithVariables(myprop, mysys, mywfdata, mypseudo, output);
}

//----------------------------------------------------------------------

/*!


 */

//MB: takes the weights inthe array of pts and finds all negatives
//MB: and cancel them out with some positive weights
//MB: only the components of energy and energy are properly adjusted
//MB: for all the other properties you will neeed to fill in

void produce_positive_weights( Array2 <Dmc_point> & pts){
  int nconfig=pts.GetDim(0);
  int npsteps=pts.GetDim(1);
  int nfunc=pts(0,0).prop.weight.GetSize();
  
  for(int p=0; p < npsteps; p++) {
    Array1 <Properties_point> prop_original(nconfig);
    Properties_point new_effective;
    doublevar ave_pos_weight=0.0;
    doublevar sum_neg_weights=0.0;
    
    vector <int>  chosen_neg_walker;
    vector <int>  all_pos_walker;
    for(int walker=0; walker < nconfig; walker++) {
      prop_original(walker)=pts(walker,p).prop;
      if( pts(walker,p).prop.weight(0)>0 ){
	ave_pos_weight+=pts(walker,p).prop.weight(0);
	all_pos_walker.push_back(walker);
      }
      else if(pts(walker,p).prop.weight(0)<0){
	sum_neg_weights+=abs(pts(walker,p).prop.weight(0));
	chosen_neg_walker.push_back(walker);
      }
    }//walker

    new_effective=prop_original(0);
    new_effective.kinetic=0;
    new_effective.potential=0;
    new_effective.nonlocal=0;
    new_effective.weight=0;

    if(all_pos_walker.size())
      ave_pos_weight/=all_pos_walker.size();
    else{
      error("All weight are negative, not able to the unbiased reweighting");
    }
    if(ave_pos_weight*all_pos_walker.size() < sum_neg_weights){
      cout << "node: "<< mpi_info.node <<" sum of possitive weights= "<<ave_pos_weight*all_pos_walker.size()<<" |sum of negative weights|= "<< sum_neg_weights<<endl;
      error("not able to the unbiased reweighting !");
    }
    //cout <<"sum of possitive weights= "<<ave_pos_weight*all_pos_walker.size()<<" |sum of negative weights|= "<< sum_neg_weights<<endl;
    if(sum_neg_weights>0){
      cout << "node: "<<mpi_info.node <<" ave_pos_weight "<<ave_pos_weight<<endl;
      cout << "node: "<<mpi_info.node <<" sum_neg_weights "<<sum_neg_weights<<endl;
    }

    vector <int>  chosen_pos_walker;
    Array1 <int> notselected(nconfig);
    notselected=1;
    doublevar residual_weight=0;
    if(sum_neg_weights > 0.0){
      
      doublevar sum_chosen_pos_weights=0;
      int tries=0;
      while(sum_chosen_pos_weights< sum_neg_weights+ave_pos_weight){
	int walker=int(nconfig*rng.ulec());
	//cout <<"randomly chosen walker: "<<walker<<" with weight "<<pts(walker,p).prop.weight(0)<<" notselected? "<<notselected(walker)<<endl;
	if(pts(walker,p).prop.weight(0)>0 && notselected(walker)){
	  //cout <<"accepted "<< walker <<endl;
	  sum_chosen_pos_weights+=pts(walker,p).prop.weight(0);
	  chosen_pos_walker.push_back(walker);
	  notselected(walker)=0;
	}
	//cout <<" tried "<<tries<<endl; 
	if(tries > 2*nconfig){
	  error("too many tries to match the negative weights");
	}
	tries++;
      }//while
      cout << "node: "<<mpi_info.node << " tries "<<tries<<endl;
     
      residual_weight=sum_chosen_pos_weights-sum_neg_weights;
      //cout <<" residual_weight "<<residual_weight<<endl;

      /*
	will need to get new effective values for these
      Array1 <doublevar> kinetic;
      Array1 <doublevar> potential;
      Array1 <doublevar> nonlocal;
      Array1 <doublevar> weight; //!< averaging weight
      Wf_return wf_val; //!< wavefunction value
      Array1 <dcomplex> z_pol; //!< =exp(i G dot sum x_j)
      Array2 <doublevar> aux_energy;
      Array2 <doublevar> aux_weight;
      Array1 <doublevar> aux_jacobian;
      Array1 <Wf_return> aux_wf_val;
      Array2 <dcomplex> aux_z_pol;
      doublevar gf_weight; //weight of green's function between this and the last point(used in DMC)
      Array1 <doublevar> aux_gf_weight;
      Array1 <Average_return> avgrets;
      */

      
      
	
      for(int k=0;k<chosen_pos_walker.size();k++){
	int w=chosen_pos_walker[k];
	doublevar weight=abs(prop_original(w).weight(0));
       	for(int f=0;f<nfunc;f++){
	  new_effective.kinetic(f)+=weight*prop_original(w).kinetic(f);
	  new_effective.potential(f)+=weight*prop_original(w).potential(f);
	  new_effective.nonlocal(f)+=weight*prop_original(w).nonlocal(f);
	}
      }
      for(int k=0;k<chosen_neg_walker.size();k++){
	int w=chosen_neg_walker[k];
	doublevar weight=abs(prop_original(w).weight(0));
	for(int f=0;f<nfunc;f++){
	  new_effective.kinetic(f)-=weight*prop_original(w).kinetic(f);
	  new_effective.potential(f)-=weight*prop_original(w).potential(f);
	  new_effective.nonlocal(f)-=weight*prop_original(w).nonlocal(f);
	}
      }
      for(int f=0;f<nfunc;f++){
	new_effective.kinetic(f)/=residual_weight;
	new_effective.potential(f)/=residual_weight;
	new_effective.nonlocal(f)/=residual_weight;
	new_effective.weight(f)=residual_weight;
      }
      //cout <<"new_effective.energy(0) "<< new_effective.energy(0)<<endl;
      //cout <<"new_effective.weight(0) "<< new_effective.weight(0)<<endl;

    }//sum_neg_weight > 0.0
    

    int counter=0;
    //coping all not chosen
    for(int walker=0; walker < nconfig; walker++) {
      //cout <<"notselected array "<<notselected(walker)<<endl;
      if(notselected(walker) && prop_original(walker).weight(0) >0 ){
	pts(counter++,p).prop=prop_original(walker);
      }
    }
    //coping 1 new effective positive walker
    if(counter<nconfig)
      pts(counter++,p).prop=new_effective;

    //making the rest all zeros
    while(counter<nconfig){
      pts(counter,p).prop.kinetic=0;
      pts(counter,p).prop.potential=0;
      pts(counter,p).prop.nonlocal=0;
      pts(counter,p).prop.weight=0;
      counter++;
    }
  }//p
}




void Rndmc_method::runWithVariables(Properties_manager & prop, 
				    System * sys, 
				    Wavefunction_data * wfdata,
				    Pseudopotential * pseudo,
				    ostream & output)
{

  allocateIntermediateVariables(sys, wfdata);

  if(!wfdata->supports(laplacian_update))
    error("RNDMC doesn't support all-electron moves..please"
          " change your wave function to use the new Jastrow");

  cout.precision(15);
  output.precision(10);

  
  prop.setSize(wf->nfunc(), nblock, nstep, nconfig, sys, 
	       wfdata);

  //MB: setting things for the myprop_unbiased and myprop_absolute
  myprop_unbiased.setSize(wf->nfunc(), nblock, nstep, nconfig, sys, 
	       wfdata);
  
  myprop_absolute.setSize(wf->nfunc(), nblock, nstep, nconfig, sys, 
	       wfdata);
  restorecheckpoint(readconfig, sys, wfdata, pseudo);

  prop.initializeLog(average_var);

  myprop_unbiased.initializeLog(average_var);

  myprop_absolute.initializeLog(average_var);

  //MB: setting initial sign for each walker and getting the average value;
  doublevar tmp_value=0.0;
  Wf_return wf_val(nwf,1);
  for(int walker=0; walker < nconfig; walker++) {
    pts(walker).config_pos.restorePos(sample);
    wf->updateVal(wfdata, sample);
    wf->getVal(wfdata,0,wf_val);
    pts(walker).sign=wf_val.sign(0);
    tmp_value+=exp(2*wf_val.amp(0,0));
    //cout <<"value "<<exp(wf_val.amp(0,0))<<endl;
  }
  //MB: rescaling the alpha for the average value
  doublevar average_value=sqrt(tmp_value)/nconfig;
  average_value=parallel_sum(average_value)/parallel_sum(1);
  if( mpi_info.node==0 )
    cout <<" average_value: " << average_value<<endl;

  alpha*=average_value;

  if( mpi_info.node==0 )
    cout <<" rescaled alpha: "<<alpha<<endl;
  //MB: seting it in the guidingwf
  guidingwf->set_alpha_for_noderelease(alpha);
  
  



  nhist=1;
  //setting the projection time for auxillary walkers to 1 a.u.
  
  doublevar teff=timestep;

  for(int block=0; block < nblock; block++) {

    int totkilled=0;  
    int totbranch=0;
    int totpoints=0;
    for(int step=0; step < nstep; ) {
      //cout <<"step: "<<step<<endl; 
      int npsteps=min(feedback_interval, nstep-step);

      Dynamics_info dinfo;
      doublevar acsum=0;
      doublevar deltar2=0;
      Array1 <doublevar> epos(3);
      
      doublevar avg_acceptance=0;
      
      doublevar rf_diffusion=0; //diffusion rate without rejection

      Wf_return wf_val(nwf,5);
      for(int walker=0; walker < nconfig; walker++) {
	
        pts(walker).config_pos.restorePos(sample);
        wf->updateLap(wfdata, sample);
	//------Do several steps without branching
        for(int p=0; p < npsteps; p++) {
          pseudo->randomize();
          
          for(int e=0; e< nelectrons; e++) {
            int acc;
            acc=dyngen->sample(e, sample, wf, wfdata, guidingwf,
                               dinfo, timestep);
            
            if(dinfo.accepted) 
              deltar2+=dinfo.diffusion_rate/(nconfig*nelectrons*npsteps);
            
            
            if(dinfo.accepted) { 
              rf_diffusion+=dinfo.diffusion_rate/(nconfig*nelectrons*npsteps);
              
              pts(walker).age(e)=0;
            }
            else { 
              pts(walker).age(e)++;
            }
            avg_acceptance+=dinfo.acceptance/(nconfig*nelectrons*npsteps);
            
            if(acc>0) acsum++;
          }
          
          totpoints++;
          Properties_point pt;
          vector <Tmove> tmov;
          doublevar subtract_out_enwt=0;
          if(tmoves) {  //------------------T-moves
            pt.setSize(nwf);
            wf->getVal(wfdata,0,pt.wf_val);
            sys->calcKinetic(wfdata,sample,wf,pt.kinetic);
            pt.potential=sys->calcLoc(sample);
            pt.weight=1.0; //this gets set later anyway
            pt.count=1;
            pseudo->calcNonlocTmove(wfdata,sys,sample,wf,pt.nonlocal,tmov);
            //cout << "choosing among " <<  tmov.size() << " tmoves " << endl;
            //Now we do the t-move
            doublevar sum=1; 
            for(vector<Tmove>::iterator mov=tmov.begin(); mov!=tmov.end(); mov++) { 
              assert(mov->vxx < 0);
              sum-=timestep*mov->vxx;  
            }
            pt.nonlocal(0)-=(sum-1)/timestep;
            subtract_out_enwt=-(sum-1)/timestep;
            //cout << "sum " << sum <<  " nonlocal " << pt.nonlocal(0) << " ratio " << sum/pt.nonlocal(0) << endl;
            assert(sum >= 0);
            doublevar rand=rng.ulec()*sum;
            sum=1; //reset to choose the move
            if(rand > sum) { 
              for(vector<Tmove>::iterator mov=tmov.begin(); mov!=tmov.end(); mov++) { 
                sum-=timestep*mov->vxx;
                if(rand < sum) { 
                  Array1 <doublevar> epos(3);
                  sample->getElectronPos(mov->e, epos);
                  //cout << "moving electron " << mov->e << " from " << epos(0) << " " << epos(1)
                  //  << " " << epos(2) << " to " << mov->pos(0) << " " << mov->pos(1) 
                  //  << " " << mov->pos(2) << endl;
                  sample->setElectronPos(mov->e,mov->pos);
                  break;
                }
              }
            }
            //wf->updateLap(wfdata, sample);
          } ///---------------------------------done with the T-moves
          else { 
	    //MB: this is where all the properties at the new step are evaluated
            mygather.gatherData(pt, pseudo, sys, wfdata, wf, 
                                sample, guidingwf);
	    	    
	  }
          
          Dmc_history new_hist;
          new_hist.main_en=pts(walker).prop.energy(0);
          pts(walker).past_energies.push_front(new_hist);
          deque<Dmc_history> & past(pts(walker).past_energies);
          if(past.size() > nhist) 
            past.erase(past.begin()+nhist, past.end());
          
	  //copy pt to pts(walker).prop
          pts(walker).prop=pt;

	  //MB: this is how we track the sign
          //MB: pts(walker).sign is before the step and the pts(walker).prop.sign(0) is after
          //MB: if they change weight will change the sign as well
	  if(pts(walker).sign!=pts(walker).prop.wf_val.sign(0)){
	    cout <<"node: "<<mpi_info.node <<" walker "<<walker<<" changed the sign "<<endl;
	    //cout <<" pts(walker).weight "<<pts(walker).weight<<" pts(walker).prop.weight(0) "<<pts(walker).prop.weight(0)<<
	    // " pts(walker).sign "<<pts(walker).sign<<" pts(walker).prop.sign "<<pts(walker).prop.sign<<endl;
	    pts(walker).sign*=-1;
	    pts(walker).weight*=-1;
	  }

	  //MB: this is part of normal DMC run
          //MB: because we do not branch adjusting the weights according 
          //MB: to etrial makes no sence, but disscuss this with lubos, I am leaving this out

	  //pts(walker).weight*=getWeight(pts(walker),teff,etrial);

          if(pts(walker).ignore_walker) {
            pts(walker).ignore_walker=0;
            pts(walker).weight=1;
            pts(walker).prop.count=0;
          }

	  //MB: this is new, since in DMC pts(walker).prop.weight(0)=1, this was not needed
          //MB: because alpha>0 it is no longer true
	  pts(walker).weight*=pts(walker).prop.weight(0);
	  //keep total weight in pts(walker).weight

	  //finally everything is put to pts(walker).prop for averaging
          pts(walker).prop.weight=pts(walker).weight;
          
	  pts(walker).prop.avgrets.Resize(1,average_var.GetDim(0));
          for(int i=0; i< average_var.GetDim(0); i++) { 
            average_var(i)->evaluate(wfdata, wf, sys, sample, pts(walker).prop.avgrets(0,i));
          }
          
	  //MB: after evaluating the properties all are stored in the prop
          prop.insertPoint(step+p, walker, pts(walker).prop);

       	  //MB: copy pts for the calculation of unbiased walker
	  pts_unbiased(walker, p).prop=pts(walker).prop;

	  //MB: copy another pts for absolute weights
	  Dmc_point  pts_absolute;
	  pts_absolute=pts(walker);
	  pts_absolute.prop.weight=fabs(pts_absolute.prop.weight(0));
	  myprop_absolute.insertPoint(step+p, walker, pts_absolute.prop);
	  
	  
          for(int i=0; i< densplt.GetDim(0); i++)
            densplt(i)->accumulate(sample,pts(walker).prop.weight(0));
	  for(int i=0; i< nldensplt.GetDim(0); i++)
	    nldensplt(i)->accumulate(sample,pts(walker).prop.weight(0),
				     wfdata,wf);
        }//p
        
        pts(walker).config_pos.savePos(sample);
      }//walker
      //---Finished moving all walkers
      
      /*
      just debug printout
      cout <<"before produce_positive_weights"<<endl;
      for(int p=0; p < npsteps; p++) {
	for(int walker=0; walker < nconfig; walker++) {
	  cout <<"energy: "<<pts_unbiased(walker,p).prop.energy(0)<<" weight "<<pts_unbiased(walker,p).prop.weight(0)<<endl;
	}
	cout <<endl;
      }
      */
      
      //MB: this makes unbiased positive  weights according the to lubos scheme
      //MB: current drawback is that happens on each node separately so meake sure you have enough walkers per node!
      produce_positive_weights(pts_unbiased);
      
      //cout <<"after produce_positive_weights"<<endl;
      //MB: store unbiased prop for averaging 
      for(int p=0; p < npsteps; p++) {
	for(int walker=0; walker < nconfig; walker++) {
	  //cout <<"energy: "<<pts_unbiased(walker,p).prop.energy(0)<<" weight "<<pts_unbiased(walker,p).prop.weight(0)<<endl;
	  myprop_unbiased.insertPoint(step+p, walker, pts_unbiased(walker,p).prop);
	 }
	//cout <<endl;
       }
      
      doublevar accept_ratio=acsum/(nconfig*nelectrons*npsteps);
      teff=timestep*accept_ratio; //deltar2/rf_diffusion; 

      updateEtrial(feedback);

      step+=npsteps;

      //MB: turnning off branching 
      int nkilled=0; //calcBranch();
      totkilled+=nkilled;
      totbranch+=nkilled;
    }//step
    //cout <<"Finished block "<<endl;
    ///----Finished block
    
    if(!low_io || block==nblock-1) {
      savecheckpoint(storeconfig,sample);
      for(int i=0; i< densplt.GetDim(0); i++)
        densplt(i)->write();
      for(int i=0; i< nldensplt.GetDim(0); i++)
        nldensplt(i)->write(log_label);
    }

    
    prop.endBlock();
    //MB: averaging over block for the unbiased and absolute weights
    myprop_unbiased.endBlock();
    myprop_absolute.endBlock();

    totbranch=parallel_sum(totbranch);
    totkilled=parallel_sum(totkilled);
    totpoints=parallel_sum(totpoints);

    Properties_final_average finavg;
    Properties_final_average finavg2;
    Properties_final_average finavg3;

    prop.getFinal(finavg);
    myprop_unbiased.getFinal(finavg2);
    myprop_absolute.getFinal(finavg3);


    Properties_block lastblock;
    Properties_block lastblock2;
    Properties_block lastblock3;

    prop.getLastBlock(lastblock);
    myprop_unbiased.getLastBlock(lastblock2);
    myprop_absolute.getLastBlock(lastblock3);
    

    //MB: finalavg: total average value over all blocks up to this time
    //MB: lastblock: average value for last block

    //MB: original DMC had eref as this:
    eref=finavg2.avg(Properties_types::total_energy,0);
    //MB: not sure this makes sence for the release node
    //MB: maybe it should be only over the last block
    //eref=lastblock2.avg(Properties_types::total_energy,0);

    updateEtrial(feedback);
   
    //MB: this is lubos's effectivity, which is ratio of averaged signed
    //MB: weights to averaged absolute weights 
    //MB: bellow is the value for each block

    doublevar effectivity;
    doublevar weight_biased=lastblock.avg(Properties_types::weight,0);
    doublevar weight_abs=lastblock3.avg(Properties_types::weight,0);

    effectivity=weight_biased/weight_abs;
        
    doublevar maxage=0;
    doublevar avgage=0;
    for(int w=0;w < nconfig; w++) {
      for(int e=0; e< nelectrons; e++) { 
        if(maxage<pts(w).age(e)) maxage=pts(w).age(e);
        avgage+=pts(w).age(e);
      }
    }
    avgage/=(nconfig*nelectrons);

    if(output) {
      //cout << "Block " << block 
      //       << " nconfig " << totconfig
      //       << " etrial " << etrial << endl;
      if(global_options::rappture ) { 
	    ofstream rapout("rappture.log");
        rapout << "=RAPPTURE-PROGRESS=>" << int(100.0*doublevar(block+1)/doublevar(nblock))
               << "  Diffusion Monte Carlo" << endl;
        cout << "=RAPPTURE-PROGRESS=>" << int(100.0*doublevar(block+1)/doublevar(nblock))
             << "  Diffusion Monte Carlo" << endl;
        rapout.close();
      }
      output << "***" << endl;
      output << "Block " << block 
             << " etrial " << etrial << endl;
      output << "maximum age " << maxage 
	     << " average age " << avgage << endl;
      dyngen->showStats(output);

      output<<"\n ---- biased weights ---- \n";
      prop.printBlockSummary(output);
      output<<"\n ---- unbiased weights ---- \n";
      myprop_unbiased.printBlockSummary(output);
     
      output<<"\n ---- absolute weights ---- \n";
      myprop_absolute.printBlockSummary(output);
      output<<"\n effectivity  "<<effectivity<<" signed weight "<<weight_biased<<" absolute weight "<<weight_abs<<endl;

      output << "Branched "
	     << totbranch << " times.  So a branch every " 
	     << doublevar(totpoints)/doublevar(totbranch)
	     << " steps " << endl;
    }

    dyngen->resetStats();

  }//block
  
  if(output) {
    output << "\n ----------Finished RN-DMC------------\n\n";
    output<<"\n ---- biased weights ---- \n";
    prop.printSummary(output,average_var);
    output<<"\n ---- unbiased weights ---- \n";
    myprop_unbiased.printSummary(output,average_var); 
    output<<"\n ---- absolute weights ---- \n";
    myprop_absolute.printSummary(output,average_var);
  }
  wfdata->clearObserver();
  deallocateIntermediateVariables();
}


//----------------------------------------------------------------------


void Rndmc_method::savecheckpoint(string & filename,                     
                                 Sample_point * config) {
  if(filename=="") return;
  ofstream checkfile(filename.c_str());
  if(!checkfile) error("Couldn't open", filename );
  checkfile.precision(15);
  
  long int is1, is2;
  rng.getseed(is1, is2);
  checkfile << "RANDNUM " << is1 << "  " << is2 << endl;

  checkfile.precision(15);
  for(int i=0; i< nconfig; i++) { 
    Dmc_point & mypt(pts(i));
    checkfile << "SAMPLE_POINT { \n";
    mypt.config_pos.restorePos(config);
    write_config(checkfile, config);
    checkfile << "   DMC { \n";
    checkfile << "DMCWEIGHT " << mypt.weight << endl;
    checkfile << "VALEN " << nwf << endl; 
    for(int w=0; w< nwf; w++) {
      checkfile << mypt.prop.wf_val.phase(w,0) << "  "
		<< mypt.prop.wf_val.amp(w,0) << "  "
		<< mypt.prop.energy(w)
		<< endl;
    }

    checkfile << "   } \n";
    checkfile << "}\n\n";
  }

  checkfile.close();
}


//----------------------------------------------------------------------



void Rndmc_method::restorecheckpoint(string & filename, System * sys,
                                    Wavefunction_data * wfdata,
                                    Pseudopotential * pseudo) {


  ifstream checkfile(filename.c_str());
  if(!checkfile) 
    error("Couldn't open ", filename);
  long int is1, is2;
  string dummy;
  checkfile >> dummy;
  if(dummy != "RANDNUM") error("Expected RANDNUM in checkfile");
  checkfile >> is1 >> is2;
  rng.seed(is1, is2);

  Array1 <Wf_return > value_temp(nconfig);
  Array2 <doublevar> energy_temp(nconfig, nwf);


  int ncread=0; //Number of configs read
  int nwread=0; //number of weights read
  while(checkfile >> dummy && 
	( ncread < nconfig && nwread < nconfig) ) {

    
    if(read_config(dummy, checkfile, sample)) { 
      pts(ncread++).config_pos.savePos(sample);
      
    }

    if(dummy=="DMC") {
      checkfile >> dummy;
      if(dummy != "{") error("Need a { after DMC");
      checkfile >> dummy >> pts(nwread).weight;
      if(dummy != "DMCWEIGHT") {
        error("expected DMCWEIGHT, got ", dummy);
      }
      int nwf_temp;
      checkfile >> dummy >> nwf_temp;
      if(dummy != "VALEN") {
        error("expected VALEN, got ", dummy);
      }
      if(nwf_temp != nwf) {
        error("Wrong number of wavefunctions in the checkpoint file");
      }
      
      //Retrieve the old values and energies from the file
      value_temp(nwread).Resize(nwf, 2);
      
      for(int w=0; w< nwf; w++) {
        checkfile >> value_temp(nwread).phase(w,0) 
		  >> value_temp(nwread).amp(w,0) 
                  >> energy_temp(nwread,w);
      }
      nwread++;

    }

  }

  //cout << "ncread " << ncread << "  nwread " << nwread << endl;
  if(nconfig!=ncread) { 
    error("nconfig doesn't match the number of walkers in the config file");
  }

 
  for(int walker=0; walker < nconfig; walker++) {
    pts(walker).config_pos.restorePos(sample);
    mygather.gatherData(pts(walker).prop, pseudo, sys,
                        wfdata, wf, sample,
                        guidingwf);
  }
  
  find_cutoffs();

  updateEtrial(start_feedback);

  if(do_cdmc) { 
    if(ncread!=nwread) {
      cout << "WARNING! do_cdmc and ncread!=nwread " << endl;
    }
    cdmcReWeight(energy_temp, value_temp);
  }

    
}


//----------------------------------------------------------------------

void Rndmc_method::cdmcReWeight(Array2 <doublevar> & energy_temp, 
                              Array1 < Wf_return > & value_temp
                              ) {
  //The following is for C-DMC(from Jeff)  It shouldn't
  //do anything unless the atomic coordinates have moved
  //It's commented out for the moment from the rewrite.
  //It's viewed as experimental code, so watch out!

  //Get the energies for the old ionic configurations and 
  //this configuration
  /*
  Array1 <doublevar> effenergy_temp(nconfig, 0.0);
  Array1 <doublevar> effoldenergy(nconfig,0.0);

  for(int i=0; i< nconfig; i++) {
    for(int w=0; w< nwf; w++) {
      effenergy_temp(i)+=(energy_temp(i,w)+offset(w))
        *guidingwf->getOperatorWeight(value_temp(i),w);

      doublevar olden=trace[0][i].energy(w);
      effoldenergy(i)+=(olden+offset(w))
        *guidingwf->getOperatorWeight(trace[0][i].wf_val, w);
    }
  }

  

  doublevar average_temp, variance_temp;

  average(0, nconfig, effenergy_temp,
          dmcweight, average_temp, variance_temp, 1);

  doublevar sigma_temp;
  if(variance_temp >0) {
    sigma_temp=sqrt(variance_temp);
  }
  else {
    sigma_temp=0;
    error("negative variance when reading in weights");
  }

  debug_write(cout, "average_temp ", average_temp, "\n");
  debug_write(cout, "sigma_temp ", sigma_temp, "\n");

  //Reweight and check whether we had a sign flip
  //(overlap will be either positive or negative

  if(do_cdmc) {
    doublevar norm1=0, norm2=0;
    for(int i=0; i< nconfig; i++) {
      //norm1+=exp(prop.trace(0,i).wf_val(0,1)*2.0);
      norm1+=exp(trace[0][i].wf_val.amp(0,0)*2.0);
      
      norm2+=exp(2.0*value_temp(i).amp(0,0));
    }
    
    int totconfig=parallel_sum(nconfig);
    norm1=parallel_sum(norm1)/totconfig;
    norm2=parallel_sum(norm2)/totconfig;
    
    //Removed the normalization(norm2/norm1), because it causes problems
    //for bigger systems
    
    doublevar overlap=0;
    for(int i=0; i < nconfig; i++) {
      doublevar ratio;
      ratio=guidingwf->getTrialRatio(trace[0][i].wf_val, value_temp(i));
      overlap+=ratio/nconfig; //check if there's an overall sign change
      
      //cout << "walker " << i << " new val " << trace[0][i].wf_val(0,1) 
      //     << " old one " << value_temp(i)(0,1) << endl;

      dmcweight(i)*= ratio*ratio//  *norm2/norm1
        *exp(-(effoldenergy(i)-effenergy_temp(i))*timestep/2);
    }
    
    single_write(cout, "overlap between old and new ", overlap);
    doublevar average_old=0.0, diff_sigma=0.0;
     
    for(int w=0; w< nconfig; w++)  {
       average_old+=effoldenergy(w)/nconfig;
       diff_sigma+=(effoldenergy(w)-effenergy_temp(w))*(effoldenergy(w)-effenergy_temp(w))/nconfig;
    }
    //cout << "difference " << average_temp-average_old <<"  +/-  " << diff_sigma/nconfig << endl;

    //We ignore walkers that either a)cross a node, or 
    //b) are outside of 10 sigmas(for the first move)
    //ofstream diffout("diff_en.dat");
    int ncross=0, nsigma=0;
    for(int i=0; i< nconfig; i++) {
      doublevar ratio;
      ratio=guidingwf->getTrialRatio(trace[0][i].wf_val, value_temp(i));
      //cout << "ratio " << i << "   " << ratio << endl;
      if(ratio*overlap < 0) {
        debug_write(cout, "crossed node ", i , "\n");
        ncross++;
        trace[0][i].count=0;
        dmcweight(i)=1;
        ignore_walker(i)=1;
      }
      else if(fabs(effenergy_temp(i)-average_temp) > 10*sigma_temp) {
        debug_write(cout, "walker out of 10 sigmas ", i, "\n");
	
        nsigma++;
        trace[0][i].count=0;
        dmcweight(i)=1;
        ignore_walker(i)=1;
      }
      // diffout << i << "   " << effenergy_temp(i)-effoldenergy(i) << endl;
    }
    //diffout.close();

    single_write(cout, "nsigma ", parallel_sum(nsigma));
    single_write(cout, "  ncross  ", parallel_sum(ncross), "\n");

  }

  */
  

}


//----------------------------------------------------------------------

void Rndmc_method::updateEtrial(doublevar feed) {
  
  doublevar totweight=0;
  for(int walker=0; walker < nconfig; walker++)
    totweight+=pts(walker).weight;

  etrial=eref-feed*log(totweight/double(nconfig));

#ifdef USE_MPI
  doublevar mpitot=0;
  MPI_Allreduce(&totweight,&mpitot, 1,
                MPI_DOUBLE, MPI_SUM, MPI_Comm_grp);
  etrial=eref-feed*log(mpitot/double(mpi_info.nprocs*nconfig));
#endif
  
  //cout << "etrial " << etrial << " total weight " << totweight << endl;

}

//----------------------------------------------------------------------

doublevar Rndmc_method::getWeight(Dmc_point & pt,
                                doublevar teff, doublevar etr) {
  doublevar teffac=teff/2.0;

  doublevar effenergy=0, effoldenergy=0;

  effenergy=pt.prop.energy(0);
  effoldenergy=pt.past_energies[0].main_en;

  doublevar fbet=max(etr-effenergy, etr-effoldenergy);

  if(fbet > branchcut_stop) {
    teffac=0;
  }
  else if(fbet > branchcut_start) {
    teffac=teffac*(1.-(fbet-branchcut_start)
                   /(branchcut_stop-branchcut_start));
  }

  doublevar return_weight;
  //if(tmoves) return_weight=exp(teffac*2.0*(etr-effenergy));
  //else 
  return_weight=exp(teffac*(etr*2-effenergy-effoldenergy));

  return return_weight;
} 
                                
//----------------------------------------------------------------------


//----------------------------------------------------------------------
struct Queue_element { 
  int from_node;
  int to_node;
  Queue_element() { } 
  Queue_element(int from, int to) { from_node=from; to_node=to; } 
};

struct Walker_sort { 
  int node;
  int index; //on the node
  int branch; //how many copies to make
  doublevar weight; 
};



int Rndmc_method::calcBranch() { 
  int totwalkers=mpi_info.nprocs*nconfig;
  Array1 <doublevar> weights(totwalkers);
  Array1 <doublevar> my_weights(nconfig);
  
  for(int walker=0; walker < nconfig; walker++)
    my_weights(walker)=pts(walker).weight;
#ifdef USE_MPI
  MPI_Allgather(my_weights.v,nconfig, MPI_DOUBLE, weights.v,nconfig,MPI_DOUBLE, MPI_Comm_grp);
#else
  weights=my_weights;
#endif
  Array1 <int> my_branch(nconfig);
  Array1 <int> nwalkers(mpi_info.nprocs);
  nwalkers=0;
  int root=0;
  if(mpi_info.node==root) {  //this if/else clause may be refactored out

    Array1 <int> branch(totwalkers);
    //----Find which walkers branch/die
    //we do it on one node since otherwise we'll have different random numbers!
    //we'll assign the weight for each copy that will be produced
    //this is the core of the branching algorithm..
    //my homegrown algo, based on Umrigar, Nightingale, and Runge
    branch=-1;
    const doublevar split_threshold=1.8;
    for(int w=0; w< totwalkers; w++) { 
      if(weights(w) > split_threshold && branch(w)==-1) { 
        //find branching partner
        doublevar smallestwt=100;
        int smallest=-1;
        for(int w2=0; w2 < totwalkers; w2++) { 
          if(branch(w2)==-1 && w2!= w && weights(w2) < smallestwt) { 
            smallest=w2;
            smallestwt=weights(w2);
          }
        }
        if(smallest != -1) { 
          doublevar weight1=weights(w)/(weights(w)+weights(smallest));
          if(weight1+rng.ulec() >= 1.0) { 
            branch(w)=2;
            branch(smallest)=0;
            weights(w)+=weights(smallest);
            weights(w)/=2.0;
          }
          else { 
            branch(w)=0;
            branch(smallest)=2;
            weights(smallest)+=weights(w);
            weights(smallest)/=2.0;
          }
        }
      }
    }
    for(int w=0; w< totwalkers; w++) {
      if(branch(w)==-1) branch(w)=1;
    }
    //----end homegrown algo
    //count how many walkers each node will have
    //without balancing
    int walk=0;
    for(int n=0; n< mpi_info.nprocs; n++) { 
      for(int i=0; i< nconfig; i++) {
        nwalkers(n)+=branch(walk);
        walk++;
      }
      //cout << "nwalkers " << n << " " << nwalkers(n) << endl;
    }
    //now send nwalkers and which to branch to all processors
    for(int i=0; i< nconfig; i++) { 
      my_branch(i)=branch(i);
      my_weights(i)=weights(i);
    }
#ifdef USE_MPI
    MPI_Bcast(nwalkers.v, mpi_info.nprocs, MPI_INT, mpi_info.node, MPI_Comm_grp);
    for(int i=1; i< mpi_info.nprocs; i++) {
      MPI_Send(branch.v+i*nconfig,nconfig,MPI_INT,i,0,MPI_Comm_grp);
      MPI_Send(weights.v+i*nconfig, nconfig, MPI_DOUBLE, i,0,MPI_Comm_grp);
    }
#endif
               
  }
  else { 
#ifdef USE_MPI
    MPI_Bcast(nwalkers.v, mpi_info.nprocs, MPI_INT, root, MPI_Comm_grp);
    MPI_Status status;
    MPI_Recv(my_branch.v,nconfig, MPI_INT,root,0,MPI_Comm_grp, &status);
    MPI_Recv(my_weights.v,nconfig, MPI_DOUBLE,root,0,MPI_Comm_grp, &status);
#endif	
  }
  //--end if/else clause

  for(int i=0; i< nconfig; i++) { 
    pts(i).weight=my_weights(i);
  }
  
  //Now we all have my_branch and nwalkers..we need to figure out who 
  //needs to send walkers to whom--after this, nwalkers should be a flat array equal to 
  //nconfig(so don't try to use it for anything useful afterwards)
  vector <Queue_element> send_queue;
  int curr_needs_walker=0;
  int nnwalkers=nwalkers(mpi_info.node); //remember how many total we should have
  for(int i=0; i< mpi_info.nprocs; i++) { 
    while(nwalkers(i) > nconfig) {
      if(nwalkers(curr_needs_walker) < nconfig) { 
        nwalkers(curr_needs_walker)++;
        nwalkers(i)--;
        send_queue.push_back(Queue_element(i,curr_needs_walker));
        //cout << mpi_info.node << ":nwalkers " << nwalkers(i) << "  " << nwalkers(curr_needs_walker) << endl;
        //cout << mpi_info.node << ":send " << i << "  " << curr_needs_walker << endl;
      }
      else { 
        curr_needs_walker++;
      }
    }
  }
  
  for(int i=0; i< mpi_info.nprocs; i++) assert(nwalkers(i)==nconfig);
  int killsize=0;
  for(int i=0; i< nconfig; i++) {
    //cout << mpi_info.node << ":branch " << i << "  " << my_branch(i) << " weight " << pts(i).weight <<  endl;
    if(my_branch(i)==0) killsize++;
  }
  //cout << mpi_info.node << ": send queue= " << send_queue.size() << endl;
  //now do branching for the walkers that we get to keep
  Array1 <Dmc_point> savepts=pts;
  int curr=0; //what walker we're currently copying from
  int curr_copy=0; //what walker we're currently copying to
  while(curr_copy < min(nnwalkers,nconfig)) { 
    if(my_branch(curr)>0) { 
      //cout << mpi_info.node << ": copying " << curr << " to " << curr_copy << " branch " << my_branch(curr) << endl;
      my_branch(curr)--;
      pts(curr_copy)=savepts(curr);
      //pts(curr_copy).weight=1;
      curr_copy++;
    }
    else curr++;
  }
  
  //Finally, send or receive spillover walkers 
  if(nnwalkers > nconfig) { 
    vector<Queue_element>::iterator queue_pos=send_queue.begin();
    while(curr < nconfig) { 
      if(my_branch(curr) > 0) { 
        my_branch(curr)--;
        while(queue_pos->from_node != mpi_info.node) { 
          queue_pos++;
        }
        //cout << mpi_info.node << ":curr " << curr << " my_branch " << my_branch(curr) << endl;
        //cout << mpi_info.node << ":sending " << queue_pos->from_node << " to " << queue_pos->to_node << endl;
        savepts(curr).mpiSend(queue_pos->to_node);
        queue_pos++;
      }
      else curr++;
    }
    
  }
  else { //if nnwalkers == nconfig, then this will just get skipped immediately
    vector <Queue_element>::iterator queue_pos=send_queue.begin();
    while(curr_copy < nconfig) { 
      while(queue_pos->to_node != mpi_info.node) queue_pos++;
      //cout << mpi_info.node <<":receiving from " << queue_pos->from_node << " to " << curr_copy << endl;
      pts(curr_copy).mpiReceive(queue_pos->from_node);
      //pts(curr_copy).weight=1;
      curr_copy++;
      queue_pos++;
    }
  }
  return killsize;
  //exit(0);
}
//----------------------------------------------------------------------

