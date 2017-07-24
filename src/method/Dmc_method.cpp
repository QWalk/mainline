/*
 
Copyright (C) 2007 Lucas K. Wagner
addition of the PURE DMC: Michal Bajdich and Fernando A. Reboredo

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



#include "Dmc_method.h"
#include "qmc_io.h"
#include "ulec.h"
#include "Program_options.h"
#include "average.h"
#include "Generate_sample.h"
#include <algorithm>
#include <ctime>
#include <cstdio>
void Dmc_method::read(vector <string> words,
                      unsigned int & pos,
                      Program_options & options)
{
  ndim=3;

  have_read_options=1;

  //vector <string> offsetwords;

  //required options
  if(!readvalue(words, pos=0, timestep, "TIMESTEP"))
    error("Need TIMESTEP in METHOD section");

  //optional options

  if(!readvalue(words, pos=0, nblock, "NBLOCK"))
    nblock=100;
 
  int target_config;
  if(!readvalue(words,pos=0,target_config,"TARGET_CONFIG"))
    target_config=2048;

  if(!readvalue(words, pos=0, nconfig, "NCONFIG"))
    nconfig=max(target_config/mpi_info.nprocs,1);

  if(!readvalue(words, pos=0, nstep, "NSTEP"))
    nstep=max(int(1.0/timestep+0.5),1);
  
  if(!readvalue(words, pos=0, readconfig, "READCONFIG"))
    readconfig=options.runid+".config";

  if(!readvalue(words, pos=0, eref, "EREF"))
    eref=0.0;

  if(!readvalue(words,pos=0, max_poss_weight, "MAX_POSS_WEIGHT")) 
    max_poss_weight=7.0;
  if(haskeyword(words, pos=0, "TMOVES")) tmoves=1;
  else tmoves=0;
  if(haskeyword(words, pos=0, "TMOVESSC")) tmoves_sizeconsistent=1;
  else tmoves_sizeconsistent=0;

  all_electron_moves=false;
  if(haskeyword(words,pos=0,"ALL_ELECTRON")) all_electron_moves=true;


  readvalue(words,pos=0,save_trace,"SAVE_TRACE");

  if(!readvalue(words, pos=0, nhist, "CORR_HIST")) 
    nhist=-1;

  if(haskeyword(words, pos=0, "PURE_DMC")) { 
    pure_dmc=1;
    if( nhist < 1 )
      error("CORR_HIST must be > 1 when PURE_DMC is on");
  }
  else pure_dmc=0;

  if(!readvalue(words, pos=0, storeconfig, "STORECONFIG"))
    storeconfig=options.runid+".config";

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
  
  allocate(dynamics_words, dyngen);
  dyngen->enforceNodes(1);

  //MB: forwark walking lengths 
  fw_length.Resize(0);
  fw_length=0;
  vector <string> fw_words;
  max_fw_length=0;
  if(readsection(words, pos=0, fw_words, "fw_length")){
    fw_length.Resize(fw_words.size());
    for(int i=0;i<fw_words.size();i++){
      fw_length(i)=atoi(fw_words[i].c_str());
      if(fw_length(i) > max_fw_length)
        max_fw_length=fw_length(i);
    }
  }
  //cout <<" Maximum forward walking history is "<<max_fw_length<<" steps"<<endl;

}

//----------------------------------------------------------------------

int Dmc_method::generateVariables(Program_options & options) {

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



int Dmc_method::allocateIntermediateVariables(System * sys,
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

  guidingwf=new Primary;

  pts.Resize(nconfig);
  for(int i=0; i < nconfig; i++) {
    pts(i).age.Resize(nelectrons);
    pts(i).age=0;
  }
  
  average_var.Resize(avg_words.size());
  average_var=NULL;
  for(int i=0; i< average_var.GetDim(0); i++) { 
    allocate(avg_words[i], sys, wfdata, average_var(i));
  }
  

  return 1;
}

//----------------------------------------------------------------------

int Dmc_method::showinfo(ostream & os)
{

  if(have_allocated_variables) {
    mysys->showinfo(os);
    mypseudo->showinfo(os);
    mywfdata->showinfo(os);
  }
  os << "###########################################################\n";
  os << "Diffusion Monte Carlo:\n";
  os << "Number of processors " <<           mpi_info.nprocs << endl;
  os << "Blocks: " <<                        nblock    << endl;
  os << "Steps per block: " <<               nstep     << endl;
  os << "Timestep: " <<                      timestep  << endl;
  if(tmoves) 
    os << "T-moves turned on" << endl;
  if(tmoves_sizeconsistent)
    os << "Size-consistent T-moves turned on" << endl;

  string indent="  ";

  if(max_fw_length){
    os << "Forward walking averaging over lengths "; 
    for(int s=0;s<fw_length.GetSize();s++)
      os<<fw_length(s)<<"  ";
    os<< endl;
  }
  dyngen->showinfo(indent, os);

  os << "###########################################################" << endl;
  return 1;
}

//----------------------------------------------------------------------

void Dmc_method::find_cutoffs() {
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

void Dmc_method::run(Program_options & options, ostream & output) {
  if(!have_allocated_variables) 
    error("Must generate variables to use Dmc_method::run");
  string logfile=options.runid+".log";
  
  if(mpi_info.node==0 ) {
    ofstream logout(logfile.c_str(), ios::app);
    logout << "#-------------------------------------------------\n";
    logout << "#DMC run: timestep " << timestep 
           << endl;
    logout << "#-------------------------------------------------\n\n\n";
    logout.close();
  }

  myprop.setLog(logfile, log_label);
  runWithVariables(myprop, mysys, mywfdata, mypseudo,output);
}

//----------------------------------------------------------------------

/*!

 */
void Dmc_method::runWithVariables(Properties_manager & prop, 
                                  System * sys, 
                                  Wavefunction_data * wfdata,
                                  Pseudopotential * pseudo,
                                  ostream & output)
{

  allocateIntermediateVariables(sys, wfdata);
  if(!wfdata->supports(laplacian_update))
    error("DMC doesn't support all-electron moves..please"
          " change your wave function to use the new Jastrow");

  cout.precision(15);
  output.precision(10);

  
  prop.setSize(wf->nfunc(), nblock, nstep, nconfig, sys, 
	       wfdata);

  restorecheckpoint(readconfig, sys, wfdata, pseudo);
  prop.initializeLog(average_var);

  //MB: new properties manager for forward walking (one per each length)
  Array1 <Properties_manager> prop_fw;
  prop_fw.Resize(fw_length.GetSize());
  if(max_fw_length){
    for(int s=0;s<fw_length.GetSize();s++){
      string logfile, label_temp;
      prop.getLog(logfile, label_temp);
      char strbuff[40];
      sprintf(strbuff, "%d", fw_length(s));
      label_temp+="_fw";
      label_temp+=strbuff;
      prop_fw(s).setLog(logfile, label_temp);
      prop_fw(s).setSize(wf->nfunc(), nblock, nstep, nconfig, sys, 
		    wfdata);
      prop_fw(s).initializeLog(average_var);
    }
  }


  if(nhist==-1)
    nhist=1;
  
  doublevar teff=timestep;
  for(int block=0; block < nblock; block++) {

    int totkilled=0;  
    int totbranch=0;
    int totpoints=0;
    for(int step=0; step < nstep; ) {
      int npsteps=min(feedback_interval, nstep-step);

      Dynamics_info dinfo;
      doublevar acsum=0;
      doublevar deltar2=0;
      Array1 <doublevar> epos(3);
      
      doublevar avg_acceptance=0;
      
      for(int walker=0; walker < nconfig; walker++) {
        pts(walker).config_pos.restorePos(sample);
        wf->updateLap(wfdata, sample);
	//------Do several steps without branching
        for(int p=0; p < npsteps; p++) {
          pseudo->randomize();
          if(all_electron_moves) move_all_electron(wfdata,wf,sample,guidingwf,pts(walker),acsum);
          else move_electron_by_electron(wfdata,wf,sample,guidingwf,pts(walker),acsum);
          totpoints++;
          Properties_point pt;
          if(tmoves or tmoves_sizeconsistent) {  //------------------T-moves
            doTmove(pt,pseudo,sys,wfdata,wf,sample,guidingwf);
          } ///---------------------------------done with the T-moves
          else {
            mygather.gatherData(pt, pseudo, sys, wfdata, wf, 
                                sample, guidingwf);
          }
          Dmc_history new_hist;
          new_hist.main_en=pts(walker).prop.energy(0);
          pts(walker).past_energies.push_front(new_hist);
          deque<Dmc_history> & past(pts(walker).past_energies);
          if(past.size() > nhist) 
            past.erase(past.begin()+nhist, past.end());
          
          pts(walker).prop=pt;
          if(!pure_dmc) { 
            pts(walker).weight*=getWeight(pts(walker),teff,etrial);
            //Introduce a small bias to avoid instability.
            if(pts(walker).weight>max_poss_weight) pts(walker).weight=max_poss_weight;
          }
          else
            pts(walker).weight=getWeightPURE_DMC(pts(walker),teff,etrial);
          
          if(pts(walker).ignore_walker) {
            pts(walker).ignore_walker=0;
            pts(walker).weight=1;
            pts(walker).prop.count=0;
          }
          pts(walker).prop.weight=pts(walker).weight;
          //This is somewhat inaccurate..will need to change it later
          //For the moment, the autocorrelation will be slightly
          //underestimated
          pts(walker).prop.parent=walker;
          pts(walker).prop.nchildren=1;
          pts(walker).prop.children(0)=walker;
          pts(walker).prop.avgrets.Resize(1,average_var.GetDim(0));
          for(int i=0; i< average_var.GetDim(0); i++) { 
            average_var(i)->randomize(wfdata,wf,sys,sample);
            average_var(i)->evaluate(wfdata, wf, sys, pseudo, sample, pts(walker).prop, pts(walker).prop.avgrets(0,i));
          }
          prop.insertPoint(step+p, walker, pts(walker).prop);
          for(int i=0; i< densplt.GetDim(0); i++)
            densplt(i)->accumulate(sample,pts(walker).prop.weight(0));
          for(int i=0; i< nldensplt.GetDim(0); i++)
            nldensplt(i)->accumulate(sample,pts(walker).prop.weight(0),
                                     wfdata,wf);
          
          
          //MB: making the history of prop.avgrets for forward walking
          if(max_fw_length){
            forwardWalking(walker, step+p,prop_fw);
          }//if FW
          
        }

        pts(walker).config_pos.savePos(sample);
      }
      //---Finished moving all walkers

      doublevar accept_ratio=acsum/(nconfig*nelectrons*npsteps);
      teff=timestep*accept_ratio; //deltar2/rf_diffusion; 
      if(output) output << "Acceptance " << accept_ratio << endl;
      
      updateEtrial(feedback);

      step+=npsteps;

      int nkilled;
      if(!pure_dmc)
        nkilled=calcBranch();
      else
        nkilled=0;
      
      totkilled+=nkilled;
      totbranch+=nkilled;
    }

    ///----Finished block
    
    if(!low_io || block==nblock-1) {
      savecheckpoint(storeconfig,sample);
      for(int i=0; i< densplt.GetDim(0); i++)
        densplt(i)->write();
      for(int i=0; i< nldensplt.GetDim(0); i++)
        nldensplt(i)->write(log_label);
    }
    if(!pure_dmc){
      prop.endBlock();
    }
    else
      prop.endBlockSHDMC();
    
    if(max_fw_length){
      for(int s=0;s<fw_length.GetSize();s++){
        //prop_fw(s).endBlock();
        prop_fw(s).endBlock_per_step();
      }
    }

    totbranch=parallel_sum(totbranch);
    totkilled=parallel_sum(totkilled);
    totpoints=parallel_sum(totpoints);

    Properties_final_average finavg;
    prop.getFinal(finavg);
    eref=finavg.avg(Properties_types::total_energy,0);
    updateEtrial(feedback);
    
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
      if(!all_electron_moves)  dyngen->showStats(output);
      prop.printBlockSummary(output);

      output << "Branched "
	     << totbranch << " times.  So a branch every " 
	     << doublevar(totpoints)/doublevar(totbranch)
	     << " steps " << endl;
    }
    dyngen->resetStats();
  }
  
  if(output) {
    output << "\n ----------Finished DMC------------\n\n";
    prop.printSummary(output,average_var);  
    
    //MB: final printout for FW
    if(max_fw_length){
      for(int s=0;s<fw_length.GetSize();s++)
        prop_fw(s).printSummary(output,average_var);
    }

  }
  wfdata->clearObserver();
  deallocateIntermediateVariables();
}


//----------------------------------------------------------------------



void Dmc_method::move_electron_by_electron(Wavefunction_data * wfdata,
                                           Wavefunction * wf, 
                                           Sample_point * sample,
                                           Guiding_function * guidingwf,
                                           Dmc_point & pt,
                                           doublevar & acsum) { 
  Dynamics_info dinfo;
  for(int e=0; e< nelectrons; e++) {
    int acc;
    acc=dyngen->sample(e, sample, wf, wfdata, guidingwf,
        dinfo, timestep);

    if(dinfo.accepted) {               
      pt.age(e)=0;
    }
    else { 
      pt.age(e)++;
    }
    if(acc>0) acsum++;
  }

}


//----------------------------------------------------------------------

void Dmc_method::move_all_electron(Wavefunction_data * wfdata,
                         Wavefunction * wf, 
                         Sample_point * sample,
                         Guiding_function * guidingwf,
                         Dmc_point & pt,
                         doublevar & acsum) { 
  wf->updateLap(wfdata,sample);
  Wf_return lap(wf->nfunc(),5),vali(wf->nfunc(),1),valf(wf->nfunc(),1);
  doublevar signi,signf;
  int ndim=sample->ndim();
  Array1 < Array1 < doublevar> > drift1(nelectrons),drift2(nelectrons),
                              translate1(nelectrons),translate2(nelectrons),
                              gauss(nelectrons),pos0(nelectrons);
  for(int e=0; e < nelectrons; e++) {
    drift1(e).Resize(ndim);
    drift2(e).Resize(ndim);
    translate1(e).Resize(ndim);
    translate2(e).Resize(ndim);
    gauss(e).Resize(ndim);
    for(int d=0; d< ndim; d++) gauss(e)(d)=rng.gasdev();
    pt.age(e)=0;
    pos0(e).Resize(ndim);
    sample->getElectronPos(e,pos0(e));
  }

  signi=sample->overallSign();
  wf->getVal(wfdata,0,vali);
  
  for(int e=0; e < nelectrons; e++) { 
    wf->getLap(wfdata,e,lap);
    guidingwf->getLap(lap,drift1(e));
    limDrift(drift1(e),timestep,drift_cyrus);
    for(int d=0; d< ndim; d++) {
      translate1(e)(d)=gauss(e)(d)*sqrt(timestep)+drift1(e)(d);
    }
  }
  
  for(int e=0; e< nelectrons; e++) { 
    sample->translateElectron(e,translate1(e));
  }
  wf->notify(all_electrons_move,0);
  wf->updateLap(wfdata,sample);


  for(int e=0; e < nelectrons; e++) { 
    wf->getLap(wfdata,e,lap);
    guidingwf->getLap(lap,drift2(e));
    limDrift(drift2(e),timestep,drift_cyrus);
    for(int d=0; d< ndim; d++) {
      translate2(e)(d)=gauss(e)(d)*sqrt(timestep)+
           0.5*(drift1(e)(d)+drift2(e)(d))-translate1(e)(d);
    }
  }
  
  for(int e=0; e< nelectrons; e++) { 
    sample->translateElectron(e,translate2(e));
  }
  wf->notify(all_electrons_move,0);
  

  wf->updateLap(wfdata,sample);
  signf=sample->overallSign();
  wf->getVal(wfdata,0,valf);

  doublevar ratio=guidingwf->getTrialRatio(valf,vali)*signi*signf;
  if(ratio <0 and false) { 
    for(int e=0; e< nelectrons; e++){
      sample->setElectronPos(e,pos0(e));
    }
    wf->notify(all_electrons_move,0);
    
  }
  else { 
    acsum+=nelectrons;
  }
  
}

//----------------------------------------------------------------------



void Dmc_method::savecheckpoint(string & filename,                     
                                 Sample_point * config) {

  if(save_trace!="") { 
    if(mpi_info.node==0) { 
      cout << "entering trace write" << endl;
      long int time_ent=clock();
      FILE * f=fopen(save_trace.c_str(),"a");
      for(int i=0;i<nconfig; i++) {
        pts(i).config_pos.writeBinary(f,pts(i).weight);
      //  fwrite(&pts(i).weight, sizeof(doublevar),1, f);
      }

#ifdef USE_MPI
      Dmc_point tmppt;
      for(int p=1; p < mpi_info.nprocs; p++) { 
        int nconfigthis;
        MPI_Recv(nconfigthis,p);
        for(int i=0; i < nconfigthis; i++) { 
          tmppt.config_pos.mpiReceive(p);
          MPI_Recv(tmppt.weight,p);
          tmppt.config_pos.writeBinary(f,tmppt.weight);
        }
      }
#endif
      long int time_b=clock();
      single_write(cout,"writing to trace : ",double(time_b-time_ent)/CLOCKS_PER_SEC,"\n");
      fclose(f);
    }
#ifdef USE_MPI
    else { 
      MPI_Send(nconfig,0);
      for(int i=0; i< nconfig; i++) { 
        pts(i).config_pos.mpiSend(0);
        MPI_Send(pts(i).weight,0);
      }
    }
#endif
  }
  if(filename=="") return;
  write_configurations(filename, pts);


  return;
}


//----------------------------------------------------------------------



void Dmc_method::restorecheckpoint(string & filename, System * sys,
                                    Wavefunction_data * wfdata,
                                    Pseudopotential * pseudo) {

  ifstream is(filename.c_str());
  if(is) { 
    is.close();
    read_configurations(filename, pts);
  }
  else { 
    Array1 <Config_save_point>  configs;
    generate_sample(sample,wf,wfdata,guidingwf,nconfig,configs);
    pts.Resize(nconfig);
    for(int i=0; i< nconfig; i++) 
      pts(i).config_pos=configs(i);
  }
  int ncread=pts.GetDim(0);
  
  //cout << "ncread " << ncread << "  nwread " << nwread << endl;
  if(nconfig < ncread) { 
    Array1 <Dmc_point> tmp_pts(nconfig);
    for(int i=0; i< nconfig; i++) tmp_pts(i)=pts(i);
    pts=tmp_pts;
  }
  else if(nconfig > ncread) { 
    error("Not enough configurations in ", filename);
  }

  for(int walker=0; walker < nconfig; walker++) {
    pts(walker).config_pos.restorePos(sample);
    mygather.gatherData(pts(walker).prop, pseudo, sys,
                        wfdata, wf, sample,
                        guidingwf);
    pts(walker).age.Resize(sys->nelectrons(0)+sys->nelectrons(1));
    pts(walker).age=0;
  }
  find_cutoffs();

  updateEtrial(start_feedback);
    
}




//----------------------------------------------------------------------

void Dmc_method::updateEtrial(doublevar feed) {
  
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

doublevar Dmc_method::getWeight(Dmc_point & pt,
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
  return_weight=exp(teffac*(etr*2-effenergy-effoldenergy));

  return return_weight;
} 
                                
//----------------------------------------------------------------------
//using the pure DMC weight from last N steps in history
doublevar Dmc_method::getWeightPURE_DMC(Dmc_point & pt,
                                doublevar teff, doublevar etr) {
  doublevar teffac=teff;

  int history=pt.past_energies.size();
  assert(history>0);
  doublevar sum=0;
  for(int i=0;i<history;i++){
    sum+=pt.past_energies[i].main_en;
  }
  sum/=history;
  
  doublevar return_weight;
  return_weight=exp(teffac*history*(etr-sum));

  //cout <<"SHDDMC weight "<<return_weight<<" energy average "<<sum<<endl;
  //cout<<"current history size "<<history<<endl;
  return return_weight;
}

//----------------------------------------------------------------------

void Dmc_method::doTmove(Properties_point & pt,Pseudopotential * pseudo, System * sys,
                         Wavefunction_data * wfdata, Wavefunction * wf, Sample_point * sample,
                         Guiding_function * guideingwf) {
  vector <Tmove> tmov;
  pt.setSize(nwf);
  wf->getVal(wfdata,0,pt.wf_val);
  sys->calcKinetic(wfdata,sample,wf,pt.kinetic);
  pt.potential=sys->calcLoc(sample);
  pt.weight=1.0; //this gets set later anyway
  pt.count=1;
  pseudo->calcNonlocTmove(wfdata,sys,sample,wf,pt.nonlocal,tmov);
  doublevar sum=1;
  for(vector<Tmove>::iterator mov=tmov.begin(); mov!=tmov.end(); mov++) {
    assert(mov->vxx < 0);
    sum-=timestep*mov->vxx;
  }
  pt.nonlocal(0)-=(sum-1)/timestep;
  //subtract_out_enwt=-(sum-1)/timestep;
  assert(sum >= 0);
  if(tmoves) { ///Non-size consistent
    doublevar rand=rng.ulec()*sum;
    sum=1; //reset to choose the move
    if(rand > sum) {
      for(vector<Tmove>::iterator mov=tmov.begin(); mov!=tmov.end(); mov++) {
        sum-=timestep*mov->vxx;
        if(rand < sum) {
          sample->translateElectron(mov->e,mov->pos);
          break;
        }
      }
    }
  }
  else { // Size-consistent
    vector < vector<Tmove> > tmv_by_e(nelectrons);
    for(vector<Tmove>::iterator mov=tmov.begin(); mov!=tmov.end(); mov++) {
      tmv_by_e[mov->e].push_back(*mov);
    }
    for(int e=0; e< nelectrons; e++) {
      doublevar sum_e=1.0;
      for(vector<Tmove>::iterator mov=tmv_by_e[e].begin(); mov!=tmv_by_e[e].end(); mov++) {
        sum_e-=timestep*mov->vxx;
      }
      doublevar rand=rng.ulec()*sum_e;
      doublevar sel_sum=1;
      if(rand > sel_sum) {
        for(vector<Tmove>::iterator mov=tmv_by_e[e].begin(); mov!=tmv_by_e[e].end(); mov++) {
          sel_sum-=timestep*mov->vxx;
          if(rand < sel_sum) {
            sample->translateElectron(e,mov->pos);
            break;
          }
        }
      }
    }
  }

}
//----------------------------------------------------------------------
void Dmc_method::forwardWalking(int walker, int step, Array1<Properties_manager> & prop_fw) {
  
  //store the observables
  Dmc_history_avgrets new_avgrets;
  new_avgrets.avgrets=pts(walker).prop.avgrets;
  new_avgrets.weight=pts(walker).prop.weight(0);
  
  
  pts(walker).past_properties.push_front(new_avgrets);
  
  //tailor the array if larger than max_fw_length
  deque<Dmc_history_avgrets> & past_avgrets(pts(walker).past_properties);
  if(past_avgrets.size() > max_fw_length)
    past_avgrets.erase(past_avgrets.begin()+max_fw_length, past_avgrets.end());
  
  int size=pts(walker).past_properties.size();
  
  //find the element from the past
  doublevar oldweight;
  for(int s=0;s<fw_length.GetSize();s++){
    if(fw_length(s)>size){
      pts(walker).prop.avgrets=pts(walker).past_properties[size-1].avgrets;
      oldweight=pts(walker).past_properties[size-1].weight;
    }
    else{
      //call the prop.avgrets from fw_length(s) steps a go
      pts(walker).prop.avgrets=pts(walker).past_properties[fw_length(s)-1].avgrets;
      oldweight=pts(walker).past_properties[fw_length(s)-1].weight;
    }
    
    //insert it into observables
    prop_fw(s).insertPoint(step, walker, pts(walker).prop);
  }
}


//----------------------------------------------------------------------
//Auxilliary functions for branching
struct Queue_element {
  int from_node;
  int to_node;
  Queue_element() { }
  Queue_element(int from, int to) { from_node=from; to_node=to; }
};

struct weight_obj {
  double w;
  int i;
};

bool operator<(const weight_obj & a,const weight_obj & b) {
  return a.w < b.w;
}
void match_walkers(Array1<double> & weights, Array1 <int> & branch) {
  const double split_threshold=1.8;
  branch=-1;
  int totwalkers=weights.GetDim(0);
  Array1<weight_obj> walkers(totwalkers);
  for(int i=0; i< totwalkers; i++) {
    walkers(i).w=weights(i);
    walkers(i).i=i;
  }
  sort(walkers.v,walkers.v+totwalkers);
  int currsmallest=0;
  for(int i=totwalkers-1; i > 1; i--) { 
    if(walkers(i).w < split_threshold) break;
    int w=walkers(i).i;
    int smallest=walkers(currsmallest).i;
    double weight1=weights(w)/(weights(w)+weights(smallest));
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
    currsmallest++;
  }
  //for(int i=0; i< totwalkers; i++) {
  //  cout << walkers(i).w << endl;
  //}

}

//----------------------------------------------------------------------

int Dmc_method::calcBranch() { 
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

    long int time_a=clock();
    match_walkers(weights,branch);
    long int time_b=clock();
    single_write(cout,"matching: ",double(time_b-time_a)/CLOCKS_PER_SEC,"\n");
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
    time_a=clock();
#ifdef USE_MPI
    MPI_Bcast(nwalkers.v, mpi_info.nprocs, MPI_INT, mpi_info.node, MPI_Comm_grp);
    for(int i=1; i< mpi_info.nprocs; i++) {
      MPI_Send(branch.v+i*nconfig,nconfig,MPI_INT,i,0,MPI_Comm_grp);
      MPI_Send(weights.v+i*nconfig, nconfig, MPI_DOUBLE, i,0,MPI_Comm_grp);
    }
#endif
    time_b=clock();
    single_write(cout,"sending branch: ",double(time_b-time_a)/CLOCKS_PER_SEC,"\n");

               
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

  long int time_a=clock();
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
  long int time_b=clock();
  single_write(cout,"Finding out where to send: ",double(time_b-time_a)/CLOCKS_PER_SEC,"\n");
  
  time_a=clock(); 
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
  time_b=clock();
  single_write(cout,"sending walkers:",double(time_b-time_a)/CLOCKS_PER_SEC,"\n");
  
  return killsize;
  //exit(0);
}
//----------------------------------------------------------------------

//----------------------------------------------------------------------


void Dmc_point::mpiSend(int node) { 
#ifdef USE_MPI
  prop.mpiSend(node);
  config_pos.mpiSend(node);
  int n=past_energies.size();
  MPI_Send(&n,1, MPI_INT,node,0,MPI_Comm_grp);
  for(deque<Dmc_history>::iterator i=past_energies.begin();
      i!= past_energies.end(); i++) { 
    i->mpiSend(node);
  }
  
  //MB: sending the past_properties for the forward walking
  int m=past_properties.size();
  MPI_Send(&m,1, MPI_INT,node,0,MPI_Comm_grp);
  for(deque<Dmc_history_avgrets>::iterator i=past_properties.begin();
      i!= past_properties.end(); i++) { 
    i->mpiSend(node);
  }

  MPI_Send(&weight,1,MPI_DOUBLE,node,0,MPI_Comm_grp);
  MPI_Send(&ignore_walker,1, MPI_INT,node,0,MPI_Comm_grp);

  int nelectrons=age.GetDim(0);
  MPI_Send(&nelectrons,1, MPI_INT,node,0,MPI_Comm_grp);
  MPI_Send(age.v,nelectrons, MPI_DOUBLE, node,0,MPI_Comm_grp);
#endif
}

//----------------------------------------------------------------------

void Dmc_point::mpiReceive(int node) {
#ifdef USE_MPI
  prop.mpiReceive(node);
  config_pos.mpiReceive(node);
  int n;
  MPI_Status status;
  MPI_Recv(&n,1,MPI_INT,node,0,MPI_Comm_grp, &status);
  Dmc_history tmp_hist;
  past_energies.clear();
  for(int i=0; i< n; i++) {
    tmp_hist.mpiReceive(node);
    past_energies.push_back(tmp_hist);
  }

  //MB: receiving the past_properties for the forward walking
  int m;
  MPI_Recv(&m,1,MPI_INT,node,0,MPI_Comm_grp, &status);
  Dmc_history_avgrets tmp_prop;
  past_properties.clear();
  for(int i=0; i< m; i++) {
    tmp_prop.mpiReceive(node);
    past_properties.push_back(tmp_prop);
  }


  MPI_Recv(&weight,1,MPI_DOUBLE,node,0,MPI_Comm_grp, &status);
  MPI_Recv(&ignore_walker,1,MPI_INT,node,0,MPI_Comm_grp, &status);
      
  int nelectrons;
  MPI_Recv(&nelectrons,1,MPI_INT,node,0,MPI_Comm_grp,&status);
  age.Resize(nelectrons);
  MPI_Recv(age.v,nelectrons,MPI_DOUBLE, node,0,MPI_Comm_grp,&status);
#endif
}
//----------------------------------------------------------------------

void Dmc_point::write(ostream & os) { 
  string indent="";
  config_pos.write(os);
  os << "DMC { \n";
  //prop.write(indent,os);
  os << "weight " << weight<< endl;
  os << "sign " << sign << endl;
  /*
  for(deque<Dmc_history>::iterator i=past_energies.begin(); 
      i!=past_energies.end(); i++) { 
    os << "past_energies { ";
    i->write(os);
    os << "}\n";    
  }
  for(deque<Dmc_history_avgrets>::iterator i=past_properties.begin(); 
      i!=past_properties.end(); i++) {
    os << "past_properties { ";
    i->write(os);
    os << "}\n";
  }
   */
  os << "}\n";
}

void Dmc_point::read(istream & is) { 
  config_pos.read(is);
  int filepos=is.tellg();
  string dum;
  is >> dum;
  if(!caseless_eq(dum, "DMC")) {
    is.seekg(filepos);
    return;
  }
  
  is >> dum; //the {
  //prop.read(is);
  is >> dum >> weight;
  is >> dum >> sign;
  //ignoring the past stuff for the moment..
}


//----------------------------------------------------------------------
void Dmc_history_avgrets::write(ostream & os) { 
  os << "weight " << weight << endl;
  for(int i=0; i< avgrets.GetDim(0); i++) { 
    os << "avgret { ";
  }
    
}

void Dmc_history_avgrets::read(istream & is) { 
}

//----------------------------------------------------------------------

void Dmc_history::mpiSend(int node) { 
#ifdef USE_MPI
  MPI_Send(&main_en,1,MPI_DOUBLE,node,0,MPI_Comm_grp);
#endif
}


void Dmc_history::mpiReceive(int node) { 
#ifdef USE_MPI
  MPI_Status status;
  int n1,n2;
  MPI_Recv(&main_en,1,MPI_DOUBLE,node,0,MPI_Comm_grp,&status);
#endif
}


void Dmc_history_avgrets::mpiSend(int node) { 
#ifdef USE_MPI
  MPI_Send(&weight,1,MPI_DOUBLE,node,0,MPI_Comm_grp);
  int n1=avgrets.GetDim(1);
  MPI_Send(&n1,1,MPI_INT,node,0,MPI_Comm_grp);
  for(int j=0;j<n1;j++){
    int n2=avgrets(0,j).vals.GetDim(0);
    MPI_Send(&n2,1,MPI_INT,node,0,MPI_Comm_grp);
    MPI_Send(avgrets(0,j).vals.v,n2,MPI_DOUBLE,node,0,MPI_Comm_grp);
  }
#endif
}

void Dmc_history_avgrets::mpiReceive(int node) { 
#ifdef USE_MPI
  MPI_Status status;
  int n1,n2;
  MPI_Recv(&weight,1,MPI_DOUBLE,node,0,MPI_Comm_grp,&status);
  MPI_Recv(&n1,1,MPI_INT,node,0,MPI_Comm_grp,&status);
  avgrets.Resize(1,n1);
  for(int j=0;j<n1;j++){
    MPI_Recv(&n2,1,MPI_INT,node,0,MPI_Comm_grp,&status);
    avgrets(0,j).vals.Resize(n2);
    MPI_Recv(avgrets(0,j).vals.v,n2,MPI_DOUBLE,node,0,MPI_Comm_grp,&status);  
  }
#endif
}
