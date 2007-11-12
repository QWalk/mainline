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



#include "Dmc_method.h"
#include "qmc_io.h"
#include "ulec.h"
#include "Program_options.h"
#include "average.h"


void Dmc_method::read(vector <string> words,
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
  
  
  vector <string> dynamics_words;
  if(!readsection(words, pos=0, dynamics_words, "DYNAMICS") ) 
    dynamics_words.push_back("SPLIT");

  low_io=0;
  if(haskeyword(words, pos=0,"LOW_IO")) low_io=1;
  
  allocate(dynamics_words, dyngen);
  dyngen->enforceNodes(1);

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
  os << "Diffusion Monte Carlo: $Date: 2006/12/05 00:17:04 $\n";
  os << "Number of processors " <<           mpi_info.nprocs << endl;
  os << "Blocks: " <<                        nblock    << endl;
  os << "Steps per block: " <<               nstep     << endl;
  os << "Timestep: " <<                      timestep  << endl;
  if(tmoves) 
    os << "T-moves turned on" << endl;
  string indent="  ";

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

  aux_converge=1;
  
  prop.setSize(wf->nfunc(), nblock, nstep, nconfig, sys, 
	       wfdata, mygather.nAux(), aux_converge);

  restorecheckpoint(readconfig, sys, wfdata, pseudo);


  int naux=mygather.nAux();

  nhist=1;
  //setting the projection time for auxillary walkers to 1 a.u.
  if(naux >0 && nhist <0) nhist=int(1.0/timestep)+1;
  if(naux > 0 && tmoves) error("Can't do t-moves and auxillary walks yet");
  //if(tmoves) pseudo->setDeterministic(1); //this may not be necessary..not sure yet.
  
  doublevar teff=timestep;
  Array1 <doublevar> aux_timestep(naux,teff);

  for(int block=0; block < nblock; block++) {

    int totkilled=0;  
    int totbranch=0;
    int totpoints=0;
    for(int step=0; step < nstep; ) {
      int npsteps=min(feedback_interval, nstep-step);

      Dynamics_info dinfo;
      doublevar acsum=0;
      doublevar deltar2=0;
      Array1 <doublevar> aux_rf_diffusion(naux,0.0);
      Array1 <doublevar> epos(3);
      
      doublevar avg_acceptance=0;
      
      doublevar rf_diffusion=0; //diffusion rate without rejection

      
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
            
            Array1 <Dynamics_info> aux_dinfo(naux);
            mygather.getGreensFunctions(dyngen, dinfo, e,sample, guidingwf,
                                        timestep, aux_dinfo, 0);
            
            if(dinfo.accepted) { 
              rf_diffusion+=dinfo.diffusion_rate/(nconfig*nelectrons*npsteps);
              for(int i=0; i< naux; i++) 
                aux_rf_diffusion(i)+=
                  aux_dinfo(i).diffusion_rate/(nconfig*nelectrons*npsteps);
              
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
            pt.setSize(nwf,0,0);
            wf->getVal(wfdata,0,pt.wf_val);
            sys->calcKinetic(wfdata,sample,wf,pt.kinetic);
            pt.potential=sys->calcLoc(sample);
            pt.weight=1.0; //this gets set later anyway
            pt.count=1;
            pseudo->calcNonlocTmove(wfdata,sample,wf,pt.nonlocal,tmov);
            getZpol(sys, sample, pt.z_pol,1); //always do the many-body zpol, because it's correct
            //cout << "choosing among " <<  tmov.size() << " tmoves " << endl;
            //Now we do the t-move
            doublevar sum=0; 
            for(vector<Tmove>::iterator mov=tmov.begin(); mov!=tmov.end(); mov++) { 
              assert(mov->vxx < 0);
              sum-=timestep*mov->vxx;  
            }
            pt.nonlocal(0)-=sum/timestep;
            subtract_out_enwt=-sum/timestep;
            //cout << "sum " << sum <<  " nonlocal " << pt.nonlocal(0) << " ratio " << sum/pt.nonlocal(0) << endl;
            assert(sum >= 0);
            doublevar rand=rng.ulec()*sum;
            sum=0; //reset to choose the move
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
            //wf->updateLap(wfdata, sample);
          } ///---------------------------------done with the T-moves
          else { 
            mygather.gatherData(pt, pseudo, sys, wfdata, wf, 
                                sample, guidingwf, aux_converge,0);
          }
          
          //the history is only used for setting the weights,
          //so we just subtract out anything that was used for t-moves..
          Dmc_history new_hist;
          new_hist.main_en=pts(walker).prop.energy(0);
          new_hist.aux_en=pts(walker).prop.aux_energy;
          //cout << " main en " << pts(walker).prop.energy(0) << endl;
          //cout << " aux en " << pts(walker).prop.aux_energy(0,0) << endl;
          pts(walker).past_energies.push_front(new_hist);
          deque<Dmc_history> & past(pts(walker).past_energies);
          if(past.size() > nhist) 
            past.erase(past.begin()+nhist, past.end());
          
          pts(walker).prop=pt;
          pts(walker).weight*=getWeight(pts(walker),teff,etrial);
          if(pts(walker).ignore_walker) {
            pts(walker).ignore_walker=0;
            pts(walker).weight=1;
            pts(walker).prop.count=0;
          }
          pts(walker).prop.weight=pts(walker).weight;
          getAuxWeight(pts(walker), teff, aux_timestep,pts(walker).prop.aux_weight);
          for(int a=0; a< naux; a++) pts(walker).prop.aux_weight(a,0)*=pts(walker).weight;
          //This is somewhat inaccurate..will need to change it later
          //For the moment, the autocorrelation will be slightly
          //underestimated
          pts(walker).prop.parent=walker;
          pts(walker).prop.nchildren=1;
          pts(walker).prop.children(0)=walker;
          
          prop.insertPoint(step+p, walker, pts(walker).prop);
          for(int i=0; i< densplt.GetDim(0); i++)
            densplt(i)->accumulate(sample,pts(walker).prop.weight(0));
        }
        
        pts(walker).config_pos.savePos(sample);
      }
      //---Finished moving all walkers

      //-----------Find the effective timesteps(of main and aux walks)
      doublevar accept_ratio=acsum/(nconfig*nelectrons*npsteps);
      teff=timestep*accept_ratio; //deltar2/rf_diffusion; 
      for(int i=0; i< naux; i++)  {
        aux_timestep(i)=teff*(aux_rf_diffusion(i)/rf_diffusion);
      }

      updateEtrial(feedback);

      step+=npsteps;

      int nkilled=calcBranch();
      totkilled+=nkilled;
      totbranch+=nkilled;
    }

    ///----Finished block
    
    if(!low_io || block==nblock-1) {
      savecheckpoint(storeconfig,sample);
      for(int i=0; i< densplt.GetDim(0); i++)
        densplt(i)->write();      
    }
    prop.endBlock();


    //loadbalance(); //not needed with global branching


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
      dyngen->showStats(output);

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
    prop.printSummary(output);  
  }
  wfdata->clearObserver();
  deallocateIntermediateVariables();
}


//----------------------------------------------------------------------


void Dmc_method::savecheckpoint(string & filename,                     
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



void Dmc_method::restorecheckpoint(string & filename, System * sys,
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

void Dmc_method::cdmcReWeight(Array2 <doublevar> & energy_temp, 
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

void Dmc_method::updateEtrial(doublevar feed) {
  
  doublevar totweight=0;
  for(int walker=0; walker < nconfig; walker++)
    totweight+=pts(walker).weight;

  etrial=eref-feed*log(totweight/double(nconfig));

#ifdef USE_MPI
  doublevar mpitot=0;
  MPI_Allreduce(&totweight,&mpitot, 1,
                MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
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
  //if(tmoves) return_weight=exp(teffac*2.0*(etr-effenergy));
  //else 
  return_weight=exp(teffac*(etr*2-effenergy-effoldenergy));

  return return_weight;
} 
                                
//----------------------------------------------------------------------


void Dmc_method::getAuxWeight(Dmc_point & pt, 
                              doublevar teff, 
                              Array1 <doublevar>&  aux_teff, 
			      Array2 <doublevar> & aux_weight) {


  int naux=pt.prop.aux_energy.GetDim(0);
  doublevar etr=eref;
  aux_weight.Resize(naux,1);
  for(int a=0; a< naux; a++) { 
    doublevar w=1;
    w*=pt.prop.aux_jacobian(a);
    doublevar ratio=guidingwf->getTrialRatio(pt.prop.aux_wf_val(a), 
					     pt.prop.wf_val);
    w*= ratio*ratio;

    doublevar fbet=max(etr-pt.prop.aux_energy(a,0), etr-pt.prop.energy(0));
    if(fbet > branchcut_stop) w=0;
    w*=exp(-.5*(aux_teff(a)*pt.prop.aux_energy(a,0)-teff*pt.prop.energy(0)));
    for(deque<Dmc_history>::iterator i=pt.past_energies.begin();
        i!=pt.past_energies.end()-1; i++) {
      doublevar fbet=max(etr-i->main_en, etr-i->aux_en(a,0));
      if(fbet < branchcut_stop) 
        w*=exp(-(aux_teff(a)*i->aux_en(a,0)-teff*i->main_en));
    }
    Dmc_history & last(*(pt.past_energies.end()-1));
     fbet=max(etr-last.main_en, etr-last.aux_en(a,0));
    if(fbet< branchcut_stop) 
      w*=exp(-(aux_teff(a)*last.aux_en(a,0)-teff*last.main_en));
    aux_weight(a,0)=w; 
  }
  
}

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



int Dmc_method::calcBranch() { 
  int totwalkers=mpi_info.nprocs*nconfig;
  Array1 <doublevar> weights(totwalkers);
  Array1 <doublevar> my_weights(nconfig);
  
  for(int walker=0; walker < nconfig; walker++)
    my_weights(walker)=pts(walker).weight;
#ifdef USE_MPI
  MPI_Allgather(my_weights.v,nconfig, MPI_DOUBLE, weights.v,nconfig,MPI_DOUBLE, MPI_COMM_WORLD);
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
    MPI_Bcast(nwalkers.v, mpi_info.nprocs, MPI_INT, mpi_info.node, MPI_COMM_WORLD);
    for(int i=1; i< mpi_info.nprocs; i++) {
      MPI_Send(branch.v+i*nconfig,nconfig,MPI_INT,i,0,MPI_COMM_WORLD);
      MPI_Send(weights.v+i*nconfig, nconfig, MPI_DOUBLE, i,0,MPI_COMM_WORLD);
    }
#endif
               
  }
  else { 
#ifdef USE_MPI
    MPI_Bcast(nwalkers.v, mpi_info.nprocs, MPI_INT, root, MPI_COMM_WORLD);
    MPI_Status status;
    MPI_Recv(my_branch.v,nconfig, MPI_INT,root,0,MPI_COMM_WORLD, &status);
    MPI_Recv(my_weights.v,nconfig, MPI_DOUBLE,root,0,MPI_COMM_WORLD, &status);
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
/*
int Dmc_method::calcBranch() {

  Array1 <int> branch(nconfig);
  branch=0;

  //Array1 <int> killwalkers(nconfig);
  //killwalkers=0;
  int killsize=0;
  //branchout << "********branching step " << bstep++ << endl;

  branch=-1;
  for(int walker=0; walker < nconfig; walker++) {
    const doublevar split_threshold=1.8;
    if(pts(walker).weight > split_threshold 
       && branch(walker)==-1) {
      //branchout << "should branch " << walker << " weight " << dmcweight(walker) << endl;
      doublevar smallestgr=100;
      int smallest=-1;
      for(int walker2=0; walker2 < nconfig; walker2++) {
        if(branch(walker2)==-1 && walker2!=walker) {
          if(pts(walker2).weight < smallestgr) {
            smallestgr=pts(walker2).weight;
            smallest=walker2;
          }
        }
      }
      if(smallest!=-1) {
        //cout  << "combining " << walker << " and " 
        //          << smallest << " (weight= " << dmcweight(smallest) << endl;
        doublevar weight1=pts(walker).weight
	  /(pts(walker).weight+pts(smallest).weight);
        
        if(weight1+rng.ulec() >=1.0) {
          branch(walker)=2;
          branch(smallest)=0;
          pts(walker).weight+=pts(smallest).weight;
          pts(walker).weight/=2;
        }
        else {
          branch(walker)=0;
          branch(smallest)=2;
          pts(smallest).weight+=pts(walker).weight;
          pts(smallest).weight/=2;
        }
        killsize++;
      }
    }
  }

  for(int walker=0; walker < nconfig; walker++) {
    if(branch(walker)==-1)
      branch(walker)=1;
  }

  int newnconfig=0;
  for(int walker=0; walker < nconfig; walker++) {
    newnconfig+=branch(walker);
  }

  if(newnconfig != nconfig) 
    error("Error in branching..newnconfig != nconfig");


  Array1 <Dmc_point> newpts;
  newpts=pts;

  int counter=0;
  for(int walker=0; walker < nconfig; walker++)  {
    for(int i=0; i< branch(walker); i++)  {
      assert(counter < newnconfig);
      pts(counter)=newpts(walker);
      //pts(walker).prop.children(pts[walker].prop.nchildren++)=counter;
      counter++;
    }
  }

  assert(counter == newnconfig);
  //Finally, reset nconfig
  nconfig=newnconfig;  
  return killsize;
}
*/
//----------------------------------------------------------------------


#include <algorithm>
struct procwt {
  doublevar totwt;
  int procnum;
};
bool operator<(const procwt & p1,const procwt & p2) {
  return p1.totwt < p2.totwt;
}

void Dmc_method::loadbalance() {
#ifdef USE_MPI
  //string ldout="load";
  //canonical_filename(ldout);
  //ofstream loadout(ldout.c_str(),ios::app);
  //ostream & loadout(cout);

  doublevar totwt=0;
  for(int i=0; i< nconfig; i++) 
    totwt+=pts(i).weight;

  double * weights=new double[mpi_info.nprocs];
  MPI_Allgather(&totwt, 1, MPI_DOUBLE, weights, 1, MPI_DOUBLE,
                MPI_COMM_WORLD);

  vector <procwt> procs(mpi_info.nprocs);
  int nprocs=mpi_info.nprocs;
  for(int i=0; i< nprocs; i++) {
    procs[i].totwt=weights[i];
    procs[i].procnum=i;
  }
  sort(procs.begin(), procs.end());
  
  int myplace=-1;
  for(int i=0; i< nprocs; i++) 
    if(procs[i].procnum==mpi_info.node) {
      myplace=i;
      break;
    }
  assert(myplace != -1);

  //Find the biggest and smallest weights
  //on this processor
  int big=-1, small=-1;
  doublevar bigw=-1, smallw=1000;
  for(int i=0; i< nconfig; i++) {
    if(pts(i).weight > bigw) {
      bigw=pts(i).weight;
      big=i;
    }
    if(pts(i).weight < smallw) {
      smallw=pts(i).weight;
      small=i;
    }
  }
      
  //loadout << mpi_info.node << " here " << myplace << endl;
//MPI_Barrier(MPI_COMM_WORLD);

  int mate=0;
  if(myplace > nprocs/2) {
    mate=procs[nprocs-myplace-1].procnum;
    pts[big].mpiSend(mate);
    //loadout << "sent " << endl;
    pts(big).mpiReceive(mate);
  }
  else if(nprocs-myplace-1 > nprocs/2) {
    mate=procs[nprocs-myplace-1].procnum;
    Dmc_point tmp_pt;
    tmp_pt.mpiReceive(mate);
    pts(small).mpiSend(mate);
    pts(small)=tmp_pt;

  }
  

  delete [] weights;


#endif
}

//----------------------------------------------------------------------


void Dmc_point::mpiSend(int node) { 
#ifdef USE_MPI
  prop.mpiSend(node);
  config_pos.mpiSend(node);
  int n=past_energies.size();
  MPI_Send(&n,1, MPI_INT,node,0,MPI_COMM_WORLD);
  for(deque<Dmc_history>::iterator i=past_energies.begin();
      i!= past_energies.end(); i++) { 
    i->mpiSend(node);
  }
  
  MPI_Send(&weight,1,MPI_DOUBLE,node,0,MPI_COMM_WORLD);
  MPI_Send(&ignore_walker,1, MPI_INT,node,0,MPI_COMM_WORLD);

  int nelectrons=age.GetDim(0);
  MPI_Send(&nelectrons,1, MPI_INT,node,0,MPI_COMM_WORLD);
  MPI_Send(age.v,nelectrons, MPI_DOUBLE, node,0,MPI_COMM_WORLD);
#endif
}

void Dmc_point::mpiReceive(int node) {
#ifdef USE_MPI
  prop.mpiReceive(node);
  config_pos.mpiReceive(node);
  int n;
  MPI_Status status;
  MPI_Recv(&n,1,MPI_INT,node,0,MPI_COMM_WORLD, &status);
  Dmc_history tmp_hist;
  past_energies.clear();
  for(int i=0; i< n; i++) {
    tmp_hist.mpiReceive(node);
    past_energies.push_back(tmp_hist);
  }
  MPI_Recv(&weight,1,MPI_DOUBLE,node,0,MPI_COMM_WORLD, &status);
  MPI_Recv(&ignore_walker,1,MPI_INT,node,0,MPI_COMM_WORLD, &status);
      
  int nelectrons;
  MPI_Recv(&nelectrons,1,MPI_INT,node,0,MPI_COMM_WORLD,&status);
  age.Resize(nelectrons);
  MPI_Recv(age.v,nelectrons,MPI_DOUBLE, node,0,MPI_COMM_WORLD,&status);
#endif
}

void Dmc_history::mpiSend(int node) { 
#ifdef USE_MPI
  MPI_Send(&main_en,1,MPI_DOUBLE,node,0,MPI_COMM_WORLD);
  int n1=aux_en.GetDim(0);
  int n2=aux_en.GetDim(1);
  MPI_Send(&n1,1,MPI_INT,node,0,MPI_COMM_WORLD);
  MPI_Send(&n2,1,MPI_INT,node,0,MPI_COMM_WORLD);
  MPI_Send(aux_en.v,n1*n2,MPI_DOUBLE,node,0,MPI_COMM_WORLD);
#endif
}

void Dmc_history::mpiReceive(int node) { 
#ifdef USE_MPI
  MPI_Status status;
  int n1,n2;
  MPI_Recv(&main_en,1,MPI_DOUBLE,node,0,MPI_COMM_WORLD,&status);
  MPI_Recv(&n1,1,MPI_INT,node,0,MPI_COMM_WORLD,&status);
  MPI_Recv(&n2,1,MPI_INT,node,0,MPI_COMM_WORLD,&status);
  aux_en.Resize(n1,n2);
  MPI_Recv(aux_en.v,n1*n2,MPI_DOUBLE,node,0,MPI_COMM_WORLD,&status);  
#endif
}
