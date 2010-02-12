/*
 
Copyright (C) 2010 Lucas K. Wagner

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



#include "Dmc_corr_method.h"
#include "qmc_io.h"
#include "ulec.h"
#include "Program_options.h"
#include "average.h"


void Dmc_corr_method::read(vector <string> words,
                      unsigned int & pos,
                      Program_options & options)
{

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

  //if(readvalue(words, pos=0, readconfig, "READCONFIG"))
  //  canonical_filename(readconfig);
  //else
  //  error("Must give READCONFIG for DMC");

  //optional options

  if(!readvalue(words, pos=0, nhist, "CORR_HIST")) 
    nhist=50;

  if(!readvalue(words, pos=0,npsteps, "PARTIAL_STEPS")) 
    npsteps=5;

  if(readvalue(words, pos=0, storeconfig, "STORECONFIG"))
    canonical_filename(storeconfig);

  if(!readvalue(words, pos=0, log_label, "LABEL"))
    log_label="dmc";
  
  vector <string> warpsec;
  if(readsection(words, pos=0, warpsec,"WARPER")) {
    warper.read(warpsec);
  }
  
  if(haskeyword(words, pos=0, "NOWARP")) 
    warper.set_warp(0);
/*
  if(!readvalue(words, pos=0, start_feedback, "START_FEEDBACK"))
    start_feedback=1;

  if(readvalue(words, pos=0, feedback_interval, "FEEDBACK_INTERVAL")) {
    if(feedback_interval < 1) 
      error("FEEDBACK_INTERVAL must be greater than or equal to 1");
  }
  else feedback_interval=5;

  if(!readvalue(words, pos=0, feedback, "FEEDBACK"))
    feedback=1.0;
*/
  //if(!readvalue(words, pos=0, branch_start_cutoff, "BRANCH_START_CUTOFF")) 
  //  branch_start_cutoff=10;
  

  //branch_stop_cutoff=branch_start_cutoff*1.5;
  
  
  vector <string> dynamics_words;
  if(!readsection(words, pos=0, dynamics_words, "DYNAMICS") ) 
    dynamics_words.push_back("SPLIT");

  allocate(dynamics_words, dyngen);
  dyngen->enforceNodes(1);


}

//----------------------------------------------------------------------

int Dmc_corr_method::generateVariables(Program_options & options) {

  if(!have_read_options) 
    error("need to call Dmc_corr_method::read before generateVariables");
  if(have_allocated_variables) 
    error("already allocated variables in Dmc_corr_method");

  have_allocated_variables=1;
  int nsys=options.systemtext.size();
  mysys.Resize(nsys);
  mywfdata.Resize(nsys);

  for(int i=0; i< nsys; i++) { 
    mysys(i)=NULL; mywfdata(i)=NULL;
    allocate(options.systemtext[i], mysys(i));
    allocate(options.twftext[i], mysys(i), mywfdata(i));
  }
  mysys(0)->generatePseudo(options.pseudotext, mypseudo);
  
  
  return 1;
}

//----------------------------------------------------------------------




int Dmc_corr_method::showinfo(ostream & os) {
  int nsys=mysys.GetDim(0);
  for(int i=0; i< nsys; i++) { 
    mysys(i)->showinfo(os);
    mywfdata(i)->showinfo(os);
  }
    
  mypseudo->showinfo(os);
  os << "###########################################################\n";
  os << "Correlated Diffusion Monte Carlo:\n";
  os << "Number of processors " <<           mpi_info.nprocs << endl;
  os << "Blocks: " <<                        nblock    << endl;
  os << "Steps per block: " <<               nstep     << endl;
  os << "Timestep: " <<                      timestep  << endl;
  os << "Number of systems: " <<             mywfdata.GetDim(0) << endl;
  string indent="  ";

  dyngen->showinfo(indent, os);

  os << "###########################################################" << endl;
  return 1;
}

//----------------------------------------------------------------------


//----------------------------------------------------------------------

void Dmc_corr_method::run(Program_options & options, ostream & output) {
  if(!have_allocated_variables) 
    error("Must generate variables to use Dmc_corr_method::run");
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
  Array1 <Average_generator*> avg_gen;
  myprop.initializeLog(avg_gen);
  int nsys=mywfdata.GetDim(0);
  assert(mysys.GetDim(0)==nsys);
  wf.Resize(nsys);
  sample.Resize(nsys);
  for(int i=0; i< nsys; i++) { 
    wf(i)=NULL;
    sample(i)=NULL;
    mywfdata(i)->generateWavefunction(wf(i));
    mysys(i)->generateSample(sample(i));
    sample(i)->attachObserver(wf(i));
  }
  guidingwf=new Primary;
  pts.Resize(nconfig);
  int nelectrons=mysys(0)->nelectrons(0)+mysys(0)->nelectrons(1);
  for(int i=0; i< nconfig; i++) { 
    pts(i).age.Resize(nelectrons);
    pts(i).age=0;
    pts(i).system= i%nsys;
    pts(i).jacobian.Resize(nsys);
    pts(i).jacobian=1;
    //cout << "walker " << i << " system " << pts(i).system << endl;
  }
  
  ///-----------------Generate VMC configs
  for(int i=0; i< nconfig; i++) {
    int currsys=pts(i).system;
    sample(currsys)->randomGuess();
    pts(i).config_pos.savePos(sample(currsys));
  }
  
  int nblock_vmc=5;
  int nstep_vmc=10;
  doublevar timestep_vmc=1.0;
  for(int block=0; block < nblock_vmc; block++) {
    for(int walker=0; walker < nconfig; walker++) { 
      int currsys=pts(walker).system;
      pts(walker).config_pos.restorePos(sample(currsys));
      wf(currsys)->updateLap(mywfdata(currsys), sample(currsys));
      for(int step=0; step < nstep_vmc; step++) { 
        Dynamics_info dinfo;
        for(int e=0; e< nelectrons; e++) {
          dyngen->sample(e,sample(currsys), wf(currsys),
                         mywfdata(currsys), guidingwf, dinfo, timestep_vmc);
        }
      }
      //cout << "finished walker " << walker << " sysy " << currsys << endl;
      pts(walker).config_pos.savePos(sample(currsys));
    }
  }
  myprop.setSize(nsys, nblock, nstep+1, nconfig, mysys(0), 
                 mywfdata(0), 0, 0);

  ///cout << "here " << endl;
  eref.Resize(nsys);
  eref=0;
  Array1 <int> nwalkers_sys(nsys);
  nwalkers_sys=0;
  
  for(int walker=0; walker < nconfig; walker++) { 
    //cout << "add_point " << endl;
    int currsys=pts(walker).system;
    pts(walker).config_pos.restorePos(sample(currsys));
    add_point(walker);
    pts(walker).initial_sign.Resize(nsys);
    for(int i=0; i< nsys; i++) { 
      Wf_return wfval(wf(i)->nfunc(), 2);
      wf(i)->getVal(mywfdata(i), 0, wfval);
      pts(walker).initial_sign(i)=wfval.sign(0);
    }
    //cout << "eref " << endl;
    eref(pts(walker).system)+=pts(walker).prop.energy(currsys);
    //cout << "energy " << pts(walker).prop.energy(pts(walker).system) << endl;
    nwalkers_sys(currsys)++;
  }

  
  for(int sys=0; sys< nsys; sys++) { 
    eref(sys)=parallel_sum(eref(sys));
    nwalkers_sys(sys)=parallel_sum(nwalkers_sys(sys));
  }
  
  for(int sys=0; sys < nsys; sys++) { 
    eref(sys)/=nwalkers_sys(sys);
    cout << mpi_info.node <<":sys " << sys << " eref " << eref(sys) 
         << " nwalkers " << nwalkers_sys(sys) << endl;
  }
  etrial=eref;
  ///cout << "done " << endl;
  //exit(0);
  //------------------------------------------
  
  for(int block=0; block < nblock; block++) { 
    for(int step=0; step < nstep; ) { 
      //cout << "block " << block << " step " << step << endl;
      Dynamics_info dinfo;
      
      doublevar avg_acceptance=0;  
      int n_partial=min(npsteps,nstep-step);
      //cout << mpi_info.node << ":move" << endl;
      int ncross=0;
      for(int walker=0; walker < nconfig; walker++) {
        int currsys=pts(walker).system;
        pts(walker).config_pos.restorePos(sample(currsys));
        //wf(currsys)->updateLap(mywfdata(currsys), sample(currsys));
        for(int i=0; i < nsys; i++) { 
          wf(i)->updateLap(mywfdata(i), sample(i));
        }
        Array2 <Config_save_point> configs(n_partial,nsys );
        Array2 <double> values(n_partial, nsys);
        //cout << "start partial " << endl;
        for(int p=0; p < n_partial; p++) { 
          mypseudo->randomize();            
          
          
          propagate_walker(walker);
          myprop.insertPoint(step+p, walker, pts(walker).prop);
          //------checking for node crossing
          
          for(int i=0; i< nsys; i++) { 
            configs(p,i).savePos(sample(i));
            cout.precision(15);
            Wf_return wfval(wf(i)->nfunc(), 2);
            wf(i)->getVal(mywfdata(i),0,wfval);
            values(p,i)=wfval.amp(0,0);
            if(wfval.sign(0) != pts(walker).initial_sign(i)) { 
              //cout << "Node crossed! system "<< i << " " << endl;
              pts(walker).initial_sign(i)*=-1;
              ncross++;
              if(p == n_partial-1) { 
                cout << "crossed node in system " << i <<  " currsys " << currsys << endl;
                for(int sy=0; sy< nsys; sy++) { 
                  cout << "----system " << sy << " amplitude " << endl;
                  for(int j=0; j<=p; j++) { 
                    cout << "amp " << values(j,i) << endl;
                  }
                }
              }
            }
          }//---finished node crossing checking 
        }
        pts(walker).config_pos.savePos(sample(currsys));

      }
      if(ncross > 0) { 
        cout << "proportion crossed " << double(ncross)/(nsys*nconfig*n_partial) 
        << endl;
      }
      step+=n_partial;
      //cout << mpi_info.node << ":branch" << endl;
      //calcBranch();
      //cout << mpi_info.node << ":etrial" << endl;

      updateEtrial();
      if(output) { 
        for(int i=0; i< nsys; i++) { 
          output << "etrial " << i << "  " << etrial(i) << endl;
        }
      }
    }
    myprop.endBlock();

    Properties_final_average finavg;
    myprop.getFinal(finavg);
    for(int i=0; i< nsys; i++) { 
      eref(i)=finavg.avg(Properties_types::total_energy,i);    
    }
    
    if(output) { 
      myprop.printBlockSummary(output);
    }
  }
  
  //loop over blocks 
  //loop over systems1
  
  deallocateIntermediateVariables();
}


//----------------------------------------------------------------------

//Move an walker by one DMC step, updating the pts(walker) variable and leaving
//the samples and wave functions at the new position.
//it's assumed that at least sample(pts(walker).system) is updated for that walker.

doublevar Dmc_corr_method::propagate_walker(int walker) { 
  int currsys=pts(walker).system;
  int nsys=mywfdata.GetDim(0);
  Dynamics_info dinfo;
  doublevar avg_acceptance=0;
  int nelectrons=mysys(currsys)->nelectrons(0)+mysys(currsys)->nelectrons(1);
  //------Do several steps without branching
  
  for(int e=0; e< nelectrons; e++) {
    int acc;
    acc=dyngen->sample(e, sample(currsys), wf(currsys),
                       mywfdata(currsys), guidingwf,dinfo, timestep);
    
    if(dinfo.accepted)  {
      pts(walker).age(e)=0;
    }
    else { 
      pts(walker).age(e)++;
    }
    avg_acceptance+=dinfo.acceptance/(nelectrons*npsteps);
    
    //if(acc>0) acsum++;
  }
  add_point(walker);
  
  
  //----------------------------------------Finally, update the weights
  
  //The weight for the primary system
  deque<Dmc_corr_history>::iterator h=pts(walker).past_energies.begin();
  //pts(walker).weight*=exp(-timestep*0.5*(h->main_en(currsys)
  //                                       +(h+1)->main_en(currsys)-2*etrial(currsys)));
  //cout << "weight " << walker << "  " << pts(walker).weight << endl;
  Array1 <doublevar> logpi(nsys); //the logarithm of the pdf for each system
  logpi=0;
  int pc_gf=1;
  for(int i=0; i< nsys; i++) { 
    for(deque<Dmc_corr_history>::iterator h=pts(walker).past_energies.begin(); 
        h!= pts(walker).past_energies.end()-1; h++) { 
      logpi(i)-= 0.5*timestep*(h->main_en(i)+(h+1)->main_en(i) - 2*eref(i));
    }
    if(i==currsys) { 
       pts(walker).weight=exp(logpi(i));
       //cout << "walker " << walker << " weight " << pts(walker).weight << " projection " 
       //     << pts(walker).past_energies.size() << endl;
    }
    
    
    if(pc_gf) { //---------
      for(deque<Dmc_corr_history>::iterator h=pts(walker).past_energies.begin(); 
          h!= pts(walker).past_energies.end()-1; h++) { 
        Array1 <doublevar> pos1(3), pos2(3), drift1(3), drift2(3);
        for(int e=0; e< nelectrons; e++) { 
          h->configs(i).getPos(e,pos1);
          (h+1)->configs(i).getPos(e,pos2);
          doublevar d12=0,d22=0;
          doublevar rdiff2=0;
          doublevar dot=0;
          for(int d=0; d< 3; d++) { 
            drift1(d)=h->wfs(i,e).amp(0,d+1);
            drift2(d)=(h+1)->wfs(i,e).amp(0,d+1);
            d12+=drift1(d)*drift1(d);
            d22+=drift2(d)*drift2(d);
            rdiff2+=(pos1(d)-pos2(d))*(pos1(d)-pos2(d));
            dot+=(pos2(d)-pos1(d))*(drift2(d)-drift1(d));
          }
          logpi(i)-=0.25*timestep*(d12+d22)+rdiff2/(2*timestep)+0.5*dot;
        }
      }
    }
    
    
    Wf_return wfval(wf(i)->nfunc(), 2);
    wf(i)->getVal(mywfdata(i),0, wfval);
    if(pc_gf) { 
      h=pts(walker).past_energies.begin();
      logpi(i)+=h->wfs(i,0).amp(0,0);
      h=pts(walker).past_energies.end()-1;
      logpi(i)+=h->wfs(i,0).amp(0,0);      
    }
    else { 
      logpi(i)+=2*wfval.amp(0,0);
    }
    //cout  << setw(18) << wfval.amp(0,0) << setw(3) << wfval.sign(0);
  }
  
  Array1 <doublevar> totjacob(nsys);
  for(int i=0; i< nsys; i++) { 
    if(pc_gf) { 
      totjacob(i)=1;
      for(deque<Dmc_corr_history>::iterator h=pts(walker).past_energies.begin(); 
          h!= pts(walker).past_energies.end()-1; h++) { 
        totjacob(i)*=h->jacobian(i);
      }
    }
    else { 
      totjacob(i)=pts(walker).jacobian(i);
    }
  }
  
  //cout << endl;
  for(int i=0; i< nsys; i++) {
    double ratio=0;
    for(int j=0; j < nsys; j++) { 

        ratio+=totjacob(j)*exp(logpi(j)-logpi(i));
    }
    pts(walker).prop.weight(i)=pts(walker).weight* totjacob(i)/ratio;
    //pts(walker).prop.weight(i)=pts(walker).jacobian(i)/ratio;
  }
  //This is somewhat inaccurate..will need to change it later
  //For the moment, the autocorrelation will be slightly
  //underestimated
  pts(walker).prop.parent=walker;
  pts(walker).prop.nchildren=1;
  pts(walker).prop.children(0)=walker;
  pts(walker).prop.count=1;
     
    

  return avg_acceptance;
}


//----------------------------------------------------------------------


//-------------------------------Evaluate the new points
void Dmc_corr_method::add_point(int walker) { 
  int currsys=pts(walker).system;
  int nsys=mywfdata.GetDim(0);
  int nelectrons=mysys(currsys)->nelectrons(0)+mysys(currsys)->nelectrons(1);
  int nrandvar=mypseudo->nTest();
  //cout << "nrandvar " << nrandvar  << endl;
  Array1 <doublevar> rand_num(nrandvar);
  for(int i=0; i< nrandvar; i++) rand_num(i)=rng.ulec();
  
  Properties_point pt;
  pt.setSize(nsys, nsys, 1);
  for(int i=0; i< nsys; i++) { 
    Array1 <doublevar> kinetic(wf(i)->nfunc());
    Array1 <doublevar> nonlocal(wf(i)->nfunc());
    //cout << "sys " << i << endl;
    doublevar tot_jacobian=1;
    if(i!=currsys) {  //Warp the non-currsys systems
      doublevar jacobian=1;
      //cout << "warping " << endl;
      for(int e=0; e< nelectrons; e++) { 
        Array1 <doublevar> ref_pos(3);
        Array1 <doublevar> warp_pos(3);
        sample(currsys)->getElectronPos(e,ref_pos);
        //cout << "warping " << endl;
        warper.space_warp(sample(currsys), sample(i), e, ref_pos, warp_pos, jacobian);
        tot_jacobian*=jacobian;
        sample(i)->setElectronPos(e,warp_pos);
        //cout << "ref pos " << ref_pos(0) << " " << ref_pos(1) << " " << ref_pos(2)
        //<< " warp pos " << warp_pos(0) << " " << warp_pos(1) << " " << warp_pos(2)
        //<< endl;
      }
      //cout << "notifying and updating " << endl;
      wf(i)->notify(all_electrons_move,0);
      //cout << "update " << endl;
    }
    //cout << "here2 " << endl;
    wf(i)->updateLap(mywfdata(i), sample(i));

    mysys(i)->calcKinetic(mywfdata(i),sample(i), 
                          wf(i), kinetic);

    pt.kinetic(i)=kinetic(0);
    //cout << " blah " << pt.kinetic.GetDim(0) << " " << mysys.GetDim(0) 
    //<< " " << sample.GetDim(0) << endl;
    pt.potential(i)=mysys(i)->calcLoc(sample(i));
    //cout <<  " nonloc " << endl;
    mypseudo->calcNonlocWithTest(mywfdata(i), sample(i), wf(i),rand_num,nonlocal);
    pt.nonlocal(i)=nonlocal(0);
    pts(walker).jacobian(i)=tot_jacobian;
  }
  //cout << "there " << endl;
  //------------------------------Now store the information for the walker
  
  Dmc_corr_history new_hist;
  new_hist.jacobian=pts(walker).jacobian;
  new_hist.main_en.Resize(nsys);
  for(int i=0; i< nsys; i++) { 
    new_hist.main_en(i)=pt.energy(i);
  }
  
  new_hist.wfs.Resize(nsys, nelectrons);
  //This is technically getLaping twice.  Might be slightly inefficient.
  //Of course, this could be rectified by making getLap a bit smarter..
  for(int i=0; i< nsys; i++) { 
    for(int j=0; j< nelectrons; j++) { 
      new_hist.wfs(i,j).Resize(1,6);
      wf(i)->getLap(mywfdata(i), j, new_hist.wfs(i,j));
    }
  }
  new_hist.configs.Resize(nsys);
  for(int i=0; i< nsys; i++) { 
    new_hist.configs(i).savePos(sample(i));
  }
  
  pts(walker).past_energies.push_front(new_hist);
  deque<Dmc_corr_history> & past(pts(walker).past_energies);
  if(past.size() > nhist) 
    past.erase(past.begin()+nhist, past.end());
  
  pts(walker).prop=pt;
}  

//----------------------------------------------------------------------
void Dmc_corr_method::updateEtrial() {
  double feed=1.0;
  for(int i=0; i< eref.GetDim(0); i++) { 
    doublevar totweight=0;
    int nc=0;
    
    for(int walker=0; walker < nconfig; walker++) { 
      if(pts(walker).system==i) { 
        totweight+=pts(walker).weight;
        nc++;
      }
    }
    totweight=parallel_sum(totweight);
    nc=parallel_sum(nc);
    etrial(i)=eref(i)-feed*log(totweight/double(nc));
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


int Dmc_corr_method::calcBranch() { 
  int totwalkers=mpi_info.nprocs*nconfig;
  Array1 <doublevar> weights(totwalkers);
  Array1 <doublevar> my_weights(nconfig);
  Array1 <int> systems(totwalkers);
  Array1 <int> my_systems(totwalkers);
  for(int walker=0; walker < nconfig; walker++) { 
    my_weights(walker)=pts(walker).weight;
    my_systems(walker)=pts(walker).system;
  }
#ifdef USE_MPI
  MPI_Allgather(my_weights.v,nconfig, MPI_DOUBLE, weights.v,nconfig,MPI_DOUBLE, MPI_Comm_grp);
  MPI_Allgather(my_systems.v,nconfig, MPI_INT, systems.v,nconfig,MPI_INT, MPI_Comm_grp);
#else
  weights=my_weights;
#endif
  Array1 <int> my_branch(nconfig);
  Array1 <int> nwalkers(mpi_info.nprocs);
  nwalkers=0;
  int root=0;
  //if(mpi_info.node==root) {
  //  for(int i=0; i< totwalkers; i++) { 
  //    cout << "walker " << i << " " << weights(i) << endl;
  //  }
  //}
  
  if(mpi_info.node==root) {  //this if/else clause may be refactored out
    
    Array1 <int> branch(totwalkers);
    //----Find which walkers branch/die
    //we do it on one node since otherwise we'll have different random numbers!
    //we'll assign the weight for each copy that will be produced
    //this is the core of the branching algorithm..
    //my homegrown algo, based on Umrigar, Nightingale, and Runge
    branch=-1;
    const doublevar split_threshold=1.1;
    for(int w=0; w< totwalkers; w++) { 
      if(weights(w) > split_threshold && branch(w)==-1) { 
        //find branching partner
        doublevar smallestwt=100;
        int smallest=-1;
        for(int w2=0; w2 < totwalkers; w2++) { 
          if(branch(w2)==-1 && w2!= w && weights(w2) < smallestwt
             && systems(w)==systems(w2)) { //restrict branching to the same system
            smallest=w2;
            smallestwt=weights(w2);
          }
        }
        //cout << "pairing " << w << " " << smallest << " systems " << systems(w)
        //     << " " << systems(smallest) << " weights " << weights(w) << "  "
        //<< weights(smallest) << endl;
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
  Array1 <Dmc_corr_point> savepts=pts;
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
#ifdef USE_MPI
  //Finally, send or receive spillover walkers 
  if(nnwalkers > nconfig) { 
    vector<Queue_element>::iterator queue_pos=send_queue.begin();
    while(curr < nconfig) { 
      if(my_branch(curr) > 0) { 
        my_branch(curr)--;
        while(queue_pos->from_node != mpi_info.node) { 
          queue_pos++;
        }
        cout << mpi_info.node << ":curr " << curr << " my_branch " << my_branch(curr) << endl;
        cout << mpi_info.node << ":sending " << queue_pos->from_node << " to " << queue_pos->to_node << endl;
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
      cout << mpi_info.node <<":receiving from " << queue_pos->from_node << " to " << curr_copy << endl;
      pts(curr_copy).mpiReceive(queue_pos->from_node);
      //pts(curr_copy).weight=1;
      curr_copy++;
      queue_pos++;
    }
  }
#endif
  return killsize;
  //exit(0);
}



//----------------------------------------------------------------------

void Dmc_corr_point::mpiSend(int node) { 
#ifdef USE_MPI
  prop.mpiSend(node);
  config_pos.mpiSend(node);
  int n=past_energies.size();
  MPI_Send(&n,1, MPI_INT,node,0,MPI_Comm_grp);
  for(deque<Dmc_corr_history>::iterator i=past_energies.begin();
      i!= past_energies.end(); i++) { 
    i->mpiSend(node);
  }
  
  MPI_Send(&weight,1,MPI_DOUBLE,node,0,MPI_Comm_grp);  
  int nelectrons=age.GetDim(0);
  MPI_Send(&nelectrons,1, MPI_INT,node,0,MPI_Comm_grp);
  MPI_Send(age.v,nelectrons, MPI_DOUBLE, node,0,MPI_Comm_grp);
  MPI_Send(&system,1, MPI_INT,node,0,MPI_Comm_grp);
  
  int njacob=jacobian.GetDim(0);
  MPI_Send(&njacob, 1, MPI_INT, node, 0, MPI_Comm_grp);
  MPI_Send(jacobian.v, njacob, MPI_DOUBLE, node, 0, MPI_Comm_grp);
  MPI_Send(initial_sign.v, njacob, MPI_INT, node, 0, MPI_Comm_grp);
#endif
}

void Dmc_corr_point::mpiReceive(int node) {
#ifdef USE_MPI
  prop.mpiReceive(node);
  config_pos.mpiReceive(node);
  int n;
  MPI_Status status;
  MPI_Recv(&n,1,MPI_INT,node,0,MPI_Comm_grp, &status);
  Dmc_corr_history tmp_hist;
  past_energies.clear();
  for(int i=0; i< n; i++) {
    tmp_hist.mpiReceive(node);
    past_energies.push_back(tmp_hist);
  }  
  
  MPI_Recv(&weight,1,MPI_DOUBLE,node,0,MPI_Comm_grp, &status);
  
  int nelectrons;
  MPI_Recv(&nelectrons,1,MPI_INT,node,0,MPI_Comm_grp,&status);
  age.Resize(nelectrons);
  MPI_Recv(age.v,nelectrons,MPI_DOUBLE, node,0,MPI_Comm_grp,&status);
  MPI_Recv(&system,1,MPI_INT,node,0,MPI_Comm_grp,&status);
  int njacob;
  MPI_Recv(&njacob,1,MPI_INT, node, 0, MPI_Comm_grp, &status);
  jacobian.Resize(njacob);
  MPI_Recv(jacobian.v, njacob, MPI_DOUBLE, node, 0, MPI_Comm_grp, &status);
  initial_sign.Resize(njacob);
  MPI_Recv(initial_sign.v, njacob, MPI_INT, node, 0, MPI_Comm_grp, &status);
#endif
}

void Dmc_corr_history::mpiSend(int node) { 
#ifdef USE_MPI
  int nen=main_en.GetDim(0);
  MPI_Send(&nen, 1,MPI_INT, node,0,MPI_Comm_grp);
  MPI_Send(main_en.v,nen,MPI_DOUBLE,node,0,MPI_Comm_grp);
#endif
}

void Dmc_corr_history::mpiReceive(int node) { 
#ifdef USE_MPI
  MPI_Status status;
  int nen;
  MPI_Recv(&nen, 1,MPI_INT, node, 0, MPI_Comm_grp, & status);
  main_en.Resize(nen);
  MPI_Recv(main_en.v,nen,MPI_DOUBLE,node,0,MPI_Comm_grp,&status);
#endif
}



//----------------------------------------------------------------------


