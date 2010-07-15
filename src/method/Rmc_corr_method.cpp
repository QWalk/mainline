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



#include "Rmc_corr_method.h"
#include "qmc_io.h"
#include "ulec.h"
#include "Program_options.h"
#include "average.h"


void Rmc_corr_method::read(vector <string> words,
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

  
  if(!readvalue(words, pos=0, nhist, "LENGTH")) 
    nhist=50;


  //optional options


  if(!readvalue(words, pos=0,npsteps, "PARTIAL_STEPS")) 
    npsteps=5;

  readvalue(words, pos=0, storeconfig, "STORECONFIG");
  readvalue(words, pos=0, readconfig, "READCONFIG");

  if(!readvalue(words, pos=0, log_label, "LABEL"))
    log_label="rmc_corr";
  
  vector <string> warpsec;
  if(readsection(words, pos=0, warpsec,"WARPER")) {
    warper.read(warpsec);
  }
  
  if(haskeyword(words, pos=0, "NOWARP")) 
    warper.set_warp(0);
  
  pc_gf=0;
  if(haskeyword(words, pos=0, "PC_GF")) pc_gf=1;
  
  dmc_gf=0;
  if(haskeyword(words, pos=0, "DMC_GF")) { 
    pc_gf=0;
    dmc_gf=1;
  }
  
  
  vector <string> avg_gen;
  pos=0;
  while(readsection(words, pos, avg_gen, "AVERAGE"))
    avgwords.push_back(avg_gen);
  pos=0;
  while(readsection(words, pos, avg_gen, "DENSITY"))
    denswords.push_back(avg_gen);
  
  vector <string> dynamics_words;
  if(!readsection(words, pos=0, dynamics_words, "DYNAMICS") ) 
    dynamics_words.push_back("SPLIT");

  allocate(dynamics_words, dyngen);
  dyngen->enforceNodes(1);


}

//----------------------------------------------------------------------

int Rmc_corr_method::generateVariables(Program_options & options) {

  if(!have_read_options) 
    error("need to call Rmc_corr_method::read before generateVariables");
  if(have_allocated_variables) 
    error("already allocated variables in Rmc_corr_method");

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




int Rmc_corr_method::showinfo(ostream & os) {
  int nsys=mysys.GetDim(0);
  for(int i=0; i< nsys; i++) { 
    mysys(i)->showinfo(os);
    mywfdata(i)->showinfo(os);
  }
    
  mypseudo->showinfo(os);
  os << "###########################################################\n";
  os << "Reptation Monte Carlo:\n";
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

void Rmc_corr_method::run(Program_options & options, ostream & output) {
  if(!have_allocated_variables) 
    error("Must generate variables to use Rmc_corr_method::run");
  string logfile=options.runid+".log";
  int nsys=mywfdata.GetDim(0);

  if(mpi_info.node==0 ) {
    ofstream logout(logfile.c_str(), ios::app);
    logout << "#-------------------------------------------------\n";
    logout << "#RMC run: timestep " << timestep 
           << endl;
    logout << "#-------------------------------------------------\n\n\n";
    logout.close();
  }
  
  //-------------Properties generators
  int navg=avgwords.size();
  avggen.Resize(nsys, navg);
  for(int sys=0; sys < nsys; sys++) {
    for(int i=0; i< navg; i++) 
      allocate(avgwords[i],mysys(sys),mywfdata(sys),avggen(sys,i));
  }
  int ndens=denswords.size();
  local_dens.Resize(nsys,ndens);
  for(int sys=0; sys< nsys; sys++) { 
    string runidsys=options.runid;
    append_number(runidsys, sys);
    for(int i=0; i< ndens; i++) 
      allocate(denswords[i], mysys(sys), options.runid, local_dens(sys,i));
  }
  
  
  myprop.setLog(logfile, log_label);
  myprop.initializeLog(avggen);
  string center_label=log_label+"_center";
  centerprop.setLog(logfile, center_label);
  centerprop.initializeLog(avggen);
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
  if(nconfig%nsys!=0) { 
    error("nconfig must be an even multiple of the number of systems");
  }
  for(int i=0; i< nconfig; i++) { 
    pts(i).system= i%nsys;
  }
  myprop.setSize(nsys, nblock, nstep+1, nconfig, mysys(0),mywfdata(0));
  centerprop.setSize(nsys, nblock,nstep+1,nconfig, mysys(0),mywfdata(0));
  
  ///-----------------Generate VMC configs if we don't have a checkfile
  if(readconfig=="") { 
    int nstep_vmc=50;
    doublevar timestep_vmc=1.0;
    for(int walker=0; walker < nconfig; walker++) { 
      int currsys=pts(walker).system;
      sample(currsys)->randomGuess();
      wf(currsys)->updateLap(mywfdata(currsys), sample(currsys));
      for(int step=0; step < nstep_vmc; step++) { 
        Dynamics_info dinfo;
        for(int e=0; e< nelectrons; e++) {
          dyngen->sample(e,sample(currsys), wf(currsys),
                         mywfdata(currsys), guidingwf, dinfo, timestep_vmc);
        }
      }
      pts(walker).path.push_back(add_point(walker));
    }
    
    //Generate first reptile  
    //cout << mpi_info.node << ":generating reptile " << endl;
    for(int walker=0; walker < nconfig; walker++) { 
      //cout << "size " << pts(walker).path.size() << endl;  
      int currsiz=pts(walker).path.size();
      int currsys=pts(walker).system;
      pts(walker).path[currsiz-1].configs(currsys).restorePos(sample(currsys));
      wf(currsys)->updateLap(mywfdata(currsys), sample(currsys));
      while(pts(walker).path.size() < nhist) { 
        Dynamics_info dinfo;
        for(int e=0; e< nelectrons; e++) {
          dyngen->sample(e,sample(currsys), wf(currsys),
                         mywfdata(currsys), guidingwf, dinfo, timestep);
        }
        pts(walker).path.push_back(add_point(walker));
      }
      //cout << "final size " << pts(walker).path.size() << endl;
    }
  }
  else readcheck();
  
  //-------------------------calculate the green's functions for the updates
  for(int walker=0; walker < nconfig; walker++) {
    pts(walker).branching.Resize(nsys);
    pts(walker).totgf.Resize(nsys);
    pts(walker).totgf=0; pts(walker).branching=0;
    for(int i=0; i< nsys; i++) { 
      for(deque<Rmc_corr_history>::iterator h=pts(walker).path.begin(); 
          h!= pts(walker).path.end()-1; h++) { 
        //cout << "h " << endl;
        doublevar btmp;
        pts(walker).totgf(i)+=get_green_weight(*h,*(h+1),i,btmp);
        pts(walker).branching(i)+=btmp;
      }
    }
    pts(walker).recalc_gf=0;
  }
  Array1 <Sample_point *> center_samp(nsys);
  center_samp=NULL;
  for(int i=0; i< nsys; i++) 
    mysys(i)->generateSample(center_samp(i));
  
  //------------------------------------------
  npsteps=nstep;
  for(int block=0; block < nblock; block++) { 
    int accept=0;
    //for(int step=0; step < nstep; ) { 
      //cout << "block " << block << " step " << step << endl;
      Dynamics_info dinfo;
      
      doublevar avg_acceptance=0;  
      //int n_partial=min(npsteps,nstep-step);
      int ncross=0;
      for(int walker=0; walker < nconfig; walker++) {
        int currsys=pts(walker).system;
        if(pts(walker).direction==1) 
          (pts(walker).path.end()-1)->configs(currsys).restorePos(sample(currsys));
        else 
          pts(walker).path.begin()->configs(currsys).restorePos(sample(currsys));
        
        
        for(int step=0; step < nstep; step++) { 
          mypseudo->randomize();            
          
          
          accept+=propagate_walker(walker);
          
          Properties_point pt=pts(walker).path.begin()->prop;
          pt.parent=walker;
          pt.nchildren=1;
          pt.children(0)=walker;
          pt.count=1;
          pt.weight=pts(walker).weight;
          cout.precision(15);
          //cout << "adding point " << pt.avgrets(0,0).vals[2] << " " << pt.avgrets(1,0).vals[2] << endl;
          myprop.insertPoint(step, walker, pt);
          //cout << "here " << endl;
          int half=pts(walker).path.size()/2;
          //cout << "half " << half << endl;
          pt=pts(walker).path[half].prop;
          pt.parent=walker; pt.children=1; pt.children(0)=walker;pt.count=1;
          pt.weight=pts(walker).weight;
          centerprop.insertPoint(step,walker,pt);
          if(local_dens.GetDim(1) > 0) { 
            for(int sys=0; sys < nsys; sys++) { 
              pts(walker).path[half].configs(sys).restorePos(center_samp(sys));
              for(int i=0; i< local_dens.GetDim(1); i++) { 
                local_dens(sys,i)->accumulate(center_samp(sys), pts(walker).weight(sys));
              }
            }
          }
                                           
        }

      }
      if(ncross > 0) { 
        cout << "proportion crossed " << double(ncross)/(nsys*nconfig*nstep) 
        << endl;
      }
    
      //step+=n_partial;
    
    storecheck();
    myprop.endBlock();
    centerprop.endBlock();
    doublevar acceptance=doublevar(parallel_sum(accept))/doublevar(parallel_sum(nconfig*nstep));
    Properties_final_average finavg;
    myprop.getFinal(finavg);
    
    if(output) { 
      output << "acceptance " << acceptance << 
        " average steps before a bounce " << 1.0/(1-acceptance) << endl;
      myprop.printBlockSummary(output);
    }
  }
  
  //loop over blocks 
  //loop over systems1
  
  deallocateIntermediateVariables();
}


//----------------------------------------------------------------------

doublevar Rmc_corr_method::get_green_weight(//deque <Rmc_corr_history>::iterator a, 
                                            //deque <Rmc_corr_history>::iterator b,
                                            Rmc_corr_history & a,
                                            Rmc_corr_history & b,
                                            int i, doublevar & branching) {
  
  branching=0;
  branching-= 0.5*timestep*(a.main_en(i)+b.main_en(i) );
  int nelectrons=mysys(i)->nelectrons(0)+mysys(i)->nelectrons(1);

  
  doublevar logpi=0;
  if(pc_gf) { //---------
    Array1 <doublevar> pos1(3), pos2(3), drift1(3), drift2(3);
    for(int e=0; e< nelectrons; e++) { 
      a.configs(i).getPos(e,pos1);
      b.configs(i).getPos(e,pos2);
      doublevar d12=0,d22=0;
      doublevar rdiff2=0;
      doublevar dot=0;
      for(int d=0; d< 3; d++) { 
        drift1(d)=a.wfs(i,e).amp(0,d+1);
        drift2(d)=b.wfs(i,e).amp(0,d+1);
        d12+=drift1(d)*drift1(d);
        d22+=drift2(d)*drift2(d);
        rdiff2+=(pos1(d)-pos2(d))*(pos1(d)-pos2(d));
        dot+=(pos2(d)-pos1(d))*(drift2(d)-drift1(d));
      }
      logpi-=0.25*timestep*(d12+d22)+rdiff2/(2*timestep)+0.5*dot;
    }
  }
  else if(dmc_gf) { 
    Array1 <doublevar> pos1(3), pos2(3), drift1(3), drift2(3);
    for(int e=0; e< nelectrons; e++) { 
      a.configs(i).getPos(e,pos1);
      b.configs(i).getPos(e,pos2);
      double atob=0, btoa=0;
      for(int d=0; d< 3; d++) { 
        drift1(d)=a.wfs(i,e).amp(0,d+1);
        drift2(d)=b.wfs(i,e).amp(0,d+1);
        btoa+=(pos1(d)-pos2(d)-drift2(d)*timestep)*(pos1(d)-pos2(d)-drift2(d)*timestep);
        atob+=(pos2(d)-pos1(d)-drift1(d)*timestep)*(pos2(d)-pos1(d)-drift1(d)*timestep);
      }
      //doublevar prob=1;
      logpi-=atob/(2*timestep);
    }
    
  }
  return logpi+branching;
}

//Move an walker by one DMC step, updating the pts(walker) variable and leaving
//the samples and wave functions at the new position.
//it's assumed that at least sample(pts(walker).system) is updated for that walker.

doublevar Rmc_corr_method::propagate_walker(int walker) { 
  int currsys=pts(walker).system;
  int nsys=mywfdata.GetDim(0);
  Dynamics_info dinfo;
  doublevar avg_acceptance=0;
  int nelectrons=mysys(currsys)->nelectrons(0)+mysys(currsys)->nelectrons(1);
  for(int e=0; e< nelectrons; e++) {
    int acc;
    acc=dyngen->sample(e, sample(currsys), wf(currsys),
                       mywfdata(currsys), guidingwf,dinfo, timestep);
    avg_acceptance+=dinfo.acceptance/(nelectrons*npsteps);
  }
  Rmc_corr_history newpt=add_point(walker);
  
  
  
  Array1 <doublevar> delbranch(nsys);
  Array1 <doublevar> deltot(nsys);
  delbranch=0;
  deltot=0;
  
  if(pts(walker).direction==1) { 
    for(int i=0; i < nsys; i++) { 
      doublevar branchupdate=0;
      deltot(i)+=get_green_weight(*(pts(walker).path.end()-1), 
                                  newpt,
                                  i, branchupdate);
      delbranch(i)+=branchupdate;
      deltot(i)+=get_green_weight(*pts(walker).path.begin(),
                                  *(pts(walker).path.begin()+1),i,branchupdate);
      delbranch(i)-=branchupdate;
    }
    
  }
  else  { 
    for(int i=0; i < nsys; i++) { 
      doublevar branchupdate=0;
      deltot(i)+=get_green_weight(newpt,
                                  *pts(walker).path.begin(),
                                  i, branchupdate);
      delbranch(i)+=branchupdate;
      deltot(i)+=get_green_weight(*(pts(walker).path.end()-1),
                                  *(pts(walker).path.end()-2),i,branchupdate);
      
      delbranch(i)-=branchupdate;
    }
  }
  int accepted=0;

  if(int(exp(delbranch(currsys))+rng.ulec())) { 
    
    accepted=1;
    if(pts(walker).direction==1) { 
      pts(walker).path.push_back(newpt);
      pts(walker).path.pop_front();
    }
    else { 
      pts(walker).path.push_front(newpt);
      pts(walker).path.pop_back();
    }
    
    for(int i=0; i < nsys; i++) { 
      pts(walker).branching(i)+=delbranch(i);
      pts(walker).totgf(i)+=deltot(i);
    }
  }
  else { 
    
    pts(walker).direction*=-1;
    if(pts(walker).direction==1) 
      (pts(walker).path.end()-1)->configs(currsys).restorePos(sample(currsys));
    else 
      pts(walker).path.begin()->configs(currsys).restorePos(sample(currsys));    
    
  }
  
  
  Array1 <doublevar> logpi(nsys); //the logarithm of the pdf for each system
  logpi=pts(walker).totgf;
     
  
  pts(walker).weight=exp(pts(walker).branching(currsys));
  deque<Rmc_corr_history>::iterator h;
  for(int i=0; i< nsys; i++) { 
    
    Wf_return wfval(wf(i)->nfunc(), 2);
    wf(i)->getVal(mywfdata(i),0, wfval);
    if(pc_gf || dmc_gf ) { 
      h=pts(walker).path.begin();
      logpi(i)+=h->wfs(i,0).amp(0,0);
      h=pts(walker).path.end()-1;
      logpi(i)+=h->wfs(i,0).amp(0,0);      
    }
    else { 
      logpi(i)+=2*wfval.amp(0,0);
    }
    //cout  << setw(18) << wfval.amp(0,0) << setw(3) << wfval.sign(0);
  }
  
  Array1 <doublevar> totjacob(nsys);
  for(int i=0; i< nsys; i++) { 
    if(pc_gf || dmc_gf) { 
      totjacob(i)=1;
      for(deque<Rmc_corr_history>::iterator h=pts(walker).path.begin(); 
          h!= pts(walker).path.end()-1; h++) { 
        totjacob(i)*=h->jacobian(i);
      }
    }
    else { 
      totjacob(i)=pts(walker).path.begin()->jacobian(i);
    }
  }
  
  //cout << endl;
  pts(walker).weight.Resize(nsys);
  for(int i=0; i< nsys; i++) {
    double ratio=0;
    for(int j=0; j < nsys; j++) { 
        ratio+=totjacob(j)*exp(logpi(j)-logpi(i));
    }
    //pts(walker).weight(i)=pts(walker).weight* totjacob(i)/ratio;
    pts(walker).weight(i)=totjacob(i)/ratio;
  }
  
  return accepted;
}


//----------------------------------------------------------------------


//-------------------------------Evaluate a new point
Rmc_corr_history Rmc_corr_method::add_point(int walker) { 
  Rmc_corr_history new_hist;

  int currsys=pts(walker).system;
  int nsys=mywfdata.GetDim(0);
  int nelectrons=mysys(currsys)->nelectrons(0)+mysys(currsys)->nelectrons(1);
  int nrandvar=mypseudo->nTest();
  //cout << "nrandvar " << nrandvar  << endl;
  Array1 <doublevar> rand_num(nrandvar);
  for(int i=0; i< nrandvar; i++) rand_num(i)=rng.ulec();
  Array1 <doublevar> jacobian_save(nsys);
  new_hist.prop.setSize(nsys);
  new_hist.prop.avgrets.Resize(avggen.GetDim(0),avggen.GetDim(1));
  for(int i=0; i< nsys; i++) { 
    Array1 <doublevar> kinetic(wf(i)->nfunc());
    Array1 <doublevar> nonlocal(wf(i)->nfunc());
    doublevar tot_jacobian=1;
    if(i!=currsys) {  //Warp the non-currsys systems
      doublevar jacobian=1;
      for(int e=0; e< nelectrons; e++) { 
        Array1 <doublevar> ref_pos(3);
        Array1 <doublevar> warp_pos(3);
        sample(currsys)->getElectronPos(e,ref_pos);
        warper.space_warp(sample(currsys), sample(i), e, ref_pos, warp_pos, jacobian);
        tot_jacobian*=jacobian;
        sample(i)->setElectronPos(e,warp_pos);
      }
      wf(i)->notify(all_electrons_move,0);
    }
    //cout << "here2 " << endl;
    wf(i)->updateLap(mywfdata(i), sample(i));

    mysys(i)->calcKinetic(mywfdata(i),sample(i), 
                          wf(i), kinetic);

    new_hist.prop.kinetic(i)=kinetic(0);
    new_hist.prop.potential(i)=mysys(i)->calcLoc(sample(i));
    mypseudo->calcNonlocWithTest(mywfdata(i), mysys(i), sample(i), wf(i),rand_num,nonlocal);
    new_hist.prop.nonlocal(i)=nonlocal(0);
    jacobian_save(i)=tot_jacobian;
    for(int a=0; a < avggen.GetDim(1); a++) { 
      avggen(i,a)->evaluate(mywfdata(i),wf(i),mysys(i),sample(i),new_hist.prop.avgrets(i,a));
    }
  }
  //------------------------------Now store the information for the walker
  
  new_hist.jacobian=jacobian_save;
  new_hist.main_en.Resize(nsys);
  for(int i=0; i< nsys; i++) { 
    new_hist.main_en(i)=new_hist.prop.energy(i);
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
  
  return new_hist;
  //pts(walker).path.push_front(new_hist);
  
  //pts(walker).prop=pt;
}  

//----------------------------------------------------------------------

void Rmc_corr_method::storecheck() { 
  if(storeconfig=="") return;
  write_configurations(storeconfig, pts);
}

//----------------------------------------------------------------------
void Rmc_corr_method::readcheck() { 
  if(readconfig=="") return;
  read_configurations(readconfig, pts);
  if(pts.GetDim(0)!=nconfig) error("Changed number of configurations in RMC_CORR");
  for(int w=0; w < pts.GetDim(0); w++) { 
    int sys=pts(w).system;
    for(deque<Rmc_corr_history>::iterator i=pts(w).path.begin(); 
        i!=pts(w).path.end(); i++) { 
      i->configs(sys).restorePos(sample(sys));
      *i=add_point(w);
    }
  }
}

//----------------------------------------------------------------------

void Rmc_corr_history::write(ostream & os ) { 
  //We'll write only the configurations and expect the rest to 
  //be recalculated.
  os << "nsys " << configs.GetDim(0) << endl;
  for(int i=0; i< configs.GetDim(0); i++) { 
    configs(i).write(os);
  }
}
//----------------------------------------------------------------------
void Rmc_corr_history::read(istream & is) { 
  string dummy;
  int nsys;
  is >> dummy >> nsys;
  configs.Resize(nsys);
  for(int i=0; i< nsys; i++) configs(i).read(is);
}
//----------------------------------------------------------------------
void Rmc_corr_point::write(ostream & os) { 
  os << "direction " << direction << endl;
  os << "system " << system << endl;
  os << "path_length " << path.size() << endl;
  for(deque<Rmc_corr_history>::iterator i=path.begin(); i!=path.end(); 
      i++) { 
    i->write(os);
  }
}
//----------------------------------------------------------------------
void Rmc_corr_point::read(istream & is) { 
  string dummy;
  is >> dummy >> direction;
  is >> dummy >> system; 
  int path_length;
  is >> dummy >> path_length;
  path.resize(path_length);
  for(deque<Rmc_corr_history>::iterator i=path.begin(); i!=path.end(); 
      i++) i->read(is);
}
//----------------------------------------------------------------------


void Rmc_corr_history::mpiSend(int node) { 
  int nsys=configs.GetDim(0);
  MPI_Send(nsys, node);
  for(int i=0; i< nsys; i++) 
    configs(i).mpiSend(node);
}


//----------------------------------------------------------------------

void Rmc_corr_history::mpiReceive(int node) { 
  int nsys;
  MPI_Recv(nsys, node);
  configs.Resize(nsys);
  for(int i=0; i < nsys; i++) 
    configs(i).mpiReceive(node);
}


//----------------------------------------------------------------------


void Rmc_corr_point::mpiSend(int node) { 
  int npath=path.size();
  MPI_Send(direction, node);
  MPI_Send(system, node);
  MPI_Send(npath, node);
  for(deque<Rmc_corr_history>::iterator i=path.begin(); i!=path.end(); i++)
    i->mpiSend(node);
  
}

//----------------------------------------------------------------------


void Rmc_corr_point::mpiReceive(int node) { 
  int npath;
  MPI_Recv(direction, node);
  MPI_Recv(system, node);
  MPI_Recv(npath, node);
  path.resize(npath);
  for(deque<Rmc_corr_history>::iterator i=path.begin(); i!=path.end(); i++)
    i->mpiReceive(node);
}

//----------------------------------------------------------------------
