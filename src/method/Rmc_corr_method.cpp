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

  if(readvalue(words, pos=0, storeconfig, "STORECONFIG"))
    canonical_filename(storeconfig);
  
  if(readvalue(words, pos=0, readconfig, "READCONFIG"))
    canonical_filename(readconfig);

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

void Rmc_corr_method::run(Program_options & options, ostream & output) {
  if(!have_allocated_variables) 
    error("Must generate variables to use Rmc_corr_method::run");
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
  if(nconfig%nsys!=0) { 
    error("nconfig must be an even multiple of the number of systems");
  }
  for(int i=0; i< nconfig; i++) { 
    pts(i).system= i%nsys;
  }
  myprop.setSize(nsys, nblock, nstep+1, nconfig, mysys(0), 
                 mywfdata(0));
  
  ///-----------------Generate VMC configs
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
      cout << "finished walker " << walker << " sysy " << currsys << endl;
      pts(walker).path.push_back(add_point(walker));
    }
    
    //Generate first reptile  
    cout << mpi_info.node << ":generating reptile " << endl;
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
  
  
  //------------------------------------------
  npsteps=nstep;
  for(int block=0; block < nblock; block++) { 
    int accept=0;
    for(int step=0; step < nstep; ) { 
      //cout << "block " << block << " step " << step << endl;
      Dynamics_info dinfo;
      
      doublevar avg_acceptance=0;  
      int n_partial=min(npsteps,nstep-step);
      //cout << mpi_info.node << ":move" << endl;
      int ncross=0;
      for(int walker=0; walker < nconfig; walker++) {
        int currsys=pts(walker).system;
        if(pts(walker).direction==1) 
          (pts(walker).path.end()-1)->configs(currsys).restorePos(sample(currsys));
        else 
          pts(walker).path.begin()->configs(currsys).restorePos(sample(currsys));
        
        
        for(int p=0; p < n_partial; p++) { 
          mypseudo->randomize();            
          
          
          accept+=propagate_walker(walker);
          
          Properties_point pt=pts(walker).path.begin()->prop;
          pt.parent=walker;
          pt.nchildren=1;
          pt.children(0)=walker;
          pt.count=1;
          pt.weight=pts(walker).weight;
          myprop.insertPoint(step+p, walker, pt);
        }
        //pts(walker).config_pos.savePos(sample(currsys));

      }
      //cout << "finished partial " << endl;
      if(ncross > 0) { 
        cout << "proportion crossed " << double(ncross)/(nsys*nconfig*n_partial) 
        << endl;
      }
    
      step+=n_partial;
    }
    myprop.endBlock();
    storecheck();
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

    new_hist.prop.kinetic(i)=kinetic(0);
    //cout << " blah " << pt.kinetic.GetDim(0) << " " << mysys.GetDim(0) 
    //<< " " << sample.GetDim(0) << endl;
    new_hist.prop.potential(i)=mysys(i)->calcLoc(sample(i));
    //cout <<  " nonloc " << endl;
    mypseudo->calcNonlocWithTest(mywfdata(i), mysys(i), sample(i), wf(i),rand_num,nonlocal);
    new_hist.prop.nonlocal(i)=nonlocal(0);
    jacobian_save(i)=tot_jacobian;
  }
  //cout << "there " << endl;
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
  ofstream os(storeconfig.c_str());
  long int is1, is2;
  rng.getseed(is1,is2);
  os << "RANDNUM " << is1 << " " << is2 << endl;
  assert(pts.GetDim(0)==nconfig);
  for(int w=0; w< nconfig; w++) { 
    os << "REPTILE { " << endl;
    os << "SYSTEM " << pts(w).system << endl;
    os << "DIRECTION " << pts(w).direction << endl;
    os << "PATH_SIZE " << pts(w).path.size() << endl;
    for(deque<Rmc_corr_history>::iterator i=pts(w).path.begin(); 
        i!=pts(w).path.end(); i++) { 
      os << "SAMPLE_POINT { \n";
      i->configs(pts(w).system).write(os);
      os << " } \n";      
    }
    os << " } \n";
  }
}

void Rmc_corr_method::readcheck() { 
  if(readconfig=="") return;
  ifstream is(readconfig.c_str());
  if(!is) {
    error("Could not open ",readconfig);
  }
  long int is1, is2;
  string dummy;
  is >> dummy >> is1 >> is2;
  if(dummy!="RANDNUM") error("expected RANDNUM");
  rng.seed(is1,is2);
  pts.Resize(nconfig);
  Config_save_point tmpsave;

  for(int w=0; w< nconfig; w++) { 
    is >> dummy >> dummy;
    is >> dummy >> pts(w).system;
    is >>  dummy >> pts(w).direction;
    if(dummy!="DIRECTION") error("expected DIRECTION");
    int path_size;
    is >> dummy >> path_size;
    if(path_size != nhist) error("Can't restart with a different path length");
    int currsys=pts(w).system;
    for(int i=0; i< nhist; i++) { 
      //cout << "walker " << w << " hist " << i << endl;
      is >> dummy ; 
      if(dummy!="SAMPLE_POINT") error("expected SAMPLE_POINT");
      is >> dummy; 
      tmpsave.read(is);
      is >> dummy;
      tmpsave.restorePos(sample(currsys));
      wf(currsys)->updateLap(mywfdata(currsys),sample(currsys));
      pts(w).path.push_back(add_point(w));
    }
    is >> dummy; //last end bracket
  }
  
}


