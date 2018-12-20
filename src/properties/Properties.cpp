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
#include "Properties.h"

#include "qmc_io.h"
#include "Wavefunction_data.h"
#include "ulec.h"
#include <iomanip>

//######################################################################

void allocate(vector<string> & words, System * sys, string & runid,
	 Local_density_accumulator *& denspt) { 
  if(words.size() < 1) error("empty density section");
  if(caseless_eq(words[0], "DENSITY"))
    denspt=new One_particle_density;
  else if(caseless_eq(words[0], "POTENTIAL"))
    denspt=new Local_potential_density;
  else
    error("Didn't understand density keyword",words[0]);
  denspt->init(words, sys, runid);
}

void allocate(vector<string> & words, System * sys, string & runid,
	 Nonlocal_density_accumulator *& nldenspt) { 
  if(words.size() < 1) error("empty non-local density section");
  //  if(caseless_eq(words[0],"OBDM"))
  //    nldenspt=new OBDM;
  //  else if(caseless_eq(words[0],"TBDM"))
  //    nldenspt=new TBDM;
  //  else
  error("Didn't understand non-local density keyword",words[0]);
  nldenspt->init(words, sys, runid);
}

//######################################################################


void One_particle_density::init(vector<string> & words, System * sys,
                                string & runid) {
  int ndim=sys->ndim();
  norm=1;
  outputfile=runid+".cube";
  nup=sys->nelectrons(0);
  unsigned int pos=0;
  start_electron=0;
  end_electron=sys->nelectrons(0)+sys->nelectrons(1);
  if(haskeyword(words, pos=0, "UP")) { 
    end_electron=sys->nelectrons(0);
    outputfile=runid+".up.cube";
  }
  
  if(haskeyword(words,pos=0, "DOWN")) { 
    start_electron=sys->nelectrons(0);
    outputfile=runid+".down.cube";
  }

  if(end_electron==start_electron)
    error("In one-particle density, UP and DOWN are incompatible");
  
  
  readvalue(words, pos=0, outputfile,"OUTPUTFILE");
  
  //the domain is chosen as follows:
  //1. If MIN or MAX is input, set the respective quantity to user's values
  //2. If we have a periodic cell, they are set so that the cell is covered.
  //3. Otherwise we set it to include each ion plus  4 bohrs to capture the tails
  min_.Resize(ndim); max.Resize(ndim);

  
  for(int d=0; d< ndim; d++) { 
    min_(d)=0.0; max(d)=0.0;
  }
  
  int nions=sys->nIons();
  atominfo.Resize(nions, 4);
  Array1 <doublevar> ionpos(3);
  for(int i=0; i< nions; i++) {
    sys->getIonPos(i,ionpos);
    for(int d=0; d< 3; d++) {
      atominfo(i,d+1)=ionpos(d);
      if(ionpos(d) < min_(d)) min_(d) = ionpos(d);
      if(ionpos(d) > max(d)) max(d)=ionpos(d);
    }
    atominfo(i,0)=sys->getIonCharge(i);
  }
  
  for(int d=0; d< 3; d++) {
    min_(d)-=4.0;
    max(d)+=4.0;
  }
  
  Array2 <doublevar> latvec;
  Array1<doublevar> origin(3);
  if(sys->getBounds(latvec,origin)) { 
    //assume origin is zero for the moment
    min_=0; max=0;
    for(int d=0; d< ndim; d++) {
      for(int i=0; i< ndim; i++) {
        if(latvec(i,d)>0) max(d)+=latvec(i,d);
        if(latvec(i,d)<0) min_(d)+=latvec(i,d);
      }
    }
    //correct for the origin
    for(int d=0; d< ndim; d++) {
      max(d)+=origin(d);
      min_(d)+=origin(d);
    }
      
  }

  
  vector <string> mintxt;
  if(readsection(words, pos=0, mintxt, "MIN")) {
    if(mintxt.size()!=3) error("MIN must have exactly 3 elements");
    for(int d=0;d  < 3; d++) {
      min_(d)=atof(mintxt[d].c_str());
    }
  }

  vector<string> maxtxt;
  if(readsection(words, pos=0, maxtxt, "MAX")) {
    if(maxtxt.size()!=3) error("MAX must have exactly 3 elements");
    for(int d=0;d  < 3; d++) {
      max(d)=atof(maxtxt[d].c_str());
    }
  }

  //------end min/max stuff

  if(!readvalue(words, pos=0, resolution, "RESOLUTION"))
    resolution=.1;
  

  
   npoints.Resize(3);
   npoints=1;
  for(int d=0; d < ndim; d++) {
    npoints(d)=int((max(d)-min_(d))/resolution)+1;
  }

  bin.Resize(npoints(0), npoints(1), npoints(2));
  bin=0;
  nsample=0;

  


  //try to recover the previous run's density
  ifstream is(outputfile.c_str());
  if(is) { 
    string dummy;
    //each processor gets an equal share, since we gather them when
    //we write the file
    is >> dummy >> dummy >> nsample;
    nsample /= (mpi_info.nprocs);
    is.ignore(180, '\n'); //finish first line
    is.ignore(180, '\n');
    doublevar dum;
    is >> dum;
    if(dum != nions) error("different number of ions in previous density file");
    for(int d=0; d< 3; d++) {
      is >> dum;
      if(fabs(dum-min_(d)) > 1e-4 )error("different min in density file");
    }
    for(int d=0; d< 3; d++) {
      is >> dum;
      if(fabs(dum-npoints(d))> 1e-4) 
        error("different number of points in density file");
      for(int i=0; i< d; i++) 
        is >> dum;
      is >> dum;
      if(fabs(dum-resolution) > 1e-4) 
        error("different resolution in density file");
      is.ignore(180, '\n');
    }
    for(int at=0; at< nions; at++) 
      is.ignore(180, '\n');
    for(int x=0; x < npoints(0); x++) {
      for(int y=0; y < npoints(1); y++) {
        for(int z=0; z< npoints(2); z++) {
          is >> bin(x,y,z);

          bin(x,y,z)*=nsample/(norm);
        }
      }
    }   
    
    is.close();
  }
 
}

//----------------------------------------------------------------------

void One_particle_density::accumulate(Sample_point * sample, double weight) { 

  //int nelectrons=sample->electronSize();
  Array1 <int> place(3);
  Array1 <doublevar> epos(3);
  for(int e=start_electron; e< end_electron; e++) {
    sample->getElectronPos(e,epos);
    int use=1;
    for(int d=0; d< 3; d++) {
      place(d)=int( (epos(d)-min_(d))/resolution+0.5);
      if(place(d) <0 || place(d) >= npoints(d))
        use=0;
    }
    if(use) { 
      nsample+=weight;
      bin(place(0), place(1), place(2))+=weight;
    }
  }
}

//----------------------------------------------------------------------

void One_particle_density::add_single(Array1 <doublevar> & pos, doublevar val, 
    doublevar weight) { 
  //int nelec=end_electron-start_electron;
  Array1 <int> place(3);
  int use=1;
  for(int d=0; d< 3; d++) {
    place(d)=int((pos(d)-min_(d))/resolution+0.5);
    if(place(d)<0 ||place(d) >=npoints(d)) use=0;
  }
  //cout << "use: " << use << " place " << place(0) << " " << place(1) << " " << place(2) << endl;
  if(use) { 
    nsample+=weight;
    bin(place(0),place(1),place(2))+=val*weight;
  }
}


//----------------------------------------------------------------------

void One_particle_density::write() { 
#ifdef USE_MPI
  Array3 <doublevar> bin_tmp(npoints(0), npoints(1), npoints(2));
  bin_tmp=0;

  //MPI defaults to at max 4MB send queue(at least in the P4 communicator),
  //so we should split the array up into chunks of about 2MB, which is about
  //260,000 doubles, to be safe.  This can be avoided by setting P4_GLOBMEMSIZE
  //to some large value, but it's really annoying to do that, and there's no
  //huge performance benefit.

  int n=npoints(0)*npoints(1)*npoints(2);
  double * ptr=bin.v;
  double * ptr2=bin_tmp.v;
  double * last=bin.v+n;
  int interval=260000;
  
  while(ptr<last) { 
    int region=min( int(last-ptr) , interval);
    MPI_Reduce(ptr, ptr2, region, MPI_DOUBLE, MPI_SUM,
	       0, MPI_Comm_grp);
    ptr+=region;
    ptr2+=region;
  }

  doublevar nsample_tmp=parallel_sum(nsample);
#else
  Array3 <doublevar> & bin_tmp(bin);
  doublevar nsample_tmp=nsample;
#endif
  
  if(mpi_info.node==0) {
    string temp_file=outputfile+".backup";
    ofstream os(temp_file.c_str());
    os << "QWalk: nsamples " << nsample_tmp << "\n";
    os << "Electron density" << endl;
    
    int nions=atominfo.GetDim(0);
    os << "  " << nions << "   " << min_(0) << "   "
       << min_(1) << "   " << min_(2) << endl;
    os << npoints(0) << "   " << resolution << "  0.0000   0.0000" << endl;
    os << npoints(1) << "   0.0000   " << resolution << "  0.0000" << endl;
    os << npoints(2) << "   0.0000    0.0000    " << resolution<< endl;
    
    for(int at=0; at < nions; at++) {
      os << "   " << atominfo(at,0) << "  0.0000  " << atominfo(at,1)
         << "   " << atominfo(at,2) << "   " << atominfo(at,3) << endl;
    }

    int counter=0;
    for(int x=0; x < npoints(0); x++) {
      for(int y=0; y < npoints(1); y++) {
        for(int z=0; z< npoints(2); z++) {
          os << norm*bin_tmp(x,y,z)/nsample_tmp << "   ";
          if((counter++)%6==5) os << endl;
        }
      }
    }
    os << endl;
    os.close();
    rename(temp_file.c_str(), outputfile.c_str());
  }
        
  
}

//######################################################################

void Local_potential_density::init(vector <string> & words, 
    System * sys_, string & runid) { 
  sys=sys_;
  string runid2b=runid+"2b";
  density_2b.init(words,sys,runid2b);
  string runid1b=runid+"1b";
  density_1b.init(words,sys,runid1b);
}

void Local_potential_density::accumulate(Sample_point *sample,doublevar weight) { 
  Array1 <doublevar> onebody;
  Array2 <doublevar> twobody;
  sys->separatedLocal(sample,onebody, twobody);
  int nelec=twobody.GetDim(0);
  Array1<doublevar> two_sum(twobody.GetDim(0));
  two_sum=0.0;
  for(int i=0; i < nelec; i++) {
    for(int j=i+1; j< nelec; j++) { 
      two_sum(i)+=twobody(i,j);
      two_sum(j)+=twobody(i,j);
    }
  }
  Array1 <doublevar> pos(3);
  //cout << "adding points\n";
  for(int i=0; i< nelec; i++) { 
    sample->getElectronPos(i,pos);
    //cout << pos(0) <<"  " << pos(1) << " " << pos(2) << "  : " << weight << " : " << two_sum(i) << endl;
    density_2b.add_single(pos,two_sum(i),weight);
    density_1b.add_single(pos,onebody(i), weight);
  }
}


void Local_potential_density::write() { 
  density_2b.write();
  density_1b.write();
}

//######################################################################
//--------------------------------------------------

Properties_manager::Properties_manager() {
  nwf=1;
  maxchildren=3;
  start_avg_step=0;
  current_block=0;
  autocorr_depth=0;
  max_autocorr_depth=10;
  log_file="";
  //maxhist=0;
}

//--------------------------------------------------

Properties_manager::~Properties_manager() {
}

//--------------------------------------------------

void Properties_manager::read(vector <string> & words, 
                              vector <string> & systxt, 
                              vector <string> & wftxt) {

  if(block_avg.GetDim(0) > 0) 
    error("Must call Properties_manager::read before Properties_manager::setSize()");

  unsigned int pos=0;

  readvalue(words, pos=0, max_autocorr_depth, "MAX_AUTOCORR");

}


//--------------------------------------------------

void Properties_manager::setSize(int nwf_, int nblocks, int nsteps, 
                                 int maxwalkers,
                                 System * sys, 
                                 Wavefunction_data * wfdata) {
  assert(nblocks >= 0 );
  assert(nsteps >= 0);
  assert(maxwalkers >= 0);
  assert(nwf_ > 0);
  nwf=nwf_;

  current_block=0;
  npoints_this_block=0;

  block_avg.Resize(nblocks);
  for(int i=0; i< nblocks; i++) { 
    block_avg(i).setSize(nwf, 0,0);
    block_avg(i).aux_size=1;
  }
  
  autocorr_depth=min(nsteps-1, max_autocorr_depth);
  //trace.Resize(nsteps, maxwalkers);
  weighted_sum.setSize(nwf);
}

//--------------------------------------------------



void update_avgvar(doublevar & avg, doublevar & var, int pt_number, doublevar nwpt) { 
  doublevar oldavg=avg;
  doublevar oldvar=var;
  avg=oldavg+(nwpt-oldavg)/(pt_number+1);
  var=oldvar+(nwpt*nwpt-oldvar)/(pt_number+1);
  //if(pt_number > 0) 
  //  var=(1.0-1.0/pt_number)*oldvar+(pt_number+1)*(avg-oldavg);
  //else
  //  var=0;
}


void Properties_manager::insertPoint(int step, 
                                     int walker, 
                                     const Properties_point & pt) {

  //int nwf=pt.avgrets.GetDim(0);
  int navg=pt.avgrets.GetDim(1);
  assert(pt.kinetic.GetDim(0)==nwf);
  assert(pt.potential.GetDim(0)==nwf);
  assert(pt.nonlocal.GetDim(0)==nwf);
  assert(pt.weight.GetDim(0)==nwf);

  if(npoints_this_block==0) { 
    weighted_sum.setSize(nwf);
    sample_avg.setSize(nwf);
    sample_var.setSize(nwf);
    sample_avg.weight=0;
    sample_var.weight=0;
    weighted_sum.avgrets=pt.avgrets;
    weighted_sum.weight=0;
    energy_avg.Resize(nwf);
    energy_var.Resize(nwf);
    energy_avg=0;
    energy_var=0;
    //sample_avg.avgrets=pt.avgrets;
    //sample_var.avgrets=pt.avgrets;
    for(int w=0; w < nwf; w++) { 
      for(int a=0; a< navg; a++) { 
        for(int j=0; j< pt.avgrets(w,a).vals.GetDim(0); j++) { 
          weighted_sum.avgrets(w,a).vals(j)=0;
          //sample_avg.avgrets(w,a).vals(j)=0;
          //sample_var.avgrets(w,a).vals(j)=0;
        }
      }
    }
  }


  weighted_sum.weighted_add(pt);
  for(int w=0; w< nwf; w++) { 
    update_avgvar(sample_avg.kinetic(w),sample_var.kinetic(w),npoints_this_block,pt.kinetic(w));
    update_avgvar(sample_avg.potential(w),sample_var.potential(w),npoints_this_block,pt.potential(w));
    update_avgvar(sample_avg.nonlocal(w),sample_var.nonlocal(w),npoints_this_block,pt.nonlocal(w));
    update_avgvar(sample_avg.weight(w),sample_var.weight(w),npoints_this_block,pt.weight(w));
    update_avgvar(energy_avg(w),energy_var(w),npoints_this_block,pt.energy(w));
  }
  
  npoints_this_block++;
  

}


//--------------------------------------------------


/*!
find the energy autocorrelation.  In general, we have a directed 
graph, so let's start at the bottom and follow parents up.
*/
void Properties_manager::autocorrelation(Array2 <doublevar> & autocorr,
                                         int depth) {
  error("autocorrelation currently unsupported");
  /*
  using namespace Properties_types;
  assert(depth >= 0);
  int nsteps=trace.GetDim(0);
  int nwalkers=trace.GetDim(1);
  int nwf=trace(0,0).kinetic.GetDim(0);
  autocorr.Resize(nwf, depth);
  autocorr=0;
  
  for(int d=1; d< depth+1; d++) {
    int npts=0;
    for(int step=d+start_avg_step; step < nsteps; step++) {
      for(int walker=0; walker < nwalkers; walker++) {
        if(trace(step, walker).count) {

          //Here we trace back up to the dth ancestor of the walker
          int parent=trace(step, walker).parent;
          for(int d1=1; d1 < d; d1++) {
            //if(trace(step-d1, parent).parent <0) {
            //  cout << "step " << step << " d1 " << d1 << "parent " << parent
            //       << endl;
            //}
            assert(trace(step-d1, parent).parent >= 0);
            
            parent=trace(step-d1, parent).parent;
          }
          
          if(trace(step-d, parent).count) {
            npts++;
            
            
            for(int w=0; w< nwf; w++) {
              doublevar en=trace(step, walker).energy(w);
              doublevar enp=trace(step-d, parent).energy(w);

              //doublevar wt_av=block_avg(current_block).weight(w);
              //doublevar wt=trace(step, walker).weight(w);
              //doublevar pwt=trace(step, walker).weight(w);
              doublevar av=block_avg(current_block).avg(total_energy,w);
              
              autocorr(w,d-1)+=(en-av)*(enp-av);
            }


          }
        }
      }
    }
    
    npts=parallel_sum(npts);
    for(int w=0; w< nwf; w++) {
      autocorr(w,d-1)=parallel_sum(autocorr(w,d-1));
      autocorr(w,d-1)/=npts;
    }
  }
  for(int w=0; w< nwf; w++) {
    for(int d=0; d< depth; d++) {
      autocorr(w,d)/= block_avg(current_block).var(total_energy,w);
    }
  }
*/
}


//--------------------------------------------------


//a and a2 are assumed not to be updated over all processors
doublevar condense_variance(doublevar a2,doublevar a,int totpts, int npoints_this) { 
  doublevar a_avg=parallel_sum(npoints_this*a/totpts);
  return (parallel_sum(npoints_this*a2/totpts)-a_avg*a_avg)*(totpts/(totpts-1));
}
void Properties_manager::endBlock() {
  using namespace Properties_types;

  
  //int nwalkers=trace.GetDim(1);
  
  //int nsteps=trace.GetDim(0);
  //For the generalized averaging, this can be broken if 
  //the calling program isn't consistent with the ordering and 
  //number of Average_returns.  I can't think of any reason why someone
  //would want to do that other than spite, though.
  int navg_gen=weighted_sum.avgrets.GetDim(1);
  block_avg(current_block).avgrets.Resize(nwf, navg_gen);
  assert(weighted_sum.avgrets.GetDim(0)==nwf);
  for(int w=0; w < nwf; w++) { 
    doublevar totweight=parallel_sum(weighted_sum.weight(w));
    int totpts=parallel_sum(npoints_this_block);
    block_avg(current_block).avg(kinetic,w)=parallel_sum(weighted_sum.kinetic(w))/totweight;
    block_avg(current_block).avg(potential,w)=parallel_sum(weighted_sum.potential(w))/totweight;
    block_avg(current_block).avg(nonlocal,w)=parallel_sum(weighted_sum.nonlocal(w))/totweight;
    block_avg(current_block).avg(total_energy,w)=parallel_sum(weighted_sum.energy(w))/totweight;
    block_avg(current_block).avg(weight,w)=totweight/totpts;
    block_avg(current_block).totweight=totpts;
    for(int i=0; i< navg_gen; i++) { 
      int nvals=weighted_sum.avgrets(w,i).vals.GetDim(0);

      block_avg(current_block).avgrets(w,i)=weighted_sum.avgrets(w,i);
      parallel_sum(weighted_sum.avgrets(w,i).vals,
                   block_avg(current_block).avgrets(w,i).vals);
      for(int j=0; j< nvals; j++) { 
      //  block_avg(current_block).avgrets(w,i).vals(j)=
      //           parallel_sum(weighted_sum.avgrets(w,i).vals(j))/totweight;
          block_avg(current_block).avgrets(w,i).vals(j)/=totweight;
      }
    }

    int nprocs=mpi_info.nprocs;
    block_avg(current_block).var(kinetic, w)=
      condense_variance(sample_var.kinetic(w),sample_avg.kinetic(w),totpts,npoints_this_block);
    block_avg(current_block).var(potential, w)=
      condense_variance(sample_var.potential(w),sample_avg.potential(w),totpts,npoints_this_block);
    block_avg(current_block).var(nonlocal, w)=
      condense_variance(sample_var.nonlocal(w),sample_avg.nonlocal(w),totpts,npoints_this_block);
    block_avg(current_block).var(weight, w)=
      condense_variance(sample_var.weight(w),sample_avg.weight(w),totpts,npoints_this_block);
    block_avg(current_block).var(total_energy, w)=
      condense_variance(energy_var(w),energy_avg(w),totpts,npoints_this_block);
  }
  block_avg(current_block).autocorr.Resize(nwf,0);
  block_avg(current_block).autocorr=0;

  if(mpi_info.node==0 && log_file != "") {
    ofstream logout(log_file.c_str(), ios::app);
    logout.precision(16);
    logout << "block { " << endl;
    string indent="   ";
    block_avg(current_block).storeToLog(indent, logout, log_label);
    logout << "} " << endl << endl;
    logout.close();
  }
  current_block++;
  npoints_this_block=0;
  final_avg.blockReduce(block_avg, 0, current_block,1);
}


//--------------------------------------------------

/*!
\todo
Use the recurrence relations from 
http://mathworld.wolfram.com/SampleVarianceComputation.html
to remove the extraneous trace in Properties_manager..
 */

void Properties_manager::endBlock_per_step() {
  error("endBlock_per_step() depreciated");
  /*
  using namespace Properties_types;

  
  int nwalkers=trace.GetDim(1);
  
  int nsteps=trace.GetDim(0);
  int totpts=0;
  //For the generalized averaging, this can be broken if 
  //the calling program isn't consistent with the ordering and 
  //number of Average_returns.  I can't think of any reason why someone
  //would want to do that other than spite, though.
  int navg_gen=trace(0,0).avgrets.GetDim(1);
  block_avg(current_block).avgrets.Resize(nwf,navg_gen);
  
  for(int w=0; w< nwf; w++) {
    doublevar totweight=0;
    doublevar avgkin=0;
    doublevar avgpot=0;
    doublevar avgnonloc=0;
    doublevar avgen=0;
    for(int i=0; i< navg_gen; i++) { 
      block_avg(current_block).avgrets(w,i).vals.Resize(trace(0,0).avgrets(w,i).vals.GetDim(0));      
      block_avg(current_block).avgrets(w,i).vals=0;
      block_avg(current_block).avgrets(w,i).type=trace(0,0).avgrets(w,i).type;
    }
    
    Array1 <doublevar> weight_per_step(nsteps-start_avg_step,0.0);
    Array1 <int> npts_per_step(nsteps-start_avg_step,0);
    int npts=0;
    for(int step=start_avg_step; step < nsteps; step++) {
      for(int walker=0; walker < nwalkers; walker++) {
        if(trace(step, walker).count) {
          npts++;
          doublevar wt=trace(step,walker).weight(w);
          
          weight_per_step(step-start_avg_step)+=wt;
          npts_per_step(step-start_avg_step)++;
          //cout <<wt<<weight_per_step(step-start_avg_step)<<npts_per_step(step-start_avg_step)<<endl;
          
          totweight+=wt;
        }
      }
    }
    //cout << " totpts " << endl;
    totpts=parallel_sum(npts);
    //cout << "donepsum" << endl;
    
    Array1 <int> totnpts_per_step(nsteps-start_avg_step,0);
    Array1 <doublevar> totweight_per_step(nsteps-start_avg_step,0.0);
    for(int step=start_avg_step; step < nsteps; step++) {
      totnpts_per_step(step-start_avg_step)=parallel_sum(npts_per_step(step-start_avg_step));
      totweight_per_step(step-start_avg_step)=parallel_sum(weight_per_step(step-start_avg_step));
      //cout <<"step"<<step<<" totnpts_per_step "<<totnpts_per_step(step-start_avg_step)<<" totweight_per_step " <<totweight_per_step(step-start_avg_step)<<endl;
    }
    doublevar weight_sum=parallel_sum(totweight)/totpts;
    
    Array1 < Array1 < Array1 <doublevar> > > avgrets_per_step(nsteps-start_avg_step);
    for(int step=start_avg_step; step < nsteps; step++) {
      if(w==0) { 
        avgrets_per_step(step-start_avg_step).Resize(navg_gen);
        for(int i=0; i< navg_gen; i++){
          avgrets_per_step(step-start_avg_step)(i).Resize( block_avg(current_block).avgrets(w,i).vals.GetDim(0));
          avgrets_per_step(step-start_avg_step)(i)=0.0;
        }
      }
      
      for(int walker=0; walker < nwalkers; walker++) {
        if(trace(step, walker).count) {
          doublevar wt=trace(step,walker).weight(w);
          avgkin+=trace(step,walker).kinetic(w)*wt;
          avgpot+=trace(step, walker).potential(w)*wt;
          avgnonloc+=trace(step, walker).nonlocal(w)*wt;
          avgen+=(trace(step, walker).kinetic(w)
                  +trace(step, walker).potential(w)
                  +trace(step, walker).nonlocal(w))*wt;
          if(w==0) { 
            for(int i=0; i< navg_gen; i++) { 
              for(int j=0; j< block_avg(current_block).avgrets(w,i).vals.GetDim(0); j++) { 
                //block_avg(current_block).avgrets(i).vals(j)+=wt*trace(step,walker).avgrets(i).vals(j);
                avgrets_per_step(step-start_avg_step)(i)(j)+=wt*trace(step,walker).avgrets(w,i).vals(j);
              }
            }
          }
        }
      }
    }
    
    weight_sum*=totpts;
    
    //block_avg(current_block).totweight=weight_sum;
    block_avg(current_block).totweight=0.0;
    for(int step=start_avg_step; step < nsteps; step++) {
      if(totnpts_per_step(step-start_avg_step))
	block_avg(current_block).totweight+=totweight_per_step(step-start_avg_step);
    }
    
    block_avg(current_block).avg(kinetic,w)=parallel_sum(avgkin)/weight_sum;
    block_avg(current_block).avg(potential,w)=parallel_sum(avgpot)/weight_sum;
    block_avg(current_block).avg(nonlocal,w)=parallel_sum(avgnonloc)/weight_sum;
    block_avg(current_block).avg(total_energy,w)=parallel_sum(avgen)/weight_sum;
    //block_avg(current_block).avg(weight,w)=weight_sum/totpts;
    block_avg(current_block).avg(weight,w)=0.0;
    for(int step=start_avg_step; step < nsteps; step++) {
      if(totnpts_per_step(step-start_avg_step))
	block_avg(current_block).avg(weight,w)+=totweight_per_step(step-start_avg_step)/totnpts_per_step(step-start_avg_step);
    }
    block_avg(current_block).avg(weight,w)/=(nsteps-start_avg_step);
    
    if(w==0) { 
      for(int i=0; i< navg_gen; i++) { 
        for(int j=0; j< block_avg(current_block).avgrets(w,i).vals.GetDim(0); j++) { 
          //block_avg(current_block).avgrets(i).vals(j)=
          //parallel_sum(block_avg(current_block).avgrets(i).vals(j))/weight_sum;
	  block_avg(current_block).avgrets(w,i).vals(j)=0.0;
	  for(int step=start_avg_step; step < nsteps; step++) {
	    if(totnpts_per_step(step-start_avg_step))
	      block_avg(current_block).avgrets(w,i).vals(j)+=parallel_sum(avgrets_per_step(step-start_avg_step)(i)(j))/totnpts_per_step(step-start_avg_step);
	  }
	  block_avg(current_block).avgrets(w,i).vals(j)/=(nsteps-start_avg_step);
        }
      }
    }
    //cout << mpi_info.node << ":npoints" << npts << endl;
    

    //Variance
    doublevar varkin=0, varpot=0, varnonloc=0, varen=0;
    doublevar varweight=0;


    for(int step=start_avg_step; step < nsteps; step++) {
      for(int walker=0; walker < nwalkers; walker++) {
        if(trace(step, walker).count) {
          doublevar kin=trace(step, walker).kinetic(w);
          doublevar pot=trace(step, walker).potential(w);
          doublevar nonloc=trace(step, walker).nonlocal(w);
          doublevar en=kin+pot+nonloc;

          doublevar wt=trace(step, walker).weight(w);

          varkin+=(kin-block_avg(current_block).avg(kinetic,w))
            *(kin-block_avg(current_block).avg(kinetic,w))*wt;
          varpot+=(pot-block_avg(current_block).avg(potential,w))
            *(pot-block_avg(current_block).avg(potential,w))*wt;
          varnonloc+=(nonloc-block_avg(current_block).avg(nonlocal,w))
           *(nonloc-block_avg(current_block).avg(nonlocal,w))*wt;
          varen+=(en-block_avg(current_block).avg(total_energy,w))
            *(en-block_avg(current_block).avg(total_energy,w))*wt;

          varweight+=(wt-block_avg(current_block).avg(weight,w))
            *(wt-block_avg(current_block).avg(weight,w));

        }
      }
    }
    block_avg(current_block).var(kinetic,w)= parallel_sum(varkin)/weight_sum;
    block_avg(current_block).var(potential,w)=parallel_sum(varpot)/weight_sum;
    block_avg(current_block).var(nonlocal, w)=parallel_sum(varnonloc)/weight_sum;
    block_avg(current_block).var(total_energy, w)=parallel_sum(varen)/weight_sum;
    block_avg(current_block).var(weight,w)=parallel_sum(varweight)/totpts;
  }
  
  
  autocorrelation(block_avg(current_block).autocorr,
                  autocorr_depth);
  
  //cout << "writing block " << endl;
  if(mpi_info.node==0 && log_file != "") {
    ofstream logout(log_file.c_str(), ios::app);
    logout.precision(16);
    logout << "block { " << endl;
    string indent="   ";
    block_avg(current_block).storeToLog(indent, logout, log_label);
    logout << "} " << endl << endl;
    logout.close();
  }
  current_block++;
  for(int step=0; step < nsteps; step++)
    for(int walker=0; walker < nwalkers; walker++) 
      trace(step, walker).reset();
  
  //cout << "average " << endl;
  final_avg.blockReduce(block_avg, 0, current_block,1);
  //cout << "done" << endl;
  */
}

//--------------------------------------------------


void Properties_manager::endBlockSHDMC() {
  error("endBlockSHDMC() depreciated");
  /*
  using namespace Properties_types;

  
  int nwalkers=trace.GetDim(1);
  
  int nsteps=trace.GetDim(0);
  int totpts=0;
  //For the generalized averaging, this can be broken if 
  //the calling program isn't consistent with the ordering and 
  //number of Average_returns.  I can't think of any reason why someone
  //would want to do that other than spite, though.
  int navg_gen=trace(0,0).avgrets.GetDim(1);
  block_avg(current_block).avgrets.Resize(nwf,navg_gen);
  
  for(int w=0; w< nwf; w++) {
    doublevar totweight=0;
    doublevar avgkin=0;
    doublevar avgpot=0;
    doublevar avgnonloc=0;
    doublevar avgen=0;
    if(w==0) { 
      for(int i=0; i< navg_gen; i++) { 
        block_avg(current_block).avgrets(w,i).vals.Resize(trace(0,0).avgrets(w,i).vals.GetDim(0));      
        block_avg(current_block).avgrets(w,i).vals=0;
        block_avg(current_block).avgrets(w,i).type=trace(0,0).avgrets(w,i).type;
      }
    }
    int npts=0;
    for(int step=start_avg_step; step < nsteps; step++) {
      for(int walker=0; walker < nwalkers; walker++) {
        if(trace(step, walker).count) {
          npts++;
          doublevar wt=trace(step,walker).weight(w);
          totweight+=wt;
        }
      }
    }
    //cout << " totpts " << endl;
    totpts=parallel_sum(npts);
    //cout << "donepsum" << endl;
    doublevar weight_sum=parallel_sum(totweight)/totpts;
    
    for(int step=start_avg_step; step < nsteps; step++) {
      for(int walker=0; walker < nwalkers; walker++) {
        if(trace(step, walker).count) {
          doublevar wt=trace(step,walker).weight(w);
          avgkin+=trace(step,walker).kinetic(w)*wt;
          avgpot+=trace(step, walker).potential(w)*wt;
          avgnonloc+=trace(step, walker).nonlocal(w)*wt;
          avgen+=(trace(step, walker).kinetic(w)
                  +trace(step, walker).potential(w)
                  +trace(step, walker).nonlocal(w))*wt;
          if(w==0) { 
	    //using the averaging with w-1 if avgrets(i).type=="linear_delta_der"
            for(int i=0; i< navg_gen; i++) { 
              for(int j=0; j< block_avg(current_block).avgrets(w,i).vals.GetDim(0); j++) { 
                if(j<block_avg(current_block).avgrets(w,i).vals.GetDim(0)-1 && block_avg(current_block).avgrets(w,i).type=="linear_delta_der"){
		    block_avg(current_block).avgrets(w,i).vals(j)+=(wt-1.0)*trace(step,walker).avgrets(w,i).vals(j);
		}
		else 
		  block_avg(current_block).avgrets(w,i).vals(j)+=wt*trace(step,walker).avgrets(w,i).vals(j);
              }
            }
          }
        }
      }
    }
    
    weight_sum*=totpts;
    
    block_avg(current_block).totweight=weight_sum;
    
    block_avg(current_block).avg(kinetic,w)=parallel_sum(avgkin)/weight_sum;
    block_avg(current_block).avg(potential,w)=parallel_sum(avgpot)/weight_sum;
    block_avg(current_block).avg(nonlocal,w)=parallel_sum(avgnonloc)/weight_sum;
    block_avg(current_block).avg(total_energy,w)=parallel_sum(avgen)/weight_sum;
    block_avg(current_block).avg(weight,w)=weight_sum/totpts;
    
    if(w==0) { 
      for(int i=0; i< navg_gen; i++) { 
        for(int j=0; j< block_avg(current_block).avgrets(w,i).vals.GetDim(0); j++) { 
          block_avg(current_block).avgrets(w,i).vals(j)=
          parallel_sum(block_avg(current_block).avgrets(w,i).vals(j))/weight_sum;
        }
      }
    }
    //cout << mpi_info.node << ":npoints" << npts << endl;
    

    //Variance
    doublevar varkin=0, varpot=0, varnonloc=0, varen=0;
    doublevar varweight=0;


    for(int step=start_avg_step; step < nsteps; step++) {
      for(int walker=0; walker < nwalkers; walker++) {
        if(trace(step, walker).count) {
          doublevar kin=trace(step, walker).kinetic(w);
          doublevar pot=trace(step, walker).potential(w);
          doublevar nonloc=trace(step, walker).nonlocal(w);
          doublevar en=kin+pot+nonloc;

          doublevar wt=trace(step, walker).weight(w);

          varkin+=(kin-block_avg(current_block).avg(kinetic,w))
            *(kin-block_avg(current_block).avg(kinetic,w))*wt;
          varpot+=(pot-block_avg(current_block).avg(potential,w))
            *(pot-block_avg(current_block).avg(potential,w))*wt;
          varnonloc+=(nonloc-block_avg(current_block).avg(nonlocal,w))
           *(nonloc-block_avg(current_block).avg(nonlocal,w))*wt;
          varen+=(en-block_avg(current_block).avg(total_energy,w))
            *(en-block_avg(current_block).avg(total_energy,w))*wt;

          varweight+=(wt-block_avg(current_block).avg(weight,w))
            *(wt-block_avg(current_block).avg(weight,w));
        }
      }
    }
    block_avg(current_block).var(kinetic,w)= parallel_sum(varkin)/weight_sum;
    block_avg(current_block).var(potential,w)=parallel_sum(varpot)/weight_sum;
    block_avg(current_block).var(nonlocal, w)=parallel_sum(varnonloc)/weight_sum;
    block_avg(current_block).var(total_energy, w)=parallel_sum(varen)/weight_sum;
    block_avg(current_block).var(weight,w)=parallel_sum(varweight)/totpts;
    
  }
  
  autocorrelation(block_avg(current_block).autocorr,
                  autocorr_depth);
  
  if(mpi_info.node==0 && log_file != "") {
    ofstream logout(log_file.c_str(), ios::app);
    logout.precision(16);
    logout << "block { " << endl;
    string indent="   ";
    block_avg(current_block).storeToLog(indent, logout, log_label);
    logout << "} " << endl << endl;
    logout.close();
  }
  current_block++;
  for(int step=0; step < nsteps; step++)
    for(int walker=0; walker < nwalkers; walker++) 
      trace(step, walker).reset();
  
  final_avg.blockReduce(block_avg, 0, current_block,1);
  */
}


//--------------------------------------------------

void Properties_manager::initializeLog(Array2 <Average_generator*> & avg_gen) {
  if(mpi_info.node==0 && log_file != ""){ 
    ofstream os(log_file.c_str(), ios::app);
    os << "init { \n";
    string indent="   ";
    os << indent << "label " << log_label << endl;
    os << indent << "nsys " << avg_gen.GetDim(0) <<  " navg " << avg_gen.GetDim(1) << endl;
    for(int i=0; i< avg_gen.GetDim(0); i++) { 
      for(int j=0; j< avg_gen.GetDim(1); j++) { 
        os << indent << "average_generator { \n";
        avg_gen(i,j)->write_init(indent,os);
        os << indent << "}\n";
      }
    }
    os << "}\n";
    os.close();
  }
}



//--------------------------------------------------

void Properties_manager::initializeLog(Array1 <Average_generator*> & avg_gen) {
  if(mpi_info.node==0 && log_file != ""){ 
    ofstream os(log_file.c_str(), ios::app);
    os << "init { \n";
    string indent="   ";
    os << indent << "label " << log_label << endl;
    for(int i=0; i< avg_gen.GetDim(0); i++) { 
      os << indent << "average_generator { \n";
      avg_gen(i)->write_init(indent,os);
      os << indent << "}\n";
    }
    os << "}\n";
    os.close();
  }
}


//--------------------------------------------------

  
  

void Properties_manager::printBlockSummary(ostream & os) {
  assert(current_block > 0);
  block_avg(current_block-1).printBlockSummary(os);
}

//--------------------------------------------------


void Properties_manager::printSummary(ostream & os, Array1 <Average_generator*> & avg_gen) {
  final_avg.showSummary(os,avg_gen);
}

//--------------------------------------------------
