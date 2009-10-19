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

//######################################################################

void allocate(vector<string> & words, System * sys, string & runid,
	 Local_density_accumulator *& denspt) { 
  if(words.size() < 1) error("empty density section");
  if(caseless_eq(words[0], "DENSITY"))
    denspt=new One_particle_density;
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
  norm=50;
  outputfile=runid+".cube";
  nup=sys->nelectrons(0);
  unsigned int pos=0;
  
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
  if(sys->getBounds(latvec)) { 
    //assume origin is zero for the moment
    min_=0; max=0;
    for(int d=0; d< ndim; d++) {
      for(int i=0; i< ndim; i++) {
        if(latvec(i,d)>0) max(d)+=latvec(i,d);
        if(latvec(i,d)<0) min_(d)+=latvec(i,d);
      }
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
  

  start_electron=0;
  end_electron=sys->nelectrons(0)+sys->nelectrons(1);
  if(haskeyword(words, pos=0, "UP"))
    end_electron=sys->nelectrons(0);
  
  if(haskeyword(words,pos=0, "DOWN"))
    start_electron=sys->nelectrons(0);

  if(end_electron==start_electron)
    error("In one-particle density, UP and DOWN are incompatible");
  
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
    ofstream os(outputfile.c_str());
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
    /*  This isn't really used by anyone, is it?  It takes up a lot of space to write the density twice..
    string outputfile2 = outputfile+".dx";
    os.clear();
    os.open(outputfile2.c_str());
    
    /////////////######### DX FILE WRITING ###########
    
    /////##### DX FILE - Probabilities ##########
    //
    //    os << "  " << nions << "   " << min_(0) << "   "
    //       << min_(1) << "   " << min_(2) << endl;
    
    
    os << "object 1 class gridpositions counts " << npoints(0) << " "
      << npoints(1) << " " << npoints(2) << "\n";
    os << "origin " << min_(0) << " "
      << min_(1) << " " << min_(2) << " " << endl;
    
    os << "delta " << resolution << "   0.0   0.0" << endl;
    os << "delta 0.0   " << resolution << "   0.0" << endl;
    os << "delta 0.0   0.0   " << resolution << endl;
    
    os << endl;
    os << "object 2 class gridconnections counts " << npoints(0) << " "
      << npoints(1) << " " << npoints(2) << "\n";
    os << "attribute \"element type\" string \"cubes\" " << endl;
    os << "attribute \"ref\" string \"positions\" " << endl;
    
    
    os << endl;
    os << "object 3 class array type float rank 0 items " << (npoints(0) * npoints(1) * npoints(2)) <<
      " data follows" << endl;
    os << endl;
    
    
    counter=0;
    for(int x=0; x < npoints(0); x++) {
      for(int y=0; y < npoints(1); y++) {
        for(int z=0; z< npoints(2); z++) {
          os << norm*bin_tmp(x,y,z)/nsample_tmp << "   ";
          if((counter++)%6==5) os << endl;
        }
      }
    }
    os << endl;
    
    os << "#attribute \"dep\" string \"positions\" " << endl;
    os << "object \"regular positions regular connections\" class field" << endl;
    os << "component \"positions\" value 1" << endl;
    os << "component \"connections\" value 2" << endl;
    os << "component \"data\" value 3" << endl;
    os << "end" << endl;
    
    
    os.close();
    */
    //////// ######### ion .dx write #############
    /* This needs to utilize the format for .dx code which specifies
      * the attributes of each point based on a non-uniform grid.
      */
    ///////////// ############## end .dx write ########    
    
  }
        
  
}

//######################################################################

//--------------------------------------------------

Properties_manager::Properties_manager() {
  nwf=1;
  maxchildren=3;
  start_avg_step=0;
  current_block=0;
  autocorr_depth=0;
  num_aux_converge=1;
  naux=0;
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
                                 Wavefunction_data * wfdata, int naux_, 
				 int n_aux_cvg) {
  assert(nblocks >= 0 );
  assert(nsteps >= 0);
  assert(maxwalkers >= 0);
  assert(nwf_ > 0);
  nwf=nwf_;

  current_block=0;
  naux=naux_;
  num_aux_converge=n_aux_cvg;
  
  block_avg.Resize(nblocks);
  for(int i=0; i< nblocks; i++) { 
    block_avg(i).setSize(nwf, naux, n_aux_cvg);
    block_avg(i).aux_size=1;
  }
  
  autocorr_depth=min(nsteps-1, max_autocorr_depth);
  trace.Resize(nsteps, maxwalkers);
}

//--------------------------------------------------




void Properties_manager::insertPoint(int step, 
                                     int walker, 
                                     const Properties_point & pt) {
  assert(walker < trace.GetDim(1));
  assert(step < trace.GetDim(0));
  //cout << "insertpt " << pt.z_pol(0) << endl;
  trace(step, walker)=pt;
  if(step >= 1 && pt.parent < 0) {
    error("Problem in insertPoint; parent not set");
  }
}


//--------------------------------------------------


/*!
find the energy autocorrelation.  In general, we have a directed 
graph, so let's start at the bottom and follow parents up.
*/
void Properties_manager::autocorrelation(Array2 <doublevar> & autocorr,
					 Array2 <doublevar> & aux_autocorr,
                                         int depth) {

  using namespace Properties_types;
  assert(depth >= 0);
  int nsteps=trace.GetDim(0);
  int nwalkers=trace.GetDim(1);
  int nwf=trace(0,0).kinetic.GetDim(0);
  autocorr.Resize(nwf, depth);
  autocorr=0;
  aux_autocorr.Resize(naux, depth);
  aux_autocorr=0;

  
  
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

	    for(int i=0; i< naux; i++) {
	      doublevar avdiff=block_avg(current_block).aux_energy(i,0)
		-block_avg(current_block).avg(total_energy,0);
	      //doublevar av_auxwt=block_avg(current_block).aux_weight(i,0);
	      //doublevar av_wt=block_avg(current_block).weight(0);
	      doublevar diff=trace(step, walker).aux_energy(i,0)
		-trace(step, walker).energy(0);
	      doublevar pdiff=trace(step-d, parent).aux_energy(i,0)
		-trace(step-d, parent).energy(0);
	      aux_autocorr(i,d-1)+=(diff-avdiff)*(pdiff-avdiff);

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
    for(int i=0; i< naux; i++) {
      aux_autocorr(i,d-1)=parallel_sum(aux_autocorr(i,d-1));
      aux_autocorr(i,d-1)/=(npts*block_avg(current_block).aux_diffvar(i,0));
    }      
  }
  for(int w=0; w< nwf; w++) {
    for(int d=0; d< depth; d++) {
      autocorr(w,d)/= block_avg(current_block).var(total_energy,w);
    }
  }

}


//--------------------------------------------------

/*!
\todo
Use the recurrence relations from 
http://mathworld.wolfram.com/SampleVarianceComputation.html
to remove the extraneous trace in Properties_manager..
 */

void Properties_manager::endBlock() {
  using namespace Properties_types;

  
  int nwalkers=trace.GetDim(1);
  int n_cvg=num_aux_converge;
  
  int nsteps=trace.GetDim(0);
  int totpts=0;
  //For the generalized averaging, this can be broken if 
  //the calling program isn't consistent with the ordering and 
  //number of Average_returns.  I can't think of any reason why someone
  //would want to do that other than spite, though.
  int navg_gen=trace(0,0).avgrets.GetDim(0);
  block_avg(current_block).avgrets.Resize(navg_gen);
  
  for(int w=0; w< nwf; w++) {
    doublevar totweight=0;
    doublevar avgkin=0;
    doublevar avgpot=0;
    doublevar avgnonloc=0;
    doublevar avgen=0;
    if(w==0) { 
      for(int i=0; i< navg_gen; i++) { 
        block_avg(current_block).avgrets(i).vals.Resize(trace(0,0).avgrets(i).vals.GetDim(0));      
        block_avg(current_block).avgrets(i).vals=0;
        block_avg(current_block).avgrets(i).type=trace(0,0).avgrets(i).type;
      }
    }
    Array2 <doublevar> auxen(naux, n_cvg,0.0);
    Array2 <doublevar> auxwt(naux,n_cvg,0.0);
    Array2 <doublevar> auxdiff(naux,n_cvg, 0.0);
    int npts=0;
    for(int step=start_avg_step; step < nsteps; step++) {
      for(int walker=0; walker < nwalkers; walker++) {
        if(trace(step, walker).count) {
          npts++;
          doublevar wt=trace(step,walker).weight(w);
          totweight+=wt;
          for(int i=0; i< naux; i++) {
            for(int c=0; c< n_cvg; c++) 
              auxwt(i,c)+=trace(step, walker).aux_weight(i,c);
          }
        }
      }
    }
    //cout << " totpts " << endl;
    totpts=parallel_sum(npts);
    //cout << "donepsum" << endl;
    Array2 <doublevar> aux_wt_sum(naux, n_cvg);
    for(int i=0; i < naux; i++) 
      for(int c=0; c< n_cvg; c++) 
        aux_wt_sum(i,c)=parallel_sum(auxwt(i,c))/totpts;
    doublevar weight_sum=parallel_sum(totweight)/totpts;
    doublevar sum_inv=1.0/weight_sum;
    Array1 <dcomplex> z_pol(3, dcomplex(0.0, 0.0));
    Array2 <dcomplex> aux_z_pol(naux, 3, dcomplex(0.0, 0.0));
    
    
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
            for(int i=0; i< navg_gen; i++) { 
              for(int j=0; j< block_avg(current_block).avgrets(i).vals.GetDim(0); j++) { 
                block_avg(current_block).avgrets(i).vals(j)+=wt*trace(step,walker).avgrets(i).vals(j);
              }
            }
          }
          for(int d=0; d< 3; d++)
            z_pol(d)+=trace(step, walker).z_pol(d)*wt;
          for(int i=0; i< naux; i++) {
            //the convention(in RMC and VMC) is that the last
            //converge point is the 'pure' estimator
            int last=n_cvg-1;
            doublevar aux_wt=trace(step,walker).aux_weight(i,last)
              /(aux_wt_sum(i,last)*totpts);
            for(int d=0; d< 3; d++) 
              aux_z_pol(i,d)+=trace(step,walker).aux_z_pol(i,d)*aux_wt;
            
            for(int c=0; c< n_cvg; c++) {
              auxen(i,c)+=trace(step, walker).aux_energy(i,c)
              *trace(step, walker).aux_weight(i,c)/(aux_wt_sum(i,c)*totpts);
              auxdiff(i,c)+=trace(step,walker).aux_energy(i,c)
                *trace(step,walker).aux_weight(i,c)/aux_wt_sum(i,c)
                -trace(step,walker).energy(w)*wt*sum_inv;
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
    for(int d=0;d < 3; d++) {
      block_avg(current_block).z_pol(d)=parallel_sum(z_pol(d))/weight_sum;
    }
    for(int i=0; i< naux; i++) {
      for(int d=0;d < 3; d++) 
        block_avg(current_block).aux_z_pol(i,d)=parallel_sum(aux_z_pol(i,d));
      for(int c=0; c< n_cvg; c++) {
        aux_wt_sum(i,c)*=totpts;
        block_avg(current_block).aux_energy(i,c)=parallel_sum(auxen(i,c));
        block_avg(current_block).aux_weight(i,c)=aux_wt_sum(i,c)/totpts;
        block_avg(current_block).aux_diff(i,c)=parallel_sum(auxdiff(i,c))/totpts;
      }
    }
    
    if(w==0) { 
      for(int i=0; i< navg_gen; i++) { 
        for(int j=0; j< block_avg(current_block).avgrets(i).vals.GetDim(0); j++) { 
          block_avg(current_block).avgrets(i).vals(j)=
          parallel_sum(block_avg(current_block).avgrets(i).vals(j))/weight_sum;
        }
      }
    }
    //cout << mpi_info.node << ":npoints" << npts << endl;
    

    //Variance
    doublevar varkin=0, varpot=0, varnonloc=0, varen=0;
    doublevar varweight=0;

    Array2 <doublevar> auxvaren(naux,n_cvg, 0.0);
    Array2 <doublevar> auxvarwt(naux,n_cvg, 0.0);

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


          for(int i=0; i< naux; i++) {
            for(int c=0; c< n_cvg; c++) {
            auxvarwt(i,c)+=(trace(step, walker).aux_weight(i,c)
                          -block_avg(current_block).aux_weight(i,c))
                         *(trace(step, walker).aux_weight(i,c)
                           -block_avg(current_block).aux_weight(i,c));
            auxvaren(i,c)+=(trace(step, walker).aux_energy(i,c)
                          -block_avg(current_block).aux_energy(i,c))
                        *(trace(step, walker).aux_energy(i,c)
                          -block_avg(current_block).aux_energy(i,c))
                         *trace(step, walker).aux_weight(i,c);
            }
          }

        }
      }
    }
    block_avg(current_block).var(kinetic,w)= parallel_sum(varkin)/weight_sum;
    block_avg(current_block).var(potential,w)=parallel_sum(varpot)/weight_sum;
    block_avg(current_block).var(nonlocal, w)=parallel_sum(varnonloc)/weight_sum;
    block_avg(current_block).var(total_energy, w)=parallel_sum(varen)/weight_sum;
    block_avg(current_block).var(weight,w)=parallel_sum(varweight)/totpts;
    
    for(int i=0; i< naux; i++) {
      for(int c=0; c< n_cvg; c++) {
        block_avg(current_block).aux_energyvar(i,c)
        =parallel_sum(auxvaren(i,c))/aux_wt_sum(i,c);
        block_avg(current_block).aux_weightvar(i,c)
          =parallel_sum(auxvarwt(i,c))/totpts;
      }
    }
    
    Array2 <doublevar> auxvardiff(naux,n_cvg, 0.0);
    for(int step=start_avg_step; step < nsteps; step++) {
      for(int walker=0; walker < nwalkers; walker++) {
        if(trace(step, walker).count) {
          for(int i=0; i< naux; i++) {
            for(int c=0; c< n_cvg; c++) {
              doublevar avdiff=block_avg(current_block).aux_energy(i,c)
              -block_avg(current_block).avg(total_energy,0);
              doublevar diff=trace(step, walker).aux_energy(i,c)
                -trace(step, walker).energy(0);
              auxvardiff(i,c)+=(diff-avdiff)*(diff-avdiff);
            }
          }
        }
      }
    }
    
    for(int i=0; i< naux; i++) {
      for(int c=0; c< n_cvg; c++) {
        block_avg(current_block).aux_diffvar(i,c)
        =parallel_sum(auxvardiff(i,c))/totpts;
      }
    }
    
    
  }
  

  
  //cout << "autocorrelation " << endl;
  if(naux==1 && n_cvg==2) { 
    
    doublevar weight_corr=0;
    for(int step=start_avg_step; step < nsteps; step++) {
      for(int walker=0; walker < nwalkers; walker++) {
        if(trace(step, walker).count) {
          
          weight_corr+=(trace(step,walker).aux_weight(0,0)-
                        block_avg(current_block).aux_weight(0,0))
          *(trace(step,walker).aux_weight(0,0)-
            block_avg(current_block).aux_weight(0,0));
        }
      }
    }
    
    block_avg(current_block).aux_weight_correlation=
      parallel_sum(weight_corr)/
      (totpts*sqrt(block_avg(current_block).aux_weightvar(0,0)
                   *block_avg(current_block).aux_weightvar(0,1)));
  }
  
  autocorrelation(block_avg(current_block).autocorr,
                  block_avg(current_block).aux_autocorr,
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
}


//--------------------------------------------------

/*!
\todo
Use the recurrence relations from 
http://mathworld.wolfram.com/SampleVarianceComputation.html
to remove the extraneous trace in Properties_manager..
 */

void Properties_manager::endBlock_per_step() {
  using namespace Properties_types;

  
  int nwalkers=trace.GetDim(1);
  int n_cvg=num_aux_converge;
  
  int nsteps=trace.GetDim(0);
  int totpts=0;
  //For the generalized averaging, this can be broken if 
  //the calling program isn't consistent with the ordering and 
  //number of Average_returns.  I can't think of any reason why someone
  //would want to do that other than spite, though.
  int navg_gen=trace(0,0).avgrets.GetDim(0);
  block_avg(current_block).avgrets.Resize(navg_gen);
  
  for(int w=0; w< nwf; w++) {
    doublevar totweight=0;
    doublevar avgkin=0;
    doublevar avgpot=0;
    doublevar avgnonloc=0;
    doublevar avgen=0;
    if(w==0) { 
      for(int i=0; i< navg_gen; i++) { 
        block_avg(current_block).avgrets(i).vals.Resize(trace(0,0).avgrets(i).vals.GetDim(0));      
        block_avg(current_block).avgrets(i).vals=0;
        block_avg(current_block).avgrets(i).type=trace(0,0).avgrets(i).type;
      }
    }
    Array2 <doublevar> auxen(naux, n_cvg,0.0);
    Array2 <doublevar> auxwt(naux,n_cvg,0.0);
    Array2 <doublevar> auxdiff(naux,n_cvg, 0.0);
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
          for(int i=0; i< naux; i++) {
            for(int c=0; c< n_cvg; c++) 
              auxwt(i,c)+=trace(step, walker).aux_weight(i,c);
          }
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
    Array2 <doublevar> aux_wt_sum(naux, n_cvg);
    for(int i=0; i < naux; i++) 
      for(int c=0; c< n_cvg; c++) 
        aux_wt_sum(i,c)=parallel_sum(auxwt(i,c))/totpts;
    doublevar weight_sum=parallel_sum(totweight)/totpts;
    doublevar sum_inv=1.0/weight_sum;
    Array1 <dcomplex> z_pol(3, dcomplex(0.0, 0.0));
    Array2 <dcomplex> aux_z_pol(naux, 3, dcomplex(0.0, 0.0));
    
    

    Array1 < Array1 < Array1 <doublevar> > > avgrets_per_step(nsteps-start_avg_step);
    for(int step=start_avg_step; step < nsteps; step++) {
      if(w==0) { 
	avgrets_per_step(step-start_avg_step).Resize(navg_gen);
	for(int i=0; i< navg_gen; i++){
	  avgrets_per_step(step-start_avg_step)(i).Resize( block_avg(current_block).avgrets(i).vals.GetDim(0));
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
              for(int j=0; j< block_avg(current_block).avgrets(i).vals.GetDim(0); j++) { 
                //block_avg(current_block).avgrets(i).vals(j)+=wt*trace(step,walker).avgrets(i).vals(j);
		avgrets_per_step(step-start_avg_step)(i)(j)+=wt*trace(step,walker).avgrets(i).vals(j);
              }
            }
          }
          for(int d=0; d< 3; d++)
            z_pol(d)+=trace(step, walker).z_pol(d)*wt;
          for(int i=0; i< naux; i++) {
            //the convention(in RMC and VMC) is that the last
            //converge point is the 'pure' estimator
            int last=n_cvg-1;
            doublevar aux_wt=trace(step,walker).aux_weight(i,last)
              /(aux_wt_sum(i,last)*totpts);
            for(int d=0; d< 3; d++) 
              aux_z_pol(i,d)+=trace(step,walker).aux_z_pol(i,d)*aux_wt;
            
            for(int c=0; c< n_cvg; c++) {
              auxen(i,c)+=trace(step, walker).aux_energy(i,c)
              *trace(step, walker).aux_weight(i,c)/(aux_wt_sum(i,c)*totpts);
              auxdiff(i,c)+=trace(step,walker).aux_energy(i,c)
                *trace(step,walker).aux_weight(i,c)/aux_wt_sum(i,c)
                -trace(step,walker).energy(w)*wt*sum_inv;
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

    
    for(int d=0;d < 3; d++) {
      block_avg(current_block).z_pol(d)=parallel_sum(z_pol(d))/weight_sum;
    }
    for(int i=0; i< naux; i++) {
      for(int d=0;d < 3; d++) 
        block_avg(current_block).aux_z_pol(i,d)=parallel_sum(aux_z_pol(i,d));
      for(int c=0; c< n_cvg; c++) {
        aux_wt_sum(i,c)*=totpts;
        block_avg(current_block).aux_energy(i,c)=parallel_sum(auxen(i,c));
        block_avg(current_block).aux_weight(i,c)=aux_wt_sum(i,c)/totpts;
        block_avg(current_block).aux_diff(i,c)=parallel_sum(auxdiff(i,c))/totpts;
      }
    }
    
    if(w==0) { 
      for(int i=0; i< navg_gen; i++) { 
        for(int j=0; j< block_avg(current_block).avgrets(i).vals.GetDim(0); j++) { 
          //block_avg(current_block).avgrets(i).vals(j)=
          //parallel_sum(block_avg(current_block).avgrets(i).vals(j))/weight_sum;
	  block_avg(current_block).avgrets(i).vals(j)=0.0;
	  for(int step=start_avg_step; step < nsteps; step++) {
	    if(totnpts_per_step(step-start_avg_step))
	      block_avg(current_block).avgrets(i).vals(j)+=parallel_sum(avgrets_per_step(step-start_avg_step)(i)(j))/totnpts_per_step(step-start_avg_step);
	  }
	  block_avg(current_block).avgrets(i).vals(j)/=(nsteps-start_avg_step);
        }
      }
    }
    //cout << mpi_info.node << ":npoints" << npts << endl;
    

    //Variance
    doublevar varkin=0, varpot=0, varnonloc=0, varen=0;
    doublevar varweight=0;

    Array2 <doublevar> auxvaren(naux,n_cvg, 0.0);
    Array2 <doublevar> auxvarwt(naux,n_cvg, 0.0);

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


          for(int i=0; i< naux; i++) {
            for(int c=0; c< n_cvg; c++) {
            auxvarwt(i,c)+=(trace(step, walker).aux_weight(i,c)
                          -block_avg(current_block).aux_weight(i,c))
                         *(trace(step, walker).aux_weight(i,c)
                           -block_avg(current_block).aux_weight(i,c));
            auxvaren(i,c)+=(trace(step, walker).aux_energy(i,c)
                          -block_avg(current_block).aux_energy(i,c))
                        *(trace(step, walker).aux_energy(i,c)
                          -block_avg(current_block).aux_energy(i,c))
                         *trace(step, walker).aux_weight(i,c);
            }
          }

        }
      }
    }
    block_avg(current_block).var(kinetic,w)= parallel_sum(varkin)/weight_sum;
    block_avg(current_block).var(potential,w)=parallel_sum(varpot)/weight_sum;
    block_avg(current_block).var(nonlocal, w)=parallel_sum(varnonloc)/weight_sum;
    block_avg(current_block).var(total_energy, w)=parallel_sum(varen)/weight_sum;
    block_avg(current_block).var(weight,w)=parallel_sum(varweight)/totpts;
    
    for(int i=0; i< naux; i++) {
      for(int c=0; c< n_cvg; c++) {
        block_avg(current_block).aux_energyvar(i,c)
        =parallel_sum(auxvaren(i,c))/aux_wt_sum(i,c);
        block_avg(current_block).aux_weightvar(i,c)
          =parallel_sum(auxvarwt(i,c))/totpts;
      }
    }
    
    Array2 <doublevar> auxvardiff(naux,n_cvg, 0.0);
    for(int step=start_avg_step; step < nsteps; step++) {
      for(int walker=0; walker < nwalkers; walker++) {
        if(trace(step, walker).count) {
          for(int i=0; i< naux; i++) {
            for(int c=0; c< n_cvg; c++) {
              doublevar avdiff=block_avg(current_block).aux_energy(i,c)
              -block_avg(current_block).avg(total_energy,0);
              doublevar diff=trace(step, walker).aux_energy(i,c)
                -trace(step, walker).energy(0);
              auxvardiff(i,c)+=(diff-avdiff)*(diff-avdiff);
            }
          }
        }
      }
    }
    
    for(int i=0; i< naux; i++) {
      for(int c=0; c< n_cvg; c++) {
        block_avg(current_block).aux_diffvar(i,c)
        =parallel_sum(auxvardiff(i,c))/totpts;
      }
    }
    
    
  }
  

  
  //cout << "autocorrelation " << endl;
  if(naux==1 && n_cvg==2) { 
    
    doublevar weight_corr=0;
    for(int step=start_avg_step; step < nsteps; step++) {
      for(int walker=0; walker < nwalkers; walker++) {
        if(trace(step, walker).count) {
          
          weight_corr+=(trace(step,walker).aux_weight(0,0)-
                        block_avg(current_block).aux_weight(0,0))
          *(trace(step,walker).aux_weight(0,0)-
            block_avg(current_block).aux_weight(0,0));
        }
      }
    }
    
    block_avg(current_block).aux_weight_correlation=
      parallel_sum(weight_corr)/
      (totpts*sqrt(block_avg(current_block).aux_weightvar(0,0)
                   *block_avg(current_block).aux_weightvar(0,1)));
  }
  
  autocorrelation(block_avg(current_block).autocorr,
                  block_avg(current_block).aux_autocorr,
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
}

//--------------------------------------------------

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
