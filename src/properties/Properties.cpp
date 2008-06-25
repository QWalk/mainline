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
  else if(caseless_eq(words[0], "SK"))
    denspt=new Structure_factor;
  else if(caseless_eq(words[0], "LM"))
    denspt=new Local_moments;
  else if(caseless_eq(words[0],"POL"))
    denspt=new Polarization;
  else if(caseless_eq(words[0],"GR"))
    denspt=new Electron_correlation;
  else
    error("Didn't understand density keyword",words[0]);
  denspt->init(words, sys, runid);
}

void allocate(vector<string> & words, System * sys, string & runid,
	 Nonlocal_density_accumulator *& nldenspt) { 
  if(words.size() < 1) error("empty non-local density section");
  if(caseless_eq(words[0],"OBDM"))
    nldenspt=new OBDM;
  else
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
    
    //////// ######### ion .dx write #############
    /* This needs to utilize the format for .dx code which specifies
      * the attributes of each point based on a non-uniform grid.
      */
    ///////////// ############## end .dx write ########    
    
  }
        
  
}

//----------------------------------------------------------------------
//######################################################################

void Electron_correlation::init(vector <string> & words, System * sys, 
                                string & runid) { 
  outputfile=runid+".gr";
  resolution=.1;
  doublevar range=10.0;
  npoints=int(range/resolution);
  gr_up.Resize(npoints); gr_down.Resize(npoints);
  gr_unlike.Resize(npoints);
  gr_up=0.0; gr_down=0.0; gr_unlike=0.0;
  nsample_up=nsample_down=nsample_unlike=0;
  nup=sys->nelectrons(0);
  
  ifstream is(outputfile.c_str());
  if(is) { 
    is.ignore(180,'\n'); //we currently can't set npoints & resolution
    string dum;
    is >> dum;
    is >> nsample_up >> nsample_down >> nsample_unlike;
    nsample_up/=mpi_info.nprocs;
    nsample_down/=mpi_info.nprocs;
    nsample_unlike/=mpi_info.nprocs;
    is.ignore(180,'\n'); //finish line
    is.ignore(180,'\n'); //info line
    double ddum;
    for(int i=0; i< npoints; i++) { 
      is >> ddum >> gr_up(i) >> gr_down(i) >> gr_unlike(i);
      gr_up(i)*=nsample_up;
      gr_down(i)*=nsample_down;
      gr_unlike(i)*=nsample_unlike;
    }
    
    is.close();
  }
}

//----------------------------------------------------------------------


void Electron_correlation::accumulate(Sample_point * sample, double weight) { 
  //cout << "acc " << endl;
  int nelectrons=sample->electronSize();
  int nions=sample->ionSize();
  sample->updateEEDist();
  sample->updateEIDist();
  Array1 <doublevar> dist(5);
  Array1 <doublevar> dist_ion(5);
  for(int e=0; e< nup; e++) {

    for(int e2=e+1; e2 < nup; e2++) { 
      sample->getEEDist(e,e2,dist);
      int place=int(dist(0)/resolution);
      if(place < npoints) { 
        gr_up(place)+=weight/(4*pi*resolution*dist(1));
        nsample_up+=weight;
      }
    }
    //cout << "unlike " << endl;
    //unlike spins
    for(int e2=nup; e2 < nelectrons; e2++) {
      
      sample->getEEDist(e,e2,dist);
      int place=int(dist(0)/resolution);
      if(place < npoints) { 
        gr_unlike(place)+=weight/(4*pi*resolution*dist(1));
        nsample_unlike+=weight;
      }
    }    
  }
  //cout << "down " << endl;
  
  //like down spins
  for(int e=nup; e < nelectrons; e++) { 
    for(int e2=e+1; e2< nelectrons; e2++) { 
      sample->getEEDist(e,e2,dist);
      int place=int(dist(0)/resolution);
      if(place < npoints) { 
        gr_down(place) += weight/(4*pi*resolution*dist(1));
        nsample_down+=weight;
      }
    }
  }
  //cout << "done " << endl;
}

//----------------------------------------------------------------------


void Electron_correlation::write() { 
  doublevar n_up_tmp=parallel_sum(nsample_up);
  doublevar n_down_tmp=parallel_sum(nsample_down);
  doublevar n_unlike_tmp=parallel_sum(nsample_unlike);
  
#ifdef USE_MPI
  Array1 <doublevar> down_tmp(npoints),up_tmp(npoints),unlike_tmp(npoints);
  MPI_Reduce(gr_down.v,down_tmp.v,npoints,MPI_DOUBLE, MPI_SUM, 0,MPI_Comm_grp);
  MPI_Reduce(gr_up.v,up_tmp.v,npoints,MPI_DOUBLE, MPI_SUM, 0,MPI_Comm_grp);
  MPI_Reduce(gr_unlike.v,unlike_tmp.v,npoints,MPI_DOUBLE, MPI_SUM, 0,MPI_Comm_grp);

#else

  Array1 <doublevar> & down_tmp(gr_down);
  Array1 <doublevar> & up_tmp(gr_up);
  Array1 <doublevar> & unlike_tmp(gr_unlike);
#endif
  
  if(mpi_info.node==0) { 
    ofstream os(outputfile.c_str());
    os << "#npoints " << npoints << " resolution " << resolution << endl;
    os << "#weights " << n_up_tmp << "  " << n_down_tmp << "  " << n_unlike_tmp << endl;
    os << "#r    g(r)(spin up)  g(r)(spin down)   g(r)(unlike spins) \n";
    for(int i=0; i< npoints; i++) { 
      os << (i+.5)*resolution << "   " << up_tmp(i)/n_up_tmp
      << "  " << down_tmp(i)/n_down_tmp << "  " << unlike_tmp(i)/n_unlike_tmp
      << endl;
    }
  }
  
}

//######################################################################
void Structure_factor::init(vector <string> & words, System * sys,
			 string & runid) { 

  single_write(cout, "Structure factor calculation!");
  outputfile=runid+".sk";

  nsample=0;
  wnsample=0.0;
  unsigned int pos=0;
  if(!readvalue(words, pos=0, np_side, "NGRID")) 
    np_side=5;
  
  npoints=np_side*np_side*np_side; 
  grid.Resize(npoints);
  grid=0.0;

  vector <string> gvec_sec;
  if(readsection(words, pos=0, gvec_sec, "GVEC")) {
    int count=0;
    gvec.Resize(3,3);
    for(int i=0; i< 3; i++) { 
      for(int j=0; j< 3; j++) { 
        gvec(i,j)=atof(gvec_sec[count++].c_str());
      }
    }
  }
  else { 
    if(!sys->getRecipLattice(gvec)) 
      error("You don't have a periodic cell and you haven't specified GVEC for S(k)");
  }
  
  kpts.Resize(npoints, 3);
  kpts=0;
  int c=0;
  for(int ix=0; ix < np_side; ix++) {
    for(int iy=0; iy < np_side; iy++) {
      for(int iz=0; iz < np_side; iz++) {
        for(int i=0; i< 3; i++)
          kpts(c,i)+=2*pi*(gvec(0,i)*ix+gvec(1,i)*iy+gvec(2,i)*iz);
        c++;
      }
    }
  }
}

//----------------------------------------------------------------------

void Structure_factor::accumulate(Sample_point * sample, doublevar weight) { 
  nelectrons=sample->electronSize();

  Array1 <doublevar> pos(3);
  Array1 <doublevar> k(3);

  for(int p=0; p < npoints; p++) { 
    doublevar sum_cos=0, sum_sin=0;
    
    for(int e=0; e< nelectrons; e++) {
      sample->getElectronPos(e,pos);
      doublevar dot=0;
      for(int d=0; d< 3; d++) dot+=pos(d)*kpts(p,d);
      sum_cos+=cos(dot);
      sum_sin+=sin(dot);
    }
    grid(p)+=weight*(sum_cos*sum_cos+sum_sin*sum_sin);
    
  }
  nsample+=1;
  wnsample+=weight;
}

//----------------------------------------------------------------------

void Structure_factor::write() { 
  int nsample_tmp=parallel_sum(nsample);
  doublevar wnsample_tmp=parallel_sum(wnsample);
#ifdef USE_MPI
  Array1 <doublevar> grid_tmp(npoints);
  grid_tmp=0;
  MPI_Reduce(grid.v, grid_tmp.v, npoints, MPI_DOUBLE,MPI_SUM,
	     0,MPI_Comm_grp);
#else
  Array1 <doublevar> & grid_tmp(grid);
  
#endif  

  if(mpi_info.node==0) { 
    ofstream out(outputfile.c_str());
    out << "#nsamples  " << nsample_tmp << " ("
	<< wnsample_tmp << " weighted)" << endl;
    out << "#k  s(k)  kx  ky  kz" << endl;
    
    for(int i=0; i< npoints; i++) {
      doublevar sk=grid_tmp(i)/wnsample_tmp;
      doublevar r=0;
      for(int d=0; d< 3; d++) r+=kpts(i,d)*kpts(i,d);
      r=sqrt(r);
      out << r << "   " << sk/nelectrons << "   "
          << kpts(i,0) << "   " << kpts(i,1) << "   " << kpts(i,2) << endl;
    }
  }
}
//######################################################################


//######################################################################
void Local_moments::init(vector <string> & words, System * sys,
			 string & runid) { 

  single_write(cout, "Local moments will be calculated.\n");
  outputfile=runid+".lm";
  outputfilelog=runid+".lm.log";

  nup=sys->nelectrons(0);
  nelectrons=nup+sys->nelectrons(1);
  natoms=sys->nIons();
  
  moment.Resize(natoms+1);  // last atom is interstitial
  moment=0;
  charge.Resize(natoms+1);
  charge=0;
  rMT.Resize(natoms+1);
  rMT=-1;
  sys->getAtomicLabels(atomnames);
  if ( words.size() == 1 ) {
    single_write(cout,"No muffin-tin radii given.\n");
  } else if ( 2*(words.size()/2) == words.size() ) {
    single_write(cout,"Malformed muffin-tin radii data.\n");
  } else {
    // let's try to decode the input
    int npair=(words.size()-1)/2;
    for (int pair=0; pair < npair; pair++ ) {
      for (unsigned int i=0; i < atomnames.size(); i++) {
	if ( caseless_eq(words[2*pair+1],atomnames[i]) ) {
	  sscanf(words[2*pair+2].c_str(),"%lf",&(rMT(i)));
	}
      }
    }
  }
  single_write(cout,"Muffin-tin radii:\n");
  for (unsigned int i=0; i < atomnames.size(); i++) {
    if ( rMT(i) > 0 ) {
      single_write(cout,atomnames[i],"  ",rMT(i),"\n");
    } else {
      rMT(i)=2;
      single_write(cout,atomnames[i],"  ",rMT(i)," (default used)\n");
    }
  }
  atomnames.push_back("interst.");

  nblock=0;
  nsample=0;
  wnsample=0.0;

}

//----------------------------------------------------------------------

void Local_moments::accumulate(Sample_point * sample, doublevar weight) { 

  Array1 <doublevar> dist(5);

  sample->updateEIDist();
  
  for(int e=0; e<nup; e++) {
    bool interstitial=true;
    for(int i=0; i<natoms; i++) {
      sample->getEIDist(e,i,dist);
      // zero element of dist is norm, second is norm^2
      if ( dist(0) < rMT(i) ) {
	moment(i)+=weight;
	charge(i)+=weight;
	interstitial=false;
      }
    }
    if ( interstitial ) {
      moment(natoms)+=weight;
      charge(natoms)+=weight;
    }
  }
  
  for(int e=nup; e<nelectrons; e++) {
    bool interstitial=true;
    for(int i=0; i<natoms; i++) {
      sample->getEIDist(e,i,dist);
      if ( dist(0) < rMT(i) ) {
	moment(i)-=weight;
	charge(i)+=weight;
	interstitial=false;
      }
    }
    if ( interstitial ) {
      moment(natoms)-=weight;
      charge(natoms)+=weight;
    }
  }

  nsample+=1;
  wnsample+=weight;

}

//----------------------------------------------------------------------

void Local_moments::write() { 

  nblock+=1;
  LM_storage this_block;
  this_block.allocate(natoms+1);            // last atom is interstitial
  this_block.nsample=parallel_sum(nsample);
  this_block.wnsample=parallel_sum(wnsample);

#ifdef USE_MPI
  MPI_Reduce(moment.v, this_block.moment.v, natoms+1, MPI_DOUBLE,MPI_SUM,
	     0,MPI_Comm_grp);
  MPI_Reduce(charge.v, this_block.charge.v, natoms+1, MPI_DOUBLE,MPI_SUM,
	     0,MPI_Comm_grp);
#else
  this_block.moment=moment;
  this_block.charge=charge;
#endif

  if(mpi_info.node==0) { 
    
    blocks.push_back(this_block);
    
    int nsample_tmp=0;
    doublevar wnsample_tmp=0.0;
    for ( list <LM_storage>::const_iterator block = blocks.begin();
	    block != blocks.end(); block++ ) {
      nsample_tmp+=block->nsample;
      wnsample_tmp+=block->wnsample;
    }

    ofstream out(outputfile.c_str());
    out << "No. samples: " << nsample_tmp << " ("
	<< wnsample_tmp << " weighted)" << endl;
    out << "No. blocks:  " << nblock << endl << endl;
    out << "    atom(rMT)           moment                  charge" << endl;
    out.flags(ios::fixed);
    
    for(int i=0; i<natoms+1; i++) {
      // average
      doublevar m_av=0, ch_av=0;
      for ( list <LM_storage>::const_iterator block = blocks.begin();
	    block != blocks.end(); block++ ) {
	m_av += block->moment(i)/block->wnsample;
	ch_av+= block->charge(i)/block->wnsample;
      }
      m_av =m_av/nblock;
      ch_av=ch_av/nblock;
      // variance
      doublevar m_var=0, ch_var=0;
      if ( nblock > 1 ) {
	for ( list <LM_storage>::const_iterator block = blocks.begin();
	      block != blocks.end(); block++ ) {
	  m_var += (block->moment(i)/block->wnsample-m_av)
	    *(block->moment(i)/block->wnsample-m_av);
	  ch_var+= (block->charge(i)/block->wnsample-ch_av)
	    *(block->charge(i)/block->wnsample-ch_av);
	}
	m_var =m_var/(nblock-1)/nblock;
	ch_var=ch_var/(nblock-1)/nblock;
      }
      // write out (at last!)
      out.precision(2); 
      out << setw(8) << atomnames[i];
      if ( i < natoms ) {
	out << "(" << setw(4) << rMT(i) << ")" ;
      } else {
	out << "      " ;
      }
      out.precision(5);
      out << "   " << setw(8) << m_av << " +/- " << setw(8) << sqrt(m_var)
	  << "   " << setw(8) << ch_av << " +/- " << setw(8) << sqrt(ch_var)
	  << endl;
    }
    
    out.close();
    
    // the out, as it is now, is fine but not perfect because it gives
    // wrong errorbars if blocks are too short and reblocking is
    // needed.
    ofstream outlog(outputfilelog.c_str(),ios_base::app);
    outlog << endl;
    outlog << "block {" << endl;
    outlog << "   totweight " << this_block.wnsample << endl;
    
    for(int i=0; i<natoms+1; i++)
      outlog << "   moment_" << atomnames[i] << "_" << i 
	     << " { " << this_block.moment(i)/this_block.wnsample << " } " << endl;
    for(int i=0; i<natoms+1; i++)
      outlog << "   charge_" << atomnames[i] << "_" << i
	     << " { " << this_block.charge(i)/this_block.wnsample << " } " << endl;
    
    outlog << "}" << endl;
    outlog.close();
    
  } // if(mpi_info.node==0)
  
  // reset per-block quantities
  nsample=0;
  wnsample=0.0;
  moment=0;
  charge=0;
  
}
//######################################################################

void Polarization::init(vector <string> & words, System * sys, string & runid) { 
  outputfile=runid+".pol";
  if(!sys->getRecipLattice(gvec))
    error("There must be a reciprocal lattice for polarization to work");
  ncvg=6;
  single_pol_cvg.Resize(ncvg,3);
  manye_pol_cvg.Resize(ncvg,3);
  single_pol_cvg=dcomplex(0.0,0.0);
  manye_pol_cvg=dcomplex(0.0,0.0);
  nsample=0;
}

//------------------------------------------------------------------------


void Polarization::accumulate(Sample_point * sample, doublevar weight) { 
  int nelectrons=sample->electronSize();
  Array1 <doublevar> sum(3,0.0);
  Array1 <doublevar> pos(3);
  for(int e=0; e< nelectrons; e++) { 
    sample->getElectronPos(e,pos);
    for(int i=0; i< 3; i++) { 
      doublevar tmp=0;
      for(int d=0; d< 3; d++) { 
        tmp+=gvec(i,d)*pos(d);
      }
      sum(i)+=tmp;
      for(int n=0; n< ncvg; n++) { 
        int nf=n+1;
        single_pol_cvg(n,i)+=weight*dcomplex(cos(nf*2*pi*tmp),sin(nf*2*pi*tmp))/doublevar(nelectrons);
      }
    }
  }
  
  for(int n=0; n < ncvg; n++) { 
    for(int i=0; i< 3; i++) { 
      int nf=n+1;
      manye_pol_cvg(n,i)+=weight*dcomplex(cos(nf*2*pi*sum(i)),sin(nf*2*pi*sum(i)));
    }
  }
  nsample+=weight;
                                       
}

//------------------------------------------------------------------------

void Polarization::write() { 
  for(int n=0; n < ncvg; n++) { 
    for(int i=0; i< 3; i++) { 
      single_pol_cvg(n,i)=parallel_sum(single_pol_cvg(n,i));
      manye_pol_cvg(n,i)=parallel_sum(manye_pol_cvg(n,i));
    }
  }
  nsample=parallel_sum(nsample);
  
  if(mpi_info.node==0) { 
    
    ofstream out(outputfile.c_str(), ios::app);
    out.precision(15);
    
    for(int n=0; n< ncvg; n++) {
      for(int i=0; i< 3; i++) { 
        out << single_pol_cvg(n,i).real()/nsample  << "  " << single_pol_cvg(n,i).imag()/nsample << "   ";
      }
      for(int i=0; i< 3; i++) { 
        out << manye_pol_cvg(n,i).real()/nsample << "  " << manye_pol_cvg(n,i).imag()/nsample << "  ";
      }
    }
    out << endl;
    out.close();

  }
  nsample=0;
  single_pol_cvg=dcomplex(0.0,0.0);
  manye_pol_cvg=dcomplex(0.0,0.0);
}


//######################################################################
void OBDM::init(vector <string> & words, System * sys,
		string & runid) { 

  single_write(cout, "One-body density matrix will be calculated.\n");
  //  outputfile=runid+".obdm";
  outputfilelog=runid+".obdm.log";

  nelectrons=sys->nelectrons(0)+sys->nelectrons(1);

  unsigned int pos=0;
  if(!readvalue(words, pos=0, np_side, "NGRID")) 
    np_side=5;
  
  // I assume spherical symmetry for now, shifts will be along the first
  // lattice vector
  npoints=np_side;          
  dmtrx.Resize(npoints);
  dmtrx=0.0;
  
  int ndim=3;
  latVec.Resize(ndim,ndim);
  if(!sys->getBounds(latVec)) 
    error("As of now, OBDM is only implemented for periodic systems.");
  //dR=sqrt(latVec(0,0)*latVec(0,0)
  //	  +latVec(0,1)*latVec(0,1)+latVec(0,2)*latVec(0,2))/npoints;

  Array2 <doublevar> crossProduct(ndim,ndim);
  crossProduct(0,0)=(latVec(1,1)*latVec(2,2)-latVec(1,2)*latVec(2,1));
  crossProduct(0,1)=(latVec(1,2)*latVec(2,0)-latVec(1,0)*latVec(2,2));
  crossProduct(0,2)=(latVec(1,0)*latVec(2,1)-latVec(1,1)*latVec(2,0));

  crossProduct(1,0)=(latVec(2,1)*latVec(0,2)-latVec(2,2)*latVec(0,1));
  crossProduct(1,1)=(latVec(2,2)*latVec(0,0)-latVec(2,0)*latVec(0,2));
  crossProduct(1,2)=(latVec(2,0)*latVec(0,1)-latVec(2,1)*latVec(0,0));

  crossProduct(2,0)=(latVec(0,1)*latVec(1,2)-latVec(0,2)*latVec(1,1));
  crossProduct(2,1)=(latVec(0,2)*latVec(1,0)-latVec(0,0)*latVec(1,2));
  crossProduct(2,2)=(latVec(0,0)*latVec(1,1)-latVec(0,1)*latVec(1,0));
  
  doublevar smallestheight=1e99;
  for(int i=0; i< ndim; i++) {
    doublevar tempheight=0;
    doublevar length=0;
    for(int j=0; j< ndim; j++) {
      tempheight+=crossProduct(i,j)*latVec(i,j);
      length+=crossProduct(i,j)*crossProduct(i,j);
    }
    tempheight=fabs(tempheight)/sqrt(length);
    if(tempheight < smallestheight ) smallestheight=tempheight;
  }

  // maximum distance for OBDM evaluation
  pos=0;
  doublevar cutoff;
  if(!readvalue(words, pos=0, cutoff, "CUTOFF"))
    dR=smallestheight/npoints/2;
  else
    dR=cutoff/npoints;
  
  // angles and weights for Gaussian quadrature (spherical averaging)
  // (the same material as in system/gesqua.cpp)

  pos=0;
  if(!readvalue(words, pos=0, np_aver, "AIP")) np_aver=1;
  wt.Resize(np_aver);
  ptc.Resize(np_aver,3);             //!< cartesian coordinates of int. points
  Array2 <doublevar> pts(np_aver,2); //!< spherical coordinates of int. points

  switch (np_aver) {
  case 1:
    // no spherical averaging
    np_aver=1;
    wt=1;
    pts(0,0)=pi/2;
    pts(0,1)=0;
    break;
  case 12:
    // icosahedron 12 point rule
    wt=1.0/12;
    pts(0,0)=0;
    pts(0,1)=0;
    pts(1,0)=pi;
    pts(1,1)=0;
    for (int i=0; i<5; i++) {
      pts(i+2,0)=atan(2.0);
      pts(i+2,1)=2*i*pi/5;
      pts(i+7,0)=pi-atan(2.0);
      pts(i+7,1)=(2*i+1)*pi/5;
    }
    break;
  case 32:
    // icosahedron 32 point rule
    { // braces needed due to declaration of th1,th2 below
      for (int i=0; i<12; i++) wt(i)=5.0/168;
      for (int i=12; i<32; i++) wt(i)=27.0/840;
      pts(0,0)=0;
      pts(0,1)=0;
      pts(1,0)=pi;
      pts(1,1)=0;
      double th1=acos((2+sqrt(5.0))/sqrt(15+6*sqrt(5.0)));
      double th2=acos(1.0/sqrt(15+6*sqrt(5.0)));
      for (int i=0; i<5; i++) {
	pts(i+2,0)=atan(2.0);
	pts(i+2,1)=2*i*pi/5;
	pts(i+7,0)=pi-atan(2.0);
	pts(i+7,1)=(2*i+1)*pi/5;
	pts(i+12,0)=th1;
	pts(i+12,1)=(2*i+1)*pi/5;
	pts(i+17,0)=th2;
	pts(i+17,1)=(2*i+1)*pi/5;
	pts(i+22,0)=pi-th1;
	pts(i+22,1)=2*i*pi/5;
	pts(i+27,0)=pi-th2;
	pts(i+27,1)=2*i*pi/5;
      }
    }
    break;
  default:
    error("Invalid number of integration points (AIP) for spherical averaging in OBDM, allowed are 1, 12 and 32.");
  }

  for (int i=0; i<np_aver; i++) {
    ptc(i,0)=sin(pts(i,0))*cos(pts(i,1));
    ptc(i,1)=sin(pts(i,0))*sin(pts(i,1));
    ptc(i,2)=cos(pts(i,0));
  }

  nsample=0;
  wnsample=0.0;
  
  if(mpi_info.node==0) { 
    ofstream outlog(outputfilelog.c_str(),ios_base::app);
    outlog << endl;
    outlog << "#----------------------------------------------------------" 
	   << endl;
    outlog << "#One body density matrix, spherical average with "
	   << np_aver << " points." << endl;
    outlog << "#----------------------------------------------------------" 
	   << endl;
  }

}

//----------------------------------------------------------------------

void OBDM::accumulate(Sample_point * sample, doublevar weight, 
		      Wavefunction_data * wfdata, Wavefunction * wf) { 

  //cout << "OBDM accumulate" << endl;

  int nwf=wf->nfunc();
  Wf_return wfval_new(nwf,2);   // this structure I just don't understand
  Wf_return wfval_old(nwf,2);
  Storage_container wfStore;
  Array1 <doublevar> oldpos(3), newpos(3), transl(3);

  Array1 <doublevar> x(3), y(3), z(3);
  Array2 <doublevar> pt(np_aver,3);
  generate_random_rotation(x,y,z);
  for (int i=0; i<np_aver; i++) {
    pt(i,0)=ptc(i,0)*x(0)+ptc(i,1)*y(0)+ptc(i,2)*z(0);
    pt(i,1)=ptc(i,0)*x(1)+ptc(i,1)*y(1)+ptc(i,2)*z(1);
    pt(i,2)=ptc(i,0)*x(2)+ptc(i,1)*y(2)+ptc(i,2)*z(2);
  }

  wf->updateVal(wfdata, sample);
  wfStore.initialize(sample, wf);
  wfStore.saveUpdate(sample,wf,0);
  wf->getVal(wfdata, 0, wfval_old);
  sample->getElectronPos(0,oldpos);

  for( int i=0; i<npoints; i++) {
    for ( int k=0; k<np_aver; k++) {
      // shift makes sense only up to half the lattice vector, we are hitting
      // periodicity for larger distances
      //for (int j=0; j<3; j++) newpos(j)=oldpos(j)+i*latVec(0,j)/2/npoints;
      //sample->setElectronPos(0,newpos);
      transl(0)=(i+1)*dR*pt(k,0);
      transl(1)=(i+1)*dR*pt(k,1);
      transl(2)=(i+1)*dR*pt(k,2);
      sample->translateElectron(0,transl);
      // TODO: check k-points /= Gamma
      wf->updateVal(wfdata, sample);
      wf->getVal(wfdata, 0, wfval_new);
      // how to handle the possibility of wf_old=0? Can it happen?
      dmtrx(i)+=wt(k)*wfval_new.sign(0)*wfval_old.sign(0)
	*exp(wfval_new.amp(0,0)-wfval_old.amp(0,0))*weight;
      sample->setElectronPos(0,oldpos);
    }
  }
  
  wfStore.restoreUpdate(sample,wf,0);
  nsample+=1;
  wnsample+=weight;
  
}

//----------------------------------------------------------------------

void OBDM::write(string & log_label) { 

  doublevar wnsample_tmp=parallel_sum(wnsample);
#ifdef USE_MPI
  Array1 <doublevar> dmtrx_tmp(npoints);
  dmtrx_tmp=0;
  MPI_Reduce(dmtrx.v, dmtrx_tmp.v, npoints, MPI_DOUBLE,MPI_SUM,
	     0,MPI_Comm_grp);
#else
  Array1 <doublevar> & dmtrx_tmp(dmtrx);
#endif

  if(mpi_info.node==0) { 
    
    // the out, as it is now, is fine but not perfect because it gives
    // wrong errorbars if blocks are too short and reblocking is
    // needed.
    ofstream outlog(outputfilelog.c_str(),ios_base::app);
    outlog << endl;
    outlog << "block {" << endl;
    outlog << "   label " << log_label << endl;
    outlog << "   totweight " << wnsample_tmp << endl;
    
    for(int i=0; i<npoints; i++)
      outlog 
	//<< "   r(" << i << ")="
	<< dR*(i+1)
	<< " { " << dmtrx_tmp(i)/wnsample_tmp
	<< " } " << endl;

    outlog << "}" << endl;
    outlog.close();

  } // if(mpi_info.node==0)

  // reset per-block quantities
  nsample=0;
  wnsample=0.0;
  dmtrx=0.0;
  
}
//######################################################################


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
