// Properties.h

class LM_storage {
 public:
  Array1 <doublevar> moment;
  Array1 <doublevar> charge;
  int nsample;
  doublevar wnsample;

  void allocate(int natoms) {
    moment.Resize(natoms);
    moment=0;
    charge.Resize(natoms);
    charge=0;
  }
};

class Local_moments:public Local_density_accumulator { 
 public:
  virtual void init(vector <string> & , System *, string & runid);
  virtual void accumulate(Sample_point *, doublevar weight);
  virtual void write();
 protected:
  string outputfile;
  string outputfilelog;
  int nsample;         // samples in a block
  doublevar wnsample;  // weighted sample count in a block
  int nblock;
  int nelectrons;
  int nup;             // number of up electrons
  int natoms;

  Array1 <doublevar> rMT;      // muffin-tin radius
  Array1 <doublevar> moment;   // hit-up - hit-down in MT sphere
  Array1 <doublevar> charge;   // hit count into MT sphere
  vector <string> atomnames;
  list <LM_storage> blocks;
  
};

// ***************************************************************************

// Properties.cpp
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
