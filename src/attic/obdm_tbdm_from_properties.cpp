// Properties.h

class OBDM:public Nonlocal_density_accumulator { 
 public:
    virtual void init(vector <string> & , System *, string & runid);
    virtual void accumulate(Sample_point *, doublevar weight,
			    Wavefunction_data *, Wavefunction *);
    virtual void write(string & log_label);
 protected:
    string outputfilelog;
    int nsample;
    doublevar wnsample;
    int nelectrons;    
    int npoints;
    int np_side;  //!< number of points on the side

    doublevar dR; //!< increment in the parameter of ODBM

    int np_aver;  //!< number of points of quadrature for spherical average
    Array1 <doublevar> wt;  //!< angles and weights for Gaussian quadrature
    Array2 <doublevar> ptc; //!< cartesian coordinates of integration points

    Array2 <doublevar> latVec;
    Array1 <doublevar> dmtrx;

};


/*!
Two-body density matrix (TBDM) implemented as Nonlocal_density_accumulator. More precisely, 
 it is only a particular element of TBDM that allows
 detection of up-down, BCS-like, condensate. Calculation of TBDM needs two-particle storage in Wavefunction,
 which is supported in Slater, Jastrow and BCS so far. Implemented by Matous Ringel.
*/
class TBDM:public Nonlocal_density_accumulator { 
 public:
    virtual void init(vector <string> & , System *, string & runid);
    virtual void accumulate(Sample_point *, doublevar weight,
			    Wavefunction_data *, Wavefunction *);
    virtual void write(string & log_label);
 protected:
    string outputfilelog;
    int nsample;
    doublevar wnsample;
    int nelectrons;    
    int npoints;
    int np_side;  //!< number of points on the side

    doublevar dR; //!< increment in the parameter of TDBM

    int np_aver;  //!< number of points of quadrature for spherical average
    Array1 <doublevar> wt;  //!< angles and weights for Gaussian quadrature
    Array2 <doublevar> ptc; //!< cartesian coordinates of integration points

    Array2 <doublevar> latVec;
    Array1 <doublevar> dmtrx;

};

// ****************************************************************************

// Properties.cpp
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
  Array1 <doublevar> x(np_aver), y(np_aver), z(np_aver);

  switch (np_aver) {
  case 1:
    // no spherical averaging
    wt=1;
    ptc(0,0)=1.0;
    ptc(0,1)=0.0;
    ptc(0,2)=0.0;
    break;
  default:
    gesqua(np_aver,x,y,z,wt);
    for (int i=0; i<np_aver; i++) {
      ptc(i,0)=x(i);
      ptc(i,1)=y(i);
      ptc(i,2)=z(i);
    }
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

// Added by Matous
//######################################################################
void TBDM::init(vector <string> & words, System * sys,
		string & runid) { 

  single_write(cout, "Two-body density matrix will be calculated.\n");
  //  outputfile=runid+".tbdm";
  outputfilelog=runid+".tbdm.log";

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
    error("As of now, TBDM is only implemented for periodic systems.");
 
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

  // maximum distance for TBDM evaluation
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
  Array1 <doublevar> x(np_aver), y(np_aver), z(np_aver);

  switch (np_aver) {
  case 1:
    // no spherical averaging
    wt=1;
    ptc(0,0)=1.0;
    ptc(0,1)=0.0;
    ptc(0,2)=0.0;
    break;
  default:
    gesqua(np_aver,x,y,z,wt);
    for (int i=0; i<np_aver; i++) {
      ptc(i,0)=x(i);
      ptc(i,1)=y(i);
      ptc(i,2)=z(i);
    }
  }

  nsample=0;
  wnsample=0.0;
  
  if(mpi_info.node==0) { 
    ofstream outlog(outputfilelog.c_str(),ios_base::app);
    outlog << endl;
    outlog << "#----------------------------------------------------------" 
	   << endl;
    outlog << "#Two-body density matrix, spherical average with "
	   << np_aver << " points." << endl;
    outlog << "#----------------------------------------------------------" 
	   << endl;
  }

}

//----------------------------------------------------------------------

void TBDM::accumulate(Sample_point * sample, doublevar weight, 
		      Wavefunction_data * wfdata, Wavefunction * wf) { 

  //cout << "TBDM accumulate" << endl;

  int nwf=wf->nfunc();
  Wf_return wfval_new(nwf,2);   // this structure I just don't understand
  Wf_return wfval_old(nwf,2);
  
  Wavefunction_storage * wfStore=NULL;
  //Config_save_point  sampleStore;
  Sample_storage * sampleStorage0=NULL, * sampleStorage1 = NULL;
 
  Array1 <doublevar> oldpos(3), newpos(3), transl(3);
  Array1 <doublevar> oldpos2(3); 

  Array1 <doublevar> x(3), y(3), z(3);
  Array2 <doublevar> pt(np_aver,3);
  generate_random_rotation(x,y,z);
  for (int i=0; i<np_aver; i++) {
    pt(i,0)=ptc(i,0)*x(0)+ptc(i,1)*y(0)+ptc(i,2)*z(0);
    pt(i,1)=ptc(i,0)*x(1)+ptc(i,1)*y(1)+ptc(i,2)*z(1);
    pt(i,2)=ptc(i,0)*x(2)+ptc(i,1)*y(2)+ptc(i,2)*z(2);
  }

  wf->updateVal(wfdata, sample);
	
  //We need to generate storage capable to store two electron moves
  wf->generateStorage( wfStore );
  sample->generateStorage( sampleStorage0 ); 
  sample->generateStorage( sampleStorage1 ); 

  //The two elecrons to be moved should have opposite spin
  int e0=0, e1=nelectrons-1;
  //assert( spin(e0) != spin(e1) );


  wf->saveUpdate(sample,e0,e1,wfStore);
  sample->saveUpdate( e0, sampleStorage0 );
  sample->saveUpdate( e1, sampleStorage1 );


  wf->getVal(wfdata, e0, wfval_old);
  sample->getElectronPos(e0,oldpos);
  sample->getElectronPos(e1,oldpos2);

  for ( int k=0; k<np_aver; k++) {
    
    transl(0)=dR*pt(k,0);
    transl(1)=dR*pt(k,1);
    transl(2)=dR*pt(k,2);

    for( int i=0; i<npoints; i++) {
      // shift makes sense only up to half the lattice vector, we are hitting
      // periodicity for larger distances
      
      sample->translateElectron(e0,transl);
      sample->translateElectron(e1,transl);
      wf->updateVal(wfdata, sample);
      
      wf->getVal(wfdata, e1, wfval_new);
      
      // how to handle the possibility of wf_old=0? Can it happen?
      dmtrx(i)+=wt(k)*wfval_new.sign(0)*wfval_old.sign(0)
	*exp(wfval_new.amp(0,0)-wfval_old.amp(0,0))*weight;
    }

    //The samples need to be restore because we're using the 
    //translateElectron function (important for k-points)
    sample->setElectronPos(e0,oldpos);
    sample->setElectronPos(e1,oldpos2);

    sample->restoreUpdate(e1, sampleStorage1 );
    sample->restoreUpdate(e0, sampleStorage0 );
    wf->restoreUpdate(sample,e0,e1,wfStore);
  }

  
  if(wfStore != NULL) delete wfStore;
  if(sampleStorage0 != NULL) delete sampleStorage0;
  if(sampleStorage1 != NULL) delete sampleStorage1;

  nsample+=1;
  wnsample+=weight;
  
}

//----------------------------------------------------------------------

void TBDM::write(string & log_label) { 

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
