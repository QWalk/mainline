#include "Average_generator.h"
#include "qmc_io.h"


//-----------------------------------------------------------------------------
int decide_averager(string & label, Average_generator *& avg) { 
  if(caseless_eq(label, "DIPOLE") ) 
    avg=new Average_dipole;
  else if(caseless_eq(label,"SK")) 
    avg=new Average_structure_factor;
  else if(caseless_eq(label, "GR"))
    avg=new Average_twobody_correlation;
  else if(caseless_eq(label, "MANYBODY_POL"))
    avg=new Average_manybody_polarization;
  else 
    error("Didn't understand ", label, " in Average_generator.");
  
  return 1;
}

//-----------------------------------------------------------------------------


int allocate(vector<string> & words, System * sys, Wavefunction_data * wfdata, Average_generator * & avg) { 
  decide_averager(words[0], avg);
  avg->read(sys, wfdata, words);
  return 1;
}
//-----------------------------------------------------------------------------


int allocate(vector<string> & words, Average_generator * & avg) { 
  decide_averager(words[0], avg);
  avg->read(words);
  return 1;
}

//############################################################################
//Average_dipole

void Average_dipole::read(System * sys, Wavefunction_data * wfdata, vector <string> & words) { 
}

void Average_dipole::evaluate(Wavefunction_data * wfdata, Wavefunction * wf,
                              System * sys, Sample_point * sample, Average_return & avg ) { 
  avg.type="dipole";
  int ndim=3;
  avg.vals.Resize(ndim);
  avg.vals=0.0;
  int nelectrons=sample->electronSize();
  Array1 <doublevar> pos(ndim);
  for(int e=0; e< nelectrons; e++) {
    sample->getElectronPos(e,pos);
    for(int d=0; d< ndim; d++) {
      avg.vals(d)-=pos(d);
    }
  }
  int nions=sample->ionSize();
  for(int at=0; at < nions; at++) {
    sample->getIonPos(at,pos);
    doublevar charge=sample->getIonCharge(at);
    for(int d=0; d< ndim; d++) {
      avg.vals(d)+=charge*pos(d);
    }
  }
  
}

void Average_dipole::write_init(string & indent, ostream & os) { 
  os << indent << "dipole \n";
}

void Average_dipole::read(vector <string> & words) { 
  
}

void Average_dipole::write_summary(Average_return & avg, Average_return & err, ostream & os) {
  int ndim=avg.vals.GetDim(0);
  assert(ndim <= err.vals.GetDim(0));
  //Could put this in Debye if we want to be nice.
  os << "Dipole moment (a.u.) \n";
  for(int d=0; d< ndim; d++) { 
    if(d==0) os << "x ";
    else if(d==1) os << "y ";
    else if(d==2) os << "z ";
    os << avg.vals(d) << " +/- " << err.vals(d) << endl;
  }
    
}
//############################################################################


void Average_structure_factor::read(System * sys, Wavefunction_data * wfdata, vector <string> & words) {
  unsigned int pos=0;
  int np_side;
  if(!readvalue(words, pos=0, np_side, "NGRID")) 
    np_side=5;
  
  Array2 <doublevar> gvec;
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
  
  int direction;
  if(readvalue(words, pos=0, direction, "DIRECTION")) { 
    npoints=np_side;
    kpts.Resize(npoints, 3);
    kpts=0;
    for(int ix=0; ix < np_side; ix++) { 
      for(int i=0; i< 3; i++) 
        kpts(ix,i)=2*pi*gvec(direction,i)*ix;
    }
  }
  else { 
    npoints=np_side*np_side*np_side; 

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
}  

//-----------------------------------------------------------------------------

void Average_structure_factor::evaluate(Wavefunction_data * wfdata, Wavefunction * wf,
                                        System * sys, Sample_point * sample, Average_return & avg) {
  avg.type="sk";
  int nelectrons=sample->electronSize();
  
  Array1 <doublevar> pos(3);
  avg.vals.Resize(npoints);
  avg.vals=0;
  for(int p=0; p < npoints; p++) { 
    doublevar sum_cos=0, sum_sin=0;
    
    for(int e=0; e< nelectrons; e++) {
      sample->getElectronPos(e,pos);
      doublevar dot=0;
      for(int d=0; d< 3; d++) dot+=pos(d)*kpts(p,d);
      sum_cos+=cos(dot);
      sum_sin+=sin(dot);
    }
    avg.vals(p)+=(sum_cos*sum_cos+sum_sin*sum_sin)/nelectrons;
    
  }
  
}

//-----------------------------------------------------------------------------

void Average_structure_factor::write_init(string & indent, ostream & os) { 
  os << indent << "sk\n";
  for(int i=0; i< npoints; i++) { 
    os << indent << "kpoint { " << kpts(i,0) << " " << kpts(i,1) << " " 
    << kpts(i,2) << " } " <<endl;
  }
}
//-----------------------------------------------------------------------------

void Average_structure_factor::read(vector <string> & words) { 
  vector <vector <string> > kpttext;
  vector <string> tmp;
  unsigned int pos=0;
  while(readsection(words, pos, tmp, "KPOINT")) kpttext.push_back(tmp);
  npoints=kpttext.size();
  kpts.Resize(npoints,3);
  for(int i=0; i< npoints; i++) { 
    for(int d=0; d< 3; d++) { 
      kpts(i,d)=atof(kpttext[i][d].c_str());
    }
  }
}
//-----------------------------------------------------------------------------

void Average_structure_factor::write_summary(Average_return & avg, Average_return & err, ostream & os) { 
  os << "Structure factor \n";
  os << "    k  s(k) s(k)err  kx  ky  kz" << endl;
  assert(avg.vals.GetDim(0) >=npoints);
  assert(err.vals.GetDim(0) >=npoints);
  
  for(int i=0; i< npoints; i++) {
    doublevar sk=avg.vals(i);
    doublevar skerr=err.vals(i);
    doublevar r=0;
    for(int d=0; d< 3; d++) r+=kpts(i,d)*kpts(i,d);
    r=sqrt(r);
    os << "sk_out " <<  r << "   " << sk << "   " << skerr 
      << "  " << kpts(i,0) << "   " << kpts(i,1) << "   " << kpts(i,2) << endl;
  }
  
}
//############################################################################

void Average_twobody_correlation::read(System * sys, Wavefunction_data * wfdata, 
                                       vector <string> & words) { 
  doublevar range;
  unsigned int pos=0;
  if(!readvalue(words,pos=0,resolution, "RESOLUTION"))
    resolution=0.1;
  if(!readvalue(words,pos=0,range,"RANGE"))
    range=10.0;
  if(readvalue(words, pos=0, direction, "DIRECTION")) { 
    direction+=2;
  }
  else
    direction=0;
  
  npoints=int(range/resolution)+1;  
  
}

//-----------------------------------------------------------------------------

void Average_twobody_correlation::evaluate(Wavefunction_data * wfdata, Wavefunction * wf,
                                        System * sys, Sample_point * sample, Average_return & avg) {
  avg.type="gr";
  int nelectrons=sample->electronSize();
  avg.vals.Resize(2*npoints); //for like and unlike
  avg.vals=0;
  int nup=sys->nelectrons(0);
  int ndown=nelectrons-nup;
  sample->updateEEDist();
  Array1 <doublevar> dist(5);
  for(int e=0; e< nup; e++) {
    for(int e2=e+1; e2 < nup; e2++) {
      sample->getEEDist(e,e2,dist);
      int place=int(fabs(dist(direction))/resolution);
      if(place < npoints) { 
        avg.vals(place)+=1;
      }
    }
    for(int e2=nup; e2 < nelectrons; e2++) { 
      sample->getEEDist(e,e2,dist);
      int place=int(fabs(dist(direction))/resolution);
      if(place < npoints) { 
        avg.vals(npoints+place)+=1;
      }
    }
  }
  
  for(int e=nup; e < nelectrons; e++) { 
    for(int e2=e+1; e2< nelectrons; e2++) { 
      sample->getEEDist(e,e2,dist);
      int place=int(fabs(dist(direction))/resolution);
      if(place < npoints) { 
        avg.vals(place) += 1.0;
      }
    }
  }
  
  if(direction==0) { 
    for(int i=0; i< npoints; i++) { 
      double r=resolution*i;
      avg.vals(i)/=4*pi*resolution*r*r;
      avg.vals(npoints+i)/=4*pi*resolution*r*r;
    }
  }
  
  
}

//-----------------------------------------------------------------------------

void Average_twobody_correlation::write_init(string & indent, ostream & os) { 
  os << indent << "gr\n";
  os << indent << "npoints " << npoints << endl;
  os << indent << "resolution " << resolution << endl;
  os << indent << "direction " << direction << endl;
}
//-----------------------------------------------------------------------------

void Average_twobody_correlation::read(vector <string> & words) { 
  unsigned int pos=0;
  readvalue(words, pos=0, resolution, "resolution");
  readvalue(words, pos=0, npoints, "npoints");
  readvalue(words, pos=0, direction, "direction");
}
//-----------------------------------------------------------------------------

void Average_twobody_correlation::write_summary(Average_return & avg, Average_return & err, ostream & os) { 
  os << "Electron correlation for like and unlike spins\n";
  os << "    r  g(r) sigma(g(r))   g(r)  sigma(g(r))" << endl;
  assert(avg.vals.GetDim(0) >=2*npoints);
  assert(err.vals.GetDim(0) >=2*npoints);
  
  for(int i=0; i< npoints; i++) {
    os << "gr_out " << i*resolution << " " << avg.vals(i) << " " << err.vals(i) 
    << "  " << avg.vals(i+npoints) << "  " << err.vals(i+npoints) << endl;
  }
  
}


//############################################################################


void Average_manybody_polarization::read(System *sys, Wavefunction_data * wfdata, vector <string> & words) { 
  gvec.Resize(3,3);
  if(!sys->getRecipLattice(gvec) ) { 
    error("The manybody polarization operator works only for periodic systems!");
  }
}

//-----------------------------------------------------------------------------
void Average_manybody_polarization::evaluate(Wavefunction_data * wfdata, Wavefunction * wf,
                                             System * sys, Sample_point * sample, Average_return & avg) {
  avg.type="manybody_pol";
  int nelectrons=sample->electronSize();
  avg.vals.Resize(6);
  Array1 <doublevar> sum(3,0.0);
  Array1 <doublevar> pos(3);
  for(int e=0; e< nelectrons; e++) { 
    sample->getElectronPos(e,pos);
    for(int i=0; i< 3; i++) { 
      for(int d=0; d< 3; d++) { 
        sum(i)+=gvec(i,d)*pos(d);
      }
    }
  }
  for(int i=0; i< 3; i++) { 
    avg.vals(2*i)=cos(2*pi*sum(i));
    avg.vals(2*i+1)=sin(2*pi*sum(i));
  }
}
//-----------------------------------------------------------------------------

void Average_manybody_polarization::write_init(string & indent, ostream & os) { 
  os << indent << "manybody_pol" << endl;
  os << indent << "gvec { ";
  for(int i=0; i< 3; i++) 
    for(int j=0; j< 3; j++) os << gvec(i,j) << "  ";
  os << "}\n";
}
//-----------------------------------------------------------------------------
void Average_manybody_polarization::read(vector <string> & words) { 
  unsigned int pos=0;
  vector <string> gvec_sec;
  readsection(words, pos=0, gvec_sec, "gvec");
  int count=0;
  gvec.Resize(3,3);
  for(int i=0; i< 3; i++)
    for(int j=0; j< 3; j++) gvec(i,j)=atof(gvec_sec[count++].c_str());
  
}
//-----------------------------------------------------------------------------


void Average_manybody_polarization::write_summary(Average_return & avg, Average_return & err, ostream & os) { 
  os << "Manybody polarization operator " << endl;
  for(int i=0; i< 3; i++) { 
    os << avg.vals(2*i) << " + i " << avg.vals(2*i+1) << " +/- " << err.vals(2*i) << " + i " << err.vals(2*i+1) << endl;
  }
}


