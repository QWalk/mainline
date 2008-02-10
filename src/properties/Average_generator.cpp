#include "Average_generator.h"
#include "qmc_io.h"


//-----------------------------------------------------------------------------
int decide_averager(string & label, Average_generator *& avg) { 
  if(caseless_eq(label, "DIPOLE") ) 
    avg=new Average_dipole;
  else if(caseless_eq(label,"SK")) 
    avg=new Average_structure_factor;
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
  if(!readvalue(words, pos=0, np_side, "NGRID")) 
    np_side=5;
  
  npoints=np_side*np_side*np_side; 
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

//-----------------------------------------------------------------------------

void Average_structure_factor::evaluate(Wavefunction_data * wfdata, Wavefunction * wf,
                                        System * sys, Sample_point * sample, Average_return & avg) {
  avg.type="sk";
  nelectrons=sample->electronSize();
  
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
      << kpts(i,0) << "   " << kpts(i,1) << "   " << kpts(i,2) << endl;
  }
  
}
//############################################################################

