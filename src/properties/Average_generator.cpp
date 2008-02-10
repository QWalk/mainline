#include "Average_generator.h"
#include "qmc_io.h"


//-----------------------------------------------------------------------------
int decide_averager(string & label, Average_generator *& avg) { 
  if(caseless_eq(label, "DIPOLE") ){ 
    avg=new Average_dipole;
  }
  else { 
    error("Didn't understand ", label, " in Average_generator.");
  }
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

