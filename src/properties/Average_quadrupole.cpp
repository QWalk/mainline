#include "Average_quadrupole.h"
#include "jsontools.h"
//----------------------------------------------------------------------

void Average_quadrupole::evaluate(Wavefunction_data * wfdata, Wavefunction * wf,
                        System * sys, Sample_point * sample, Average_return & avg) 
{ 
  int ndim=3;
  int nelectrons=sample->electronSize();
  int nions=sample->ionSize();
  Array1 <doublevar> charges(nelectrons+nions);
  Array2 <doublevar> positions(nelectrons+nions,ndim);
  Array1<doublevar> pos(ndim);

  for(int e=0; e< nelectrons; e++) {
    sample->getElectronPos(e,pos);
    charges(e)=-1;
    for(int d=0; d< 3; d++) 
      positions(e,d)=pos(d);
  }

  for(int at=0; at < nions; at++) {
    sample->getIonPos(at,pos);
    charges(at+nelectrons)=sample->getIonCharge(at);
    for(int d=0; d< 3; d++) 
      positions(at+nelectrons,d)=pos(d);
  }

  
  //Serialization:
  // quad(i,j)= vals[i*n+j]
  avg.type="quadrupole";
  avg.vals.Resize(ndim*ndim*3);
  avg.vals=0.0;
  for(int j=0; j < ndim; j++) { 
    for(int k=0; k< ndim; k++) { 
      for(doublevar shift=0; shift < 0.2; shift+=0.05) { 
      doublevar sum=0;      
      for(int p=0; p < nelectrons+nions; p++) { 
        avg.vals(j*ndim+k)+=charges(p)*positions(p,j)*positions(p,k);

        for(int a=0; a< ndim; a++) { 
          for(int b=0; b< ndim; b++) { 
            sum+=charges(p)*gvecs(j,a)*gvecs(k,b)
                       *(positions(p,a)+shift)
                       *(positions(p,b)+shift);
          }
        }
      }
      cout << "shift " << shift << " sum " << sum << " cos(sum) " 
        << cos(2*pi*sum) << endl;
      }
      //avg.vals(ndim*ndim+j*ndim+k)=cos(2*pi*sum);
      //avg.vals(2*ndim*ndim+j*ndim+k)=sin(2*pi*sum);
    }
  }

  
} 
//----------------------------------------------------------------------

void Average_quadrupole::read(System * sys, Wavefunction_data * wfdata, vector
                    <string> & words) { 
  int ndim=3;
  unsigned int pos=0;
  gvecs.Resize(3,3);
  if(!sys->getRecipLattice(gvecs)) { 
    doublevar len=10.0;
    readvalue(words,pos=0,len,"LENGTH");
    gvecs=0.0;
    for(int d=0; d< ndim; d++) gvecs(d,d)=1./len;
  }
  
}
//----------------------------------------------------------------------

void Average_quadrupole::write_init(string & indent, ostream & os) { 
  os << indent << "quadrupole" << endl;
  os << indent << "GVECS { ";
  for(int i=0 ;i< 3; i++) { 
    for(int j=0; j< 3; j++) {
      os << gvecs(i,j) << " ";
    }
  }
  os << " } " << endl;
  
} 

//----------------------------------------------------------------------

void Average_quadrupole::read(vector <string> & words) { 
  unsigned int pos=0;
  vector<string> gtext;
  readsection(words,pos,gtext,"GVECS");
  if(gtext.size()!=9) error("Expected 9 terms in GVECS");
  gvecs.Resize(3,3);
  int count=0;
  for(int i=0; i < 3; i++) { 
    for(int j=0; j< 3; j++) 
          gvecs(i,j)=atof(gtext[count++].c_str());
  }
  
} 

//----------------------------------------------------------------------

void Average_quadrupole::write_summary(Average_return &avg,
                                       Average_return &err, 
                                       ostream & os) {
  int ndim=3;
  for(int k=0; k < 3; k++) { 
    os << "Quadrupole moment" << endl;
    for(int i=0; i < ndim; i++) { 
      for(int j=0; j< ndim; j++) { 
        os << avg.vals(k*ndim*ndim+i*ndim+j) << " ";
      }
      os << endl;
    }
    os << "Error" << endl;
    for(int i=0; i < ndim; i++) { 
      for(int j=0; j< ndim; j++) { 
        os << err.vals(k*ndim*ndim+i*ndim+j) << " ";
      }
      os << endl;
    }
  }

} 
//----------------------------------------------------------------------

void Average_quadrupole::jsonOutput(Average_return &avg,
                                    Average_return &err, 
                                    ostream & os) { 
  os << "\"" << avg.type << "\":{" << endl;

  vector <string> nms;
  nms.push_back("quadrupole"); 
  nms.push_back("realpart"); 
  nms.push_back("imagpart");
  int ndim=3;
  os << "\"gvecs\":";
  jsonarray(os,gvecs);
  os << ",\n";

  Array2 <doublevar> tmp(ndim,ndim);
  
  for(int k=0; k < 3; k++) { 

    for(int i=0; i < ndim; i++) { 
      for(int j=0; j< ndim; j++) { 
        tmp(i,j)=avg.vals(k*ndim*ndim+i*ndim+j);
      }
    }
    os << "\"" << nms[k] << "\":";
    jsonarray(os,tmp);
    os << ",";

    for(int i=0; i < ndim; i++) { 
      for(int j=0; j< ndim; j++) { 
        tmp(i,j)=err.vals(k*ndim*ndim+i*ndim+j);
      }
    }
    os << "\"" << nms[k] << "_err" << "\":";
    jsonarray(os,tmp);
    if(k!=2) os << ",";
    
  }
  os << "}" << endl;  

} 


//----------------------------------------------------------------------


