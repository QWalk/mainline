#include "Average_quadrupole.h"
#include "jsontools.h"
//----------------------------------------------------------------------

void Average_quadrupole::evaluate(Wavefunction_data * wfdata, Wavefunction * wf,
                        System * sys, Sample_point * sample, Average_return & avg) 
{ 
  int ndim=3;
  Array2 <doublevar> quad(ndim,ndim,0.0);
  
  int nelectrons=sample->electronSize();
  Array1<doublevar> pos(ndim);
  for(int e=0; e< nelectrons; e++) {
    sample->getElectronPos(e,pos);
    doublevar r2=0;
    for(int d=0; d< ndim; d++) 
      r2+=pos(d)*pos(d);

    for(int d1=0; d1< ndim; d1++) { 
      for(int d2=0; d2< ndim; d2++) { 
        quad(d1,d2)+=-1.*pos(d1)*pos(d2)*3.;
      }
      quad(d1,d1)-=-1.*r2;
    }
  }

  int nions=sample->ionSize();
  for(int at=0; at < nions; at++) {
    sample->getIonPos(at,pos);
    doublevar charge=sample->getIonCharge(at);
    doublevar r2=0.0;
    for(int d=0; d< ndim; d++) 
      r2+=pos(d)*pos(d);
    
    for(int d1=0; d1< ndim; d1++) {
      for(int d2=0; d2< ndim; d2++) { 
        quad(d1,d2)+=charge*pos(d1)*pos(d2)*3.;
      }
      quad(d1,d1)-=charge*r2;
    }
  }
  
  //Serialization:
  // quad(i,j)= vals[i*n+j]
  avg.type="quadrupole";
  avg.vals.Resize(ndim*ndim*3);
  for(int i=0; i < ndim; i++) { 
    for(int j=0; j< ndim; j++) { 
      avg.vals(i*ndim+j)=quad(i,j);
      avg.vals(ndim*ndim+i*ndim+j)=cos(2*pi*gvecs(i,j)*quad(i,j));
      avg.vals(2*ndim*ndim+i*ndim+j)=sin(2*pi*gvecs(i,j)*quad(i,j));
    }
  }
  
} 
//----------------------------------------------------------------------

void Average_quadrupole::read(System * sys, Wavefunction_data * wfdata, vector
                    <string> & words) { 
  int ndim=3;
  unsigned int pos=0;
  Array2<doublevar> gveclat(ndim,ndim);
  if(!sys->getRecipLattice(gveclat)) { 
    doublevar len=10.0;
    readvalue(words,pos=0,len,"LENGTH");
    gveclat=0.0;
    for(int d=0; d< ndim; d++) gveclat(d,d)=1./len;
  }
  gvecs.Resize(ndim,ndim);
  gvecs=0.0;
  for(int d1=0; d1 < ndim; d1++) { 
    for(int d2=0; d2 < ndim; d2++) { 
// This is just for orthorhombic cells. 
      doublevar factor=0.5;
      if(d1==d2) factor=1.0;
      gvecs(d1,d2)=factor*gveclat(d1,d1)*gveclat(d2,d2);
      //for(int d3=0; d3< ndim; d3++) { 
      //  gvecs(d1,d2)+=gveclat(d1,d3)*gveclat(d2,d3);
      //}
    }
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
  for(int i=0; i < 3; i++) { 
    for(int j=0; j< 3; j++) gvecs(i,j)=atof(gtext[i*3+j].c_str());
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

