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
      doublevar sum=0;      
      for(int p=0; p < nelectrons+nions; p++) { 
        avg.vals(j*ndim+k)+=charges(p)*positions(p,j)*positions(p,k);

        for(int a=0; a< ndim; a++) { 
          for(int b=0; b< ndim; b++) { 
            sum+=charges(p)*gvecs(j,a)*gvecs(k,b)
                       *(positions(p,a))
                       *(positions(p,b));
          }
        }
      }
      avg.vals(ndim*ndim+j*ndim+k)=cos(2*pi*sum);
      avg.vals(2*ndim*ndim+j*ndim+k)=sin(2*pi*sum);
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
//######################################################################
//

void Average_joint_pbc::evaluate(Wavefunction_data * wfdata, Wavefunction * wf,
                        System * sys, Sample_point * sample, Average_return & avg) 
{ 
  int ndim=3;  
  int nelectrons=sample->electronSize();  
  Array2 <doublevar> pol(nelectrons, ndim,0.0);
  Array3 <doublevar> quadplus(nelectrons,ndim,ndim,0.0),quadminus(nelectrons,ndim,ndim,0.0);
  Array1 <doublevar> pos(ndim);
  for(int e=0; e< nelectrons; e++) {
    sample->getElectronPos(e,pos);
    Array1 <doublevar> reduced(ndim,0.0);
    for(int i=0; i< ndim; i++ ) { 
      for(int j=0; j< ndim; j++) { 
        reduced(i)+=gvecs(i,j)*pos(j);
      }
    }

    for(int i=0; i< ndim; i++) { 
      pol(e,i)+=reduced(i);
      for(int j=i; j< ndim; j++) { 
        quadplus(e,i,j)=reduced(i)+reduced(j);
        quadminus(e,i,j)=reduced(i)-reduced(j);
      }
    }
  }

  avg.type="joint_pbc";
  int ntot=2*(nelectrons*ndim*(ndim+1)+nelectrons*ndim);
  avg.vals.Resize(ntot);

  int count=0;
  for(int e=0; e < nelectrons; e++) { 
    for(int i=0; i < ndim; i++) {
      avg.vals(count++)=cos(pol(e,i));
      avg.vals(count++)=sin(pol(e,i));
      
    }
  }
  
  for(int e=0; e < nelectrons; e++) { 
    for(int i=0; i < ndim; i++) {
      for(int j=i; j< ndim; j++) { 
        avg.vals(count++)=cos(quadplus(e,i,j));
        avg.vals(count++)=sin(quadplus(e,i,j));
        
      }
    }
  }
  for(int e=0; e < nelectrons; e++) { 
    for(int i=0; i < ndim; i++) {
      for(int j=i; j< ndim; j++) { 
        avg.vals(count++)=cos(quadminus(e,i,j));
        avg.vals(count++)=sin(quadminus(e,i,j));
        
      }
    }
  }
  
} 

//----------------------------------------------------------------------

void Average_joint_pbc::read(System * sys, Wavefunction_data * wfdata, vector
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

void Average_joint_pbc::write_init(string & indent, ostream & os) { 
  os << indent << "joint_pbc" << endl;
  os << indent << "GVECS { ";
  for(int i=0 ;i< 3; i++) { 
    for(int j=0; j< 3; j++) {
      os << gvecs(i,j) << " ";
    }
  }
  os << " } " << endl;
  
} 

//----------------------------------------------------------------------

void Average_joint_pbc::read(vector <string> & words) { 
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

void Average_joint_pbc::write_summary(Average_return &avg,
                                       Average_return &err, 
                                       ostream & os) {
} 
//----------------------------------------------------------------------

void Average_joint_pbc::jsonOutput(Average_return &avg,
                                    Average_return &err, 
                                    ostream & os) { 


  int ndim=3;
  int ntot=avg.vals.GetDim(0);
  int nelectrons=ntot/(ndim*(ndim+1)+ndim)/2;
  int nd=ndim*(ndim+1)/2;
  Array2 <doublevar> pol(nelectrons, ndim,0.0),poli(nelectrons, ndim,0.0);
  Array2 <doublevar> quadplus(nelectrons,nd,0.0),quadminus(nelectrons,nd,0.0);
  Array2 <doublevar> quadplusi(nelectrons,nd,0.0),quadminusi(nelectrons,nd,0.0);
  
  Array2 <doublevar> pol_err(nelectrons, ndim,0.0);
  Array2 <doublevar> quadplus_err(nelectrons,nd,0.0),quadminus_err(nelectrons,nd,0.0);
  
  int count=0;
  for(int e=0; e < nelectrons; e++) { 
    for(int i=0; i < ndim; i++) {
      pol(e,i)=avg.vals(count);
      pol_err(e,i)=err.vals(count);
      count++;
      poli(e,i)=avg.vals(count);
      count++;
    }
  }
  
  for(int e=0; e < nelectrons; e++) { 
    int cij=0;
    for(int i=0; i < ndim; i++) {
      for(int j=i; j< ndim; j++) { 
        quadplus(e,cij)=avg.vals(count);
        quadplus_err(e,cij)=err.vals(count);
        count++;
        quadplusi(e,cij)=avg.vals(count);

        count++;
        cij++;

      }
    }
  }
  for(int e=0; e < nelectrons; e++) { 
    int cij=0;
    for(int i=0; i < ndim; i++) {
      for(int j=i; j< ndim; j++) { 
        quadminus(e,cij)=avg.vals(count);
        quadminus_err(e,cij)=err.vals(count);
        count++;
        quadminusi(e,cij)=avg.vals(count);
        count++;
        cij++;
      }
    }
  }
  os << "\"" << avg.type << "\":{" << endl;
  
  os << "\"" << "pol" << "\":";  
  jsonarray(os,pol);
  os << "," << endl;
  os << "\"" << "pol_err" << "\":";  
  jsonarray(os,pol_err);
  os << "}" << endl;  
  os << "\"" << "poli" << "\":";  
  jsonarray(os,poli);
  os << "," << endl;
  
  os << "\"" << "quadplus" << "\":";  
  jsonarray(os,quadplus);
  os << "," << endl;
  os << "\"" << "quadplus_err" << "\":";  
  jsonarray(os,quadplus_err);
  os << "}" << endl;  
  os << "\"" << "quadplusi" << "\":";  
  jsonarray(os,quadplusi);
  os << "," << endl;
  
  
  os << "\"" << "quadminus=" << "\":";  
  jsonarray(os,quadminus);
  os << "," << endl;

                       
  os << "\"" << "quadminusi" << "\":";  
  jsonarray(os,quadminusi);
  os << "," << endl;
  os << "\"" << "quadminus_err" << "\":";  
  jsonarray(os,quadplus_err);
  os << "}" << endl;  
}
