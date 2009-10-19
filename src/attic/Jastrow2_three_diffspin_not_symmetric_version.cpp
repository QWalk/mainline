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

#include "Jastrow2_wf.h"
#include "qmc_io.h"
#include "System.h"
#include "Sample_point.h"
#include "Jastrow2_three_diffspin.h"
#include "Jastrow2_one.h"

//--------------------------------------------------------  

int Jastrow_threebody_piece_diffspin::make_default_list() { 
  klm.Resize(12,3);
  klm(0,0)=1; klm(0,1)=1; klm(0,2)=0;
  klm(1,0)=1; klm(1,1)=0; klm(1,2)=1;
  klm(2,0)=1; klm(2,1)=1; klm(2,2)=1;
  klm(3,0)=2; klm(3,1)=2; klm(3,2)=0;
  klm(4,0)=2; klm(4,1)=0; klm(4,2)=1;
  klm(5,0)=2; klm(5,1)=0; klm(5,2)=2;
  klm(6,0)=2; klm(6,1)=2; klm(6,2)=2;
  klm(7,0)=3; klm(7,1)=3; klm(7,2)=0;
  klm(8,0)=3; klm(8,1)=0; klm(8,2)=2;
  klm(9,0)=3; klm(9,1)=3; klm(9,2)=2;
  klm(10,0)=1; klm(10,1)=2; klm(10,2)=2;
  klm(11,0)=2; klm(11,1)=3; klm(11,2)=2;
  return 12;

}

void get_onebody_parms_diffspin(vector<string> & words, vector<string> & atomnames,
				Array2 <doublevar> & unique_parameters,
				int s,
				Array1 <int> & nparms,
				vector <string> & parm_labels,
				Array2 <int> & linear_parms,
                                int & counter,
				Array1 <int> & parm_centers) {
  vector < vector < string> > parmtxt;
  vector <string> parmtmp;
  unsigned int pos=0;
  if(s==0){
    while(readsection(words, pos,parmtmp, "LIKE_COEFFICIENTS"))
      parmtxt.push_back(parmtmp);
    if(parmtxt.size()==0)
      error("Didn't find LIKE COEFFICIENTS in one of threebody spin Jastrow section");
  }
  else{
    while(readsection(words, pos,parmtmp, "UNLIKE_COEFFICIENTS"))
      parmtxt.push_back(parmtmp);
    if(parmtxt.size()==0)
      error("Didn't find UNLIKE COEFFICIENTS in one of threebody spin Jastrow section");
  }
    
  int nsec=parmtxt.size();
  int maxsize=0;
  nparms.Resize(nsec);
  for(int i=0; i< nsec; i++)
    nparms(i)=parmtxt[i].size()-1;

  for(int i=0; i< nsec; i++)
    if(maxsize < nparms(i)) maxsize=nparms(i);

  unique_parameters.Resize(nsec, maxsize);
  linear_parms.Resize(nsec, maxsize);
  //cout <<"spin "<<s<<" counter "<<counter<<endl;
  for(int i=0; i< nsec; i++) {
    for(int j=0; j< nparms(i); j++) {
      
      unique_parameters(i,j)=atof(parmtxt[i][j+1].c_str());
      linear_parms(i,j)=counter++;
    }
  }

  parm_labels.resize(nsec);
  for(int i=0; i< nsec; i++)
    parm_labels[i]=parmtxt[i][0];

  int natoms=atomnames.size();
  parm_centers.Resize(natoms);
  parm_centers=-1;

  for(int i=0; i< nsec; i++) {
    int found=0;
    for(int j=0; j< natoms; j++) {
      if(atomnames[j]==parm_labels[i]) {
        parm_centers(j)=i;
        found=1;
      }
    }
    if(!found)
      error("Couldn't find a matching center for ", parm_labels[i]);
  }

  doublevar scale;
  if(readvalue(words, pos=0, scale, "SCALE")) {
    for(int i=0; i< unique_parameters.GetDim(0); i++) {
      for(int j=0; j< nparms(i); j++) {
        unique_parameters(i,j)*=scale;
      }
    }
  }
}

//--------------------------------------------------------  




void Jastrow_threebody_piece_diffspin::set_up(vector <string> & words, 
                                   System * sys) {
  //cout <<"Start Jastrow_threebody_piece_diffspin::set_up"<<endl;                                 
  vector<string> atomnames;
  sys->getAtomicLabels(atomnames);
  unique_parameters_spin.Resize(2);
  linear_parms.Resize(2);
  int tmpcounter=0;
  for(int s=0;s<unique_parameters_spin.GetSize();s++){
    get_onebody_parms_diffspin(words, atomnames, unique_parameters_spin(s),s,
			       _nparms, parm_labels,
			       linear_parms(s), tmpcounter, parm_centers);
  }

  


  if(unique_parameters_spin(0).GetDim(0)!=unique_parameters_spin(1).GetDim(0))
    error("Need same amount of LIKE_COEFFICIENTS and UNLIKE_COEFFICIENTS sections in three body");
  
  if(unique_parameters_spin(0).GetDim(1)!=unique_parameters_spin(1).GetDim(1))
    error("Need same amount of LIKE_COEFFICIENTS and UNLIKE_COEFFICIENTS parameters in three body");
  

  unsigned int pos=0;  
  int maxnparms=0;
  vector <string> klm_list;
  if(readsection(words,pos=0,klm_list,"SM_TERMS")) { 
    if(klm_list.size()%3 !=0) 
      error("SM_TERMS must contain a multiple of 3 elements");
    maxnparms=klm_list.size()/3;
    int counter=0;
    klm.Resize(maxnparms,3);
    for(int i=0; i< maxnparms; i++) { 
      for(int j=0; j< 3; j++) {
	klm(i,j)=atoi(klm_list[counter++].c_str());
      }
    }
  }
  else maxnparms=make_default_list();

  
  


  int natoms=atomnames.size();
  eebasis_max=0;
  eibasis_max.Resize(natoms);
  eibasis_max=0;
  for(int at=0; at< natoms; at++) {
    if(nparms(at) > maxnparms) 
      error("too many three-body spin parameters on ", atomnames[at]);
    for(int p=0; p < nparms(at); p++) { 
      if(eibasis_max(at) <= klm(p,0)) eibasis_max(at)=klm(p,0)+1;
      if(eibasis_max(at) <= klm(p,1)) eibasis_max(at)=klm(p,1)+1;
      if(eebasis_max <= klm(p,2)) eebasis_max=klm(p,2)+1;
    }
  }
  
  freeze=haskeyword(words, pos=0, "FREEZE");
  nspin_up=sys->nelectrons(0);
  nelectrons=sys->nelectrons(0)+sys->nelectrons(1);
  //cout <<"End Jastrow_threebody_piece_diffspin::set_up"<<endl;  

}
//--------------------------------------------------------------------------

int Jastrow_threebody_piece_diffspin::writeinput(string & indent, ostream & os) {
  if(freeze) os << indent << "FREEZE" << endl;

  
  for(int i=0; i< unique_parameters_spin(0).GetDim(0); i++) {
    os << indent << "LIKE_COEFFICIENTS { " << parm_labels[i] << "  ";
    for(int j=0; j< _nparms(i); j++)
      os << unique_parameters_spin(0)(i,j) << "  ";
    os << " } " << endl;
  }
  
  for(int i=0; i< unique_parameters_spin(1).GetDim(0); i++) {
    os << indent << "UNLIKE_COEFFICIENTS { " << parm_labels[i] << "  ";
    for(int j=0; j< _nparms(i); j++)
      os << unique_parameters_spin(1)(i,j) << "  ";
    os << " } " << endl;
  }
  
  os << indent << "SM_TERMS {\n";
  for(int i=0; i< klm.GetDim(0); i++) {
    os << indent << "  ";
    for(int j=0; j< klm.GetDim(1); j++) { 
      os << klm(i,j) << "  ";
    }
    os << endl;
  }
  os << indent << "}\n";

  return 1;

}

//--------------------------------------------------------------------------

int Jastrow_threebody_piece_diffspin::showinfo(string & indent, ostream & os) {
  // cout <<"Showinfo:start"<<endl;
  os << indent << "Atom   Like Coefficients " << endl;
  for(int i=0; i< unique_parameters_spin(0).GetDim(0); i++) {
    os << indent << parm_labels[i] << "  ";
    for(int j=0; j< _nparms(i); j++)
      os << unique_parameters_spin(0)(i,j) << "  ";
    os << endl;
  }
  os << indent << "Atom   Unlike Coefficients " << endl;
  for(int i=0; i< unique_parameters_spin(1).GetDim(0); i++) {
    os << indent << parm_labels[i] << "  ";
    for(int j=0; j< _nparms(i); j++)
      os << unique_parameters_spin(1)(i,j) << "  ";
    os << endl;
  }
  return 1;
  //cout <<"Showinfo:end"<<endl;
}


//--------------------------------------------------------------------------

void Jastrow_threebody_piece_diffspin::updateLap(int e,
                 const Array4 <doublevar> & eibasis,
                 const Array3 <doublevar> & eebasis,
                 Array3 <doublevar> & lap) {
  //cout << "Jastrow_threebody_piece::updateLap" << endl;
  assert(lap.GetDim(2) >= 5);
  assert(lap.GetDim(0) >=2);
  
  assert(eibasis.GetDim(1) >= parm_centers.GetDim(0));
  int natoms=parm_centers.GetDim(0);
  int nelectrons=eebasis.GetDim(0);
  assert(lap.GetDim(1) >= nelectrons);
  Array1 <doublevar> vkl_e(5); //derivative wrt e of a_e a_j + a_j a_e
  Array1 <doublevar> vkl_j(5); //wrt j
  
  int s;
  const doublevar tiny=1e-14;
  //cout << "updateLap " << endl;
  for(int at=0; at < natoms; at++) {
    int p=parm_centers(at);
    
    for(int i=0; i< _nparms(p); i++) {
      doublevar parm;
      int k=klm(i,0), el=klm(i,1), m=klm(i,2);
      if(fabs(eibasis(e,at,k,0)) > tiny
	 || fabs(eibasis(e,at,el,0)) > tiny) { 
	for(int j=0; j< e; j++) {
	  s=(j < nspin_up) != (e < nspin_up);
	  parm=unique_parameters_spin(s)(p,i);
	  doublevar vkl=parm*(eibasis(e,at,k,0)*eibasis(j,at,el,0)
			      //+eibasis(j,at,k,0)*eibasis(e,at,el,0)
			      );
	  
	  
	  for(int d=1; d< 5; d++) {
	    vkl_e[d]=parm*(eibasis(e,at,k,d)*eibasis(j,at,el,0)
			   //+eibasis(j,at,k,0)*eibasis(e,at,el,d)
			   );
	    vkl_j[d]=parm*(eibasis(e,at,k,0)*eibasis(j,at,el,d)
			   //+eibasis(j,at,k,d)*eibasis(e,at,el,0)
			   );
	  } 
	  lap(0,j,0)+=vkl*eebasis(j,m,0);
	  doublevar dot_e=0, dot_j=0;
	  for(int d=1; d< 4; d++) {
	    lap(0,j,d)+=vkl_e[d]*eebasis(j,m,0);
	    lap(0,j,d)-=vkl*eebasis(j,m,d);
	    lap(1,j,d)+=vkl_j[d]*eebasis(j,m,0);
	    lap(1,j,d)+=vkl*eebasis(j,m,d);
	    
	    dot_e+=vkl_e[d]*eebasis(j,m,d);
	    dot_j+=vkl_j[d]*eebasis(j,m,d);
	  }
	  
	  
	  lap(0,j,4)+=vkl*eebasis(j,m,4);
	  lap(1,j,4)+=vkl*eebasis(j,m,4);
	  
	  lap(0,j,4)+=vkl_e[4]*eebasis(j,m,0);
	  lap(1,j,4)+=vkl_j[4]*eebasis(j,m,0);
	  
	  lap(0,j,4)-=2*dot_e;
	  lap(1,j,4)+=2*dot_j;
	}
	for(int j=e+1; j< nelectrons; j++) {
	  s=(j < nspin_up) != (e < nspin_up);
	  parm=unique_parameters_spin(s)(p,i);
	  doublevar vkl=parm*(//eibasis(e,at,k,0)*eibasis(j,at,el,0)
			      +eibasis(j,at,k,0)*eibasis(e,at,el,0)
			      );
	  for(int d=1; d< 5; d++) {
	    vkl_e[d]=parm*(//eibasis(e,at,k,d)*eibasis(j,at,el,0)
			   +eibasis(j,at,k,0)*eibasis(e,at,el,d)
			   );
	    vkl_j[d]=parm*(//eibasis(e,at,k,0)*eibasis(j,at,el,d)
			   +eibasis(j,at,k,d)*eibasis(e,at,el,0)
			   );
	  } 
	  
	  lap(0,j,0)+=vkl*eebasis(j,m,0);
	  doublevar dot_e=0, dot_j=0;
	  for(int d=1; d< 4; d++) {
	    lap(0,j,d)+=vkl_e[d]*eebasis(j,m,0);
	    lap(0,j,d)+=vkl*eebasis(j,m,d);
	    lap(1,j,d)+=vkl_j[d]*eebasis(j,m,0);
	    lap(1,j,d)-=vkl*eebasis(j,m,d);
	    
	    dot_e+=vkl_e[d]*eebasis(j,m,d);
	    dot_j+=vkl_j[d]*eebasis(j,m,d);
	  }
	  lap(0,j,4)+=vkl*eebasis(j,m,4);
	  lap(1,j,4)+=vkl*eebasis(j,m,4);
	  
	  lap(0,j,4)+=vkl_e[4]*eebasis(j,m,0);
	  lap(1,j,4)+=vkl_j[4]*eebasis(j,m,0);
	  
	  lap(0,j,4)+=2*dot_e;
	  lap(1,j,4)-=2*dot_j;
	}
      }
    }
  }
  
  
}

//-----------------------------------------------------------
void Jastrow_threebody_piece_diffspin::updateVal(int e,
                                         const Array4 <doublevar> & eibasis,
                                         const Array3 <doublevar> & eebasis,
                                         Array1 <doublevar> & updated_val) {
  
  assert(eibasis.GetDim(1) >= parm_centers.GetDim(0));
  int natoms=parm_centers.GetDim(0);
  int nelectrons=eebasis.GetDim(0);
  int s;

  const doublevar tiny=1e-14;
  //cout << "updateLap " << endl;
  for(int at=0; at < natoms; at++) {
    int p=parm_centers(at);
    for(int i=0; i< _nparms(p); i++) {
      int k=klm(i,0), el=klm(i,1), m=klm(i,2);
      if(fabs(eibasis(e,at,k,0)) > tiny
	 || fabs(eibasis(e,at,el,0)) > tiny) { 
	for(int j=0; j< e; j++) {
	  s=(j < nspin_up) != (e < nspin_up); 
	  doublevar vkl=unique_parameters_spin(s)(p,i)*(eibasis(e,at,k,0)*eibasis(j,at,el,0)
							//+eibasis(j,at,k,0)*eibasis(e,at,el,0)
							);
	  updated_val(j)+=vkl*eebasis(j,m,0);
	}
	for(int j=e+1; j< nelectrons; j++) {
	  s=(j < nspin_up) != (e < nspin_up);
	  doublevar vkl=unique_parameters_spin(s)(p,i)*(//eibasis(e,at,k,0)*eibasis(j,at,el,0)
							+eibasis(j,at,k,0)*eibasis(e,at,el,0));
	  updated_val(j)+=vkl*eebasis(j,m,0);
	}
      }
    }
  } 
}



//--------------------------------------------------------------------------
//start: added for ei back-flow
void Jastrow_threebody_piece_diffspin::updateVal_E_I(int e,
						     const Array4 <doublevar> & eibasis,
						     const Array3 <doublevar> & eebasis,
						     Array3 <doublevar> & lap) {
  //cout << "Jastrow_threebody_piece::updateLap_E_I" << endl;
   if(eebasis.GetDim(1)<eebasis_max)
     error("Need more eebasis functions in three-body spin");
  int natoms=parm_centers.GetDim(0);
  int nelectrons=eebasis.GetDim(0);

  assert(lap.GetDim(2) >= 5);
  assert(lap.GetDim(0) >= nelectrons);
  assert(eibasis.GetDim(1) >= parm_centers.GetDim(0));
  assert(lap.GetDim(1) >= natoms);

   
  int s;
  const doublevar tiny=1e-14;
  //cout << "updateLap " << endl;
  for(int at=0; at < natoms; at++) {
    int p=parm_centers(at);
    
    for(int i=0; i< _nparms(p); i++) {
      doublevar parm;
      int k=klm(i,0), el=klm(i,1), m=klm(i,2);
      if(fabs(eibasis(e,at,k,0)) > tiny
	 || fabs(eibasis(e,at,el,0)) > tiny) { 
	for(int j=0; j< e; j++) {
	  s=(j < nspin_up) != (e < nspin_up);
	  parm=unique_parameters_spin(s)(p,i);
	  doublevar vkl=parm*(eibasis(e,at,k,0)*eibasis(j,at,el,0)
			      //+eibasis(j,at,k,0)*eibasis(e,at,el,0)
			      );
	  lap(e,at,0)+=vkl*eebasis(j,m,0);
	}
	for(int j=e+1; j< nelectrons; j++) {
	  s=(j < nspin_up) != (e < nspin_up);
	  parm=unique_parameters_spin(s)(p,i);
	  doublevar vkl=parm*(eibasis(e,at,k,0)*eibasis(j,at,el,0)
			      //+eibasis(j,at,k,0)*eibasis(e,at,el,0)
			      );
	  lap(e,at,0)+=vkl*eebasis(j,m,0);
	}
      }
    }
  }
}



//--------------------------------------------------------------------------

void Jastrow_threebody_piece_diffspin::updateLap_E_I(int e,
						     const Array4 <doublevar> & eibasis,
						     const Array3 <doublevar> & eebasis,
						     Array3 <doublevar> & lap) {
  //cout << "Jastrow_threebody_piece::updateLap_E_I" << endl;
 
  int natoms=parm_centers.GetDim(0);
  int nelectrons=eebasis.GetDim(0);
  if(eebasis.GetDim(1)<eebasis_max)
    error("Need more eebasis functions in three-body spin");

  assert(lap.GetDim(2) >= 5);
  assert(lap.GetDim(0) >= nelectrons);
  assert(eibasis.GetDim(1) >= parm_centers.GetDim(0));
  assert(lap.GetDim(1) >= natoms);

  Array1 <doublevar> vkl_e(5); //derivative wrt e of a_e a_j + a_j a_e
  Array1 <doublevar> vkl_j(5); //wrt j
  
  int s;
  const doublevar tiny=1e-14;
  //cout << "updateLap " << endl;
  for(int at=0; at < natoms; at++) {
    int p=parm_centers(at);
    for(int i=0; i< _nparms(p); i++) {
      doublevar parm;
      int k=klm(i,0), el=klm(i,1), m=klm(i,2);
      if(fabs(eibasis(e,at,k,0)) > tiny
	 || fabs(eibasis(e,at,el,0)) > tiny) { 
	for(int j=0; j< e; j++) {
	  s=(j < nspin_up) != (e < nspin_up);
	  parm=unique_parameters_spin(s)(p,i);
	  doublevar vkl=parm*(eibasis(e,at,k,0)*eibasis(j,at,el,0)
			      //+eibasis(j,at,k,0)*eibasis(e,at,el,0)
			      );
	  
	  
	  for(int d=1; d< 5; d++) {
	    vkl_e[d]=parm*(eibasis(e,at,k,d)*eibasis(j,at,el,0)
			   //+eibasis(j,at,k,0)*eibasis(e,at,el,d)
			   );
	    vkl_j[d]=parm*(eibasis(e,at,k,0)*eibasis(j,at,el,d)
			   //+eibasis(j,at,k,d)*eibasis(e,at,el,0)
			   );
	  } 
	  lap(e,at,0)+=vkl*eebasis(j,m,0);
	  doublevar dot_e=0, dot_j=0;
	  for(int d=1; d< 4; d++) {
	    lap(e,at,d)+=vkl_e[d]*eebasis(j,m,0);
	    lap(e,at,d)-=vkl*eebasis(j,m,d);
	    lap(j,at,d)+=vkl_j[d]*eebasis(j,m,0);
	    lap(j,at,d)+=vkl*eebasis(j,m,d);
	    
	    dot_e+=vkl_e[d]*eebasis(j,m,d);
	    dot_j+=vkl_j[d]*eebasis(j,m,d);
	  }
	  
	  
	  lap(e,at,4)+=vkl*eebasis(j,m,4);
	  lap(j,at,4)+=vkl*eebasis(j,m,4);
	  
	  lap(e,at,4)+=vkl_e[4]*eebasis(j,m,0);
	  lap(j,at,4)+=vkl_j[4]*eebasis(j,m,0);
	  
	  lap(e,at,4)-=2*dot_e;
	  lap(j,at,4)+=2*dot_j;
	}
	for(int j=e+1; j< nelectrons; j++) {
	  s=(j < nspin_up) != (e < nspin_up);
	  parm=unique_parameters_spin(s)(p,i);
	  doublevar vkl=parm*(eibasis(e,at,k,0)*eibasis(j,at,el,0)
			      //+eibasis(j,at,k,0)*eibasis(e,at,el,0)
			      );
	  for(int d=1; d< 5; d++) {
	    vkl_e[d]=parm*(eibasis(e,at,k,d)*eibasis(j,at,el,0)
			   //+eibasis(j,at,k,0)*eibasis(e,at,el,d)
			   );
	    vkl_j[d]=parm*(eibasis(e,at,k,0)*eibasis(j,at,el,d)
			   //+eibasis(j,at,k,d)*eibasis(e,at,el,0)
			   );
	  } 
	  
	  lap(e,at,0)+=vkl*eebasis(j,m,0);
	  doublevar dot_e=0, dot_j=0;
	  for(int d=1; d< 4; d++) {
	    lap(e,at,d)+=vkl_e[d]*eebasis(j,m,0);
	    lap(e,at,d)+=vkl*eebasis(j,m,d);
	    lap(j,at,d)+=vkl_j[d]*eebasis(j,m,0);
	    lap(j,at,d)-=vkl*eebasis(j,m,d);
	    
	    dot_e+=vkl_e[d]*eebasis(j,m,d);
	    dot_j+=vkl_j[d]*eebasis(j,m,d);
	  }
	  lap(e,at,4)+=vkl*eebasis(j,m,4);
	  lap(j,at,4)+=vkl*eebasis(j,m,4);
	  
	  lap(e,at,4)+=vkl_e[4]*eebasis(j,m,0);
	  lap(j,at,4)+=vkl_j[4]*eebasis(j,m,0);
	  
	  lap(e,at,4)+=2*dot_e;
	  lap(j,at,4)-=2*dot_j;
	}
      }
    }
  }
  
  //cout<<"Done"<<endl;
}
//end: added for ei back-flow

//-----------------------------------------------------------

void Jastrow_threebody_piece_diffspin::getParmDeriv(const Array3 <doublevar> & eibasis,
                                        const Array3 <doublevar> & eebasis,
                                       Parm_deriv_return & deriv) {
  //cout << "Jastrow_threebody_piece_diffspin::getParmDeriv" << endl;
  assert(eibasis.GetDim(1) >= parm_centers.GetDim(0));
  int natoms=parm_centers.GetDim(0);
  int nelectrons=eebasis.GetDim(0);
  
  const doublevar tiny=1e-14;
  //cout << "updateLap " << endl;
  for(int at=0; at < natoms; at++) {
    int p=parm_centers(at);
    for(int i=0; i< _nparms(p); i++) {
      int k=klm(i,0), el=klm(i,1), m=klm(i,2);
      for(int e=0; e< nelectrons; e++) { 
        if(fabs(eibasis(e,at,k)) > tiny
           || fabs(eibasis(e,at,el)) > tiny) { 
          for(int j=e+1; j< nelectrons; j++) {
            int s=(j < nspin_up) != (e < nspin_up);
            int index=linear_parms(s)(p,i);
            doublevar vkl=(//eibasis(e,at,k)*eibasis(j,at,el)
			   +eibasis(j,at,k)*eibasis(e,at,el));
            deriv.gradient(index)+=vkl*eebasis(e,j,m);
          }
        }
      }
    }
  } 
}

//-----------------------------------------------------------

int Jastrow_threebody_piece_diffspin::nparms() {
  int tot=0;
  int nsec;
  if(freeze) return 0;
  for(int s=0; s< unique_parameters_spin.GetSize();s++){
    nsec=unique_parameters_spin(s).GetDim(0);
    assert(nsec==_nparms.GetDim(0));
    for(int i=0; i<nsec; i++)
      tot+=_nparms(i);
  }
  return tot;
}
//-----------------------------------------------------------

void Jastrow_threebody_piece_diffspin::getParms(Array1 <doublevar> & parms) {
  if(freeze) {
    parms.Resize(0);
  }
  else {
    parms.Resize(nparms());
    int counter=0;
    int natoms=parm_centers.GetDim(0);
    for(int s=0; s< unique_parameters_spin.GetSize();s++){
      for(int i=0; i< unique_parameters_spin(s).GetDim(0); i++) {
	for(int j=0; j< _nparms(i); j++) {
	  parms(counter++) = unique_parameters_spin(s)(i,j);//*natoms;
	}
      }
    }
  }

}
//-----------------------------------------------------------

void Jastrow_threebody_piece_diffspin::setParms(Array1 <doublevar> & parms) {
  assert(parms.GetDim(0)==nparms());

  if(freeze) return;
  int counter=0;
  int natoms=parm_centers.GetDim(0);
  for(int s=0; s< unique_parameters_spin.GetSize();s++){
    for(int i=0; i< unique_parameters_spin(s).GetDim(0); i++) {
      for(int j=0; j< _nparms(i); j++) {
	unique_parameters_spin(s)(i,j)=parms(counter++);//natoms;
	//cout << "set one-body " << unique_parameters(i,j) << endl;
      }
    }
  }
}

//--------------------------------------------------------------------------


