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

void get_onebody_parms(vector<string> & words, vector<string> & atomnames,
                        Array2 <doublevar> & unique_parameters,
                        Array1 <int> & nparms,
                        vector <string> & parm_labels,
                        Array2 <int> & linear_parms,
                        Array1 <int> & parm_centers) {
  vector < vector < string> > parmtxt;
  vector <string> parmtmp;
  unsigned int pos=0;
  while(readsection(words, pos,parmtmp, "COEFFICIENTS"))
    parmtxt.push_back(parmtmp);

  if(parmtxt.size()==0) 
    error("Didn't find COEFFICIENTS in one of threebody or onebody Jastrow sections");

  int nsec=parmtxt.size();
  int maxsize=0;
  nparms.Resize(nsec);
  for(int i=0; i< nsec; i++)
    nparms(i)=parmtxt[i].size()-1;

  for(int i=0; i< nsec; i++)
    if(maxsize < nparms(i)) maxsize=nparms(i);

  unique_parameters.Resize(nsec, maxsize);
  linear_parms.Resize(nsec, maxsize);
  linear_parms=0; unique_parameters=0;
  int counter=0;
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

void Jastrow_onebody_piece::set_up(vector <string> & words, 
                                   vector <string> & atomnames) {
                                   
   get_onebody_parms(words, atomnames, unique_parameters,
                        _nparms, parm_labels,
                        linear_parms, parm_centers);
  unsigned int pos=0;

  if(haskeyword(words, pos=0, "FREEZE"))
    freeze=1;
  else freeze=0;
}
//--------------------------------------------------------------------------

int Jastrow_onebody_piece::writeinput(string & indent, ostream & os) {
  if(freeze) os << indent << "FREEZE" << endl;
  for(int i=0; i< unique_parameters.GetDim(0); i++) {
    os << indent << "COEFFICIENTS { " << parm_labels[i] << "  ";
    for(int j=0; j< _nparms(i); j++)
      os << unique_parameters(i,j) << "  ";
    os << " } " << endl;
  }
  return 1;

}

//--------------------------------------------------------------------------

int Jastrow_onebody_piece::showinfo(string & indent, ostream & os) {
  os << indent << "Atom   Coefficients " << endl;
  for(int i=0; i< unique_parameters.GetDim(0); i++) {
    os << indent << parm_labels[i] << "  ";
    for(int j=0; j< _nparms(i); j++)
      os << unique_parameters(i,j) << "  ";
    os << endl;
  }
  return 1;
  
}


//--------------------------------------------------------------------------


void Jastrow_onebody_piece::updateLap(int e,
                 const Array3 <doublevar> & eibasis,
                 Array1 <doublevar> & lap) {
  //cout << "Jastrow_onebody_piece::updateLap" << endl;
  assert(lap.GetDim(0) >= 5);
  assert(eibasis.GetDim(0) >= parm_centers.GetDim(0));
  int natoms=parm_centers.GetDim(0);
  for(int at=0; at < natoms; at++) {
    int p=parm_centers(at);
    for(int i=0; i< _nparms(p); i++) {
      for(int d=0; d< 5; d++) {
        //lap(d)+=unique_parameters(p,i)*eibasis(at,i+1,d);
	lap(d)+=unique_parameters(p,i)*eibasis(at,i,d);
      }
    }
  }
  //cout << "done " << endl;
}

//----------------------------------------------------------------------


void Jastrow_onebody_piece::updateLap_ion(int e,
                 const Array3 <doublevar> & eibasis,
                 Array3 <doublevar> & lap) {
  //cout << "Jastrow_onebody_piece::updateLap" << endl;
  assert(lap.GetDim(2) >= 5);
  assert(eibasis.GetDim(0) >= parm_centers.GetDim(0));
  int natoms=parm_centers.GetDim(0);
  for(int at=0; at < natoms; at++) {
    int p=parm_centers(at);
    for(int i=0; i< _nparms(p); i++) {
      for(int d=0; d< 5; d++) {
        lap(e,at,d)+=unique_parameters(p,i)*eibasis(at,i,d);
      }
    }
  }
  //cout << "done " << endl;
}



//-----------------------------------------------------------

void Jastrow_onebody_piece::updateVal(int e,
               const Array3 <doublevar> & eibasis,
               doublevar & val) {

  assert(eibasis.GetDim(0) >= parm_centers.GetDim(0));
  int natoms=parm_centers.GetDim(0);
  for(int at=0; at < natoms; at++) {
    int p=parm_centers(at);
    for(int i=0; i< _nparms(p); i++) {
      val+=unique_parameters(p,i)*eibasis(at,i,0);
    }
  }

}
//-----------------------------------------------------------

//will just add to the parm_deriv, since that'll make it easy and quick
//to loop through all the electrons
void Jastrow_onebody_piece::getParmDeriv(int e, const Array3 <doublevar> & eibasis,
                  Parm_deriv_return & parm_deriv) { 
  assert(parm_deriv.gradient.GetDim(0)==nparms() || freeze);
  assert(eibasis.GetDim(0) >= parm_centers.GetDim(0));
  if(freeze) return;
  int natoms=parm_centers.GetDim(0);
  for(int at=0; at < natoms; at++) {
    int p=parm_centers(at);
    for(int i=0; i< _nparms(p); i++) {
      int index=linear_parms(p,i);
      parm_deriv.gradient(index)+=eibasis(at,i,0);
    }
  }
}
//-----------------------------------------------------------

int Jastrow_onebody_piece::nparms() {
  int nsec=unique_parameters.GetDim(0);
  assert(nsec==_nparms.GetDim(0));

  if(freeze) return 0;
  int tot=0;
  for(int i=0; i<nsec; i++)
    tot+=_nparms(i);
  return tot;
}
//-----------------------------------------------------------

void Jastrow_onebody_piece::getParms(Array1 <doublevar> & parms) {
  if(freeze) {
    parms.Resize(0);
  }
  else {
    parms.Resize(nparms());
    int counter=0;
    int natoms=parm_centers.GetDim(0);
    for(int i=0; i< unique_parameters.GetDim(0); i++) {
      for(int j=0; j< _nparms(i); j++) {
        parms(counter++) = unique_parameters(i,j);//*natoms;
        
      }
    }
  }

}
//-----------------------------------------------------------

void Jastrow_onebody_piece::setParms(Array1 <doublevar> & parms) {
  assert(parms.GetDim(0)==nparms());

  if(freeze) return;
  int counter=0;
  int natoms=parm_centers.GetDim(0);
  for(int i=0; i< unique_parameters.GetDim(0); i++) {
    for(int j=0; j< _nparms(i); j++) {
      unique_parameters(i,j)=parms(counter++);///natoms;
      //cout << "set one-body " << unique_parameters(i,j) << endl;
    }
  }
}
//-----------------------------------------------------------


void Jastrow_onebody_piece::parameterSaveVal(int e,
                      const Array3 <doublevar> & eibasis,
                      Array1 <doublevar> & save, int begin_fill) {
  assert(eibasis.GetDim(0) >= parm_centers.GetDim(0));
  assert(save.GetDim(0) >= begin_fill+nparms());

  int natoms=parm_centers.GetDim(0);

  int totparm=nparms();
  for(int i=begin_fill; i< begin_fill+totparm; i++) {
    save(i)=0;
  }

  for(int at=0; at < natoms; at++) {
    int p=parm_centers(at);
    for(int n=0; n < _nparms(p); n++) {
      int index=begin_fill+linear_parms(p,n);
      save(index)+=eibasis(at, n,0);
    }
  }

}

//--------------------------------------------------------------------------
