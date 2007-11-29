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



//----------------------------------------------------------------------

void Jastrow_twobody_piece::set_up(vector <string> & words, System * sys) {
  vector <string> coefftxt;
  unsigned int pos;
  //cout << "readsection " << endl;
  if(haskeyword(words, pos=0, "FREEZE"))
    freeze=1;
  else freeze=0;


  if(!readsection(words, pos=0, coefftxt, "COEFFICIENTS"))
    error("need COEFFICIENTS in jastrow twobody");
  int np=coefftxt.size();
  parameters.Resize(np);
  //cout << np << " two-body  parameters " << endl;
  for(int p=0; p< np; p++) {
    parameters(p)=atof(coefftxt[p].c_str());
  }
  nelectrons=sys->nelectrons(0)+sys->nelectrons(1);
  
}

//----------------------------------------------------------------------

int Jastrow_twobody_piece::writeinput(string & indent, ostream & os) {
  if(freeze)
    os << indent << "FREEZE" << endl;

  os << indent << "COEFFICIENTS { ";
  for(int i=0; i< parameters.GetDim(0); i++) {
    os << parameters(i) << "  ";
  }
  os << " } " << endl;

  return 1;
}

//----------------------------------------------------------------------

int Jastrow_twobody_piece::showinfo(string & indent, ostream & os) {
  if(freeze)
    os << indent << "Coefficients frozen" << endl;

  os << indent << "Coefficients ";
  for(int i=0; i< parameters.GetDim(0); i++) {
    os << parameters(i) << "  ";
  }
  os <<  endl;

  return 1;
}

//----------------------------------------------------------------------

void Jastrow_twobody_piece::updateLap(int e, const Array3 <doublevar> & eebasis,
                                      Array2 <doublevar> & lap) {

  //cout << "Jastrow_twobody_piece::updateLap" << endl;
  assert(lap.GetDim(0) == eebasis.GetDim(0));
  assert(parameters.GetDim(0)<= eebasis.GetDim(1));
  int nelectrons=eebasis.GetDim(0);
  int np=parameters.GetDim(0);

  //lap=0;

  for(int i=0; i < e; i++) {
    for(int p=0; p < np; p++) {
      for(int d=0; d< 5; d++)
        //lap(i,d)+=parameters(p)*eebasis(i,p+1,d);
	lap(i,d)+=parameters(p)*eebasis(i,p,d);
    }
  }

  for(int j=e+1; j < nelectrons; j++) {
    for(int p=0; p < np; p++) {
      for(int d=0; d< 5; d++)
        lap(j,d)+=parameters(p)*eebasis(j,p,d);
    }
  }
  //cout << "done" << endl;
}

//----------------------------------------------------------------------

void Jastrow_twobody_piece::updateVal(int e, const Array3 <doublevar> & eebasis,
                                      Array1 <doublevar> & val) {
  assert(val.GetDim(0) == eebasis.GetDim(0));
  int nelectrons=eebasis.GetDim(0);
  int np=parameters.GetDim(0);

  //val=0;

  for(int i=0; i < e; i++) {
    for(int p=0; p < np; p++) {
      val(i)+=parameters(p)*eebasis(i,p,0);
    }
  }

  for(int j=e+1; j < nelectrons; j++) {
    for(int p=0; p < np; p++) {
      val(j)+=parameters(p)*eebasis(j,p,0);
    }
  }

}
//----------------------------------------------------------------------


void Jastrow_twobody_piece::getParmDeriv(const Array3 <doublevar> & eebasis,
                                         Parm_deriv_return & parm_ret) { 
  assert(parm_ret.gradient.GetDim(0)>=nparms());
  int nelectrons=eebasis.GetDim(0);
  int np=parameters.GetDim(0);
  for(int i=0; i< nelectrons; i++) {
    for(int j=i+1; j < nelectrons; j++) { 
      for(int p=0; p < np; p++) { 
        parm_ret.gradient(p)+=eebasis(i,j,p);
      }
    }
  }
}
//----------------------------------------------------------------------

void Jastrow_twobody_piece::getParms(Array1 <doublevar> & parms) {
  if(freeze) parms.Resize(0);
  else {
    parms.Resize(parameters.GetDim(0));
    int nparms=parameters.GetDim(0);
    for(int p=0; p< nparms; p++)
      parms(p)=parameters(p);//*nelectrons;
  }
}

void Jastrow_twobody_piece::setParms(Array1 <doublevar> & parms) {
  if(freeze) { assert(parms.GetDim(0)==0); }
  else {
    assert(parms.GetDim(0)==parameters.GetDim(0));
    int nparms=parameters.GetDim(0);
    for(int p=0; p< nparms; p++) {
      parameters(p)=parms(p);///nelectrons;
      //cout << "set two-body " << parameters(p) << endl;
    }
    
  }
}



//----------------------------------------------------------------------
void Jastrow_twobody_piece::parameterSaveVal(int e,
                        const Array3 <doublevar> & eebasis,
                        Array1 <doublevar> & save, int begin_fill){}
//----------------------------------------------------------------------

void Jastrow_twobody_piece::parameterSaveLap(const Array4 <doublevar> & eebasis,
                        Array2 <doublevar> & save_lap, int begin_fill) {}

//----------------------------------------------------------------------
//######################################################################
void Jastrow_twobody_piece_diffspin::set_up(vector <string> & words, System * sys) {
  vector <string> coeff_like;
  vector <string> coeff_unlike;
  unsigned int pos;
  //cout << "readsection " << endl;
  if(haskeyword(words, pos=0, "FREEZE"))
    freeze=1;
  else freeze=0;

  if(!readsection(words, pos=0, coeff_like, "LIKE_COEFFICIENTS"))
    error("need LIKE_COEFFICIENTS in jastrow twobody");
  if(!readsection(words, pos=0, coeff_unlike, "UNLIKE_COEFFICIENTS"))
    error("need UNLIKE_COEFFICIENTS in jastrow twobody");

  if(coeff_like.size() != coeff_unlike.size() )
    error("in Jastrow twobody, must have the same number of like and unlike"
          "coefficients.");

  int np=coeff_like.size();
  spin_parms.Resize(2,np);
  //cout << 2*np << " two-body  parameters " << endl;
  for(int p=0; p< np; p++) {
    spin_parms(0,p)=atof(coeff_like[p].c_str());
    spin_parms(1,p)=atof(coeff_unlike[p].c_str());
  }

  nspin_up=sys->nelectrons(0);

}
//----------------------------------------------------------------------

int Jastrow_twobody_piece_diffspin::writeinput(string & indent, ostream & os) {

  if(freeze)
    os << indent << "FREEZE" << endl;

  os << indent << "LIKE_COEFFICIENTS { ";
  for(int i=0; i< spin_parms.GetDim(1); i++) {
    os << spin_parms(0,i) << "  ";
  }

  os << " } " <<endl;
  os << indent << "UNLIKE_COEFFICIENTS { ";
  for(int i=0; i< spin_parms.GetDim(1); i++) {
    os << spin_parms(1,i) << "  ";
  }
  os << " }  " << endl;
  return 1;
}

//----------------------------------------------------------------------

int Jastrow_twobody_piece_diffspin::showinfo(string & indent, ostream & os) {

  if(freeze)
    os << indent << "Coefficients frozen" << endl;

  os << indent << "Like coefficients ";
  for(int i=0; i< spin_parms.GetDim(1); i++) {
    os << spin_parms(0,i) << "  ";
  }
  os <<endl;

  os << indent << "Unlike Coefficients ";
  for(int i=0; i< spin_parms.GetDim(1); i++) {
    os << spin_parms(1,i) << "  ";
  }
  os << endl;
  return 1;
}


//----------------------------------------------------------------------

void Jastrow_twobody_piece_diffspin::updateLap(int e, const Array3 <doublevar> & eebasis,
                                      Array2 <doublevar> & lap) {

  //cout << "Jastrow_twobody_piece::updateLap" << endl;
  assert(lap.GetDim(0) == eebasis.GetDim(0));
  assert(spin_parms.GetDim(1) <= eebasis.GetDim(1));
  int nelectrons=eebasis.GetDim(0);
  int np=spin_parms.GetDim(1);

  //lap=0;

  int s;

  for(int i=0; i < e; i++) {
    s=(i < nspin_up) != (e<  nspin_up); //0 if they're the same, 1 otherwise
    for(int p=0; p < np; p++) {
      for(int d=0; d< 5; d++)
        lap(i,d)+=spin_parms(s,p)*eebasis(i,p,d);
    }
  }

  for(int j=e+1; j < nelectrons; j++) {
    s=(j< nspin_up) != (e < nspin_up);
    for(int p=0; p < np; p++) {
      for(int d=0; d< 5; d++)
        lap(j,d)+=spin_parms(s,p)*eebasis(j,p,d);
    }
  }
  //cout << "done" << endl;
}

//----------------------------------------------------------------------

void Jastrow_twobody_piece_diffspin::updateVal(int e, const Array3 <doublevar> & eebasis,
                                      Array1 <doublevar> & val) {
  assert(val.GetDim(0) == eebasis.GetDim(0));
  int nelectrons=eebasis.GetDim(0);
  int np=spin_parms.GetDim(1);

  int s;
  for(int i=0; i < e; i++) {
    s=(i < nspin_up) != (e<  nspin_up); //0 if they're the same, 1 otherwise
    for(int p=0; p < np; p++) {
      val(i)+=spin_parms(s,p)*eebasis(i,p,0);
    }
  }

  for(int j=e+1; j < nelectrons; j++) {
    s=(j< nspin_up) != (e < nspin_up);
    for(int p=0; p < np; p++) {
      val(j)+=spin_parms(s,p)*eebasis(j,p,0);
    }
  }

}

//----------------------------------------------------------------------

void Jastrow_twobody_piece_diffspin::getParms(Array1 <doublevar> & parms) {
  if(freeze) parms.Resize(0);
  else {
    parms.Resize(2*spin_parms.GetDim(1));
    int counter=0;
    for(int i=0; i< 2; i++) {
      for(int j=0; j< spin_parms.GetDim(1); j++)
        parms(counter++)=spin_parms(i,j);
    }
  }
}

void Jastrow_twobody_piece_diffspin::setParms(Array1 <doublevar> & parms) {
  if(freeze) assert(parms.GetDim(0)==0);
  else {
    assert(parms.GetDim(0)==2*spin_parms.GetDim(1));

    int counter=0;
    for(int i=0; i< 2; i++) {
      for(int j=0; j< spin_parms.GetDim(1); j++)
        spin_parms(i,j)=parms(counter++);
    }
  }


}

//----------------------------------------------------------------------

void Jastrow_twobody_piece_diffspin::getParmDeriv(const Array3 <doublevar> & eebasis,
                                                  Parm_deriv_return & parm_ret) { 
  assert(parm_ret.gradient.GetDim(0)>=2*spin_parms.GetDim(1));
  Array2 <doublevar> grad_spin_parms(2,spin_parms.GetDim(1));
  grad_spin_parms=0;

  int nelectrons=eebasis.GetDim(0);
  int np=spin_parms.GetDim(1);
  for(int i=0; i< nelectrons; i++) {
    for(int j=i+1; j < nelectrons; j++) { 
      int s=(j< nspin_up) != (i < nspin_up);
      for(int p=0; p < np; p++) { 
        grad_spin_parms(s,p)+=eebasis(i,j,p);
      }
    }
  }
  
  int counter=0;
  for(int i=0; i< 2; i++) {
    for(int j=0; j< spin_parms.GetDim(1); j++)
      parm_ret.gradient(counter++)=grad_spin_parms(i,j);
  }

}
//----------------------------------------------------------------------



//----------------------------------------------------------------------
void Jastrow_twobody_piece_diffspin::parameterSaveVal(int e,
                        const Array3 <doublevar> & eebasis,
                        Array1 <doublevar> & save, int begin_fill){}
//----------------------------------------------------------------------

void Jastrow_twobody_piece_diffspin::parameterSaveLap(const Array4 <doublevar> & eebasis,
                        Array2 <doublevar> & save_lap, int begin_fill) {}




//--------------------------------------------------------------------------
