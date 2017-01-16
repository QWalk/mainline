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
//------------------------------------------------------------------------
//src/Ext_slat_data.cpp

#include "Qmc_std.h"
#include "qmc_io.h"
#include "Ext_slat_data.h"
#include "Wavefunction_data.h"
#include "Ext_slat.h"
#include <algorithm>



/*!
*/
void Ext_slat_data::read(vector <string> & words, unsigned int & pos,
                        System * sys)
{
  vector <string> motxt;
  if(!readsection(words, pos=0, motxt, "ORBITALS") ) 
    error("Need section ORBITALS in EXT_SLAT");

  vector <string> pairtxt;
  if(!readsection(words, pos=0, pairtxt, "PAIRS"))
    error("Need section PAIRS in EXT_SLAT");
  if(pairtxt.size()%2 !=0)
    error("PAIRS needs a set of doubles ");
  
  vector <string> coefftxt;
  if(!readsection(words, pos=0, coefftxt, "COEFF")) 
    error("Need COEFF");

  vector <string> amptxt;
  if(!readsection(words, pos=0, amptxt, "AMPLITUDE")) 
    error("Need AMPLITUDE");

  if(coefftxt.size() != pairtxt.size()/2) 
    error("coefftxt != pairtxt/2");

  if(amptxt.size() != pairtxt.size()/2) 
    error("amptxt != pairtxt/2");



  allocate(motxt, sys, molecorb);
  
  int nup=sys->nelectrons(0);
  int ndown=sys->nelectrons(1);
  
  int nelectrons=nup+ndown;

  pairs.Resize(2,nelectrons, 2);
  pair_coeff.Resize(2,nelectrons, 2);
  amplitude.Resize(2,nelectrons);

  //------------------------------------------------------
  //This part is specific to total spin=0..

  if(nup != ndown) 
    error("open-shell not supported yet");

  cout << "nelectrons " << nelectrons 
       << " pair size " << pairtxt.size() << endl;
  if(pairtxt.size() != 2*nelectrons) 
    error("pairs must contain a double for every orbital");

  //a,b
  int bigmo=0;
  int count=0;
  for(int e=0; e< nelectrons/2; e++) {
    cout << "electron " << e << " count " << count << endl;
    pairs(0,e,0)=atoi(pairtxt[count++].c_str())-1;
    pairs(0,e,1)=atoi(pairtxt[count++].c_str())-1;

    if(pairs(0,e,0) > bigmo) bigmo=pairs(0,e,0);
    if(pairs(0,e,1) > bigmo) bigmo=pairs(0,e,1);

    pair_coeff(0,e,1)=atof(coefftxt[e].c_str());
    //pair_coeff(0,e,0)=sqrt(1-pair_coeff(0,e,1)*pair_coeff(0,e,1));
    pair_coeff(0,e,0)=1;
    amplitude(0,e)=atof(amptxt[e].c_str());
  }

  for(int e=nelectrons-1; e >= nelectrons/2; e--) {
    cout << "electron " << e << " count " << count << endl;
    pairs(1,e,0)=atoi(pairtxt[count++].c_str())-1;
    pairs(1,e,1)=atoi(pairtxt[count++].c_str())-1;

    if(pairs(1,e,0) > bigmo) bigmo=pairs(0,e,0);
    if(pairs(1,e,1) > bigmo) bigmo=pairs(0,e,1);

    pair_coeff(1,e,1)=atof(coefftxt[e].c_str());
    //pair_coeff(1,e,0)=sqrt(1-pair_coeff(1,e,1)*pair_coeff(1,e,1));
    pair_coeff(1,e,0)=1;
    amplitude(1,e)=atof(amptxt[e].c_str());
  }


  bigmo+=1;
  propagate_irreducible();
  nmo=bigmo;

  for(int s=0; s< 2; s++) {
    cout << "****For spin channel " << s << endl;
    for(int e=0; e< nelectrons; e++) {
      cout << "pair " << e << " : ";
      for(int p=0; p < 2; p++) {
        cout << pairs(s,e,p) << "  ";
      }
      cout << "   amplitude " << amplitude(s,e);
      cout << "  coeff ";
      for(int p=0; p < 2; p++) {
        cout << pair_coeff(s,e,p) << "  ";
      }      
      cout << endl;
    }
  }

  

  //----------------------------------------------------

  if(nmo > molecorb->getNmo() ) 
    error("Asked for more MO's than the MO object has.");

  Array1 < Array1 <int > > totoccupation(1);
  totoccupation(0).Resize(bigmo);
  for(int i=0; i< bigmo; i++) {
    totoccupation(0)(i)=i;
  }
  molecorb->buildLists(totoccupation);
  
  spin.Resize(nelectrons);
  for(int e=0; e < nup; e++) 
    spin(e)=0;
  for(int e=nup; e< nelectrons; e++)
    spin(e)=1;


}


//----------------------------------------------------------------------

void Ext_slat_data::propagate_irreducible() {

  //BUG: we should pad with zeros if we have a spin-polarized case

  int nelectrons=pairs.GetDim(1);
  //Get the coefficients for the c,d,..

  //cout << "doing same spin" << endl;

  for(int e=nelectrons/2; e< nelectrons; e++) {
    int oppe=e-nelectrons/2;
    pairs(0,e,0)=pairs(0,oppe,0);
    pairs(0,e,1)=pairs(0,oppe,1);
    pair_coeff(0,e,0)=pair_coeff(0,oppe,1);
    pair_coeff(0,e,1)=-pair_coeff(0,oppe,0);
    amplitude(0,e)=1-amplitude(0,oppe);
  }
                                 
  //cout << "doing opposite spin " << endl;
  //Now for the opposite spin, we reverse the order
  for(int e=0; e < nelectrons/2; e++) {
    int oppe=e+nelectrons/2;
    
    pairs(1,e,0)=pairs(1,oppe,0);
    pairs(1,e,1)=pairs(1,oppe,1);
    pair_coeff(1,e,0)=pair_coeff(1,oppe,1);
    pair_coeff(1,e,1)=-pair_coeff(1,oppe,0);
    amplitude(1,e)=1-amplitude(1,oppe);
  }
}
  


//----------------------------------------------------------------------

int Ext_slat_data::supports(wf_support_type support) {
  switch(support) {
  case laplacian_update:
    return 1;
  case density:
    return 1;
  default:
    return 0;
  }
}

//----------------------------------------------------------------------

//----------------------------------------------------------------------

void Ext_slat_data::generateWavefunction(Wavefunction *& wf)
{
  assert(wf==NULL);
  wf=new Ext_slat;
  Ext_slat * slatwf;
  recast(wf, slatwf);
  slatwf->init(this);
  attachObserver(slatwf);
}

int Ext_slat_data::showinfo(ostream & os)
{
  os << "Extended Slater function" << endl;
  

  os << "Molecular Orbital object : ";
  molecorb->showinfo(os);

  return 1;
}

//----------------------------------------------------------------------

int Ext_slat_data::writeinput(string & indent, ostream & os)
{
  os << indent << "EXT_SLAT" << endl;

  os << indent << "PAIRS { " << endl;
  int npair=pairs.GetDim(1)/2;
  int nelectrons=npair*2;

  cout << "nelectrons " << nelectrons << endl;

  for(int i=0; i< npair; i++) {
    os << indent << "  " << pairs(0,i,0)+1 << "  " << pairs(0,i,1)+1 << endl;
  }

  for(int e=nelectrons-1; e >= nelectrons/2; e--) {
    os << indent << "   " << pairs(1,e,0)+1
       << "  " << pairs(1,e,1)+1 << endl;
  }


  os << indent << "}\n";

  os << indent << "COEFF { ";
  for(int i=0; i< npair; i++) {
    os << pair_coeff(0,i,1) << "  ";
  }
  os <<endl << indent << "     ";
  for(int e=nelectrons/2; e < nelectrons; e++) {
    os << pair_coeff(1,e,1) << "  ";
  }
  os << "}" <<  endl;

  os << indent << "AMPLITUDE { ";
  for(int i=0; i < npair; i++) {
    os << amplitude(0,i) << "  ";
  }
  os << endl << indent << "             ";
  for(int e=nelectrons-1; e >= nelectrons/2; e--) {
    os << amplitude(1,e) << "   ";
  }
  os << "}" <<  endl;
  

  

  os << indent << "ORBITALS {\n";
  string indent2=indent+"  ";
  molecorb->writeinput(indent2, os);
  os << indent << "}\n";

  return 1;
}

//------------------------------------------------------------------------
void Ext_slat_data::getVarParms(Array1 <doublevar> & parms)
{
  parms.Resize(nparms());
  int npair=amplitude.GetDim(1)/2;
  int counter=0;
  int nelectrons=amplitude.GetDim(1);
  for(int s=0; s< 2; s++) {
    for(int i=0; i< npair; i++) {
      int f=i+s*nelectrons/2;
      //parms(counter++)=tan((pi-.2)*(amplitude(s,f)-.5));
      //cout << "amplitude " << amplitude(s,f) << " parm " << parms(counter-1) << endl;
      
    }
    
    for(int i=0; i< npair; i++) {
      int f=i+s*nelectrons/2;
      parms(counter++)=tan(.5*(pi-.2)*(pair_coeff(s,f,1)));
    }
  }  
}

void Ext_slat_data::setVarParms(Array1 <doublevar> & parms)
{
  assert(parms.GetDim(0)==nparms());

  int counter=0;
  int npair=amplitude.GetDim(1)/2;
  int nelectrons=amplitude.GetDim(1);
  for(int s=0; s< 2; s++) {
    for(int i=0; i< npair; i++) {
      int f=i+s*nelectrons/2;      
      //amplitude(s,f)=(atan(parms(counter++)))/(pi-.2)+.5;
    }
    
    for(int i=0; i< npair; i++) {
      int f=i+s*nelectrons/2;
      pair_coeff(s,f,1)=(atan(parms(counter++)))/(.5*(pi-.2)) ;
      //pair_coeff(s,f,0)=sqrt(1-pair_coeff(s,f,1)*pair_coeff(s,f,1));
    }
  }

  propagate_irreducible();
  
  for(int s=0; s< 2; s++) {
    cout << "****For spin channel " << s << endl;
    for(int e=0; e< nelectrons; e++) {
      cout << "pair " << e << " : ";
      for(int p=0; p < 2; p++) {
        cout << pairs(s,e,p) << "  ";
      }
      cout << "   amplitude " << amplitude(s,e);
      cout << "  coeff ";
      for(int p=0; p < 2; p++) {
        cout << pair_coeff(s,e,p) << "  ";
      }      
      cout << endl;
    }
  }



  int max=wfObserver.size();
  for(int i=0; i< max; i++) {
    wfObserver[i]->notify(all_wf_parms_change, 0);
  }
}
