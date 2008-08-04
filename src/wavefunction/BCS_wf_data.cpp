/*
 
Copyright (C) 2007 Lucas K. Wagner, 2008 Jindrich Kolorenc

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

#include "Qmc_std.h"
#include "qmc_io.h"
#include "Wavefunction_data.h"
#include "BCS_wf_data.h"
#include "BCS_wf.h"

/*!
*/
void BCS_wf_data::read(vector <string> & words, unsigned int & pos,
                        System * sys)
{

  /* 

  // In a polarized system we would need also one-body orbitals

  vector <string> mowords;
  if(!readsection(words, pos=0, mowords, "ORBITALS") ) 
    error("Need ORBITALS");

  allocate(mowords,sys,molecorb);

  int nmo=molecorb->getNmo();
  totoccupation(0).Resize(nmo);
  for(int i=0; i< nmo; i++) { 
    totoccupation(0)(i)=i;
  }

  occupation(0,0)=totoccupation(0);
  molecorb->buildLists(totoccupation);

  */

  if(!readvalue(words, pos=0, magnification_factor, "MAGNIFY")) {
    magnification_factor=1;
  }

  vector <string> correlation;
  if(!readsection(words, pos=0, correlation, "PAIR_ORBITAL"))
    error("Need PAIR_ORBITAL");
  jastdata.read(correlation,pos=0,sys);

  nelectrons.Resize(2);
  nelectrons(0)=sys->nelectrons(0);
  nelectrons(1)=sys->nelectrons(1);

  int tote=nelectrons(0)+nelectrons(1);
  spin.Resize(tote);
  rede.Resize(tote);
  int eup=nelectrons(0);
  for(int e=0; e<eup; e++) {
    spin(e)=0;
    rede(e)=e;
  }
  for(int e=eup; e<tote; e++) {
    spin(e)=1;
    rede(e)=e-eup;
  }
  


}

//----------------------------------------------------------------------

int BCS_wf_data::supports(wf_support_type support) {
  switch(support) {
  case laplacian_update:
    return 1;
  case density:
    return 0;
  case parameter_derivatives:
    return jastdata.supports(support);
    //return 0;
  default:
    return 0;
  }
}

//----------------------------------------------------------------------



void BCS_wf_data::generateWavefunction(Wavefunction *& wf)
{
  assert(wf==NULL);

  wf=new BCS_wf;
  BCS_wf * slatwf;
  recast(wf, slatwf);
  slatwf->init(this);
  attachObserver(slatwf);
}

int BCS_wf_data::showinfo(ostream & os)
{

  os << "BCS determinant with pair orbital of Jastrow form." << endl;
  os << "START of pair orbital definition..." << endl << endl;
  jastdata.showinfo(os);
  os << "...END of pair orbital definition." << endl << endl;

  return 1;
}

//----------------------------------------------------------------------

int BCS_wf_data::writeinput(string & indent, ostream & os)
{


  os << indent << "BCS" << endl;

  string indent2=indent+"  ";
  os << indent << "PAIR_ORBITAL { " << endl;
  jastdata.writeinput(indent2,os);
  os << indent << "}" << endl;
 
  /*
  os << indent << "ORBITALS { " << endl;
  molecorb->writeinput(indent2,os);
  os << indent << "}" << endl;
  */
  return 1;
}

//------------------------------------------------------------------------
void BCS_wf_data::getVarParms(Array1 <doublevar> & parms)
{

  parms.Resize(nparms());

  Array1 <doublevar> jastparms;
  jastdata.getVarParms(jastparms);

  int njparm=jastdata.nparms();
  for(int i=0; i< njparm; i++) parms(i)=jastparms(i);

}


void BCS_wf_data::setVarParms(Array1 <doublevar> & parms)
{

  int njparm=jastdata.nparms();
  Array1 <doublevar> jastparms(njparm);
  for(int i=0; i< njparm; i++) jastparms(i)=parms(i);
  jastdata.setVarParms(jastparms);

  int max=wfObserver.size();
  for(int i=0; i< max; i++) {
    wfObserver[i]->notify(all_wf_parms_change, 0);
  }

}

//----------------------------------------------------------------------
