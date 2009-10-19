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

#include "Qmc_std.h"
#include "qmc_io.h"
#include "Slat_Jastrow_data.h"
#include "Slat_Jastrow.h"

void Slat_Jastrow_data::attachObserver(Wavefunction * wf)
{
  Slat_Jastrow * newwf;
  recast(wf, newwf);
  slater->attachObserver(newwf->slater_wf);
  jastrow->attachObserver(newwf->jastrow_wf);
}

/*!
*/
void Slat_Jastrow_data::read(vector <string> & words, unsigned int & pos,
                             System * sys)
{
  vector <string> slater_section;
  vector <string> jastrow_section;
  unsigned int startpos=pos;
  if(!readsection(words, pos, slater_section, "WF1"))
    error("Need WF1 section in SLATER-JASTROW");
  pos=startpos;
  if(!readsection(words, pos, jastrow_section, "WF2"))
    error("Need WF2 section in SLATER-JASTROW");

  allocate(slater_section, sys, slater);
  allocate(jastrow_section, sys, jastrow);
}

//------------------------------------------------------------------------

int Slat_Jastrow_data::supports(wf_support_type support) {
  return slater->supports(support) && jastrow->supports(support);
}


//------------------------------------------------------------------------
void Slat_Jastrow_data::generateWavefunction(Wavefunction *& wf)
{
  assert(wf == NULL);

  wf=new Slat_Jastrow;
  Slat_Jastrow * sjwf;
  recast(wf, sjwf);
  sjwf->init(this);
  attachObserver(sjwf);
}

int Slat_Jastrow_data::valSize()
{
  return slater->valSize()+jastrow->valSize();
}

void Slat_Jastrow_data::getVarParms(Array1 <doublevar> & parms)
{
  Array1 <doublevar> wf1parms;
  Array1 <doublevar> wf2parms;

  slater->getVarParms(wf1parms);
  jastrow->getVarParms(wf2parms);

  parms.Resize(wf1parms.GetDim(0)+wf2parms.GetDim(0));

  for(int i=0; i< wf1parms.GetDim(0); i++)
  {
    parms(i)=wf1parms(i);
  }

  for(int i=wf1parms.GetDim(0); i< parms.GetDim(0); i++)
  {
    parms(i)=wf2parms(i-wf1parms.GetDim(0));
  }
}

void Slat_Jastrow_data::setVarParms(Array1 <doublevar> & parms)
{
  int nslatparms=slater->nparms();
  int njasparms=jastrow->nparms();
  
  assert(nslatparms+njasparms==parms.GetDim(0));

  Array1 <doublevar> slatparms(nslatparms);
  Array1 <doublevar> jasparms(njasparms);

  for(int i=0; i< nslatparms; i++)
  {
    slatparms(i)=parms(i);
  }

  for(int i=0; i< njasparms; i++)
  {
    jasparms(i)=parms(i+nslatparms);
  }

  slater->setVarParms(slatparms);
  jastrow->setVarParms(jasparms);
}

int Slat_Jastrow_data::showinfo(ostream & os)
{
  slater->showinfo(os);
  os << "\n";
  jastrow->showinfo(os);
  return 1;
}

int Slat_Jastrow_data::writeinput(string & indent, ostream & os)
{
  os << indent << "SLATER-JASTROW" << endl;
  string indent2=indent+"  ";
  os << indent << "WF1 { " << endl;
  slater->writeinput(indent2, os);
  os << indent << "}" << endl;

  os << indent << "WF2 { " << endl;
  jastrow->writeinput(indent2, os);
  os << indent << "}" << endl;
  return 1;
}
//------------------------------------------------------------------------
