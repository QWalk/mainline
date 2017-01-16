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


#include "ulec.h"
#include "SHO_sample.h"

#include "Wavefunction.h"




void SHO_sample::generateStorage(Sample_storage * & store) {
  //int nions=parent->ions.size();
  store=new Sample_storage;
  //store->iondist_temp.Resize(5,nions);
  store->pointdist_temp.Resize(5,nelectrons);
  store->pos_temp.Resize(3);
}

void SHO_sample::saveUpdate(int e, Sample_storage * store) {
  //cout << "saveUpdate" << endl;
  getElectronPos(e, store->pos_temp);
  //int nions=parent->ions.size();
  for(int d=0; d< 5; d++)
  {
    for(int i=0; i<e; i++)
    {
      store->pointdist_temp(d,i)=pointdist(d,i,e);
    }
    for(int j=e; j < nelectrons; j++)
    {
      store->pointdist_temp(d,j)=pointdist(d,e,j);
    }
  }
}

void SHO_sample::restoreUpdate(int e, Sample_storage * store) {
  //cout << "restoreUpdate" << endl;
  for(int i=0; i<3; i++)
  {
    elecpos(e,i)=store->pos_temp(i);
  }
  //int nions=parent->ions.size();
  for(int d=0; d< 5; d++)
  {
    for(int i=0; i<e; i++)
    {
      pointdist(d,i,e)=store->pointdist_temp(d,i);
    }
    for(int j=e; j < nelectrons; j++)
    {
      pointdist(d,e,j)=store->pointdist_temp(d,j);
    }
  }
}

void SHO_sample::updateEIDist(){}

void SHO_sample::updateEEDist(){}

void SHO_sample::init(System * sys) {
  assert(sys != NULL);
  recast(sys, parent); //assign sys to parent

  nd=parent->ndim();

  nelectrons=parent->nelectrons(0)+parent->nelectrons(1);

  elecpos.Resize(nelectrons, 3);
  elecpos=0.0;
  pointdist.Resize(5, nelectrons, nelectrons);
  pointdist=0.0;
  elecDistStale.Resize(nelectrons);

  elecDistStale=1;

}



void SHO_sample::randomGuess()
{
  //cout << "randomGuess " << endl;

  for(int e=0; e< nelectrons; e++) {
    for(int d=0; d< nd; d++) {
      elecpos(e,d)=rng.gasdev();
    }
    for(int d=nd; d< 3; d++) {
      elecpos(e,d)=0;
    }
  }

  //ionDistStale=1;
  elecDistStale=1;

  if(wfObserver) {
    wfObserver->notify(all_electrons_move, 0);
  }

}



void SHO_sample::setElectronPos(const int e,
                                    const Array1 <doublevar> & position)
{
  //cout << "setElectronPos" << endl;
  assert( position.GetDim(0) == 3 );

  for(int d=0; d< nd; d++) {
    elecpos(e,d)=position(d);
  }

  //ionDistStale(e)=1;
  elecDistStale(e)=1;


  if(wfObserver)
    wfObserver->notify(electron_move, e);
  

}


void SHO_sample::moveIon(const int ion, const Array1 <doublevar> & r)
{
  error("Don't support moveIon in SHO");
}

//----------------------------------------------------------------------


#include <iomanip>

void SHO_sample::rawOutput(ostream & os)
{
  os << "SHO_sample\n";

  os << "   numElectrons   " << nelectrons << endl;
  for(int e=0; e< nelectrons; e++) {
    for(int d=0; d< 3; d++) {
      os << setw(20) << setprecision(12) << elecpos(e,d);
    }
    os << "\n";
  }
}

//----------------------------------------------------------------------


void SHO_sample::rawInput(istream & is)
{
  string text;
  is >> text;
  if(text != "SHO_sample")
    error("expected SHO_sample, got ", text);

  is >> text;
  if(text != "numElectrons")
    error("expected numElectrons, got ", text);

  is >> nelectrons;
  for(int e=0; e< nelectrons; e++)
    for(int d=0; d< 3; d++)
      is >> elecpos(e,d);


}

//-------------------------------------------------------------------------
