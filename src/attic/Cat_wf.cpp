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

//--------------------------------------------------------------------------
// Cat_wf.cpp
//
//
#include "Cat_wf.h"
#include "qmc_io.h"
#include "Qmc_std.h"
#include "Sample_point.h"


void Cat_wf::notify(change_type change , int num )
{
  for(int i=0; i< nUniqueFunc; i++)
  {
    wfarray(i)->notify(change, num);
  }
}



void Cat_wf::storeParmIndVal(Wavefunction_data * wfdata, Sample_point * sample,
                             int e, Array1 <doublevar> & vals )
{
  Cat_wf_data * dataptr;
  recast(wfdata, dataptr);
  assert(dataptr != NULL);
  Array1 <Array1 <doublevar> > temparray(nUniqueFunc);

  for(int i=0; i< nUniqueFunc; i++)
  {
    temparray(i).Resize(dataptr->data_array(i)->valSize());
    wfarray(i)->storeParmIndVal(dataptr->data_array(i),sample,  e,temparray(i));
  }
  int counter=0;
  for(int i=0; i<nUniqueFunc; i++)
  {
    for(int j=0; j< temparray(i).GetDim(0); j++)
    {
      vals(counter)=temparray(i)(j);
      counter++;
    }
  }
}



void Cat_wf::getParmDepVal(Wavefunction_data * wfdata,
                           Sample_point * sample,
                           int e,
                           Array1 <doublevar> & oldval,
                           Wf_return & newval)
{
  Cat_wf_data * dataptr;
  recast(wfdata, dataptr);
  assert(dataptr != NULL);
  assert(oldval.GetDim(0) >= dataptr->valSize());
  Array1 < Array1 <doublevar> > temparray(nUniqueFunc);
  Array1 <Wf_return > tempval(nUniqueFunc);
  int counter=0;

  for(int i=0; i< nUniqueFunc; i++)
  {
    tempval(i).Resize(wfarray(i)->nfunc(),2);
    temparray(i).Resize(dataptr->data_array(i)->valSize());
    for(int j=0; j< temparray(i).GetDim(0); j++)
    {
      temparray(i)(j)=oldval(counter);
      counter++;
    }
    wfarray(i)->getParmDepVal(dataptr->data_array(i), sample, e,
                              temparray(i),tempval(i));
  }

  counter=0;
  for(int i=0; i< ordering.GetDim(0); i++)
  {
    int index=ordering(i);
    for(int j=0; j< tempval(index).amp.GetDim(0); j++)
    {
      for(int k=0; k<2; k++)
      {
        newval.amp(counter,k)=tempval(index).amp(j,k);
        newval.phase(counter,k)=tempval(index).phase(j,k);
      }
      counter++;
      assert(counter <= nFunc);
    }
  }

  counter=0;
  for(int i=0; i<nUniqueFunc; i++)
  {
    for(int j=0; j< temparray(i).GetDim(0); j++)
    {
      oldval(counter)=temparray(i)(j);
      counter++;
    }
  }

}




void Cat_wf::generateStorage(Wavefunction_storage * & wfstore)
{
  wfstore=new Cat_wf_storage;
  Cat_wf_storage * store;
  recast(wfstore, store);
  store->store_array.Resize(nUniqueFunc);
  for(int i=0; i< nUniqueFunc; i++)
  {
    wfarray(i)->generateStorage(store->store_array(i));
  }
}

void Cat_wf::saveUpdate(Sample_point * sample, int e, 
                        Wavefunction_storage * wfstore)
{
  //cout << "saveUpdate\n";
  Cat_wf_storage * store;
  recast(wfstore, store);

  for(int i=0; i< nUniqueFunc; i++)
  {
    wfarray(i)->saveUpdate(sample, e, store->store_array(i));
  }
}

void Cat_wf::restoreUpdate(Sample_point * sample, int e, 
                           Wavefunction_storage * wfstore)
{
  //cout << "restoreUpdate\n";
  Cat_wf_storage * store;
  recast(wfstore, store);

  for(int i=0; i< nUniqueFunc; i++)
  {
    wfarray(i)->restoreUpdate(sample, e, store->store_array(i));
  }

}

void Cat_wf::getDensity(Wavefunction_data * wfdata, int e,
                        Array2 <doublevar> & dens)
{
  assert(dens.GetDim(0) >= nFunc);
  Cat_wf_data * dataptr;
  recast(wfdata, dataptr);
  Array1 <Array2 <doublevar> > tempval(nUniqueFunc);
  for(int i=0; i< nUniqueFunc; i++)
  {
    tempval(i).Resize(wfarray(i)->nfunc(),2);
    wfarray(i)->getDensity(dataptr->data_array(i),e, tempval(i));
  }
  int counter=0;
  for(int i=0; i< ordering.GetDim(0); i++)
  {
    int index=ordering(i);
    for(int j=0; j< tempval(index).GetDim(0); j++)
    {
      for(int d=0; d< 2; d++)
      {
        dens(counter,d)=tempval(index)(j,d);
      }
      counter++;
      assert(counter <= nFunc);
    }
  }

}


void Cat_wf::getVal(Wavefunction_data * wfdata, int e,
                    Wf_return & val, int startpos)
{
  assert(val.amp.GetDim(0) >= nFunc+startpos);
  Array1 <Wf_return > tempval(nUniqueFunc);
  Cat_wf_data * dataptr;
  recast(wfdata, dataptr);
  assert(dataptr != NULL);

  for(int i=0; i< nUniqueFunc; i++)
  {
    tempval(i).Resize(wfarray(i)->nfunc(),2);
    wfarray(i)->getVal(dataptr->data_array(i), e, tempval(i));
  }
  int counter=startpos;
  for(int i=0; i< ordering.GetDim(0); i++)
  {
    int index=ordering(i);
    for(int j=0; j< tempval(index).amp.GetDim(0); j++)
    {
      for(int d=0; d< 2; d++)
      {
        val.amp(counter,d)=tempval(index).amp(j,d);
        val.phase(counter,d)=tempval(index).phase(j,d);
      }
      counter++;
      assert(counter <= nFunc+startpos);
    }
  }

}

void Cat_wf::getForceBias(Wavefunction_data * wfdata, int e,
                          Wf_return & bias, int startpos)
{


  assert(bias.amp.GetDim(0) >= nFunc+startpos);
  assert(bias.amp.GetDim(1) >= 5);
  Cat_wf_data * dataptr;
  recast(wfdata, dataptr);

  assert(dataptr != NULL);

  Array1 <Wf_return > tempbias(nUniqueFunc);

  for(int i=0; i< nUniqueFunc; i++)
  {
    tempbias(i).Resize(wfarray(i)->nfunc(), 5);
    wfarray(i)->getForceBias(dataptr->data_array(i), e, tempbias(i));
  }

  int counter=startpos;
  for(int i=0; i< ordering.GetDim(0); i++)
  {
    int index=ordering(i);
    for(int j=0; j< tempbias(index).amp.GetDim(0); j++)
    {
      for(int d=0; d< 5; d++)
      {
        bias.amp(counter,d)=tempbias(index).amp(j,d);
        bias.phase(counter,d)=tempbias(index).phase(j,d);
      }
      counter++;
      assert(counter <= nFunc);
    }
  }


}


void Cat_wf::getLap(Wavefunction_data * wfdata, int e,
                    Wf_return & lap, int startpos)
{

  assert(lap.amp.GetDim(0) >= nFunc+startpos);
  assert(lap.amp.GetDim(1) >= 6);
  Cat_wf_data * dataptr;
  recast(wfdata, dataptr);
  assert(dataptr != NULL);

  Array1 <Wf_return> templap(nUniqueFunc);

  for(int i=0; i< nUniqueFunc; i++)
  {
    templap(i).Resize(wfarray(i)->nfunc(), 6);
    wfarray(i)->getLap(dataptr->data_array(i), e, templap(i));
  }



  int counter=startpos;
  for(int i=0; i< ordering.GetDim(0); i++)
  {
    int index=ordering(i);
    for(int j=0; j< templap(index).amp.GetDim(0); j++)
    {
      for(int d=0; d< 6; d++)
      {

        lap.amp(counter,d)=templap(index).amp(j,d);
        lap.phase(counter,d)=templap(index).phase(j,d);

      }
      counter++;
      assert(counter <= nFunc);
    }
  }

}

void Cat_wf::updateVal(Wavefunction_data * wfdata, Sample_point * sample)
{
  Cat_wf_data * dataptr;
  recast(wfdata, dataptr);
  assert(dataptr != NULL);
  for(int i=0; i< nUniqueFunc; i++)
  {
    wfarray(i)->updateVal(dataptr->data_array(i), sample);
  }
}

void Cat_wf::updateLap(Wavefunction_data * wfdata, Sample_point * sample)
{
  //cout << "updateLap " << endl;
  Cat_wf_data * dataptr;
  recast(wfdata, dataptr);

  assert(dataptr != NULL);
  for(int i=0; i< nUniqueFunc; i++)
  {
    wfarray(i)->updateLap(dataptr->data_array(i), sample);
  }
  //cout << "done Lap" << endl;
}

void Cat_wf::updateForceBias(Wavefunction_data * wfdata, Sample_point * sample)
{
  Cat_wf_data * dataptr;
  recast(wfdata, dataptr);
  assert(dataptr != NULL);
  for(int i=0; i<nUniqueFunc; i++)
  {
    wfarray(i)->updateForceBias(dataptr->data_array(i), sample);
  }
}
