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
// src/Ext_slat_calc.cpp
//
//
#include "Qmc_std.h"
#include "Ext_slat.h"
#include "MatrixAlgebra.h"
#include "Sample_point.h"
#include "Ext_slat_data.h"

//----------------------------------------------------------------------

/*!
*/
void Ext_slat::notify(change_type change, int num)
{
  switch(change)
  {
  case electron_move:
    electronIsStaleVal(num)=1;
    electronIsStaleLap(num)=1;
    break;
  case all_electrons_move:
    updateEverythingVal=1;
    updateEverythingLap=1;
    break;
  case wf_parm_change:  
  case all_wf_parms_change:
    updateEverythingVal=1;
    updateEverythingLap=1;
    break;
  case sample_attach:
    sampleAttached=1;
    updateEverythingVal=1;
    updateEverythingLap=1;
    break;
  case data_attach:
    dataAttached=1;
    updateEverythingVal=1;
    updateEverythingLap=1;
    break;
  case sample_static:
    staticSample=1;
    break;
  case sample_dynamic:
    staticSample=0;
    break;
  default:
    updateEverythingVal=1;
    updateEverythingLap=1;
  }
}

//----------------------------------------------------------------------

void Ext_slat::init(Wavefunction_data * wfdata)
{
  recast(wfdata, parent);

  int nmo=parent->nmo;


  nelectrons=parent->spin.GetDim(0);


  //Properties and intermediate calculation storage.


  moVal.Resize(nelectrons, nmo, 5);
  updatedMoVal.Resize(nmo,5);

  detVal=1;
  inverse.Resize(nelectrons, nelectrons);

  inverse=0;
  for(int e=0; e< nelectrons; e++) {
    inverse(e,e)=1;
  }


  electronIsStaleVal.Resize(nelectrons);
  electronIsStaleLap.Resize(nelectrons);

  electronIsStaleVal=0;
  electronIsStaleLap=0;
  updateEverythingVal=1;
  updateEverythingLap=1;
  sampleAttached=0;
  dataAttached=0;
  staticSample=0;
}
//----------------------------------------------------------------------

void Ext_slat::generateStorage(Wavefunction_storage * & wfstore)
{
  wfstore=new Ext_slat_storage;
  Ext_slat_storage * store;
  recast(wfstore, store);
  nelectrons=parent->spin.GetDim(0);
  store->moVal_temp.Resize (parent->nmo,5);
  store->inverse_temp.Resize(nelectrons, nelectrons);

}


//----------------------------------------------------------------------

void Ext_slat::saveUpdate(Sample_point * sample, int e,
                         Wavefunction_storage * wfstore)
{
  Ext_slat_storage * store;
  recast(wfstore, store);
  for(int m=0; m < parent->nmo; m++) {
    for(int d=0; d< 5; d++) {
      store->moVal_temp(m,d)=moVal(e,m,d);
    }
  }

  array_cp(store->inverse_temp, inverse);
  store->detVal=detVal;
}

//----------------------------------------------------------------------

void Ext_slat::restoreUpdate(Sample_point * sample, int e,
                            Wavefunction_storage * wfstore)
{
  Ext_slat_storage * store;
  recast(wfstore, store);
  for(int m=0; m < parent->nmo; m++) {
    for(int d=0; d< 5; d++) {
      moVal(e,m,d)=store->moVal_temp(m,d);
    }
  }

  array_cp(inverse, store->inverse_temp);
  detVal=store->detVal;

}

//----------------------------------------------------------------------

void Ext_slat::updateVal(Wavefunction_data * wfdata,
                        Sample_point * sample)
{
  assert(sampleAttached);
  assert(dataAttached);

  Ext_slat_data * slatdata;
  recast(wfdata, slatdata);

  if(updateEverythingVal==1) {
    calcVal(slatdata, sample);
    updateEverythingVal=0;
    electronIsStaleVal=0;
  }
  else
    {
      for(int e=0; e< nelectrons; e++)
      {
        if(electronIsStaleVal(e))
        {
          updateVal(slatdata, sample, e);
          electronIsStaleVal(e)=0;
        }
      }
    }
  
}

//----------------------------------------------------------------------

void Ext_slat::updateForceBias(Wavefunction_data * wfdata,
                              Sample_point * sample)
{
  //We just use the laplacian, since the speed diff is minimal.
  assert(sampleAttached);
  assert(dataAttached);

  updateLap(wfdata, sample);
}


//----------------------------------------------------------------------
/*!

 */
void Ext_slat::updateLap( Wavefunction_data * wfdata,
                        Sample_point * sample)
{
  assert(sampleAttached);
  assert(dataAttached);

  //cout << "updateLap\n";

  Ext_slat_data * slatdata;
  recast(wfdata, slatdata);
  if(updateEverythingLap==1) {
    calcLap(slatdata, sample);
    updateEverythingVal=0;
    updateEverythingLap=0;
    electronIsStaleLap=0;
    electronIsStaleVal=0;
  }
  else {
    for(int e=0; e< nelectrons; e++) {
      if(electronIsStaleLap(e)) {
          updateLap(slatdata, sample, e);
          electronIsStaleLap(e)=0;
          electronIsStaleVal(e)=0;
      }
    }
  }


}

//----------------------------------------------------------------------


void Ext_slat::storeParmIndVal(Wavefunction_data * wfdata, Sample_point * sample,
                              int e, Array1 <doublevar> & vals )
{
  
  assert(vals.GetDim(0)==0);
}

//----------------------------------------------------------------------

void Ext_slat::getParmDepVal(Wavefunction_data * wfdata,
                            Sample_point * sample,
                            int e,
                            Array1 <doublevar> & oldval,
                            Array2 <doublevar> & newval)
{

  updateVal(wfdata, sample);
  getVal(wfdata, e, newval);

}


//------------------------------------------------------------------------


void Ext_slat::calcVal(Ext_slat_data * dataptr, Sample_point * sample)
{

  calcLap(dataptr, sample);

}

//------------------------------------------------------------------------
/*!

*/
void Ext_slat::updateVal( Ext_slat_data * dataptr, Sample_point * sample,int e)
{
  updateLap(dataptr, sample, e);
}

//------------------------------------------------------------------------


void Ext_slat::getVal(Wavefunction_data * wfdata, int e,
                     Array2 <doublevar> & val, int startpos)
{
  val(startpos,0)=sign(detVal);
  val(startpos,1)=log(fabs(detVal));
}

//----------------------------------------------------------------------

void Ext_slat::getDensity(Wavefunction_data * wfdata, int e,
                         Array2 <doublevar> & dens)
{

  error("Ext_slat::No density");

}

//----------------------------------------------------------------------------


void Ext_slat::getForceBias(Wavefunction_data * wfdata, int e,
                           Array2 <doublevar> & bias, int startpos)
{
  //cout << "getForceBias" << endl;
  assert(bias.GetDim(0) >=startpos+1);
  assert(bias.GetDim(1) >=5);

  Array2 <doublevar> bias_temp(1,6);
  getLap(wfdata, e, bias_temp);
  int shiftf=1+startpos;
  for(int i=0; i< 5; i++)
      bias(shiftf,i)=bias_temp(0,i);
}

//------------------------------------------------------------------------

void Ext_slat::calcLap(Ext_slat_data * dataptr, Sample_point * sample)
{
  //cout << "calcLap " << endl;
  for(int e=0; e< nelectrons; e++) {
    parent->molecorb->updateLap(sample, e, 0,updatedMoVal);
    for(int d=0; d< 5; d++) {
      for(int i=0; i< updatedMoVal.GetDim(0); i++) {
        moVal(e,i,d)=updatedMoVal(i,d);
      }
    }
  }

  Array2 <doublevar> modet(nelectrons, nelectrons);
 
  for(int i=0; i< nelectrons; i++) {
    //cout << "modet: ";
    int s=parent->spin(i);
    for(int j=0; j< nelectrons; j++) {

      modet(i,j)=parent->amplitude(s,j)*(
                                         parent->pair_coeff(s,j,0)
                                         *moVal(i,parent->pairs(s,j,0),0) 
                                         +parent->pair_coeff(s,j,1)
                                         *moVal(i,parent->pairs(s,j,1),0));
      //cout << modet(i,j) << "  ";
    }
    //cout << endl;
  }

  detVal=TransposeInverseMatrix(modet, inverse, nelectrons);
  //cout << "detVal " << detVal << endl;


}

//------------------------------------------------------------------------


/*!
*/

void Ext_slat::getLap(Wavefunction_data * wfdata,
                     int e, Array2 <doublevar> & lap, int startpos)
{
  //cout << "---getLap---" << endl;

  assert(lap.GetDim(1) >= 5);
  assert(lap.GetDim(0) >= startpos);
  getVal(wfdata, e, lap, startpos);
  int s=parent->spin(e);
  for(int d=1; d< 5; d++) {
    doublevar temp=0;
    for(int i=0; i< nelectrons; i++) {

      temp+=inverse(e,i)
        *parent->amplitude(s,i)*(
                                 parent->pair_coeff(s,i,0)
                                 *moVal(e,parent->pairs(s,i,0),d) 
                                 +parent->pair_coeff(s,i,1)
                                 *moVal(e,parent->pairs(s,i,1),d));
      //cout << "amplitude " << parent->amplitude(s,i) << endl;
      
      //cout << "temp " << temp 
      //     << " pair " <<  parent->pair_coeff(s,i,0)
      //     << "  " << parent->pair_coeff(s,i,1)
      //     << "  vals " << moVal(e,parent->pairs(s,i,0),d) 
      //   << "  " << moVal(e,parent->pairs(s,i,1),d) << endl;
      //cout << "inverse " << inverse(e,i) << endl;
    }
    lap(startpos,d+1)=temp;
  }

  //cout << "lap  ";
  //for(int d=0; d< 6; d++) {
  //  cout << lap(startpos, d) << "   ";
  //}
  //cout << endl;

}

//-------------------------------------------------------------------------

/*!
*/
void Ext_slat::updateLap(
  Ext_slat_data * dataptr,
  Sample_point * sample,
  int e
)
{
  assert(dataptr != NULL);


  int s=dataptr->spin(e);
  sample->updateEIDist();
  doublevar ratio;

  int maxmatsize=nelectrons;
  static Array1 <doublevar> modet(maxmatsize);
  //cout << "electron " << e << " detVal " << detVal << endl;

  //check to make sure the determinant isn't zero
  if(detVal==0) {
    cout << "WARNING: Determinant zero" << endl;
    //error("Determinant zero ");
    calcLap(dataptr, sample);
    return;
  }


  //cout << "mo update\n";
  //update all the mo's that we will be using.
  dataptr->molecorb->updateLap(sample, e,
                              0,
                              updatedMoVal);

  //for(int i=0; i< nelectrons; i++) {
  //  cout << "updatedMoVal " << i << "  " << updatedMoVal(i,0) << endl;
  //}

  for(int i = 0; i < nelectrons; i++) {
    modet(i)=parent->amplitude(s,i)*(
                                     parent->pair_coeff(s,i,0)
                                     *updatedMoVal(parent->pairs(s,i,0),0) 
                                     +parent->pair_coeff(s,i,1)
                                     *updatedMoVal(parent->pairs(s,i,1), 0));
    //cout << "modet " << i << "  "<< modet(i) << endl;
  }

  //for(int i=0; i< nelectrons; i++) {
  //  cout << "inverse : ";
  //  for(int j=0; j< nelectrons; j++) {
  //    cout << inverse(i,j) << "  ";
  //  }
  //  cout << endl;
  //}
  
  doublevar tmpratio=InverseUpdateColumn(inverse,
                                         modet,e,
                                         nelectrons);

  ratio=1./tmpratio;
  if(tmpratio==0)
    ratio=0;

  //cout << "detval before " << detVal << "  ratio " << ratio << endl;
  detVal=ratio*detVal;


  for(int d=0; d< 5; d++)
  {
    for(int i=0; i< updatedMoVal.GetDim(0); i++)
    {
      moVal(e,i,d)=updatedMoVal(i,d);
    }
  }

}

//-------------------------------------------------------------------------

