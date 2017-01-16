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
// src/Jastrow_wf.cpp
//
//
#include "Qmc_std.h"
#include "Jastrow_wf.h"
#include "Sample_point.h"


//----------------------------------------------------------------------
void Jastrow_wf::notify(change_type change, int num)
{
  switch(change)
  {
  case electron_move:
    electronIsStaleVal(num)=1;
    updateEverythingLap=1;
    break;
  case all_electrons_move:
    updateEverythingVal=1;
    updateEverythingLap=1;
    break;
  case wf_parm_change:
  case all_wf_parms_change:
    parmChanged=1;
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
    updateStatic=1;
    break;
  case sample_dynamic:
    staticSample=0;
    updateStatic=0;
    break;
  default:
    updateEverythingVal=1;
    updateEverythingLap=1;
  }
}


//----------------------------------------------------------------------

void Jastrow_wf::init(Wavefunction_data * wfdata)
{
  Jastrow_wf_data * dataptr;
  recast(wfdata, dataptr);
  //cout << "Jastrow_wf" << endl;
  nelectrons=dataptr->nelectrons;
  nions=dataptr->nions;

  derivatives.Resize(nelectrons, 5);
  valPartialSum.Resize(nelectrons);
  ionPartialSum.Resize(nelectrons);

  elecPartialSum.Resize(nelectrons, nelectrons);
  a_kval.Resize(nions, nelectrons);
  for(int i=0; i< nions; i++)
  {
    for(int j=0; j < nelectrons; j++)
    {
      a_kval(i,j).Resize(dataptr->elecIonBasis->nfunc()+1);
      a_kval(i,j)(0)=1;
    }
  }
  electronIsStaleVal.Resize(nelectrons);
  electronIsStaleVal=0;
  updateEverythingVal=1;
  updateEverythingLap=1;
  sampleAttached=0;
  dataAttached=0;
  staticSample=0;
  updateStatic=0;
  parmChanged=0;
  // cout << "done" << endl;

}

//----------------------------------------------------------------------

void Jastrow_wf::updateVal(Wavefunction_data * wfdata,
                           Sample_point * sample)
{
  assert(sampleAttached);
  assert(dataAttached);

  Jastrow_wf_data * jastdata;
  recast(wfdata, jastdata);


  if(updateEverythingVal==1)
  {
    calcVal(jastdata, sample);
    updateEverythingVal=0;
    electronIsStaleVal=0;
  }
  else
  {
    for(int e=0; e< nelectrons; e++)
    {
      if(electronIsStaleVal(e)==1)
      {
        updateVal(jastdata, sample, e);
        //calcLap(jastdata, sample);
        electronIsStaleVal(e)=0;
      }
    }
  }

}

//----------------------------------------------------------------------

void Jastrow_wf::updateForceBias(Wavefunction_data * wfdata,
                                 Sample_point * sample)
{
  assert(sampleAttached);
  assert(dataAttached);

  updateVal(wfdata, sample);
}


//----------------------------------------------------------------------

void Jastrow_wf::updateLap(Wavefunction_data * wfdata, Sample_point * sample)
{
  assert(sampleAttached);
  assert(dataAttached);

  Jastrow_wf_data * jastdata;
  recast(wfdata, jastdata);

  //cout << "updateEverythingLap " << updateEverythingLap << endl;

  //For a static sample, we trust that only the parameters
  //are changing
  //If we're optimizing the basis, we have to recalculate everything
  if(staticSample && parmChanged && !(jastdata->optimize_basis))
  {
    if(updateStatic)
    {  //if we haven't saved since set static
      jastdata->makeStaticSave(sample, values_static,
                               partialValues_static, deriv_static);
      updateStatic=0;
    }
    jastdata->updateParms(sample, values_static,partialValues_static,
                          deriv_static, value,valPartialSum, derivatives);
    parmChanged=0;
  }
  else
  { //Otherwise do a regular update
    if(updateEverythingLap==1 || parmChanged==1)
    {

      calcLap(jastdata, sample);
      updateEverythingVal=0;
      updateEverythingLap=0;
      electronIsStaleVal=0;
      parmChanged=0;
    }
    else
    {
      assert(updateEverythingVal==0);
      for(int e=0; e< nelectrons; e++)
      {
        assert(electronIsStaleVal(e) == 0);
      }
    }
  }

}

//----------------------------------------------------------------------


void Jastrow_wf::storeParmIndVal(Wavefunction_data * wfdata,
                                 Sample_point * sample,
                                 int e, Array1 <doublevar> & vals )
{
  Jastrow_wf_data * dataptr;
  recast(wfdata, dataptr);
  assert(vals.GetDim(0) >= dataptr->valSize());

  if(dataptr->optimize_basis) {
    //do nothing.
  }
  else {
    dataptr->fillParmInd(sample, e, vals);
  }
}




//----------------------------------------------------------------------
/*!

 */
void Jastrow_wf::getParmDepVal(Wavefunction_data * wfdata,
                               Sample_point * sample,
                               int e,
                               Array1 <doublevar> & oldval,
                               Array2 <doublevar> & newval)
{
  assert(staticSample);
  Jastrow_wf_data * dataptr;
  recast(wfdata, dataptr);
  assert(staticSample==1);

  if(dataptr->optimize_basis) {
    //cout << "value before  " << value << endl;
    updateVal(dataptr, sample, e);
    //cout << "value after " << value << endl;
    getVal(dataptr, e, newval);
  }
  else {
    assert(updateStatic==0);

    doublevar newvalue=0;
    doublevar oldvalue=valPartialSum(e);
    //cout << "oldvalue " << oldvalue << endl;

    dataptr->updateParmsVal(sample, e, newvalue, oldval);
    //cout << "updateParmsnewvalue " << newvalue << endl;
    newval(0,0)=1;
    newval(0,1)=value+newvalue-oldvalue-dataptr->normalization;
    //cout << "value " << value << endl;
  }

  //cout << "newvalue " << newval(0,0) << "  " << newval(0,1) << endl;

}
//----------------------------------------------------------------------


int Jastrow_wf::getParmDeriv(Wavefunction_data *  wfdata, 
			  Sample_point * sample ,
			  Parm_deriv_return & derivatives){

  error("Not implemented for this wavefunction");
  /*
  assert(staticSample);
  Jastrow_wf_data * dataptr;
  recast(wfdata, dataptr);
  assert(staticSample==1);
  
  derivatives.gradient.Resize(dataptr->nparms());
  derivatives.hessian.Resize(dataptr->nparms(),dataptr->nparms());

  if(dataptr->nparms()) {
    error("Not implemented yet!");
    return 0;
  }
  else { 
    derivatives.gradient=0;
    derivatives.hessian=0;
    return 1;
  }
  */
  return 0;
}

//----------------------------------------------------------------------

void Jastrow_wf::generateStorage(Wavefunction_storage * & wfstore)
{
  wfstore=new Jastrow_wf_storage;
  Jastrow_wf_storage * store;
  recast(wfstore, store);
  store->elecPartialSum_temp.Resize(nelectrons);
  store->a_kval_temp.Resize(nions);
  for(int i=0; i< nions; i++)
  {
    store->a_kval_temp(i).Resize(a_kval(i,0).GetDim(0));
  }
}

//----------------------------------------------------------------------

void Jastrow_wf::saveUpdate(Sample_point * sample, int e, 
                            Wavefunction_storage * wfstore)
{
  //cout << "Jastrow::save" << endl;
  Jastrow_wf_storage * store;
  recast(wfstore, store);

  store->value_temp=value;

  store->valPartialSum_temp=valPartialSum(e);

  for(int i=0; i< e; i++)
  {
    store->elecPartialSum_temp(i)=elecPartialSum(i,e);
  }
  for(int j=e+1; j< nelectrons; j++)
  {
    store->elecPartialSum_temp(j)=elecPartialSum(e,j);
  }

  for(int i=0; i< nions; i++) {
    store->a_kval_temp(i)=a_kval(i,e);
  }

  store->ionPartialSum_temp=ionPartialSum(e);

}

//----------------------------------------------------------------------

void Jastrow_wf::restoreUpdate(Sample_point * sample, int e, 
                               Wavefunction_storage * wfstore)
{

  //cout << "Jastrow::restore" << endl;
  Jastrow_wf_storage * store;

  recast(wfstore, store);
  //cout << "ionPartialSum " << ionPartialSum(e) << "    "
  //     << store->ionPartialSum_temp << endl;

  ionPartialSum(e)=store->ionPartialSum_temp;

  for(int i=0; i< e; i++)
  {
    //cout << "elecPartialSum " << elecPartialSum(i,e)
    //   << "   " << store->elecPartialSum_temp(i) << endl;
    elecPartialSum(i,e)=store->elecPartialSum_temp(i);
  }
  for(int j=e+1; j< nelectrons; j++)
  {
    //cout << "elecPartialSum " << elecPartialSum(e,j)
    //   << "   " << store->elecPartialSum_temp(j) << endl;
    elecPartialSum(e,j)=store->elecPartialSum_temp(j);
  }


  for(int i=0; i< nions; i++)
  {
    a_kval(i,e)=store->a_kval_temp(i);
  }
  //cout << "value " << value << "   " << store->value_temp << endl;
  value=store->value_temp;
  valPartialSum(e)=store->valPartialSum_temp;
  updateEverythingVal=0;
  electronIsStaleVal=0;
}

//------------------------------------------------------------------------

void Jastrow_wf::calcVal(Jastrow_wf_data * dataptr, Sample_point * sample)
{

  for(int e=0; e< nelectrons; e++)
  {
    updateVal(dataptr, sample, e);
  }
}

//------------------------------------------------------------------------

/*!

 */
void Jastrow_wf::updateVal(Jastrow_wf_data * dataptr, Sample_point * sample,int e)
{



  Array1 <doublevar> elecPartialSum_temp(nelectrons);
  doublevar ionPartialSum_temp;
  //Save old values
  for(int i=0; i< e; i++)
  {
    elecPartialSum_temp(i)=elecPartialSum(i,e);
    //cout << "elecPartialSum_temp(" << i << ")"
    // << elecPartialSum_temp(i) << endl;
  }
  for(int j=e+1; j< nelectrons; j++)
  {
    elecPartialSum_temp(j)=elecPartialSum(e,j);
    //cout << "elecPartialSum_temp(" << j << ")"
    // << elecPartialSum_temp(j) << endl;

  }
  ionPartialSum_temp=ionPartialSum(e);
  //cout << "ionPartialSum_temp  " << ionPartialSum_temp << endl;



  dataptr->updateVal(a_kval, sample, e,ionPartialSum, elecPartialSum);

  //Assign the value
  value+=ionPartialSum(e)-ionPartialSum_temp;
  //cout << "ionpartialSum " << ionPartialSum(e) << endl;
  for(int i=0; i< e; i++)
  {
    value+=elecPartialSum(i,e)-elecPartialSum_temp(i);
    //cout << "elecPartialSum(" << i << ") " << elecPartialSum(i,e)
    //   <<endl;
  }

  for(int j=e+1; j<nelectrons; j++)
  {
    value+=elecPartialSum(e,j)-elecPartialSum_temp(j);
    // cout << "elecPartialSum(" << j << ") " << elecPartialSum(e,j)
    // <<endl;
  }




}

//------------------------------------------------------------------------

void Jastrow_wf::getVal(Wavefunction_data * wfdata, int e, Array2 <doublevar> & val, int startpos)
{
  assert(val.GetDim(0) >= startpos+1);
  Jastrow_wf_data * jastdata;
  recast(wfdata, jastdata);
  //cout << "getVal " << value << endl;
  val(startpos,0)=1;  //always positive
  val(startpos,1)=value-jastdata->normalization;
}

//----------------------------------------------------------------------

void Jastrow_wf::getDensity(Wavefunction_data * wfdata,
                            int e,
                            Array2 <doublevar> & dens)
{
  dens(0,0)=1;
}

//----------------------------------------------------------------------------

void Jastrow_wf::getForceBias(Wavefunction_data * wfdata, int e, 
                              Array2 <doublevar> & bias, int startpos)
{
  bias=0;
  getVal(wfdata, e,bias, startpos);
}

//------------------------------------------------------------------------
/*!

*/
void Jastrow_wf::calcLap(Jastrow_wf_data * dataptr, Sample_point * sample)
{
  //cout << "calcLap" << endl;;
  assert(dataptr != NULL);



  //Get the necessary basis functions

  //Electron-electron
  //cout << "Electron-electron\n";
  Array1 <doublevar> eedist(5);
  Array1 <doublevar> temp(5);
  
  sample->updateEIDist();
  sample->updateEEDist();

  for(int i=0; i< nelectrons; i++)
  {

    for(int j=i+1; j< nelectrons; j++)
    {
      //cout << "getDist" << endl;
      sample->getEEDist(i,j, eedist);
      //cout << "calcLap " << endl;

      dataptr->elecElecBasis->
      calcLap(eedist, dataptr->b_m(i,j), dataptr->basisoffset);
      //cout << "cusp " << endl;
      //cout << "eedist " << eedist(0) << endl;
      dataptr->calcEECuspLap(i,j,eedist);
      //cout << "done2 " << endl;
    }
  }

  //cout << "Electron-ion\n";
  //Electron-ion
  Array1 <doublevar> eidist(5);
  for(int j=0; j< nelectrons; j++)
  {
    for(int i=0; i< nions; i++)
    {
      sample->getEIDist(j,i,eidist);
      //cout << "eidist " << eidist(0) << endl;
      dataptr->elecIonBasis->
      calcLap(eidist,dataptr->a_k(i,j), dataptr->basisoffset);
      int max=dataptr->a_k(i,j).GetDim(0);
      for(int k=0; k< max; k++) {
        a_kval(i,j)(k)=dataptr->a_k(i,j)(k,0);
      }
    }
  }

  //Initialize sum variables
  derivatives=0;
  ionPartialSum=0;
  elecPartialSum=0;
  //Calculate the value, gradient, and laplacian

  //   cout << "sumElecIon\n";
  //Electron-ion correlation

  dataptr->sumElecIon(ionPartialSum, derivatives);
  //cout << "sumElecElec\n";
  //electron-electron and electron-electron-ion correlations
  dataptr->sumElecElec(elecPartialSum, derivatives);

  value=0;
  for(int i=0; i< nelectrons; i++) {
    value+=ionPartialSum(i);
  }
  //cout << "u from ion " << value << endl;
  doublevar tot=0;
  for(int i=0; i< nelectrons; i++)
  {

    //cout << "ionPartialSum " << ionPartialSum(i) << endl;
    for(int j=i+1; j< nelectrons; j++)
    {
      tot+=elecPartialSum(i,j);
      //cout << "elecPartialSum " << elecPartialSum(i,j) << endl;
      //cout << "calcLap " << value << endl;
    }
  }
  value+=tot;
  //cout << "u from electron " << tot << endl;



}


//------------------------------------------------------------------------


void Jastrow_wf::getLap(Wavefunction_data * wfdata, int e,
                        Array2 <doublevar> & lap, int startpos)
{
  assert(lap.GetDim(0) >=startpos+1);
  assert(lap.GetDim(1) >=6 );
  Jastrow_wf_data * jastdata;
  recast(wfdata, jastdata);

  lap(startpos,0)=1;
  lap(startpos,1)=value-jastdata->normalization;

  doublevar dotproduct=0;
  for(int d=1; d< 4; d++)
  {
    lap(startpos,d+1)=derivatives(e,d);
    dotproduct+=derivatives(e,d)*derivatives(e,d);
    //cout << "derivative " << e << "  " << d+1 << "  " << derivatives(e,d)
    //     << endl;
  }
  lap(startpos,5)=derivatives(e,4)+dotproduct;
  //cout << "derivative " << e << "  " << "5" << "  " << lap(0,5)
  // << endl;
}

//-------------------------------------------------------------------------



