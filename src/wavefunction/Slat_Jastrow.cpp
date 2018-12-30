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

#include "Slat_Jastrow.h"
#include "qmc_io.h"
#include "Qmc_std.h"
#include "Sample_point.h"

void Slat_Jastrow::notify(change_type change , int num )
{
  slater_wf->notify(change, num);
  jastrow_wf->notify(change, num);
}



void Slat_Jastrow::init(Wavefunction_data * wfdata)
{
  Slat_Jastrow_data * dataptr;
  recast(wfdata, dataptr);
  dataptr->slater->generateWavefunction(slater_wf);
  dataptr->jastrow->generateWavefunction(jastrow_wf);
  if(jastrow_wf->nfunc() != slater_wf->nfunc())
  {
    error("You must put the same number of functions in the "
          "functions to be multiplied in SLAT-JASTROW.");
  }
  nfunc_=slater_wf->nfunc();
}



void Slat_Jastrow::getVal(Wavefunction_data * wfdata,
                          int e, Wf_return & val)
{
  assert(val.amp.GetDim(0) >= nfunc_);
  Wf_return slat_val(nfunc_, 2);
  Wf_return jast_val(nfunc_, 2);
  Slat_Jastrow_data * dataptr;
  recast(wfdata, dataptr);
  assert(dataptr != NULL);
  
  slater_wf->getVal(dataptr->slater, e, slat_val);
  jastrow_wf->getVal(dataptr->jastrow, e, jast_val);

  if ( slat_val.is_complex==1 || jast_val.is_complex==1 )
    val.is_complex=1;


  for(int i=0; i< nfunc_; i++) {
    val.phase(i,0)=slat_val.phase(i,0)+jast_val.phase(i,0); 
    val.amp(i,0)=slat_val.amp(i,0)+jast_val.amp(i,0);  //add the logarithm
  }
}


void Slat_Jastrow::evalTestPos(Array1 <doublevar> & pos, Sample_point * sample,Array1 <Wf_return> & wf) {
  Array1 <Wf_return> slat_val,jast_val;
  slater_wf->evalTestPos(pos,sample,slat_val);
  jastrow_wf->evalTestPos(pos,sample,jast_val);
  assert(slat_val.GetDim(0)==jast_val.GetDim(0));
  int n=slat_val.GetDim(0);
  wf.Resize(n);
  for(int i=0; i< n; i++) { 
    wf(i).Resize(nfunc_,2);
    if ( slat_val(i).is_complex==1 || jast_val(i).is_complex==1 )
      wf(i).is_complex=1;
    
    for(int j=0; j< nfunc_; j++) { 
      wf(i).phase(j,0)=slat_val(i).phase(j,0)+jast_val(i).phase(j,0); 
      wf(i).amp(j,0)=slat_val(i).amp(j,0)+jast_val(i).amp(j,0);  //add the logarithm
    }
  }

}



void Slat_Jastrow::getSymmetricVal(Wavefunction_data * wfdata,
			     int e, Wf_return & val){

  Wf_return slat_val(nfunc_, 2);
  Wf_return jast_val(nfunc_, 2);
 
  Slat_Jastrow_data * dataptr;
  recast(wfdata, dataptr);
  assert(dataptr != NULL);
  
  slater_wf->getSymmetricVal(dataptr->slater, e, slat_val);
  jastrow_wf->getSymmetricVal(dataptr->jastrow, e, jast_val);

  //cout <<"slat_val.amp(0,0) "<<slat_val.amp(0,0)<<"  jast_val.amp(0,0) "<<jast_val.amp(0,0)<<endl;

  for(int i=0; i< nfunc_; i++) {
    val.amp(i,0)=slat_val.amp(i,0)+jast_val.amp(i,0);  //add the logarithm
  }
}





int Slat_Jastrow::getParmDeriv(Wavefunction_data *  wfdata, 
			       Sample_point * sample ,
			       Parm_deriv_return & derivatives){

  //cout <<"Slat_Jastrow::getParmDeriv"<<endl;
  Slat_Jastrow_data * dataptr;
  recast(wfdata, dataptr);
  assert(dataptr != NULL);
   
  Parm_deriv_return slaterval;
  Parm_deriv_return jastval;
  slaterval.need_hessian=jastval.need_hessian=derivatives.need_hessian;
  
  int nslater=dataptr->slater->nparms();
  int njast=dataptr->jastrow->nparms();
  int nparms=nslater + njast;
  derivatives.val_gradient.Resize(sample->electronSize(),5);
  slater_wf->getParmDeriv(dataptr->slater,sample, derivatives);
  jastrow_wf->getParmDeriv(dataptr->jastrow,sample, jastval);
  extend_parm_deriv(derivatives,jastval);
 
  return 1;
}



// JK: complex wavefunctions still behave weird, let's try to rewrite even
// this one, not touching cvals now, though
// The most important is setting is_complex, but I do differ at other places
// too.
void Slat_Jastrow::getLap(Wavefunction_data * wfdata, int e,
                          Wf_return & lap)
{
  //cout << "getLap\n";
  assert(lap.amp.GetDim(0) >= nfunc_);
  assert(lap.amp.GetDim(1) >= 5);
  Slat_Jastrow_data * dataptr;
  recast(wfdata, dataptr);
  assert(dataptr != NULL);

  Wf_return slat_lap(nfunc_,5);
  Wf_return jast_lap(nfunc_,5);

  slater_wf->getLap(dataptr->slater, e, slat_lap);
  jastrow_wf->getLap(dataptr->jastrow, e, jast_lap);

  if ( slat_lap.is_complex==1 || jast_lap.is_complex==1 )
    lap.is_complex=1;

  for(int i=0; i< nfunc_; i++)
  {

    lap.phase(i,0)=slat_lap.phase(i,0)+jast_lap.phase(i,0);
    lap.amp(i,0)=slat_lap.amp(i,0)+jast_lap.amp(i,0);
    doublevar dotproduct=0;
    dcomplex dot(0.0, 0.0);
    for(int d=1; d<4; d++)
    {
      lap.amp(i,d)=slat_lap.amp(i,d)+jast_lap.amp(i,d);
      lap.phase(i,d)=slat_lap.phase(i,d)+jast_lap.phase(i,d);
      dotproduct+=slat_lap.amp(i,d)*jast_lap.amp(i,d);

      lap.cvals(i,d)=slat_lap.cvals(i,d)+jast_lap.cvals(i,d);
      dot+=slat_lap.cvals(i,d)*jast_lap.cvals(i,d);
    }

    lap.amp(i,4)=slat_lap.amp(i,4)+jast_lap.amp(i,4)
      +2*dotproduct;
    lap.phase(i,4)=slat_lap.phase(i,4)+jast_lap.phase(i,4);

    lap.cvals(i,4)=slat_lap.cvals(i,4)+jast_lap.cvals(i,4)
      +2.0*dot;

  }
}


//----------------------------------------------------------------------

void Slat_Jastrow::updateVal(Wavefunction_data * wfdata, Sample_point * sample)
{
  Slat_Jastrow_data * dataptr;
  recast(wfdata, dataptr);
  assert(dataptr != NULL);
  slater_wf->updateVal(dataptr->slater, sample);
  jastrow_wf->updateVal(dataptr->jastrow, sample);
}

void Slat_Jastrow::updateLap(Wavefunction_data * wfdata, Sample_point * sample)
{
  Slat_Jastrow_data * dataptr;
  recast(wfdata, dataptr);
  assert(dataptr != NULL);
  slater_wf->updateLap(dataptr->slater, sample);
  jastrow_wf->updateLap(dataptr->jastrow, sample);
}


//----------------------------------------------------------------------


void Slat_Jastrow::plot1DInternals(Array1 <doublevar> & xdata,
				   vector <Array1 <doublevar> > & data,
				   vector <string> & desc,
				   string desc0) {

  slater_wf->plot1DInternals(xdata,data,desc,desc0);
  jastrow_wf->plot1DInternals(xdata,data,desc,desc0);

}


//----------------------------------------------------------------------
