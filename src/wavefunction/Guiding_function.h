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

#ifndef GUIDING_FUNCTION_H_INCLUDED
#define GUIDING_FUNCTION_H_INCLUDED

#include "Array.h"
#include "Qmc_std.h"



//----------------------------------------------------------------------

class Guiding_function {
 public:
  /*!
    \brief
    get weight of function, val(w)/guide
  */
  virtual doublevar getWeight(Wf_return & val, 
                              Wf_return & guide, 
                              int w)=0;

  /*!
    \brief
    condense to drift
   */
  virtual void getLap(Wf_return & lap, Array1 <doublevar> &)=0;

  /*!
    \brief
    get the ratio guide1/guide2
   */
  virtual doublevar getTrialRatio(Wf_return & wf1, 
                                  Wf_return & wf2)=0;

  virtual ~Guiding_function() {};
};


//----------------------------------------------------------------------

class Dmc_guiding_function:public Guiding_function {
 public:
  virtual void set_alpha_for_noderelease(doublevar & a){ }

  virtual doublevar getWeight(Wf_return & val, 
                              Wf_return & guide, 
                              int w)=0;
  virtual void getLap(Wf_return & lap, 
                      Array1 <doublevar> &)=0;
  virtual doublevar getTrialRatio(Wf_return & wf1, 
                                  Wf_return & wf2)=0;

  /*!
    \brief
    get the part for a linear operator 
  */
  virtual doublevar getOperatorWeight(Wf_return & lap, int w)=0;

};  

//----------------------------------------------------------------------

/*!
\brief
Represents the function \f$ \Psi_{guide}=\sum_i \Psi_i \f$
 */
/*
class Sum_guide:public Dmc_guiding_function
{
private:
  Array1 <doublevar> parm;
public:

  Sum_guide() { parm.Resize(1); parm(0)=1.0; }
  Sum_guide(Array1 <doublevar> & p) {
    parm.Resize(p.GetDim(0));
    parm=p;
  }

  virtual doublevar getOperatorWeight(Wf_return & vals, int w)
  {
    int nwf=vals.GetDim(0);
    assert(parm.GetDim(0)==nwf);

    doublevar ratio=0;
    doublevar total=0;
    doublevar wfval=0;
    for(int i=0; i<nwf; i++)
    {
      ratio=parm(i)*vals(0,0)*vals(i,0)*exp(vals(i,1)-vals(0,1));
      total+=ratio;
      if(i==w)
      {
        wfval=ratio;
      }
    }
    // cout << " weight " << w << "   " << wfval/total << endl;
    if(total==0)
    {
      cout << "WARNING: guiding function is zero \n";
      for(int i=0; i< nwf; i++)
      {
        cout << i << "   " << vals(i,0) << "  " << vals(i,1) << endl;
      }
      return 0;
    }
    return wfval/total;
  }


  virtual doublevar getWeight(Wf_return & vals, 
                              Wf_return & guide, int w)
  {
    int nwf=vals.GetDim(0);
    assert(parm.GetDim(0)==nwf);

    doublevar ratio=0;
    doublevar total=0;
    doublevar wfval=0;
    for(int i=0; i<nwf; i++)
    {
      ratio=guide(0,0)*guide(i,0)*exp(guide(i,1)-guide(0,1));
      total+=parm(i)*ratio;
      //if(i==w)
      //{
      //  wfval=ratio;
      //}
    }

    wfval=guide(0,0)*vals(w,0)*exp(vals(w,1)-guide(0,1));

    if(total==0)
    {
      cout << "WARNING: guiding function is zero \n";
      for(int i=0; i< nwf; i++)
      {
        cout << i << "   " << vals(i,0) << "  " << vals(i,1) << endl;
      }
      return 0;
    }
    return wfval/total;
  }


  void getLap(Wf_return & laps, Array1 <doublevar> & ret )
  {
    ret=0;
    int nwf=laps.GetDim(0);
    assert(parm.GetDim(0)==nwf);
    for(int w=0; w< nwf; w++)
    {
      doublevar weight=getOperatorWeight(laps,w);

      for(int d=2; d< laps.GetDim(1); d++)
      {
        ret(d)+=laps(w,d)*weight;
      }
    }
  }

  doublevar getTrialRatio(Wf_return & newfunc,
                          Wf_return & oldfunc)
  {
    int nwf=newfunc.GetDim(0);
    assert(newfunc.GetDim(0)==oldfunc.GetDim(0));
    Array1 <doublevar> ratio1(nwf);
    Array1 <doublevar> ratio2(nwf);

    doublevar sum1=0, sum2=0;
    //Sums relative to the first wave functions in each.
    for(int i=0; i<nwf; i++)
    {
      sum1+=parm(i)*newfunc(0,0)*newfunc(i,0)*exp(newfunc(i,1)-newfunc(0,1));
      sum2+=parm(i)*oldfunc(0,0)*oldfunc(i,0)*exp(oldfunc(i,1)-oldfunc(0,1));
    }
    //ratio of first wave functions.
    doublevar firstratio=newfunc(0,0)*oldfunc(0,0)
                         *exp(newfunc(0,1)-oldfunc(0,1));

    return firstratio*sum1/sum2;
  }
};
*/

//----------------------------------------------------------------------


/*!
\brief
Represents the function \f$ \Psi_{guide}^2 =\sum_i \Psi_i^2 \f$
 */
class Vmc_sum_squares:public Guiding_function
{
public:

  virtual doublevar getWeight(Wf_return & vals,
                              Wf_return & guide, int w) {
    doublevar denom=0;
    for(int i=0; i< guide.amp.GetDim(0) && i < 2; i++) {
      denom+=exp(2.0*(guide.amp(i,0)-vals.amp(w,0)));
    }
    return 1/sqrt(denom);
  }


  virtual void getLap(Wf_return & laps, Array1 <doublevar> & ret)
  {

    assert(ret.GetDim(0)>=3);
    int max=min(ret.GetDim(0)+1, laps.amp.GetDim(1));

    for(int d=1; d< max; d++) {
      ret(d-1)=laps.amp(0,d);
    }
  }

  virtual doublevar getTrialRatio(Wf_return & newfunc,
                                  Wf_return & oldfunc)
  {
    doublevar wfratio1=0;
    doublevar wfratio2=0;
    for(int i=0; i< newfunc.amp.GetDim(0) && i<2; i++)
    {
      doublevar tempratio1=exp(newfunc.amp(i,0)-newfunc.amp(0,0));
      doublevar tempratio2=exp(oldfunc.amp(i,0)-oldfunc.amp(0,0));
      wfratio1+=tempratio1*tempratio1;
      wfratio2+=tempratio2*tempratio2;
    }

    return exp((newfunc.amp(0,0)-oldfunc.amp(0,0)))*sqrt(wfratio1/wfratio2);
  }


};




//----------------------------------------------------------------------
/*!
\brief
Represents the function \f$ \Psi_{guide}^2 = \Psi_1^2 \f$
*/
class Primary:public Dmc_guiding_function {
public:

  virtual doublevar getWeight(Wf_return & vals,
                              Wf_return & guide, int w) {
    return exp((vals.amp(w,0)-guide.amp(0,0)));
  }


  virtual void getLap(Wf_return & laps, 
                      Array1 <doublevar> & ret) {
    assert(ret.GetDim(0) >= 3);
    assert(laps.amp.GetDim(0) >=1);
    int max=min(ret.GetDim(0)+1, laps.amp.GetDim(1));
    for(int d=1; d< max; d++)
      ret(d-1)=laps.amp(0,d);
  }

  virtual doublevar getOperatorWeight(Wf_return & lap, int w) {
    if(w==0) return 1;
    else return 0;
  }


  virtual doublevar getTrialRatio(Wf_return & newfunc,
                                  Wf_return & oldfunc) {
    return newfunc.sign(0)*oldfunc.sign(0)*exp((newfunc.amp(0,0)-oldfunc.amp(0,0)));
  }
};


//---------------------------------------------------------------------------
/*!
\brief
Represents the function \f$ \Psi_{guide} = \sqrt(\Psi_1^2+\alpha) \f$
*/
class Primary_noderelease:public Dmc_guiding_function {
 private:
  doublevar alpha;
 public:
  void set_alpha_for_noderelease(doublevar & a){
    alpha=a;
  }
  

  virtual doublevar getWeight(Wf_return & vals,
                              Wf_return & guide, int w) {
    
    //cout <<"alpha "<<alpha<<" guide.amp(0,0) "<<guide.amp(0,0)<<endl; 

    doublevar factor;
    if(alpha>0)
      factor= sqrt(1.0/(1.0+exp(log(alpha)-2.0*guide.amp(0,0))));
    else
      factor=1;

    //cout <<" getWeight:  "<<exp((vals.amp(w,0)-guide.amp(0,0)))*factor<<endl;

    return exp((vals.amp(w,0)-guide.amp(0,0)))*factor;
  }


  virtual void getLap(Wf_return & laps, 
                      Array1 <doublevar> & ret) {
    assert(ret.GetDim(0) >= 3);
    assert(laps.amp.GetDim(0) >=1);
    int max=min(ret.GetDim(0)+1, laps.amp.GetDim(1));

    doublevar factor;
    if(alpha>0)
      factor= 1.0/(1.0+exp(log(alpha)-2.0*laps.amp(0,0)));
    else
      factor=1;

    //cout <<" getLap: factor "<<factor<<endl;

    for(int d=1; d< max; d++)
      ret(d-1)=laps.amp(0,d)*factor;
  }

  virtual doublevar getOperatorWeight(Wf_return & lap, int w) {
    if(w==0) return 1;
    else return 0;
  }


  virtual doublevar getTrialRatio(Wf_return & newfunc,
                                  Wf_return & oldfunc) {
    doublevar factor;
      if(alpha>0)
	factor = sqrt((1.0+exp(log(alpha)-2.0*newfunc.amp(0,0)))/(1.0+exp(log(alpha)-2.0*oldfunc.amp(0,0))));
      else
	factor=1;

      //cout <<" getTrialRatio: factor "<<factor<<endl; 
    return factor*newfunc.sign(0)*oldfunc.sign(0)*exp((newfunc.amp(0,0)-oldfunc.amp(0,0)));
  }
};


#endif //GUIDING_FUNCTION_H_INCLUDED
//------------------------------------------------------------------------
