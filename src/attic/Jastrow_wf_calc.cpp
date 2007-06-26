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
//----------------------------------------------------------------------
//src/Jastrow_wf_calc.cpp

#include "Qmc_std.h"
#include "Jastrow_wf_data.h"
#include "Sample_point.h"


void Jastrow_wf_data::calcEECuspVal(int i, int j, Array1 <doublevar> & eedist)
{
  assert(i<j);
  //cout << "cuspval  " << endl;

  int jspin=spin(j);
  int ispin=spin(i);

  //Calculate cusp
  static Array1 <doublevar> tempval(1);
  if(jspin==ispin) {
    cuspBasis(0)->calcVal(eedist, tempval);
  }
  else {
    cuspBasis(1)->calcVal(eedist, tempval);
  }

  eecusp(0,i,j)=tempval(0);
  //cout << "done " << endl;
}


void Jastrow_wf_data::calcEECuspLap(int i, int j, Array1 <doublevar> & eedist)
{
  //cout << "cusplap " << endl;
  assert(i<j);


  int jspin=spin(j);
  int ispin=spin(i);

  static Array2 <doublevar> tempvals(1,5);
  if(jspin==ispin) {
    cuspBasis(0)->calcLap(eedist, tempvals);
  }
  else {
    cuspBasis(1)->calcLap(eedist, tempvals);
  }

  for(int k=0; k< 5; k++) {
    eecusp(k,i,j)=tempvals(0,k);
  }


}

//----------------------------------------------------------------------

void Jastrow_wf_data::sumElecIon(Array1 <doublevar> & ionPartialSum,
                                 Array2 <doublevar> & derivatives)
{


  for(int I=0; I< nions; I++)
  {
    int k=0;

    for(int seq=eiCorr_start(I);
        seq < eiCorr_end(I); seq++)
    {
      //cout << "seq " << seq << endl;
      doublevar coeff=facco*eiCorrelation(seq);
      //cout << "coeff " << coeff << endl;
      for(int e=0; e< nelectrons; e++)
      {
        ionPartialSum(e)+=coeff*a_k(I,e)(k+1,0);
        //cout << "ionPartialSum("<< e<<") " << ionPartialSum(e)
        //    << " coeff " << coeff << " basis " << a_k(I,e)(k+1, 0) << endl;
        for(int d=1; d< 5; d++)
        {
          derivatives(e,d)+=coeff*a_k(I,e)(k+1,d);
          //cout << "deriv " << e << "   " << d << "  "
          // << derivatives(e,d) << endl;
        }
        //cout << "ecor " << derivatives(e,4)*(-.5) << endl;

      }
      k++;  //Count the function on this center.
    }
  }



}

//----------------------------------------------------------------------

/*!
adds to value, valPartialSum, and derivatives.
*/
void Jastrow_wf_data::sumElecElec(Array2 <doublevar> & elecPartialSum,
                                  Array2 <doublevar> & derivatives)
{




  //cout << "eecusp " << endl;
  for(int i=0; i< nelectrons; i++)
  {
    for(int j=i+1; j< nelectrons; j++)
    {
      elecPartialSum(i,j)+=eecusp(0,i,j);
      //cout << "eecusp " << eecusp(0,i,j);
      //cout << "elecPartialSum(" << i << "," << j
      //     << ")  " << elecPartialSum(i,j) << endl;
    }
  }
  for(int d=1; d< 4; d++)
  {
    for(int i=0; i< nelectrons; i++)
    {
      for(int j=i+1; j< nelectrons; j++)
      {
        derivatives(i,d)+=eecusp(d,i,j);
        derivatives(j,d)-=eecusp(d,i,j);
      }
    }
  }
  for(int i=0; i< nelectrons; i++)
  {
    for(int j=i+1; j< nelectrons; j++)
    {
      derivatives(i,4)+=eecusp(4,i,j);
      derivatives(j,4)+=eecusp(4,i,j);
    }
  }

  Array3 <doublevar> driftij(nelectrons, nelectrons,5);
  Array3 <doublevar> driftion_i(nelectrons, nelectrons,5);
  Array3 <doublevar> driftion_j(nelectrons, nelectrons,5);

  driftij=0;
  driftion_i=0;
  driftion_j=0;

  //----------e-e correlation----------
  //cout << "eecorr" << endl;
  for(int seq=0;
      seq< eeCorrelation.GetDim(0);
      seq++)
  {

    doublevar coeff=eeCorrelation(seq)*faccp;
    for(int i=0; i< nelectrons; i++)
    {
      for(int j=i+1; j< nelectrons; j++)
      {
        elecPartialSum(i,j)+=coeff*b_m(i,j)(seq+1,0);

        //cout << "elecPartialSum(" << i << "," << j
        //     << ")  " << elecPartialSum(i,j)
        //     << " coeff " << coeff
        //     << " basis " << b_m(i,j)(seq+1,0)
        //     << endl;
        for(int d=1; d<5; d++)
        {
          driftij(i, j,d)+=coeff*b_m(i,j)(seq+1,d);
        }
        //cout << "add to laplacian " << -.5*coeff*b_m(i,j)(seq+1, 4) << endl;
      }
    }
  }


  //-----electron electron ion---------
  //Array1 <doublevar> v_kl_i(5);
  //Array1 <doublevar> v_kl_j(5);
  //cout << "eeicorr" << endl;
  for(int i=0; i< nelectrons; i++)
  {
    int istride=i*nelectrons;
    doublevar dotproduct_i, dotproduct_j;
    int ijd;
    doublevar v_klval;
    doublevar v_kl_i[5], v_kl_j[5];
    doublevar bmij[5];
    doublevar akj_temp, alj_temp;
    doublevar aki_drift[5], ali_drift[5];

    for(int j=i+1; j< nelectrons; j++)
    {
      for(int I=0; I< nions; I++)
      {
        int ionstride=I*a_k.GetDim(1);
        int ionstridei=ionstride+i;
        int ionstridej=ionstride+j;
        int istridej=istride+j;

        for(int seq=eeiCorr_start(I);
            seq < eeiCorr_end(I);
            seq++)
        {
          int c_kt=c_k(seq); //temporary variable lookups
          int c_lt=c_l(seq);
          int c_mt=c_m(seq);
          doublevar coeff=faccp*eeiCorrelation(seq);
          c_kt*=5;
          c_mt*=5;
          c_lt*=5; //account for the stride.

          doublevar aki_temp=a_k.v[ionstridei].v[c_kt];
          doublevar ali_temp=a_k.v[ionstridei].v[c_lt];

          for(int d=1; d< 5; d++)
          {
            aki_drift[d]=a_k.v[ionstridei].v[c_kt+d];
            ali_drift[d]=a_k.v[ionstridei].v[c_lt+d];
          }


          for(int d=0; d<5; d++)
          {
            bmij[d]=b_m.v[istridej].v[c_mt+d];
          }


          akj_temp=a_k.v[ionstridej].v[c_kt];
          alj_temp=a_k.v[ionstridej].v[c_lt];


          v_klval=coeff*(aki_temp*alj_temp+ali_temp*akj_temp);
          elecPartialSum(i,j)+=v_klval*bmij[0];
          //cout << "elecPartialSum(" << i << "," << j
          //   << ")  " << elecPartialSum(i,j) << endl;
          for(int d=1; d< 5; d++)
          {
            v_kl_i[d]=coeff*(aki_drift[d]*alj_temp
                             +ali_drift[d]*akj_temp);

            v_kl_j[d]=coeff*(aki_temp*a_k.v[ionstridej].v[c_lt+d]
                             +ali_temp*a_k.v[ionstridej].v[c_kt+d]);
          }

          ijd=(istridej)*5;
          for(int d=1; d< 5; d++)
          {
            //ijd=(istridej)*5+d;
            ijd++;

            driftij.v[ijd]+=v_klval*bmij[d];
            driftion_i.v[ijd]+=v_kl_i[d]*bmij[0];
            driftion_j.v[ijd]+=v_kl_j[d]*bmij[0];

          }

          //cout << "done " << endl;

          //Cross terms in the laplacian
          dotproduct_i=0;
          dotproduct_j=0;
          for(int d=1; d<4; d++)
          {
            dotproduct_i +=v_kl_i[d]*bmij[d];//*b_m.v[istride+j].v[c_mt+d];
            dotproduct_j +=v_kl_j[d]*bmij[d];//b_m.v[istride+j].v[c_mt+d];

          }
          //since we don't change the
          //sign for the lap later,
          //we need to do it here.
          //driftion_i(i,j,4)+=2.0*dotproduct_i;
          //driftion_j(i,j,4)-=2.0*dotproduct_j;
          int driftoffset=(istride+j)*5+4;
          driftion_i.v[driftoffset]+=2.0*dotproduct_i;
          driftion_j.v[driftoffset]-=2.0*dotproduct_j;
        }
      }
    }
  }


  //cout << "gradient adding\n";
  //add up the gradient for i and j
  //Now sum up the derivatives.
  for(int i=0; i< nelectrons; i++)
  {
    for(int j=i+1; j< nelectrons; j++)
    {
      for(int d=1; d<4; d++)
      {
        derivatives(i,d)+=driftion_i(i,j,d)+driftij(i,j,d);
        derivatives(j,d)+=driftion_j(i,j,d)-driftij(i,j,d);
      }

      //Laplacian is added for both i and j
      derivatives(i,4)+=driftion_i(i,j,4)+driftij(i,j,4);
      derivatives(j,4)+=driftion_j(i,j,4)+driftij(i,j,4);
      //cout << "pair " << i << "   " << j
      //     << " driftion " << driftion_i(i,j,4)
      //     << " driftij " << driftij(i,j,4) << endl;
    }
  }
  // cout << "done" << endl;

}

//----------------------------------------------------------------------

/*!
 */
void Jastrow_wf_data::updateVal(Array2 <Array1 <doublevar> > & a_kval,
                                Sample_point * sample, int e,
                                Array1 <doublevar> & ionPartialSum,
                                Array2 <doublevar> & elecPartialSum
                               )
{
  //cout << "updateVal\n";
  Array1 <doublevar> eedist(5);
  Array1 <doublevar> temp(5);

  sample->updateEIDist();
  sample->updateEEDist();

  //Calculate the updated basis functions for electron-electron
  for(int i=0; i< e; i++)
  {
    sample->getEEDist(i,e,eedist);
    //for(int k=0; k< 5; k++) {
    //  cout << "eedist " << k << "  " << eedist(k) << endl;
    //}
    calcEECuspVal( i,e,eedist);
    elecElecBasis->calcVal(eedist, bupdate(i),basisoffset );
  }
  for(int j=e+1; j< nelectrons; j++)
  {
    sample->getEEDist(e,j,eedist);
    calcEECuspVal(e,j,eedist);
    elecElecBasis->calcVal(eedist,  bupdate(j), basisoffset);
  }

  //Calculate the updated basis functions for electron-ion
  Array1 <doublevar> eidist(5);
  Array1 <doublevar> tempa(elecIonBasis->nfunc());
  for(int i=0; i< nions; i++)
  {
    //for(int j=0; j< nelectrons; j++)
    //{
    //  sample->getEIDist(j,i,eidist);

    // elecIonBasis->calcVal(eidist, a_kval(i,j), basisoffset);
    //elecIonBasis->calcLap(eidist, temp, a_k(i,j));
    //}
    sample->getEIDist(e,i,eidist);
    elecIonBasis->calcVal(eidist, a_kval(i,e), basisoffset);
  }

  //---------electron-ion----------


  ionPartialSum(e)=0;
  for(int I=0; I< nions; I++)
  {
    int k=0;
    for(int seq=eiCorr_start(I);
        seq < eiCorr_end(I); seq++)
    {
      doublevar coeff=facco*eiCorrelation(seq);
      //ionPartialSum(e)+=coeff*a_k(I,e)(k+1,0);
      ionPartialSum(e)+=coeff*a_kval(I,e)(k+1);
      k++;
    }
  }

  //int estride=a_k.GetDim(1);
  int estride=a_kval.GetDim(1);

  //Electron-electron interaction



  //-----electron electron cusp----------
  for(int i=0; i< e; i++)
  {
    elecPartialSum(i,e)=eecusp(0,i,e); //initialize and assign
  }
  for(int j=e+1; j< nelectrons; j++)
  {
    elecPartialSum(e,j)=eecusp(0,e,j);
  }


  //----electron electron correlation------

  for(int seq=0;
      seq< eeCorrelation.GetDim(0);
      seq++)
  {
    doublevar coeff=eeCorrelation(seq)*faccp;
    for(int i=0; i< e; i++)
    {
      elecPartialSum(i,e)+=coeff*bupdate(i)(seq+1);
    }

    for(int j=e+1; j< nelectrons; j++)
    {
      elecPartialSum(e,j)+=coeff*bupdate(j)(seq+1);
    }

  }


  //--------electron electron ion------------

  //Array1 <doublevar> elecPartial_temp(nelectrons);
  elecPartial_temp=0;

  for(int I=0; I< nions; I++)
  {
    //cout << "center " << I << endl;

    doublevar v_klval;
    int ionoffset=I*estride;
    int eoffset;
    //cout << "ionoffset " << ionoffset << endl;
    for(int seq=eeiCorr_start(I);
        seq < eeiCorr_end(I);
        seq++)
    {
      //doublevar total=0;
      int c_kt=c_k(seq);
      int c_lt=c_l(seq);
      int c_mt=c_m(seq);
      doublevar coeff=faccp*eeiCorrelation(seq);

      doublevar tmpa_k=a_kval.v[ionoffset+e].v[c_kt];
      doublevar tmpa_l=a_kval.v[ionoffset+e].v[c_lt];
      for(int i=0; i< e; i++)
      {
        // v_klval=coeff*
        //  (a_k(I,i)(c_k,0)
        //   *a_k(I,e)(c_l,0)
        //   +a_k(I,i)(c_l,0)
        //   *a_k(I,e)(c_k,0) );
        //  elecPartialSum(i,e)+=v_klval*b_m(i,e)(c_m,0);
        //eoffset=i*estride;
        v_klval=coeff*(
                  a_kval.v[ionoffset+i].v[c_kt]
                  *tmpa_l
                  +a_kval.v[ionoffset+i].v[c_lt]
                  *tmpa_k);
        elecPartial_temp.v[i]+=v_klval*bupdate.v[i].v[c_mt];
      }
      //cout << "total-mideei " << total << endl;
      eoffset=e*estride;
      for(int j=e+1; j< nelectrons; j++)
      {
        //v_klval=coeff*(a_k(I,j)(c_k,0)*a_k(I,e)(c_l,0)
        //	 +a_k(I,j)(c_l,0)*a_k(I,e)(c_k,0));
        v_klval=coeff*(
                  a_kval.v[ionoffset+j].v[c_kt]
                  *tmpa_l
                  +a_kval.v[ionoffset+j].v[c_lt]
                  *tmpa_k);
        elecPartial_temp.v[j]+=v_klval*bupdate.v[j].v[c_mt];
      }
      //cout << "total-eei " << total << endl;
    }
  }

  for(int i=0; i< e; i++)
  {
    elecPartialSum(i,e)+=elecPartial_temp(i);
  }
  for(int j=e+1; j< nelectrons; j++)
  {
    elecPartialSum(e,j)+=elecPartial_temp(j);
  }

}

//----------------------------------------------------------------------

void Jastrow_wf_data::makeStaticSave(Sample_point * sample,
                                     Array1 <doublevar> & values,
                                     Array2 <doublevar> & partialValues,
                                     Array3 <doublevar> & derivatives)
{

  //cout << "makeStatic " << endl;
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

      elecElecBasis->calcLap(eedist, b_m(i,j), basisoffset);
      //cout << "cusp " << endl;

      calcEECuspLap(i,j,eedist);
      //cout << "done " << endl;
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
      elecIonBasis->calcLap(eidist, a_k(i,j), basisoffset);
    }
  }

  //now fill the partial sums.
  //cout << "valSize() " << valSize() << endl;
  values.Resize(valSize());
  partialValues.Resize(valSize(), nelectrons);
  derivatives.Resize(valSize(), nelectrons, 5);
  values=0;
  partialValues=0;
  derivatives=0;
  int counter=0;
  for(int I=0; I< nions; I++)
  {
    int k=0;

    for(int seq=eiCorr_start(I);
        seq < eiCorr_end(I); seq++)
    {
      //doublevar total=0;
      int place=seq+counter;

      for(int e=0; e< nelectrons; e++)
      {
        values(place)+=a_k(I,e)(k+1,0);
        partialValues(place, e)+=a_k(I,e)(k+1,0);
        for(int d=1; d< 5; d++)
        {
          derivatives(place,e,d)+=a_k(I,e)(k+1,d);
        }
      }
      k++;
    }
  }

  counter+=eiCorrelation.GetDim(0);

  //cout << "counter " << counter << endl;

  //----------e-e correlation----------

  for(int seq=0;
      seq< eeCorrelation.GetDim(0);
      seq++)
  {

    for(int i=0; i< nelectrons; i++)
    {
      for(int j=i+1; j< nelectrons; j++)
      {
        values(counter)+=b_m(i,j)(seq+1, 0);
        partialValues(counter, i)+=b_m(i,j)(seq+1, 0);
        partialValues(counter, j)+=b_m(i,j)(seq+1, 0);
        for(int d=1; d<4; d++)
        {
          derivatives(counter, i, d)+=b_m(i,j)(seq+1,d);
          derivatives(counter, j, d)-=b_m(i,j)(seq+1,d);
        }
        derivatives(counter, i, 4)+=b_m(i,j)(seq+1,4);
        derivatives(counter, j, 4)+=b_m(i,j)(seq+1,4);
      }
    }
    counter++;
  }

  //cout << "done ee " << endl;
  //------e-e-i correlation---------

  for(int I=0; I< nions; I++)
  {
    for(int seq=eeiCorr_start(I);
        seq < eeiCorr_end(I);
        seq++)
    {
      int place=seq+counter;
      int c_kt=c_k(seq); //temporary variable lookups
      int c_lt=c_l(seq);
      int c_mt=c_m(seq);

      doublevar dotproduct_i, dotproduct_j;

      int ionstride=I*a_k.GetDim(1);

      c_kt*=5;
      c_mt*=5;
      c_lt*=5;
      doublevar v_klval;
      doublevar v_kl_i[5], v_kl_j[5];

      for(int i=0; i< nelectrons; i++)
      {
        int istride=i*nelectrons;
        int ionstridei=ionstride+i;
        doublevar aki_temp=a_k.v[ionstridei].v[c_kt];
        doublevar ali_temp=a_k.v[ionstridei].v[c_lt];
        doublevar akj_temp, alj_temp;
        doublevar aki_drift[5], ali_drift[5];

        for(int d=1; d< 5; d++)
        {
          aki_drift[d]=a_k.v[ionstridei].v[c_kt+d];
          ali_drift[d]=a_k.v[ionstridei].v[c_lt+d];
        }
        doublevar bmij[5];
        for(int j=i+1; j< nelectrons; j++)
        {
          int ionstridej=ionstride+j;
          int istridej=istride+j;

          for(int d=0; d<5; d++)
          {
            bmij[d]=b_m.v[istridej].v[c_mt+d];
          }


          akj_temp=a_k.v[ionstridej].v[c_kt];
          alj_temp=a_k.v[ionstridej].v[c_lt];


          v_klval=(aki_temp*alj_temp+ali_temp*akj_temp);
          //elecPartialSum.v[istridej]+=v_klval*bmij[0];
          values(place)+=v_klval*bmij[0];
          partialValues(place,i)+=v_klval*bmij[0];
          partialValues(place,j)+=v_klval*bmij[0];

          for(int d=1; d< 5; d++)
          {
            //v_kl_i(d)=coeff*(a_k(I,i)(c_kt,d)*a_k(I,j)(c_lt,0)
            //	       +a_k(I,i)(c_lt,d)*a_k(I,j)(c_kt,0));
            //v_kl_j(d)=coeff*(a_k(I,i)(c_kt,0)*a_k(I,j)(c_lt,d)
            //	       +a_k(I,i)(c_lt,0)*a_k(I,j)(c_kt,d));
            //driftij(i,j,d)+=v_klval*b_m(i,j)(c_mt,d);
            //driftion_i(i,j,d)+=v_kl_i(d)*b_m(i,j)(c_mt, 0);
            //driftion_j(i,j,d)+=v_kl_j(d)*b_m(i,j)(c_mt, 0);
            //cout << "d " << d << endl;

            v_kl_i[d]=(aki_drift[d]*alj_temp
                       +ali_drift[d]*akj_temp);

            v_kl_j[d]=(aki_temp*a_k.v[ionstridej].v[c_lt+d]
                       +ali_temp*a_k.v[ionstridej].v[c_kt+d]);
          }



          for(int d=1; d< 4; d++)
          {
            derivatives(place, i,d) += v_klval*bmij[d];
            derivatives(place, j,d) -= v_klval*bmij[d];
            derivatives(place, i,d) += v_kl_i[d]*bmij[0];
            derivatives(place, j,d) += v_kl_j[d]*bmij[0];

          }

          derivatives(place, i,4) += v_klval*bmij[4];
          derivatives(place, j,4) += v_klval*bmij[4];
          derivatives(place, i,4) += v_kl_i[4]*bmij[0];
          derivatives(place, j,4) += v_kl_j[4]*bmij[0];
          //cout << "done " << endl;

          //Cross terms in the laplacian
          dotproduct_i=0;
          dotproduct_j=0;
          for(int d=1; d<4; d++)
          {
            //dotproduct_i +=v_kl_i(d)*b_m(i,j)(c_mt,d);
            //dotproduct_j +=v_kl_j(d)*b_m(i,j)(c_mt,d);
            dotproduct_i +=v_kl_i[d]*bmij[d];//*b_m.v[istride+j].v[c_mt+d];
            dotproduct_j +=v_kl_j[d]*bmij[d];//b_m.v[istride+j].v[c_mt+d];

          }
          //since we don't change the
          //sign for the lap later,
          //we need to do it here.
          //driftion_i(i,j,4)+=2.0*dotproduct_i;
          //driftion_j(i,j,4)-=2.0*dotproduct_j;
          //int driftoffset=(istride+j)*5+4;
          //driftion_i.v[driftoffset]+=2.0*dotproduct_i;
          //driftion_j.v[driftoffset]-=2.0*dotproduct_j;
          derivatives(place, i, 4)+=2.0*dotproduct_i;
          derivatives(place, j, 4)-=2.0*dotproduct_j;
        }
      }
    }
  }
  counter += eeiCorrelation.GetDim(0);

  //cout << "done " << endl;
}

//----------------------------------------------------------------------


void Jastrow_wf_data::updateParms(Sample_point * sample,
                                  Array1 <doublevar> & oldvalue,
                                  Array2 <doublevar> & oldValPartial,
                                  Array3 <doublevar> & oldderiv,
                                  doublevar & value,
                                  Array1 <doublevar> & valPartialSum,
                                  Array2 <doublevar> & derivatives)
{
  //cout << "updateParms" << endl;
  value=0;
  valPartialSum=0;
  derivatives=0;
  Array1 <doublevar> eedist(5);

  Array1 <doublevar> temp(5);
  sample->updateEIDist();
  sample->updateEEDist();

  for(int i=0; i< nelectrons; i++)
  {
    for(int j=i+1; j< nelectrons; j++)
    {
      sample->getEEDist(i,j, eedist);
      calcEECuspLap(i,j,eedist);
    }
  }

  for(int i=0; i< nelectrons; i++)
  {
    for(int j=i+1; j< nelectrons; j++)
    {
      value+=eecusp(0,i,j);
      valPartialSum(i)+=eecusp(0,i,j);
      valPartialSum(j)+=eecusp(0,i,j);
    }
  }
  for(int d=1; d< 4; d++)
  {
    for(int i=0; i< nelectrons; i++)
    {
      for(int j=i+1; j< nelectrons; j++)
      {
        derivatives(i,d)+=eecusp(d,i,j);
        derivatives(j,d)-=eecusp(d,i,j);
      }
    }
  }
  for(int i=0; i< nelectrons; i++)
  {
    for(int j=i+1; j< nelectrons; j++)
    {
      derivatives(i,4)+=eecusp(4,i,j);
      derivatives(j,4)+=eecusp(4,i,j);
    }
  }



  int counter=0;
  for(int seq=0; seq < eiCorrelation.GetDim(0); seq++)
  {
    doublevar coeff=facco*eiCorrelation(seq);
    int place=seq+counter;
    //cout << "coeff " << coeff << endl;
    //cout << "counter " << counter << "  oldval " << oldvalue(counter)
    //   << endl;
    value+=coeff*oldvalue(place);
    for(int e=0; e< nelectrons; e++)
    {
      valPartialSum(e)+=coeff*oldValPartial(place,e);
      for(int d=1; d< 5; d++)
      {
        derivatives(e, d)+=coeff*oldderiv(place,e,d);
        //cout << "deriv " << e << "   " << d << "  "
        //   << derivatives(e,d) << endl;
      }
    }
  }

  counter+=eiCorrelation.GetDim(0);

  //cout << "eecorr" << endl;
  //----electron electron correlation------

  for(int seq=0;
      seq< eeCorrelation.GetDim(0);
      seq++)
  {
    doublevar coeff=eeCorrelation(seq)*faccp;
    value+=coeff*oldvalue(counter);
    for(int e=0; e< nelectrons; e++)
    {
      valPartialSum(e)+=coeff*oldValPartial(counter,e);
      for(int d=1; d< 5; d++)
      {
        derivatives(e, d)+=coeff*oldderiv(counter, e,d);
      }
    }
    counter++;
  }



  for(int seq=0; seq < eeiCorrelation.GetDim(0); seq++)
  {
    doublevar coeff=faccp*eeiCorrelation(seq);
    int place=seq+counter;

    value+=coeff*oldvalue(place);
    for(int e=0; e< nelectrons; e++)
    {
      valPartialSum(e)+=coeff*oldValPartial(place,e);
      for(int d=1; d< 5; d++)
      {
        derivatives(e, d)+=coeff*oldderiv(place, e,d);
      }
    }
  }
  counter+=eeiCorrelation.GetDim(0);

  //cout << "eeiCorr" << endl;
  //for(int i=0; i< nelectrons; i++) {
  //for(int d=1; d< 5; d++) {
  // cout << "derivatives(" << i << "," << d
  //   <<") " << derivatives(i,d) << endl;
  //}
  //}
  //cout << "done updateParms " << endl;
}

//----------------------------------------------------------------------


/*!
 
 */
void Jastrow_wf_data::fillParmInd(Sample_point * sample,
                                  int e,
                                  Array1 <doublevar> & vals)
{

  //assert(vals.GetDim(0) >= totalsize);
  //cout << "fillparmInd" << endl;
  Array1 <doublevar> eedist(5);
  Array1 <doublevar> temp(5);


  Array1 <doublevar> eiCorr_save(eiCorrelation.GetDim(0));
  Array1 <doublevar> eeCorr_save(eeCorrelation.GetDim(0));
  Array1 <doublevar> eeiCorr_save(eeiCorrelation.GetDim(0));
  eiCorr_save=0;
  eeCorr_save=0;
  eeiCorr_save=0;


  sample->updateEIDist();
  sample->updateEEDist();


  //Calculate the updated basis functions for electron-electron
  for(int i=0; i< e; i++)
  {
    sample->getEEDist(i,e,eedist);
    elecElecBasis->calcVal(eedist, bupdate(i), basisoffset);
  }
  for(int j=e+1; j< nelectrons; j++)
  {
    sample->getEEDist(e,j,eedist);
    elecElecBasis->calcVal(eedist,  bupdate(j), basisoffset);
  }

  //Calculate the updated basis functions for electron-ion
  Array1 <doublevar> eidist(5);
  for(int i=0; i< nions; i++)
  {
    for(int j=0; j< nelectrons; j++)
    {
      sample->getEIDist(j,i,eidist);
      elecIonBasis->calcLap(eidist, a_k(i,j), basisoffset);
    }
  }


  //Fill the saved array

  //---------electron-ion----------

  for(int I=0; I< nions; I++)
  {
    int k=0;
    for(int seq=eiCorr_start(I);
        seq < eiCorr_end(I); seq++)
    {
      eiCorr_save(seq)+=a_k(I,e)(k+1,0);
      k++;
    }
  }


  //----electron electron correlation------
  for(int seq=0;
      seq< eeCorrelation.GetDim(0);
      seq++)
  {
    doublevar total=0;
    for(int i=0; i< e; i++)
    {
      total+=bupdate(i)(seq+1);
    }

    for(int j=e+1; j< nelectrons; j++)
    {
      total+=bupdate(j)(seq+1);
    }
    eeCorr_save(seq)+=total;
    //cout << "total-fill " << total << endl;
  }

  int estride=a_k.GetDim(1);

  //--------electron electron ion------------

  for(int I=0; I< nions; I++)
  {
    //cout << "center " << I << endl;

    doublevar v_klval;
    int ionoffset=I*estride;
    int eoffset;
    //cout << "ionoffset " << ionoffset << endl;
    for(int seq=eeiCorr_start(I);
        seq < eeiCorr_end(I);
        seq++)
    {
      doublevar total=0;
      int c_kt=c_k(seq);
      int c_lt=c_l(seq);
      int c_mt=c_m(seq);
      c_kt*=5;
      c_lt*=5; //account for stride.

      doublevar tmpa_k=a_k.v[ionoffset+e].v[c_kt];
      doublevar tmpa_l=a_k.v[ionoffset+e].v[c_lt];
      for(int i=0; i< e; i++)
      {
        // v_klval=coeff*
        //  (a_k(I,i)(c_k,0)
        //   *a_k(I,e)(c_l,0)
        //   +a_k(I,i)(c_l,0)
        //   *a_k(I,e)(c_k,0) );
        //  elecPartialSum(i,e)+=v_klval*b_m(i,e)(c_m,0);
        v_klval=(
                  a_k.v[ionoffset+i].v[c_kt]
                  *tmpa_l
                  +a_k.v[ionoffset+i].v[c_lt]
                  *tmpa_k);
        total+=v_klval*bupdate.v[i].v[c_mt];
      }

      eoffset=e*estride;
      for(int j=e+1; j< nelectrons; j++)
      {
        //v_klval=coeff*(a_k(I,j)(c_k,0)*a_k(I,e)(c_l,0)
        //	 +a_k(I,j)(c_l,0)*a_k(I,e)(c_k,0));
        //elecPartialSum(e,j)+=v_klval*b_m(e,j)(c_m,0);
        v_klval=(
                  a_k.v[ionoffset+j].v[c_kt]
                  *tmpa_l
                  +a_k.v[ionoffset+j].v[c_lt]
                  *tmpa_k);
        total+=v_klval*bupdate.v[j].v[c_mt];
      }
      eeiCorr_save(seq)+=total;
    }
  }

  int counter=0;
  for(int i=0; i< eiCorr_save.GetDim(0); i++)
  {
    vals(counter)=eiCorr_save(i);
    counter++;
  }
  for(int i=0; i< eeCorr_save.GetDim(0); i++)
  {
    vals(counter)=eeCorr_save(i);
    counter++;
  }
  for(int i=0; i< eeiCorr_save.GetDim(0); i++)
  {
    vals(counter)=eeiCorr_save(i);
    counter++;
  }

  //for(int i=0; i< counter; i++) {
  //  cout << "vals " << i << "   " << vals(i) << endl;
  //}
  //cout << "done " << endl;

}

//----------------------------------------------------------------------


void Jastrow_wf_data::updateParmsVal(Sample_point * sample,
                                     int e, doublevar & value,
                                     Array1 <doublevar> & oldvalue
                                    )
{
  //cout << "updateParmsVal " << endl;
  Array1 <doublevar> eedist(5);

  sample->updateEEDist();
  value=0;

  //Calculate the updated basis functions for electron-electron
  for(int i=0; i< e; i++)
  {
    sample->getEEDist(i,e,eedist);
    calcEECuspVal( i,e,eedist);
  }
  for(int j=e+1; j< nelectrons; j++)
  {
    sample->getEEDist(e,j,eedist);
    calcEECuspVal(e,j,eedist);

  }
  int counter=0;
  for(int seq=0; seq < eiCorrelation.GetDim(0); seq++)
  {
    doublevar coeff=facco*eiCorrelation(seq);
    value+=coeff*oldvalue(counter);
    counter++;
  }
  //cout << "e-ion " << value << endl;
  //Electron-electron interaction
  //-----electron electron cusp----------
  for(int i=0; i< e; i++)
  {
    value+= eecusp(0, i,e);
  }
  for(int j=e+1; j< nelectrons; j++)
  {
    value+= eecusp(0,e,j);
  }

  //cout << "eecusp " << value << endl;
  //----electron electron correlation------
  for(int seq=0;
      seq< eeCorrelation.GetDim(0);
      seq++)
  {
    doublevar coeff=eeCorrelation(seq)*faccp;
    value+=coeff*oldvalue(counter);
    counter++;
  }

  //cout << "eeCorr " << value << endl;

  for(int seq=0; seq < eeiCorrelation.GetDim(0); seq++)
  {
    doublevar coeff=faccp*eeiCorrelation(seq);
    value+=coeff*oldvalue(counter);
    counter++;
  }

  //cout << "done " << endl;
  //cout << "eeiCorr " << value << endl;
}

//----------------------------------------------------------------------
