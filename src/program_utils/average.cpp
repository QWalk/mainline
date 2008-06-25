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

#include "average.h"
#include "qmc_io.h"




/*
Averaging without weights.
 */
void average(int start, int end,
             Array2 <doublevar> & vals,
             Array1 <doublevar> & average,
             Array1 <doublevar> & variance,
             int paravg
            )
{

  int navg=vals.GetDim(0);  //number of averaging variables
  int npoints=end-start;//vals.GetDim(1); //Number of points

  assert(vals.GetDim(1) >= npoints);

  average.Resize(navg);
  variance.Resize(navg);

  average=0;
  variance=0;
  int totpoints=0;


#ifdef USE_MPI

  if(paravg)
  {
    MPI_Allreduce(&npoints, &totpoints, 1,
                  MPI_INT, MPI_SUM, MPI_Comm_grp);
  }
  else
  {
    totpoints=npoints;
  }
#else
  totpoints=npoints;
#endif
  //----------------------------------------
  //take average
  for(int n=0; n< navg; n++)
  {
    doublevar average_temp=0;
    for(int point=start; point < end; point++)
    {
      average_temp+=vals(n,point);
    }

#ifdef USE_MPI
    if(paravg)
    {
      MPI_Allreduce(&average_temp, &(average(n)), 1,
                    MPI_DOUBLE, MPI_SUM, MPI_Comm_grp);
    }
    else
    {
      average(n)=average_temp;
    }
#else
    average(n)=average_temp;
#endif

    average(n)/=totpoints;
  }
  //----------------------------------------
  //get variance
  for(int n=0; n< navg; n++)
  {
    doublevar variance_temp=0;
    for(int point=start; point < end; point++)
    {
      variance_temp+=
        (vals(n,point)-average(n))
        *(vals(n,point)-average(n));
    }

#ifdef USE_MPI
    if(paravg)
    {
      MPI_Allreduce(&variance_temp, &(variance(n)), 1,
                    MPI_DOUBLE, MPI_SUM, MPI_Comm_grp);
    }
    else
    {
      variance(n)=variance_temp;
    }
#else
    variance(n)=variance_temp;
#endif

    variance(n)/=totpoints;
  }
}

//----------------------------------------------------------------------
void average(int start, int end,
             Array1 <doublevar> & vals, Array1 <doublevar> & weights,
             doublevar & average, doublevar & variance,
             int paravg
            )
{

  //int npoints=end-start;

  assert(end-start<= weights.GetDim(0));
  assert(end-start <= vals.GetDim(0));


  average=0;
  variance=0;

  doublevar weightsum; //sum of weights
  doublevar negweightsum; //sum of negative weights
  negweightsum=0;
  weightsum=0;
  //take average
  doublevar weight_temp=0;
  doublevar avg_temp=0;
  doublevar negweight_temp=0;
  for(int point=start; point < end; point++)
  {
    weight_temp+=weights(point);
    avg_temp+=vals(point)*weights(point);
    if(weights(point) < 0)
    {
      negweight_temp+=weights(point);
    }
  }
#ifdef USE_MPI
  if(paravg)
  {

    MPI_Allreduce(&weight_temp, &(weightsum), 1,
                  MPI_DOUBLE, MPI_SUM, MPI_Comm_grp);
    MPI_Allreduce(&avg_temp, &(average), 1,
                  MPI_DOUBLE, MPI_SUM, MPI_Comm_grp);
    MPI_Allreduce(&negweight_temp, &(negweightsum), 1,
                  MPI_DOUBLE, MPI_SUM, MPI_Comm_grp);

    average/=weightsum;
  }
  else
  {
    negweightsum=negweight_temp;
    weightsum=weight_temp;
    average=avg_temp/weight_temp;
  }
#else
  negweightsum=negweight_temp;
  weightsum=weight_temp;
  average=avg_temp/weight_temp;
#endif
  //cout << "weightsum(" << n << ") " << weightsum(n) << endl;
  //average(n)/=weightsum(n);
  if(mpi_info.node ==0 )
  {
    if( negweightsum/weightsum > 1e-10) 
    debug_write(cout, "negative proportion of weight : ",
                fabs(negweightsum/weightsum), "\n");
    if(fabs(negweightsum/weightsum) > .5) {
      cout << "Warning!  Negative proportion of weight is "
      << fabs(negweightsum/weightsum) << endl;
    }
    if(weightsum==0)
      cout << "WARNING: sum of weights for is zero!" << endl;
  }


  //get variance
  doublevar var_temp=0;
  for(int point=start; point < end; point++)
  {
    var_temp+=
      (vals(point)-average)
      *(vals(point)-average)*weights(point);
  }
#ifdef USE_MPI
  if(paravg)
  {
    MPI_Allreduce(&var_temp, &(variance), 1,
                  MPI_DOUBLE, MPI_SUM, MPI_Comm_grp);
    //cout << mpi_info.node << " var_temp " << var_temp << endl;
    //cout << mpi_info.node << " variance " << variance(n) << endl;
  }
  else
  {
    variance=var_temp;
  }
#else
  variance=var_temp;
#endif
  //cout << "variance(" << n << ") " << variance(n) << endl;
  variance/=weightsum;

}

//----------------------------------------------------------------------

/*
Averaging with weights
 */
void average(int start, int end,
             Array2 <doublevar> & vals, Array2 <doublevar> & weights,
             Array1 <doublevar> & average, Array1 <doublevar> & variance,
             int paravg
            )
{
  int navg=vals.GetDim(0);  //number of averaging variables
  //int npoints=end-start;
  assert(navg==weights.GetDim(0));

  assert(end-start<= weights.GetDim(1));
  assert(end-start <= vals.GetDim(1));

  average.Resize(navg);
  variance.Resize(navg);
  average=0;
  variance=0;

  Array1 <doublevar> weightsum(navg); //sum of weights
  Array1 <doublevar> negweightsum(navg); //sum of negative weights
  negweightsum=0;
  weightsum=0;
  //take average
  for(int n=0; n< navg; n++)
  {
    doublevar weight_temp=0;
    doublevar avg_temp=0;
    doublevar negweight_temp=0;
    for(int point=start; point < end; point++)
    {
      weight_temp+=weights(n,point);
      avg_temp+=vals(n,point)*weights(n,point);
      if(weights(n,point) < 0)
      {
        negweight_temp+=weights(n,point);
      }
    }
#ifdef USE_MPI
    if(paravg)
    {

      MPI_Allreduce(&weight_temp, &(weightsum(n)), 1,
                    MPI_DOUBLE, MPI_SUM, MPI_Comm_grp);
      MPI_Allreduce(&avg_temp, &(average(n)), 1,
                    MPI_DOUBLE, MPI_SUM, MPI_Comm_grp);
      MPI_Allreduce(&negweight_temp, &(negweightsum(n)), 1,
                    MPI_DOUBLE, MPI_SUM, MPI_Comm_grp);

      average(n)/=weightsum(n);
    }
    else
    {
      negweightsum(n)=negweight_temp;
      weightsum(n)=weight_temp;
      average(n)=avg_temp/weight_temp;
    }
#else
    negweightsum(n)=negweight_temp;
    weightsum(n)=weight_temp;
    average(n)=avg_temp/weight_temp;
#endif
    //cout << "weightsum(" << n << ") " << weightsum(n) << endl;
    //average(n)/=weightsum(n);
    if(mpi_info.node ==0 )
    {
      doublevar negprop=fabs(negweightsum(n)/weightsum(n));
      if(negprop > 1e-10) {
        debug_write(cout , "negative proportion of weight for " , n , " : ");
        debug_write(cout, negprop, "\n");
      }
      if(negprop > .5) {
        cout << "Warning! Negative proportion of weight for " << n
        << " is " << negprop << endl;
      }
      if(weightsum(n)==0)
        cout << "WARNING: sum of weights for " << n << " is zero!" << endl;
    }
  }

  //get variance

  for(int n=0; n< navg; n++)
  {
    doublevar var_temp=0;
    for(int point=start; point < end; point++)
    {
      var_temp+=
        (vals(n,point)-average(n))
        *(vals(n,point)-average(n))*weights(n,point);
    }
#ifdef USE_MPI
    if(paravg)
    {
      MPI_Allreduce(&var_temp, &(variance(n)), 1,
                    MPI_DOUBLE, MPI_SUM, MPI_Comm_grp);
      //cout << mpi_info.node << " var_temp " << var_temp << endl;
      //cout << mpi_info.node << " variance " << variance(n) << endl;
    }
    else
    {
      variance(n)=var_temp;
    }
#else
    variance(n)=var_temp;
#endif
    //cout << "variance(" << n << ") " << variance(n) << endl;
    variance(n)/=weightsum(n);

  }

}



//----------------------------------------------------------------------

void autocorrelation(int start,
		     int end,
		     const Array2 <doublevar> & vals, //!< values
		     const Array1 <doublevar> & avg,
		     const Array1 <doublevar> & var, //!< variance
		     Array2 <doublevar> & autocorr, //!< autocorrelation in (wf,n-1 steps)
		     int depth,
		     int paravg //!< sync across processes
		     ) {
  int nwf=vals.GetDim(0);
  assert(end-start >0);
  assert(end-start <= vals.GetDim(1));
  //int depth=4;
  autocorr.Resize(nwf, depth);
  autocorr=0;
  for(int w=0; w< nwf; w++) {
    for(int d=1; d< depth+1; d++) {
      int npoints=end-start-d;

      for(int i=start; i< end-d; i++) {
	autocorr(w,d-1)+=(vals(w,i)-avg(w))*(vals(w,i+d)-avg(w));
      }

      if(paravg) {
        npoints=parallel_sum(npoints);
        autocorr(w,d-1)=parallel_sum(autocorr(w,d-1));
      }
      autocorr(w,d-1)/=npoints;

    }
        
    for(int d=0; d< depth; d++) {
      
      autocorr(w,d)/= sqrt(var(w)*var(w));
    }
  }
}

//----------------------------------------------------------------------
