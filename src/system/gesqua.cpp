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

#include "Array.h"
#include "Qmc_std.h"

/*!
For these points, see
Mitas, Shirley, Ceperly, J. Chem. Phys.  95, p 3467  (1991)
 
\todo add the possibility for more points
 
*/

void gesqua (int & nq,
             Array1<doublevar> &xq, //x points
             Array1<doublevar> &yq, //y points
             Array1<doublevar> &zq, //z points
             Array1<doublevar> &wq) //weights
{
  if( nq < 4)
  {
    error("In gesqua, AIP must be greater than or equal to 4.");
  }
  else if(nq >= 4 && nq < 6)
  {
    nq=4;
    for(int i=0; i< nq; i++)
    {
      wq(i)=0.25;
    }
    doublevar q=1./sqrt(3.);
    xq(0)=   q;
    xq(1)=   q;
    xq(2)=  -q;
    xq(3)=  -q;

    yq(0)=   q;
    yq(1)=  -q;
    yq(2)=   q;
    yq(3)=  -q;

    zq(0)=   q;
    zq(1)=  -q;
    zq(2)=  -q;
    zq(3)=   q;
  }

  else if(nq >=6 && nq < 12)
  {
    nq=6;
    for(int i=0; i< nq; i++)
    {
      wq(i)=1./6.;
    }

    xq=yq=zq=0;

    xq(0)=  1;
    xq(1)= -1;

    yq(2)=  1;
    yq(3)= -1;

    zq(4)=  1;
    zq(5)= -1;
  }

  
  else if( nq >= 12 && nq < 18)
  {
    nq=12;
    for(int i=0; i< nq; i++)
    {
      wq(i)=1./12.;
    }

    doublevar fi0=pi/5.;
    doublevar cstha=1./sqrt(5.0);
    doublevar sntha=2./sqrt(5.0);
    xq=yq=zq=0;

    zq(0)=1;
    zq(1)=-1;

    for(int i=0; i< 5; i++)
    {
      doublevar rk2=2*i;
      doublevar crk2=cos(rk2*fi0);
      doublevar srk2=sin(rk2*fi0);

      xq(i+2)=sntha*crk2;
      yq(i+2)=sntha*srk2;
      zq(i+2)=cstha;

      doublevar crk21=cos((rk2+1)*fi0);
      doublevar srk21=sin((rk2+1)*fi0);
      xq(i+7)=sntha*crk21;
      yq(i+7)=sntha*srk21;
      zq(i+7)=-cstha;
      
    }
  }
  else
  {
    error("Bad number of AIP's: ", nq);
  }
}

//------------------------------------------------------------------------
