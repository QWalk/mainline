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
  else if( nq >= 18 && nq < 26)
  {
    //see Handbook of mathematical functions, Abramowitz, pg. 894
    nq=18;
    for(int i=0; i< nq; i++)
    {
      if(i<6)
	wq(i)=1./30.;
      else
	wq(i)=1./15.;
    }
    xq=yq=zq=0;
    xq(0)=  1;
    xq(1)= -1;

    yq(2)=  1;
    yq(3)= -1;

    zq(4)=  1;
    zq(5)= -1;

    doublevar cstha=1./sqrt(2.);
    xq(6)=cstha;
    yq(6)=cstha;

    xq(7)=-cstha;
    yq(7)=cstha;

    xq(8)=cstha;
    yq(8)=-cstha;

    xq(9)=-cstha;
    yq(9)=-cstha;

    xq(10)=cstha;
    zq(10)=cstha;

    xq(11)=-cstha;
    zq(11)=cstha;

    xq(12)=cstha;
    zq(12)=-cstha;

    xq(13)=-cstha;
    zq(13)=-cstha;

    zq(14)=cstha;
    yq(14)=cstha;

    zq(15)=-cstha;
    yq(15)=cstha;

    zq(16)=cstha;
    yq(16)=-cstha;

    zq(17)=-cstha;
    yq(17)=-cstha;
  }
  else if( nq >= 26 && nq < 32 )
  {
    for(int i=0; i< nq; i++)
    {
      if(i<6)
	wq(i)=40./840.;
      else if(i<18)
	wq(i)=32./840.;
      else
	wq(i)=27./840.;
    }
    xq=yq=zq=0;
    xq(0)=  1;
    xq(1)= -1;

    yq(2)=  1;
    yq(3)= -1;

    zq(4)=  1;
    zq(5)= -1;

    doublevar cstha=1./sqrt(2.);
    xq(6)=cstha;
    yq(6)=cstha;

    xq(7)=-cstha;
    yq(7)=cstha;

    xq(8)=cstha;
    yq(8)=-cstha;

    xq(9)=-cstha;
    yq(9)=-cstha;

    xq(10)=cstha;
    zq(10)=cstha;

    xq(11)=-cstha;
    zq(11)=cstha;

    xq(12)=cstha;
    zq(12)=-cstha;

    xq(13)=-cstha;
    zq(13)=-cstha;

    zq(14)=cstha;
    yq(14)=cstha;

    zq(15)=-cstha;
    yq(15)=cstha;

    zq(16)=cstha;
    yq(16)=-cstha;

    zq(17)=-cstha;
    yq(17)=-cstha;

    doublevar csthb=1./sqrt(3.);
    xq(18)=csthb;
    yq(18)=csthb;
    zq(18)=csthb;

    xq(19)=csthb;
    yq(19)=-csthb;
    zq(19)=csthb;

    xq(20)=csthb;
    yq(20)=csthb;
    zq(20)=-csthb;

    xq(21)=csthb;
    yq(21)=-csthb;
    zq(21)=-csthb;

    xq(22)=-csthb;
    yq(22)=csthb;
    zq(22)=csthb;

    xq(23)=-csthb;
    yq(23)=-csthb;
    zq(23)=csthb;

    xq(24)=-csthb;
    yq(24)=csthb;
    zq(24)=-csthb;

    xq(25)=-csthb;
    yq(25)=-csthb;
    zq(25)=-csthb;
  }
  else if( nq == 32 ) {
    // icosahedron 32 point rule
    for (int i=0; i<12; i++) wq(i)=5.0/168;
    for (int i=12; i<32; i++) wq(i)=27.0/840;
    xq(0)=0.0;
    yq(0)=0.0;
    zq(0)=1.0;

    xq(1)=0.0;
    yq(1)=0.0;
    zq(1)=-1.0;
    double th1=acos((2+sqrt(5.0))/sqrt(15+6*sqrt(5.0)));
    double th2=acos(1.0/sqrt(15+6*sqrt(5.0)));
    for (int i=0; i<5; i++) {
      doublevar th=atan(2.0);
      doublevar phi=2*i*pi/5;
      xq(i+2)=sin(th)*cos(phi);
      yq(i+2)=sin(th)*sin(phi);
      zq(i+2)=cos(th);

      th=pi-atan(2.0);
      phi=(2*i+1)*pi/5;
      xq(i+7)=sin(th)*cos(phi);
      yq(i+7)=sin(th)*sin(phi);
      zq(i+7)=cos(th);     

      th=th1;
      phi=(2*i+1)*pi/5;
      xq(i+12)=sin(th)*cos(phi);
      yq(i+12)=sin(th)*sin(phi);
      zq(i+12)=cos(th);

      th=th2;
      phi=(2*i+1)*pi/5;
      xq(i+17)=sin(th)*cos(phi);
      yq(i+17)=sin(th)*sin(phi);
      zq(i+17)=cos(th);

      th=pi-th1;
      phi=2*i*pi/5;
      xq(i+22)=sin(th)*cos(phi);
      yq(i+22)=sin(th)*sin(phi);
      zq(i+22)=cos(th);

      th=pi-th2;
      phi=2*i*pi/5;
      xq(i+27)=sin(th)*cos(phi);
      yq(i+27)=sin(th)*sin(phi);
      zq(i+27)=cos(th);
    }
  }
  else
  {
    error("Bad number of AIP's: ", nq);
  }
}

//------------------------------------------------------------------------
