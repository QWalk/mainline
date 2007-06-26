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

Random_generator rng;


double unif()
{
  static int ix=1234567;
  int k1=ix/127773;
  ix=16807*(ix-k1*127773)-k1*2836;
  if(ix < 0)
    ix+=2147483647;
  return 4.656612875e-10 * ix;
}


doublevar ranr2exponential() {
  while(1) { 
    doublevar ex=-log(rng.ulec())/.4;
    //cout << "ex " << ex  << "  prob " << (.5/.8)*ex*ex*exp(-ex+0.4*ex) << endl;
    if(rng.ulec()< (.5/.8)*ex*ex*exp(-ex+0.4*ex) ) { 
      //cout << "acc " << endl;
       return ex;
    }
  }
}

doublevar rancos() {
  return acos(1-2*rng.ulec());
}

//----------------------------------------------------------------------

void generate_random_rotation(Array1 <doublevar> & x,
                              Array1 <doublevar> & y,
                              Array1 <doublevar> & z) {
  doublevar usum, u2,uu,yu2,yu1, zz1, yy3, zz2, u1, zz3, yu3, zu2,zu1, zu3,
  cfi, sthet, sfi, yy2, yy1, xsum,  theta,  xsum2;
  assert(x.GetDim(0)==3);
  assert(y.GetDim(0)==3);
  assert(z.GetDim(0)==3);

  do
  {
    x(0)=1.-2.*unif();
    x(1)=1.-2.*unif();
    xsum=(x(0)*x(0)+x(1)*x(1));
  }
  while (xsum >= 1.);

  //cout << "xsum  " << xsum << endl;

  xsum2=2.*sqrt(fabs(1.-xsum));
  x(0)=x(0)*xsum2;
  x(1)=x(1)*xsum2;
  x(2)=1.-2.*xsum;

  theta=acos(x(2));
  sthet=sin(theta);
  if (sthet < 1.e-05)
  {
    cfi=1.;
    sfi=0.;
    x(0)=0;
    x(1)=0.;
    x(2)=1;
    sthet=0.;
  }
  else
  {
    cfi=x(0)/sthet;
    sfi=x(1)/sthet;
  }

  yy1=x(2)*cfi;
  yy2=x(2)*sfi;
  yy3=-sthet;
  zz1=-sfi;
  zz2=cfi;
  zz3=0.;
  do
  {
    u1=unif()*2.-1.;
    u2=unif()*2.-1.;
    usum=u1*u1+u2*u2;
  }
  while(usum >= 1.0);

  //cout << "usum " << usum << endl;

  uu=sqrt(usum);
  u1=u1/uu;
  u2=u2/uu;
  yu1=yy1*u1;
  yu2=yy2*u1;
  yu3=yy3*u1;
  zu1=zz1*u2;
  zu2=zz2*u2;
  zu3=zz3*u2;
  y(0)=yy1*u1+zz1*u2;
  y(1)=yy2*u1+zz2*u2;
  y(2)=yy3*u1+zz3*u2;
  z(0)=yy1*u2-zz1*u1;
  z(1)=yy2*u2-zz2*u1;
  z(2)=yy3*u2-zz3*u1;


}

//----------------------------------------------------------------------
