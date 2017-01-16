/*
 
Copyright (C) 2007 Michal Bajdich

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

#ifndef INTEGRALS_H_INCLUDED
#define INTEGRALS_H_INCLUDED

#include "Array.h"
#include <iostream>

using namespace std;



inline doublevar Hermite(int n, doublevar x)
{

  
  switch(n) {
  case 0:
    return 1.0;
  case 1:
    return 2.0*x;
  case 2:
    return 4.0*x*x-2.0;
  case 3:
    return 8.0*x*x*x-12.0*x;
  case 4:
    return x*x*(16.0*x*x-48.0)+12;
  default:
    error("order greater than 4 or less than zero in Hermite");
    return 0;
  }

}

inline doublevar Integral1D(doublevar G, int n, doublevar alpha)
  //only real part
{

  doublevar ret= 0.88622692545275801365*Hermite(n,G/(2.0*sqrt(alpha)))*exp(-G*G/(4.0*alpha))
    /pow(sqrt(alpha),n+1);
  //0.88622692545=sqrt(pi)/2
  if(n==0) return 2.0*ret;
  if(n==2) return 0.5*ret;
  return ret;


}


//------------------------------------------------------------------



inline doublevar pwIntegral1D_f(doublevar G, int n,
                               doublevar alpha,
                               doublevar sqrtalpha) {
  const doublevar sqrtpi=1.77245385091;
  doublevar ex=sqrtpi;//*exp(-G*G/(4.0*alpha));
  doublevar alfin=1.0/sqrtalpha;
  switch(n) {
    case 0:
      return ex*alfin;
    case 1:
      return G*0.5*ex*alfin/alpha;
    case 2:
      return 0.5*ex*alfin/alpha * (G*G/(2.0*alpha) - 1);
    case 3:
      return G*ex*alfin/(alpha*alpha) * (0.75*G - 0.125*G*G/alpha); //check this one
    default:
      error("pwIntegral1D_f doesn't support n higher than 3");
      return 0;
  }

}



//--------------------------------------------------------------------

inline doublevar P_x_c(Array1 <int> & l, Array1 <doublevar> & center,
		       doublevar alpha, Array2 <doublevar> & coef, 
		       Array2 <doublevar> & G, 
		       Array1 <doublevar> & sin_arr, Array1 <doublevar> & cos_arr)
{
  // Gives the specified gaussian value of the element P_matrix times c_vector.
  // Later has to be normalized by S11.
  // Array coef(i,j) has i=0..M-1 and j=0,1 (0 for real part, 1 for imaginary)
  // and G(i,j) i=0..M-1, j=0,1,2 for x,y,z.
  // where M is number of independent G vectors with respect to G=-G symmetry including zero!
  // Gaussian, as everywhere else is defined by 3 integers l(3) for powers of x,y,z,
  // exponent alpha and vector center(3).


  //Array1 <doublevar> temp_I(3),scalarproduct(G.GetDim(0));
  doublevar scalarproduct;
  doublevar tmp, A,B,A0;
  int sum;
  sum=0;
  scalarproduct=0.0;

  for(int i=0;i<3;i++)
    sum+=l(i);
  doublevar sqrtalpha=sqrt(alpha);

  A0=coef(0,0)*Integral1D(0.0,l(0),alpha)
              *Integral1D(0.0,l(1),alpha)
              *Integral1D(0.0,l(2),alpha);
  A=B=0.0;
  int ngvec=G.GetDim(0);

  //doublevar t_cos, t_sin;
  int i_scal;

  if(sum %2 == 0 ) {
    for(int i=1; i < ngvec; i++) {
      //tmp= pwIntegral1D_f(G(i,0),l(0),alpha, sqrtalpha)
      //    *pwIntegral1D_f(G(i,1),l(1),alpha, sqrtalpha)
      //    *pwIntegral1D_f(G(i,2),l(2),alpha, sqrtalpha);
      //A+=(coef(i,0)*cos_arr(i)-coef(i,1)*sin_arr(i))*tmp;
      i_scal=i*3;
      tmp= pwIntegral1D_f(G.v[i_scal],l.v[0],alpha, sqrtalpha)
          *pwIntegral1D_f(G.v[i_scal+1],l.v[1],alpha, sqrtalpha)
          *pwIntegral1D_f(G.v[i_scal+2],l.v[2],alpha, sqrtalpha)
          *exp( -(G.v[i_scal]*G.v[i_scal]
                +G.v[i_scal+1]*G.v[i_scal+1]
                +G.v[i_scal+2]*G.v[i_scal+2] )/(4*alpha));
      A+=(coef.v[2*i]*cos_arr.v[i]-coef.v[2*i+1]*sin_arr.v[i])*tmp;

    }
  }
  else {
    for(int i=1; i < ngvec; i++) {
      //tmp= pwIntegral1D_f(G(i,0),l(0),alpha, sqrtalpha)
      //    *pwIntegral1D_f(G(i,1),l(1),alpha, sqrtalpha)
      //     *pwIntegral1D_f(G(i,2),l(2),alpha, sqrtalpha);
      //B+=(coef(i,1)*cos_arr(i)+coef(i,0)*sin_arr(i))*tmp;
      i_scal=i*3;
      tmp= pwIntegral1D_f(G.v[i_scal],l.v[0],alpha, sqrtalpha)
          *pwIntegral1D_f(G.v[i_scal+1],l.v[1],alpha, sqrtalpha)
          *pwIntegral1D_f(G.v[i_scal+2],l.v[2],alpha, sqrtalpha)
	  *exp( -(G.v[i_scal]*G.v[i_scal]
		 +G.v[i_scal+1]*G.v[i_scal+1]
		 +G.v[i_scal+2]*G.v[i_scal+2] )/(4*alpha));

      B+=(coef.v[2*i+1]*cos_arr.v[i]+coef.v[2*i]*sin_arr.v[i])*tmp;
    }
  }



  switch (sum%4) {
  case 0:
    return  A0+2.0*A;
    break;
  case 1:
    return -2.0*B;
    break;
  case 2:
    return -A0-2.0*A;
    break;
  case 3:
    return  2.0*B;
    break;
  default:
    cout <<" something wrong in function P_x_c "<<endl;
    return 0.0;
  }



}


inline int factorial(int n)
{
  assert(n>=0);
  //if (n<0) error ("factorial from negative number requested");
  int dummy;  
  dummy=1;
  for (int i=2;i<=n;i++) dummy*=i;
  return dummy;
}

inline int doublefactorial(int n)
{
  //if (n<-1) error ("doublefactorial from negative number requested");
  assert(n >= -1);
  // if(n==-1) return 1;
  int dummy=1;
  for (int i=n;i>1;i-=2)
    dummy*=i;
  return dummy;
}



inline int binomial2(const int n,const int k) {
  assert(k<=n);
  int num=1;
  int den=1;
  for(int i=n-k+1; i<= n; i++) {
    num*=i;
  }
  
  for(int i=2; i<=k; i++) {
    den*=i;
  }

  return num/den;
}

inline int binomial(int n, int k)
{
  if (k>n) error("n needs to be grater than k");
  return factorial(n)/(factorial(k)*factorial(n-k));
}


inline doublevar fsum(int k, int l1, int l2, doublevar x, doublevar y)
{
  doublevar tmp;
  int i,j;
  tmp=0.0;
  int qstart=max(-k, k-2*l2);
  int qend=min(k,2*l1-k);
  //for (int q=max(-k,k-2*l2);q<=min(k,2*l1-k);q+=2){
  for(int q=qstart; q <= qend; q+=2) {
    i=(k+q)/2;
    j=(k-q)/2;
    tmp+=binomial2(l1,i)*binomial2(l2,j)*pow(x,l1-i)*pow(y,l2-j);
  }
  return tmp;
}

inline doublevar Integ(int l1,int l2, doublevar x, doublevar y, doublevar exponent)
{
  doublevar tmp;
  tmp=0.0;
  for (int i=0;i<=(l1+l2)/2.0;i++)
    tmp+=fsum(2*i,l1,l2,x,y)*doublefactorial(2*i-1)*pow(2.0*exponent,-i);
  tmp*=sqrt(pi/exponent);
  return tmp;
}



inline doublevar S11(Array1 <int> & l, doublevar alpha)//, Array1 <doublevar> & center)
{
  //same as S12 but in the case when 1=2.
  return pow(sqrt(1.57079632679489661925/alpha),3)*pow(4.0*alpha,-l(0)-l(1)-l(2))
    *doublefactorial(2*l(0)-1)*doublefactorial(2*l(1)-1)*doublefactorial(2*l(2)-1);
  // pi/2=1.570796
}



inline doublevar S12(Array1 <int> & l1,Array1 <int> & l2, doublevar & alpha1, doublevar & alpha2,
	      Array1 <doublevar> & center1, Array1 <doublevar> & center2)
{
  //gives overlap matrix element for gaussian 1 and 2
  doublevar gamma,AB2,dummyI;
  //Array1 <doublevar> PA(3),PB(3);
  doublevar PA, PB;
  gamma=alpha1+alpha2;
  AB2=0.0;
  dummyI=1.0;
  for (int i=0;i<3;i++) {
    AB2+=(center1(i)-center2(i))*(center1(i)-center2(i));
    PA=(center2(i)-center1(i))*alpha2/gamma;
    PB=(center1(i)-center2(i))*alpha1/gamma;
    dummyI*=Integ(l1(i),l2(i),PA,PB,gamma);
  }
  return exp(-alpha1*alpha2*AB2/gamma)*dummyI;

}


//inline doublevar Integ(int l1,int l2, doublevar x, doublevar y, 
//                       doublevar exponent)
//{
//doublevar tmp;
//  tmp=0.0;
//  for (int i=0;i<=(l1+l2)/2.0;i++)
//    tmp+=fsum(2*i,l1,l2,x,y)*doublefactorial(2*i-1)*pow(2.0*exponent,-i);
//  tmp*=sqrt(pi/exponent);
//  return tmp;
//}



inline doublevar S12_f(Array1 <int> & l1, Array1 <int> & l2, 
		       doublevar alpha1, doublevar alpha2,
		       Array1 <doublevar> & center1, 
		       Array1 <doublevar> & center2) {
  doublevar gamma=alpha1+alpha2;
  doublevar invgamma=1.0/gamma;
  doublevar PA, PB, AB2=0.0, dummyI=1.0;
  doublevar sqrtgamma=sqrt(pi/gamma);
  doublevar diff;
  doublevar tmpint;
  for(int i=0; i< 3; i++) {
    diff=center2(i)-center1(i);
    AB2+=diff*diff;
    PA=diff*alpha2*invgamma;
    PB=-diff*alpha1*invgamma;
    tmpint=0.0;
    for(int j=0; j<= 0.5*(l1(i)+l2(i)); j++) {
      tmpint+=fsum(2*j,l1(i),l2(i),PA, PB)
              *doublefactorial(2*j-1)*pow(2.0*gamma, -j);
    }
    dummyI*=tmpint*sqrtgamma;
    //dummyI*=Integ(l1(i), l2(i), PA, PB, gamma);
  }
  return exp(-alpha1*alpha2*AB2*invgamma)*dummyI;
    
}


#endif //INTEGRALS_H_INCLUDED


//----------------------------------------------------------------------
