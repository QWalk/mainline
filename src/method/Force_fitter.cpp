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

#include "Force_fitter.h"
#include "MatrixAlgebra.h"

void Force_fitter::setup(doublevar  maxcut, int nexpansion ) {
  m=2;
  cut=maxcut;
  nexp=nexpansion;
  coeff.Resize(nexp);
  Array2 <doublevar> S(nexp, nexp);
  Array1 <doublevar> h(nexp);
  for(int i=0; i< nexp; i++) {
    for(int j=0; j< nexp; j++) {
      S(i,j)=pow(cut,m+i+j+1)/(m+i+j+1);
    }
    h(i)=pow(cut,i+1)/(i+1);
  }
  
  Array2 <doublevar> Sinv(nexp, nexp);
  InvertMatrix(S, Sinv, nexp);
  
  coeff=0;
  for(int i=0; i< nexp; i++) 
    for(int j=0;j< nexp;j++) 
      coeff(i)+=Sinv(i,j)*h(j);
}


void Force_fitter::fit_force(Array1 <doublevar> & r, 
                             Array1 <doublevar> & bare_force,
                             Array1 <doublevar> & fit_force) {
  int ndim=bare_force.GetDim(0);

  fit_force=bare_force; return;

  fit_force.Resize(ndim);
  fit_force=0;

  if(r(0) < cut) {
    for(int d=0; d< ndim; d++) {
      doublevar basis=0;
      for(int i=0; i< nexp; i++) {
        basis+=coeff(i)*pow(r(0),i+m);
      }
      fit_force(d)=bare_force(d)*basis;
    }
  }
  else {
    fit_force=bare_force;
  }


}
    
