/*
 
Copyright (C) 2007 Burkhard Militzer
 with some modifications by Lucas K. Wagner
 and extensions to Pfaffian algebra by Michal Bajdich

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
#ifndef MATRIXALGEBRA_H_INCLUDED
#define MATRIXALGEBRA_H_INCLUDED

#include "Array.h"
#include "Qmc_std.h"

doublevar Determinant(const Array2 <doublevar> & a, const int n);
dcomplex Determinant(const Array2 <dcomplex> & a, const int n);

void MultiplyMatrices(const Array2 <doublevar> & a, const Array2 <doublevar> & b,
                      Array2 <doublevar> & c, int n);

void InvertMatrix(const Array2 <doublevar> & a, Array2 <doublevar> & a1, const int n);
doublevar InverseUpdateRow(Array2 <doublevar> & a1, const Array2 <doublevar> & a,
                           const int lRow, const int n);
doublevar InverseUpdateRow(Array2 <doublevar> & a1, const Array1 <doublevar> & newRow,
                           const int lRow, const int n);
doublevar InverseUpdateColumn(Array2 <doublevar> & a1, const Array2 <doublevar> & a,
                              const int lCol, const int n);
doublevar InverseUpdateColumn(Array2 <doublevar> & a1, const Array1 <doublevar> & newCol,
                              const int lCol, const int n);

void TransposeMatrix(Array2 <doublevar> & a, const int n);
doublevar TransposeInverseMatrix(const Array2 <doublevar> & a, Array2 <doublevar> & a1, const int n);
doublevar TransposeInverseUpdateRow(Array2 <doublevar> & a1, const Array2 <doublevar> & a,
                                    const int lCol, const int n);
doublevar TransposeInverseUpdateColumn(Array2 <doublevar> & a1, const Array2 <doublevar> & a,
                                       const int lRow, const int n);
doublevar Pfaffian_nopivot(const Array2 <doublevar> & a);

doublevar Pfaffian_partialpivot(const Array2 <doublevar> & a);
     
doublevar UpdateInversePfaffianMatrix(Array2 <doublevar> & a, Array1 <doublevar> & row,
                              Array1 <doublevar> & column, int n);
doublevar GetUpdatedPfaffianValue(Array2 <doublevar> & in, 
                                      Array1 <doublevar> & row, 
                                      int e);
doublevar PfaffianInverseMatrix(const Array2 <doublevar> & a, Array2 <doublevar> & a1);


inline doublevar cabs(dcomplex & a) {
  return sqrt(a.real()*a.real()+a.imag()*a.imag());
}

int ludcmp(Array2 <doublevar> & a, const int n, Array1 <int> & indx,
           doublevar & d);
void lubksb(Array2 <doublevar> & a, int n, Array1 <int> & indx,
            Array1 <doublevar> & b);
void Jacobi(const Array2 < doublevar > & Ain, Array1 < doublevar> & evals, Array2 < doublevar> & evecs);

//----------------------------------------------------------------------
//complex versions of determinant stuff

int ludcmp(Array2 <dcomplex > & a, const int n, 
           Array1 <int> & indx, doublevar & d);
void lubksb(Array2 < dcomplex > & a, int n, 
            Array1 <int> & indx, Array1 <dcomplex > & b);

dcomplex
TransposeInverseMatrix(const Array2 <dcomplex > & a, 
                       Array2 <dcomplex> & a1, 
                       const int n);
dcomplex
InverseUpdateColumn(Array2 <dcomplex > & a1, 
                    const Array1 <dcomplex > & newCol,
                    const int lCol, const int n);
void Jacobi(const Array2 < dcomplex > & Ain, Array1 <doublevar> & evals, Array2 < dcomplex > & evecs);

#endif // MATRIXALGEBRA_H_INCLUDED
//--------------------------------------------------------------------------
