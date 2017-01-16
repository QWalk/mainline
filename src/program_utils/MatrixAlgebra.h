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


/*!
  This is basically a replacement for doubles, but always with a logarithmic representation.
  I've tried to make it so it automatically does the right thing with standard 
  multiplication. Automatic downconversion to double/complex is possible, but I've excluded it to 
  prevent accidental loss of precision.  Use the explicit val() function instead.
  Note also that there may be performance problems with using these automatic conversions, so if performance is critical, one may need to write some specialized code.
  */

template <class T> struct log_value { 
  T logval;
  int sign;
  T val() { return doublevar(sign)*exp(logval); }  //doublevar() a bit ugly, but it works ok
  log_value() { logval=0; sign=1; } 
  log_value(T t) {  } 
  log_value<T> & operator *=(const log_value<T> & right) {
    this->logval+=right.logval;
    this->sign*=right.sign;
    return *this;
  } 
};

template<> inline log_value<doublevar>::log_value(doublevar t) {
    doublevar ft=fabs(t);
    if(ft > 0) logval=log(ft); 
    else logval=-1e201;
    sign=t<0?-1:1; 
}

template<> inline log_value<dcomplex>::log_value(dcomplex t) { 
  sign=1;
  doublevar ft=abs(t);
  if(ft > 0) { 
    logval=dcomplex(log(abs(t)),arg(t));
  }
  else logval=dcomplex(-1e201,0.0);
}
typedef log_value<doublevar> log_real_value;
typedef log_value<dcomplex> log_complex_value;

//--------Finished with the class definitions
//Here are a few helper functions.
template <class T> inline log_value<T> operator*(T t, log_value<T> u) { 
  log_value<T> v(t);
  v*=u;
  return v;
}

template <class T> inline log_value<T>  operator*(log_value<T> u,double t) { return t*u; }
template <class T> inline log_value<T> operator*(log_value<T> t, log_value<T> u) { 
  log_value<T> v;
  v.logval=t.logval+u.logval;
  v.sign=t.sign*u.sign;
  return v;
}

inline log_value<dcomplex> operator*(log_value<doublevar> & t, log_value<dcomplex> & u) {
  log_value<dcomplex> v;
  v.logval=t.logval+u.logval;
  v.sign=t.sign*u.sign;
  return v;
}

//Try to safely sum a series of log_values
template <class T> inline log_value<T> sum(const Array1 <log_value<T> > & vec) { 
  T s=0;
  int n=vec.GetDim(0);
  int piv=0;
  for(int i=0; i< n; i++) 
    s+=T(vec(i).sign*vec(piv).sign)*exp(vec(i).logval-vec(piv).logval);
  return s*vec(piv);
}

//--------------------------------------------------------------

doublevar dot(const Array1 <doublevar> & a,  const Array1 <doublevar> & b);
dcomplex  dot(const Array1 <doublevar> & a,  const Array1 <dcomplex> & b);
dcomplex  dot(const Array1 <dcomplex> & a,  const Array1 <doublevar> & b);
dcomplex  dot(const Array1 <dcomplex> & a,  const Array1 <dcomplex> & b);


doublevar Determinant(const Array2 <doublevar> & a, const int n);
dcomplex Determinant(const Array2 <dcomplex> & a, const int n);

void MultiplyMatrices(const Array2 <doublevar> & a, const Array2 <doublevar> & b,
                      Array2 <doublevar> & c, int n);

void InvertMatrix(const Array2 <doublevar> & a, Array2 <doublevar> & a1, const int n);
doublevar InverseUpdateRow(Array2 <doublevar> & a1, const Array2 <doublevar> & a,
                           const int lRow, const int n);
doublevar InverseUpdateRow(Array2 <doublevar> & a1, const Array1 <doublevar> & newRow,
                           const int lRow, const int n);
//doublevar InverseGetNewRatio(const Array2 <doublevar> & a1, const Array1 <doublevar> & newCol,
//                             const int lCol, const int n);
//Get the new ratio without updating the inverse..

template <class T> T InverseGetNewRatio(const Array2 <T> & a1, 
     const Array1 <T> & newCol,
                             const int lCol, const int n) { 
  T f=T(0.0);  
  for(int i=0;i<n;++i) {
    f += a1(lCol,i)*newCol[i];
  }
  f =1.0/f;
  
  return f;
   
}


doublevar InverseGetNewRatioRow(const Array2 <doublevar> & a1, const Array1 <doublevar> & newRow,
                             const int lRow, const int n);
doublevar InverseUpdateColumn(Array2 <doublevar> & a1, const Array2 <doublevar> & a,
                              const int lCol, const int n);
doublevar InverseUpdateColumn(Array2 <doublevar> & a1, const Array1 <doublevar> & newCol,
                              const int lCol, const int n);

void TransposeMatrix(Array2 <doublevar> & a, const int n);
log_real_value TransposeInverseMatrix(const Array2 <doublevar> & a, Array2 <doublevar> & a1, const int n);
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

void EigenSystemSolverRealSymmetricMatrix(const Array2 < doublevar > & Ain, Array1 < doublevar> & evals, Array2 < doublevar> & evecs);

void DGGEV(Array2 <doublevar> & A, Array2 <doublevar> & B,
           Array1 <dcomplex> & alpha, 
           Array2 <doublevar> & VL, Array2 <doublevar> & VR);

void GeneralizedEigenSystemSolverRealSymmetricMatrices(const Array2 < doublevar > & Ain, const Array2 < doublevar> & Bin, Array1 < doublevar> & evals, Array2 < doublevar> & evecs);
void GeneralizedEigenSystemSolverComplexGeneralMatrices(Array2 < dcomplex > & Ain, 
         Array1 <dcomplex> & W, Array2 <dcomplex> & VL, Array2 <dcomplex> & VR); 
void GeneralizedEigenSystemSolverRealGeneralMatrices(Array2 < doublevar > & Ain, 
         Array1 <dcomplex> & W, Array2 <doublevar> & VL, Array2 <doublevar> & VR); 


void sort_abs_values_descending(const Array1 <doublevar> & vals, Array1 <doublevar> & newvals, Array1 <int> & list);
void sort_according_list(const Array1 <doublevar> & vals, Array1 <doublevar> & newvals, const Array1 <int> & list);
//------------------------------------------------------------------------
//IF Lapack was used
#ifdef USE_LAPACK
extern "C" { 
  void dsyevr_(char *JOBZp, char *RANGEp, char *UPLOp, int *Np,
               double *A, int *LDAp, double *VLp, double *VUp,
               int *ILp, int *IUp, double *ABSTOLp, int *Mp,
               double *W, double *Z, int *LDZp, int *ISUPPZ,
               double *WORK, int *LWORKp, int *IWORK, int *LIWORKp,
               int *INFOp);
  double dlamch_(char *CMACHp);
  void dsygv_(int *itype, char *jobz, char *uplo, int *n, double *a, int *lda, 
              double *b, int *ldb, double *w, double *WORK, int *IWORK, int *info);
  void dgetrf_(int * m, int * n, double * A, int * lda, int * ipiv, int *info);
  void dgetrs_(char * trans, int *n, int * nrhs, double * A, int * lda, 
               int *ipiv, double * B, int * ldb, int * info);

  void dgeev_(char * JOBVL, char * JOBVR,int * N,double * A,int * lda,
      double * WR, double * WI,double *VL, int *LDVL,double * VR, int *LDVR,
      double * WORK, int * LWORK, int * INFO);
  
  void zgeev_(char * JOBVL, char * JOBVR,int * N,void * A,int * lda,
      void * W, void *VL, int *LDVL,void * VR, int *LDVR,
      void * WORK, int * LWORK,void * RWORK, int * INFO);

  void dggev_(char * JOBVL, char * JOBVR, int * N, void * A, int * lda,
      void * B, int * ldb, void * alphar, void * alphai, void * beta,
      void* vl, int * ldvl, void * vr, int * ldvr, void * work, int * lwork, int * info);
};
#endif

//----------------------------------------------------------------------
//complex versions of determinant stuff

int ludcmp(Array2 <dcomplex > & a, const int n, 
           Array1 <int> & indx, doublevar & d);
void lubksb(Array2 < dcomplex > & a, int n, 
            Array1 <int> & indx, Array1 <dcomplex > & b);

log_value <dcomplex>
TransposeInverseMatrix(const Array2 <dcomplex > & a, 
                       Array2 <dcomplex> & a1, 
                       const int n);
dcomplex
InverseUpdateColumn(Array2 <dcomplex > & a1, 
                    const Array1 <dcomplex > & newCol,
                    const int lCol, const int n);
dcomplex InverseUpdateRow(Array2 <dcomplex> & a1, const Array1 <dcomplex> & newRow,
                           const int lRow, const int n);
void Jacobi(const Array2 < dcomplex > & Ain, Array1 <doublevar> & evals, Array2 < dcomplex > & evecs);

#endif // MATRIXALGEBRA_H_INCLUDED
//--------------------------------------------------------------------------
