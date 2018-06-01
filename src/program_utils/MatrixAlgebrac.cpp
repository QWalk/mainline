/*
 
Copyright (C) 2007 Burkhard Militzer
 with modifications by Lucas K. Wagner
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
//--------------------------------------------------------------------------
// src/MatrixAlgebrac.cpp
//
//
// Matrix operations
//
// Burkhard Militzer                                    Urbana 4-1-99
//
#include <math.h>
#include "Array.h"
#include "Qmc_std.h"
#include "MatrixAlgebra.h"

#ifdef USE_LAPACK
//C++ wrappwr of Lapack routines
doublevar dlamch(char CMACH)
{
  return dlamch_(&CMACH);
}

int dsyevr(char JOBZ, char RANGE, char UPLO, int N,
           double *A, int LDA, double VL, double VU,
           int IL, int IU, double ABSTOL, int *M,
           double *W, double *Z, int LDZ, int *ISUPPZ,
           double *WORK, int LWORK, int *IWORK, int LIWORK){
  
  int INFO;
  dsyevr_(&JOBZ, &RANGE, &UPLO, &N, A, &LDA, &VL, &VU,
          &IL, &IU, &ABSTOL, M, W, Z, &LDZ, ISUPPZ,
          WORK, &LWORK, IWORK, &LIWORK, &INFO);
  return INFO;
}

int dsygv(int itype, char jobz, char uplo,  int  n,  
          double *a,  int  lda,  double *b, int ldb, double *w, 
          double *WORK, int IWORK){
  int INFO;
  dsygv_(&itype, &jobz, &uplo, &n, a, &lda, b, &ldb, w, WORK, &IWORK, &INFO);
  return INFO;
}


int dgeev(char JOBVL, char JOBVR,int N,doublevar * A,int lda,
    doublevar * WR, doublevar * WI, doublevar* VL, int LDVL,doublevar * VR, int LDVR,
    doublevar * WORK, int LWORK) { 
  int INFO;
  dgeev_(&JOBVL,&JOBVR, &N,A,&lda,WR,WI,VL,&LDVL,VR,&LDVR,WORK, &LWORK,&INFO);
  return INFO;
}



int zgeev(char JOBVL, char JOBVR,int N,dcomplex * A,int lda,
    dcomplex * W, dcomplex* VL, int LDVL,dcomplex * VR, int LDVR,
    dcomplex * WORK, int LWORK, dcomplex * RWORK) { 
  int INFO;
  zgeev_(&JOBVL,&JOBVR, &N,A,&lda,W,VL,&LDVL,VR,&LDVR,WORK, &LWORK,RWORK,&INFO);
  return INFO;
}

    


int dgetrf(int m, int n, double * A, int lda, int * ipiv) { 
  int info;
  dgetrf_(&m, &n, A, &lda, ipiv, &info);
  return info;
}

int dgetrs(char trans, int n, int nrhs, double * A, int lda, int * ipiv, 
           double * B, int ldb) { 
  int info;
  dgetrs_(&trans, &n, &nrhs, A, &lda, ipiv, B, &ldb, &info);
  return info;
}
#endif


const doublevar TINY=1.0e-20;

Array2 <doublevar> tmp2;
Array1 <doublevar> tmp11,tmp12;
Array1 <int> itmp1;

// LU decomposition of matrix a, which is overwritten
//Modified to return 1 if successful, 0 if not.
int ludcmp(Array2 <doublevar> & a, const int n, Array1 <int> & indx, doublevar & d)
{
  //	cout << "ludcmp" << endl;
  Array1 <doublevar>& vv(tmp11);
  vv.Resize(n);

  //cout << "start " << endl;
  d=1.0;
  for (int i=0;i<n;++i)
  {
    doublevar big=0.0;
    for (int j=0;j<n;++j)
    {
      //doublevar temp;
      //if ((temp=fabs(a(i,j))) > big) big=temp;
      doublevar temp=fabs(a(i,j));
      if(temp>big)
        big=temp;
    }
    if (big == 0.0)
      return 0;//error("Singular matrix in routine ludcmp");
    vv[i]=1.0/big;
  }

  //cout << "middle " << endl;
  for (int j=0;j<n;++j)
  {
    int imax;
    //cout << "j " << j << endl;
    for (int i=0;i<j;++i)
    {
      //cout << "1i " << i << endl;
      doublevar sum=a(i,j);
      for (int k=0;k<i;++k)
      {
        sum -= a(i,k)*a(k,j);
      }
      a(i,j)=sum;
    }
    //cout << "imax " << imax << endl;
    doublevar big(0.0);
    for (int i=j;i<n;++i)
    {
      //cout << " i " << i << endl;
      doublevar sum=a(i,j);
      for (int k=0;k<j;++k)
        sum -= a(i,k)*a(k,j);
      a(i,j)=sum;
      //doublevar dum;
      //cout << "dum part " << endl;
      //if ( (dum=vv[i]*fabs(sum)) >= big) {
      doublevar dum=vv[i]*fabs(sum);
      if(dum >=big)
      {
        big=dum;
        imax=i;
      }
    }
    //cout << "3 " << endl;
    if (j != imax)
    {
      for (int k=0;k<n;++k)
      {
        doublevar dum(a(imax,k));
        a(imax,k)=a(j,k);
        a(j,k)=dum;
      }
      d = -d;
      vv[imax]=vv[j];
    }
    indx[j]=imax;
    if (a(j,j) == 0.0)
    {
      a(j,j)=TINY;
      //Write(n);
      for (int q=0;q<n;++q)
        for (int qq=0;qq<n;++qq)
          //Write3(q,qq,a(q,qq));
          //error("Singular matrix in routine ludcmp II.");
          return 0;
    }
    //cout << "5 " << endl;
    if (j != n-1)
    {
      doublevar dum=1.0/(a(j,j));
      for (int i=j+1;i<n;++i)
        a(i,j) *= dum;
    }
  }
  return 1;
}

void lubksb(Array2 <doublevar> & a, int n, Array1 <int> & indx, Array1 <doublevar> & b)
{
  int ii(-1);
  for (int i=0;i<n;++i)
  {
    int ip=indx[i];
    doublevar sum(b[ip]);
    b[ip]=b[i];
    if (ii>=0)
      for (int j=ii;j<i;++j)
        sum -= a(i,j)*b[j];
    else
      if (sum!=0.0)
        ii=i;
    b[i]=sum;
  }

  for (int i=n-1;i>=0;--i)
  {
    doublevar sum(b[i]);
    for (int j=i+1;j<n;++j)
      sum -= a(i,j)*b[j];
    b[i]=sum/a(i,i);
  }
}

// Invert matrix a, but do not change a
// Inverse in a1
void InvertMatrix(const Array2 <doublevar> & a, Array2 <doublevar> & a1, const int n)
{
  Array2 <doublevar>& temp(tmp2);
  temp.Resize(n,n);
  Array1 <int>& indx(itmp1);
  indx.Resize(n);
  doublevar d;
  // a(i,j) first index i is row index (convention)
  // elements of column vectors are stored contiguous in memory in C style arrays
  // a(i) refers to a column vector

  // calculate the inverse of the transposed matrix because this
  // allows to pass a column vector to lubksb() instead of a row

  // put the transposed matrix in temp
  for(int i=0;i<n;++i)
  {
    for(int j=0;j<n;++j)
    {
      temp(i,j)=a(j,i);
      a1(i,j)=0.0;
    }
    a1(i,i)=1.0;
  }
  

  if(!ludcmp(temp,n,indx,d)) error("singular matrix in inversion");

  
  
  for(int j=0;j<n;++j)
  {
    // get column vector
    Array1 <doublevar> yy;//(a1(j));
    yy.refer(a1(j));
    lubksb(temp,n,indx,yy);
  }

  
  
}


int InvertPfaffianMatrix(const Array2 <doublevar> & a, Array2 <doublevar> & a1, const int n)
{
  Array2 <doublevar>& temp(tmp2);
  temp.Resize(n,n);
  Array1 <int>& indx(itmp1);
  indx.Resize(n);
  doublevar d;

  // a(i,j) first index i is row index (convention)
  // elements of column vectors are stored contiguous in memory in C style arrays
  // a(i) refers to a column vector

  // calculate the inverse of the transposed matrix because this
  // allows to pass a column vector to lubksb() instead of a row

  // put the transposed matrix in temp
  for(int i=0;i<n;++i)
  {
    for(int j=0;j<n;++j)
    {
      temp(i,j)=a(j,i);
      a1(i,j)=0.0;
    }
    a1(i,i)=1.0;
  }

  if(!ludcmp(temp,n,indx,d)) {
    cout <<"ERROR: singular matrix in inversion\n";
    return 0;
  }

  for(int j=0;j<n;++j)
  {
    // get column vector
    Array1 <doublevar> yy;//(a1(j));
    yy.refer(a1(j));
    lubksb(temp,n,indx,yy);
  }
  return 1;
}



// Calculate the transpose inverse of matrix a
// and return the determinant
log_real_value TransposeInverseMatrix(const Array2 <doublevar> & a, Array2 <doublevar> & a1, const int n)
{
  Array2 <doublevar> &temp(tmp2);
  temp.Resize(n,n);
  Array1 <int>& indx(itmp1);
  indx.Resize(n);
  doublevar d=1;
 
  log_real_value logdet;
  logdet.logval=0; logdet.sign=1;
  //for(int i=0; i< n; i++) { 
  //  cout << "matrix ";
  //  for(int j=0; j< n; j++) cout << a(i,j) << " ";
  //  cout << endl;
  //}
  
#ifdef USE_LAPACK
  //LAPACK routines don't handle n==1 case??
  if(n==1) { 
    a1(0,0)=1.0/a(0,0);
    logdet.logval=log(fabs(a(0,0)));
    logdet.sign=a(0,0)<0?-1:1;
    return logdet;
  }
  else { 
  
    for(int i=0; i < n;++i) {
      for(int j=0; j< n; ++j) { 
        temp(j,i)=a(i,j);
        a1(i,j)=0.0;
      }
      a1(i,i)=1.0;
    }
    if(dgetrf(n, n, temp.v, n, indx.v)> 0) { 
      return 0.0;
    }
    for(int j=0; j< n; ++j) { 
      dgetrs('N',n,1,temp.v,n,indx.v,a1.v+j*n,n);
    }
  }

  for(int i=0; i< n; i++) { 
    if(indx(i)!=i+1) logdet.sign*=-1;
    logdet.logval+=log(fabs(temp(i,i)));
    if(temp(i,i) <0) logdet.sign*=-1;
  }

  //cout << " det " << det << " logval " << logdet.val() << endl;
  //return det;
  return logdet;
//#endif  
#else 
  
  // a(i,j) first index i is row index (convention)
  // elements of column vectors are stored contiguous in memory in C style arrays
  // a(i) refers to a column vector

  // calculate the inverse of the transposed matrix because this
  // allows to pass a column vector to lubksb() instead of a row

  // put the transposed matrix in temp
  //cout << "temp " << endl;
  for(int i=0;i<n;++i)
  {
    for(int j=0;j<n;++j)
    {
      temp(i,j)=a(i,j);
      a1(i,j)=0.0;
    }
    a1(i,i)=1.0;
  }

  //cout << "ludcmp" << endl;
  //if the matrix is singular, the determinant is zero.
  d=1;
  if(ludcmp(temp,n,indx,d)==0)
    return 0;

  //cout << "lubksb" << endl;

  for(int j=0;j<n;++j)
  {
    // get column vector
    Array1 <doublevar> yy;//(a1(j));
    yy.refer(a1(j));
    lubksb(temp,n,indx,yy);
  }
  
  //for(int j=0;j<n;++j) {
  //  d *= temp(j,j);
  //}

  logdet.logval=0;
  logdet.sign=1;
  for(int i=0; i< n; i++) { 
    if(indx(i)!=i) logdet.sign*=-1;
    logdet.logval+=log(fabs(temp(i,i)));
    if(temp(i,i) <0) logdet.sign*=-1;
  }

  //cout << " det " << d << " logval " << logdet.val() << endl;
  

  return logdet;
#endif
}

// Return det|a| and leave matrix a constant
doublevar Determinant(const Array2 <doublevar> & a, const int n)
{
  Array2 <doublevar>& temp(tmp2);
  temp.Resize(n,n);
  Array1 <int>& indx(itmp1);
  indx.Resize(n);
  doublevar d;

  // a(i,j) first index i is row index (convention)
  // elements of column vectors are stored contiguous in memory in C style arrays
  // a(i) refers to a column vector

  // calculate the inverse of the transposed matrix because this
  // allows to pass a column vector to lubksb() instead of a row

  // put the transposed matrix in temp
//#ifdef USE_LAPACK
/*  if(n==0) return 1;
  for(int i=0; i < n;++i) {
    for(int j=0; j< n; ++j) { 
      temp(i,j)=a(i,j);
    }
  }
  dgetrf(n, n, temp.v, n, indx.v);
  double det=1;
  for(int i=0; i< n; i++) { 
    if(indx(i) != i+1)
      det*= -temp(i,i);
    else det*=temp(i,i);
  }
  return det;
  */
/*
  for(int i=0; i< n; i++) { 
    if(indx(i)!=i+1) logdet.sign*=-1;
    logdet.logval+=log(fabs(temp(i,i)));
    if(temp(i,i) <0) logdet.sign*=-1;
  }

  //cout << " det " << det << " logval " << logdet.val() << endl;
  //return det;
  return logdet;
  */
//#endif  
//#else 
 
  for(int i=0;i<n;++i) {
    for(int j=0;j<n;++j) {
      temp(i,j)=a(i,j);
    }
  }

  ludcmp(temp,n,indx,d);

  for(int j=0;j<n;++j) {
    d *= temp(j,j);
  }
  return d;
//#endif
}


// Return det|a| and leave matrix a constant(complex version)
dcomplex Determinant(const Array2 <dcomplex> & a, const int n)
{
  Array2 <dcomplex> temp(n,n);
  temp.Resize(n,n);
  Array1 <int>& indx(itmp1);
  indx.Resize(n);
  dcomplex d;

  // a(i,j) first index i is row index (convention)
  // elements of column vectors are stored contiguous in memory in C style arrays
  // a(i) refers to a column vector

  // calculate the inverse of the transposed matrix because this
  // allows to pass a column vector to lubksb() instead of a row

  // put the transposed matrix in temp
  for(int i=0;i<n;++i)
  {
    for(int j=0;j<n;++j)
    {
      temp(i,j)=a(i,j);
    }
  }

  doublevar det_sign;
  ludcmp(temp,n,indx,det_sign);

  d=det_sign;
  for(int j=0;j<n;++j)
  {
    d *= temp(j,j);
  }
  return d;
}


// transpose matrix a
void TransposeMatrix(Array2 <doublevar> & a, const int n)
{
  for(int i=0;i<n;++i)
  {
    for(int j=i+1;j<n;++j)
    {
      doublevar & r(a(i,j));
      doublevar & rt(a(j,i));
      doublevar x=r;
      r=rt;
      rt=x;
    }
  }
}

// multiply c=a*b
void MultiplyMatrices(const Array2 <doublevar> & a, const Array2 <doublevar> & b,
                      Array2 <doublevar> & c, int n)
{
  for(int i=0;i<n;++i)
  {
    for(int j=0;j<n;++j)
    {
      doublevar & cc(c(i,j));
      cc=0.0;
      for(int k=0;k<n;++k)
      {
        cc +=a(i,k)*b(k,j);
      }
    }
  }
}

//////////////////////// Update Inverse /////////////////////


// Update inverse a1 after row in matrix a has changed
// get new row  out of matrix a
//
// a1 = old inverse
// a  = new matrix with new row lRow
// returns Det(a_old)/Det(a_new)
doublevar InverseUpdateRow(Array2 <doublevar> & a1, const Array2 <doublevar> & a,
                           const int lRow, const int n)
{
  Array1 <doublevar> & tmpColL(tmp11);
  tmpColL.Resize(n);
  Array1 <doublevar> & prod(tmp12);
  prod.Resize(n);

  doublevar f=0.0;
  for(int i=0;i<n;++i)
  {
    f += a(lRow,i)*a1(i,lRow);
  }
  f =-1.0/f;

  for(int j=0;j<n;++j)
  {
    prod[j]   =0.0;
    tmpColL[j]=a1(j,lRow);
    for(int i=0;i<n;++i)
    {
      prod[j] += a(lRow,i)*a1(i,j);
    }
    prod[j] *= f;
  }

  for(int ii=0;ii<n;++ii)
  {
    doublevar & t(tmpColL[ii]);
    for(int j=0;j<n;++j)
    {
      a1(ii,j) += t*prod[j];
    }
  }

  f = -f;
  for(int i=0;i<n;++i)
  {
    a1(i,lRow) = f*tmpColL[i];
  }
  return f;
}

// Update inverse a1 after row in matrix a has changed
// get new row out of array1 newRow
//
// a1 = old inverse
// newRow  = new row lRow in new matrix a_new
// returns Det(a_old)/Det(a_new)
doublevar InverseUpdateRow(Array2 <doublevar> & a1, const Array1 <doublevar> & newRow,
                           const int lRow, const int n)
{
  Array1 <doublevar> & tmpColL(tmp11);
  tmpColL.Resize(n);
  Array1 <doublevar> & prod(tmp12);
  prod.Resize(n);

  doublevar f=0.0;
  for(int i=0;i<n;++i) {
    f += newRow[i]*a1(i,lRow);
  }
  f =-1.0/f;

  for(int j=0;j<n;++j){
    prod[j]   =0.0;
    tmpColL[j]=a1(j,lRow);
    for(int i=0;i<n;++i) {
      prod[j] += newRow[i]*a1(i,j);
    }
    prod[j] *= f;
  }

  for(int i=0;i<n;++i) {
    doublevar & t(tmpColL[i]);
    for(int j=0;j<n;++j) {
      a1(i,j) += t*prod[j];
    }
  }

  f = -f;
  for(int i=0;i<n;++i) {
    a1(i,lRow) = f*tmpColL[i];
  }
  return f;
}

//We write different versions in case we want to use BLAS routines, which don't template well..
dcomplex InverseUpdateRow(Array2 <dcomplex> & a1, const Array1 <dcomplex> & newRow,
                           const int lRow, const int n)
{
  Array1 <dcomplex>  tmpColL;
  tmpColL.Resize(n);
  Array1 <dcomplex>  prod;
  prod.Resize(n);

  dcomplex f=0.0;
  for(int i=0;i<n;++i) {
    f += newRow[i]*a1(i,lRow);
  }
  f =-1.0/f;

  for(int j=0;j<n;++j){
    prod[j]   =0.0;
    tmpColL[j]=a1(j,lRow);
    for(int i=0;i<n;++i) {
      prod[j] += newRow[i]*a1(i,j);
    }
    prod[j] *= f;
  }

  for(int i=0;i<n;++i) {
    dcomplex t(tmpColL[i]);
    for(int j=0;j<n;++j) {
      a1(i,j) += t*prod[j];
    }
  }

  f = -f;
  for(int i=0;i<n;++i) {
    a1(i,lRow) = f*tmpColL[i];
  }
  return f;
}




// Update inverse a1 after column in matrix a has changed
// get new column out of matrix a
//
// a1= old inverse
// a = new matrix with new column lCol
// returns Det(a_old)/Det(a_new)
doublevar InverseUpdateColumn(Array2 <doublevar> & a1, const Array2 <doublevar> & a,
                              const int lCol, const int n)
{
  Array1 <doublevar> & tmpColL(tmp11);
  tmpColL.Resize(n);
  Array1 <doublevar> & prod(tmp12);
  prod.Resize(n);

  doublevar f=0.0;
  for(int i=0;i<n;++i)
  {
    f += a1(lCol,i)*a(i,lCol);
  }
  f =-1.0/f;

  for(int j=0;j<n;++j)
  {
    tmpColL[j]=a1(lCol,j);
    prod[j]   =0.0;
    for(int i=0;i<n;++i)
    {
      prod[j] += a1(j,i)*a(i,lCol);
    }
    prod[j] *= f;
  }

  for(int i=0;i<n;++i)
  {
    doublevar & p(prod[i]);
    for(int j=0;j<n;++j)
    {
      a1(i,j) += tmpColL[j]*p;
    }
  }

  f = -f;
  for(int j=0;j<n;++j)
  {
    a1(lCol,j) = f*tmpColL[j];
  }
  return f;
}




//Get the new ratio without updating the inverse.. for a row change (transpose of InverseGetNewRatio)

doublevar InverseGetNewRatioRow(const Array2 <doublevar> & a1, const Array1 <doublevar> & newRow,
                             const int lRow, const int n) { 
  doublevar f=0.0;  
  for(int i=0;i<n;++i) {
    f += a1(i,lRow)*newRow[i];
  }
  f =1.0/f;
  
  return f;
   
}


// Update inverse a1 after column in matrix a has changed
// get new column out of array1 newCol
//
// a1= old inverse
// newCol = new column lCol in the new matrix a_new
// returns Det(a_old)/Det(a_new)
doublevar InverseUpdateColumn(Array2 <doublevar> & a1, const Array1 <doublevar> & newCol,
                              const int lCol, const int n)
{
  Array1 <doublevar> & tmpColL(tmp11);
  tmpColL.Resize(n);
  Array1 <doublevar> & prod(tmp12);
  prod.Resize(n);

  doublevar f=0.0;

#ifdef USE_BLAS
  int a1size=a1.GetDim(1);

  doublevar * a1col=a1.v+lCol*a1size;

  f=cblas_ddot(n,a1col, 1, newCol.v, 1);
  f=-1.0/f;

  cblas_dcopy(n,a1col,1,tmpColL.v,1);
  
  cblas_dgemv(CblasRowMajor,CblasNoTrans,n,n,
              1.0,a1.v,a1size,
              newCol.v,1,
              0.0,prod.v,1);

  cblas_dscal(n,f,prod.v,1);

  cblas_dger(CblasRowMajor, n,n,1.0,
             prod.v,1,
             tmpColL.v,1,
             a1.v,a1size);
  f=-f;
  cblas_dcopy(n,tmpColL.v,1,a1col,1);
  cblas_dscal(n,f,a1col,1);

#else 

  for(int i=0;i<n;++i)
  {
    f += a1(lCol,i)*newCol[i];
  }
  f =-1.0/f;

  for(int j=0;j<n;++j)
  {
    tmpColL[j]=a1(lCol,j);
    prod[j]   =0.0;
    for(int i=0;i<n;++i)
    {
      prod[j] += a1(j,i)*newCol[i];
    }
    prod[j] *= f;
  }

  for(int i=0;i<n;++i)
  {
    doublevar & p(prod[i]);
    for(int j=0;j<n;++j)
    {
      a1(i,j) += tmpColL[j]*p;
    }
  }

  f = -f;
  for(int j=0;j<n;++j)
  {
    a1(lCol,j) = f*tmpColL[j];
  }

#endif

  return f;
}

///////////////// Update Transpose Inverse /////////////////////


// Update transpose inverse a1 after row in matrix a has changed
// get new row out of matrix a
// (This is actually the modified routine InverseUpdateColumn)
//
// a1= old inverse
// a = new matrix with new column lCol
// returns Det(a_old)/Det(a_new)
doublevar TransposeInverseUpdateRow(Array2 <doublevar> & a1, const Array2 <doublevar> & a,
                                    const int lCol, const int n)
{
  Array1 <doublevar> & tmpColL(tmp11);
  tmpColL.Resize(n);
  Array1 <doublevar> & prod(tmp12);
  prod.Resize(n);

  doublevar f=0.0;
  for(int i=0;i<n;++i)
  {
    f += a1(lCol,i)*a(lCol,i);
  }
  f =-1.0/f;

  for(int j=0;j<n;++j)
  {
    tmpColL[j]=a1(lCol,j);
    prod[j]   =0.0;
    for(int i=0;i<n;++i)
    {
      prod[j] += a1(j,i)*a(lCol,i);
    }
    prod[j] *= f;
  }

  for(int i=0;i<n;++i)
  {
    doublevar & p(prod[i]);
    for(int j=0;j<n;++j)
    {
      a1(i,j) += tmpColL[j]*p;
    }
  }

  f = -f;
  for(int j=0;j<n;++j)
  {
    a1(lCol,j) = f*tmpColL[j];
  }
  return f;
}

// Update inverse a1 after column in matrix a has changed
// get new column out of matrix a
// (This is actually the modified routine InverseUpdateRow)
//
// a1 = old inverse
// a  = new matrix with new row lRow
// returns Det(a_old)/Det(a_new)
doublevar TransposeInverseUpdateColumn(Array2 <doublevar> & a1, const Array2 <doublevar> & a,
                                       const int lRow, const int n)
{
  Array1 <doublevar> & tmpColL(tmp11);
  tmpColL.Resize(n);
  Array1 <doublevar> & prod(tmp12);
  prod.Resize(n);

  doublevar f=0.0;
  for(int i=0;i<n;++i)
  {
    f += a(i,lRow)*a1(i,lRow);
  }
  f =-1.0/f;

  for(int j=0;j<n;++j)
  {
    prod[j]   =0.0;
    tmpColL[j]=a1(j,lRow);
    for(int i=0;i<n;++i)
    {
      prod[j] += a(i,lRow)*a1(i,j);
    }
    prod[j] *= f;
  }

  for(int i=0;i<n;++i)
  {
    doublevar & t(tmpColL[i]);
    for(int j=0;j<n;++j)
    {
      a1(i,j) += t*prod[j];
    }
  }

  f = -f;
  for(int i=0;i<n;++i)
  {
    a1(i,lRow) = f*tmpColL[i];
  }
  return f;
}
/*
doublevar ulec() {
  static long int is1=12345;
  static long int is2=56789;
  long int k,iz;
  k=is1/53668;
  is1=is1-k*53668;
  is1=40014*is1-k*12211;
  if (is1<0) is1 +=2147483563;
  k=is2/52774;
  is2=is2-k*52774;
  is2=40692*is2-k*3791;
  if (is2<0) is2 +=2147483399;
  iz=is1-is2;
  if (iz<0) iz +=2147483562;
  return (4.656613057e-10*iz);
}
*/
//--------------------------------------------------------------------------


doublevar Pfaffian_nopivot(const Array2 <doublevar> & in){
  //  calculate the pfaffian of skew-symmetric matrix
  //  no pivoting
  assert(in.dim[0]==in.dim[1]);
  if (in.dim[0]%2!=0) return 0.0;
  int n=in.dim[0]/2;
  Array2 <doublevar> tmp(2*n,2*n);
  doublevar PF=1.0;
  doublevar fac;
  for (int i=0;i<2*n;i++)
    for (int l=0;l<2*n;l++)
      tmp(i,l)=in(i,l);
 
  
  for (int i=0;i<2*n;i=i+2){
    for (int j=i+2;j<2*n;j++){
      fac=-tmp(i,j)/tmp(i,i+1);
      for (int k=i+1;k<2*n;k++){
        tmp(k,j)=tmp(k,j)+fac*tmp(k,i+1);
        tmp(j,k)=tmp(j,k)+fac*tmp(i+1,k);
      }
    }
    PF=PF*tmp(i,i+1);
  }
  return PF;
}



int RowPivoting(Array2 <doublevar> & tmp, int i, int n){
  doublevar big;
  doublevar temp;
  doublevar TINY=1e-20;
  Array1 <doublevar> backup(2*n) ;
  int d=1;
  int k=0;
  big=0.0;
  for (int j=i+1;j<2*n;j++){
    temp=fabs(tmp(i,j));
    if(temp > big){
      big=temp;
      k=j;
      // cout <<"found one";
    }
  }
  if (big<TINY){
    cout <<"Singular row in matrix!!! "<<endl;
    tmp(i,i+1)=TINY;
  }
  
  // cout <<"row: "<<i+1<<"   element: "<<k+1<<"    value: "<<big<<endl;
  
  if (k!=i+1){
    //exchange k-th column with 2-nd column;
     for (int j=i;j<2*n;j++){
       backup(j)=tmp(j,i+1);
       tmp(j,i+1)=tmp(j,k);
       tmp(j,k)=backup(j);
     }
     //exchange k-th row with 2-nd row;
     for (int j=i;j<2*n;j++){
       backup(j)=tmp(i+1,j);
       tmp(i+1,j)=tmp(k,j);
       tmp(k,j)=backup(j);
     }
     
     d*=-1; //sign change of pfaffian
  }
  return d;
}


doublevar Pfaffian_partialpivot(const Array2 <doublevar> & in){
  //cout << "Pfaffian_partialpivot start"<<endl;
  //  calculate the pfaffian of skew-symmetric matrix
  //  partial pivoting
  assert(in.dim[0]==in.dim[1]);
  if (in.dim[0]%2!=0) return 0.0;
  int n=in.dim[0]/2;
  Array2 <doublevar> tmp(2*n,2*n);
  doublevar PF=1.0;
  doublevar fac;
  for (int i=0;i<2*n;i++)
    for (int l=0;l<2*n;l++)
      tmp(i,l)=in(i,l);

  
  
 
  int d=1;
  for (int i=0;i<2*n;i=i+2){
    //for given row look for pivoting element
    //exchange if needed
    d*=RowPivoting(tmp, i, n);
    // printout_matrix(tmp);

    for (int j=i+2;j<2*n;j++){
      fac=-tmp(i,j)/tmp(i,i+1);
      for (int k=i+1;k<2*n;k++){
        tmp(k,j)=tmp(k,j)+fac*tmp(k,i+1);
        tmp(j,k)=tmp(j,k)+fac*tmp(i+1,k);
      }
    }
    PF=PF*tmp(i,i+1);
  }
  //cout << "Pfaffian_partialpivot end"<<endl;
  return PF*d;
  
}

doublevar UpdateInversePfaffianMatrix(Array2 <doublevar> & in, 
                                      Array1 <doublevar> & row, 
                                      Array1 <doublevar> & column, 
                                      int e)
{
  //update row and collum of skew-symmetric inverse matrix
  assert(in.dim[0]==in.dim[1]);
  int n=in.dim[0]/2;

  assert(e<=2*n);

  for (int i=0;i<2*n;i++){
    column(i)=0.0;
    for (int j=0;j<2*n;j++)
      column(i)+=row(j)*in(j,i);
  }

  //to avoid the catastrophe in later division
  if (column(e)==0)
    column(e)=1e-20;

  for(int i=0;i<2*n;i++){
    if (i==e){
      in(i,i)=0.0;
    }
    else {
      in(e,i)=+in(e,i)/column(e);
      in(i,e)=-in(e,i);
    }
  }

  for(int j=0;j<2*n;j++){
    if (j!=e){
      for (int k=0;k<2*n;k++){
        if (k==j) {
          in(k,k)=0.0;
        }
        else {
          in(k,j)-=column(j)*in(k,e);
          in(j,k)=-in(k,j);
        }
      }
    }
  }

  /*
  for (int i=0;i<2*n;i++){
    for (int j=0;j<2*n;j++)
      cout << in(i,j)<< " ";
    cout <<endl;
  }
  */
  return  column(e);
}

doublevar GetUpdatedPfaffianValue(Array2 <doublevar> & in, 
                                      Array1 <doublevar> & row, 
                                      int e)
{
  //update row and collum of skew-symmetric inverse matrix
  assert(in.dim[0]==in.dim[1]);
  int n=in.dim[0]/2;

  assert(e<=2*n);

  doublevar new_val_tmp=0.0;
  for (int j=0;j<2*n;j++)
      new_val_tmp+=row(j)*in(j,e);
  

  return  new_val_tmp;
}


/*
doublevar UpdateInversePfaffianMatrix(Array2 <doublevar> & in, 
                                      Array1 <doublevar> & row, 
                                      Array1 <doublevar> & column, 
                                      int ii)
{
  //update row and collum of skew-symmetric matrix
  assert(in.dim[0]==in.dim[1]);
  int n=in.dim[0]/2;

  for (int i=0;i<2*n;i++){
    column(i)=0.0;
    for (int j=0;j<2*n;j++)
      column(i)+=row(j)*in(j,i);
  }

  if (fabs(column(ii))<1e-10) {
    cout <<"Ratio of pfaffians is essentially zero!!!\n";
    return 0.0;
  }
  else {

    for(int i=0;i<2*n;i++){
      in(ii,i)=-in(ii,i)/column(ii);
      in(i,ii)=-in(i,ii)/column(ii);
    }

    for(int j=0;j<2*n;j++){
      if (j!=ii){
        for (int k=0;k<2*n;k++){
          in(j,k)=in(j,k)+column(j)*in(ii,k);
          in(k,j)=in(k,j)+column(j)*in(k,ii);
        }
      }
    }
    return column(ii);
  }
}
*/


doublevar PfaffianInverseMatrix(const Array2 <doublevar> & a, Array2 <doublevar> & a1){
  if(!InvertPfaffianMatrix(a,a1,a.dim[0]))
    return 0.0;
  else
    return  Pfaffian_partialpivot(a);
}



//----------------------------------------------------------------------

int ludcmp(Array2 <complex  < doublevar> > & a, const int n, 
           Array1 <int> & indx, doublevar & d)
{

  const doublevar TINY=1.0e-20;
  //	cout << "ludcmp" << endl;
  Array1 <doublevar> vv(n);

  //cout << "start " << endl;
  d=1.0;
  for (int i=0;i<n;++i)
  {
    doublevar big=0.0;
    for (int j=0;j<n;++j)
    {
      //doublevar temp;
      //if ((temp=fabs(a(i,j))) > big) big=temp;
      doublevar temp=cabs(a(i,j));
      if(temp>big)
        big=temp;
    }
    if (big == 0.0)
      return 0;//error("Singular matrix in routine ludcmp");
    vv[i]=1.0/big;
  }

  //cout << "middle " << endl;
  for (int j=0;j<n;++j)
  {
    int imax;
    //cout << "j " << j << endl;
    for (int i=0;i<j;++i)
    {
      //cout << "1i " << i << endl;
      complex <double>  sum=a(i,j);
      for (int k=0;k<i;++k)
      {
        sum -= a(i,k)*a(k,j);
      }
      a(i,j)=sum;
    }
    //cout << "imax " << imax << endl;
    doublevar big(0.0);
    for (int i=j;i<n;++i)
    {
      //cout << " i " << i << endl;
      complex < double>  sum=a(i,j);
      for (int k=0;k<j;++k)
        sum -= a(i,k)*a(k,j);
      a(i,j)=sum;
      //doublevar dum;
      //cout << "dum part " << endl;
      //if ( (dum=vv[i]*fabs(sum)) >= big) {
      doublevar dum=vv[i]*cabs(sum);
      if(dum >=big)
      {
        big=dum;
        imax=i;
      }
    }
    //cout << "3 " << endl;
    if (j != imax)
    {
      for (int k=0;k<n;++k)
      {
        complex <double> dum(a(imax,k));
        a(imax,k)=a(j,k);
        a(j,k)=dum;
      }
      d = -d;
      vv[imax]=vv[j];
    }
    indx[j]=imax;
    if (a(j,j) == 0.0)
    {
      a(j,j)=TINY;
      //Write(n);
      for (int q=0;q<n;++q)
        for (int qq=0;qq<n;++qq)
          //Write3(q,qq,a(q,qq));
          //error("Singular matrix in routine ludcmp II.");
          return 0;
    }
    //cout << "5 " << endl;
    if (j != n-1)
    {
      complex < double>  dum=1.0/(a(j,j));
      for (int i=j+1;i<n;++i)
        a(i,j) *= dum;
    }
  }
  return 1;
}

//----------------------------------------------------------------------

void lubksb(Array2 < complex <doublevar> > & a, int n, 
            Array1 <int> & indx, Array1 <complex <doublevar> > & b) {
  int ii(-1);
  for (int i=0;i<n;++i)
  {
    int ip=indx[i];
    complex < doublevar > sum(b[ip]);
    b[ip]=b[i];
    if (ii>=0)
      for (int j=ii;j<i;++j)
        sum -= a(i,j)*b[j];
    else
      if (sum!=0.0)
        ii=i;
    b[i]=sum;
  }

  for (int i=n-1;i>=0;--i)
  {
    complex < doublevar > sum(b[i]);
    for (int j=i+1;j<n;++j)
      sum -= a(i,j)*b[j];
    b[i]=sum/a(i,i);
  }
}



//----------------------------------------------------------------------

log_value<dcomplex>
TransposeInverseMatrix(const Array2 < complex <doublevar> > & a, 
                       Array2 < complex <doublevar> > & a1, 
                       const int n)
{
  Array2 <complex <doublevar> >  temp(n,n);
  Array1 <int> indx(n);
  doublevar d;

  // a(i,j) first index i is row index (convention)
  // elements of column vectors are stored contiguous in memory in C style arrays
  // a(i) refers to a column vector

  // calculate the inverse of the transposed matrix because this
  // allows to pass a column vector to lubksb() instead of a row

  // put the transposed matrix in temp
  //cout << "temp " << endl;
  for(int i=0;i<n;++i)
  {
    for(int j=0;j<n;++j)
    {
      temp(i,j)=a(i,j);
      a1(i,j)=complex <doublevar> (0.0,0.0);
    }
    a1(i,i)=complex <doublevar> (1.0,0.0);
  }

  //cout << "ludcmp" << endl;
  //if the matrix is singular, the determinant is zero.
  if(ludcmp(temp,n,indx,d)==0) { 
#ifdef SUPERDEBUG
    cout << "log_value<dcomplex>TransposeInverseMatrix:zero determinant " << endl;
#endif
    return dcomplex(0.0,0.0);
  }

  //cout << "lubksb" << endl;

  for(int j=0;j<n;++j)
  {
    // get column vector
    Array1 <complex <doublevar> > yy;//(a1(j));
    yy.refer(a1(j));
    lubksb(temp,n,indx,yy);
  }


  //complex <doublevar> det(d,0);
  log_value<dcomplex> det=dcomplex(d,0);
  for(int j=0;j<n;++j) {
    det *= temp(j,j);
  }
  return det;
}

//----------------------------------------------------------------------
// Update inverse a1 after column in matrix a has changed
// get new column out of array1 newCol
//
// a1= old inverse
// newCol = new column lCol in the new matrix a_new
// returns Det(a_old)/Det(a_new)
dcomplex
InverseUpdateColumn(Array2 <complex <doublevar> > & a1, 
                    const Array1 <complex <doublevar> > & newCol,
                              const int lCol, const int n)
{
  Array1 <complex <doublevar> >  tmpColL(n);
  Array1 <complex <doublevar> >  prod(n);

  complex <doublevar> f=0.0;
  for(int i=0;i<n;++i)
  {
    f += a1(lCol,i)*newCol[i];
  }
  f =-1.0/f;

  for(int j=0;j<n;++j)
  {
    tmpColL[j]=a1(lCol,j);
    prod[j]   =0.0;
    for(int i=0;i<n;++i)
    {
      prod[j] += a1(j,i)*newCol[i];
    }
    prod[j] *= f;
  }

  for(int i=0;i<n;++i)
  {
    complex <doublevar> & p(prod[i]);
    for(int j=0;j<n;++j)
    {
      a1(i,j) += tmpColL[j]*p;
    }
  }

  f = -f;
  for(int j=0;j<n;++j)
  {
    a1(lCol,j) = f*tmpColL[j];
  }
  return f;
}


//----------------------------------------------------------------------

//----------------------------------------------------------------------



void EigenSystemSolverRealSymmetricMatrix(const Array2 <doublevar > & Ain, Array1 <doublevar> & evals, Array2 <doublevar> & evecs){
  //returns eigenvalues from largest to lowest and
  //eigenvectors, where for i-th eigenvalue, the eigenvector components are evecs(*,i)
#ifdef USE_LAPACK //if LAPACK
  int N=Ain.dim[0];
  Array2 <doublevar> Ain2(N,N);
  //need to copy the array!!
  Ain2=Ain;
  /* allocate and initialise the matrix */
  Array1 <doublevar> W, Z, WORK;
  Array1 <int> ISUPPZ, IWORK;
  int  M;
  
  /* allocate space for the output parameters and workspace arrays */
  W.Resize(N);
  Z.Resize(N*N);
  ISUPPZ.Resize(2*N);
  WORK.Resize(26*N);
  IWORK.Resize(10*N);

  int info;
  /* get the eigenvalues and eigenvectors */
  info=dsyevr('V', 'A', 'L', N, Ain2.v, N, 0.0, 0.0, 0, 0, dlamch('S'), &M,
         W.v, Z.v, N, ISUPPZ.v, WORK.v, 26*N, IWORK.v, 10*N);
  if(info>0)
    error("Internal error in the LAPACK routine dsyevr");
  if(info<0)
    error("Problem with the input parameter of LAPACK routine dsyevr in position "-info);

  for (int i=0; i<N; i++)
    evals(i)=W[N-1-i];
  for (int i=0; i<N; i++) {
    for (int j=0; j<N; j++) {
      evecs(j,i)=Z[j+(N-1-i)*N];
    }
  }
 //END OF LAPACK 
#else //IF NO LAPACK
  const int n = Ain.dim[0];
  Array2 < dcomplex > Ain_complex(n,n);
  Array2 <dcomplex> evecs_complex(n,n);
  for (int i=0; i < n; i++)
    for (int j=0; j < n; j++) {
      Ain_complex(i,j)=dcomplex(Ain(i,j),0.0);
    }
  Jacobi(Ain_complex, evals, evecs_complex);
   for (int i=0; i < n; i++)
     for (int j=0; j < n; j++){
       evecs(i,j)=real(evecs_complex(i,j));
     }
#endif //END OF NO LAPACK
}

void GeneralizedEigenSystemSolverRealSymmetricMatrices(const Array2 < doublevar > & Ain, const Array2 < doublevar> & Bin, Array1 < doublevar> & evals, Array2 < doublevar> & evecs){
  //solves generalized eigensystem problem A.x=labda*B.x
  //returns eigenvalues from largest to lowest and
  //eigenvectors, where for i-th eigenvalue, the eigenvector components are evecs(*,i)
  //eigenvectors are normalized such that: evecs**T*B*evecs = I;
#ifdef USE_LAPACK //if LAPACK
  int N=Ain.dim[0];
  
  /* allocate and initialise the matrix */
  Array2 <doublevar> A_temp(N,N), B_temp(N,N);
  Array1 <doublevar>  W,WORK;
  
  
  /* allocate space for the output parameters and workspace arrays */
  W.Resize(N);
  A_temp=Ain;
  B_temp=Bin;
  
  int info;
  int NB=64;
  int NMAX=N;
  int lda=NMAX;
  int ldb=NMAX;
  int LWORK=(NB+2)*NMAX;
  WORK.Resize(LWORK);

  /* get the eigenvalues and eigenvectors */
  info=dsygv(1, 'V', 'U' , N,  A_temp.v,  lda,  B_temp.v, ldb, W.v, WORK.v, LWORK);

  if(info>0)
    error("Internal error in the LAPACK routine dsyevr");
  if(info<0)
    error("Problem with the input parameter of LAPACK routine dsyevr in position "-info);

  for (int i=0; i<N; i++)
    evals(i)=W[N-1-i];

  for (int i=0; i<N; i++) {
    for (int j=0; j<N; j++) {
      evecs(j,i)=A_temp(N-1-i,j);
    }
  }
 //END OF LAPACK 
#else //IF NO LAPACK
  //for now we will solve it only approximatively
   int N=Ain.dim[0];
   Array2 <doublevar> B_inverse(N,N),A_renorm(N,N);
   InvertMatrix(Bin,B_inverse,N);
   MultiplyMatrices(B_inverse,Ain,A_renorm,N);
   //note A_renorm is not explicitly symmetric
   EigenSystemSolverRealSymmetricMatrix(A_renorm,evals,evecs);
#endif //END OF NO LAPACK
}

//----------------------------------------------------------------------

void GeneralizedEigenSystemSolverComplexGeneralMatrices(Array2 < dcomplex > & Ain, 
         Array1 <dcomplex> & W, Array2 <dcomplex> & VL, Array2 <dcomplex> & VR) { 
#ifdef USE_LAPACK //if LAPACK
  int N=Ain.dim[0];
  
  Array2 <dcomplex> A_temp=Ain; //,VL(N,N),VR(N,N);
  Array1 <dcomplex> WORK,RWORK(2*N);
  W.Resize(N);
  VL.Resize(N,N);
  VR.Resize(N,N);

  int info;
  int NB=64;
  int NMAX=N;
  int lda=NMAX;
  int ldb=NMAX;
  int LWORK=(NB+2)*NMAX;
  WORK.Resize(LWORK);

  //info=dsygv(1, 'V', 'U' , N,  A_temp.v,  lda,  B_temp.v, ldb, W.v, WORK.v, LWORK);

  info=  zgeev('V','V',N,A_temp.v, lda,W.v,VL.v,lda,VR.v,lda,WORK.v, LWORK,RWORK.v);

  if(info>0)
    error("Internal error in the LAPACK routine zgeev");
  if(info<0)
    error("Problem with the input parameter of LAPACK routine dsyevr in position ",-info);

//  for (int i=0; i<N; i++)
//    evals(i)=W[N-1-i];

 // for (int i=0; i<N; i++) {
 //   for (int j=0; j<N; j++) {
 //     evecs(j,i)=A_temp(N-1-i,j);
 //   }
 // }
 //END OF LAPACK 
#else //IF NO LAPACK
  error("need LAPACK for eigensystem solver for general matrices");
#endif //END OF NO LAPACK
}




void GeneralizedEigenSystemSolverRealGeneralMatrices(Array2 < doublevar > & Ain, 
         Array1 <dcomplex> & W, Array2 <doublevar> & VL, Array2 <doublevar> & VR) { 
#ifdef USE_LAPACK //if LAPACK
  int N=Ain.dim[0];
  
  Array2 <doublevar> A_temp=Ain; //,VL(N,N),VR(N,N);
  Array1 <doublevar> WORK,RWORK(2*N),WI(N),WR(N);
  WI.Resize(N);
  VL.Resize(N,N);
  VR.Resize(N,N);

  int info;
  int NB=64;
  int NMAX=N;
  int lda=NMAX;
  int ldb=NMAX;
  int LWORK=5*NMAX;
  WORK.Resize(LWORK);


  info=  dgeev('V','V',N,A_temp.v, lda,WR.v,WI.v,VL.v,lda,VR.v,lda,WORK.v,LWORK);

  if(info>0)
    error("Internal error in the LAPACK routine dgeev",info);
  if(info<0)
    error("Problem with the input parameter of LAPACK routine dgeev in position ",-info);
  W.Resize(N);
  for(int i=0; i< N; i++) { 
    W(i)=dcomplex(WR(i),WI(i));
  }

//  for (int i=0; i<N; i++)
//    evals(i)=W[N-1-i];

 // for (int i=0; i<N; i++) {
 //   for (int j=0; j<N; j++) {
 //     evecs(j,i)=A_temp(N-1-i,j);
 //   }
 // }
 //END OF LAPACK 
#else //IF NO LAPACK
  error("need LAPACK for eigensystem solver for general matrices");
#endif //END OF NO LAPACK
}
//-----------------------------------------------------------

void DGGEV(Array2 <doublevar> & A, Array2 <doublevar> & B,
           Array1 <dcomplex> & alpha, 
           Array2 <doublevar> & VL, Array2 <doublevar> & VR) {
  int N=A.dim[0];
  Array2 <doublevar> A_tmp=A;
  Array2 <doublevar> B_tmp=B;
  Array1<doublevar> alphar(N),alphai(N),beta(N);

  //Translate between row and column ordering (Fortran)
  for(int i=0; i< N; i++) { 
    for(int j=0; j< N; j++) { 
      A_tmp(i,j)=A(j,i);
      B_tmp(i,j)=B(j,i);
    }
  }
  VL.Resize(N,N);
  VR.Resize(N,N);
  alpha.Resize(N);
  int lda=N, ldb=N,LWORK=10*N,ldvl=N,ldvr=N;
  Array1 <doublevar> WORK(LWORK);
  int info;
#ifdef USE_LAPACK
  char * optionN="N";
  char * optionV="V";
  dggev_(optionV,optionV,&N,A_tmp.v,&lda,B_tmp.v,&ldb,
         alphar.v,alphai.v,beta.v,
         VL.v,&ldvl,VR.v,&ldvr,
         WORK.v,&LWORK,&info);
  alpha.Resize(N);
  for(int i=0; i< N; i++) { 
    alpha(i)=dcomplex(alphar(i)/beta(i),alphai(i)/beta(i));
  }
  
  //cout << "A " << endl;
  //for(int i=0; i < N; i++) { 
  //  for(int j=0; j< N; j++) { 
  //    cout << A(i,j) << " ";
  //  }
  //  cout << endl;
  //}
  //cout << "B " << endl;
  //for(int i=0; i < N; i++) { 
  //  for(int j=0; j< N; j++) { 
  //    cout << B(i,j) << " ";
  //  }
  //  cout << endl;
  //}
  //cout << "VR " << endl;
  //for(int i=0; i < N; i++) { 
  //  for(int j=0; j< N; j++) { 
  //    cout << VR(i,j) << " ";
  //  }
  //  cout << endl;
  //}
  
#else //IF NO LAPACK
  error("need LAPACK for eigensystem solver for general matrices");
#endif //END OF NO LAPACK

}


void DGESVD(Array2 <doublevar> & A, Array1 <doublevar> & s, 
           Array2 <doublevar> & U, Array2 <doublevar> & VT) { 

  int N=A.dim[0];
  Array2 <doublevar> A_tmp=A;

  //Translate between row and column ordering (Fortran)
  for(int i=0; i< N; i++) { 
    for(int j=0; j< N; j++) { 
      A_tmp(i,j)=A(j,i);
    }
  }
  U.Resize(N,N);
  VT.Resize(N,N);
  s.Resize(N);
  int lda=N, LWORK=10*N,ldvl=N,ldvr=N;
  Array1 <doublevar> WORK(LWORK);
  int info;
#ifdef USE_LAPACK
  char * jobU="A";
  char * jobVT="A";
  dgesvd_(jobU,jobU,&N,&N,A_tmp.v,&lda,
         s.v,U.v,&lda,
         VT.v,&lda,
         WORK.v,&LWORK,&info);

  doublevar tmp=0;
  for(int i=0;i < N; i++) { 
    for(int j=0; j< N; j++) { 
      tmp=U(i,j);
      U(i,j)=U(j,i);
      U(j,i)=tmp;
      tmp=VT(i,j);
      VT(i,j)=VT(j,i);
      VT(j,i)=tmp;
    }
  }
      
#else //IF NO LAPACK
  error("need LAPACK for eigensystem solver for general matrices");
#endif

}

//-----------------------------------------------------------


void Jacobi(const Array2 < dcomplex > & Ain, Array1 <doublevar> & evals, Array2 < dcomplex > & evecs)
{
  //assert(Ain.nr == Ain.nc);
  //assert(Ain.hermitian);

  
  //numerical recepies Jacobi transfomation eigen value solver
  const int n = Ain.dim[0];
  
  int p,q;
  doublevar fabsApq,fabsApp,fabsAqq;
  doublevar a,c,blen,sitheta;
  dcomplex b;
  doublevar t,zeta;
  doublevar alpha,stau;
  dcomplex beta,betastar;
  dcomplex  zp,zq;

  const int MAXSWEEP=100;

  /* local copy to work on */
  Array2 < dcomplex > A(n,n);

  /* Copy Ain to A */
  for (int i=0; i < n; i++)
    for (int j=0; j < n; j++) {
      A(i,j) = Ain(i,j);
      //cout <<  A(i,j).imag()<<endl;
    }

  /* Set evecs to identity */
  for (int i=0; i < n; i++)
    for (int j=0; j < n; j++) {
      if (i==j) evecs(i,j) = 1.0;
      else evecs(i,j) = 0.0;
    }

  int sweep;
  for (sweep=0; sweep < MAXSWEEP; sweep++) {
    /* S = sum of abs of real and imag of off diag terms */
    /* maxoffdiag is maximum size of off-diagonal element, and (p,q)
     * will contain its coordinates in the matrix A */
    doublevar maxoffdiag = 0.0;
    doublevar S = 0.0;
    for (int i=0; i < n; i++)
      for (int j=i+1; j < n; j++) {
        const doublevar temp = fabs(real(A(i,j)))+fabs(imag(A(i,j)));
        S += temp;
        if (temp >= maxoffdiag) {
          maxoffdiag = temp;
          p = i;
          q = j;
        }
      }
    if (maxoffdiag == 0.0) break;
    // If we're done to machine precision, go home!
    fabsApp = fabs(real(A(p,p)))+fabs(imag(A(p,p)));
    fabsAqq = fabs(real(A(q,q)))+fabs(imag(A(q,q)));
    if ( (maxoffdiag+fabsApp==fabsApp) && (maxoffdiag+fabsAqq==fabsAqq) )
      break;

    // set threshold
    doublevar thresh = 0.0;
    if (sweep < 5) thresh = 0.4*S/(n*n);

    // Loop over off diagonal terms of A in upper triangle: p < q
    for (int p=0; p < n; p++)
      for (int q=p+1; q < n; q++) {
        // If A(p,q) is too small compared to A(p,p) and A(q,q),
        // skip the rotation
        fabsApq = fabs(real(A(p,q)))+fabs(imag(A(p,q)));
        fabsApp = fabs(real(A(p,p)))+fabs(imag(A(p,p)));
        fabsAqq = fabs(real(A(q,q)))+fabs(imag(A(q,q)));
        if ( (fabsApq+fabsApp==fabsApp) && (fabsApq+fabsAqq==fabsAqq) )
          continue;
      
        // If A(p,q) is smaller than the threshold, then also skip
        // the rotation
        if (fabsApq <= thresh)
          continue;

        // the 2x2 matrix we diagonalize is [ [a b] ; [conj(b) a] ]
        a = real(A(p,p));
        c = real(A(q,q));
        b = A(p,q);
        blen = abs(b);
        zeta = (c-a)/(2.0*blen);

        // t = sgn(zeta)/(|zeta|+sqrt(zeta^2+1)), but if zeta is too
        // huge, then we set t = 1/(2*zeta)
        if ( fabs(zeta)>1.0e200 )
          t = 1/(2.0*zeta);
        else {
          t = 1.0/(fabs(zeta)+sqrt(zeta*zeta+1.0));
          if (zeta<0.0) t = -t;
        }

        /* The matrix we use to diagonalize the 2x2 block above is
         * [ [alpha beta] ; [-conj(beta) alpha] ] 
         * where alpha is real and positive and alpha = cos(theta)
         * and beta = sin(theta)*b/|b|.
         * The angle theta is chosen to diagonalize the 2x2 block.
         * The relevant formula are sin(theta)=cos(theta)*t and
         * cos(theta)=1/sqrt(1+t^2).
         * stau = (1-alpha) cleverly written. */
        alpha = 1.0/sqrt(t*t+1.0);
        sitheta = t*alpha;
        stau = sitheta*sitheta/(1.0+alpha);
        beta = b*sitheta/blen;
        betastar = conj(beta);

        /* Now we update the elements of A: */
        /* This involves chaning the p'th and q'th column of A */
        for (int i=0; i < n; i++) {
          if (i==p) {
            A(i,p) -= blen*t;
            A(i,q) = 0.0;
          }
          else if (i==q) {
            A(i,p) = 0.0;
            A(i,q) += blen*t;
          }
          else {
            zp = A(i,p);
            zq = A(i,q);
            A(i,p) -= stau*zp + betastar*zq;
            A(i,q) += beta*zp - stau*zq;
            A(p,i) = conj(A(i,p));
            A(q,i) = conj(A(i,q));
          }
        }
        /*
        for (int i=0; i < n; i++){
          for (int j=i; j < n; j++) {
            if(fabs(A(i,j).real())>1e-6)
              cout <<  A(i,j).real()<<"  ";
            else 
              cout <<"0  ";
          }
          cout << endl;
        }
        */

        /* Now we must update the eigenvector matrix with this
         * rotation:  evecs <- evecs*P_pq.
         * Update p'th and q'th column of evecs */
        for (int i=0; i < n; i++) {
          zp = evecs(i,p);
          zq = evecs(i,q);
          evecs(i,p) = alpha*zp - betastar*zq;
          evecs(i,q) = beta*zp  + alpha*zq;
        }
      } /* (p,q) rotation loop */
  } /* end of sweep loop */
    
  if (sweep == MAXSWEEP) {
    cout << endl <<"Warning:  Jacobi() needs more than "<<MAXSWEEP <<" sweeps "<<endl;
  }
  
  for (int i=0; i < n; i++) {
    evals[i] = real(A(i,i));
    //cout << evals[i]<<endl;
  }

  // sort eigs and evecs by ascending eigenvalue
  Array1 <int> list(n); 
  for (int i=0; i < n; i++) 
    list[i] = i;
  
  for (int i=1; i < n; i++) {
    const double temp = evals[i];
    int j;
    for (j=i-1; j>=0 && evals[j]<temp; j--) {
      evals[j+1] = evals[j];
      list[j+1] = list[j];
    }
    evals[j+1] = temp;
    list[j+1] = i;
  }
  
  A = evecs;
  for (int i=0; i < n; i++)
    for (int j=0; j < n; j++)
      evecs(i,j) = A(i,list[j]);
  
}
//----------------------------------------------------------------------

void sort_abs_values_descending(const Array1 <doublevar> & vals, Array1 <doublevar> & newvals, Array1 <int> & list){
  int n=vals.GetSize();
  list.Resize(n); 
  for (int i=0; i < n; i++) {
    list(i) = i;
    newvals(i)=vals(i);
  }
  
  for (int i=1; i < n; i++) {
    const double temp = newvals[i];
    const double abstemp =  fabs(newvals[i]);
    int j;
    for (j=i-1; j>=0 && fabs(newvals[j])<abstemp; j--) {
      newvals[j+1] = newvals[j];
      list[j+1] = list[j];
    }
    newvals[j+1] = temp;
    list[j+1] = i;
  }
  //for (int i=0; i < n; i++) 
  //cout << list[i]<<endl;
}
//----------------------------------------------------------------------

void sort_according_list(const Array1 <doublevar> & vals, Array1 <doublevar> & newvals, const Array1 <int> & list){
  int n=vals.GetSize();
  for (int i=0; i < n; i++) {
    newvals(i)=vals(list(i));
  }
}

//----------------------------------------------------------------------
doublevar dot(const Array1 <doublevar> & a,  const Array1 <doublevar> & b){
  int dim=a.GetSize();
  assert(a.GetSize()==b.GetSize());
  doublevar sum=0;
  for(int i=0;i<dim;i++)
    sum+=a(i)*b(i);
  return sum;
}
//----------------------------------------------------------------------

dcomplex  dot(const Array1 <doublevar> & a,  const Array1 <dcomplex> & b){
  int dim=a.GetSize();
  assert(a.GetSize()==b.GetSize());
  dcomplex sum(0.0,0.0);
  for(int i=0;i<dim;i++)
    sum+=a(i)*b(i);
  return sum;
}
//----------------------------------------------------------------------

dcomplex  dot(const Array1 <dcomplex> & a,  const Array1 <doublevar> & b){
  int dim=a.GetSize();
  assert(a.GetSize()==b.GetSize());
  dcomplex sum(0.0,0.0);
  for(int i=0;i<dim;i++)
    sum+=a(i)*b(i);
  return sum;
}

//----------------------------------------------------------------------

dcomplex  dot(const Array1 <dcomplex> & a,  const Array1 <dcomplex> & b){
  int dim=a.GetSize();
  assert(a.GetSize()==b.GetSize());
  dcomplex sum(0.0,0.0);
  for(int i=0;i<dim;i++)
    sum+=a(i)*b(i);
  return sum;
}
