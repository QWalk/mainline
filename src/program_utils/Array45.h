/*
 
Copyright (C) 2007 Lucas K. Wagner
 extension of Arrays written by Burkhard Militzer

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
// include/Array45.h
//
//
#ifndef ARRAY45_H_INCLUDED
#define ARRAY45_H_INCLUDED

//This is somewhat hacked.  Use with caution.
#include "Array.h"

template <class T>
class Array4 {
public:
  int step1, step2, step3;
  int size;
  const int rank;
  int dim[4];
  T* v;

  Array4():step1(0),step2(0),step3(0),size(0),rank(4),v(0)
  {
    dim[0]=dim[1]=dim[2]=dim[3]=0;
  }

  ~Array4()
  {
    if (v)
      delete[] v;
  }

  Array4(int n1, int n2, int n3, int n4):           // Create without initialization
      step1(n2*n3*n4),
      step2(n3*n4),
      step3(n4),
      size(n1*step1),
      rank(4),
      v(new T[size])
  {
    dim[0] = n1;
    dim[1] = n2;
    dim[2] = n3;
    dim[3]=n4;
  }

  Array4(int n1, int n2, int n3, int n4, const T vv):size(0),rank(4),v(0)
  {    // Create with initialization
    Resize(n1,n2,n3,n4);
    for(int i=0;i<size;++i)
      v[i]=vv;
  }

  const T& operator() (int n1, int n2, int n3, int n4) const
  { // Read element
    Limits(n1,dim[0]);
    Limits(n2,dim[1]);
    Limits(n3,dim[2]);
    Limits(n4,dim[3]);
    return v[n1*step1+n2*step2+n3*step3+n4];
  }

  T& operator() (int n1, int n2, int n3, int n4)
  {       // Access element
    Limits(n1,dim[0]);
    Limits(n2,dim[1]);
    Limits(n3,dim[2]);
    Limits(n4,dim[3]);
    return v[n1*step1+n2*step2+n3*step3+n4];
  }

  /*
  Array1 <T> operator() (const int n1,const int n2, const int n3)
  {        // Access 1D sub-array
    //    cout << "Reference sub-array\n";
    Limits(n1,dim[0]);
    Limits(n2,dim[1]);
    Limits(n3,dim[2]);
    return Array1 <T> (dim[2],v+n1*step1+n2*step2+n3*step3);
  };

  Array2 <T> operator
  () (const int n1, const int n2)
  {        // Access 2D sub-array
    //    cout << "Reference sub-array\n";
    Limits(n1,dim[0]);
    Limits(n2,dim[1]);
    return Array2 <T> (dim[1],dim[2],v+n1*step1);
  };

  Array3 <T> operator
  () (const int n1)
  {        // Access 3D sub-array
    //    cout << "Reference sub-array\n";
    Limits(n1,dim[0]);
    return Array3 <T> (dim[1],dim[2],dim[3],v+n1*step1);
  };
*/

  /*
  T* operator() (const int n1, const int n2) {   // Access sub-array
    Limits(n1,dim[0]); Limits(n2,dim[1]); 
    return v+n1*step1+n2*step2;
  }
  */

  Array4 operator=(const Array4 x)
  {     // Copy array like 1D
    Resize(x.dim[0],x.dim[1],x.dim[2],x.dim[3]);
    for(int i=0;i<size;++i)
      v[i]=x.v[i];
    return *this;
  }

  Array4 & operator=(const T x)
  {        // Fill with one element
    for(int i=0;i<size;++i)
      v[i]=x;
    return *this;
  }

  void Resize(int n1, int n2, int n3, int n4)
  {
    int sizen= n1*n2*n3*n4;
    
    if (sizen>size)
    {
      if (v)
        delete[] v;
      size=sizen;
      v = new T[size];
    }
    dim[0] = n1;
    dim[1] = n2;
    dim[2] = n3;
    dim[3]=n4;
    step3=n4;
    step2=step3*n3;
    step1=step2*n2;
  }

  virtual int GetDim(const int d) const
  {
    Limits(d,rank);
    return dim[d];
  }
  virtual int* GetDims() const
  {
    return (int*)dim;
  }
  virtual int GetRank() const
  {
    return rank;
  }
  inline int GetSize() const
  {
    return size;
  }
  inline int Size() const
  {
    return size;
  }
  virtual void Read(istream& is)
  {
    error("Cannot read standard Array4\nNeed specialization of template");
  }
};

template <class T>
class Array5  {
public:
  int step1, step2, step3, step4;
  int size;
  const int rank;
  int dim[5];
  T* v;

  Array5():step1(0),step2(0),step3(0),step4(0),size(0),rank(4),v(0)
  {
    dim[0]=dim[1]=dim[2]=dim[3]=0;
  }

  ~Array5()
  {
    if (v)
      delete[] v;
  }

  Array5(int n1, int n2, int n3, int n4, int n5):           // Create without initialization
      step1(n2*n3*n4*n5),
      step2(n3*n4*n5),
      step3(n4*n5),
      step4(n5),
      size(n1*step1),
      rank(5),
      v(new T[size])
  {
    dim[0] = n1;
    dim[1] = n2;
    dim[2] = n3;
    dim[3]=n4;
    dim[4]=n5;
  }

  Array5(int n1, int n2, int n3, int n4, int n5, const T vv):size(0),rank(5),v(0)
  {    // Create with initialization
    Resize(n1,n2,n3,n4,n5);
    for(int i=0;i<size;++i)
      v[i]=vv;
  }

  const T& operator() (int n1, int n2, int n3, int n4, int n5) const
  { // Read element
    Limits(n1,dim[0]);
    Limits(n2,dim[1]);
    Limits(n3,dim[2]);
    Limits(n4,dim[3]);
    Limits(n5, dim[4]);
    return v[n1*step1+n2*step2+n3*step3+n4*step4+n5];
  }

  T& operator() (int n1, int n2, int n3, int n4, int n5)
  {       // Access element
    Limits(n1,dim[0]);
    Limits(n2,dim[1]);
    Limits(n3,dim[2]);
    Limits(n4,dim[3]);
    Limits(n5, dim[4]);
    return v[n1*step1+n2*step2+n3*step3+n4*step4+n5];
  }

  Array1 <T> operator() (const int n1,const int n2, const int n3, const int n4)
  {        // Access 1D sub-array
    //    cout << "Reference sub-array\n";
    Limits(n1,dim[0]);
    Limits(n2,dim[1]);
    Limits(n3,dim[2]);
    Limits(n4, dim[3]);
    return Array1 <T> (dim[4],v+n1*step1+n2*step2+n3*step3+n4*step4);
  };

  Array2 <T> operator
  () (const int n1, const int n2, const int n3)
  {        // Access 2D sub-array
    //    cout << "Reference sub-array\n";
    Limits(n1,dim[0]);
    Limits(n2,dim[1]);
    Limits(n3,dim[2]);
    return Array2 <T> (dim[3],dim[4],v+n1*step1+n2*step2+n3*step3);
  };

  Array3 <T>  operator() (const int n1, const int n2)
  {
    Limits(n1, dim[0]);
    Limits(n2, dim[1]);
    return Array3 <T> (dim[2], dim[3], dim[4], v+n1*step1+n2*step2);
  }

  /*
  T* operator() (const int n1, const int n2) {   // Access sub-array
    Limits(n1,dim[0]); Limits(n2,dim[1]); 
    return v+n1*step1+n2*step2;
  }
  */

  Array5 operator=(const Array5 x)
  {     // Copy array like 1D
    Resize(x.dim[0],x.dim[1],x.dim[2],x.dim[3], x.dim[4]);
    for(int i=0;i<size;++i)
      v[i]=x.v[i];
    return *this;
  }

  Array5 & operator=(const T x)
  {        // Fill with one element
    for(int i=0;i<size;++i)
      v[i]=x;
    return *this;
  }

  void Resize(int n1, int n2, int n3, int n4, int n5)
  {
    int sizen= n1*n2*n3*n4*n5;
    //    cout << "Resize dd " << sizen << " " << size << endl;
    if (sizen>size)
    {
      if (v)
        delete[] v;
      size=sizen;
      v = new T[size];
    }
    dim[0] = n1;
    dim[1] = n2;
    dim[2] = n3;
    dim[3]=n4;
    dim[4]=n5;
    step4=n5;
    step3=step4*n4;
    step2=step3*n3;
    step1 = step2*n2;
  }

  virtual int GetDim(const int d) const
  {
    Limits(d,rank);
    return dim[d];
  }
  virtual int* GetDims() const
  {
    return (int*)dim;
  }
  virtual int GetRank() const
  {
    return rank;
  }
  virtual int GetSize() const
  {
    return size;
  }
  inline int Size() const
  {
    return size;
  }
  virtual void Read(istream& is)
  {
    error("Cannot read standard Array5\nNeed specialization of template");
  }
};


#endif // ARRAY45_H_INCLUDED
//--------------------------------------------------------------------------
