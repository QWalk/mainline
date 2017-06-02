/*
 
Copyright (C) 2007 Burkhard Militzer
 LKW--Changed RANGE_CHECKING conditional to whether or not we _call_ the
 limits functions, not what the functions do.  The other way seemed to
 confuse some complilers, notably Red Hat's gcc 2.96.

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


#ifndef ARRAY_H_INCLUDED
#define ARRAY_H_INCLUDED
#include "Qmc_std.h"
/** These are templates for matrices. We define Array1, Array2, and
    Array3, for 1, 2, and 3 dimensional arrays. */


// Limits is an inline function that checks indices
inline void Limits(const int n, const int max)
{
#ifdef RANGE_CHECKING
  if ((n<0) || (n>=max))
  {
    //    error("Array Error: Index out of range ",n,max);
    cerr << "Array error: Index out of range:  0<= "
    << n << " < " << max << "\n" ;
    Terminate();
  }
#endif // RANGE_CHECKING
}

// Limits is an inline function that checks indices
// n==max is allowed
inline void LimitsInclusive(const int n, const int max)
{
#ifdef RANGE_CHECKING
  if ((n<0) || (n>max))
  {
    //    error("Array Error: Index out of range ",n,max);
    cerr << "Array error: Upper limit for index out of range:  0<= "
    << n << " <= " << max << "\n" ;
    Terminate();
  }
#endif // RANGE_CHECKING
}

inline void EqualLimits(const int max1, const int max2)
{
#ifdef RANGE_CHECKING
  if (max1!=max2)
  {
    cerr << "Array copy error: array sizes not equal:"
    << max1 << "," << max2 << endl;
    Terminate();
  }
#endif // RANGE_CHECKING
}

inline void BiggerLimit(const int lower, const int upper)
{
#ifdef RANGE_CHECKING
  if (lower>=upper)
  {
    cerr << "Sub-array limits error: lower limit not lower "
    << lower << "," << upper << endl;
    Terminate();
  }
#endif // RANGE_CHECKING
}


// Array1 can accessed by () or []. [] keeps
// compatibility with <vector.h> form the STL
// Read is empty function since many elements do not provide
// the << operator
template <class T>
class Array1  {
public:
  int size,mSize;
  const int rank;
  T* v;
  bool b;
  Array1():size(0),mSize(0),rank(1),v(0),b(true)
  {}
  ;


  void clear() {
    if(v && b) delete [] v;
    size=0;
    mSize=0;
    v=NULL;
    b=true;
  }
  
  // Make a copy of the array and copy elements
  Array1(const Array1& a):size(a.size),mSize(a.mSize), rank(1),v(new T[mSize]),b(a.b) {
    //cout << "1";
    //cout << "  Copying array1 size="
  //<< size << "mSize=" << mSize << " " << b << endl;
    for(int i=0;i<size;++i) v[i]=a[i];
  }
  

  void refer(const Array1 & a) {
    size=a.size;
    v=a.v;
    b=false;
    mSize=a.mSize;    
  }

  // Create without initialization
  Array1(int n):size(n),mSize(n),rank(1),v(new T[n]),b(true)
  {
    //    cout << "  Array [" << n << "] of size "
    //	 << size << " x " << sizeof(T) << " allocated.\n";
  }

  // With extra memory and initialization
  Array1(int n, int mn, const T vv):size(0),mSize(0),rank(1),v(0),b(true)
  {
    cout << "  Constructing array1\n";
    Resize(mn);
    size=n;
    for(int i=0;i<size;++i)
      v[i]=vv;
  }


  // Create with initialization
  Array1(int n, const T vv):size(0),mSize(0),rank(1),v(0),b(true)
  {
    Resize(n);                        // This conflicts with Array1(int n, int mn) if T=int
    for(int i=0;i<size;++i)
      v[i]=vv;  // use Array1(int n, int mn, const T vv) instead
  }

  // Create if memory is already allocated, do not deallocate memory later
  Array1(int n, T* vv):size(n),mSize(n),rank(1),v(vv),b(false)
  {
    //    cout << "Subarray created at" << vv << endl;
  }

  ~Array1()
  {
    //cout << "Destructing";
    //      cout << "Destructor " << b << " " << v << endl;
    if (v && b)
      delete[] v;
  };

  const T& operator() (int n) const
  {     // Read element
#ifdef RANGE_CHECKING
    Limits(n,size);
#endif
    return v[n];
  };

  const T& operator[] (int n) const
  {     // Read element
#ifdef RANGE_CHECKING
    Limits(n,size);
#endif
    return v[n];
  };

  T& operator() (int n)
  {                 // Access element
#ifdef RANGE_CHECKING
    Limits(n,size);
#endif
    return v[n];
  };

  T& operator[] (int n)
  {                 // Access element
#ifdef RANGE_CHECKING
    Limits(n,size);
#endif
    return v[n];
  };

  
  Array1 operator=(const Array1 & x) {      // Copy array like 1D
    
    Resize(x.mSize);
    size=x.size;
    for(int i=0;i<size;++i) v[i]=x.v[i];
    return *this;
  }
  
  
  //Array1 & operator=(const Array1 & x)
  //{      // Copy into predefined array
//#ifdef RANGE_CHECKING
//    EqualLimits(size,x.size);
//#endif
    //    cout << "Copy array1 elements\n";
//    for(int i=0;i<size;++i)
//      v[i]=x.v[i];
//    return *this;
//  }
  

  Array1 & operator=(const T x)
  {           // Fill with one element
    for(int i=0;i<size;++i)
      v[i]=x;
    return *this;
  }


  void Resize(int n)
  {                    // Resize and do not worry
    //    cout << "Resizing " << n << " " << mSize << " " << v << endl;
    if (!b)
      error ("Cannot resize sub-array");
    if (n>mSize)
    {
      if (v)
        delete[] v;
        v = new T[n];
#ifdef NO_EXCEPTIONS
        if(v==NULL) error("Couldn't allocate Array");
#endif

      mSize=n;
    }
    size = n;
  }

  void Set(const T x)
  {                   // Set all elements, equivalent to
    for(int i=0;i<size;++i)
      v[i]=x;       // Array1=x; but no copying here!!!
  }

  void CopyResize(int n)
  {                   // Resize and copy elements
    if (!b)
      error ("Cannot copy and resize sub-array");
    if (n>mSize)
    {
      T* vv = new T[n];
      for(int i=0;i<size;i++)
        vv[i] = v[i];
      if (v)
        delete[] v;
      v = vv;
      mSize=n;
    }
    size = n;
  }

  // Resize and copy elements and reserve memory (mn) for later CopyResize
  void CopyResize(int n, int mn)
  {
    if (!b)
      error ("Cannot copy and resize sub-array");
    if (mn<n)
      mn=n;
    if (mn>mSize)
    {
      T* vv = new T[n];
#ifdef NO_EXCEPTIONS
        if(vv==NULL) error("Couldn't allocate storage to resize array");
#endif
      for(int i=0;i<size;i++)
        vv[i] = v[i];
      if (v)
        delete[] v;
      v = vv;
      mSize=mn;
    }
    size = n;
  }

   int GetDim(const int d) const
  {
#ifdef RANGE_CHECKING
    Limits(d,rank);
#endif
    return size;
  }

   int* GetDims() const
  {
    return (int*)&size;
  }
   int GetRank() const
  {
    return rank;
  }
   int GetSize() const
  {
    return size;
  }
  inline int Size() const
  {
    return size;
  }
   void Read(istream& is)
  {
    error("Cannot read standard Array1\nNeed specialization of template");
  }
  /*
  friend ostream& operator<<(ostream &os, const Array1 &a) {
      os << "Array1 [0..." << a.size-1 << "(" << a.mSize-1 << ")]= (";
      for (int i=0;i<a.size;++i) {
        if (i>0) os << ",";
        os << a[i];
      }
      os << ")";
      return os;
  }
  */
  
  
};



inline void MPI_Send(Array1 <doublevar> & arr, int node) {
#ifdef USE_MPI
  int s=arr.GetDim(0);
  MPI_Send(s, node);
  MPI_Send(arr.v,s, MPI_DOUBLE, node, 0,MPI_Comm_grp);
#endif
}
inline void MPI_Recv(Array1 <doublevar> & arr, int node) { 
#ifdef USE_MPI
  int s;
  MPI_Recv(s,node);
  arr.Resize(s);
  MPI_Status status;
  MPI_Recv(arr.v, s, MPI_DOUBLE,node,0, MPI_Comm_grp, &status);    
#endif
}


template < class T> 
void write_array(ostream & os, const Array1 <T> & a) {
  int size=a.GetDim(0);
  for(int i=0; i < size; i++) {
    os << a(i) << "   ";
  }
}

template <class T> 
void read_array(istream & is, int n, Array1 <T> & a) {
  a.Resize(n);
  for(int i=0; i< n; i++) 
    is >> a(i);
}

template <class T>
void write_json(ostream & os, const Array1 <T> & a) {
  int size=a.GetDim(0);
  os << "[";
  for(int i=0; i< size-1; i++) {
    os << a[i] << ",";
  }
  os << a[size-1] << "]";
}

//#####################################################################

template <class T>
class Array2  {
public:
  int step1;
  int size,mSize;
  const int rank;
  int dim[2];
  T* v;
  bool b;
  Array2():step1(0),size(0),mSize(0),rank(2),v(0),b(true)
  {
    dim[0]=dim[1]=0;
  };

  void clear() {
    if(v && b) delete [] v;
    size=0;
    mSize=0;
    step1=0;
    dim[0]=dim[1]=0;
    v=NULL;
    b=true;
  }


  // Create without initialization
  Array2(int n1, int n2):
      step1(n2),
      size(n1*n2),
      mSize(size),
      rank(2),
      v(new T[size]),
      b(true)
  {
   
    dim[0] = n1;
    step1=dim[1]=n2;
#ifdef NO_EXCEPTIONS
     if(v==NULL) error("Couldn't allocate Array");
#endif
    //    cout << "  Array [" << n1 << "," << n2 << "] of size "
    //	 << size << " x " << sizeof(T) << " allocated.\n";
  }
  // Create with initialization
  Array2(int n1, int n2, const T vv):size(0),mSize(0),rank(2),v(0),b(true)
  {

    Resize(n1,n2);
    for(int i=0;i<size;++i)
      v[i]=vv;
  }
  // Create if memory is already allocated, do not deallocate memory later
  Array2(int n1, int n2, T* vv):size(n1*n2),mSize(size),rank(2),v(vv),b(false)
  {

    dim[0] = n1;
    step1=dim[1]=n2;
  }
  
  Array2(const Array2 & x):step1(0),size(0),mSize(0),rank(2),v(0),b(true){

    Resize(x.dim[0], x.dim[1]);
    for(int i=0; i< size; i++) v[i]=x.v[i];
  }
  
  ~Array2() {

    if (v && b)
      delete[] v;
    //if(!b) cout << "didn't delete subarray" << endl;
  };
  const T& operator() (int n1, int n2) const
  {  // Read element
#ifdef RANGE_CHECKING
    Limits(n1,dim[0]);
    Limits(n2,dim[1]);
#endif
    return v[n1*step1+n2];
  };
  T& operator() (int n1, int n2)
  {              // Access element
#ifdef RANGE_CHECKING
    Limits(n1,dim[0]);
    Limits(n2,dim[1]);
#endif
    return v[n1*step1+n2];
  };

  /* 
  T* operator() (const int n1) {                // Access sub-array
    Limits(n1,dim[0]);
    return v+n1*step1;
  };
  */

  // Create a sub-array Array1(n1,0..<dim[1]) without allocating memory
  // returning Array1 <T> requires making a copy !!!
  Array1 <T> operator() (const int n1)
  {
    //    cout << "Reference sub-array\n";
  #ifdef RANGE_CHECKING
    Limits(n1,dim[0]);
  #endif
    return Array1 <T> (dim[1],v+n1*step1);
   }

  // Create a sub-array Array1(n1,n21..<n22) new array without allocating memory
  // returning it requires making a copy !!!
  Array1 <T> operator() (const int n1, const int n21, const int n22)
  {
    //    cout << "Reference sub-array\n";
#ifdef RANGE_CHECKING
    Limits(n1,dim[0]);
    Limits(n21,dim[1]);
    LimitsInclusive(n22,dim[1]);
    BiggerLimit(n21,n22);
#endif
    return Array1 <T> (n22-n21,v+n1*step1+n21);
  };

  // Copy array like 1D
  
  Array2 & operator=(const Array2 & x) {
    
    Resize(x.dim[0],x.dim[1]);
    for(int i=0;i<size;++i) v[i]=x.v[i];
    return *this;
  }
  

  
  // Copy array like 1D in predefined array
  //Array2 & operator=(const Array2 & x)
  //{
//#ifdef RANGE_CHECKING
 //   EqualLimits(x.dim[0],dim[0]);
 //   EqualLimits(x.dim[1],dim[1]);
//#endif
    //    cout << "Copy array2 elements " << v << endl;
 //   for(int i=0;i<size;++i)
 //     v[i]=x.v[i];
 //   return *this;
 // }
  

  Array2 & operator=(const T x)
  {                 // Fill with one element
    for(int i=0;i<size;++i)
      v[i]=x;
    return *this;
  }

  void Resize(int n1, int n2)
  {
    if (!b)
      error ("Cannot resize sub-array");
    int nn= n1*n2;
    if (nn>mSize)
    {
      if (v)
        delete[] v;
      v = new T[nn];
#ifdef NO_EXCEPTIONS
        if(v==NULL) error("Couldn't allocate Array");
#endif
      mSize=nn;
      //      cout << "Resizing 2d array (" << dim[0] << "," << dim[1] << ")->("
      //	   << n1 << "," << n2 << ")\n";
    }
    size=nn;
    dim[0] = n1;
    dim[1] = n2;
    step1 = dim[1];
  }



  void CopyResize(int n1, int n2)
  {
    if (!b)
      error ("Cannot resize sub-array");
    int nn= n1*n2;
    if ((n1!=dim[0]) || (nn>mSize))
    {
      T* vv=new T[nn];
#ifdef NO_EXCEPTIONS
        if(vv==NULL) error("Couldn't allocate Array");
#endif
      for(int i1=0;i1<n1;++i1)
      {
        for(int i2=0;i2<n2;++i2)
        {
          vv[i1*n1+i2] = v[i1*step1+i2];
        }
      }
      if (v)
        delete[] v;
      v =vv;
      mSize=nn;
      step1=dim[0]=n1;
    }
    size=nn;
    dim[1]=n2;
  }

   int GetDim(const int d) const
  {
#ifdef RANGE_CHECKING
    Limits(d,rank);
#endif
    return dim[d];
  }

   int* GetDims() const
  {
    return (int*)dim;
  }
   int GetRank() const
  {
    return rank;
  }
   int GetSize() const
  {
    return size;
  }
   void Read(istream& is)
  {
    error("Cannot read standard Array2\nNeed specialization of template");
  }
  /*
   friend ostream& operator<<(ostream &os, const Array2 &a) {
    os << "Array2 (0..." << a.dim[0]-1 << ",0..." << a.dim[1]-1 << ")= (\n";
    for (int i=0;i<a.dim[0];++i) {
      os << "   " << i << " (";
      for (int j=0;j<a.dim[1];++j) {
  if (j>0) { os << ","; }
  os << a(i,j);
      }
      os << ")\n";
    }
    os << ")\n";
    return os;
    }
  */
};


template < class T> 
void write_array(ostream & os, const Array2 <T> & a) {
  int size1=a.GetDim(0); int size2=a.GetDim(1);
  for(int i=0; i < size1; i++) {
    for(int j=0; j< size2; j++)         
      os << a(i,j) << "   ";
    os << endl;
  }
}

template <class T> 
void read_array(istream & is, int n, int m, Array2 <T> & a) {
  a.Resize(n,m);
  for(int i=0; i< n; i++) 
    for(int j=0; j< m; j++)
      is >> a(i,j);
  
}


inline void MPI_Send(Array2 <doublevar> & arr, int node) {
#ifdef USE_MPI
  int s=arr.GetDim(0);
  int v=arr.GetDim(1);
  MPI_Send(s, node);
  MPI_Send(v,node);
  MPI_Send(arr.v,s*v, MPI_DOUBLE, node, 0,MPI_Comm_grp);
#endif
}
inline void MPI_Recv(Array2 <doublevar> & arr, int node) { 
#ifdef USE_MPI
  int s,v;
  MPI_Recv(s,node);
  MPI_Recv(v,node);
  arr.Resize(s,v);
  MPI_Status status;
  MPI_Recv(arr.v, s*v, MPI_DOUBLE,node,0, MPI_Comm_grp, &status);    
#endif
}


//#####################################################################

template <class T>
class Array3  {
public:
  int step1, step2;
  int size;
  const int rank;
  int dim[3];
  T* v;
  bool b;
  Array3():step1(0),step2(0),size(0),rank(3),v(0),b(true)
  {
    dim[0]=dim[1]=dim[2]=0;
  }

  void clear() {
    if(v && b) delete [] v;
    size=0;
    step1=step2=0;
    dim[0]=dim[1]=dim[2]=0;
    v=NULL;
    b=true;
  }


  ~Array3()
  {
    if (v && b)
      delete[] v;
    //if(!b) cout << "didn't delete sub-array" << endl;
  }

  Array3(int n1, int n2, int n3):                        // Create without initialization
      step1(n2*n3),
      step2(n3),
      size(n1*step1),
      rank(3),
      v(new T[size]),
      b(true)
  {
#ifdef NO_EXCEPTIONS
        if(v==NULL) error("Couldn't allocate Array");
#endif
    dim[0] = n1;
    dim[1] = n2;
    dim[2] = n3;
  }

  Array3(int n1, int n2, int n3, const T vv):size(0),rank(3),v(0), b(true)
  {    // Create with initialization
    Resize(n1,n2,n3);
    for(int i=0;i<size;++i)
      v[i]=vv;
  }
  
  Array3(const Array3 & x):size(0), rank(3), v(0), b(true) {
    Resize(x.dim[0], x.dim[1], x.dim[2]);
    for(int i=0; i< size; i++) v[i]=x.v[i];
  }

  Array3(int n1, int n2, int n3, T* vv):size(n1*n2*n3), rank(3),
      v(vv), b(false)
  {
    dim[0]=n1;
    step1=dim[1]=n2*n3;
    step2=dim[2]=n3;
  }

  const T& operator() (int n1, int n2, int n3) const
  { // Read element
#ifdef RANGE_CHECKING
    Limits(n1,dim[0]);
    Limits(n2,dim[1]);
    Limits(n3,dim[2]);
#endif
    return v[n1*step1+n2*step2+n3];
  }

  T& operator() (int n1, int n2, int n3)
  {       // Access element
#ifdef RANGE_CHECKING
    Limits(n1,dim[0]);
    Limits(n2,dim[1]);
    Limits(n3,dim[2]);
#endif
    return v[n1*step1+n2*step2+n3];
  }

  Array1 <T> operator() (const int n1,const int n2)
  {        // Access 1D sub-array
    //    cout << "Reference sub-array\n";
#ifdef RANGE_CHECKING
    Limits(n1,dim[0]);
    Limits(n2,dim[1]);
#endif
    return Array1 <T> (dim[2],v+n1*step1+n2*step2);
  };

  Array2 <T> operator
  () (const int n1)
  {        // Access 2D sub-array
    //    cout << "Reference sub-array\n";
#ifdef RANGE_CHECKING
    Limits(n1,dim[0]);
#endif
    return Array2 <T> (dim[1],dim[2],v+n1*step1);
  };

  /*
    T* operator() (const int n1, const int n2) {   // Access sub-array
      Limits(n1,dim[0]); Limits(n2,dim[1]);
      return v+n1*step1+n2*step2;
    }
    */

  
  Array3 operator=(const Array3 & x)
  {     // Copy array like 1D
   Resize(x.dim[0],x.dim[1],x.dim[2]);
//#ifdef RANGE_CHECKING
//   assert(x.dim[0]=dim[0]);
//   assert(x.dim[1]=dim[1]);
//   assert(x.dim[2]=dim[2]);
//#endif    
   //cout << "3";
    for(int i=0;i<size;++i) {
	//cout << "i " << i << endl;
      v[i]=x.v[i];
    }
    return *this;
  }

  Array3 & operator=(const T x)
  {        // Fill with one element
    for(int i=0;i<size;++i)
      v[i]=x;
    return *this;
  }

  void Resize(int n1, int n2, int n3)
  {
    int sizen= n1*n2*n3;
    assert(sizen >= 0);
    if (sizen>size)
    {
      if (v)
        delete[] v;
      size=sizen;
      v = new T[size];
#ifdef NO_EXCEPTIONS
        if(v==NULL) error("Couldn't allocate Array");
#endif
    }
    dim[0] = n1;
    dim[1] = n2;
    dim[2] = n3;
    step2=dim[2]; 
    step1 = dim[1]*(step2);
  }

   int GetDim(const int d) const
  {
#ifdef RANGE_CHECKING
    Limits(d,rank);
#endif
    return dim[d];
  }
   int* GetDims() const
  {
    return (int*)dim;
  }
   int GetRank() const
  {
    return rank;
  }
   int GetSize() const
  {
    return size;
  }
  inline int Size() const
  {
    return size;
  }

   void Read(istream& is)
  {
    double dummy;
    is >> dummy;
    error("Cannot read standard Array3\nNeed specialization of template");
  }
};



template < class T > 
inline void array_cp(Array1 <T> & x, const Array1 <T> & y) {
  x.Resize(y.GetDim(0));
  x=y;
}

template < class T > 
inline void array_cp(Array2 <T> & x, const Array2 <T> & y) {
  x.Resize(y.GetDim(0), y.GetDim(1));
  x=y;
}

/*!
These replace the arrays with the sum over all processors.
 */

inline void parallel_sum(Array1 <doublevar> & arr) { 
#ifdef USE_MPI
  int n=arr.GetDim(0);
  Array1 <doublevar> a(n);
  MPI_Allreduce(arr.v,a.v,n,MPI_DOUBLE, MPI_SUM, MPI_Comm_grp);
  arr=a;
#endif
}

inline void parallel_sum(Array1 <doublevar> & in,Array1 <doublevar> & out) { 
#ifdef USE_MPI
  int n=in.GetDim(0);
  out.Resize(n);
  MPI_Allreduce(in.v,out.v,n,MPI_DOUBLE, MPI_SUM, MPI_Comm_grp);
#endif
}


inline void parallel_sum(Array2 <doublevar> & arr) { 
#ifdef USE_MPI
  int n=arr.GetDim(0);
  int m=arr.GetDim(1);
  Array2 <doublevar> a(n,m);
  MPI_Allreduce(arr.v,a.v,n*m,MPI_DOUBLE, MPI_SUM, MPI_Comm_grp);
  arr=a;
#endif
}



#endif // ARRAY_H_INCLUDED
//--------------------------------------------------------------------------
