#include <iostream>
#include "Array.h"
#include "Array45.h"

#ifndef JSONTOOLS_H_INCLUDED
#define JSONTOOLS_H_INCLUDED

template <class T> void jsonarray(std::ostream & os, Array1 <T> & arr) { 
  os << "[";
  int npts=arr.GetDim(0);
  for(int i=0; i< npts-1; i++)
    os << arr[i] << ",";
  os << arr[npts-1] << "]";
}

template <class T> void jsonarray(std::ostream & os, Array2 <T> & arr) { 
  int nptsi=arr.GetDim(0);
  int nptsj=arr.GetDim(1);
  os << "[";
  
  for(int i=0; i< nptsi; i++) { 
    os << "[";
    for(int j=0; j< nptsj-1; j++)
      os << arr(i,j) << ",";
    os << arr(i,nptsj-1) << "]";
    if(i< nptsi-1) os << ",\n";
  }
  os << "\n ] ";
}

#endif
