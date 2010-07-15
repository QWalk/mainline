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

#ifndef QMC_STD_INCLUDED
#define QMC_STD_INCLUDED

#ifdef USE_MPI
#include "mpi.h"
#endif

#include <assert.h>
#include <vector>
#include <string>
#include <cstdio>
#include <complex>

#ifdef USE_BLAS
extern "C" {
#include "cblas.h"
}
#endif

using namespace std;

typedef double doublevar;
typedef complex <doublevar> dcomplex;
const doublevar zero_exp=-1e3; //!< effective zero for an exponent



class mpi_info_struct
{
public:
  mpi_info_struct()
  {
    node=0;
    nprocs=1;
  }
  int node;
  int nprocs;
};
extern mpi_info_struct mpi_info;

#ifdef USE_MPI
extern MPI_Comm MPI_Comm_grp;  // communicator for each independent process
#endif

int parallel_sum(int inp);
doublevar parallel_sum(doublevar inp);
dcomplex parallel_sum(dcomplex inp);

int MPI_Send_complex(dcomplex & , int node);
int MPI_Recv_complex(dcomplex &, int node);
int MPI_Send(double &, int node);
int MPI_Recv(double &, int node);
int MPI_Send(int &, int node);
int MPI_Recv(int &, int node);



namespace global_options {
  extern int rappture;
}

class Program_options;

int qmcgethostname(char *, size_t);


class Qmc_error {

};


//put wait_turn() and finish_turn() around a block
//to make the processes go in order from 0,..,n-1
//Useful for I/O through NFS servers(especially writing configs)
void wait_turn();
void finish_turn();



template <class T>
void exchange(T & a, T & b) {
  T tmp;
  tmp=a;
  a=b;
  b=tmp;
}


#include <iostream>
#include <fstream>
#include <stdlib.h>


#include <math.h>
using namespace std;

#ifndef pi
const double pi=3.1415926535897932385;
#endif

#ifndef I
const dcomplex I(0.0,1.0);
#endif

//--------------------------------------------------------------------------
// include/Burk_standard.h
//
//
// Standard definitions and functions
//
// Burkhard Militzer                                    Urbana 4-9-99
// edited by Lucas Wagner, NCSU, 2002


inline double sign(double x)
{
  return (x>0.0) ? 1.0 : ((x<0.0) ? -1.0 : 0.0);
}

inline int sign(int x)
{
  return (x>0) ? 1 : ((x<0) ? -1 : 0);
}

void Terminate();

inline void error(const char* m)
{
  cout<<"Error   " << m << "\n";
  Terminate();
}

template<class T>
inline
void error(const char* m, const T& n)
{
  cout.precision(16);
  cout<<"Error   "<< m << n <<"\n";
  Terminate();
}

template<class T, class U>
inline
void error(const char* m, const T& t, const U& u)
{
  cout.precision(16);
  cout<<"Error   "<< m << t << u <<"\n";
  Terminate();
}

template<class T, class U, class V>
inline
void error(const char* m, const T& t, const U& u, const V& v)
{
  cout.precision(16);
  cout<<"Error   "<< m << t << u << v <<"\n";
  Terminate();
}

template<class T, class U, class V, class W>
inline
void error(const char* m, const T& t, const U& u, const V& v, const W& w)
{
  cout.precision(16);
  cout<<"Error   "<< m << t << u << v << w << "\n";
  Terminate();
}

int roundoff(double x);

//print out a progress bar
void banner(doublevar percentage, int length, ostream & os);


/*!
This is a slightly tricky bit that I'd like to get rid of if possible.
When we have separate data and calculator wavefunction objects, the
calculator needs to use a pointer to the derived class to access the
data members of the data class, but of course the higher level parts
of the program shouldn't know about it.  So we have to recast the
base pointer to a derived pointer.
*/
template <class T, class U>
void recast(T * & baseptr, U * & derivedptr)
{
  derivedptr = dynamic_cast <U *> (baseptr);
  if(derivedptr == 0)
  {
    error("Died while attempting to downcast a variable.  There's probably"
          "something wrong with the code.  Make sure you're using recast() on"
          "a pointer to a base class and a derived pointer, and that everything"
          "is getting generated correctly.");
  }
}
#include "Array.h"



#endif //QMC_STD_INCLUDED

//--------------------------------------------------------------------------
