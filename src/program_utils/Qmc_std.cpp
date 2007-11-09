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


#include "Qmc_std.h"

namespace global_options {
  int rappture=0;
}

mpi_info_struct mpi_info;

int parallel_sum(int inp) {
#ifdef USE_MPI
  int ret;
  MPI_Allreduce(&inp, &ret, 1,MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  return ret;
#endif
  return inp;
}

doublevar parallel_sum(doublevar inp) {
#ifdef USE_MPI
  doublevar ret;
  MPI_Allreduce(&inp, &ret, 1,MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  return ret;
#endif
  return inp;
}


dcomplex parallel_sum(dcomplex inp) { 
#ifdef USE_MPI
 doublevar real, imag, dum;
 dum=inp.real();
 MPI_Allreduce(&dum, &real,1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
 dum=inp.imag();
 MPI_Allreduce(&dum, &imag,1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
 return dcomplex(real, imag);
#endif
  return inp;
} 


int MPI_Send_complex(dcomplex & c, int node) {
#ifdef USE_MPI
  doublevar tmp=c.real();
  MPI_Send(&tmp, 1, MPI_DOUBLE, node, 0, MPI_COMM_WORLD);
  tmp=c.imag();
  MPI_Send(&tmp, 1, MPI_DOUBLE, node, 0, MPI_COMM_WORLD);
#endif
  return 1;
}
int MPI_Recv_complex(dcomplex & c , int node) { 
#ifdef USE_MPI
  doublevar real, imag;
  MPI_Status status;
  MPI_Recv(&real, 1, MPI_DOUBLE, node, 0, MPI_COMM_WORLD, &status);
  MPI_Recv(&imag, 1, MPI_DOUBLE, node, 0, MPI_COMM_WORLD, &status);
  c=dcomplex(real, imag);
#endif
  return 1; 
}



//----------------------------------------------------------------------

void wait_turn() {
#ifdef USE_MPI
  if(mpi_info.node==0) return;
  else {
    int rec=0;
    int recnode=mpi_info.node-1;
    MPI_Status status;
    MPI_Recv(&rec, 1, MPI_INT, recnode, 0, MPI_COMM_WORLD, & status);
  }
#endif
}

void finish_turn() {
#ifdef USE_MPI
  int node=mpi_info.node+1;
  if(node < mpi_info.nprocs) {
    int send=1; 
    MPI_Send(&send, 1, MPI_INT, node, 0, MPI_COMM_WORLD);
  }
#endif
}

//----------------------------------------------------------------------


//This works for unix, at least.  May not be completely portable
//to windows..can remove it if it isn't.
#include <unistd.h>
int qmcgethostname(char * name, size_t size)
{
  return gethostname(name, size);
}

//----------------------------------------------------------------------

void Terminate()
{
#ifdef USE_MPI
  cout << "Node " << mpi_info.node << " terminating " << endl;
#endif

#ifndef NO_EXCEPTIONS
     Qmc_error err_obj;
     throw err_obj;
#else
  cerr <<"Exiting now\n";
#ifdef USE_MPI
  MPI_Abort(MPI_COMM_WORLD, 75);
#endif

  exit(1);
#endif
}

//----------------------------------------------------------------------


void banner(doublevar percentage, int length, ostream & os){
  string processbar;
  processbar.insert(0, "[");
  for(int i=1;i<length+1;i++)
    processbar.insert(i, " ");
  processbar.insert(length+1, "]");
   
  int pos=int(length*percentage);
  if(pos==0){
    os <<"Starting processing..."<<endl;
  }
  for(int i=1;i<pos+1;i++)
    processbar.replace(i,1, "=");
  processbar.replace(pos+1,1, ">");
  if(pos==length){
    processbar.replace(pos+1,1, "]");
    os<<processbar<<"\n";
  }
  else
    os<<processbar<<"\r";
  os.flush();
}


//----------------------------------------------------------------------
int roundoff(double x)
{
  double y;
  y=floor(x);
  if (x>=y+0.5) return int(y+1);
  else return int(y);
}
