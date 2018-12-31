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
#include "Properties_point.h"

void Properties_point::setSize(int nwf) {
  kinetic.Resize(nwf);
  potential.Resize(nwf);
  nonlocal.Resize(nwf);
  weight.Resize(nwf);
  wf_val.Resize(nwf, 1);
  weight=1;
  count=0;
  kinetic=0;
  potential=0;
  nonlocal=0;

}


//---------------------------------------------------------------------

void Properties_point::mpiSend(int node) {
#ifdef USE_MPI
  int nwf=kinetic.GetDim(0);

  MPI_Send(&nwf, 1, MPI_INT, node, 0, MPI_Comm_grp);
  MPI_Send(&count, 1, MPI_INT, node, 0, MPI_Comm_grp);
  MPI_Send(kinetic.v, kinetic.GetDim(0), MPI_DOUBLE, 
           node, 0, MPI_Comm_grp);
  MPI_Send(potential.v, potential.GetDim(0), MPI_DOUBLE,
           node, 0, MPI_Comm_grp);
  MPI_Send(nonlocal.v, nonlocal.GetDim(0), MPI_DOUBLE, 
           node, 0, MPI_Comm_grp);
  MPI_Send(weight.v, weight.GetDim(0), MPI_DOUBLE, node, 
           0, MPI_Comm_grp);
  MPI_Send(wf_val.amp.v, wf_val.amp.GetDim(0)*wf_val.amp.GetDim(1),
           MPI_DOUBLE, node, 0, MPI_Comm_grp);
  MPI_Send(wf_val.phase.v, wf_val.phase.GetDim(0)*wf_val.phase.GetDim(1),
           MPI_DOUBLE, node, 0, MPI_Comm_grp);

  
  int ni=avgrets.GetDim(0);
  int nj=avgrets.GetDim(1);
  MPI_Send(ni,node);
  MPI_Send(nj,node);
  for(int i=0; i< ni; i++) { 
    for(int j=0; j< nj; j++) { 
      MPI_Send(avgrets(i,j).vals,node);
      MPI_Send(avgrets(i,j).type,node);
    }
  }
#else
    error("Properties_point::mpi_send: not using MPI,"
          " this is most likely a bug");
#endif
}

//---------------------------------------------------------------------

void Properties_point::mpiReceive(int node) {
#ifdef USE_MPI
  MPI_Status status;

  int nwf;
  MPI_Recv(&nwf, 1, MPI_INT, node, 0, MPI_Comm_grp, &status);
  //cout << mpi_info.node << "  " << nwf << endl;

  setSize(nwf);
  MPI_Recv(&count, 1, MPI_INT, node, 0, MPI_Comm_grp, &status);
  MPI_Recv(kinetic.v, kinetic.GetDim(0), MPI_DOUBLE, 
           node, 0, MPI_Comm_grp, & status);
  MPI_Recv(potential.v, potential.GetDim(0), MPI_DOUBLE,
           node, 0, MPI_Comm_grp, & status);

  MPI_Recv(nonlocal.v, nonlocal.GetDim(0), MPI_DOUBLE, 
           node, 0, MPI_Comm_grp, & status);
  MPI_Recv(weight.v, weight.GetDim(0), MPI_DOUBLE, node, 
           0, MPI_Comm_grp, & status);
  MPI_Recv(wf_val.amp.v, wf_val.amp.GetDim(0)*wf_val.amp.GetDim(1),
           MPI_DOUBLE,node, 0, MPI_Comm_grp, & status);
  MPI_Recv(wf_val.phase.v, wf_val.phase.GetDim(0)*wf_val.phase.GetDim(1),
           MPI_DOUBLE,node, 0, MPI_Comm_grp, & status);
  int ni,nj;
  MPI_Recv(ni,node);
  MPI_Recv(nj,node);
  avgrets.Resize(ni,nj);
  for(int i=0; i< ni; i++) { 
    for(int j=0; j< nj; j++) { 
      MPI_Recv(avgrets(i,j).vals,node);
      MPI_Recv(avgrets(i,j).type,node);
    }
  }

#else
  
  error("Properties_point::mpiRecieve: not using MPI,"
        " this is most likely a bug");
#endif
}

//----------------------------------------------------------------------

#include "qmc_io.h"

void Properties_point::read(istream & is) { 
  int nwf;
  string dummy;
  const string errmsg="Properties error in checkpoint read";
  is >> dummy >> nwf;
  if(dummy != "nwf") error("expected nwf, got ", dummy);
  is >> dummy;
  if(dummy=="naux") {
    debug_write(cout,"Trying to read from old Properties_point..");
    int naux;
    is >> naux;
    if(naux > 0) error("Don't support auxiliary wave functions any more!");
    is >> dummy >> dummy;
  }
  setSize(nwf);
  //is >> dummy;
  if(dummy!="kinetic") error("expected kinetic, got ",dummy);
  read_array(is, nwf, kinetic);
  is >> dummy;
  read_array(is,nwf, potential);
  is >> dummy;
  read_array(is, nwf, nonlocal);
  is >> dummy;
  if(dummy!="weight") error("expected weight, got ",dummy);
  read_array(is, nwf, weight);
  is >> dummy >> dummy;  //Wf_val { 
  wf_val.read(is);
  is >> dummy; //}
}

//----------------------------------------------------------------------


void Properties_point::write(string & indent, ostream & os) { 
  int nwf=kinetic.GetDim(0);
  os << indent << "nwf " << nwf << endl;
  os << indent << "kinetic ";
  write_array(os, kinetic);
  os << endl << indent << "potential ";
  write_array(os, potential);
  os << endl << indent << "nonlocal ";
  write_array(os, nonlocal);
  os << endl << indent << "weight ";
  write_array(os, weight);
  os << endl;

  os << indent << "Wf_val { \n";
  string indent2=indent+"  ";
  wf_val.write(indent2, os);
  os << indent << "}\n";
}

//----------------------------------------------------------------------
void Properties_point::weighted_add(const Properties_point & pt) { 
  int nwf=kinetic.GetDim(0);
  assert(nwf==pt.kinetic.GetDim(0));
  for(int w=0;w < nwf; w++) {
    kinetic(w)+=pt.weight(w)*pt.kinetic(w);
    potential(w)+=pt.weight(w)*pt.potential(w);
    nonlocal(w)+=pt.weight(w)*pt.nonlocal(w);
    weight(w)+=pt.weight(w);
    int navg=avgrets.GetDim(1);
    assert(navg==pt.avgrets.GetDim(1));
    for(int a=0; a< navg; a++) { 
      for(int j=0; j< pt.avgrets(w,a).vals.GetDim(0); j++) { 
        avgrets(w,a).vals(j)+=pt.weight(w)*pt.avgrets(w,a).vals(j);
      }
    }
  }
}


//----------------------------------------------------------------------

void Properties_point::unweighted_add(const Properties_point & pt,doublevar pre) { 
  int nwf=kinetic.GetDim(0);
  assert(nwf==pt.kinetic.GetDim(0));
  for(int w=0;w < nwf; w++) {
    kinetic(w)+=pre*pt.kinetic(w);
    potential(w)+=pre*pt.potential(w);
    nonlocal(w)+=pre*pt.nonlocal(w);
    weight(w)+=pre*pt.weight(w);
    int navg=avgrets.GetDim(1);
    assert(navg==pt.avgrets.GetDim(1));
    for(int a=0; a< navg; a++) { 
      for(int j=0; j< pt.avgrets(w,a).vals.GetDim(0); j++) { 
        avgrets(w,a).vals(j)+=pre*pt.avgrets(w,a).vals(j);
      }
    }
  }
}
//----------------------------------------------------------------------

