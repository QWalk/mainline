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

void Properties_point::setSize(int nwf, 
                               int n_aux, int n_aux_cvg) {

  
  kinetic.Resize(1);
  potential.Resize(nwf);
  nonlocal.Resize(nwf);
  weight.Resize(nwf);
  wf_val.Resize(nwf, 1);
  
  moved=0;
  weight=1;
  count=0;
  
  aux_energy.Resize(n_aux, n_aux_cvg);
  aux_weight.Resize(n_aux, n_aux_cvg);
  aux_energy=0;
  aux_weight=0;
  
  aux_jacobian.Resize(n_aux);
  aux_wf_val.Resize(n_aux);
  for(int i=0; i< n_aux; i++) {
    aux_wf_val(i).Resize(nwf,1);
  }
  aux_gf_weight.Resize(n_aux);
   z_pol.Resize(3);

  aux_z_pol.Resize(n_aux,3);
}


//---------------------------------------------------------------------

void Properties_point::mpiSend(int node) {
#ifdef USE_MPI
  int nwf=kinetic.GetDim(0);
  int naux=aux_energy.GetDim(0);
  int n_aux_cvg=aux_energy.GetDim(1);
  //cout << mpi_info.node << "  " << nwf << "  " << naux << endl;

  MPI_Send(&nwf, 1, MPI_INT, node, 0, MPI_COMM_WORLD);
  MPI_Send(&naux, 1, MPI_INT, node, 0, MPI_COMM_WORLD);
  MPI_Send(&n_aux_cvg, 1, MPI_INT, node, 0, MPI_COMM_WORLD);

  MPI_Send(children.v, children.GetDim(0), MPI_INT,
           node, 0,MPI_COMM_WORLD);
  MPI_Send(&nchildren, 1, MPI_INT, node, 0, MPI_COMM_WORLD);
  MPI_Send(&parent, 1, MPI_INT, node, 0, MPI_COMM_WORLD);
  MPI_Send(&count, 1, MPI_INT, node, 0, MPI_COMM_WORLD);
  MPI_Send(kinetic.v, kinetic.GetDim(0), MPI_DOUBLE, 
           node, 0, MPI_COMM_WORLD);
  MPI_Send(potential.v, potential.GetDim(0), MPI_DOUBLE,
           node, 0, MPI_COMM_WORLD);
  MPI_Send(nonlocal.v, nonlocal.GetDim(0), MPI_DOUBLE, 
           node, 0, MPI_COMM_WORLD);
  MPI_Send(weight.v, weight.GetDim(0), MPI_DOUBLE, node, 
           0, MPI_COMM_WORLD);
  MPI_Send(wf_val.amp.v, wf_val.amp.GetDim(0)*wf_val.amp.GetDim(1),
           MPI_DOUBLE, node, 0, MPI_COMM_WORLD);
  MPI_Send(wf_val.phase.v, wf_val.phase.GetDim(0)*wf_val.phase.GetDim(1),
           MPI_DOUBLE, node, 0, MPI_COMM_WORLD);
  
  // cout << mpi_info.node << ":sending aux_energy " << endl;
  MPI_Send(&gf_weight, 1, MPI_DOUBLE, node, 0, MPI_COMM_WORLD);  
  MPI_Send(aux_energy.v, aux_energy.GetDim(0)*aux_energy.GetDim(1),
           MPI_DOUBLE, node, 0, MPI_COMM_WORLD);
  MPI_Send(aux_weight.v, aux_weight.GetDim(0)*aux_weight.GetDim(1),
           MPI_DOUBLE, node, 0, MPI_COMM_WORLD);
  MPI_Send(aux_gf_weight.v, aux_gf_weight.GetDim(0),
           MPI_DOUBLE, node, 0, MPI_COMM_WORLD);
  MPI_Send(aux_jacobian.v, aux_jacobian.GetDim(0), MPI_DOUBLE,
           node, 0, MPI_COMM_WORLD); 
  for(int i=0; i< naux; i++) {
    aux_wf_val(i).mpiSend(node);
  }


  for(int d=0; d< 3; d++) 
    MPI_Send_complex(z_pol(d), node);

  for(int a=0; a< naux; a++) 
    for(int d=0; d< 3; d++) 
      MPI_Send_complex(aux_z_pol(a,d),node);

#else
    error("Properties_point::mpi_send: not using MPI,"
          " this is most likely a bug");
#endif
}

//--------------------------------------------------


void Properties_point::mpiReceive(int node) {
#ifdef USE_MPI
  MPI_Status status;

  int naux, nwf;
  int n_aux_cvg;
  MPI_Recv(&nwf, 1, MPI_INT, node, 0, MPI_COMM_WORLD, &status);
  MPI_Recv(&naux, 1, MPI_INT, node, 0, MPI_COMM_WORLD, &status);
  MPI_Recv(&n_aux_cvg, 1, MPI_INT, node, 0, MPI_COMM_WORLD, &status);
  //cout << mpi_info.node << "  " << nwf << "  " << naux << endl;

  setSize(nwf, naux, n_aux_cvg);
  MPI_Recv(children.v, children.GetDim(0), MPI_INT,
           node, 0,MPI_COMM_WORLD, &status);
  MPI_Recv(&nchildren, 1, MPI_INT, node, 0, MPI_COMM_WORLD, &status);
  MPI_Recv(&parent, 1, MPI_INT, node, 0, MPI_COMM_WORLD, &status);

  MPI_Recv(&count, 1, MPI_INT, node, 0, MPI_COMM_WORLD, &status);
  MPI_Recv(kinetic.v, kinetic.GetDim(0), MPI_DOUBLE, 
           node, 0, MPI_COMM_WORLD, & status);
  MPI_Recv(potential.v, potential.GetDim(0), MPI_DOUBLE,
           node, 0, MPI_COMM_WORLD, & status);

  MPI_Recv(nonlocal.v, nonlocal.GetDim(0), MPI_DOUBLE, 
           node, 0, MPI_COMM_WORLD, & status);
  MPI_Recv(weight.v, weight.GetDim(0), MPI_DOUBLE, node, 
           0, MPI_COMM_WORLD, & status);
  MPI_Recv(wf_val.amp.v, wf_val.amp.GetDim(0)*wf_val.amp.GetDim(1),
           MPI_DOUBLE,node, 0, MPI_COMM_WORLD, & status);
  MPI_Recv(wf_val.phase.v, wf_val.phase.GetDim(0)*wf_val.phase.GetDim(1),
           MPI_DOUBLE,node, 0, MPI_COMM_WORLD, & status);

  MPI_Recv(&gf_weight, 1, MPI_DOUBLE, node, 0, MPI_COMM_WORLD, & status);
  MPI_Recv(aux_energy.v, aux_energy.GetDim(0)*aux_energy.GetDim(1),
           MPI_DOUBLE, node, 0, MPI_COMM_WORLD, & status);
  MPI_Recv(aux_weight.v, aux_weight.GetDim(0)*aux_weight.GetDim(1),
           MPI_DOUBLE, node, 0, MPI_COMM_WORLD, & status);
  MPI_Recv(aux_gf_weight.v, aux_gf_weight.GetDim(0),
           MPI_DOUBLE, node, 0, MPI_COMM_WORLD, & status);
  MPI_Recv(aux_jacobian.v, aux_jacobian.GetDim(0),
           MPI_DOUBLE, node, 0, MPI_COMM_WORLD, & status);
  for(int i=0; i< naux; i++) {
    aux_wf_val(i).mpiRecieve(node);
  } 


  for(int d=0; d< 3; d++) 
    MPI_Recv_complex(z_pol(d), node);

  for(int a=0; a< naux; a++) 
    for(int d=0; d< 3; d++) 
      MPI_Recv_complex(aux_z_pol(a,d), node);

#else
  
  error("Properties_point::mpiRecieve: not using MPI,"
        " this is most likely a bug");
#endif
}

//----------------------------------------------------------------------


void Properties_point::read(istream & is) { 
  int nwf, naux, n_aux_cvg;
  string dummy;
  const string errmsg="Properties error in checkpoint read";
  is >> dummy >> nwf;
  if(dummy != "nwf") error("expected nwf, got ", dummy);
  is >> dummy >> naux;
  is >> dummy >> n_aux_cvg;
  setSize(nwf, naux, n_aux_cvg);
  //is >> dummy >> nchildren;
  //read_array(is,nchildren, children);
  //is >> dummy >> parent;
  //is >> dummy >> moved;
  //is >> dummy >> count;
  is >> dummy;
  read_array(is, nwf, kinetic);
  is >> dummy;
  read_array(is,nwf, potential);
  is >> dummy;
  read_array(is, nwf, nonlocal);
  is >> dummy;
  read_array(is, nwf, weight);
  is >> dummy >> dummy;  //Wf_val { 
  wf_val.read(is);
  is >> dummy; //}
  
  is >> dummy; //z_pol
  read_array(is, 3, z_pol);
  if(naux > 0) { 
    is >> dummy; //aux_energy
    read_array(is, naux, n_aux_cvg, aux_energy);
    is >> dummy; //aux_weight
    read_array(is, naux, n_aux_cvg, aux_weight);
    is >> dummy; //aux_jacobian
    read_array(is, naux, aux_jacobian);
    for(int a=0; a< naux; a++) {
      is >> dummy >> dummy; //aux_wf_val { 
      aux_wf_val(a).read(is);
      is >> dummy; // } 
    }
    is >> dummy >> gf_weight;
    is >> dummy; //aux_gf_weight
    read_array(is, naux, aux_gf_weight);

    basic_istream<char>::pos_type place=is.tellg();    
    is >> dummy; // aux_z_pol
    if(dummy != "aux_z_pol") is.seekg(place);
    else read_array(is, naux, 3, aux_z_pol);
  }

}

//----------------------------------------------------------------------


void Properties_point::write(string & indent, ostream & os) { 
  int nwf=kinetic.GetDim(0);
  int naux=aux_energy.GetDim(0);
  int n_aux_cvg=aux_energy.GetDim(1);
  os << indent << "nwf " << nwf << endl;
  os << indent << "naux " << naux << endl;
  os << indent << "n_aux_cvg " << n_aux_cvg << endl;
  
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
  
  
  os << indent << " z_pol "; write_array(os, z_pol); os << endl;
  
  if(naux > 0) { 
    os << indent << " aux_energy "; write_array(os, aux_energy); os << endl;
    os << indent << " aux_weight "; write_array(os, aux_weight); os << endl;
    os << indent << " aux_jacobian "; write_array(os, aux_jacobian); os << endl;
    
    for(int a=0; a< naux; a++) {
      os << indent << "aux_wf_val { \n";
      aux_wf_val(a).write(indent2, os);
      os << indent << "}\n";
    }
    
    os << indent << "gf_weight " << gf_weight << endl;
    os << indent << "aux_gf_weight "; write_array(os,aux_gf_weight); os << endl;
    os << indent << "aux_z_pol " ; write_array(os, aux_z_pol); os << endl;

  }
  
}

//----------------------------------------------------------------------
