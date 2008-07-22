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

#include "Sample_point.h"
#include "Wavefunction.h"
void Sample_point::attachObserver(Wavefunction * wfptr)
{
  assert(wfptr != NULL);
  wfObserver=wfptr;
  wfptr->notify(sample_attach, 0);
}


int Sample_point::ndim() { return 3; }

//----------------------------------------------------------------------

int read_config(string & last_read, istream & is, 
                Sample_point * sample) {
  if(last_read=="CONFIGS") {
    is >> last_read;
    if(last_read != "{") error("Need a { after CONFIGS");
    sample->rawInput(is);
    is >> last_read;
    if(last_read != "}") error("Need a closing }  for CONFIGS");
    return 1;
  }
  else 
    return 0;
}


void write_config(ostream & os, 
                  Sample_point * sample) {
  os << "   CONFIGS { \n";
  sample->rawOutput(os);
  os << "   }\n";
}


//-------------------------------------------------------------------------


void Config_save_point::savePos(Sample_point * sample) {
  
  int nelectrons=sample->electronSize();
  
  if(electronpos.GetDim(0)!=nelectrons) {
    electronpos.Resize(nelectrons);
    for(int i=0; i< nelectrons; i++) {
      electronpos(i).Resize(3);
  
    }
  }
  
  
  for(int i=0; i< nelectrons; i++) {
    sample->getElectronPos(i,electronpos(i));
  }
  
}

void Config_save_point::restorePos(Sample_point * sample) {
  int nelectrons=sample->electronSize();
  assert(nelectrons==electronpos.GetDim(0));
  for(int i=0; i< nelectrons; i++) 
    sample->setElectronPos(i,electronpos(i));
}

void Config_save_point::mpiReceive(int node) {
#ifdef USE_MPI
  MPI_Status status;
  int nelectrons;
  MPI_Recv(&nelectrons,1, MPI_INT, node, 0, MPI_Comm_grp,
           &status);
  electronpos.Resize(nelectrons);
  for(int e=0; e< nelectrons; e++) {
    electronpos(e).Resize(3);
    for(int d=0; d< 3; d++) {
      MPI_Recv(&(electronpos(e)(d)), 1, MPI_DOUBLE,
               node, 0, MPI_Comm_grp, &status);
    }
  }
#endif
}

void Config_save_point::mpiSend(int node) {
#ifdef USE_MPI
  int nelectrons=electronpos.GetDim(0);
  MPI_Send(&nelectrons,1, MPI_INT, node, 0, MPI_Comm_grp);
  for(int e=0; e< nelectrons; e++) {
    for(int d=0; d< 3; d++) {
      MPI_Send(&(electronpos(e)(d)), 1, MPI_DOUBLE,
                 node, 0, MPI_Comm_grp);
    }
  }
#endif
}

//----------------------------------------------------------------------
