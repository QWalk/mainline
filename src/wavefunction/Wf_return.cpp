/*
 
 Copyright (C) 2016 Lucas K. Wagner
 
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

#include "Wf_return.h"

//----------------------------------------------------------------------

void Wf_return::read(istream & is){
  int nfunc, nst;
  string dummy;
  is >> dummy >> nfunc;
  if(dummy != "nfunc") error("expected nfunc, got ", dummy);
  is >> dummy >> nst;
  Resize(nfunc, nst);
  is >> dummy >> is_complex;
  
  is >> dummy; read_array(is, nfunc, nst, amp);
  is >> dummy; read_array(is, nfunc, nst, phase);
  is >> dummy; read_array(is, nfunc, nst, cvals);
}
//----------------------------------------------------------------------


void Wf_return::write(string & indent, ostream & os) {
  //  int is_complex;
  //  Array2 <doublevar> amp;//!< ln( |psi| ), grad ln( |psi| ), grad^2 |psi|/|psi|
  //  Array2 <doublevar> phase; //!< phase and derivatives
  //  Array2 <dcomplex> cvals; //!< (null), grad ln(|psi|), grad^2 psi/psi  for
  os << indent << "nfunc " << amp.GetDim(0) << " nst " << amp.GetDim(1) << endl;
  os << indent << "is_complex " << is_complex << endl;
  os << indent << "amp "; write_array(os, amp) ;
  os << indent << "phase "; write_array(os, phase) ;
  os << indent << "cvals "; write_array(os, cvals) ;
  
  
  
}

//----------------------------------------------------------------------




void Wf_return::setVals(Array2 <dcomplex> & vals ,
                        Array1 <doublevar> &  p) {
  
  // here we extract amplitude and phase, and their gradients and laplacians,
  // from the complex value and its gradient and derivative
  
  is_complex=1;
  cvals=vals;
  int ntype=vals.GetDim(1);
  int nwf=vals.GetDim(0);
  for (int w=0; w< nwf; w++) {
    amp(w,0)=vals(w,0).real();
    phase(w,0)=p(w);
    
    doublevar sum_ii=0;
    doublevar sum_ri=0;
    if (ntype>=4) {
      for (int i=1; i<4; i++) {
        amp(w,i)=vals(w,i).real();
        phase(w,i)=vals(w,i).imag();
        sum_ii+=phase(w,i)*phase(w,i);
        sum_ri+=amp(w,i)*phase(w,i);
      }
    }
    
    if (ntype>=5) {
      amp(w,4)=vals(w,4).real()+sum_ii;
      phase(w,4)=vals(w,4).imag()-2*sum_ri;
    }
  }
  
}

//------------------------------------------------------------------------

void Wf_return::setVals(Array2 <doublevar> & vals, Array1 <doublevar> & sign) {
  is_complex=0;
  for(int w=0; w< vals.GetDim(0); w++) {
    for(int i=0; i< vals.GetDim(1); i++) {
      cvals(w,i)=vals(w,i);
    }
  }
  amp=vals;
  
  phase=0;
  
  for(int w=0; w< sign.GetDim(0); w++) {
    phase(w,0)=.5*pi*(1-sign(w));
  }
  
}
//----------------------------------------------------------------------

void Wf_return::setVals(Array2 <log_value<doublevar> > & v ) {
  is_complex=0;
  int nfunc=v.GetDim(0);
  int nst=v.GetDim(1);
  Resize(nfunc,nst);
  for(int f=0; f< nfunc; f++) {
    amp(f,0)=v(f,0).logval;
    phase(f,0)=v(f,0).sign<0?pi:0.0;
    for(int s=1; s< nst; s++) {
      amp(f,s)=v(f,s).val();
      phase(f,s)=0.0;
    }
  }
}
//----------------------------------------------------------------------
void Wf_return::setVals(Array2 <log_value<dcomplex> > & v ) {
  is_complex=1;
  int nfunc=v.GetDim(0);
  int nst=v.GetDim(1);
  Resize(nfunc,nst);
  for(int f=0; f< nfunc; f++) {
    amp(f,0)=v(f,0).logval.real();
    phase(f,0)=v(f,0).logval.imag();
    for(int s=1; s< nst; s++) {
      amp(f,s)=v(f,s).val().real();
      phase(f,s)=v(f,s).val().imag();
    }
    if(nst > 4) {
      doublevar sum_ii=0,sum_ri=0;
      for(int s=1; s< 4; s++) {
        sum_ii+=phase(f,s)*phase(f,s);
        sum_ri+=amp(f,s)*phase(f,s);
      }
      phase(f,4)-=2*sum_ri;
      amp(f,4)+=sum_ii;
    }
    
  }
}
//----------------------------------------------------------------------


void Wf_return::mpiSend(int node) {
#ifdef USE_MPI
  int nwf, nst;
  nwf=amp.GetDim(0); nst=amp.GetDim(1);
  MPI_Send(&nwf, 1, MPI_INT, node, 0, MPI_Comm_grp);
  MPI_Send(&nst, 1, MPI_INT, node, 0, MPI_Comm_grp);
  MPI_Send(&is_complex, 1, MPI_INT, node, 0, MPI_Comm_grp);
  
  MPI_Send(amp.v, nwf*nst, MPI_DOUBLE, node, 0, MPI_Comm_grp);
  MPI_Send(phase.v, nwf*nst, MPI_DOUBLE, node, 0, MPI_Comm_grp);
  
  if(is_complex) {
    for(int w=0; w < nwf; w++) {
      for(int i=0; i < nst; i++) {
        doublevar tmp=cvals(w,i).real();
        MPI_Send(&tmp, 1, MPI_DOUBLE, node, 0, MPI_Comm_grp);
        tmp=cvals(w,i).imag();
        MPI_Send(&tmp, 1, MPI_DOUBLE, node, 0, MPI_Comm_grp);
      }
    }
  }
  
#endif
  
}

//----------------------------------------------------------------------


void Wf_return::mpiRecieve(int node) {
#ifdef USE_MPI
  int nwf, nst;
  MPI_Status status;
  
  MPI_Recv(&nwf, 1, MPI_INT, node, 0, MPI_Comm_grp, &status);
  MPI_Recv(&nst, 1, MPI_INT, node, 0, MPI_Comm_grp, &status);
  MPI_Recv(&is_complex, 1, MPI_INT, node, 0, MPI_Comm_grp, &status);
  
  Resize(nwf, nst);
  MPI_Recv(amp.v, nwf*nst, MPI_DOUBLE, node, 0, MPI_Comm_grp, & status);
  MPI_Recv(phase.v, nwf*nst, MPI_DOUBLE, node, 0, MPI_Comm_grp, &status);
  if(is_complex) {
    for(int w=0; w < nwf; w++) {
      for(int i=0; i < nst; i++) {
        doublevar tmp1, tmp2;
        MPI_Recv(&tmp1, 1, MPI_DOUBLE, node, 0, MPI_Comm_grp, &status);
        MPI_Recv(&tmp2, 1, MPI_DOUBLE, node, 0, MPI_Comm_grp, &status);
        cvals(w,i)=dcomplex(tmp1, tmp2);
      }
    }
  }  
  
#endif
}

//----------------------------------------------------------------------
