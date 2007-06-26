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

#include "Wavefunction.h"
#include "Program_options.h"
#include "MatrixAlgebra.h"

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


void Wf_return::setVals(Array2 <dcomplex> & vals, 
                        Array1 <doublevar> &  p) {

  is_complex=1;
  //cout << "wf_return:: setvals " << endl;

  //cout << "vals " << vals(0,0) << "   " << vals(0,1) << endl;
  cvals=vals;
  int ntype=vals.GetDim(1);

  Array1 <dcomplex> funcval(ntype);


  int nwf=vals.GetDim(0);
  for(int w=0; w< nwf; w++) {
    amp(w,0)=vals(w,0).real();
    phase(w,0)=p(w);

    doublevar a=exp(amp(w,0));
    
    funcval(0)=dcomplex(a*cos(p(w)), a*sin(p(w)));
    for(int i=1; i< ntype; i++) 
      funcval(i)=vals(w,i)*funcval(0);


    
    //cout << "svfuncval(0) " << funcval(0) << "  1 " << funcval(1) 
    //     << "   amp " << a << endl;
         

    doublevar dot=0;
    if(ntype>=4) {
      for(int i=1; i< 4; i++) {
        doublevar sum=funcval(i).real()+funcval(i).imag();
        amp(w,i)=sum/(a*a);
        //cout << "amp " << amp(w,i) << "  sum " << sum << endl;
        dot+=sum*sum;
      }
      //now phase:
      if(fabs(funcval(0).imag()) > 1e-8) {
        for(int i=1; i< 4; i++) {
          phase(w,i)=(amp(w,i)*funcval(0).real()/a
                               -funcval(i).real())/funcval(0).imag();
        }
      }
      else { //then we're on the real line
        for(int i=1; i< 4; i++) {
          phase(w,i)=funcval(i).imag()/funcval(0).real();
        }        
      }
    }

    if(ntype>=5) {
      amp(w,4)=(funcval(4).real()+funcval(4).imag()
        -dot/(a*a))*1/(a*a);
      phase(w,4)=0;
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


void Wf_return::mpiSend(int node) {
#ifdef USE_MPI
  int nwf, nst;
  nwf=amp.GetDim(0); nst=amp.GetDim(1);
  MPI_Send(&nwf, 1, MPI_INT, node, 0, MPI_COMM_WORLD);
  MPI_Send(&nst, 1, MPI_INT, node, 0, MPI_COMM_WORLD);
  MPI_Send(&is_complex, 1, MPI_INT, node, 0, MPI_COMM_WORLD);
  
  MPI_Send(amp.v, nwf*nst, MPI_DOUBLE, node, 0, MPI_COMM_WORLD);
  MPI_Send(phase.v, nwf*nst, MPI_DOUBLE, node, 0, MPI_COMM_WORLD);

  if(is_complex) {  
    for(int w=0; w < nwf; w++) {
      for(int i=0; i < nst; i++) {
        doublevar tmp=cvals(w,i).real();
        MPI_Send(&tmp, 1, MPI_DOUBLE, node, 0, MPI_COMM_WORLD);
        tmp=cvals(w,i).imag();
        MPI_Send(&tmp, 1, MPI_DOUBLE, node, 0, MPI_COMM_WORLD);
      }
    }
  }

#endif
  
}

void Wf_return::mpiRecieve(int node) {
#ifdef USE_MPI
  int nwf, nst;
  MPI_Status status;
  
  MPI_Recv(&nwf, 1, MPI_INT, node, 0, MPI_COMM_WORLD, &status);
  MPI_Recv(&nst, 1, MPI_INT, node, 0, MPI_COMM_WORLD, &status);
  MPI_Recv(&is_complex, 1, MPI_INT, node, 0, MPI_COMM_WORLD, &status);
  
  Resize(nwf, nst);
  MPI_Recv(amp.v, nwf*nst, MPI_DOUBLE, node, 0, MPI_COMM_WORLD, & status);
  MPI_Recv(phase.v, nwf*nst, MPI_DOUBLE, node, 0, MPI_COMM_WORLD, &status);
  if(is_complex) {  
    for(int w=0; w < nwf; w++) {
      for(int i=0; i < nst; i++) {
        doublevar tmp1, tmp2;
        MPI_Recv(&tmp1, 1, MPI_DOUBLE, node, 0, MPI_COMM_WORLD, &status);
        MPI_Recv(&tmp2, 1, MPI_DOUBLE, node, 0, MPI_COMM_WORLD, &status);
        cvals(w,i)=dcomplex(tmp1, tmp2);
      }
    }
  }  
  
#endif
}


//----------------------------------------------------------------------


//Default ForceBias to Lap methods, since that's the most
//common behavior.

void Wavefunction::updateForceBias(Wavefunction_data * wfdata,
                              Sample_point * sample) {
  updateLap(wfdata, sample);
}

void Wavefunction::getForceBias(Wavefunction_data * wfdata, int e,
                           Wf_return & bias)
{
  assert(bias.amp.GetDim(0) >=nfunc());
  assert(bias.amp.GetDim(1) >=4);

  Wf_return bias_temp(nfunc(),5);
  getLap(wfdata, e, bias_temp);
  for(int f=0; f< nfunc(); f++)  {
    for(int i=0; i< 4; i++) {
      bias.amp(f,i)=bias_temp.amp(f,i);
      bias.phase(f,i)=bias_temp.phase(f,i);
    }
  }
}

//----------------------------------------------------------------------



int deallocate(Wavefunction * & wfptr)
{
  if(wfptr == NULL)
    return 0;

  delete wfptr;
  wfptr=NULL;
  return 1;
}

int deallocate(Wavefunction_storage * & wfptr)
{
  if(wfptr == NULL)
    return 0;

  delete wfptr;
  wfptr=NULL;
  return 1;
}


//------------------------------------------------------------------------
