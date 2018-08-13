
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

#ifndef WF_RETURN_H_INCLUDED
#define WF_RETURN_H_INCLUDED
#include "Qmc_std.h"
#include "qmc_io.h"
#include "Array.h"
#include "MatrixAlgebra.h"

/*!
 \brief
 Represent the value and derivatives of one or several many-body wave functions.
 
 This is a struct, so access the member variables amp and phase  directly.
 See documentation of those member variables for the array locations.

 The values are (somewhat confusingly) as follows:

amplitude: |Psi| Re(grad Psi/Psi) Re(lap Psi/Psi) 
phase: arg(Psi) Imag(grad Psi/Psi) Imag(lap Psi/Psi)

 */
struct Wf_return {
  Wf_return() {is_complex=0;}
  
  void write(string & indent, ostream & os);
  void read(istream & is);
  void mpiSend(int node);
  void mpiRecieve(int node);
  
  Wf_return(int nfunc, int nst) {
    is_complex=0;
    Resize(nfunc,nst);
  }
  
  void Resize(int nfunc, int nst) {
    amp.Resize(nfunc, nst);
    phase.Resize(nfunc, nst);
    cvals.Resize(nfunc, nst);
  }
  
  /*!
   Note: returns 1 if complex.
   */
  int sign(int w) {
    if(is_complex) return 1;
    
    if(fabs(phase(w,0)) < 1e-6) return 1;
    else return -1;
  }
  
  /*!
   Use log_value to set wave function values for real functions
   */
  void setVals(Array2 <log_value<doublevar> > & v );
  
  /*!
   Use log_value to set wave function values for complex functions.
   */
  void setVals(Array2 <log_value<dcomplex> > & v );
  
  /*!
   used for complex functions
   vals= [ln|psi|, grad psi/psi, grad^2 psi/psi
   p=phase
   */
  void setVals(Array2 <dcomplex> & vals, Array1 <doublevar> & p);
  
  /*!
   used for real functions
   vals= [ln|psi|, grad psi/psi, grad^2 psi/psi
   sign= sign of wave function
   */
  void setVals(Array2 <doublevar> & vals, Array1 <doublevar> &  sign);
  
  int is_complex;
  Array2 <doublevar> amp;//!< ln( |psi| ), grad ln( |psi| ), grad^2 |psi|/|psi|
  Array2 <doublevar> phase; //!< phase and derivatives
  Array2 <dcomplex> cvals; //!< (null), grad ln(|psi|), grad^2 psi/psi  for debugging purposes.
};

#endif //WF_RETURN_H_INCLUDED
