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
#include "basis_writer.h"
#include <cassert>
#include <iomanip>
#include <cmath>

using namespace std;
//###############################################
//Basis functions
//###############################################

//----------------------------------------------------------------------

int Gaussian_basis_set::nfunc() { //number of functions that this set has.
  int tot=0;
  for(unsigned int i=0; i< types.size(); i++) {
    if     (types[i]=="S") tot+=1;
    else if(types[i]=="P") tot+=3;
    //else if(types[i]=="L") tot+=4;
    else if(types[i]=="5D") tot+=5;
    else if(types[i]=="6D") tot+=6;
    else if(types[i]=="7F") tot+=7;
    else if(types[i]=="10F") tot+=10;
    else if(types[i]=="15G") tot+=15;
    else std::cout << "WARNING! Unknown type '" << types[i] << "'" << std::endl;
  }
  return tot;
}

//----------------------------------------------------------------------

void Gaussian_basis_set::print_basis(ostream & os) {
    os << label << "\nAOSPLINE\n";
    os << endl;
    os << options;
    os << endl;
    if(fabs(cutoff) > 1e-8)
      os << "CUTOFF " << cutoff << endl;
    os << "  GAMESS { \n";
    int nmax=types.size();
    assert(nmax==exponents.size());
    assert(nmax==coefficients.size());
    for(int j=0; j< nmax; j++) {
      int kmax=exponents[j].size();
        os << types[j];
        os << "   " << kmax << endl;
        assert(kmax==coefficients[j].size());
        for(int k=0; k< kmax; k++) {
          os << "      " << setw(10) << k+1;
          os << setw(15) << exponents[j][k];
          os << setw(15) << coefficients[j][k] << endl;
        }
    }
    os << "}\n";
}

//######################################################################

int Pade_molecular_basis::nfunc(){
  return 1;
}

void Pade_molecular_basis::print_basis(ostream & os) {
  os << "LIKECUSP {\n"
        " NULL\n"
        " EXPONENTIAL_CUSP\n"
        " GAMMA 1.0\n"
        "  CUSP .25\n"
        "}\n\n"

        "UNLIKECUSP {\n"
        " NULL\n"
        " EXPONENTIAL_CUSP\n"
        " GAMMA 1.0\n"
        " CUSP .5\n"
        "}\n\n"

        "EEBASIS {\n"
        "  NULL\n"
        "  PADE\n"
        "  ALPHA0 1.2\n"
        "  NFUNC 2\n"
        "}\n\n"
        "EIBASIS {\n"
        "  NULL\n"
        "  PADE\n"
        "  ALPHA0 1.6\n"
        "  NFUNC 3\n"
        "}\n\n";
}


//####################################################################

int Exponential_cusp_basis::nfunc() {
  return 1;
}

void Exponential_cusp_basis::print_basis(ostream & os) {
  os << label << endl;
  os << " EXPONENTIAL_CUSP\n"
        " GAMMA 1.0\n"
        " CUSP 1.0\n";
}

//####################################################################

int Cutoff_cusp_basis::nfunc() {
  return 1;
}

void Cutoff_cusp_basis::print_basis(ostream & os) {
  os << label << endl;
  os << " CUTOFF_CUSP\n"
        " GAMMA 24.0\n"
        " CUSP 1.0\n"
        " CUTOFF " << cutoff << endl;

}


//####################################################################

int Poly_pade_basis::nfunc() {
  return nfunc_;
}

void Poly_pade_basis::print_basis(ostream & os) {
  os << label << endl;
  os << " POLYPADE\n"
     << " BETA0 " << beta0 << "\n"
     << " NFUNC " << nfunc_ << endl
     << " RCUT " << cutoff << endl;

}

//###################################################################
