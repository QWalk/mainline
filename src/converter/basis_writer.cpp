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
int nfunction_symmetry(string & str) { 
  if     (str=="S") return 1;
  else if(str=="P"|| str=="P_siesta") return  3;
  //else if(types[i]=="L") tot+=4;
  else if(str=="5D" || str=="5D_siesta") return 5;
  else if(str=="6D") return 6;
  else if(str=="7F" || str=="7F_siesta") return 7;
  else if(str=="10F") return 10;
  else if(str=="15G") return 15;
   std::cout << "WARNING! Unknown type '" << str << "'" << std::endl;
  return 1;
}



int Gaussian_basis_set::nfunc() { //number of functions that this set has.
  int tot=0;
  for(unsigned int i=0; i< types.size(); i++) {
    tot+=nfunction_symmetry(types[i]);
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
void Spline_basis_writer::print_basis(ostream & os) { 
  assert(vals.size()==rad.size());
  assert(vals.size()==types.size());
  
  os << label << "\nAOSPLINE\n";
  os << endl;
  os << endl;
  int nf=vals.size();
  for(int i=0; i< nf; i++) { 
    os << "SPLINE { \n";
    os << types[i] << endl;
    int ntot=rad[i].size();
    assert(ntot==vals[i].size());
    for(int j=0; j< ntot; j++) { 
      os << rad[i][j] << "  " << vals[i][j] << endl;
    }
    os << " } ";
  }
  
}

int Spline_basis_writer::nfunc() { 
  int tot=0;
  for(unsigned int i=0; i< types.size(); i++) {
    tot+=nfunction_symmetry(types[i]);
  }
  return tot;
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
