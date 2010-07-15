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
#include "Pseudo_writer.h"
#include <fstream>

using namespace std;

void Gaussian_pseudo_writer::print_pseudo(ostream & os) {
    string indent="";
  os.precision(32);
    int numL=exponents.size();
    int naip=6;
    if(numL > 2) naip=12;
    os << indent << "PSEUDO { \n";
    os << indent << "  " << label << endl;
    os << indent << "  AIP " << naip << " \n";
    os << indent << "  BASIS { " << label << endl;
    os << indent << "   RGAUSSIAN \n";
    os << indent << "   OLDQMC { " << endl;
    os << indent << "0.0 " << numL << endl;
    os << indent;
    for(int i=0; i< numL ; i++) {
      int nthisl=exponents[i].size();
      os << nthisl << "  ";
    }
    os << endl;
    for(int i=0; i< numL ; i++) {
      int nthisl=exponents[i].size();
      for(int j=0; j< nthisl; j++) {
        os << indent << "    " <<  nvalue[i][j] +2 << "  "
        << exponents[i][j] << "  "
        << coefficients[i][j] << endl;
      }
    }
    os << indent << "   }\n";
    os << indent << "  }\n";
    os <<indent << "}\n\n";
}

//----------------------------------------------------------------------

void Spline_pseudo_writer::print_pseudo(ostream & os) {

  string indent = "    ";
  int numL=psp_pos.size();
  int naip=6;
  if(numL > 2) naip=12;
  if(spin_dep) numL/=2;
  
  os << indent << "PSEUDO { \n";
  os << indent << "  " << label << endl;
  os << indent << "  AIP "<< naip << " \n";
  os << indent << "  ADD_ZEFF  #add zeff over r to the local function\n";
  string filename_base=label+".qmcpseudo";

  if(spin_dep) { 
    os << indent << " SPIN_DEP \n";
    os << indent << " BASIS_DN { \n";
    os << indent << " " << label << endl;

    os << indent << "  AOSPLINE \n";
    os << indent << "  NORENORMALIZE \n";
    string filename=filename_base+"_dn";
    os << indent << "  INCLUDE " << filename << "\n";
    os << indent << "  }\n";
    ofstream pseudofile(filename.c_str());
    pseudofile.precision(15);
    
    for(int i=0; i< numL; i++) {
      pseudofile << "SPLINE { \n";
      int npoints=psp_pos[i].size();
      pseudofile << "S"  << endl;
      for(int p=0; p < npoints; p++) {
        pseudofile << "     " << psp_pos[i][p] << "  "
        << psp_val[i][p] << endl;
      }
      pseudofile << "} \n";
      //pseudofile << endl;
    }
    pseudofile.close();
    
    filename=filename_base+"_up";
    os << indent << " BASIS_UP  { \n";
    os << indent << " " << label << endl;
    os << indent << "  AOSPLINE \n";
    os << indent << "  NORENORMALIZE \n";
    os << indent << "  INCLUDE " << filename << "\n";
    os << indent << "  }\n";
    pseudofile.open(filename.c_str());
    pseudofile.precision(15);
    
    for(int i=numL; i< 2*numL; i++) {
      pseudofile << "SPLINE { \n";
      int npoints=psp_pos[i].size();
      pseudofile << "S"  << endl;
      for(int p=0; p < npoints; p++) {
        pseudofile << "     " << psp_pos[i][p] << "  "
        << psp_val[i][p] << endl;
      }
      pseudofile << "} \n";
      //pseudofile << endl;
    }
    pseudofile.close();
    
  }
  else { 
    os << indent << "  BASIS { \n";
    os << indent << "  " << label << endl;
    os << indent << "  AOSPLINE \n";
    os << indent << "  NORENORMALIZE \n";
    os << indent << "  INCLUDE " << filename_base << "\n";
    os << indent << "  }\n";
    ofstream pseudofile(filename_base.c_str());
    pseudofile.precision(15);
    
    for(int i=0; i< numL; i++) {
      pseudofile << "SPLINE { \n";
      int npoints=psp_pos[i].size();
      pseudofile << "S"  << endl;
      for(int p=0; p < npoints; p++) {
        pseudofile << "     " << psp_pos[i][p] << "  "
        << psp_val[i][p] << endl;
      }
      pseudofile << "} \n";
    }
    
  }
  os << indent << "}\n\n";
  
  
}

//----------------------------------------------------------------------
