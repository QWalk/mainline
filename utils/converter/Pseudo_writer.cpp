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
    os << indent << "PSEUDO { \n";
    os << indent << "  " << label << endl;
    os << indent << "  AIP 6\n";
    os << indent << "  BASIS { " << label << endl;
    os << indent << "   RGAUSSIAN \n";
    os << indent << "   OLDQMC { " << endl;
    int numL=exponents.size();
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
    os << indent << "PSEUDO { \n";
    os << indent << "  " << label << endl;
    os << indent << "  AIP 6\n";
    os << indent << "  ADD_ZEFF  #add zeff over r to the local function\n";
    os << indent << "  BASIS { \n";
    os << indent << "  " << label << endl;
    os << indent << "  AOSPLINE \n";
    string pseudofilename=label + ".qmcpseudo";
    os << indent << "  INCLUDE " << pseudofilename << "\n";
    os << indent << "  }\n";
    os << indent << "}\n\n";


    ofstream pseudofile(pseudofilename.c_str());
    pseudofile.precision(15);
    int numL=psp_pos.size();
    for(int i=0; i< numL; i++) {
      pseudofile << "SPLINE { # l= " << i << "\n";
      int npoints=psp_pos[i].size();
      pseudofile << npoints << endl;
      for(int p=0; p < npoints; p++) {
        pseudofile << "     " << psp_pos[i][p] << "  "
        << psp_val[i][p] << endl;
      }
      pseudofile << "} \n";
   }
}

//----------------------------------------------------------------------
