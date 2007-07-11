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

#ifndef BASIS_WRITER_H_INCLUDED
#define BASIS_WRITER_H_INCLUDED
#include <string>
#include <vector>
#include <iostream>


class Basis_writer {
 public:
    std::string label;
    virtual void print_basis(std::ostream & inputfile)=0;
    virtual int nfunc()=0;
    virtual ~Basis_writer() {}
};

class Gaussian_basis_set:public Basis_writer {
 public:
  //std::string label;
  std::vector <std::string> types;
  std::vector < std::vector < double > > exponents;
  std::vector < std::vector < double > > coefficients;
  std::string options;//!< any options to pass to the basis function
  double cutoff; //!< cutoff distance to enforce
  virtual void print_basis(std::ostream & inputfile);
  int nfunc();
  Gaussian_basis_set():cutoff(0.0) {}

};


class Spline_basis_writer: public Basis_writer { 
public:
  //basis function #, radial values
  std::vector < std::vector <double > > vals;
  //basis function #, r-points at which it's defined
  std::vector < std::vector <double> > rad; 
  //S,P, D, etc of the basis functions
  virtual void print_basis(std::ostream & inputfile);
  int nfunc();
  std::vector <std::string> types;
  
};

class Pade_molecular_basis:public Basis_writer {
 public:
  virtual void print_basis(std::ostream & os);
  int nfunc();
};


class Exponential_cusp_basis:public Basis_writer {
 public:
  virtual void print_basis(std::ostream & os);
  int nfunc();
};


class Cutoff_cusp_basis:public Basis_writer {
 public:
  virtual void print_basis(std::ostream & os);
  int nfunc();
  double cutoff;
  Cutoff_cusp_basis() { cutoff=5.0; label="EE"; }
};

class Poly_pade_basis:public Basis_writer {
 public:
  virtual void print_basis(std::ostream & os);
  int nfunc();
  double cutoff;
  double beta0;
  int nfunc_;
  Poly_pade_basis() { cutoff=5.0; label="EE"; }
};


#endif
