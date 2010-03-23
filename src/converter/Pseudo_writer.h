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

#ifndef PSEUDO_WRITER_H_INCLUDED
#define PSEUDO_WRITER_H_INCLUDED

#include <iostream>
#include <string>
#include <vector>


class Pseudo_writer {
  public:
    virtual void print_pseudo(std::ostream & os)=0;
    virtual ~Pseudo_writer() {}
};

class Gaussian_pseudo_writer {
public:
  std::string label;
  int atomnum;
  std::vector <std::vector <double > > exponents; //l-value, then exponents per l.
  std::vector <std::vector <double > > coefficients;
  std::vector <std::vector <int> > nvalue;  //the coefficient in r^n
  virtual void print_pseudo(std::ostream & os);
};


class Spline_pseudo_writer {
public:
  //If spin_dep is true, then the first half of L-values is the down spin,
  //and the second is the up spin.
  
  // Pseudopotential positions.  First is for l-value, second for vector..
  std::vector < std::vector <double> > psp_pos;

  // Pseudopotential values.  First for l-value.
  std::vector <std::vector <double> > psp_val;

  std::string label;
  bool spin_dep;
  virtual void print_pseudo(std::ostream & os);
  Spline_pseudo_writer() { 
    spin_dep=false;
  }
};

#endif
