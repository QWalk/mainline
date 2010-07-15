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

#ifndef CONVERTER_H_INCLUDED
#define CONVERTER_H_INCLUDED
#include <string>
#include <vector>
#include <iostream>
#include <cassert>
#include <complex>
#include <cstdio>

typedef double doublevar;
typedef std::complex <doublevar> dcomplex;

//http://www.mindcracker.com/mindcracker/c_cafe/stl/PartingSTLStrings.asp
void split
(std::string & text, std::string & separators, std::vector<std::string> & words);


void append_number(std::string & str, int num);


class Atom {
public:
  std::vector <double> pos;
  std::string name;
  double charge;
  int basis;  //basis set number
  Atom() {
    pos.push_back(0);
    pos.push_back(0);
    pos.push_back(0);
  }
  void print_atom(std::ostream & os) {
    os.precision(10);
    os << "ATOM { " << name << "  "
       << charge << "  COOR "
       << pos[0] << "   " << pos[1]
       << "   " << pos[2] << " } " << std::endl;
  }
};

void find_unique_atoms(const std::vector<Atom> & atoms, std::vector<std::string> & unique_atoms);

class Center {
public:
  std::vector <double> pos;
  std::string name;
  int equiv_atom; //equivalent atom
  int basis;  //basis set#
  Center() {
    for(int i=0; i< 3; i++) pos.push_back(0);
  }
  void print_center(std::ostream & os) {
    os << name << "   ";
    for(int i=0; i< 3; i++) os << pos[i] << "  ";
    os << std::endl;
  }
};


void print_orbitals(std::ostream & os,
                    std::vector <Center> & centers,
                    std::vector <int> & nbasis,
                    std::vector < std::vector < double > > & mo_coeff);

void print_orbitals(std::ostream & os,
                    std::vector <Center> & centers,
                    std::vector <int> & nbasis,
                    std::vector < std::vector < dcomplex > > & mo_coeff);



void print_centers(std::ostream & os, std::vector <Center> & centers);

void print_vmc_section(std::ostream & os,
                       std::string & outputname, double eref);

void print_opt_section(std::ostream & os,
                       std::string & outputname, double eref);

void print_dmc_section(std::ostream & os,
                       std::string & outputname, double eref);


double find_centers(std::vector <double> & origin,
                    std::vector <std::vector <double> > & latvec,
                    std::vector <Atom> & atoms,
                    std::vector <Center> & centers );


bool compare_mo(std::vector <std::vector < double> > & oldMOCoeff,
                std::vector <std::vector < double> > & moCoeff,
                std::vector <int> & compare_list  );

double find_basis_cutoff(std::vector <std::vector <double> > & latvec);

const double pi=3.1415926535897932385;


#endif

//#ifndef uint
//#define uint unsigned int
//#endif



