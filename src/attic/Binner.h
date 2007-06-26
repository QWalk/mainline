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
//----------------------------------------------------------------------
//include/Binner.h
//Author: Lucas Wagner
#ifndef BINNER_H_INCLUDED
#define BINNER_H_INCLUDED

#include "Qmc_std.h"
#include "Array.h"

class Sample_point;
class System;

class Binner {
 public:
  virtual void read(vector <string> &, System * )=0;
  virtual void accumulate(Sample_point *)=0;
  virtual void write()=0;
  virtual ~Binner() {}
};

class EE_bin : public Binner {
 public:
  EE_bin():nelec(0){}
  virtual void read(vector <string> &, System *  );
  virtual void accumulate(Sample_point *);
  virtual void write();
 private:
  Array1 <unsigned long int> eebin; //A count of the distances that fall in each bin
  Array1 <unsigned long int> eebin_like; //!< For electrons of like spin
  Array1 <unsigned long int> eebin_unlike; //!< For electrons of unlike spin..
  doublevar range; //The maximum range of the bins
  doublevar width; //The width of each bin
  string output_file;
  unsigned long int npoints; //!< number of points that have been accumulated
  int nelec;   //!< number of electrons
  int high_up_spin; //!< number of up spin
  
};


void allocate(vector < string> &, Binner * &, System *);

#endif //BINNER_H_INCLUDED
//----------------------------------------------------------------------
