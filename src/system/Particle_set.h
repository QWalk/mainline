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

#ifndef PARTICLE_SET_H_INCLUDED
#define PARTICLE_SET_H_INCLUDED
#include "Array.h"
#include "Program_options.h"



//enum particle_type {electron, ion, positron};

/*!
A set of particles, which could be electrons, ions, etc
*/
class Particle_set
{
public:

  Array1 <doublevar> charge;
  //!< Charge of each of the particles
  Array2 <doublevar> r;
  //!< position, of form (dimension, particle)

  Particle_set()
  {}
 
  Particle_set(Particle_set & x) {
    charge.Resize(x.charge.GetDim(0));
    charge=x.charge;
    r.Resize(x.r.GetDim(0), x.r.GetDim(1));
    r=x.r;
    labels=x.labels;
  }


  Particle_set & operator=(Particle_set & x) {
    charge.Resize(x.charge.GetDim(0));
    charge=x.charge;
    r.Resize(x.r.GetDim(0), x.r.GetDim(1));
    r=x.r;
    labels=x.labels;
    return *this;
  }


  int size() const
  {
    return r.GetDim(1);
  }

  //void build(Program_options & options, particle_type t);
  int read(vector <string> & words, unsigned int pos);

  string getLabel(unsigned int at) {
    assert(at < labels.size());
    return labels[at];
  }



  ~Particle_set()
  {}

  int showinfo(ostream & os);
private:
  vector <string> labels;
};

#endif // PARTICLE_SET_H_INCLUDED
//--------------------------------------------------------------------------
