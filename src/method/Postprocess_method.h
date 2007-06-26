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


#ifndef POSTPROCESS_METHOD_H_INCLUDED
#define POSTPROCESS_METHOD_H_INCLUDED

#include <vector>
#include <fstream>
#include <math.h>
#include <stdlib.h>
#include "Qmc_method.h"
#include "qmc_io.h"
#include "Array.h"
#include "MO_matrix.h"
#include "Sample_point.h"
class System;
class Wavefunction;
class Program_options;

/*!
 */
class Postprocess_method : public Qmc_method
{
public:
  void read(vector <string> words,
            unsigned int & pos,
            Program_options & options);
  void run(Program_options & options, ostream & output);
  int showinfo(ostream & os);
  Postprocess_method():sample(NULL),sys(NULL) {}
  ~Postprocess_method()
  {
    if(sample != NULL) delete sample;
    if(sys != NULL) delete sys;
  }

private:
  Sample_point * sample; //a single configuration/walker
  System * sys;
  string configfile;
};

#endif //POSTPROCESS_METHOD_H_INCLUDED
//------------------------------------------------------------------------
