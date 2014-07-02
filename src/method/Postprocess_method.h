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
#include "Properties.h"
class System;
class Wavefunction;
class Program_options;

/*!
In DMC, use the option SAVE_TRACE filename to save the walker positions and their weights to a binary file (called 'filename') once every block.
These walkers should be decorrelated, so they can now be treated independently.

Now, you can use 
method { POSTPROCESS
  READCONFIG filename
  AVERAGE { ... } 
  DENSITY { ... }
  }
  include sys
  trialfunc { .. }

to take any averages that you'd like. 

There are some tradeoffs with this method. First, since we only save the walkers occasionally, the error bars will b
e larger than if you had put AVERAGE into the DMC section directly. However, if the AVERAGE takes some time to evalu
ate, such as the 2-RDM, then you will also not waste time evaluating the object on highly correlated walkers. You sh
ould use this scheme when one of the following is true:
a) You're not sure what averaging variables you want when the DMC calculation is run.
b) The averaging variables you want dominate the cost of the DMC calculation (2-RDM!)

Otherwise, evaluation on the fly (AVERAGE in the DMC section) is better.
 
 */
class Postprocess_method : public Qmc_method
{
public:
  void read(vector <string> words,
            unsigned int & pos,
            Program_options & options);
  void run(Program_options & options, ostream & output);
  int showinfo(ostream & os);
  Postprocess_method():sys(NULL),pseudo(NULL),wfdata(NULL) {}
  ~Postprocess_method()
  {
    if(sys != NULL) delete sys;
  }

private:
  System * sys;
  Pseudopotential * pseudo;
  Wavefunction_data * wfdata;
  Array1 < Local_density_accumulator *> densplt;
  Array1 < Average_generator * > average_var;
  string configfile;
};

#endif //POSTPROCESS_METHOD_H_INCLUDED
//------------------------------------------------------------------------
