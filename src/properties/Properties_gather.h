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
//include/Properties_gather.h

#ifndef PROPERTIES_GATHER_H_INCLUDED
#define PROPERTIES_GATHER_H_INCLUDED

#include "Space_warper.h"
#include "Wavefunction.h"
#include "System.h"
#include "Sample_point.h"
#include "Properties_point.h"
#include "Pseudopotential.h"
#include "Guiding_function.h"
class Dynamics_generator;
class Dynamics_info;

class Properties_gather {
 public:

  Properties_gather() {
  }

  ~Properties_gather();
  void read(vector <string> & words);

  void gatherData(Properties_point &,
                  Pseudopotential *, System *, Wavefunction_data *, 
                  Wavefunction *, Sample_point *, Guiding_function *); 
 private:
};



#endif //PROPERTIES_GATHER_H_INCLUDED

//----------------------------------------------------------------------
