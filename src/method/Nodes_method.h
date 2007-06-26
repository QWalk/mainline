/*
 
Copyright (C) 2007 Michal Bajdich

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


#ifndef NODES_METHOD_H_INCLUDED
#define NODES_METHOD_H_INCLUDED

#include "Array.h"
#include "MO_matrix.h"
#include "Qmc_std.h"
#include "Qmc_method.h"
#include "Wavefunction.h"
#include "Wavefunction_data.h"
#include "Sample_point.h"
#include "System.h"
class System;
class Wavefunction;
class Program_options;

class Nodes_method : public Qmc_method
{
public:

  void read(vector <string> words,
            unsigned int & pos,
            Program_options & options);
  void run(Program_options & options, ostream & output);
  int showinfo(ostream & os);
  Nodes_method() {
    mywalker=NULL;
    wfdata=NULL;
    wf=NULL;
    sysprop=NULL;
  }
  ~Nodes_method()
  {
    if(mywalker != NULL) delete mywalker;
    deallocate(wf);
    deallocate(wfdata);
    //  for(int i=0; i< nconfigsread; i++) { delete electrons(i); }
  }

private:
  Array1 <int> plots;	//3D projections of WF for particular e's  to be made
  Array1 <doublevar> minmax;	//xmin xmax ymin ymax zmin zmax
  Array2 <doublevar> dxyz;//shift of second electron to be moved
  bool doublemove; //switch for doublemove part
  doublevar resolution;	//grid coarsness: 10=coarser 0.1=finer
  string readconfig;
  int nconfigsread;
  System * sysprop;
  Array1 <Sample_point *>  tmp_configs;
  Wavefunction * wf;
  Wavefunction_data * wfdata;
  Sample_point * mywalker; //a single configuration/walker
  
 
};

#endif //NODES_METHOD_H_INCLUDED
//------------------------------------------------------------------------
