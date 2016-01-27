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

#ifndef MAXIMIZE_METHOD_H_INCLUDED
#define MAXIMIZE_METHOD_H_INCLUDED

#include "Qmc_std.h"
#include "Qmc_method.h"
#include "Wavefunction.h"
#include "Wavefunction_data.h"
#include "Sample_point.h"
#include "Guiding_function.h"
#include "System.h"
#include "Pseudopotential.h"
#include "Split_sample.h"
#include "Space_warper.h"
class Program_options;
#include "Properties.h"

/*!
\brief
Find the maxima in either the wave function or the local energy.
*/
class Maximize_method :public Qmc_method {
public:

  Maximize_method() {
    guidewf=NULL;
    pseudo=NULL;
    sys=NULL;
    wfdata=NULL;
  }
  void read(vector <string> words,
            unsigned int & pos,
            Program_options & options);
  void run(Program_options & options, ostream & output);
  
  int showinfo(ostream & os);

  ~Maximize_method()  {
    if(guidewf) delete guidewf;
    if(pseudo) delete pseudo;
    if(sys) delete sys;
    if(wfdata) delete wfdata;
  }

private:
  void maximize(Sample_point * sample,Wavefunction * wf,Config_save_point & pt);
  int nconfig;
  string guidetype;
  Guiding_function * guidewf;
  Pseudopotential * pseudo;
  System * sys;
  Wavefunction_data * wfdata;
};

#endif //MAXIMIZE_METHOD_H_INCLUDED
//------------------------------------------------------------------------
