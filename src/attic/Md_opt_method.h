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
//------------------------------------------------------------------------
//include/Md__opt_method.h

#ifndef MD_OPT_METHOD_H_INCLUDED
#define MD_OPT_METHOD_H_INCLUDED

#include "Qmc_std.h"
#include "Qmc_method.h"
#include "Wavefunction.h"
#include "Wavefunction_data.h"
#include "Sample_point.h"
#include "Guiding_function.h"
#include "System.h"
#include "Properties.h"

class Program_options;


/*!
\brief
 Keyword: MD 
*/
class Md_opt_method : public Qmc_method
{
public:

  Md_opt_method() {
    sys=NULL;
    wfdata=NULL;
    pseudo=NULL;
    qmc_avg=NULL;
  }
  void read(vector <string> words,
            unsigned int & pos,
            Program_options & options);
  void run(Program_options & options, ostream & output);
  
  int showinfo(ostream & os);

  int generateVariables(Program_options & options);


  ~Md_opt_method()
  {
    if(sys) delete sys;
    if(wfdata) delete wfdata;

    if(pseudo) delete pseudo;

  }

private:

  void read_check(Array3 <doublevar> & last_pos);
  void write_check(Array3 <doublevar> & pos, int);

  Array1 <doublevar> atomic_weights;

  string readcheckfile, writecheckfile;

  int nstep;
  doublevar tstep;
  doublevar damp;

  Pseudopotential * pseudo;
  System * sys;
  Wavefunction_data * wfdata;

  Qmc_avg_method * qmc_avg;

  Properties_manager prop;

};

#endif //MD_OPT_METHOD_H_INCLUDED
//------------------------------------------------------------------------
