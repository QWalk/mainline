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
//include/Fitbasis_method.h

#ifndef FITBASIS_METHOD_H_INCLUDED
#define FITBASIS_METHOD_H_INCLUDED

#include "Qmc_std.h"
#include "Qmc_method.h"
#include "Basis_function.h"
#include "Center_set.h"
class Program_options;
/*!

 */
class Fitbasis_method : public Qmc_method
{
public:
  void read(vector <string> words,
            unsigned int & pos,
            Program_options & options);
  void run(Program_options & options, ostream & output);
  Fitbasis_method() {
    sys=NULL;
  }
  ~Fitbasis_method()
  {
    if(sys) delete sys;
    for(int i=0; i< basis.GetDim(0); i++) {
      if(basis(i)) delete basis(i);
    }
  }

  int showinfo(ostream & os);
private:

  Array1 <Basis_function * > basis;
  Center_set centers;
  vector <string> file_list;

  string localize_output;
  string orb_out_file;
  int totbasis;

  System * sys;

};

#endif //FITBASIS_METHOD_H_INCLUDED
//------------------------------------------------------------------------
