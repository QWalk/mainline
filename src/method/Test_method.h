/*
 
Copyright (C) 2007 Lucas K. Wagner
 with modifications by Michal Bajdich

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


#ifndef TEST_METHOD_H_INCLUDED
#define TEST_METHOD_H_INCLUDED

#include "Qmc_std.h"
#include "Qmc_method.h"
#include "Wavefunction.h"
#include "Wavefunction_data.h"
#include "Sample_point.h"
#include "Guiding_function.h"
class Program_options;
#include "Basis_function.h"
#include "Backflow_wf_data.h"
#include "System.h"
/*!
\brief
Do some tests. Used mostly as a development tool.  Keyword: TEST
 
 */
class Test_method : public Qmc_method
{
public:
  void read(vector <string> words,
            unsigned int & pos,
            Program_options & options);
  void run(Program_options & options, ostream & output);
  ~Test_method()
  {

    deallocate(wfdata);
    if(sysprop) delete sysprop;
  }

  int showinfo(ostream & os);
 private:
  void testBackflow();
  void plotCusp(Wavefunction * mywf, Sample_point * sample);
  void testParmDeriv(Wavefunction * mywf, Sample_point * sample);
  int nelectrons; //!< Number of electrons
  string wfoutputfile;
  System * sysprop;
  string readconfig;

  Wavefunction_data * wfdata;
  Pseudopotential * psp;
  Basis_function * basis;

  Backflow_wrapper backflow;
  int test_backflow;
  int plot_cusp;
  int parms_ders;
  int testhessian;
};

#endif //TEST_METHOD_H_INCLUDED
//------------------------------------------------------------------------
