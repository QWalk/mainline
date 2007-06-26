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


#ifndef QMC_METHOD_H_INCLUDED
#define QMC_METHOD_H_INCLUDED

#include "Qmc_std.h"
class Program_options;

class Qmc_method
{
public:

  /*!
    Read the input section from words.
  */
  virtual void read(vector <string> words,
                    unsigned int & pos,
                    Program_options & options)=0;

  /*!
    Generate any large temporary variables needed for a run(some classes provide 
    an interface for using prepared data)
  */
  virtual int generateVariables(Program_options & options) {return 1 ; }
  
  /*!
    Do the calculation, spooling output to output.
  */
  virtual void run(Program_options & options, ostream & output)=0;
  
  /*!
    Print some information about the calculation setup to os.
  */
  virtual int showinfo(ostream & os)
  {
    os << "A Qmc_method  hasn't instantiated showinfo\n";
    return 1;
  }


  virtual ~Qmc_method()
  {}
}
;


class Properties_manager;
class System;
class Wavefunction_data;
class Pseudopotential;

class Qmc_avg_method: public Qmc_method
{
public:

  virtual void read(vector <string> words,
                    unsigned int & pos,
                    Program_options & options)=0;
  virtual int generateVariables(Program_options & options) {return 1 ; }
  
  virtual void run(Program_options & options, ostream & output)=0;

  virtual void runWithVariables(Properties_manager & prop, 
                                System * sys,
                                Wavefunction_data * wfdata,
                                Pseudopotential *,
                                ostream & output)=0;

  virtual int showinfo(ostream & os)
  {
    os << "A Qmc_method  hasn't instantiated showinfo\n";
    return 1;
  }


  virtual ~Qmc_avg_method()
  {}
}
;


int allocate(vector <string> & words,
             Program_options & options,
             Qmc_method * & methptr);
int allocate(vector <string> & words, 
	     Program_options & options,
	     Qmc_avg_method *& methptr);

int deallocate(Qmc_method * & methptr);

#endif //QMC_METHOD_H_INCLUDED
//------------------------------------------------------------------------
