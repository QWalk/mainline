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

#ifndef PROGRAM_OPTIONS_H_INCLUDED
#define PROGRAM_OPTIONS_H_INCLUDED

#include "Qmc_std.h"

/*!
This is just to hold an atom for the input.

*/
class Atom
{
public:

  vector <doublevar> coor;
  //!< Coordinates of the atom
  string name;
  //!< An identifier that's used to match with pseudopotentials, basis sets, etc
  doublevar charge;
  //!< atomic number minus psp removed elecs

  void reset()
  {
    coor[0]=0;
    coor[1]=0;
    coor[2]=0;
    name=" ";
    charge=0;
  }
  Atom()
  {
    coor.resize(3);
    reset();
  }

};


/*!
This class is designed to hold all the options given by the user as they
were given, so that various parts of the program can read the options.
*/
class Program_options
{
public:

  string runid;
  //!< name of the run

  vector < vector < string > > methodtext;
  //!<METHOD text

  vector <vector <string> > pseudotext;
  vector <vector <string> > systemtext;
  //!< SYSTEM text
  vector < vector <string> > twftext;
  //!< TRIALFUNC text.  Since there could be several sections, this is a vector.
  int verbose;
  //!< Possibly: 0-minimal output, 1-medium output, 2-debug output

  Program_options()
  {
    verbose=0;
  }
};


#endif //PROGRAM_OPTIONS_H_INCLUDED
//--------------------------------------------------------------------------
