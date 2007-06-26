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

#include "Step_function.h"
#include "qmc_io.h"

/*
 
*/
int Step_function::read(
  vector <string> & words,
  unsigned int & pos)
{


  centername=words[0];


  if(!readvalue(words, pos=0, rcut, "CUTOFF") ) {
    error("Need CUTOFF in step_function");
  }

  return 1;
}

void Step_function::getVarParms(Array1 <doublevar> & parms) {
  //cout << "getVarParms " << endl;
  parms.Resize(0);
}

void Step_function::setVarParms(Array1 <doublevar> & parms) {
  assert(parms.GetDim(0)==0);
}

int Step_function::nfunc()
{
  return 1;
}

int Step_function::showinfo(string & indent, ostream & os)
{
  os << indent << "Step function\n";
  os << indent << "Cutoff " << rcut << endl;

  return 1;
}

int Step_function::writeinput(string & indent, ostream & os)
{
  os << indent << centername << endl;
  os << indent << "STEP\n";
  os << indent << "CUTOFF " << rcut << endl;
  return 1;
}

void Step_function::raw_input(ifstream & input)
{error("Raw input not supported by Step_function");}

void Step_function::calcVal(const Array1 <doublevar> & r,
                                Array1 <doublevar> & symvals,
                                const int startfill)
{
  //cout << "calcVal " << endl;
  if(r(0) < rcut)
    symvals(startfill)=1;
  else
    symvals(startfill)=0;
}

void Step_function::calcLap(
  const Array1 <doublevar> & r,
  Array2 <doublevar> & symvals,
  const int startfill
)
{
  if(r(0) < rcut)
    symvals(startfill,0)=1;
  else
    symvals(startfill,0)=0;

  for(int d=1; d< 5; d++)
    symvals(startfill, d)=0;

}

//------------------------------------------------------------------------
