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

#include "Exponent_cusp.h"
#include "qmc_io.h"

/*
 
*/
int Exponent_cusp::read(
  vector <string> & words,
  unsigned int & pos)
{
  //cout << "gauss read " << endl;
  unsigned int startpos=pos;

  centername=words[0];

  vector <string> alphatxt;


  if(!readvalue(words, pos, cusp, "CUSP")) {
    error("Need CUSP section in Exponent_cusp");
  }
  pos=startpos;

  if(!readvalue(words, pos, gamma, "GAMMA") ) {
    error("Need GAMMA in Exponent_cusp");
  }
  rcut=1e99;


  return 0;
}

void Exponent_cusp::getVarParms(Array1 <doublevar> & parms) {
  //cout << "getVarParms " << endl;
  parms.Resize(1);
  parms(0)=log(gamma);
  //cout << "parms out " << parms(0) << endl;
}

void Exponent_cusp::setVarParms(Array1 <doublevar> & parms) {

  assert(parms.GetDim(0)==1);

  gamma=exp(parms(0));
  //cout << "parms in " << parms(0) << endl;
}

int Exponent_cusp::nfunc()
{
  return 1;
}

int Exponent_cusp::showinfo(string & indent, ostream & os)
{
  os << indent << "Exponential cusp\n";
  os << indent << "Gamma " << gamma;
  os << "   Cusp " << cusp << endl;
  return 1;
}

int Exponent_cusp::writeinput(string & indent, ostream & os)
{
  os << indent << centername << endl;
  os << indent << "EXPONENTIAL_CUSP\n";
  os << indent << "GAMMA " << gamma << endl;
  os << indent << "CUSP " << cusp << endl;
  return 1;
}

void Exponent_cusp::raw_input(ifstream & input)
{error("Raw input not supported by Exponent_cusp");}

void Exponent_cusp::calcVal(const Array1 <doublevar> & r,
                            Array1 <doublevar> & symvals,
                            const int startfill)
{
  //cout << "calcVal " << endl;
  int index=startfill;
  doublevar reducedexp=min(gamma*r(0), (doublevar) 22.0);

  symvals(index)=-cusp*exp(-reducedexp)/gamma;

}

void Exponent_cusp::calcLap(
  const Array1 <doublevar> & r,
  Array2 <doublevar> & symvals,
  const int startfill
)
{
  doublevar reducedexp=min(gamma*r(0), (doublevar) 22.0);
  doublevar exponent=cusp*exp(-reducedexp);
  symvals(startfill,0)=-exponent/gamma;
  for(int d=1; d< 4; d++) {
    symvals(startfill, d)=exponent*r(d+1)/r(0);
  }
  symvals(startfill, 4)=(2/r(0)-gamma)*exponent;
}

//------------------------------------------------------------------------
