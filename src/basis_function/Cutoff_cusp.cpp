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

#include "Cutoff_cusp.h"
#include "qmc_io.h"

/*
 
*/
int Cutoff_cusp::read(
  vector <string> & words,
  unsigned int & pos)
{
  //cout << "gauss read " << endl;
  unsigned int startpos=pos;

  centername=words[0];

  vector <string> alphatxt;


  if(!readvalue(words, pos, cusp, "CUSP")) {
    error("Need CUSP section in Cutoff_cusp");
  }

  pos=startpos;
  if(!readvalue(words, pos, gamma, "GAMMA") ) {
    error("Need GAMMA in Cutoff_cusp");
  }
  pos=startpos;
  if(!readvalue(words, pos, rcut, "CUTOFF") ) {
    error("Need CUTOFF in Cutoff_cusp");
  }
  rcutinv=1.0/rcut;

  if(rcut < 0) {
    error("Cutoff must be greater than zero in Cutoff_cusp.  I read ", rcut);
  }
  pade0=1./(3.+gamma);

  return 0;
}

void Cutoff_cusp::getVarParms(Array1 <doublevar> & parms) {
  //cout << "getVarParms " << endl;
  parms.Resize(1);
  parms(0)=log(gamma);
  //cout << "parms out " << parms(0) << endl;
}

void Cutoff_cusp::setVarParms(Array1 <doublevar> & parms) {

  assert(parms.GetDim(0)==1);

  gamma=exp(parms(0));
  pade0=1./(3+gamma);
  //cout << "parms in " << parms(0) << endl;
}

int Cutoff_cusp::nfunc()
{
  return 1;
}

int Cutoff_cusp::showinfo(string & indent, ostream & os)
{
  os << indent << "Cutoff cusp\n";
  os << indent << "Gamma " << gamma;
  os << "   Cusp " << cusp << endl;
  return 1;
}

int Cutoff_cusp::writeinput(string & indent, ostream & os)
{
  os << indent << centername << endl;
  os << indent << "CUTOFF_CUSP\n";
  os << indent << "GAMMA " << gamma << endl;
  os << indent << "CUSP " << cusp << endl;
  os << indent << "CUTOFF " << rcut << endl;
  return 1;
}

void Cutoff_cusp::raw_input(ifstream & input)
{error("Raw input not supported by Cutoff_cusp");}

void Cutoff_cusp::calcVal(const Array1 <doublevar> & r,
                          Array1 <doublevar> & symvals,
                          const int startfill)
{
  assert(r.GetDim(0) >= 5);
  assert(symvals.GetDim(0) >= 1+startfill);
  if(r(0) < rcut) {
    doublevar zz=r(0)/rcut;
    doublevar zz2=zz*zz;
    doublevar pp=zz-zz2+zz*zz2/3;
    doublevar pade=1./(1+gamma*pp);
    symvals(startfill)=cusp*rcut*(pp*pade-pade0);
  }
  else {
    symvals(startfill)=0;
  }

}

void Cutoff_cusp::calcLap(
  const Array1 <doublevar> & r,
  Array2 <doublevar> & symvals,
  const int startfill
)
{
  assert(r.GetDim(0) >=5);
  assert(symvals.GetDim(0) >= 1+startfill);
  assert(symvals.GetDim(1) >= 5);
  if(r(0) < rcut) {
    doublevar zz=r(0)*rcutinv;
    doublevar zz2=zz*zz;
    doublevar pp=zz-zz2+zz*zz2/3;
    doublevar pade=1./(1+gamma*pp);
    symvals(startfill,0)=cusp*rcut*(pp*pade-pade0);
    doublevar pade2=pade*pade;
    doublevar ppd=1.-2.*zz+zz2;
    doublevar ppdd=-2.+2.*zz;
    doublevar dadr=ppd*pade2/r(0);
    doublevar dadr2=ppdd*pade2*rcutinv
                    -2.*gamma*ppd*ppd*pade2*pade*rcutinv
                    +2.*dadr;

    for(int i=1; i< 4; i++) {
      symvals(startfill, i)=cusp*dadr*r(i+1);
    }
    symvals(startfill, 4)=cusp*dadr2;
  }
  else {
    for(int i=0; i< 5; i++)
      symvals(startfill, i)=0;
  }
}

//------------------------------------------------------------------------
