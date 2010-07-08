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

#include "Spherical_bessel_function.h"
#include "qmc_io.h"

/*
 
*/
int Spherical_bessel_function::read(
  vector <string> & words,
  unsigned int & pos
)
{
  unsigned int startpos=pos;
  centername=words[0];
  doublevar beta0;
  pos=startpos;
  if(!readvalue(words, pos, rcut, "RCUT"))
  {
    error("Need RCUT in SBESSEL");
  }
  pos=startpos;
  if(!readvalue(words, pos=0, nmax, "NFUNC")) {
    error("Need NFUNC in SBESSEL");
  }
  if(nmax > 100)
    error ("Too many functions requested in SBESSEL");

  return 0;
}

void Spherical_bessel_function::getVarParms(Array1 <doublevar> & parms) {
}
void Spherical_bessel_function::setVarParms(Array1 <doublevar> & parms) {
}


int Spherical_bessel_function::nfunc()
{
  return nmax;
}

int Spherical_bessel_function::showinfo(string & indent, ostream & os)
{
  os << indent << "Spherical_bessel function\n";
  os << indent << "cutoff distance : " << rcut << endl;
  os << indent << "highest order function : " << nmax << endl;
  return 1;
}

int Spherical_bessel_function::writeinput(string & indent, ostream & os)
{
  os << indent << centername << endl;
  os << indent << "SBESSEL\n";
  os << indent << "RCUT " << rcut << endl;
  os << indent << "NFUNC " << nmax << endl;
  return 1;
}

void Spherical_bessel_function::raw_input(ifstream & input)
{
  error("raw_input doesn't work for Poly_pade function");
}

void Spherical_bessel_function::calcVal(const Array1 <doublevar> & r,
                                 Array1 <doublevar> & symvals,
                                 const int startfill)
{
  //cout << "Calculating Poly_pade Function\n";
  if( r(0) > rcut)
  {
    int end=startfill+nmax;
    for(int i=startfill; i< end; i++)
    {
      symvals(i)=0;
    }
  }
  else if (r(0)<1e-6){
    int index=startfill;
    doublevar zz=r(0)/rcut;
    doublevar rcutsqrt=sqrt(rcut);
    doublevar rcutt=0.70710678118654752440*rcutsqrt*rcutsqrt*rcutsqrt;
    for(int i=0; i< nmax; i++){
       doublevar zzz=(i+1)*pi*zz;
       doublevar norm=rcutt/((i+1)*pi);
       doublevar zzz2=zzz*zzz;
       symvals(index++)=(1-zzz2/6.0+zzz2*zzz2/120.0)/norm;
    }
  }
  else
  {
    doublevar zz=r(0)/rcut; //reduced radius
    int index=startfill;
    doublevar rcutsqrt=sqrt(rcut);
    doublevar rcutt=0.70710678118654752440*rcutsqrt*rcutsqrt*rcutsqrt;
    for(int i=0; i< nmax; i++)
      {
	doublevar zzz=(i+1)*pi*zz;
	doublevar norm=rcutt/((i+1)*pi);
	symvals(index++)=sin(zzz)/(norm*zzz);
      }
  }
}

void Spherical_bessel_function::calcLap(
  const Array1 <doublevar> & r,
  Array2 <doublevar> & symvals,
  const int startfill
)
{
}

//------------------------------------------------------------------------
