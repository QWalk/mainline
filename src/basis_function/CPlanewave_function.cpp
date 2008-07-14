/*
 
Copyright (C) 2008 Jindrich Kolorenc

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

#include "CPlanewave_function.h"
#include "qmc_io.h"

/*!
   Note:  We may need to change this to read directly from a file, for memory
   constraints..
*/
int CPlanewave_function::read(
  vector <string> & words,
  unsigned int & pos)
{
  //cout << "making planewave basis " << endl;
  unsigned int startpos=pos;

  centername=words[0];

  pos=startpos;
  vector <string> gvectxt;
  if(!readsection(words, pos, gvectxt, "GVECTOR")) {
    error("Need GVECTOR in CPLANEWAVE input");
  }
  if(gvectxt.size() % 3 != 0) {
    error("Bad count on GVECTOR."
          "  There must be sets of three numbers.");
  }
  nmax=gvectxt.size()/3;
  //cout << nmax << " g-vectors " << endl;
  g_vector.Resize(nmax, 3);
  int counter=0;
  for(int i=0; i< nmax; i++) {
    for(int d=0; d< 3; d++) {
      g_vector(i,d)=atof(gvectxt[counter].c_str());
      counter++;
    }
  }
  //cout << "done " << endl;
  return 0;
}

void CPlanewave_function::getVarParms(Array1 <doublevar> & parms) {
  //cout << "getVarParms " << endl;
  parms.Resize(0);
}

void CPlanewave_function::setVarParms(Array1 <doublevar> & parms) {
  //cout << "setVArParms " << endl;
  assert(parms.GetDim(0)==0);
}

int CPlanewave_function::nfunc()
{
  //cout << "Planewave::nfunc() " << nmax << endl;
  return nmax;
}

int CPlanewave_function::showinfo(string & indent, ostream & os)
{
  os << indent << "Planewave function\n";
  os << indent << nmax << " plane waves" << endl;
  return 1;
}

int CPlanewave_function::writeinput(string & indent, ostream & os)
{
  os << indent << centername << endl;
  os << indent << "CPLANEWAVE\n";
  os << indent << "GVECTOR { \n";
  for(int i=0; i< nmax; i++) {
    os << indent << "  " << g_vector(i,0) << "   " << g_vector(i,1)
        << "  " << g_vector(i,2) << endl;
  }
  os << indent << "}" << endl;
  return 1;
}


void CPlanewave_function::calcVal(const Array1 <doublevar> & r,
                                 Array1 <dcomplex> & symvals,
                                 const int startfill)
{
  const dcomplex I(0.0,1.0);
  //cout << "calcVal " << endl;
  assert(r.GetDim(0) >= 5);
  assert(symvals.GetDim(0) >= nmax+startfill);
  int index=startfill;
  doublevar gdotr;
  for(int i=0; i< nmax; i++) {
    gdotr=g_vector(i,0)*r(2)+g_vector(i,1)*r(3)+g_vector(i,2)*r(4);
    symvals(index++)=exp(I*gdotr);
  }
  //cout << "done" << endl;
}

void CPlanewave_function::calcLap(
  const Array1 <doublevar> & r,
  Array2 <dcomplex> & symvals,
  const int startfill
)
{
  const dcomplex I(0.0,1.0);
  //cout << "calcLap " << endl;
  assert(r.GetDim(0) >= 5);
  assert(symvals.GetDim(0) >= nmax+startfill);
  assert(symvals.GetDim(1) >= 5);

  int index=startfill;
  doublevar gdotr;
  doublevar gsquared;
  dcomplex t_exp;
  for(int fn=0; fn< nmax; fn++) {
    gdotr=g_vector(fn,0)*r(2)+g_vector(fn,1)*r(3)+g_vector(fn,2)*r(4);
    t_exp=exp(I*gdotr);

    //Should probably store this one..
    gsquared=g_vector(fn,0)*g_vector(fn,0)
             +g_vector(fn,1)*g_vector(fn,1)
             +g_vector(fn,2)*g_vector(fn,2);
    
    symvals(index,0)=t_exp;
    for(int i=1; i< 4; i++) {
      symvals(index, i)=I*g_vector(fn, i-1)*t_exp;
    }
    symvals(index, 4)=-gsquared*t_exp;
    index++;
  }
  //cout << "done" << endl;
}



//------------------------------------------------------------------------
