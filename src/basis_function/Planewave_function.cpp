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

#include "Planewave_function.h"
#include "qmc_io.h"

/*!
   Note:  We may need to change this to read directly from a file, for memory
   constraints..
*/
int Planewave_function::read(
  vector <string> & words,
  unsigned int & pos)
{
  //cout << "making planewave basis " << endl;
  unsigned int startpos=pos;

  centername=words[0];

  pos=startpos;
  vector <string> gvectxt;
  if(!readsection(words, pos, gvectxt, "GVECTOR")) {
    error("Need GVECTOR in PLANEWAVE input");
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

void Planewave_function::getVarParms(Array1 <doublevar> & parms) {
  //cout << "getVarParms " << endl;
  parms.Resize(0);
}

void Planewave_function::setVarParms(Array1 <doublevar> & parms) {
  //cout << "setVArParms " << endl;
  assert(parms.GetDim(0)==0);
}

int Planewave_function::nfunc()
{
  //cout << "Planewave::nfunc() " << nmax << endl;
  return 2*nmax;
}

int Planewave_function::showinfo(string & indent, ostream & os)
{
  os << indent << "Planewave function\n";
  os << indent << nmax << " plane waves" << endl;
  return 1;
}

int Planewave_function::writeinput(string & indent, ostream & os)
{
  os << indent << centername << endl;
  os << indent << "PLANEWAVE\n";
  os << indent << "GVECTOR { \n";
  for(int i=0; i< nmax; i++) {
    os << indent << "  " << g_vector(i,0) << "   " << g_vector(i,1)
        << "  " << g_vector(i,2) << endl;
  }
  os << indent << "}" << endl;
  return 1;
}

void Planewave_function::raw_input(ifstream & input)
{error("Raw input not supported by Planewave_function");}

void Planewave_function::calcVal(const Array1 <doublevar> & r,
                                 Array1 <doublevar> & symvals,
                                 const int startfill)
{
  //cout << "calcVal " << endl;
  assert(r.GetDim(0) >= 5);
  assert(symvals.GetDim(0) >= nmax*2+startfill);
  int index=startfill;
  doublevar gdotr, t_cos, t_sin;
  for(int i=0; i< nmax; i++) {
    gdotr=g_vector(i,0)*r(2)+g_vector(i,1)*r(3)+g_vector(i,2)*r(4);
#ifdef __USE_GNU
    sincos(gdotr, &t_sin, &t_cos);
#else
    t_cos=cos(gdotr);
    t_sin=sin(gdotr);
#endif
    symvals(index++)=t_cos;
    symvals(index++)=t_sin;
  }
  //cout << "done" << endl;
}

void Planewave_function::calcLap(
  const Array1 <doublevar> & r,
  Array2 <doublevar> & symvals,
  const int startfill
)
{
  //cout << "calcLap " << endl;
  assert(r.GetDim(0) >= 5);
  assert(symvals.GetDim(0) >= nmax*2+startfill);
  assert(symvals.GetDim(1) >= 5);

  int index=startfill;
  doublevar gdotr, t_cos, t_sin;
  doublevar gsquared;
  for(int fn=0; fn< nmax; fn++) {
    gdotr=g_vector(fn,0)*r(2)+g_vector(fn,1)*r(3)+g_vector(fn,2)*r(4);

    //Should probably store this one..
    gsquared=g_vector(fn,0)*g_vector(fn,0)
             +g_vector(fn,1)*g_vector(fn,1)
             +g_vector(fn,2)*g_vector(fn,2);
#ifdef __USE_GNU
    sincos(gdotr, &t_sin, &t_cos);
#else
    t_cos=cos(gdotr);
    t_sin=sin(gdotr);
#endif
    //cos function
    symvals(index, 0)=t_cos;
    for(int i=1; i< 4; i++) {
      symvals(index, i)=-g_vector(fn, i-1)*t_sin;
    }
    symvals(index, 4)=-gsquared*t_cos;
    index++;

    //sin function
    symvals(index, 0)=t_sin;
    for(int i=1; i< 4; i++) {
      symvals(index, i)=g_vector(fn,i-1)*t_cos;
    }
    symvals(index, 4)=-gsquared*t_sin;
    index++;
  }
  //cout << "done" << endl;
}



void Planewave_function::calcHessian(
  const Array1 <doublevar> & r,
  Array2 <doublevar> & symvals,
  const int startfill
)
{
  //cout << "calcLap " << endl;
  assert(r.GetDim(0) >= 5);
  assert(symvals.GetDim(0) >= nmax*2+startfill);
  assert(symvals.GetDim(1) >= 10);

  int index=startfill;
  doublevar gdotr, t_cos, t_sin;
  doublevar gx, gy, gz;
  //  doublevar gsquared;
  for(int fn=0; fn< nmax; fn++) {
    gx=g_vector(fn, 0);
    gy=g_vector(fn, 1);
    gz=g_vector(fn, 2);
    gdotr=gx*r(2)+gy*r(3)+gz*r(4);

#ifdef __USE_GNU
    sincos(gdotr, &t_sin, &t_cos);
#else
    t_cos=cos(gdotr);
    t_sin=sin(gdotr);
#endif
    
    //cos function
    symvals(index, 0)=t_cos;
    for(int i=1; i< 4; i++) {
      symvals(index, i)=-g_vector(fn, i-1)*t_sin;
    }

    symvals(index, 4)=-gx*gx*t_cos;
    symvals(index, 5)=-gy*gy*t_cos;
    symvals(index, 6)=-gz*gz*t_cos;
    symvals(index, 7)=-gx*gy*t_cos;
    symvals(index, 8)=-gx*gz*t_cos;
    symvals(index, 9)=-gy*gz*t_cos;
    
    index++;

    //sin function
    symvals(index, 0)=t_sin;
    for(int i=1; i< 4; i++) {
      symvals(index, i)=g_vector(fn,i-1)*t_cos;
    }
    symvals(index, 4)=-gx*gx*t_sin;
    symvals(index, 5)=-gy*gy*t_sin;
    symvals(index, 6)=-gz*gz*t_sin;
    symvals(index, 7)=-gx*gy*t_sin;
    symvals(index, 8)=-gx*gz*t_sin;
    symvals(index, 9)=-gy*gz*t_sin;
    
    index++;
  }
  //cout << "done" << endl;
}







//------------------------------------------------------------------------
