/*
 
Copyright (C) 2007 Lucas K. Wagner, 2009 J. Kolorenc

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

#include "Cosine_function.h"
#include "qmc_io.h"

int Cosine_function::read(
  vector <string> & words,
  unsigned int & pos)
{
  unsigned int startpos=pos;

  centername=words[0];

  pos=startpos;
  vector <string> gvectxt;
  if(!readsection(words, pos, gvectxt, "GVECTOR")) {
    error("Need GVECTOR in COSINE input");
  }
  if(gvectxt.size() % 3 != 0) {
    error("Bad count on GVECTOR."
          "  There must be sets of three numbers.");
  }
  nmax=gvectxt.size()/3;
  //cout << nmax << " g-vectors " << endl;
  g_vector.Resize(nmax, 3);
  g_vec_sqrd.Resize(nmax);
  int counter=0;
  for(int i=0; i< nmax; i++) {
    g_vec_sqrd(i)=0.0;
    for(int d=0; d< 3; d++) {
      g_vector(i,d)=atof(gvectxt[counter].c_str());
      g_vec_sqrd(i)+=g_vector(i,d)*g_vector(i,d);
      counter++;
    }
  }
  return 0;
}

void Cosine_function::getVarParms(Array1 <doublevar> & parms) {
  //cout << "getVarParms " << endl;
  parms.Resize(0);
}

void Cosine_function::setVarParms(Array1 <doublevar> & parms) {
  //cout << "setVArParms " << endl;
  assert(parms.GetDim(0)==0);
}

int Cosine_function::nfunc()
{
  //cout << "Cosine_function::nfunc() " << nmax << endl;
  return nmax;
}

int Cosine_function::showinfo(string & indent, ostream & os)
{
  os << indent << "Cosine function" << endl;
  os << indent << nmax << " functions" << endl;
  return 1;
}

int Cosine_function::writeinput(string & indent, ostream & os)
{
  os << indent << centername << endl;
  os << indent << "COSINE" << endl;
  os << indent << "GVECTOR {" << endl;
  for(int i=0; i< nmax; i++) {
    os << indent << "  " << g_vector(i,0) << "   " << g_vector(i,1)
        << "  " << g_vector(i,2) << endl;
  }
  os << indent << "}" << endl;
  return 1;
}

void Cosine_function::raw_input(ifstream & input)
{
  error("Raw input not supported by Cosine_function");
}

void Cosine_function::calcVal(const Array1 <doublevar> & r,
			      Array1 <doublevar> & symvals,
			      const int startfill)
{
  //cout << "calcVal " << endl;
  assert(r.GetDim(0) >= 5);
  assert(symvals.GetDim(0) >= nmax+startfill);
  int index=startfill;
  doublevar gdotr, t_cos;
  for(int i=0; i< nmax; i++) {
    gdotr=g_vector(i,0)*r(2)+g_vector(i,1)*r(3)+g_vector(i,2)*r(4);
    t_cos=cos(gdotr);
    symvals(index++)=t_cos;
  }
  //cout << "done" << endl;
}

void Cosine_function::calcLap(
  const Array1 <doublevar> & r,
  Array2 <doublevar> & symvals,
  const int startfill
)
{
  //cout << "calcLap " << endl;
  assert(r.GetDim(0) >= 5);
  assert(symvals.GetDim(0) >= nmax+startfill);
  assert(symvals.GetDim(1) >= 5);

  int index=startfill;
  doublevar gdotr, t_cos, t_sin;
  doublevar gsquared;
  for(int fn=0; fn< nmax; fn++) {
    gdotr=g_vector(fn,0)*r(2)+g_vector(fn,1)*r(3)+g_vector(fn,2)*r(4);

#ifdef __USE_GNU
    sincos(gdotr, &t_sin, &t_cos);
#else
    t_cos=cos(gdotr);
    t_sin=sin(gdotr);
#endif
    symvals(index, 0)=t_cos;
    for(int i=1; i< 4; i++) {
      symvals(index, i)=-g_vector(fn, i-1)*t_sin;
    }
    symvals(index, 4)=-g_vec_sqrd(fn)*t_cos;
    index++;
  }
  //cout << "done" << endl;
}



void Cosine_function::calcHessian(
  const Array1 <doublevar> & r,
  Array2 <doublevar> & symvals,
  const int startfill
)
{
  //cout << "calcLap " << endl;
  assert(r.GetDim(0) >= 5);
  assert(symvals.GetDim(0) >= nmax+startfill);
  assert(symvals.GetDim(1) >= 10);

  int index=startfill;
  doublevar gdotr, t_cos, t_sin;
  doublevar gx, gy, gz;
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
  }
  //cout << "done" << endl;
}







//------------------------------------------------------------------------
