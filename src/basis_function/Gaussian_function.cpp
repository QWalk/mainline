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

#include "Gaussian_function.h"
#include "qmc_io.h"

/*
 
*/
int Gaussian_function::read(
  vector <string> & words,
  unsigned int & pos)
{
  //cout << "gauss read " << endl;
  unsigned int startpos=pos;

  centername=words[0];

  vector <string> alphatxt;


  if(!readsection(words, pos, alphatxt, "ALPHA")) {
    error("Need ALPHA section in Gaussian_function");
  }
  pos=startpos;

  if(!readvalue(words, pos, rcut, "CUTOFF") ) {
    rcut=1e99;
  }

  if(!readvalue(words, pos, cut_smooth, "SMOOTHING")) {
    cut_smooth=1.2;
  }
  rcut_start=rcut-cut_smooth;


  nmax=alphatxt.size();
  alpha.Resize(nmax);
  for(int i=0; i< nmax; i++) {
    alpha(i)=atof(alphatxt[i].c_str());
    //cout << "alpha(i) " << alpha(i) << " " << alphatxt[i]<<  endl;
    if( alpha(i) <= 0) {
      error("Alpha must be greater than zero in Gaussian function.");
    }
  }

  //cout << "done " << endl;
  return 0;
}

void Gaussian_function::getVarParms(Array1 <doublevar> & parms) {
  //cout << "getVarParms " << endl;
  parms.Resize(nmax);
  for(int i=0; i< nmax; i++) {
    parms(i)=log(alpha(i));
  }
  //cout << "done" << endl;
}

void Gaussian_function::setVarParms(Array1 <doublevar> & parms) {
  //cout << "setVArParms " << endl;
  assert(parms.GetDim(0)==nmax);

  for(int i=0; i< nmax; i++) {
    alpha(i)=exp(parms(i));
  }
  //cout << "done" << endl;
}

int Gaussian_function::nfunc()
{
  return nmax;
}

int Gaussian_function::showinfo(string & indent, ostream & os)
{
  os << indent << "Gaussian function\n";
  os << indent << "Number of functions " << nmax << endl;
  os << indent << "Alpha ";
  for(int i=0; i< nmax; i++) os << alpha(i) << "  ";
  os << endl;
  return 1;
}

int Gaussian_function::writeinput(string & indent, ostream & os)
{
  os << indent << centername << endl;
  os << indent << "GAUSSIAN\n";
  os << indent << "CUTOFF " << rcut << endl;
  os << indent << "SMOOTHING " << cut_smooth << endl;
  os << indent << "ALPHA  { ";
  for(int i=0; i< nmax; i++) {
    os << alpha(i) << "   ";
  }
  os << " } \n";
  return 1;
}

void Gaussian_function::raw_input(ifstream & input)
{error("Raw input not supported by Gaussian_function");}

void Gaussian_function::calcVal(const Array1 <doublevar> & r,
                                Array1 <doublevar> & symvals,
                                const int startfill)
{
  //cout << "calcVal " << endl;
  int index=startfill;

  if(r(0) > rcut) {
    for(int i=0; i< nmax; i++) {
      symvals(index++) = 0;
    }
  }
  else if(r(0) > rcut_start) {
    for(int i=0; i< nmax; i++) {
      doublevar norm=sqrt(alpha(i)/pi);
      doublevar z=(r(0)-rcut_start)/cut_smooth; //cut_smooth=SMOOTHING rcut_start=CUTOFF-SMOOTHING
      symvals(index++) = norm*exp(-alpha(i)*r(1))*(1-z*z*z*(6*z*z-15*z+10));
    }
  }
  else {
    for(int i=0; i< nmax; i++)
    {
      doublevar norm=sqrt(alpha(i)/pi);
      symvals(index++)=norm*exp(-alpha(i)*r(1));
    }
  }
  //cout << "done" << endl;

}

void Gaussian_function::calcLap(
  const Array1 <doublevar> & r,
  Array2 <doublevar> & symvals,
  const int startfill
)
{
  //cout << "calcLap " << endl;
  int index=startfill;
  doublevar exponent, dadr, d2adr2;
  if(r(0) > rcut) {
    for(int i=0; i< nmax; i++) {
      for(int d=0; d< 5; d++) {
        symvals(index,d)=0;
      }
      index++;
    }
  }
  else if(r(0) > rcut_start) {
    doublevar f, dfdr, d2fdr2; //cutoff function
    for(int i=0; i< nmax; i++) {
      doublevar norm=sqrt(alpha(i)/pi);
      exponent=exp(-alpha(i)*r(1));
      dadr=-2*alpha(i)*exponent;
      d2adr2=(4*alpha(i)*alpha(i)*r(1)-6*alpha(i))*exponent;

      //find cutoff function
      doublevar z=(r(0)-rcut_start)/cut_smooth;
      f=(1-z*z*z*(6*z*z-15*z+10));
      dfdr=z*z*30*(2*z-z*z-1)/(cut_smooth*r(0));
      d2fdr2=z*(180*z-120*z*z-60)/(cut_smooth*cut_smooth);

      symvals(index,0)=norm*f*exponent;
      //doublevar dotproduct=0;
      doublevar r_derivative=f*dadr+exponent*dfdr;
      for(int j=1; j<4; j++) {
        symvals(index,j)=norm*(r_derivative)*r(j+1);
        //dotproduct+=dadr*dfdr*r(j+1)*r(j+1);
      }
      symvals(index,4)=norm*(f*d2adr2+2*dadr*dfdr*r(1)+exponent*d2fdr2+2*r_derivative);


      index++;
    }
  }
  else {
    for(int i=0; i< nmax; i++)
    {
      doublevar norm=sqrt(alpha(i)/pi);
      exponent=exp(-alpha(i)*r(1));
      symvals(index,0)=norm*exponent;

      dadr=-2*alpha(i)*exponent;

      d2adr2=(4*alpha(i)*alpha(i)*r(1)-6*alpha(i))*exponent;


      for(int j=1; j<4; j++)
      {
        symvals(index,j)=norm*dadr*r(j+1);
      }
      symvals(index,4)=norm*d2adr2;
      index++;
    }
  }

  //cout << "done" << endl;
}

//------------------------------------------------------------------------
