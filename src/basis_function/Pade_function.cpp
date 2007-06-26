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

#include "Pade_function.h"
#include "qmc_io.h"

/*

*/
int Pade_function::read(
  vector <string> & words,
  unsigned int & pos
)
{

  unsigned int startpos=pos;
  centername=words[0];
  rcut=10;   //shouldn't be used at the moment.
  
  if(!readvalue(words, pos=startpos, nmax, "NFUNC")) {
    error("Need NMAX in PADE function");
  }
  alpha.Resize(nmax);


  if(readvalue(words, pos=startpos, alpha0, "ALPHA0")) {
    for(int i=0; i< nmax; i++)
    {
      alpha(i)=alpha0/pow(2.0, i);
    }
  }
  else {
    error("Need ALPHA0 in PADE function");
  }

  return 0;
}

void Pade_function::getVarParms(Array1 <doublevar> & parms) {
  parms.Resize(1);
  parms(0)=log(alpha0);
}

void Pade_function::setVarParms(Array1 <doublevar> & parms) {
  assert(parms.GetDim(0)==1);

  alpha0=exp(parms(0));
  for(int i=0; i< nmax; i++) {
    alpha(i)=alpha0/pow(2.0,i);
  }

}

int Pade_function::nfunc()
{
  return nmax;
}

int Pade_function::showinfo(string & indent, ostream & os)
{
  os << indent << "Pade function\n";
  os << indent << "Number of functions " << nmax << endl;
  os << indent << "Alpha_0  " << alpha(0) << endl;
  return 1;
}

int Pade_function::writeinput(string & indent, ostream & os)
{
  os << indent << centername << endl;
  os << indent << "PADE\n";
  os << indent << "NFUNC " << nmax << endl;
  os << indent << "ALPHA0 " << alpha(0) << endl;
  return 1;
}

void Pade_function::raw_input(ifstream & input)
{}

void Pade_function::calcVal(const Array1 <doublevar> & r,
                            Array1 <doublevar> & symvals,
                            const int startfill)
{
  //Removing the cutoff radius for the moment
  //cout << "Calculating Pade Function\n";
  //if( r(0) > rcut) {
  //symvals=0;
  //}
  //else {
  //if(startfill != 1) cout << "wrong startfill" << endl;
  doublevar pade;   //


  int index=startfill;
  for(int i=0; i< nmax; i++)
  {
    //This can probably be optimized with a little effort.
    //cout << "i=" << i << endl;
    //cout << "alpha " << alpha(i) << endl;
    pade=alpha(i)/(1.+alpha(i)*r(0));
    symvals(index++)=pade*pade*r(1);
  }

  //cout << "Leaving Pade function\n";

  //for(int i=0; i< symvals.GetDim(1); i++) {
  //  cout << "val:symvals " << i << symvals(i) << endl;
  //}

}

void Pade_function::calcLap(
  const Array1 <doublevar> & r,
  Array2 <doublevar> & symvals,
  const int startfill
)
{

  //if( r(0) > rcut) {
  // symvals=0;
  // }
  //else {
  //if(startfill != 1) cout << "wrong startfill" << endl;
  doublevar pade,   //
  app,  //alpha(i)*pade*pade
  dadr,       //the first derivative wrt r
  d2adr2;     //second derivative wrt r
  //symvals.Resize(nmax+1, 5);
  //symvals(0,0)=1;
  //for(int d=1; d< 5; d++) {
  //    symvals(0,d)=0;
  //}
  int index=startfill;
  for(int i=0; i< nmax; i++)
  {
    //This can probably be optimized with a little effort.
    //cout << "i=" << i << endl;
    //cout << "alpha " << alpha(i) << endl;
    pade=1.0/(1.+alpha(i)*r(0));
    app=alpha(i)*alpha(i)*pade*pade;
    symvals(index,0)=app*r(1);

    dadr=app*pade;
    dadr+=dadr;

    d2adr2=3.0*dadr*pade;

    //cout << "a/b_m  " << symvals(i+1,0) << endl;
    //cout << "a/b_md  " << dadr << endl;
    //cout << "d2adr2  " << d2adr2 << endl;
    for(int j=1; j<4; j++)
    {
      symvals(index,j)=dadr*r(j+1);
    }
    symvals(index,4)=d2adr2;
    index++;
  }
  //for(int i=0; i< symvals.GetDim(1); i++) {
  //  for(int d=0; d< 5; d++) {
  //    cout << "lap:symvals(" << i<<"," << d << ") "
  //	 << symvals(i,d) << endl;
  //  }
  //}

  //}


}

//------------------------------------------------------------------------
