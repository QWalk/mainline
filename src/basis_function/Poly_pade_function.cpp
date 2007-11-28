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

#include "Poly_pade_function.h"
#include "qmc_io.h"

/*
 
*/
int Poly_pade_function::read(
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
// rcut=10;
     error("Need RCUT in POLYPADE");
  }
  pos=startpos;

  vector <string> betatxt;
  if(readsection(words, pos, betatxt, "BETA"))
  {

    cout << "BETA has been replaced with BETA0/NFUNC." << endl;
    if(betatxt.size()==1) {
      cout << "You can replace BETA { " << betatxt[0] << " } ";
      cout << " with  BETA0 " << betatxt[0] << "  NFUNC  1 " << endl;

    }
    else {
      cout << "there are several betas, which actually doesn't make sense."
           << endl;
    }

    error("done");
    //beta.Resize(betatxt.size());
    //nmax=betatxt.size();
    //for(int i=0; i< nmax; i++)
    //{
    //  beta(i)=atof(betatxt[i].c_str());
    //  if(beta(i) <= -1)
    //  {
    //    error("In POLYPADE, beta must be greater than -1, but I read ",
    //          beta(i));
    //  }
    //}
  }

  if(readvalue(words, pos=0, beta0, "BETA0")) {
    if(!readvalue(words, pos=0, nmax, "NFUNC")) {
      error("Need NFUNC with BETA0 in POLY_PADE");
    }
    beta.Resize(nmax);
    beta(0)=beta0;
    double beta1=log(beta0+1.00001);
    for(int i=1; i< nmax; i++) {
      beta(i)=exp(beta1+1.6*i)-1;
    }
  }
  else {
    error("Couldn't find BETA in POLY_PADE function " );
  }

  return 0;
}

void Poly_pade_function::getVarParms(Array1 <doublevar> & parms) {
  parms.Resize(1);
  parms(0)=log(beta(0)+.9999999);

/*
  int parmmax=nparms();
  parms.Resize(parmmax);
  for(int i=0; i< parmmax; i++) {
    parms(i)=log(beta(i)+1.0 );
  }
*/
}
void Poly_pade_function::setVarParms(Array1 <doublevar> & parms) {
  assert(parms.GetDim(0)==1);
  double beta0=exp(parms(0))-.9999999;
  beta(0)=beta0;
  
  double beta1=log(beta0+1.00001);
  for(int i=1; i< nmax; i++) {
    beta(i)=exp(beta1+1.6*i)-1;
    //cout << "betas " << beta(i) << endl;
  }
}


int Poly_pade_function::nfunc()
{
  return nmax;
}

int Poly_pade_function::showinfo(string & indent, ostream & os)
{
  os << indent << "Poly_pade function\n";
  os << indent << "beta : ";
  for(int i=0; i< nmax; i++) 
    os << beta(i) << "  ";
  os << endl;
  os << indent << "cutoff distance : " << rcut << endl;
  return 1;
}

int Poly_pade_function::writeinput(string & indent, ostream & os)
{
  os << indent << centername << endl;
  os << indent << "POLYPADE\n";
  os << indent << "RCUT " << rcut << endl;
  //os << indent << "BETA { ";
  //for(int i=0; i< nmax; i++)
  //{
  //  os << beta(i) << "   ";
  //}
  //os << " }\n";
  os << indent << "BETA0 " << beta(0) << endl;
  os << indent << "NFUNC " << nmax << endl;
  return 1;
}

void Poly_pade_function::raw_input(ifstream & input)
{
  error("raw_input doesn't work for Poly_pade function");
}

void Poly_pade_function::calcVal(const Array1 <doublevar> & r,
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
  else
  {
    doublevar zz=r(0)/rcut; //reduced radius
    doublevar zz2=zz*zz; //zz squared
    doublevar zpp=zz2*(6.-8.*zz+3.*zz2);
    doublevar zpp1=1-zpp;
    doublevar zpade,crs;
    int index=startfill;
    for(int i=0; i< nmax; i++)
    {
      zpade=1./(1+beta(i)*zpp);
      crs=zpp1*zpade;
      symvals(index++)=crs;
    }
  }
}

void Poly_pade_function::calcLap(
  const Array1 <doublevar> & r,
  Array2 <doublevar> & symvals,
  const int startfill
)
{

  //cout << "Calculating Poly_pade Function\n";
  if( r(0) >= rcut)
  {
    //cout << "cutting off " << r(0) << endl;
    int end=startfill+nmax;
    for(int i=startfill; i< end; i++)
    {
      for(int d=0; d<5; d++)
      {
        symvals(i,d)=0;
      }
    }
  }
  else
  {
    doublevar zz=r(0)/rcut; //reduced radius
    doublevar zz2=zz*zz; //zz squared
    doublevar zpp=zz2*(6.-8.*zz+3.*zz2);
    doublevar zpp1=1-zpp;
    doublevar zppd=12.*zz*(1.-2.*zz+zz2);
    doublevar zppdd2=6.-24.*zz+18*zz2; //dzpp/dzz^2 divided by two
    doublevar zpade, zpade2, crs, crsd, crsdd2;
    int index=startfill;
    //cout << "-----------\n";
    //cout << this << " r(0)  " << r(0) << endl;
    for(int i=0; i< nmax; i++)
    {
      //cout << "beta " << beta(i) << endl;
      zpade=1./(1+beta(i)*zpp);
      zpade2=zpade*zpade;

      crs=zpp1*zpade;
      crsd=-zppd*zpade2*(beta(i)+1)/(r(0)*rcut);
      crsdd2=-zpade2*(zppdd2-beta(i)*zppd*zppd*zpade)
             *(beta(i)+1)/(rcut*rcut); //d^2a/dr^2 divided by 2
      //cout << "crs " << crs << " crsd " << crsd
      //   << "  crsdd2 " << crsdd2 << endl;
      symvals(index,0)=crs;
      for(int d=1; d< 4; d++)
      {
        symvals(index,d)=crsd*r(d+1);
      }
      symvals(index,4)=2*(crsdd2+crsd);
      //cout << " laplacian times onehalf " << 2*(crsdd2+crsd) << endl;
      index++;
    }
    }
}

//------------------------------------------------------------------------
