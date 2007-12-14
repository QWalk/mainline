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

#include "Rgaussian_function.h"
#include "qmc_io.h"

int Rgaussian_function::read(
  vector <string> & words,
  unsigned int & pos)
{
  //cout << "gauss read " << endl;
  unsigned int startpos=pos;

  centername=words[0];


  rcut=1e99;
  vector < vector <doublevar> >  exptemp;
  vector < vector <doublevar> >  coefftemp;
  vector < vector <doublevar> >  radtemp;
  vector <string> oldqmctxt;

  int maxL=0, maxExpansion=0;
  if(readsection(words, pos=0, oldqmctxt, "OLDQMC")) {
    double zeff;
    int nlval; //number of l values that we'll read in
    if(oldqmctxt.size() < 2)
      error("badly formed OLDQMC section: size less than 2");
    int it=0; //position iterator
    zeff=atof(oldqmctxt[it++].c_str());
    nlval=atoi(oldqmctxt[it++].c_str());

    if(maxL < nlval) maxL=nlval;
    Array1 <int> nPerL(nlval); //number per l-value

    for(int l=0; l < nlval; l++) {
      nPerL(l)=atoi(oldqmctxt[it++].c_str());
      if(maxExpansion < nPerL(l))
        maxExpansion=nPerL(l);
    }
    for(int l=0; l < nlval; l++) {
      vector <doublevar> radtempl;
      vector <doublevar> exptempl;
      vector <doublevar> coefftempl;

      for(int j=0; j< nPerL(l); j++) {
        doublevar temp;
        temp=atof(oldqmctxt[it++].c_str());
        temp-=2;
        radtempl.push_back(temp);
        temp=atof(oldqmctxt[it++].c_str());
        exptempl.push_back(temp);

        temp=atof(oldqmctxt[it++].c_str());
        coefftempl.push_back(temp);
      }
      radtemp.push_back(radtempl);
      exptemp.push_back(exptempl);
      coefftemp.push_back(coefftempl);
    }
  }
  else {
    error("Try OLDQMC in RGAUSSIAN.");
  }

  gaussexp.Resize(maxL,  maxExpansion);
  gausscoeff.Resize(maxL, maxExpansion);
  radiusexp.Resize(maxL, maxExpansion);
  numExpansion.Resize(maxL);
  nmax=radtemp.size();

  //cout << "maxL " << maxL << "  maxExpansion " << maxExpansion << endl;

  assert( radtemp.size() == exptemp.size() );


  for(int l=0; l< nmax; l++)
  {

    numExpansion(l) = radtemp[l].size();

    for(int e=0; e< numExpansion(l); e++)
    {
      gaussexp(l, e) = exptemp[l][e];
      gausscoeff(l, e) = coefftemp[l][e];
      radiusexp(l, e) = radtemp[l][e];
    }
  }
  return 1;
}

void Rgaussian_function::getVarParms(Array1 <doublevar> & parms) {
  //cout << "getVarParms " << endl;
  parms.Resize(nparms());

  int count=0;
  for(int l=0; l < nmax; l++) { 
    for(int e=0; e< numExpansion(l); e++) { 
      parms(count++)=log(gaussexp(l,e));
      //parms(count++)=gausscoeff(l,e);
    }
  }
}

void Rgaussian_function::setVarParms(Array1 <doublevar> & parms) {
  //cout << "setVArParms " << endl;
  assert(parms.GetDim(0)==nparms());
  int count=0;
  for(int l=0; l < nmax; l++) { 
    for(int e=0; e< numExpansion(l); e++) { 
      gaussexp(l,e)=exp(parms(count++));
      //gausscoeff(l,e)=parms(count++);
    }
  }
}

int Rgaussian_function::nfunc()
{
  return nmax;
}

int Rgaussian_function::showinfo(string & indent, ostream & os)
{
  os << indent <<  "R^n *Gaussian" << endl;
  
  for(int i=0; i< nmax; i++) {
    os << indent << "Expansion " << i << endl;
    os << indent << setw(10) << "coeff" 
       << setw(10) << "n" << setw(10) << "exponent" << endl;
    for(int j=0; j< numExpansion(i); j++) {
      os << indent << setw(10) << gausscoeff(i,j)  
         << setw(10) << radiusexp(i,j) 
         << setw(10) << gaussexp(i,j) << endl;
    }
  }
  return 1;
}

int Rgaussian_function::writeinput(string & indent, ostream & os)
{
  os << indent << centername << endl;
  os << indent << "RGAUSSIAN\n";
  os << indent << "OLDQMC { \n";
  os << indent << "  0.  " << nmax << endl;
  os << indent << "  ";
  for(int i=0; i< nmax; i++) os << numExpansion(i) << "  ";
  os << endl;
  for (int i=0; i< nmax; i++) { 
    for(int j=0; j< numExpansion(i); j++) { 
      os << indent << "   " 
	 << "  " << radiusexp(i,j)+2
	 << "   " << gaussexp(i,j) 
	 << "    " << gausscoeff(i,j) << endl;
    }
  }
  os << indent << "}\n";
  return 1;
}

void Rgaussian_function::raw_input(ifstream & input)
{error("Raw input not supported by Rgaussian_function");}

void Rgaussian_function::calcVal(const Array1 <doublevar> & r,
                                 Array1 <doublevar> & symvals,
                                 const int startfill)
{
  //cout << "Rgaussian::calcVal " << endl;
  assert(symvals.GetDim(0) >= nmax+startfill);
  assert(r.GetDim(0) >= 5);
  doublevar exponent;
  //cout << "v_l for l=" << l << "   atom=" << at << endl;
  int index=startfill;

  for(int i=0; i< nmax; i++) {
    doublevar v_l=0;
    for(int j=0; j< numExpansion(i); j++) {
      exponent=gaussexp(i, j)*r(1);
      exponent=min(exponent,(doublevar) 60.0);
      v_l+= gausscoeff(i, j)
            *pow(r(0), radiusexp(i, j) )
            *exp(-exponent);
    }
    symvals(index++) = v_l;
  }
}

void Rgaussian_function::calcLap(
  const Array1 <doublevar> & r,
  Array2 <doublevar> & symvals,
  const int startfill){

  assert(symvals.GetDim(0) >= nmax+startfill);
  assert(r.GetDim(0) >= 5);
  doublevar exponent;
  //cout << "v_l for l=" << l << "   atom=" << at << endl;
  int index=startfill;

  for(int i=0; i< nmax; i++) {
    doublevar v_l=0, dv_l=0, d2v_l=0;

    for(int j=0; j< numExpansion(i); j++) {
      exponent=gaussexp(i, j)*r(1);
      exponent=min(exponent,(doublevar) 60.0);
      exponent=exp(-exponent);
      doublevar n=radiusexp(i,j);
      v_l+= gausscoeff(i, j)
            *pow(r(0),n )
            *exponent;
      dv_l+=gausscoeff(i,j)*exponent*pow(r(0), n-1)*( n
                                                      -2*r(1)*gaussexp(i,j));

      //note that the laplacian is untested..
      d2v_l+=gausscoeff(i,j)*pow(r(0), n-2)*( n*(n-1)-2*gaussexp(i,j)*n*r(1)
                                              -2*(n+1)*gaussexp(i,j)*r(1)
                                              +4*gaussexp(i,j)*gaussexp(i,j)
                                              *r(1)*r(1))*exponent;

    }
    dv_l/=r(0);

    symvals(index,0) = v_l;
    symvals(index,1) = dv_l*r(2);
    symvals(index,2) = dv_l*r(3);
    symvals(index,3) = dv_l*r(4);
    symvals(index,4) = d2v_l+2*dv_l;
    index++;
  }
  
}

//------------------------------------------------------------------------
