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

#include "Group_function.h"
#include "qmc_io.h"

int Group_function::read(
  vector <string> & words,
  unsigned int & pos)
{

  unsigned int startpos=pos;
  pos=startpos;

  centername=words[0];

  vector < vector < string> > groupstxt;
  vector <string> grouptmp;
  while(readsection(words, pos, grouptmp, "BASIS_GROUP")) {
    groupstxt.push_back(grouptmp);
  }
  nmax=groupstxt.size();
  //cout << nmax;
  if ( nmax < 1 ) {
    error("No group defined in 'BASIS_GROUPS' basis function/construct.");
  }

  groups.Resize(nmax);
  nfunc_groups.Resize(nmax);
  groups=NULL;
  largest_cutoff=0.0;
  nfunc_max=0;
  for(int b=0; b< nmax; b++) {
    allocate(groupstxt[b],groups(b));
    nfunc_groups(b)=groups(b)->nfunc();
    if ( nfunc_groups(b) > nfunc_max ) nfunc_max=nfunc_groups(b);
    for (int i=0; i< nfunc_groups(b); i++) {
      if ( groups(b)->cutoff(i) > largest_cutoff ) {
	largest_cutoff=groups(b)->cutoff(i);
      }
    }
  }

  return 0;
}

void Group_function::getVarParms(Array1 <doublevar> & parms) {
  //cout << "getVarParms " << endl;
  parms.Resize(0);
}

void Group_function::setVarParms(Array1 <doublevar> & parms) {
  //cout << "setVArParms " << endl;
  assert(parms.GetDim(0)==0);
}

int Group_function::nfunc()
{
  return nmax;
}

int Group_function::showinfo(string & indent, ostream & os)
{
  os << indent << "Groups of basis functions" << endl;
  string indent2=indent+"    ";
  for ( int b=0; b<nmax; b++ ) {
    os << indent << "  " << "group " << b+1 << ":" << endl;
    groups(b)->showinfo(indent2,os);
  }
  return 1;
}

int Group_function::writeinput(string & indent, ostream & os)
{
  os << indent << centername << endl;
  os << indent << "BASIS_GROUPS" << endl;
  string indent2=indent+"  ";
  for(int b=0; b< nmax; b++) {
    os << indent << "BASIS_GROUP { " << endl;
    groups(b)->writeinput(indent,os);
    os << indent << "}" << endl;
  }
  return 1;
}

void Group_function::raw_input(ifstream & input)
{
  error("Raw input not supported by Group_function.");
}

void Group_function::calcVal(const Array1 <doublevar> & r,
			     Array1 <doublevar> & symvals,
			     const int startfill)
{
  assert(r.GetDim(0) >= 5);
  assert(symvals.GetDim(0) >= nmax+startfill);
  Array1<doublevar> tmpvals(nfunc_max);
  int index=startfill;
  for(int b=0; b< nmax; b++) {
    groups(b)->calcVal(r,tmpvals,0);
    doublevar val=0.0;
    for (int i=0; i< nfunc_groups(b); i++) val+=tmpvals(i);
    symvals(index++)=val;
  }
}


void Group_function::calcLap(const Array1 <doublevar> & r,
			     Array2 <doublevar> & symvals,
			     const int startfill)
{
  assert(r.GetDim(0) >= 5);
  assert(symvals.GetDim(0) >= nmax+startfill);
  assert(symvals.GetDim(1) >= 5);
  Array2<doublevar> tmpvals(nfunc_max,5);
  int index=startfill;
  for(int b=0; b< nmax; b++) {
    groups(b)->calcLap(r,tmpvals,0);
    Array1<doublevar> lap(5);
    lap=0.0;
    for (int j=0; j< 5; j++) {
      for (int i=0; i< nfunc_groups(b); i++) {
	lap(j)+=tmpvals(i,j);
      }
      symvals(index,j)=lap(j);
    }
    index++;
  }
}



void Group_function::calcHessian(const Array1 <doublevar> & r,
				 Array2 <doublevar> & symvals,
				 const int startfill)
{
  assert(r.GetDim(0) >= 5);
  assert(symvals.GetDim(0) >= nmax+startfill);
  assert(symvals.GetDim(1) >= 10);
  Array2<doublevar> tmpvals(nfunc_max,10);
  int index=startfill;
  for(int b=0; b< nmax; b++) {
    groups(b)->calcHessian(r,tmpvals,0);
    Array1<doublevar> hess(10);
    hess=0.0;
    for (int j=0; j< 10; j++) {
      for (int i=0; i< nfunc_groups(b); i++) {
	hess(j)+=tmpvals(i,j);
      }
      symvals(index,j)=hess(j);
    }
    index++;
  }
}







//------------------------------------------------------------------------
