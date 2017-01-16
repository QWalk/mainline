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

#include "Gen_pade_function.h"
#include "qmc_io.h"

/*

*/
int Gen_pade_function::read(
  vector <string> & words,
  unsigned int & pos
)
{

  cout << "creating genpade" << endl;
  unsigned int startpos=pos;
  centername=words[0];
  rcut=10;   //shouldn't be used at the moment.


  vector <string> alphatxt;
  if(readsection(words, pos=startpos, alphatxt, "ALPHA")) {
    nmax=alphatxt.size();
    alpha.Resize(nmax);
    for(int i=0; i< nmax; i++)
    {
      alpha(i)=atof(alphatxt[i].c_str());
    }
  }
  else {
    error("Need ALPHA in GENPADE function");
  }

  symmetry.Resize(nmax);
  vector <string> symmtxt;
  totfunc=0;
  if(readsection(words, pos=startpos, symmtxt, "SYMMETRY")) {
    if(symmtxt.size() != alphatxt.size()) {
      error("SYMMETRY and ALPHA must have the same size in GENPADE");
    }
    for(int i=0; i< nmax; i++) {
      if(symmtxt[i]=="S") {
        symmetry(i)=0;
        totfunc++;
      }
      else if(symmtxt[i]=="P") {
        totfunc+=3;
        symmetry(i)=1;
      }
      else
        error("Unknown symmetry in SYMMETRY of GENPADE: ", symmtxt[i]);
    }
  }
  else {
    symmetry=0;
    totfunc=nmax;
  }

  return 1;
}

//----------------------------------------------------------------------

void Gen_pade_function::getVarParms(Array1 <doublevar> & parms) {
  parms.Resize(nmax);
  for(int i=0; i< nmax; i++) {
    parms(i)=log(alpha(i));
  }
}

//----------------------------------------------------------------------

void Gen_pade_function::setVarParms(Array1 <doublevar> & parms) {
  assert(parms.GetDim(0)==nmax);

  for(int i=0; i< nmax; i++) {
    alpha(i)=exp(parms(i));
  }

}

//----------------------------------------------------------------------

int Gen_pade_function::nfunc()
{
  return totfunc;
}

//----------------------------------------------------------------------

int Gen_pade_function::showinfo(string & indent, ostream & os)
{
  os << indent << "Generalized Pade function\n";
  os << indent << "Number of functions " << nmax << endl;
  return 1;
}

//----------------------------------------------------------------------

int Gen_pade_function::writeinput(string & indent, ostream & os)
{
  os << indent << centername << endl;
  os << indent << "GENPADE\n";
  os << indent << "ALPHA { ";
  for(int i=0; i< nmax; i++) {
    os << alpha(i) << "   ";
  }
  os << " }\n";

  os << indent << "SYMMETRY { ";
  for(int i=0; i< nmax; i++ ) {
    if(symmetry(i)==0)
      os << "S ";
    else if(symmetry(i)==1)
      os << "P ";
    else error("unknown symmetry in Gen_function::writeinput");
  }
  os << " }\n";
  return 1;
}

//----------------------------------------------------------------------

void Gen_pade_function::raw_input(ifstream & input)
{}


//----------------------------------------------------------------------

void Gen_pade_function::calcVal(const Array1 <doublevar> & r,
                            Array1 <doublevar> & symvals,
                            const int startfill)
{

  //cout << "calcVal " << endl;

  doublevar pade;   //


  int index=startfill;
  doublevar f;
  for(int i=0; i< nmax; i++)
  {
    pade=alpha(i)/(1.+alpha(i)*r(0));
    f=pade*pade*r(1);
    switch(symmetry(i)) {
      case 0:
        symvals(index++)=f;
        break;
      case 1:
        symvals(index++)=f*r(2);
        symvals(index++)=f*r(3);
        symvals(index++)=f*r(4);
        break;
      default:
        error("unknown symmetry in Gen_pade_function::calcVal");
    }

  }

  
  //cout << "done" << endl;
}

//----------------------------------------------------------------------

void Gen_pade_function::calcLap(
  const Array1 <doublevar> & r,
  Array2 <doublevar> & symvals,
  const int startfill
)
{

  doublevar pade,   //
            app,        //alpha(i)*pade*pade
            dadr,       //the first derivative wrt r
            d2adr2;     //laplacian of the radial part(d2f/dr2 + 2/r df/dr)
  doublevar f;
  int index=startfill;
  for(int i=0; i< nmax; i++)
  {
    pade=1.0/(1.+alpha(i)*r(0));
    app=alpha(i)*alpha(i)*pade*pade;
    dadr=app*pade;
    dadr+=dadr;
    d2adr2=3.0*dadr*pade;
    f=app*r(1);
    
    
    switch(symmetry(i)) {
    case 0:
      symvals(index,0)=f;
      for(int j=1; j<4; j++)
      {
        symvals(index,j)=dadr*r(j+1);
      }
      symvals(index,4)=d2adr2;
      index++;
      break;
    case 1:
      doublevar rdp;

      symvals(index,0)=f*r(2); //px
      rdp=dadr*r(2);
      symvals(index,1)=rdp*r(2)+f;
      symvals(index,2)=rdp*r(3);
      symvals(index,3)=rdp*r(4);
      symvals(index,4)=d2adr2*r(2)+2.*rdp;
      index++;

      symvals(index,0)=f*r(3); //py
      rdp=dadr*r(3);
      symvals(index,1)=rdp*r(2);
      symvals(index,2)=rdp*r(3)+f;
      symvals(index,3)=rdp*r(4);
      symvals(index,4)=d2adr2*r(3)+2.*rdp;
      index++;

      symvals(index,0)=f*r(4); //pz
      rdp=dadr*r(4);
      symvals(index,1)=rdp*r(2);
      symvals(index,2)=rdp*r(3);
      symvals(index,3)=rdp*r(4)+f;
      symvals(index,4)=d2adr2*r(4)+2.*rdp;
      index++;
      break;
    default:
      error("unknown symmetry in Gen_pade_function::calcLap: ", symmetry(i));
    }
  }

}

//------------------------------------------------------------------------
