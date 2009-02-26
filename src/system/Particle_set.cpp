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


#include "qmc_io.h"
#include <iomanip>
#include "Particle_set.h"

using namespace std;


int Particle_set::showinfo(ostream & os)
{
  int field=os.precision()+10;
  os << "Atomic positions \n";
  os << setw(6) << "Label";
  os << setw(field) << "Charge"
  << setw(field) << "   x   "
  << setw(field) << "   y   "
  << setw(field) << "   z   " << endl;
  for(int i=0; i< r.GetDim(1); i++)
  {
    os << setw(6) << labels[i];
    os << setw(field) << charge(i);
    for(int j=0; j<r.GetDim(0); j++)
    {
      os << setw(field) << r(j,i);
    }
    os << endl;
  }

  return 1;
}


int Particle_set::read(vector <string> & words, unsigned int pos) {
  unsigned int startpos=0;
  vector <string> atomsec;
  vector <Atom> atoms;
  pos=startpos;
  while(readsection(words,pos,atomsec, "ATOM") != 0)
  {
    Atom tempatom;
    tempatom.name=atomsec[0];
    tempatom.charge=atof(atomsec[1].c_str());
    unsigned int pos2=0;
    if(!readvalue(atomsec,pos2, tempatom.coor[0], "COOR"))
      error("Looking for COOR in ATOM section; didn't find it");
    readnext(atomsec, pos2, tempatom.coor[1]);
    readnext(atomsec, pos2, tempatom.coor[2]);
    pos2=0;
    atoms.push_back(tempatom);
  }

  int nions=atoms.size();
  r.Resize(3, nions);
  charge.Resize(nions);
  for(int i=0; i<nions; i++)
  {
    //cout << "ion " <<  i << endl;
    for(int j=0; j<3;j++)
    {
      r(j,i)=atoms[i].coor[j];
      //cout << r(j,i) << "   ";
    }
    charge(i)=atoms[i].charge;
    labels.push_back(atoms[i].name);
    //cout << "\n charge " << charge(i) << endl;
  }
  return 1;

}
//--------------------------------------------------------------------------
