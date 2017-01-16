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

#include "Center_set.h"
#include "qmc_io.h"
#include "Basis_function.h"
#include "CBasis_function.h"
#include "Sample_point.h"

//----------------------------------------------------------------------

void Center_set::writeinput(string & indent, ostream & os)
{
  if(usingatoms)
  {
    os << indent << "USEATOMS\n";
  }
  else if(usingsampcenters) 
  {
    os << indent << "USEGLOBAL\n";
  }
  else
  {
    os << indent << "READ " << centerfile << endl;
  }
}

/*!
\todo
Allow center specification inline in the input file, rather than
having to read an external file.
 */
void Center_set::read(vector <string> & words, unsigned int pos,
                      System * sys)
{

  usingsampcenters=usingatoms=0;
  int nelectrons=sys->nelectrons(0)+sys->nelectrons(1);
  unsigned int startpos=pos;

  if(readvalue(words,pos, centerfile, "READ"))
  {

    readcenterfile(centerfile);
  }
  else
  {
    
    pos=startpos;
    if(haskeyword(words, pos,"USEATOMS"))
    {
      //cout << "using atoms " << endl;
      usingatoms=1;
      sys->getAtomicLabels(labels);
      ncenters=labels.size();
      //cout << "ncenters " << ncenters << endl;
    }
    else if(haskeyword(words, pos=0, "USEGLOBAL"))
    {
      usingsampcenters=1;
      sys->getCenterLabels(labels);
      ncenters=labels.size();

      //cout << ncenters << " centers from system " << endl;
    }
    else
    {
      error("Couldn't find anything in center section");
    }

  }

  if(usingsampcenters) {
    sys->getEquivalentCenters(equiv_centers, ncenters_atom, centers_displacement);
  }
  else {
    equiv_centers.Resize(ncenters,1);
    ncenters_atom.Resize(ncenters);
    centers_displacement.Resize(ncenters,3);
    centers_displacement=0;
    for(int cen=0; cen< ncenters; cen++) {
      equiv_centers(cen,0)=cen;
      ncenters_atom(cen)=1;
    }
  }


  edist.Resize(nelectrons,ncenters, 5);
  nbasis.Resize(ncenters);
  nbasis=0;
}

//------------------------------------------------------------

void Center_set::assignBasis(Array1 <Basis_function *> basisfunc)
{
  int totbasis=basisfunc.GetDim(0);
  basis.Resize(ncenters, totbasis);

  for(int i=0; i< ncenters; i++)
  {
    int found=0;
    for(int j=0; j < totbasis; j++)
    {
      if(basisfunc(j)->label() == labels[i])
      {
        appendbasis(j,i);
        found=1;
        //break;
      }

    }
    if(!found)
    {
      cout << "\n****WARNING*****\n"
      << "Couldn't find basis for atom " << labels[i]
      << endl;
    }
  }
}

//------------------------------------------------------------

void Center_set::assignBasis(Array1 <CBasis_function *> basisfunc)
{
  int totbasis=basisfunc.GetDim(0);
  basis.Resize(ncenters, totbasis);

  for(int i=0; i< ncenters; i++)
  {
    int found=0;
    for(int j=0; j < totbasis; j++)
    {
      if(basisfunc(j)->label() == labels[i])
      {
        appendbasis(j,i);
        found=1;
        //break;
      }

    }
    if(!found)
    {
      cout << "\n****WARNING*****\n"
      << "Couldn't find basis for atom " << labels[i]
      << endl;
    }
  }
}


//------------------------------------------------------------

void Center_set::updateDistance(int e, Sample_point * sample)
{
  if(usingatoms)
  {
    Array1 <doublevar> R(5);
    sample->updateEIDist();
    for(int i=0; i< ncenters; i++ )
    {
      sample->getEIDist(e, i, R);
      for(int d=0; d< 5; d++)
      {
        edist(e,i,d)=R(d);
      }
    }
  }
  else if(usingsampcenters) {
    Array1 <doublevar> R(5);
    sample->updateECDist();
    for(int i=0; i< ncenters; i++ )
    {
      sample->getECDist(e, i, R);
      for(int d=0; d< 5; d++)
      {
        edist(e,i,d)=R(d);
      }
    }
  }
  else
  {
    Array1 <doublevar> r(3);
    sample->getElectronPos(e, r);
    for(int i=0; i< ncenters; i++)
    {
      edist(e,i,1)=0;
      for(int d=0; d< 3; d++)
      {
        edist(e,i,d+2)=r(d)-position(i,d);
	//cout << "positions " << position(i,d) << " dist " << edist(e,i,d+2) <<  endl;
        edist(e,i,1)+=edist(e,i,d+2)*edist(e,i,d+2);
      }
    }
    for(int i=0; i< ncenters; i++)
    {
      edist(e,i,0)=sqrt(edist(e,i,1));
    }

  }
}

/*!
Uses flat file of form:
ncenters
Label x y z
Label x y z
...
 */
void Center_set::readcenterfile(string & filename)
{
  assert(!usingatoms);
  ifstream centin(filename.c_str());
  if(!centin) error("Couldn't open ", filename);
  centin >> ncenters;
  position.Resize(ncenters,3);
  string labeltemp;
  for(int i=0; i< ncenters; i++)
  {
    centin >> labeltemp;
    labels.push_back(labeltemp);
    for(int d=0; d< 3; d++)
    {
      centin >> position(i,d);
    }
  }
  centin.close();
}



//--------------------------------------------------------------------------
