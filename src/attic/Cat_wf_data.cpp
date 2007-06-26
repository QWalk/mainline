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

//------------------------------------------------------------------------
//Cat_wf_data.cpp

#include "Qmc_std.h"
#include "qmc_io.h"
#include "Cat_wf_data.h"
#include "Cat_wf.h"

/*!
*/
void Cat_wf_data::read(vector <string> & words, unsigned int & pos,
                       System * sys)
{
  vector <vector <string> > wfsection;
  unsigned int startpos=pos;
  vector <string> indsection;
  while(readsection(words, pos, indsection, "WF"))
  {
    wfsection.push_back(indsection);
  }

  pos=startpos;
  vector <string> orderingtext;
  if( ! readsection(words, pos, orderingtext, "ORDER"))
  {
    int imax=wfsection.size();
    ordering.Resize(imax);
    for(int i=0; i< imax; i++)
    {
      ordering(i)=i;
    }
  }
  else
  {
    int imax=orderingtext.size();
    ordering.Resize(imax);
    for(int i=0; i< imax; i++)
    {
      ordering(i)=atoi(orderingtext[i].c_str());
    }
  }

  nFunc=ordering.GetDim(0);
  nUniqueFunc=wfsection.size();

  data_array.Resize(nUniqueFunc);
  data_array=NULL;
  for(int i=0; i< nUniqueFunc; i++)
  {
    allocate(wfsection[i],sys, data_array(i));
  }

}
//----------------------------------------------------------------------

int Cat_wf_data::supports(wf_support_type support) {
  int ret=1;
  for(int i=0; i< nUniqueFunc; i++)
    ret=ret && data_array(i)->supports(support);
  return ret;

}

//----------------------------------------------------------------------

void Cat_wf_data::generateWavefunction(Wavefunction *& wf)
{
  assert(wf==NULL);
  wf=new Cat_wf;
  Cat_wf * catwf;
  recast(wf, catwf);

  catwf->nUniqueFunc=nUniqueFunc;
  catwf->ordering.Resize(ordering.GetDim(0));
  catwf->ordering=ordering;
  catwf->nFunc=nFunc;
  catwf->wfarray.Resize(nUniqueFunc);
  catwf->wfarray=NULL;
  for(int i=0; i< nUniqueFunc; i++)
  {
    data_array(i)->generateWavefunction(catwf->wfarray(i));
  }
  attachObserver(catwf);
}



void Cat_wf_data::attachObserver(Wavefunction * wf)
{
  Cat_wf * wfptr;
  //cout << "foo " << endl;
  recast(wf, wfptr);
  //cout << "bar " << endl;
  assert(wfptr != NULL);
  for(int i=0; i< nUniqueFunc; i++)
  {
    data_array(i)->attachObserver(wfptr->wfarray(i));
  }
}


int Cat_wf_data::nparms()
{
  int sum=0;
  for(int i=0; i< nUniqueFunc; i++)
  {
    sum+=data_array(i)->nparms();
  }
  return sum;
}

int Cat_wf_data::valSize()
{
  int total=0;
  for(int i=0; i< nUniqueFunc; i++)
  {
    total+=data_array(i)->valSize();
  }
  return total;
}



/*!
 
 */
void Cat_wf_data::getVarParms(Array1 <doublevar> & parms)
{
  Array1 <Array1 <doublevar> > temp_parms(nUniqueFunc);

  int totalsize=0;
  for(int i=0; i< nUniqueFunc; i++)
  {
    data_array(i)->getVarParms(temp_parms(i));
    totalsize+=temp_parms(i).GetDim(0);
  }
  temp_parms.Resize(totalsize);

  int par=0;
  for(int i=0; i< nUniqueFunc; i++)
  {
    for(int j=0; j< temp_parms(i).GetDim(0); j++)
    {
      parms(par)=temp_parms(i)(j);
      par++;
    }
  }

  assert(par==totalsize);

}

void Cat_wf_data::setVarParms(Array1 <doublevar> & parms)
{
  Array1 <Array1 <doublevar> > temp_parms(nUniqueFunc);
  assert( nparms() ==parms.GetDim(0));

  int par=0;
  for(int i=0; i< nUniqueFunc; i++)
  {
    int thisNParms=data_array(i)->nparms();
    temp_parms(i).Resize(thisNParms);
    for(int j=0; j< thisNParms; j++)
    {
      temp_parms(i)(j)=parms(par);
    }
    data_array(i)->setVarParms(temp_parms(i));
  }

  assert(par==nparms());
}

int Cat_wf_data::showinfo(ostream & os)
{
  os << "###############Concatenated function##############\n\n";
  for(int i=0; i< nUniqueFunc; i++)
  {
    data_array(i)->showinfo(os);
  }
  os << endl;
  os << "#############Concatenation Done####################\n";
  return 1;
}

int Cat_wf_data::writeinput(string & indent, ostream & os)
{
  os << indent << "CONCATENATE" << endl;
  os << indent << "ORDER { ";
  for(int i=0; i< nFunc; i++)
  {
    os << ordering(i) << "  ";
  }
  os << " } " << endl;

  string indent2=indent+"  ";

  for(int i=0; i< nUniqueFunc; i++)
  {
    os << indent << "WF { " << endl;
    data_array(i)->writeinput(indent2, os);
    os << indent << "}" << endl;
  }
  return 1;
}


//------------------------------------------------------------------------
