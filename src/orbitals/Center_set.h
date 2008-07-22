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

#ifndef CENTER_SET_H_INCLUDED
#define CENTER_SET_H_INCLUDED
#include "Array.h"
#include "System.h"
class Basis_function;
class CBasis_function;
class Sample_point;
/*!
\brief
Holds a set of centers for basis function evaluation.

Input specification:

Required one of:

<b> USEATOMS : </b> Sets the centers to the positions of the
atoms given in a System

--or--

<b> READ : </b> takes a single string as argument, attempts to
read center positions from there.  The format is:

ncenters                 <br>
label1  x1 y1 z1         <br>
label2  x2 y2 z2         <br>
.                        <br>
.                        <br>
.                        <br>

*/
class Center_set
{
public:

  Array2 <int> basis;
  //!< Holds the index of the basis functions associated with each particle, of form (particle, function on that particle)

  Array1 <int> nbasis;
  //!< Number of basis functions on each particle

  Center_set()
  { usingatoms=0;}

  void read(vector <string> & words, unsigned int pos,
            System * sys);
  void assignBasis(Array1 <Basis_function *>);
  void assignBasis(Array1 <CBasis_function *>);

  void updateDistance(int, Sample_point *);
  void getDistance(const int e, const int cent,
                   Array1 <doublevar> & distance)
  {
    assert(e< edist.GetDim(0));
    assert(cent < ncenters);
    assert(distance.GetDim(0) >=5);
    for(int d=0; d< 5; d++)
    {
      distance(d)=edist(e,cent,d);
    }
  }


  void writeinput(string &, ostream & );

  void Resize(int mdim, int mpar)
  {
    //r.Resize(mdim, mpar);
    nbasis.Resize(mpar);
  }


  int size() const
  {
    return nbasis.GetDim(0);
  }

  void readcenterfile(string &);
  //------------------------------------------------------------
  //Basis function methods

  void initbasis(int m)
  {
    basis.Resize(nbasis.GetDim(0),m);
    nbasis=0;
  }


  //assigns a basis function to a specific particle in the
  //set
  void assignbasis(int b, int n, int m)
  {
    basis(n,m)=b;
    if(nbasis(n) < m)
      nbasis(n)=m;
  }
  void appendbasis(int b,int n)
  {
    nbasis(n)++;
    basis(n,nbasis(n)-1)=b;
  }
  
  Array2 <int> equiv_centers; //!< for a given atom, what centers correspond to it
  Array1 <int> ncenters_atom; //!< number of centers on a given atom
  Array2 <int> centers_displacement; //!< displacement vector of a center(for periodic systems and k-points)
private:
  string centerfile;
  int ncenters;
  int usingatoms;
  int usingsampcenters;
  Array2 <doublevar> position;
  Array3 <doublevar> edist;
  
  vector <string> labels;


};

#endif //CENTER_SET_H_INCLUDED

//--------------------------------------------------------------------------
