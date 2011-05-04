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

#ifndef MO_MATRIX_H_INCLUDED
#define MO_MATRIX_H_INCLUDED

#include "Array.h"
#include "Qmc_std.h"
#include "Basis_function.h"
#include "Center_set.h"

class System;
class Sample_point;

class General_MO_matrix {
 public:
  /*!
    Build several sets of MO's to be evaluated in updateVal and updateLap.
    Each element in occupations should be a list of the MO's that should
    be evaluated.  For example, one can create a list of spin up and spin
    down MO's, and only evaluate up when an up electron is moved.
   */
  virtual void buildLists(Array1 <Array1 <int> > & occupations)=0;

  /*!
    get the number of molecular orbitals
   */
  virtual int getNmo()=0;
  virtual int showinfo(ostream & os)=0;
  virtual int writeinput(string &, ostream &)=0;
  virtual void read(vector <string> & words, unsigned int & startpos, 
                    System * sys)=0;

  virtual ~General_MO_matrix()  {  }

};

//----------------------------------------------------------------------


template <class T> class Templated_MO_matrix: public General_MO_matrix {
protected:
  Center_set centers;
  Array1 <Basis_function *> basis;
  
  int nmo;
  doublevar magnification_factor;
  string orbfile;
  int totbasis;
  int maxbasis;
  virtual void init() { } ;

  string oldsofile;

  Array1 <doublevar> kpoint; //!< the k-point of the orbitals in fractional units(1 0 0) is X, etc..
  //Array1 <doublevar> mo_counter;//!< Count how many basis functions are evaluated per MO
  //unsigned int n_calls; //!< number of calls
public:

  /*!
    Build several sets of MO's to be evaluated in updateVal and updateLap.
    Each element in occupations should be a list of the MO's that should
    be evaluated.  For example, one can create a list of spin up and spin
    down MO's, and only evaluate up when an up electron is moved.
   */
  virtual void buildLists(Array1 <Array1 <int> > & occupations)=0;
  virtual void setOrbfile(string & x) {
    orbfile=x;
  }
  /*!
    get the number of molecular orbitals
   */
  int getNmo() {
    return nmo;
  }

  virtual int showinfo(ostream & os)=0;

  virtual int writeinput(string &, ostream &)=0;
  virtual void read(vector <string> & words, unsigned int & startpos, 
                    System * sys);

  virtual void writeorb(ostream &, Array2 <doublevar> & rotation, Array1 <int> &)
  { error("writeorb not implemented"); } 

  /*!
    Get the molecular orbital coefficients
   */
  virtual void getMoCoeff(Array2 <T> & coeff) { 
    error("getMoCoeff not implemented");
  }
  
  /*!
    
   */
  virtual void setMoCoeff(Array2 <T> & coeff){
    error("setMoCoeff not implemented");
  }
  virtual int nMoCoeff() { 
    error("nMoCoeff not implemented");
  }


  virtual void updateVal(
    Sample_point * sample,
    int e,
    //!< electron number
    int listnum,
    Array2 <T> & newvals
    //!< The return: in form (MO)
  )=0;

  virtual void getBasisVal(
    Sample_point * sample,
    int e,
    //!< electron number
    Array1 <T> & newvals
    ) { 
    error("getBasisVal not implemented");
  }

  virtual void updateLap(
    Sample_point * sample,
    int e,
    //!< electron number
    int listnum,
    //!< Choose the list that was built in buildLists
    Array2 <T> & newvals
    //!< The return: in form (MO,[value gradient lap])
  )=0;

  virtual void updateHessian(Sample_point * sample,
			     int e,
			     int listnum,
			     Array2 <T>& newvals
			     //!< in form (MO, [value gradient, dxx,dyy,dzz,dxy,dxz,dyz])
			     ) { 
    error("this MO_matrix doesn't support Hessians");
  }

  Templated_MO_matrix()
  {}

  virtual ~Templated_MO_matrix()
  {
    //doublevar totcalls=0;
    //for(int i=0; i< nmo; i++) {
    //  cout << "mo_counter " << mo_counter(i)/n_calls << endl;
    //  totcalls+=mo_counter(i)/n_calls;
    //}
    //cout << " average # of basis functions evaluated " << totcalls/nmo << endl;
    for(int i=0; i< basis.GetDim(0); i++)
      deallocate(basis(i));
  }

};

typedef  Templated_MO_matrix<doublevar> MO_matrix;
typedef  Templated_MO_matrix<dcomplex> Complex_MO_matrix;
//----------------------------------------------------------------------------


int allocate(vector <string> & words, System * sys, MO_matrix *& moptr);
int allocate(vector <string> & words, System * sys, 
             Complex_MO_matrix *& moptr);

void rotate_orb(istream & orbin, ostream & orbout,
                Array2 <doublevar> & rotation,
                Array1 <int>  & moList, int nfunctions);
void rotate_Corb(istream & orbin, ostream & orbout,
		 Array2 <doublevar> & rotation,
		 Array1 <int>  & moList, int nfunctions);


int readorb(istream & input, Center_set & centers, 
            int nmo, int maxbasis, Array1 <doublevar> & kpoint, 
	    Array3 <int > & coeffmat, Array1 <doublevar> & coeff);
int readorb(istream & input, Center_set & centers, 
            int nmo, int maxbasis, Array1 <doublevar> & kpoint, 
	    Array3 <int > & coeffmat, Array1 <dcomplex> & coeff);

//----------------------------------------------------------------------------
#include "qmc_io.h"
template<> inline void Complex_MO_matrix::read(vector <string> & words,
                     unsigned int & startpos,
                     System * sys) { 
}

template<> inline void MO_matrix::read(vector <string> & words,
                     unsigned int & startpos,
                     System * sys)
{


  unsigned int pos=startpos;

  if(!readvalue(words, pos, nmo, "NMO"))
  {
    error("Need NMO in molecular orbital section");
  }

  if(nmo > 40000) 
    error("You have entered more than 40,000 for NMO.  This seems a bit big; we most likely"
        " can't handle it.");


  pos=0;
  if(!readvalue(words, pos, magnification_factor, "MAGNIFY")) {
    magnification_factor=1;
  }


  //Basis functions
  vector <vector <string> > basistext;
  vector < string > basissec;
  pos=0;
  while( readsection(words, pos, basissec, "BASIS") != 0)
  {
    basistext.insert(basistext.end(), basissec);
  }
  basis.Resize(basistext.size());
  basis=NULL;

  if(basistext.size() == 0 )
  {
    error("Didn't find a BASIS section");
  }
  for(unsigned int i=0; i<basistext.size(); i++)
  {
    allocate(basistext[i], basis(i));
  }
  
  sys->kpoint(kpoint);
  /*
  vector <string> kpt;
  if(readsection(words, pos=0, kpt, "KPOINT")) {
    if(kpt.size()!=3) error("KPOINT must be a section of size 3");
    kpoint.Resize(3);
    for(int i=0; i< 3; i++) 
      kpoint(i)=atof(kpt[i].c_str());
  }
  else {
    kpoint.Resize(3);
    kpoint=0;
  }
  */

  //OLDSOFILE is for compatibility to the old qmc f90
  //version, which stores the spline coefficients in
  //a bit of a strange format.  You still have to
  //put in the regular input, since the SO file doesn't
  //contain the symmetry information.
  pos=0;
  if(readvalue(words, pos, oldsofile, "OLDSOFILE"))
  {
    ifstream sofile(oldsofile.c_str());
    if(sofile)
    {
      single_write(cout,
                   "Overwriting the basis functions with the information ",
                   "from ",oldsofile, "\n");
      int dummy;
      sofile >> dummy;
      for(unsigned int i=0; i<basistext.size(); i++)
      {
        basis(i)->raw_input(sofile);
      }
      sofile.close();
    }

    else
    {
      error("Didn't find compatibility SO file in ", oldsofile);
    }
  }


  //------------------------------Centers
  vector <string> centertext;
  pos=startpos;
  if(readsection(words, pos, centertext, "CENTERS"))
  {}
  else
  {
    single_write(cout, "Defaulting to using the atoms as centers\n");
    string temp="USEATOMS";
    centertext.push_back(temp);
  }


  unsigned int newpos=0;
  centers.read(centertext, newpos, sys);
  centers.assignBasis(basis);

  //cout << "number of centers " << centers.size() << endl;
  totbasis=0;
  maxbasis=0;
  for(int i=0; i< centers.size(); i++)
  {
    int basiscent=0;
    for(int j=0; j< centers.nbasis(i); j++)
    {
      basiscent+=basis(centers.basis(i,j))->nfunc();
      //cout << "basiscent " << basiscent << endl;
    }
    totbasis+=basiscent;
    if(maxbasis < basiscent)
      maxbasis=basiscent;
  }

  //cout << "maxbasis " << maxbasis << endl;
  //single_write(cout, nmo, " molecular orbitals requested.\n");

  pos=0;
  if(! readvalue(words, pos, orbfile, "ORBFILE"))
  {
    error("Must specify ORBFILE for MO matrix");
  }
  init();

}
//----------------------------------------------------------------------------

#endif // MO_MATRIX_H_INCLUDED

//----------------------------------------------------------------------------
