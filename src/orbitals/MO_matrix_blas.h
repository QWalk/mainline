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


#ifndef MO_MATRIX_BLAS_H_INCLUDED
#define MO_MATRIX_BLAS_H_INCLUDED

#include "Array.h"
#include "Qmc_std.h"
#include "Basis_function.h"
#include "Center_set.h"
#include "MO_matrix.h"

class System;
class Sample_point;
//----------------------------------------------------------------------------

template <class T> class MO_matrix_blas: public Templated_MO_matrix <T> {
protected:
  void init();
  using Templated_MO_matrix<T>::centers;
  using Templated_MO_matrix<T>::basis;
  using Templated_MO_matrix<T>::nmo;
  using Templated_MO_matrix<T>::magnification_factor;
  using Templated_MO_matrix<T>::orbfile;
  using Templated_MO_matrix<T>::totbasis;
  using Templated_MO_matrix<T>::maxbasis;
  using Templated_MO_matrix<T>::kpoint;
  
private:

  Array1 <Array2 <T> > moCoeff_list;
  Array2 <T> moCoeff;

  Array1 <doublevar> obj_cutoff; //!< cutoff for each basis object
  Array1 <doublevar> cutoff;  //!< Cutoff for individual basis functions
  Array1 <int> nfunctions; //!< number of functions in each basis

public:

  /*!
    Build several sets of MO's to be evaluated in updateVal and updateLap.
    Each element in occupations should be a list of the MO's that should
    be evaluated.  For example, one can create a list of spin up and spin
    down MO's, and only evaluate up when an up electron is moved.
   */
  virtual void buildLists(Array1 <Array1 <int> > & occupations);

  /*!
    get the number of molecular orbitals
   */
  //virtual int getNmo()
  //{
  //  return nmo;
  //}

  virtual int showinfo(ostream & os);

  virtual int writeinput(string &, ostream &);


  virtual void writeorb(ostream &, Array2 <doublevar> & rotation, Array1 <int> &);

  virtual void updateVal(
    Sample_point * sample,
    int e,
    //!< electron number
    int listnum,
    Array2 <T> & newvals
    //!< The return: in form (MO)
  );
  
  virtual void updateLap(
    Sample_point * sample,
    int e,
    //!< electron number
    int listnum,
    //!< Choose the list that was built in buildLists
    Array2 <T> & newvals
    //!< The return: in form ([value gradient lap], MO)
  );

  MO_matrix_blas()
  {}

};


template <class T> struct MOBLAS_CalcObjVal { 
  T * moplace;
  doublevar sval;
};

template <class T> struct MOBLAS_CalcObjLap { 
  T * moplace;
  doublevar sval[5];
};


inline void product_kernel(int n, doublevar & a, doublevar * x, doublevar * y) { 
#ifdef USE_BLAS
  cblas_daxpy(n,a,x,1,y,1);
#else
  for(int i=0; i< n; i++) { 
    *(y+i)+=a* (*(x+i));
  }
#endif
}


inline void product_kernel(int n, doublevar & a, dcomplex * x, dcomplex * y) { 
#ifdef USE_BLAS
  dcomplex b(a);
  cblas_zaxpy(n,&b,x,1,y,1);
#else
  for(int i=0; i< n; i++) { 
    *(y+i)+=a* (*(x+i));
  }  
#endif 
}



template <class T> void MO_matrix_blas<T>::init() {
  //Determine where to cut off the basis functions
  
  cutoff.Resize(totbasis);
  int basiscounter=0;
  for(int i=0; i< centers.size(); i++) {
    for(int j=0; j< centers.nbasis(i); j++) {
      Basis_function* tempbasis=basis(centers.basis(i,j));
      for(int n=0; n< tempbasis->nfunc(); n++) {
        cutoff(basiscounter)=tempbasis->cutoff(n);
        basiscounter++;
      }
    }
  }

  obj_cutoff.Resize(basis.GetDim(0));
  for(int b=0; b< basis.GetDim(0); b++) {
    int nf=basis(b)->nfunc();
    doublevar maxcut=basis(b)->cutoff(0);
    for(int n=1; n< nf; n++) {
      doublevar cut=basis(b)->cutoff(n);
      if(cut > maxcut) maxcut=cut;
    }
    obj_cutoff(b)=maxcut;
  }
  
  
  nfunctions.Resize(basis.GetDim(0));
  for(int b=0; b< basis.GetDim(0); b++) {
    nfunctions(b)=basis(b)->nfunc();
  }
  

  moCoeff.Resize(totbasis, nmo);


  ifstream ORB(orbfile.c_str());
  if(!ORB) error("couldn't find orb file ", orbfile);

  Array3 <int> coeffmat;
  Array1 <T> coeff;
  readorb(ORB,centers, nmo, maxbasis, kpoint,coeffmat, coeff);
  string in;
  ORB.close();

  //Find the cutoffs

  int totfunc=0;

  for(int ion=0; ion<centers.size(); ion++) {
    int f=0;
    doublevar dot=0;
    for(int d=0; d<3; d++) dot+=centers.centers_displacement(ion,d)*kpoint(d);
    
    //cout << "kptfac " << cos(dot*pi) << "  displacement " 
    //    << centers.centers_displacement(ion,0) << "   "
    //    << endl;            
    doublevar kptfac=cos(dot*pi);
        
    for(int n=0; n< centers.nbasis(ion); n++) {

      int fnum=centers.basis(ion,n);
      int imax=basis(fnum)->nfunc();

      for(int i=0; i<imax; i++){ 
        for(int mo=0; mo<nmo; mo++) {	   

          if(coeffmat(mo,ion, f) == -1) {
            moCoeff(totfunc,mo)=0.0;
            //cout << "missing MO pointer: mo# " << mo << " ion # " << ion
            //<< " function on ion: " << f << endl;
            //error("In the orb file, there is a missing pointer. It might "
            //      "be a badly structured file.");
          }
          else moCoeff(totfunc, mo)=magnification_factor*kptfac*coeff(coeffmat(mo,ion,f));
        }//mo
        f++;  //keep a total of functions on center
        totfunc++;
      } //i
    } //n
  }  //ion



}

//----------------------------------------------------------------------


template <class T> void MO_matrix_blas<T>::writeorb(ostream & os, Array2 <doublevar> & rotation, Array1 <int>  &moList) {


  int nmo_write=moList.GetDim(0);
  assert(rotation.GetDim(0)==nmo_write);
  assert(rotation.GetDim(1)==nmo_write);
  os.precision(15);
  int counter=0;
  for(int m=0; m < nmo_write; m++) {
    int mo=moList(m);
    for(int ion=0; ion<centers.size(); ion++) {
      int f=0;
      for(int n=0; n< centers.nbasis(ion); n++) {
        int fnum=centers.basis(ion,n);
        int imax=basis(fnum)->nfunc();

        for(int i=0; i<imax; i++){
          os << mo+1 << "  "   << f+1 << "   " << ion+1 << "   " << counter+1 << endl;
          f++;  //keep a total of functions on center
          counter++;
        } //i
      } //n
    }  //ion
  }
  os << "COEFFICIENTS\n";
  ifstream orbin(orbfile.c_str());
  rotate_orb(orbin, os, rotation, moList, totbasis);
  orbin.close();


}
//---------------------------------------------------------------------


template <class T> void MO_matrix_blas<T>::buildLists(Array1 < Array1 <int> > & occupations)
{
  int numlists=occupations.GetDim(0);
  moCoeff_list.Resize(numlists);
  for(int lis=0; lis < numlists; lis++) {
    int nmo_list=occupations(lis).GetDim(0);
    moCoeff_list(lis).Resize(totbasis, nmo_list);
    for(int i=0; i < nmo_list; i++)  {
      int mo=occupations(lis)(i);
      for(int bas=0; bas < totbasis; bas++) {
        moCoeff_list(lis)(bas,i)=moCoeff(bas,mo);
      }
    }
  }
}


//----------------------------------------------------------------------


template <class T> int MO_matrix_blas<T>::showinfo(ostream & os)
{
  os << "Blas MO " << endl;
  os << "Number of molecular orbitals: " << nmo << endl;
  string indent="  ";
  os << "Basis functions: \n";
  for(int i=0; i< basis.GetDim(0); i++)
  {
    basis(i)->showinfo(indent, os);
  }
  return 1;
}


template <class T> int MO_matrix_blas<T>::writeinput(string & indent, ostream & os)
{
  os << indent << "BLAS_MO" << endl;
  os << indent << "NMO " << nmo << endl;
  os << indent << "ORBFILE " << orbfile << endl;
  os << indent << "MAGNIFY " << magnification_factor << endl;
  string indent2=indent+"  ";
  for(int i=0; i< basis.GetDim(0); i++)
  {
    os << indent << "BASIS { " << endl;
    basis(i)->writeinput(indent2, os);
    os << indent << "}" << endl;
  }

  os << indent << "CENTERS { " << endl;
  centers.writeinput(indent2, os);
  os << indent << "}" << endl;
  return 1;
}
//------------------------------------------------------------------------


template <class T> void MO_matrix_blas<T>::updateVal(Sample_point * sample,
                    int e, int listnum, Array2 <T> & newvals) {

  int centermax=centers.size();

  assert(e < sample->electronSize());
  assert(newvals.GetDim(1) >=1);


  Array1 <doublevar> R(5);
  static Array1 <doublevar> symmvals_temp(maxbasis);
  static Array1 <T> newvals_T;
  Array2 <T> & moCoefftmp(moCoeff_list(listnum));
  int totbasis=moCoefftmp.GetDim(0);
  int nmo_list=moCoefftmp.GetDim(1);
  newvals_T.Resize(nmo_list);
  newvals_T=T(0.0);

  int b;
  int totfunc=0;
  
  centers.updateDistance(e, sample);
  
  static Array1 <MOBLAS_CalcObjVal<T> > calcobjs(totbasis);
  int ncalcobj=0;

  for(int ion=0; ion < centermax; ion++)  {
    centers.getDistance(e, ion, R);
    for(int n=0; n< centers.nbasis(ion); n++) {
      b=centers.basis(ion, n);
      if(R(0) < obj_cutoff(b)) {

        basis(b)->calcVal(R, symmvals_temp);
        int imax=nfunctions(b);
        for(int i=0; i< imax; i++) {
          if(R(0) < cutoff(totfunc)) {
            calcobjs(ncalcobj).sval=symmvals_temp(i);
            calcobjs(ncalcobj).moplace=moCoefftmp.v+totfunc*nmo_list;
            ncalcobj++;
          }
          totfunc++;
        }
      }
      else {
        totfunc+=nfunctions(b);
      }
    }
  }


  for(int i=0; i< ncalcobj; i++) { 
    product_kernel(nmo_list,calcobjs(i).sval,calcobjs(i).moplace,newvals_T.v);
    
  }


  for(int m=0; m < nmo_list; m++) {
    newvals(m,0)=newvals_T(m);

  }  


}

//------------------------------------------------------------------------


/*!
*/

template <class T> void MO_matrix_blas<T>::updateLap(
  Sample_point * sample, 
  int e, 
  int listnum,
  Array2 <T> & newvals) {

  int centermax=centers.size();

  assert(e < sample->electronSize());
  assert(newvals.GetDim(1) >=5);


  Array1 <doublevar> R(5);
  static Array2 <doublevar> symmvals_temp(maxbasis,5);
  static Array2 <T> newvals_T;
  Array2 <T> & moCoefftmp(moCoeff_list(listnum));
  int totbasis=moCoefftmp.GetDim(0);
  int nmo_list=moCoefftmp.GetDim(1);
  newvals_T.Resize(5, nmo_list);
  newvals_T=0.0;

  int b;
  int totfunc=0;
  
  centers.updateDistance(e, sample);

  static Array1 <MOBLAS_CalcObjLap<T> > calcobjs(totbasis);
  int ncalcobj=0;
  

  for(int ion=0; ion < centermax; ion++)  {
    centers.getDistance(e, ion, R);
    for(int n=0; n< centers.nbasis(ion); n++) {
      b=centers.basis(ion, n);
      if(R(0) < obj_cutoff(b)) {
        basis(b)->calcLap(R, symmvals_temp);
        int imax=nfunctions(b);
        for(int i=0; i< imax; i++) {
          if(R(0) < cutoff(totfunc)) {
            calcobjs(ncalcobj).moplace=moCoefftmp.v+totfunc*nmo_list;
            
            for(int j=0; j< 5; j++) {
              calcobjs(ncalcobj).sval[j]=symmvals_temp(i,j);
            }           
            ncalcobj++;
          }
          totfunc++;
        }
      }
      else {
        totfunc+=nfunctions(b);
      }
    }
  }

  for(int i=0; i< ncalcobj; i++) { 
    for(int j=0; j< 5; j++) { 
      product_kernel(nmo_list,calcobjs(i).sval[j],calcobjs(i).moplace,
          newvals_T.v+j*nmo_list);
    }
  }
  

  for(int m=0; m < nmo_list; m++) {
    for(int j=0; j< 5; j++) {
      newvals(m,j)=newvals_T(j,m);
    }
  }

}


#endif // MO_MATRIX_BLAS_H_INCLUDED

//--------------------------------------------------------------------------
