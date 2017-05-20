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

#ifndef MO_MATRIX_CUTOFF_H_INCLUDED
#define MO_MATRIX_CUTOFF_H_INCLUDED

#include "Array.h"
#include "Qmc_std.h"
#include "Basis_function.h"
#include "Center_set.h"
#include "MO_matrix.h"

class System;
class Sample_point;
//----------------------------------------------------------------------------

template <class T> class MO_matrix_cutoff: public Templated_MO_matrix <T> {
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
  //Center_set centers;
  //Array1 <Basis_function *> basis;
  //int nmo;
  //int totbasis;
  //int maxbasis;
 // doublevar magnification_factor;
  //string orbfile;
  Array2 <int> mofill;
  Array2 <T> moCoeff2;
  Array1 <int> nbasis;

  Array1 <doublevar> obj_cutoff; //!< cutoff for each basis object
  Array1 <doublevar> cutoff;  //!< Cutoff for individual basis functions
  Array1 <int> nfunctions; //!< number of functions in each basis
  //Array1 <int> basismo;
  //Array2 <doublevar> moCoeff;
  //Array2 <int> basisfill;

  Array1 < Array2 <int> > basisfill_list;
  Array1 < Array2 <T> > moCoeff_list;
  Array1 < Array1 <int> > basismo_list;

 Array1 <doublevar> symmvals_temp1d;
 Array2 <doublevar> symmvals_temp2d;



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


  //virtual void read(vector <string> & words, unsigned int & startpos, System * sys);

  //! Takes an ORB file and inserts all the coefficients.
  //virtual int readorb(istream &);


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
    int listnum,
    Array2 <T> & newvals
  );
  virtual void updateHessian(Sample_point * sample,
			     int e,
			     int listnum,
			     Array2 <T>& newvals
			     //!< in form ([value gradient, dxx,dyy,dzz,dxy,dxz,dyz], MO)
			     );
  MO_matrix_cutoff()
  {}

};

//######################################################################

#include "Qmc_std.h"
#include "MO_matrix_cutoff.h"
#include "Sample_point.h"
#include "qmc_io.h"



template <class T> void MO_matrix_cutoff<T>::init() {

  
  //Determine where to cut off the basis functions
  
  cutoff.Resize(totbasis);
  int basiscounter=0;
  for(int i=0; i< centers.size(); i++)
  {
    for(int j=0; j< centers.nbasis(i); j++)
    {
      Basis_function* tempbasis=basis(centers.basis(i,j));
      for(int n=0; n< tempbasis->nfunc(); n++)
      {
        //cout << "cutoff " << endl;
        cutoff(basiscounter)=tempbasis->cutoff(n);
        //cout << "rcut(" << basiscounter << ") "
        //     << cutoff(basiscounter) << endl;
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
  

  nbasis.Resize(nmo);
  mofill.Resize(nmo, totbasis);
  moCoeff2.Resize(nmo, totbasis);

  ifstream ORB(orbfile.c_str());

  if(!ORB) {
    error("couldn't find orb file ", orbfile);
  }

  Array3 <int> coeffmat;
  Array1 <T> coeff;
  
  readorb(ORB,centers, nmo, maxbasis,kpoint, coeffmat, coeff);
  string in;
  ORB.close();

  
  //Find the cutoffs

  int totfunc=0;
  nbasis=0;
  //basismo=0;
  const doublevar threshold=1e-12;
  for(int ion=0; ion<centers.size(); ion++)
  {
    int f=0;

    doublevar dot=0;
    for(int d=0; d<3; d++) dot+=centers.centers_displacement(ion,d)*kpoint(d);
    
    //T kptfac=T(exp(dcomplex(0,1.0)*dot*pi));
    T kptfac=eval_kpoint_fac<T>(dot);
    
    for(int n=0; n< centers.nbasis(ion); n++) {
      
      int fnum=centers.basis(ion,n);
      int imax=basis(fnum)->nfunc();

      for(int i=0; i<imax; i++) { //sum over the symmetries
        for(int mo=0; mo<nmo; mo++) {      //and the MO's
          T temp;
          if(coeffmat(mo,ion, f) == -1) {
            temp=T(0.0);
            //cout << "missing MO pointer: mo# " << mo << " ion # " << ion
            //<< " function on ion: " << f << endl;
            //error("In the orb file, there is a missing pointer. It might "
            //      "be a badly structured file.");
          }
          else temp=coeff(coeffmat(mo,ion,f));
          if(abs(temp) > threshold) {
            mofill(mo, nbasis(mo))=totfunc;
            moCoeff2(mo, nbasis(mo))=kptfac*magnification_factor*temp;
            nbasis(mo)++;
          }

        }//mo
        f++;  //keep a total of functions on center
        totfunc++;
      } //i
    } //n
  }  //ion
  symmvals_temp1d.Resize(maxbasis);
  symmvals_temp2d.Resize(maxbasis,10);

}

//---------------------------------------------------------------------------------------------

template<class T> void MO_matrix_cutoff<T>::writeorb(ostream & os, 
    Array2 <doublevar> & rotation, Array1 <int>  &moList) {


  int nmo_write=moList.GetDim(0);
  assert(rotation.GetDim(0)==nmo_write);
  assert(rotation.GetDim(1)==nmo_write);

  os.precision(15);
  int counter=0;
  int totfuncs=0;
  for(int ion=0; ion<centers.equiv_centers.GetDim(0); ion++) {
    for(int n=0; n< centers.nbasis(centers.equiv_centers(ion,0)); n++) {
      int fnum=centers.basis(centers.equiv_centers(ion,0),n);
      totfuncs+=basis(fnum)->nfunc();
    }
  }
  
  for(int m=0; m < nmo_write; m++) {
    for(int ion=0; ion<centers.equiv_centers.GetDim(0); ion++) {
      int f=0;
      for(int n=0; n< centers.nbasis(centers.equiv_centers(ion,0)); n++) {
        int fnum=centers.basis(centers.equiv_centers(ion,0),n);
        int imax=basis(fnum)->nfunc();

        for(int i=0; i<imax; i++) {
          os << m+1 << "  "   << f+1 << "   " << ion+1 << "   " << counter+1 << endl;
          f++;  //keep a total of functions on center
          counter++;
        } //i
      } //n
    }  //ion
  }
  os << "COEFFICIENTS\n";
  ifstream orbin(orbfile.c_str());

  //cout << "orbfile " << orbfile << endl;
  Array1 <T> coeff;
  Array3 <int> coeffmat;
  readorb(orbin,centers, nmo, maxbasis,kpoint, coeffmat, coeff);
  orbin.close();

  Array2 <T> moCoeff(nmo,totbasis);
  moCoeff=0.0;
  int nmo_read=coeffmat.GetDim(0);
  int ncenter=coeffmat.GetDim(1);
  int maxfunc=coeffmat.GetDim(2);
  int currfunc;
  for(int mo=0; mo < nmo_read; mo++) { 
    currfunc=0;
    for(int ion=0; ion<centers.equiv_centers.GetDim(0); ion++) {
      int f=0;
      int equiv_center=centers.equiv_centers(ion,0);
      for(int n=0; n< centers.nbasis(centers.equiv_centers(ion,0)); n++) {
        int fnum=centers.basis(centers.equiv_centers(ion,0),n);
        int imax=basis(fnum)->nfunc();
        for(int i=0; i<imax; i++) {
          if(coeffmat(mo,equiv_center,f)!=-1)
            moCoeff(mo,currfunc)=coeff(coeffmat(mo,equiv_center,f));
          f++;
          currfunc++;
        }
      }
    }
  }
  // Now we rotate and write out the coefficients
  Array2 <T> rotatedMO(nmo_write,currfunc);
  rotatedMO=0.0;
  for(int mo=0; mo < nmo_write; mo++) { 
    for(int f=0; f< currfunc; f++) { 
      for(int mo2=0; mo2 < nmo_write; mo2++) { 
        int realmo=moList(mo2);
        rotatedMO(mo,f)+=rotation(mo,mo2)*moCoeff(realmo,f);
        //cout << "mo1 " << mo << " mo2 " << mo2 << " f " << f 
        //  << " realmo " << realmo << " rotation " << rotation(mo,mo2) 
        //     << " rotatedmo " << rotatedMO(mo,f) << " coeff " << moCoeff(realmo,f) << endl;
      }
    }
  }
  int counter2=1;
  for(int m=0; m < nmo_write; m++) {
    for(int f=0; f< currfunc; f++) {
      os << rotatedMO(m, f) << "   ";
      if(counter2 % 5 ==0) os << endl;
      counter2++;
    }
  }
}
//---------------------------------------------------------------------

template <class T> void MO_matrix_cutoff<T>::buildLists(Array1 < Array1 <int> > & occupations){
  int numlists=occupations.GetDim(0);
  basisfill_list.Resize(numlists);
  moCoeff_list.Resize(numlists);
  basismo_list.Resize(numlists);
  for(int lis=0; lis < numlists; lis++)
  {
    int nmo_list=occupations(lis).GetDim(0);
    basisfill_list(lis).Resize(totbasis, nmo_list);
    moCoeff_list(lis).Resize(totbasis, nmo_list);
    basismo_list(lis).Resize(totbasis);
    basismo_list(lis)=0;
    for(int i=0; i < nmo_list; i++)
    {
      int mo=occupations(lis)(i);
      for(int bas=0; bas < nbasis(mo); bas++)
      {
        int func=mofill(mo, bas);

        //basisfill_list(lis)(func, basismo_list(lis)(func))=mo;
        basisfill_list(lis)(func, basismo_list(lis)(func))=i;
        moCoeff_list(lis)(func, basismo_list(lis)(func))=moCoeff2(mo, bas);
        //cout << "basisfill_list " << 2 << "  f  "
        //     << func << "  mo " <<  mo;
        //cout << "  basis coeff " << moCoeff2(mo, bas) << endl;
        //cout << "real basisfill " << basisfill(func, basismo_list(lis)(func))
        //     << " coeff  " << moCoeff(func, basismo_list(lis)(func)) << endl;
        basismo_list(lis)(func)++;
      }
    }
  }
}


//----------------------------------------------------------------------

template <class T> int MO_matrix_cutoff<T>::showinfo(ostream & os)
{
  os << "Cutoff MO " << endl;
  os << "Number of molecular orbitals: " << nmo << endl;
  string indent="  ";
  os << "Basis functions: \n";
  for(int i=0; i< basis.GetDim(0); i++)
  {
    basis(i)->showinfo(indent, os);
  }
  return 1;
}

template <class T> int MO_matrix_cutoff<T>::writeinput(string & indent, ostream & os)
{
  os << indent << "CUTOFF_MO" << endl;
  os << indent << "NMO " << nmo << endl;
  os << indent << "ORBFILE " << orbfile << endl;
  //if(oldsofile!="") 
  //  os << indent << "OLDSOFILE " << oldsofile << endl;
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

template <class T> void MO_matrix_cutoff<T>::updateVal(
  Sample_point * sample,  int e,  int listnum,  Array2 <T> & newvals) {
  //cout << "start updateval " << endl;
  int centermax=centers.size();
  static Array1 <doublevar> R(5);
  //Array1 <doublevar> symmvals_temp(maxbasis);

  //Make references for easier access to the list variables.
  Array1 <int> & basismotmp(basismo_list(listnum));
  Array2 <int> & basisfilltmp(basisfill_list(listnum));
  Array2 <T> & moCoefftmp(moCoeff_list(listnum));
  assert(newvals.GetDim(1) >= 1);

  newvals=0;
  Basis_function * tempbasis;

  //int fn;
  T c;
  int mo=0;
  int scalebasis=basisfill_list(listnum).GetDim(1);
  int totfunc=0;
  int b; //basis
  //cout << "here " << endl;
  centers.updateDistance(e, sample);
  //int retscale=newvals.GetDim(1);
  for(int ion=0; ion < centermax; ion++) {
    //sample->getECDist(e, ion, R);
    centers.getDistance(e,ion,R);
    for(int n=0; n< centers.nbasis(ion); n++) {
      b=centers.basis(ion,n);
      tempbasis=basis(b);
      if(obj_cutoff(b) > R(0)) {
        tempbasis->calcVal(R, symmvals_temp1d);
        //cout << "ion " << ion << "b " << b << " mo "<< mo << endl;
        int imax=nfunctions(b);
        for(int i=0; i< imax; i++) {
          int reducedbasis=scalebasis*totfunc;
          if(R(0) < cutoff(totfunc)) {
            for(int basmo=0; basmo < basismotmp.v[totfunc]; basmo++) {
              //mo=basisfill(totfunc, basmo);
              //c=moCoeff(totfunc, basmo);
              //cout << "basisfill reducedbasis "<< reducedbasis
              // << "  basmo " << basmo << endl;
              mo=basisfilltmp.v[reducedbasis+basmo];
              //cout << "mocoeff (mo=" << mo <<  endl;
              //mo_counter(mo)++;
              c=moCoefftmp.v[reducedbasis+basmo];

              newvals(mo, 0)+=c*symmvals_temp1d(i);
              //newvals.v[retscale*mo]+=c*symmvals_temp.v[i];

            }
          }
          totfunc++;
        }
      }
      else {
        totfunc+=nfunctions(b);
      }
    }
  }
  //n_calls++;
    //cout << "done updateVal " << endl;
}

//------------------------------------------------------------------------

template <class T>void MO_matrix_cutoff<T>::updateLap( Sample_point * sample,
  int e, int listnum, Array2 <T> & newvals) {

  //cout << "updateLap" << endl;
  int centermax=centers.size();
  //int momax=occupation.GetDim(0);
  newvals=0;
  //assert(momax <= nmo);
  assert(e < sample->electronSize());
  assert(newvals.GetDim(1) >=5);

  // cout << "array " << endl;
  Array1 <doublevar> R(5);
  // cout << "symvals " << endl;
   //cout << "arrayref " << endl;

  //References to make the code easier to read and slightly faster.
  Array1 <int> & basismotmp(basismo_list(listnum));
  Array2 <int> & basisfilltmp(basisfill_list(listnum));
  Array2 <T> & moCoefftmp(moCoeff_list(listnum));

  Basis_function * tempbasis;

  T c;
  int scaleval=0, scalesymm=0;
  int mo=0;
  int scalebasis=basisfilltmp.GetDim(1);
  centers.updateDistance(e, sample);
  int totfunc=0;
  int b;
  int symmvals_stride=symmvals_temp2d.GetDim(1);
  for(int ion=0; ion < centermax; ion++) {
    centers.getDistance(e, ion, R);
    for(int n=0; n< centers.nbasis(ion); n++) {
      b=centers.basis(ion, n);
      tempbasis=basis(b);
      if(R(0) < obj_cutoff(b)) {
       // cout << "basis " << endl;
      tempbasis->calcLap(R, symmvals_temp2d);

      int imax=nfunctions(b);
      for(int i=0; i< imax; i++) {
        //cout << "i " << i << endl;

        int reducedbasis=scalebasis*totfunc;
        scalesymm=i*symmvals_stride;
        if(R(0) < cutoff(totfunc)) {
          for(int basmo=0; basmo < basismotmp.v[totfunc]; basmo++) {
            //mo=basisfill(basis, basmo);
            //c=moCoeff(basis, basmo);

            mo=basisfilltmp.v[reducedbasis+basmo];
            c=moCoefftmp.v[reducedbasis+basmo];
            //cout << c << "   ";

            //mo_counter(mo)++;

            scaleval=mo*5;
            //cout << "coeff " << c << endl;
            //cout << "reducedbasis " << reducedbasis << endl;
            // cout << mo << "   " << basmo << "   " << c << endl;


            for(int j=0; j< 5; j++) {
              //newvals(mo,j)+=c*symmvals(fn,j);
              newvals.v[scaleval+j]+=c*symmvals_temp2d.v[scalesymm+j];
              //cout << "newvals(" << mo << "," << j << ")  " << newvals(mo,j) << endl;
              //cout << "symmvals(" << j << ")  " << symmvals(fn,j) << endl;
            }
          }
          //cout << endl;

        }

        totfunc++;
      }
    }
    else {
      totfunc+=nfunctions(b);
    }
    }
  }


  //n_calls++;
  
  //cout << "newvals " << endl;
  //output_array(newvals);

}

//--------------------------------------------------------------------------

template <class T>void MO_matrix_cutoff<T>::updateHessian(
  Sample_point * sample,
  int e,
  int listnum,
  //const Array1 <int> & occupation,
  //!<A list of the MO's to evaluate
  Array2 <T> & newvals
  //!< The return: in form (MO, [val, grad, dxx,dyy,...])
)
{

  int centermax=centers.size();
  newvals=0;
  assert(e < sample->electronSize());
  assert(newvals.GetDim(1)==10);
  

  Array1 <doublevar> R(5);
 // static Array2 <doublevar> symmvals_temp(maxbasis,10);

  //References to make the code easier to read and slightly faster.
  Array1 <int> & basismotmp(basismo_list(listnum));
  Array2 <int> & basisfilltmp(basisfill_list(listnum));
  Array2 <T> & moCoefftmp(moCoeff_list(listnum));

  Basis_function * tempbasis;

  T c;
  int scaleval=0, scalesymm=0;
  int mo=0;
  int scalebasis=basisfilltmp.GetDim(1);
  centers.updateDistance(e, sample);
  int totfunc=0;
  int b;
  int symmvals_stride=symmvals_temp2d.GetDim(1);
  for(int ion=0; ion < centermax; ion++)  {
    centers.getDistance(e, ion, R);
    for(int n=0; n< centers.nbasis(ion); n++)  {
      b=centers.basis(ion, n);
      tempbasis=basis(b);
      if(R(0) < obj_cutoff(b)) {
        tempbasis->calcHessian(R, symmvals_temp2d);

        int imax=nfunctions(b);
        for(int i=0; i< imax; i++)  {
          int reducedbasis=scalebasis*totfunc;
          scalesymm=i*10;
          if(R(0) < cutoff(totfunc))  {
            for(int basmo=0; basmo < basismotmp.v[totfunc]; basmo++)   {
              mo=basisfilltmp.v[reducedbasis+basmo];
              c=moCoefftmp.v[reducedbasis+basmo];
              scaleval=mo*symmvals_stride;
              for(int j=0; j< 10; j++) {
                newvals.v[scaleval+j]+=c*symmvals_temp2d.v[scalesymm+j];
              }
            }
          }

          totfunc++;
        }
      }
      else {
        totfunc+=nfunctions(b);
      }
    }
  }
}

//--------------------------------------------------------------------------



#endif // MO_MATRIX_CUTOFF_H_INCLUDED

//--------------------------------------------------------------------------
