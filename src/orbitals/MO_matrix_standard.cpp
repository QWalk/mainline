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

#include "Qmc_std.h"
#include "MO_matrix_standard.h"
#include "Sample_point.h"
#include "qmc_io.h"



void MO_matrix_standard::init() {
  moCoeff.Resize(nmo, totbasis);

  single_write(cout, "Standard MO\n");
  ifstream ORB(orbfile.c_str());
  if(!ORB)
  {
    error("couldn't find orb file ", orbfile);
  }

  single_write(cout,"Reading orbitals from ",orbfile, "\n");
  Array3 <int> coeffmat;
  Array1 <doublevar> coeff;
  int tempint=readorb(ORB,centers,nmo, maxbasis,kpoint, coeffmat, coeff);
  ORB.close();
  single_write(cout, tempint," unique MO coefficients found.\n\n");

  //Fill moCoeff

  int totfunc=0;
  for(int ion=0; ion<centers.size(); ion++)
  {
    int f=0;
    doublevar dot=0;
    for(int d=0; d<3; d++) dot+=centers.centers_displacement(ion,d)*kpoint(d);
    
    cout << "kptfac " << cos(dot*pi) << "  displacement " 
        << centers.centers_displacement(ion,0) << "   "
        << endl;            
    doublevar kptfac=cos(dot*pi);
    
    for(int n=0; n< centers.nbasis(ion); n++)
    {

      //evaluate atomic orbital
      int fnum=centers.basis(ion,n);

      int imax=basis(fnum)->nfunc();

      for(int i=0; i<imax; i++)
      { //sum over the symmetries
        for(int mo=0; mo<nmo; mo++)
        {	      //and the MO's
          //cout << "ion " << ion;
          //cout << "  i " << i << " mo " << mo << "  fnum " << fnum << endl;
          //cout << "coeffmat " << coeffmat(mo,ion, f) << endl;
          //moCoeff(ion, f, mo)=coeff(coeffmat(mo,ion,f));
          if(coeffmat(mo,ion, f) == -1) {
            cout << "missing MO pointer: mo# " << mo << " ion # " << ion
            << " function on ion: " << f << endl;
            error("In the orb file, there is a missing pointer. It might "
                  "be a badly structured file.");
          }
          doublevar temp=coeff(coeffmat(mo,ion,f));
          moCoeff(mo, totfunc)=magnification_factor*kptfac*temp;

        }//mo
        f++;  //keep a total of functions on center
        totfunc++;
      } //i
    } //n
  }  //ion

}

//----------------------------------------------------------------------

void MO_matrix_standard::writeorb(ostream & os, 
                                  Array2 <doublevar> & rotation, 
                                  Array1 <int>  &moList) {

  int nmo_write=moList.GetDim(0);
  assert(rotation.GetDim(0)==nmo_write);
  assert(rotation.GetDim(1)==nmo_write);
  cout << "nmo_write " << nmo_write << endl;
 
  int counter=0;
  for(int m=0; m < nmo_write; m++) {
    int mo=moList(m);
    for(int ion=0; ion<centers.size(); ion++)
    {
      int f=0;

      for(int n=0; n< centers.nbasis(ion); n++)
      {
        int fnum=centers.basis(ion,n);
        int imax=basis(fnum)->nfunc();

        for(int i=0; i<imax; i++)
        {
          os << mo+1 << "  "   << f+1 << "   " << ion+1 << "   " 
             << counter+1 << endl;
          f++;  //keep a total of functions on center
          counter++;
        } //i
      } //n
    }  //ion
  }
  os << "COEFFICIENTS\n";

  Array2 <doublevar> rotated_mo(nmo_write, totbasis);
  rotated_mo=0;
  for(int m=0; m< nmo_write; m++) {
    for(int m2=0; m2 < nmo_write; m2++) {
      for(int f=0; f< totbasis; f++) {
        rotated_mo(m, f)+=rotation(m,m2)*moCoeff(m2,f);
      }
    }
  }

  int counter2=1;
  for(int m=0; m < nmo_write; m++) {
    for(int f=0; f< totbasis; f++) {
      os << rotated_mo(m, f) << "   ";
      if(counter2 % 5 ==0) os << endl;
    }
  }

}
//---------------------------------------------------------------------

void MO_matrix_standard::buildLists(Array1 < Array1 <int> > & occupations) {
  int numlists=occupations.GetDim(0);
  moLists.Resize(numlists);
  for(int lis=0; lis < numlists; lis++) {
    int nmo_list=occupations(lis).GetDim(0);
    moLists(lis).Resize(nmo_list);
    for(int mo=0; mo < nmo_list; mo++) {
      moLists(lis)(mo)=occupations(lis)(mo);
    }
  }
}


//----------------------------------------------------------------------

int MO_matrix_standard::showinfo(ostream & os)
{
  os << "Standard molecular orbital\n";
  string indent="  ";
  os << "Basis functions: \n";
  for(int i=0; i< basis.GetDim(0); i++) {
    basis(i)->showinfo(indent, os);
  }

  os << "Number of molecular orbitals: " << nmo << endl;

  return 1;
}

int MO_matrix_standard::writeinput(string & indent, ostream & os)
{
  os << indent << "STANDARD_MO" << endl;
  os << indent << "NMO " << nmo << endl;
  os << indent << "ORBFILE " << orbfile << endl;
  os << indent << "MAGNIFY " << magnification_factor << endl;
  string indent2=indent+"  ";
  for(int i=0; i< basis.GetDim(0); i++) {
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

void MO_matrix_standard::updateVal(Sample_point * sample, int e,
                                   int listnum,
                                   //!<which list to use
                                   Array2 <doublevar> & newvals
                                   //!< The return: in form (MO, val)
) {

  Array1 <doublevar> basisvals(totbasis);
  int centermax=centers.size();
  centers.updateDistance(e, sample);
  Basis_function * tempbasis;
  Array1 <doublevar> R(5);
  int currfunc=0;
  newvals=0;
  for(int ion=0; ion < centermax; ion++) {
    centers.getDistance(e, ion, R);
    for(int n=0; n< centers.nbasis(ion); n++)
    {

      tempbasis=basis(centers.basis(ion,n));
      tempbasis->calcVal(R, basisvals, currfunc);
      currfunc+=tempbasis->nfunc();
    }
  }

  int nmo_list=moLists(listnum).GetDim(0);
  //int valscale=newvals.GetDim(1);
  int mocoeffscale=moCoeff.GetDim(1);
  doublevar * tempnewvals=new doublevar[moCoeff.GetDim(0)];
  for(int m=0; m < nmo_list; m++) {
    int mo=moLists(listnum)(m);
    int effmo=mo*mocoeffscale;
    //int effm=m*valscale;
    tempnewvals[m]=0;
    for(int f=0; f< totbasis; f++) {
      //c=moCoeff(mo,f);
      //newvals(m,0)+=c*basisvals(f);
      //newvals.v[effm]+=moCoeff.v[effmo+f]*basisvals.v[f];
      tempnewvals[m]+=moCoeff.v[effmo+f]*basisvals.v[f];
    }
  }

  int nfill=moCoeff.GetDim(0);
  for(int m=0; m<nfill; m++) {
    newvals(m,0)=tempnewvals[m];
  }

  delete [] tempnewvals;
}

//------------------------------------------------------------------------

void MO_matrix_standard::getBasisVal(Sample_point * sample, int e,
				     Array1 <doublevar> & newvals
				     ) {

  newvals.Resize(totbasis);
  int centermax=centers.size();
  centers.updateDistance(e, sample);
  Basis_function * tempbasis;
  Array1 <doublevar> R(5);
  int currfunc=0;
  newvals=0;
  for(int ion=0; ion < centermax; ion++) {
    centers.getDistance(e, ion, R);
    for(int n=0; n< centers.nbasis(ion); n++)
    {

      tempbasis=basis(centers.basis(ion,n));
      tempbasis->calcVal(R, newvals, currfunc);
      currfunc+=tempbasis->nfunc();
    }
  }
}

//------------------------------------------------------------------------


/*!
*/

void MO_matrix_standard::updateLap(
  Sample_point * sample, int e,  int listnum,
  Array2 <doublevar> & newvals
  //!< The return: in form (MO, [val, grad, lap])
){
  Array2 <doublevar> basisvals(totbasis, 5);
  int centermax=centers.size();
  centers.updateDistance(e, sample);
  Basis_function * tempbasis;
  Array1 <doublevar> R(5);
  int currfunc=0;
  newvals=0;
  //cout << "centermax " << centermax << endl;
  for(int ion=0; ion < centermax; ion++) {
    centers.getDistance(e, ion, R);
    //cout << "nbasis " << centers.nbasis(ion) << endl;
    for(int n=0; n< centers.nbasis(ion); n++)
    {

      tempbasis=basis(centers.basis(ion,n));
      tempbasis->calcLap(R, basisvals, currfunc);
      currfunc+=tempbasis->nfunc();
    }
  }

  doublevar c;
  int nmo_list=moLists(listnum).GetDim(0);
  int valscale=newvals.GetDim(1);
  int basisscale=basisvals.GetDim(1);
  int mocoeffscale=moCoeff.GetDim(1);
  for(int m=0; m < nmo_list; m++) {
    int mo=moLists(listnum)(m);
    int effmo=mo*mocoeffscale;
    int effm=m*valscale;
    int efff;
    for(int f=0; f< totbasis; f++) {
      //c=moCoeff(mo,f);
      //cout << " c " << c << endl;
      c=moCoeff.v[effmo+f];
      //cout << "c2 " << c << endl;
      efff=f*basisscale;
      for(int d=0; d< 5; d++) {
        //cout << "basisvals " << basisvals(f,d) << endl;
        //cout << "basisvals2 " << basisvals.v[efff+d] << endl;
        //
        //newvals(m,d)+=c*basisvals(f,d);
        //newvals.v[effm+d]+=c*basisvals(f,d);
        newvals.v[effm+d]+=c*basisvals.v[efff+d];
      }
    }
  }
}

//--------------------------------------------------------------------------
