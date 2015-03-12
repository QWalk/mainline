/*
 
Copyright (C) 2007 Zachary Helms
 with further modifications by Lucas K. Wagner

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

#include "Program_options.h"
#include "Lowdin_method.h"
#include "qmc_io.h"
#include "System.h"
#include "MatrixAlgebra.h"
#include "ulec.h"
/*!
Read the "words" from the method section in the input file
via doinput() parsing, and store section information in private
variables orbs, resolution, and minmax.
Set up MO_matrix and Sample_point objects for wavefunction from input
*/
void Lowdin_method::read(vector <string> words,
                       unsigned int & pos,
                       Program_options & options)
{

  
  sys=NULL;
  allocate(options.systemtext[0],  sys);

  if(! readvalue(words,pos=0,resolution,"RESOLUTION"))
    resolution=.1;
  if(!readvalue(words,pos=0, out_orbs, "OUT_ORB"))
    out_orbs=options.runid+".orb";

  vector <vector < string> > orbgroup_txt;
  pos=0;
  vector < string> dummy;
  while( readsection(words,pos,dummy,"ORB_GROUP")) { 
    orbgroup_txt.push_back(dummy);
  }

  vector <string> orbtext;
  if(!readsection(words, pos=0, orbtext, "ORBITALS")){
      error("Need ORBITALS");
  }
  allocate(orbtext, sys, mymomat);

  mywalker=NULL;
  sys->generateSample(mywalker);

  origin.Resize(3);
  origin=0;
  LatticeVec.Resize(3,3);
   
  if(!sys->getBounds(LatticeVec,origin))
    error("need getBounds");

}

//----------------------------------------------------------------------
int Lowdin_method::showinfo(ostream & os) { 
  os << "Lowdin localization method " << endl;
  return 1;
}

//----------------------------------------------------------------------

void Lowdin_method::run(Program_options & options, ostream & output) {
  int nmo=mymomat->getNmo();
  doublevar threshold=1e-10; //ignore overlaps smaller than this

  Array2 <doublevar> S;
  cout << "calc overlap " << endl;
  Array2 <doublevar>  Rtot(nmo,nmo);

  calculate_overlap(S);

  Array1 <doublevar> eigenvals(nmo);
  Array2 <doublevar> eigenvecs(nmo,nmo);
  EigenSystemSolverRealSymmetricMatrix(S,eigenvals,eigenvecs);
  for(int i=0; i < nmo; i++) { 
    for(int j=0; j< nmo; j++) { 
      if(fabs(S(i,j)) < threshold) S(i,j)=0.0;
      cout << setw(15) << S(i,j);
    }
    cout << endl;
  }
      
  Array2<doublevar> eigenvec_inverse(nmo,nmo);
  InvertMatrix(eigenvecs,eigenvec_inverse,nmo);
  cout << "Eigenvectors" << endl;
  for(int i=0; i < nmo; i++) { 
    for(int j=0; j< nmo; j++) { 
      cout << setw(15) << eigenvecs(i,j);
    }
    cout << endl;
  }
  cout << "Eigenvalues" << endl;
  for(int i=0; i< nmo; i++) cout << eigenvals(i) << " " ;
  cout << endl;

  for(int i=0; i< nmo; i++) { 
    if(eigenvals(i) < 0) 
      error("Negative eigenval: cannot perform Lowdin orthogonalization");
  }


  Rtot=0.;
  for(int i=0; i < nmo; i++) { 
    for(int j=0; j< nmo;j++) { 
      eigenvec_inverse(i,j)/=sqrt(eigenvals(i));
    }
  }
  MultiplyMatrices(eigenvecs,eigenvec_inverse,Rtot,nmo);
  Array1 <int> allorbs(nmo);
  for(int i=0; i< nmo; i++) allorbs(i)=i;
  ofstream testorb(out_orbs.c_str());
  mymomat->writeorb(testorb, Rtot,allorbs);
  testorb.close();
  


}

//----------------------------------------------------------------------
void Lowdin_method::calculate_overlap(Array2 <doublevar> & S) { 


  Array1 <Array1 <int> > tmp_list(1);
  int nmo=mymomat->getNmo();
  tmp_list(0).Resize(nmo);
  for(int i=0; i< nmo; i++) 
    tmp_list(0)(i)=i;
  mymomat->buildLists(tmp_list);
  int list=0;

  S.Resize(nmo,nmo);
  S=0.0;
  Array1 <doublevar> prop(3);
  Array1 <doublevar> n_lat(3);
  n_lat=0.0;
  for(int d=0;d < 3; d++) { 
    for(int d1=0; d1 < 3; d1++) 
      n_lat(d)+=LatticeVec(d,d1)*LatticeVec(d,d1);
    n_lat(d)=sqrt(n_lat(d));
    prop(d)=resolution/double(n_lat(d));
    //cout << "prop " << prop(d) <<  " res " << resolution << " norm " << n_lat(d) << endl;
  }

  cout << "calculating...." << endl;
  Array2 <doublevar> mymovals(nmo,1);
  Array1 <doublevar> r(3),x(3);
  int totpts=0;
  for(r(0)=0; r(0) < 1.0; r(0)+=prop(0)) {
    cout << "percent done: " << r(0)*100 << endl;
    for(r(1)=0; r(1) < 1.0; r(1)+=prop(1)) {
      for(r(2)=0; r(2) < 1.0; r(2)+=prop(2)) {
        totpts++;
        for(int d=0; d< 3; d++) { 
          x(d)=0.0;
          for(int d1=0; d1 < 3; d1++) { 
            x(d)+=r(d1)*LatticeVec(d1,d);
          }
          x(d)+=origin(d);
        }
        
        mywalker->setElectronPos(0,x);
        mymomat->updateVal(mywalker,0,list,mymovals);
        for(int i=0; i< nmo; i++) { 
          for(int j=0; j< nmo; j++) { 
            S(i,j)+=mymovals(i,0)*mymovals(j,0);
          }
        }
      }
    }
  }
  doublevar volume=Determinant(LatticeVec,3);
  doublevar norm=volume/totpts;
  for(int i=0; i< nmo; i++) { 
    for(int j=0; j< nmo; j++) { 
      S(i,j)*=norm;
    }
  }



}
//----------------------------------------------------------------------
