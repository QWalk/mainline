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
//----------------------------------------------------------------------

#include "Space_warper.h"

#include "MatrixAlgebra.h"
#include "Sample_point.h"
#include "qmc_io.h"

Space_warper::Space_warper() { 
  weight_basis=NULL; 
    vector <string> basistxt;
    basistxt.push_back("NULL");
    basistxt.push_back("POLYPADE");
    basistxt.push_back("RCUT");
    basistxt.push_back("2.0");
    basistxt.push_back("BETA0");
    basistxt.push_back("-.9");
    basistxt.push_back("NFUNC");
    basistxt.push_back("1");
    allocate(basistxt, weight_basis);
    warp_on=1;
    ex=1;
}

void Space_warper::read(vector <string> & words) {
  unsigned int pos=0;
  readvalue(words, pos=0, ex, "EXPONENT");
}


//----------------------------------------------------------------------
void transformation_matrix (const Array2 <doublevar> &ellePrime, 
			    const Array2 <doublevar> &elle, 
			    Array2 <doublevar> &ratio) {
  
  assert(elle.GetDim(0)==elle.GetDim(1)==3);
  assert(ellePrime.GetDim(0)==ellePrime.GetDim(1)==3);
  ratio.Resize(3,3);
  Array2 <doublevar> elleInverse(3,3,0.0);

  InvertMatrix(elle, elleInverse, 3);
  MultiplyMatrices(ellePrime, elleInverse, ratio, 3);

}

//----------------------------------------------------------------------
int Space_warper::space_warp (Sample_point * refsample, Sample_point * sample,
                              int e, Array1 <doublevar> & R_old,
                              Array1 <doublevar> & R_new,
                              doublevar & jacobian ) {

  int natoms=refsample->ionSize();
  int moved=0;
  refsample->updateEIDist();
  Array1 <doublevar> dist(5);

  Array2 <doublevar>  jacob_matrix(3,3);
  jacob_matrix=0;
  assert(R_old.GetDim(0)==3);
  R_new.Resize(3);
  R_new=R_old;

  
  if(!warp_on)
  { jacobian=1; return 0; }

  //Here we stretch or shrink the cell appropriately to match the 
  //different lattice parameters.  This is a guaranteed way to 
  //make sure the two spaces actually have a 1-1 relationship in a PBC calculation.
  Array2 <doublevar> reflatvec;
  Array1 <doublevar> reforigin;
  if(refsample->getBounds(reflatvec, reforigin)) {
    Array2 <doublevar> latvec;
    Array1 <doublevar> origin;
    sample->getBounds(latvec, origin);
    assert(reflatvec.GetDim(0)==3 && reflatvec.GetDim(1)==3 && reforigin.GetDim(0)==3);
    assert(latvec.GetDim(0)==3 && latvec.GetDim(1)==3 && origin.GetDim(0)==3);

    transformation_matrix(latvec, reflatvec, jacob_matrix);
    R_new=0;
    //cout << " jacobian " << endl;
    for(int i=0; i< 3; i++) {
      for(int j=0; j< 3; j++) {
	R_new(i)+=jacob_matrix(i,j)*R_old(j);
	//cout << jacob_matrix(i,j) << "   ";
      }
      //cout << endl;
    }

    //cout << "rnew   rold \n";
    //for(int i=0; i< 3; i++) cout << R_new(i) <<"   " << R_old(i) << endl;
    //cout << endl;
    jacobian=Determinant(jacob_matrix,3);
    return 1;
  }  
  
  
  //If it is not a periodic calculation, then space-warp using
  //either a general basis function or the 1/r warping suggested
  //by Filippi and Umrigar(using 1/r by default because the full
  //estimator works better with it)
  
  //doublevar cutoff=weight_basis->cutoff(0);
  doublevar cutoff=1e99;
  
  //Array2 <doublevar> basisval(1,5);

  Array1 <doublevar> displace(3, 0.0);
  //if(0) {
  doublevar norm=0;
  Array1 <doublevar> der_norm(3,0.0);
  Array1 <doublevar> weights(natoms,0.0);
  Array2 <doublevar> weight_der(natoms, 3,0.0);
  Array2 <doublevar> deltar(natoms, 3,0.0);

  Array1 <doublevar> oldrefpos(3);
  refsample->getElectronPos(e,oldrefpos);

  
  
  
  if(natoms >0) {
    refsample->setElectronPosNoNotify(e,R_old);
    refsample->updateEIDist();
  }

  for(int at=0; at < natoms; at++) {
    
    refsample->getEIDist(e,at, dist);
    //cout << "atom " << at << "dist " << dist(0) << endl;
    if(dist(0) < cutoff) {
      //if(moved==1) cout << "Warning--warped twice " << endl;
    
    moved=1;
    Array1 <doublevar> refionpos(3);
    refsample->getIonPos(at, refionpos);
    Array1 <doublevar> newionpos(3);
    sample->getIonPos(at, newionpos);
    
    //weight_basis->calcLap(dist, basisval);
    doublevar z=1.0;//refsample->getIonCharge(at);
    
   
    
    //doublevar func=z/(dist(1)*dist(1));
    //doublevar func=z/(dist(1));
    doublevar func=z/(pow(dist(0),ex));
    //doublevar func=basisval(0,0);
    norm+=func;
    
    for(int d=0; d< 3; d++) {
      displace(d) += func*(newionpos(d)-refionpos(d));
    }
    
    weights(at)=func;
    for(int d=0; d< 3; d++) {
      deltar(at,d)=newionpos(d)-refionpos(d);
      //weight_der(at,d)= -4*z*dist(d+2)/(dist(1)*dist(1)*dist(1));
      weight_der(at,d)= -ex*z*dist(d+2)/pow(dist(0),ex+2);
      //weight_der(at,d)=basisval(0,d+1);
    }
    for(int d=0; d< 3; d++) der_norm(d)+=weight_der(at,d);
    }
  }
  

  if(natoms >0) {
    refsample->setElectronPosNoNotify(e,oldrefpos);
    refsample->updateEIDist();
  }

  //cout << "norm " << norm << endl;
  if(moved && norm > 1e-14) {
    for(int d=0; d< 3; d++) {
      R_new(d)+=displace(d)/norm;
    }

    jacob_matrix=0;
    for(int i=0; i< 3; i++) 
      for(int j=0; j< 3; j++) 
        for(int at=0; at < natoms; at++) {
          jacob_matrix(i,j)+=deltar(at,i)*(weight_der(at,j)/norm
                                           -weights(at)*der_norm(j)/(norm*norm));
          //cout << " at " << at << "der " << weight_der(at,j) << " weights " 
          //     << weights(at) << "  der_norm " << der_norm(j)
          //     << endl;
          }
  }
  //}


  for(int i=0; i< 3; i++) jacob_matrix(i,i)+=1;
  
  //cout << "jacobian matrix " << endl;
  //for(int i=0; i< 3; i++) {
  //  for(int j=0; j< 3; j++) 
  //    cout << jacob_matrix(i,j) << "   ";
  //  cout << endl;
  //}
  jacobian=Determinant(jacob_matrix, 3);
  //cout << "jacobian " << jacobian << endl;
  return moved;
}

//----------------------------------------------------------------------
