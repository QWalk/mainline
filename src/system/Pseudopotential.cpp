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


#include "Pseudopotential.h"
#include "qmc_io.h"
#include "Wavefunction.h"
#include "Sample_point.h"
#include <iomanip>
#include "ulec.h"
#include "Wavefunction_data.h"
#include "Basis_function.h"
using namespace std;
/*!
Some legendre polynomials..
*/
doublevar legendre(doublevar x, int n)
{
  switch(n)
  {
  case 0:
    return 1;
  case 1:
    return x;
  case 2:
    return .5*(3*x*x-1);
  case 3:
    return .5*(5*x*x*x - 3*x);
  case 4:
    return 0.125*(35*x*x*x*x -30*x*x +3);
  case 5:
    return 0.125*(63*x*x*x*x*x - 70*x*x*x + 15*x);
  default:
    error("Do not have legendre polynomial of order ", n);
    return 0; //shouldn't get here, but gets rid of a compiler message
  }
}

doublevar legendre_der(doublevar x, int n) {
  switch(n) {
  case 0: 
    return 0;
  case 1:
    return 1;
  case 2:
    return 3*x;
  case 3:
    return .5*(15*x*x-3);
  case 4:
    return 35*x*x*x-15*x;
  case 5:
    return .125*(315*x*x*x*x-210*x*x+15);
  default:
    error("Do not have legendre polynomial of order ", n);
    return 0;
  }
}
//----------------------------------------------------------------------


Pseudopotential::~Pseudopotential()  {
  for(int i=0; i< radial_basis.GetDim(0); i++)
    if(radial_basis(i)) delete radial_basis(i);
}






int Pseudopotential::initializeStatic(Wavefunction_data *wfdata,
                                      Sample_point * sample, Wavefunction * wf,
                                      Pseudo_buffer & output)
{

  int natoms=sample->ionSize();

  Array1 <doublevar> ionpos(3), oldpos(3), olddist(5), newpos(3);
  Array1 <doublevar> staticvals(wfdata->valSize());
  //Storage_container wfStore;
  if(! wfStore.isInitialized())
  {
    wfStore.initialize(sample, wf);
  }

  for(int at=0; at< natoms; at++)
  {
    if(numL(at) != 0)
    {
      sample->getIonPos(at, ionpos);
      for(int e=0; e < sample->electronSize(); e++)
      {
        sample->getElectronPos(e, oldpos);
        //cout << "electron position  " << oldpos(0)
        // << "   " << oldpos(1) << "   " << oldpos(2) << endl;
        sample->updateEIDist();//kind of inefficient..
        sample->getEIDist(e,at, olddist);
	
        //----------------------------------------
        //Start integral

        doublevar ratio=olddist(0)/cutoff(at);
        int shouldCalculate= ratio < 1.0;

        if(shouldCalculate)
        {
          wfStore.saveUpdate(sample, wf, e);

          for(int i=0; i< aip(at); i++)
          {
            sample->setElectronPos(e,oldpos);
            for(int d=0; d < 3; d++) 
              newpos(d)=integralpt(at,i,d)*olddist(0)-olddist(d+2);

            //for(int d=0; d < 3; d++)
            //{
            //  newpos(d)=ionpos(d)+integralpt(at, i, d)*olddist(0);
            //}
            sample->translateElectron(e, newpos);
            //cout << "wf " << endl;
            wf->storeParmIndVal(wfdata,sample, e, staticvals);
            for(int i=0; i< staticvals.GetDim(0); i++)
            {
              //fwrite(&staticvals(i), sizeof(doublevar), 1, output);
              output.push_value(staticvals(i));
            }
            sample->setElectronPos(e,oldpos);

          } //integral done


          wfStore.restoreUpdate(sample, wf, e);
        }  //if we should calculate


      }  //electron loop

    }  //if atom has any psp's
  }  //atom loop
  //delete wfStore;
  return 1;
}

//----------------------------------------------------------------------

void Pseudopotential::calcNonlocWithFile(Wavefunction_data * wfdata,
    Sample_point * sample,
    Wavefunction * wf,
    Array1 <doublevar> & totalv,
    Pseudo_buffer & input)
{
  int natoms=sample->ionSize();
  int nwf=wf->nfunc();

  assert(totalv.GetDim(0) >= nwf);
  assert(nelectrons == sample->electronSize());

  Array1 <doublevar> ionpos(3), oldpos(3), newpos(3);
  Array1 <doublevar> newdist(5), olddist(5);
  Wf_return val(nwf,2);
  Array1 <doublevar> staticval(wfdata->valSize());

  totalv=0;

  doublevar accum_local=0;
  doublevar accum_nonlocal=0;

  wf->updateVal(wfdata, sample);
  if( ! wfStore.isInitialized() ) {
    wfStore.initialize(sample, wf);
  }

  for(int at=0; at< natoms; at++)
  {
    if(numL(at) != 0)
    {
      sample->getIonPos(at, ionpos);
      for(int e=0; e < sample->electronSize(); e++)
      {

        sample->getElectronPos(e, oldpos);
        sample->updateEIDist();//kind of inefficient..
        sample->getEIDist(e,at, olddist);

        Array1 <doublevar>  nonlocal(nwf);
        nonlocal=0;

        Array1 <doublevar> v_l(numL(at));
        //cout << "getRadial " << endl;
        getRadial(at, sample, olddist, v_l);


        //----------------------------------------
        //For this one, we decide based on a cutoff, not randomly
        if(olddist(0) < cutoff(at))  {
          wfStore.saveUpdate(sample, wf, e);
          Wf_return oldWfVal(nwf,2);
          wf->getVal(wfdata, e,oldWfVal);


          Array2 <doublevar> integralpts(nwf, aip(at));
          Array1 <doublevar> rDotR(aip(at));
          
          //Gather the points for the integral
          for(int i=0; i< aip(at); i++) {
            
            
            sample->setElectronPos(e, oldpos);
            doublevar base_sign=sample->overallSign();
            //Make sure to move the electron relative to the nearest neighbor
            //in a periodic calculation(so subtract the distance rather than
            //adding to the ionic position).  This actually only matters 
            //when we're doing non-zero k-points.
            for(int d=0; d < 3; d++) 
              //newpos(d)=ionpos(d)+integralpt(at, i, d)*olddist(0)-oldpos(d);
              newpos(d)=integralpt(at,i,d)*olddist(0)-olddist(d+2);
                        

            sample->translateElectron(e, newpos);
            if(!sample->getEIDist_temp(e,at,newdist)) {
              sample->updateEIDist();
              sample->getEIDist(e,at,newdist);
            }
            

            rDotR(i)=0;
            for(int d=0; d < 3; d++) {
              rDotR(i)+=newdist(d+2)*olddist(d+2);
            }
            rDotR(i)/=(newdist(0)*olddist(0));  
            double new_sign=sample->overallSign();
            for(int j=0; j< staticval.GetDim(0); j++) {
              //fread(&staticval(j), sizeof(doublevar), 1, input);
              staticval(j)=input.next_value();
            }
            
            wf->getParmDepVal(wfdata, sample, e, staticval, val);
            for(int w=0; w< nwf; w++) {
              integralpts(w,i)=val.sign(w)*oldWfVal.sign(w)*new_sign*base_sign
                               *exp(val.amp(w,0)-oldWfVal.amp(w,0))
                               *integralweight(at, i);
            }
            //sample->setElectronPos(e, oldpos);
          }
          
          //Now do the integral
          for(int w=0; w< nwf; w++) {
            for(int i=0; i< aip(at); i++) {
              doublevar tempsum=0;
              for(int l=0; l< numL(at)-1; l++) {
                tempsum+=(2*l+1)*v_l(l)*legendre(rDotR(i), l);
              }
              nonlocal(w)+=tempsum*integralpts(w,i);
            }
          }
          wfStore.restoreUpdate(sample, wf, e);
        }


        //----------------------------------------------
        //now do the local part
        doublevar vLocal=0;
        int localL=numL(at)-1; //The l-value of the local part is
                               //the last part.

        vLocal=v_l(localL);

        accum_local+=vLocal;
        accum_nonlocal+=nonlocal(0);

        //cout << "vLocal  " << accum_local
        // << "    nonlocal   " << accum_nonlocal
        // << endl;
        for(int w=0; w< nwf; w++)
        {
          totalv(w)+=vLocal+nonlocal(w);
        }

        //cout << "totalv " << totalv << endl;

      }  //electron loop
    }  //if atom has any psp's
  }  //atom loop

  //cout << "psp: local part " << accum_local
  // << "  nonlocal part " << accum_nonlocal << endl;

}

//----------------------------------------------------------------------


/*!
Evaluates the radial part of the pseudopotential for a given
atom and l-value.

 */

void Pseudopotential::getRadial(int at, Sample_point * sample,
                                Array1 <doublevar> & r, Array1 <doublevar> & v_l) {
  assert(radial_basis(at) != NULL);
  radial_basis(at)->calcVal(r, v_l);
  if(addzeff(at)) {
    //for(int l=0; l < numL(at); l++) {
    int l=numL(at)-1;
    doublevar cutoff_rad=radial_basis(at)->cutoff(l);
    if(r(0) < cutoff_rad) {
      v_l(l) += sample->getIonCharge(at)/r(0);
    }
    //}
  }
}


void Pseudopotential::getRadial(int at, Sample_point * sample,
                                Array1 <doublevar> & r, 
                                Array2 <doublevar> & v_l){
  assert(radial_basis(at) != NULL);
  radial_basis(at)->calcLap(r, v_l);
  if(addzeff(at)) {
    int l=numL(at)-1;
    doublevar cutoff_rad=radial_basis(at)->cutoff(l);
    if(r(0) < cutoff_rad) {
      v_l(l,0) += sample->getIonCharge(at)/r(0);
      for(int d=0; d< 3; d++) {
        v_l(l,d+1) -= sample->getIonCharge(at)*r(d+2)/(r(0)*r(1));
      }
    }
  }
}



//----------------------------------------------------------------------


int Pseudopotential::nTest() {
  int tot=0;
  int natoms=numL.GetDim(0);
  for(int at=0; at < natoms; at++) {
    if(numL(at) != 0) {
      tot+=nelectrons;
    }
  }
  return tot;
}


//----------------------------------------------------------------------

void Pseudopotential::calcNonloc(Wavefunction_data * wfdata,
                                 Sample_point * sample, Wavefunction * wf,
                                 Force_fitter & fitter,
                                 Array1 <doublevar> & totalv, 
                                 Array2 <doublevar> & deriv) {
  int tot=nTest();
  Array1 <doublevar> test(tot);
  for(int i=0; i< tot; i++) {
    test(i)=rng.ulec();
  }
  calcNonlocWithTest(wfdata, sample, wf,fitter, test, totalv, deriv);
}
//----------------------------------------------------------------------

void Pseudopotential::calcNonlocWithTest(Wavefunction_data * wfdata,
                                         Sample_point * sample,
                                         Wavefunction * wf,
                                         Force_fitter & fitter,
                                         const Array1 <doublevar> & accept_var,
                                         Array1 <doublevar> & totalv, 
                                         Array2 <doublevar> & deriv)
{

  int natoms=sample->ionSize();
  int nwf=wf->nfunc();

  assert(accept_var.GetDim(0) >= nTest());
  assert(totalv.GetDim(0) >= nwf);
  assert(nelectrons == sample->electronSize());

  Array1 <doublevar> ionpos(3), oldpos(3), newpos(3);
  Array1 <doublevar> newdist(5), olddist(5);
  Wf_return val(nwf,2);

  totalv=0;
  doublevar accum_local=0;
  doublevar accum_nonlocal=0;

  wf->updateVal(wfdata, sample);
  wfStore.initialize(sample, wf);
  
  int accept_counter=0;
  deriv.Resize(natoms, 3);
  deriv=0;
  
  for(int at=0; at< natoms; at++){
    if(numL(at) != 0) {
      sample->getIonPos(at, ionpos);

      for(int e=0; e < sample->electronSize(); e++)  {
        sample->getElectronPos(e, oldpos);
        
        Array1 <doublevar> der_this_e(3);
        der_this_e=0;
        //note: this updateEIDist might become inefficient..
        //depends on how often we're rejecting/how much it costs
        //to do it.  If needed, we should add an interface to 
        //Sample_point
        sample->updateEIDist();
        sample->getEIDist(e,at, olddist);
        Array1 <doublevar>  nonlocal(nwf);
        nonlocal=0;

        Array2 <doublevar> v_l(numL(at),5);
        getRadial(at, sample, olddist, v_l);

        //----------------------------------------
        //Start integral

        int accept;
        if(deterministic) {
          accept= olddist(0) < cutoff(at);
        }
        else {
          doublevar strength=0;
          const doublevar calculate_threshold=10;
  
          for(int l=0; l<numL(at)-1; l++)
            strength+=calculate_threshold*(2*l+1)*fabs(v_l(l,0));
   
          strength=min((doublevar) 1.0, strength);
  
          doublevar rand=accept_var(accept_counter++);
         //cout << at <<"  random number  " << rand
        //   << "  p_eval  " << strength  << endl;
          for(int l=0; l<numL(at)-1; l++)
            for(int d=0; d< 5; d++) 
              v_l(l,d)/=strength;
          accept=strength>rand;
        }
      

        if(accept)  {
          wfStore.saveUpdate(sample, wf, e);
          Wf_return  oldWfVal(nwf,2);
          wf->getVal(wfdata, e,oldWfVal);

          Array2 <doublevar> integralpts(nwf, aip(at));
          Array1 <doublevar> rDotR(aip(at));
          
          for(int i=0; i< aip(at); i++) {
            sample->setElectronPos(e, oldpos);
            doublevar base_sign=sample->overallSign();
            //Make sure to move the electron relative to the nearest neighbor
            //in a periodic calculation(so subtract the distance rather than
            //adding to the ionic position).  This actually only matters 
            //when we're doing non-zero k-points.
            for(int d=0; d < 3; d++) 
	      newpos(d)=integralpt(at,i,d)*olddist(0)-olddist(d+2);
	    
            //cout << "translation " << newpos(0) << "   " 
             //   << newpos(1) << "   " << newpos(2) << endl;
            
            sample->translateElectron(e, newpos);
            sample->updateEIDist();
            sample->getEIDist(e,at,newdist);

            rDotR(i)=0;
            for(int d=0; d < 3; d++)
              rDotR(i)+=newdist(d+2)*olddist(d+2);
            double new_sign=sample->overallSign();
            
            rDotR(i)/=(newdist(0)*olddist(0));  //divide by the magnitudes
            wf->updateVal(wfdata, sample);
            wf->getVal(wfdata, e, val); 
            //cout << "signs " << base_sign << "  " << new_sign << endl;;
            for(int w=0; w< nwf; w++) {
              integralpts(w,i)=val.sign(w)*oldWfVal.sign(w) *base_sign*new_sign
                               *exp(val.amp(w,0)-oldWfVal.amp(w,0))
                  *integralweight(at, i);
            }
            sample->setElectronPos(e, oldpos);
          }

          //--------------------

          for(int w=0; w< nwf; w++)  {
            for(int i=0; i< aip(at); i++)   {
              doublevar tempsum=0;
              for(int l=0; l< numL(at)-1; l++) {
                tempsum+=(2*l+1)*v_l(l,0)*legendre(rDotR(i), l);
              }
              nonlocal(w)+=tempsum*integralpts(w,i);
              //wavefunc in (modified) f90
              //cout << "wavefunction "
              //     << integralpts(w,i)/integralweight(at,i) << endl;
            }
          }

          //putting in the part for s-nonlocality only..p and higher to follow.
          //if(olddist(0) > fitter.cutoff()) {
         
          for(int d=0; d< 3; d++) {
            for(int i=0; i< aip(at); i++) {
              doublevar tempsum=0;
              for(int l=0; l< numL(at)-1; l++) 
                tempsum+=(2*l+1)*v_l(l,d+1)*legendre(rDotR(i),l);
              der_this_e(d)+= tempsum*integralpts(0,i);
            }
          }
          
          //}

          wfStore.restoreUpdate(sample, wf, e);
        }

        //----------------------------------------------
        //now do the local part
        doublevar vLocal=0;
        int localL=numL(at)-1; //The l-value of the local part is
                               //the last part.

        vLocal=v_l(localL,0);
        accum_local+=vLocal;
        accum_nonlocal+=nonlocal(0);

        //cout << "vLocal  " << accum_local
        // << "    nonlocal   " << accum_nonlocal
        // << endl;
        for(int w=0; w< nwf; w++) {
          totalv(w)+=vLocal+nonlocal(w);
        }


        for(int d=0; d< 3; d++) 
          der_this_e(d)+=-v_l(localL,d+1);

        Array1 <doublevar> fitted(3);
        fitter.fit_force(olddist,der_this_e, fitted);
        for(int d=0; d< 3; d++) 
          deriv(at,d)+=fitted(d);

      }  //electron loop

    }  //if atom has any psp's

  }  //atom loop

  //cout << "psp: local part " << accum_local
  // << "  nonlocal part " << accum_nonlocal << endl;

}

//------------------------------------------------------------------------


void Pseudopotential::getLocalDerivative(Sample_point * sample, 
                                         Force_fitter & fitter,
                                         Array2 <doublevar> & der) {
  
  int natoms=sample->ionSize();
  der.Resize(natoms, 3);
  der=0;
  Array1 <doublevar> olddist(5);
  Array1 <doublevar> tmpder(3), tmpder_fit(3);

  for(int at=0; at< natoms; at++){
    if(numL(at) != 0){
      Array2 <doublevar> v_l(numL(at), 5);

      for(int e=0; e < sample->electronSize(); e++) {
        sample->getEIDist(e,at, olddist);
        getRadial(at, sample, olddist, v_l);
        int lastL=numL(at)-1;
        for(int d=0; d< 3; d++) 
          //der(at,d)-=v_l(lastL, d+1);
          tmpder(d)=-v_l(lastL,d+1);
        fitter.fit_force(olddist,tmpder, tmpder_fit);
        for(int d=0; d< 3; d++) 
          der(at,d)+=tmpder_fit(d);
      }
    }
  }

}

//----------------------------------------------------------------------
void Pseudopotential::randomize() {
  Array1 <doublevar> x(3), y(3), z(3);
  generate_random_rotation(x,y,z);
  rotateQuadrature(x,y,z);
}


//----------------------------------------------------------------------
void Pseudopotential::rotateQuadrature(Array1 <doublevar> & x,
                                       Array1 <doublevar> & y,
                                       Array1 <doublevar> & z)
{


  int natoms=aip.GetDim(0);

  //cout << "x1, x2, x3" << x1 << "   " << x2 << "   " << x3 << endl;
  //cout << "y1, y2, y3" << y1 << "   " << y2 << "   " << y3 << endl;
  //cout << "z1, z2, z3" << z1 << "   " << z2 << "   " << z3 << endl;

  for(int at=0; at<natoms; at++)
  {
    //cout << "-----------atom   " << at << endl;
    for(int i=0; i< aip(at); i++)
    {
      //cout << "quadrature points before:\n"
      //   << integralpt(at, i, 0) << "   "
      //   <<integralpt(at, i, 1) << "    "
      //   <<integralpt(at, i, 2) << "    \n";
      integralpt(at,i,0)=integralpt_orig(at,i,0)*x(0)
                         +integralpt_orig(at,i,1)*y(0)
                         +integralpt_orig(at,i,2)*z(0);
      integralpt(at,i,1)=integralpt_orig(at,i,0)*x(1)
                         +integralpt_orig(at,i,1)*y(1)
                         +integralpt_orig(at,i,2)*z(1);
      integralpt(at,i,2)=integralpt_orig(at,i,0)*x(2)
                         +integralpt_orig(at,i,1)*y(2)
                         +integralpt_orig(at,i,2)*z(2);
      //cout << "quadrature points after:\n"
      //   << integralpt(at, i, 0) << "   "
      //	   <<integralpt(at, i, 1) << "    "
      //   <<integralpt(at, i, 2) << "    \n";
    }
  }
}



//------------------------------------------------------------------------

#include "System.h"

/*!
*/
void Pseudopotential::read(vector <vector <string> > & pseudotext,
                           System * sys)
{

  sys->getAtomicLabels(atomnames);
  int natoms=atomnames.size();
  nelectrons=sys->nelectrons(0)+sys->nelectrons(1);
  //Get the atomic integration points
  aip.Resize(natoms);
  aip=6; //default value for atomic integration points
  integralpt.Resize(natoms, maxaip, 3);
  integralpt_orig.Resize(natoms, maxaip, 3);
  integralweight.Resize(natoms, maxaip);
  addzeff.Resize(natoms);
  addzeff=false;
  Array1 <int> atom_has_psp(natoms);
  atom_has_psp=0;
  int maxL=0;

  radial_basis.Resize(natoms);
  radial_basis=NULL;
  vector < vector < string > > basistxt;
  for(unsigned int i=0; i< pseudotext.size(); i++) {
    unsigned int pos=0;
    vector <string> basistmp;
    if(!readsection(pseudotext[i], pos=0, basistmp, "BASIS"))
      error("Need Basis section in pseudopotential for ", pseudotext[i][0]);
    vector <string> oldqmctxt;
    if(readsection(pseudotext[i], pos=0, oldqmctxt, "OLDQMC")) {
      error("OLDQMC is out of date.  Please enclose it in:\nBASIS { \n",
            pseudotext[i][0], "\n RGAUSSIAN \n  ....   \n}");
    }    
    basistxt.push_back(basistmp);
  }

  for(unsigned int i=0; i< pseudotext.size(); i++ )
  {
    for(int at=0; at<natoms; at++)
    {
      unsigned int pos=0;

      if( pseudotext[i][0] == atomnames[at] )
      {
        if(atom_has_psp(at))
          error("There are two PSEUDO sections for ", atomnames[at]);
        atom_has_psp(at)=1;
        string type;
        //vector <string> oldqmctxt;
        //vector <string> basistxt;
        //if(readsection(pseudotext[i], pos=0, basistxt, "BASIS")) {
        vector <string> tmpbasis(basistxt[i]);
          allocate(tmpbasis, radial_basis(at));
          int nlval=radial_basis(at)->nfunc();
          if(maxL < nlval) maxL=nlval;

          //}
          //else if(readsection(pseudotext[i], pos=0, oldqmctxt, "OLDQMC")) {
          //error("OLDQMC is out of date.  Please enclose it in:\nBASIS { \n",
          //     atomnames[at], "\n RGAUSSIAN \n  ....   \n}");
          // }
          // else {
          // error("Need BASIS section in pseudopotential for ", atomnames[at]);
          //}

        pos=0;
        readvalue(pseudotext[i], pos, aip(at), "AIP");
        if(haskeyword(pseudotext[i], pos=0, "ADD_ZEFF")) {
          addzeff(at)=true;
        }

      }
    }
  }



  for(int at=0; at<natoms; at++)
  {
    int aiptemp=aip(at);
    Array1 <doublevar> xpt(maxaip);
    Array1 <doublevar> ypt(maxaip);
    Array1 <doublevar> zpt(maxaip);
    Array1 <doublevar> weight(maxaip);


    gesqua(aiptemp, xpt, ypt, zpt, weight);

    aip(at)=aiptemp;

    for(int i=0; i< aip(at); i++)
    {
      integralpt(at,i,0)=integralpt_orig(at, i, 0)=xpt(i);
      integralpt(at,i,1)=integralpt_orig(at, i, 1)=ypt(i);
      integralpt(at,i,2)=integralpt_orig(at, i, 2)=zpt(i);
      integralweight(at, i)=weight(i);
    }

  }


  //allocate the storage

  numL.Resize(natoms);

  numL=0;
  for(int at=0; at < natoms; at++ )
  {
    if(radial_basis(at) != NULL)
      numL(at)=radial_basis(at)->nfunc();
  }


  //find the cutoff radius for a static calculation

  Sample_point * tempsample=NULL;
  sys->generateSample(tempsample);
  cutoff.Resize(natoms);
  const doublevar cutoff_threshold=1e-5;
  const doublevar cutoff_max=20.0;
  cutoff=0.0;
  const doublevar cutoff_interval=.05;
  for(int at=0; at< natoms; at++)
  {
    if(numL(at) > 0) {
    Array1 <doublevar> cutoffL(numL(at));
    Array1 <doublevar> v_l(numL(at));
    Array1 <doublevar> tempr(5);
    Array1 <int> foundcutoff(numL(at));
    foundcutoff=0;
    tempr=0;
    cutoffL=cutoff_max;
    for(doublevar r=cutoff_max; r> 0; r-=cutoff_interval)
    {
      //going along the z axis; watch this if we ever try non-spherical
      //potentials

      tempr(0)=r; tempr(1)=r*r; tempr(4)=r;
      getRadial(at, tempsample, tempr, v_l);
      for(int l=0; l< numL(at)-1; l++)
      {
        if(fabs(v_l(l)) > cutoff_threshold && !foundcutoff(l))
        {
          //It seems to me that it should be r+cutoff_interval,
          //but for compatability with the f90 code, we'll put it
          //at r-cutoff_interval.  Shouldn't be too critical, anyway.
          cutoffL(l)=r-cutoff_interval;
          foundcutoff(l)=1;
        }
      }
    }

    for(int l=0; l< numL(at)-1; l++)
    {
      cutoff(at)=max(cutoff(at), cutoffL(l));
    }
    }
  }
  delete tempsample;

}


//----------------------------------------------------------------------

int Pseudopotential::showinfo(ostream & os)
{
  os << "Pseudopotential " << endl;
  string indent="  ";
  vector <string> uniquenames;
  int natoms=aip.GetDim(0);
  for(int at=0; at < natoms; at++) {
    int unique=1;
    for(unsigned int i=0; i< uniquenames.size(); i++) {
      if(uniquenames[i]==atomnames[at]) {
        unique=0;
        break;
      }
    }
    if(unique) {
      uniquenames.push_back(atomnames[at]);
      
      os << "atom " << atomnames[at] << endl;
      if(numL(at)==0) { os << "No pseudopotential" << endl; }
      else {
        os << "Integration points " << aip(at) << endl;
//os << setw(10) << "x" 
//   << setw(10) << "y" 
//   << setw(10) << "z"
//   << setw(10) << "weight" << endl;
//	for(int i=0; i< aip(at); i++)
//	  {
//	    os << setw(10) << integralpt(at,i,0)
//	       << setw(10) << integralpt(at,i,1)
//	       << setw(10) << integralpt(at,i,2)
//	       << setw(10)  << integralweight(at, i) << endl;
//	  }	
        radial_basis(at)->showinfo(indent, os);
      }
    }
  }

  return 1;
}

//------------------------------------------------------------------------
