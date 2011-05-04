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
#include "MO_matrix_Ccutoff.h"
#include "Sample_point.h"
#include "qmc_io.h"



void MO_matrix_Ccutoff::init() {

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
  
  //--

  nbasis.Resize(nmo);
  mofill.Resize(nmo, totbasis);
  moCoeff2.Resize(nmo, totbasis);

  //---


  ifstream ORB(orbfile.c_str());

  if(!ORB)
  {
    error("couldn't find orb file ", orbfile);
  }

  Array3 <int> coeffmat;
  Array1 <dcomplex> coeff;
  
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
    
    dcomplex kptfac=exp(dot*pi*dcomplex(0.0,1.0));
    //cout << "kptfac " << kptfac << "  displacement " 
    //    << centers.centers_displacement(ion,0) << "   "
    //    << endl;            
    
    for(int n=0; n< centers.nbasis(ion); n++) {
      
      int fnum=centers.basis(ion,n);
      int imax=basis(fnum)->nfunc();

      for(int i=0; i<imax; i++) { //sum over the symmetries
        for(int mo=0; mo<nmo; mo++) {      //and the MO's
          if(coeffmat(mo,ion, f) == -1) {
            cout << "missing MO pointer: mo# " << mo << " ion # " << ion
            << " function on ion: " << f << endl;
            error("In the orb file, there is a missing pointer. It might "
                  "be a badly structured file.");
          }
          dcomplex temp=coeff(coeffmat(mo,ion,f));
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

  /*
  single_write(cout, "total MO's ", nmo, "\n");
  for(int basis=0; basis < totfunc; basis++) {
    single_write(cout, "basismo(", basis, ")  ");
    single_write(cout, basismo(basis), "\n");
  }
  */


  //single_write(cout,  "total functions ", totfunc, "\n") ;
  //for(int mo=0; mo < nmo; mo++) {
  //  single_write(cout,"nbasis(",mo, ")   ");
  //  single_write(cout, nbasis(mo),"\n");
  //}

}

//------------------------------------------------------------------------------

void MO_matrix_Ccutoff::writeorb(ostream & os, Array2 <doublevar> & rotation,
				 Array1 <int>  &moList) {


  int nmo_write=moList.GetDim(0);
  assert(rotation.GetDim(0)==nmo_write);
  assert(rotation.GetDim(1)==nmo_write);

  os.precision(15);
  int counter=0;
  for(int m=0; m < nmo_write; m++) {
    int mo=moList(m);
    for(int ion=0; ion<centers.equiv_centers.GetDim(0); ion++)
    {
      int f=0;

      for(int n=0; n< centers.nbasis(centers.equiv_centers(ion,0)); n++)
      {
        int fnum=centers.basis(ion,n);
        int imax=basis(fnum)->nfunc();

        for(int i=0; i<imax; i++)
        {
          os << mo+1 << "  "   << f+1 << "   " << ion+1 << "   " << counter+1 << endl;
          f++;  //keep a total of functions on center
          counter++;
        } //i
      } //n
    }  //ion
  }
  os << "COEFFICIENTS\n";
  ifstream orbin(orbfile.c_str());
  rotate_orb(orbin, os, rotation, moList,
	     int((totbasis/centers.size())*centers.equiv_centers.GetDim(0)));
  orbin.close();


}
//---------------------------------------------------------------------

void MO_matrix_Ccutoff::buildLists(Array1 < Array1 <int> > & occupations)
{
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

int MO_matrix_Ccutoff::showinfo(ostream & os)
{
  os << "Complex cutoff MO " << endl;
  os << "Number of molecular orbitals: " << nmo << endl;
  string indent="  ";
  os << "Basis functions: \n";
  for(int i=0; i< basis.GetDim(0); i++)
  {
    basis(i)->showinfo(indent, os);
  }
  return 1;
}

int MO_matrix_Ccutoff::writeinput(string & indent, ostream & os)
{
  os << indent << "CCUTOFF_MO" << endl;
  os << indent << "NMO " << nmo << endl;
  os << indent << "ORBFILE " << orbfile << endl;
  if(oldsofile!="") 
    os << indent << "OLDSOFILE " << oldsofile << endl;
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

//----------------------------------------------------------------------

/*!
Basically just a copy of MO_matrix::read. Certainly not the most elegant
way to do it but the present form of complex support here is quite
constraining - the 1D orbital needs different read(), unfortunately.
Although there is a slight difference even compared to MO_matrix_Cbasfunc,
since here we have real basis functions and there are complex basis
functions.
*/
/*
void MO_matrix_Ccutoff::read(vector <string> & words,
			     unsigned int & startpos,
			     System * sys)
{


  unsigned int pos=startpos;

  if(!readvalue(words, pos, nmo, "NMO"))
  {
    error("Need NMO in molecular orbital section");
  }

  if(nmo > 40000) 
    error("You have entered more than 40,000 for NMO. This seems a bit big; we most likely"
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
*/

//------------------------------------------------------------------------

void MO_matrix_Ccutoff::updateVal(Sample_point * sample, int e,
				  int listnum,
				  //!<which list to use
				  //const Array1 <int> & occupation,
				  //!<A list of the MO's to evaluate
				  Array2 <dcomplex> & newvals
				  //!< The return: in form (MO, val)
				  )
{
  //cout << "updateVal " << listnum << endl;
  int centermax=centers.size();
  //int momax=occupation.GetDim(0);
  static Array1 <doublevar> R(5);
  //static Array1 <doublevar> symmvals(totbasis);
  static Array1 <doublevar> symmvals_temp(maxbasis);

  //Make references for easier access to the list variables.
  Array1 <int> & basismotmp(basismo_list(listnum));
  Array2 <int> & basisfilltmp(basisfill_list(listnum));
  Array2 <dcomplex> & moCoefftmp(moCoeff_list(listnum));
  assert(newvals.GetDim(1) >= 1);

  newvals=0;
  Basis_function * tempbasis;

  dcomplex c;
  int mo=0;
  int scalebasis=basisfill_list(listnum).GetDim(1);
  int totfunc=0;
  int b; //basis
  
  centers.updateDistance(e, sample);
  //int retscale=newvals.GetDim(1);
  for(int ion=0; ion < centermax; ion++)
  {
    //sample->getECDist(e, ion, R);
    centers.getDistance(e,ion,R);
    for(int n=0; n< centers.nbasis(ion); n++)
    {
      b=centers.basis(ion,n);
      tempbasis=basis(b);
      if(obj_cutoff(b) > R(0)) {
      tempbasis->calcVal(R, symmvals_temp);
      int imax=nfunctions(b);
      for(int i=0; i< imax; i++)
      {
        //cout << "i " << i << endl;
        int reducedbasis=scalebasis*totfunc;
        if(R(0) < cutoff(totfunc))
        {
          //  cout << "basmo " << basismotmp.v[totfunc] <<  endl;
          for(int basmo=0; basmo < basismotmp.v[totfunc]; basmo++)
          {
            //mo=basisfill(totfunc, basmo);
            //c=moCoeff(totfunc, basmo);
            //cout << "basisfill reducedbasis "<< reducedbasis
            // << "  basmo " << basmo << endl;
            mo=basisfilltmp.v[reducedbasis+basmo];
            //cout << "mocoeff (mo=" << mo <<  endl;
            //mo_counter(mo)++;

            c=moCoefftmp.v[reducedbasis+basmo];
            newvals(mo, 0)+=c*symmvals_temp(i);
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
  //  cout << "done updateVal " << endl;
}

//------------------------------------------------------------------------
inline void output_array(Array2 <doublevar> & arr) {
  for(int i=0; i< arr.GetDim(0); i++) {
    for(int j=0; j < arr.GetDim(1); j++) {
      cout << arr(i,j) << "  ";
    }
    cout << endl;
  }
}

inline void output_array(Array2 <dcomplex> & arr) {
  for(int i=0; i< arr.GetDim(0); i++) {
    for(int j=0; j < arr.GetDim(1); j++) {
      cout << arr(i,j) << "  ";
    }
    cout << endl;
  }
}

/*!
*/
void MO_matrix_Ccutoff::updateLap(Sample_point * sample, int e,
				  int listnum,
				  //const Array1 <int> & occupation,
				  //!<A list of the MO's to evaluate
				  Array2 <dcomplex> & newvals
				  //!< The return: in form (MO, [val, grad, lap])
				  )
{

  //cout << "updateLap" << endl;
  int centermax=centers.size();
  //int momax=occupation.GetDim(0);
  newvals=0;
  //assert(momax <= nmo);
  assert(e < sample->electronSize());
  assert(newvals.GetDim(1) >=5);

  // cout << "array " << endl;
  Array1 <doublevar> fval(3);
  Array1 <doublevar> R(5);
  // cout << "symvals " << endl;
  static Array2 <doublevar> symmvals_temp(maxbasis,5);

  //cout << "arrayref " << endl;
  //References to make the code easier to read and slightly faster.
  Array1 <int> & basismotmp(basismo_list(listnum));
  Array2 <int> & basisfilltmp(basisfill_list(listnum));
  Array2 <dcomplex> & moCoefftmp(moCoeff_list(listnum));

  Basis_function * tempbasis;

  dcomplex c;
  int scaleval=0, scalesymm=0;
  int mo=0;
  int scalebasis=basisfilltmp.GetDim(1);
  centers.updateDistance(e, sample);
  int totfunc=0;
  int b;
  for(int ion=0; ion < centermax; ion++)
  {
    centers.getDistance(e, ion, R);
    //cout << "R  " << R(0) << endl;
    for(int n=0; n< centers.nbasis(ion); n++)
    {
      b=centers.basis(ion, n);
      tempbasis=basis(b);
      if(R(0) < obj_cutoff(b)) {
      tempbasis->calcLap(R, symmvals_temp);

      //cout << "symmvals: " << endl;
      //output_array(symmvals_temp);


      int imax=nfunctions(b);
      for(int i=0; i< imax; i++)
      {
        //cout << "i " << i << endl;

        int reducedbasis=scalebasis*totfunc;
        scalesymm=i*5;
	//if(e==60) cout << "R " << R(0) << "  cutoff " << cutoff(totfunc) <<  endl;
        if(R(0) < cutoff(totfunc))
        {
          //cout << "*****" << endl;
          //cout << "center " << ion << endl;
          //cout << "distance " << R(0) << endl;
          //cout << "function " << symmvals_temp(i,0) << endl;

          //cout << "newvals before " << totfunc << endl;
          //output_array(newvals);
          //cout << "coefficients " << endl;
          for(int basmo=0; basmo < basismotmp.v[totfunc]; basmo++)
          {
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


            for(int j=0; j< 5; j++)
            {
              //newvals(mo,j)+=c*symmvals(fn,j);
              newvals.v[scaleval+j]+=c*symmvals_temp.v[scalesymm+j];
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

  //cout << newvals(0,0) << endl;

}

//--------------------------------------------------------------------------

void MO_matrix_Ccutoff::updateHessian(Sample_point * sample, int e,
				      int listnum,
				      //const Array1 <int> & occupation,
				      //!<A list of the MO's to evaluate
				      Array2 <dcomplex> & newvals
				      //!< The return: in form (MO, [val, grad, dxx,dyy,...])
				      )
{

  int centermax=centers.size();
  newvals=0;
  assert(e < sample->electronSize());
  assert(newvals.GetDim(1)==10);
  
  Array1 <doublevar> fval(3);
  Array1 <doublevar> R(5);
  static Array2 <doublevar> symmvals_temp(maxbasis,10);

  //References to make the code easier to read and slightly faster.
  Array1 <int> & basismotmp(basismo_list(listnum));
  Array2 <int> & basisfilltmp(basisfill_list(listnum));
  Array2 <dcomplex> & moCoefftmp(moCoeff_list(listnum));

  Basis_function * tempbasis;

  dcomplex c;
  int scaleval=0, scalesymm=0;
  int mo=0;
  int scalebasis=basisfilltmp.GetDim(1);
  centers.updateDistance(e, sample);
  int totfunc=0;
  int b;
  for(int ion=0; ion < centermax; ion++)  {
    centers.getDistance(e, ion, R);
    for(int n=0; n< centers.nbasis(ion); n++)  {
      b=centers.basis(ion, n);
      tempbasis=basis(b);
      if(R(0) < obj_cutoff(b)) {
	tempbasis->calcHessian(R, symmvals_temp);
	
	int imax=nfunctions(b);
	for(int i=0; i< imax; i++)  {
	  int reducedbasis=scalebasis*totfunc;
	  scalesymm=i*10;
	  if(R(0) < cutoff(totfunc))  {
	    for(int basmo=0; basmo < basismotmp.v[totfunc]; basmo++)   {
	      mo=basisfilltmp.v[reducedbasis+basmo];
	      c=moCoefftmp.v[reducedbasis+basmo];
	      scaleval=mo*10;
	      for(int j=0; j< 10; j++) {
		newvals.v[scaleval+j]+=c*symmvals_temp.v[scalesymm+j];
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
