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
#include "qmc_io.h"
#include <iomanip>
#include "Cubic_spline.h"

string Cubic_spline::symmetry_lookup(symmetry_type s) { 
   switch(s)
   {
     case sym_S:
       return string("S");
       break;
     case sym_P:
       return string("P");
       break;
     case sym_5D:
       return string("5D");
       break;
     case sym_6D:
       return string("6D");
       break;
     case sym_7F:
       return string("7F");
       break;
     case sym_7F_crystal://CRYSTAL F orbital
       return string("7F_crystal"); 
       break; 
     case sym_10F:
       return string("10F");
       break;
     case sym_9G:
       return string("9G");
       break;
     case sym_15G:
       return string("15G");
       break;
     case sym_P_siesta:
       return string("P_siesta");
       break;
     case sym_D_siesta:
       return string("5D_siesta");
       break;
     case sym_F_siesta:
       return string("7F_siesta");
       break;
   default:
       error("Cubic_spline::symmetry_lookup found unknown symmetry");
   }
}


int Cubic_spline::symmetry_lvalue(symmetry_type s) { 
  switch(s)
  {
    case sym_S:
      return 0;
      break;
    case sym_P:
    case sym_P_siesta:
      return 1;
      break;
    case sym_5D:
    case sym_6D:
    case sym_D_siesta:
      return 2;
      break;
    case sym_7F:
    case sym_10F:
    case sym_7F_crystal://CRYSTAL F orbital
    case sym_F_siesta:
      return 3;
      break;
    case sym_9G:
    case sym_15G:
      return 4;
      break;
    default:
      error("Cubic_spline::symmetry_lookup found unknown symmetry");
  }
}



Cubic_spline::symmetry_type Cubic_spline::symmetry_lookup(string & s) { 
  if(caseless_eq(s, "S"))
    return sym_S;
  if(caseless_eq(s, "P"))
    return sym_P;
  if(caseless_eq(s, "D")|| caseless_eq(s,"6D"))
    return sym_6D;
  if(caseless_eq(s,"F") || caseless_eq(s,"10F")) 
    return sym_10F;
  if(caseless_eq(s, "G")|| caseless_eq(s,"15G"))
    return sym_15G;
  if(caseless_eq(s,"5D")) 
    return sym_5D;
  if(caseless_eq(s,"7F"))
    return sym_7F;
  if(caseless_eq(s,"9G"))
    return sym_9G;
  if (caseless_eq(s, "7F_crystal"))//CRYSTAL F orbital
    return sym_7F_crystal; //CRYSTAL F orbital
  if(caseless_eq(s,"P_siesta"))
    return sym_P_siesta;
  if(caseless_eq(s,"5D_siesta"))
    return sym_D_siesta;
  if(caseless_eq(s,"7F_siesta"))
    return sym_F_siesta;  
  error("I don't understand symmetry type ", s, "\nShould be S,P,D,etc.");
  return sym_S;
}


//-------------------------------------------------------------------------
int Cubic_spline::read(
  vector <string> & words,
  unsigned int & pos
)
{

  Array1 <double> basisparms(4);
  atomname=words[0];
  basisparms(0)=.02;  //spacing
  basisparms(1)=2500; //number of points
  basisparms(2)=1e31;
  basisparms(3)=1e31;
  if(!readvalue(words, pos=0, norm_type, "NORMTYPE")) {
    norm_type="GAMESSNORM";
  }

  if(haskeyword(words, pos=0, "NORENORMALIZE")) {
    renormalize=false;
  }
  else {
    renormalize=true;
  }

  doublevar cutmax;
  if(readvalue(words, pos=0, cutmax, "CUTOFF")) {
    requested_cutoff=cutmax;
  }
  else {
    requested_cutoff=-1;
  }

  enforce_cusp=false;
  if(readvalue(words, pos=0, cusp, "CUSP")) enforce_cusp=true;

  match_sto=false;
  if(readvalue(words, pos=0, cusp_matching, "CUSP_MATCHING")) match_sto=true;
  
  zero_derivative=false;
  if(haskeyword(words, pos=0, "ZERO_DERIVATIVE")) zero_derivative=true;

  customspacing=basisparms(0);
  if(readvalue(words, pos=0, customspacing, "SPACING")) {
    basisparms(0)=customspacing;
    basisparms(1)=50.0/customspacing;
  }

  vector <string> basisspec;
  string read_file; //read positions from a file
  if(readsection(words,pos=0, basisspec, "GAMESS"))
  {
    unsigned int newpos=0;
    return readbasis(basisspec, newpos, basisparms);
  }
  else {
    pos=0;
    if(readspline(words)) {
      return 1;
    }
    else {
      error("couldn't find proper basis section in AOSPLINE.  "
            "Try GAMESS { .. }.");
      return 0;
    }
  }

  //We assume that all splines use the same interval,
  //which saves a few floating divisions, so check to make
  //sure the assumption is good.
  for(int i=0; i< splines.GetDim(0); i++) {
    if(!splines(0).match(splines(i))) 
      error("spline ", i, " doesn't match the first one ");
  }
}


//-------------------------------------------------------------------------
int Cubic_spline::readspline(vector <string> & words) {
  unsigned int pos=0;
  vector <vector <string> > splinefits;
  vector <string> onesection;
  while(readsection(words, pos, onesection, "SPLINE")) {
    splinefits.push_back(onesection);
  }
  if(splinefits.size()==0) {
    return 0;
  }

  nsplines=splinefits.size();
  splines.Resize(nsplines);
  symmetry.Resize(nsplines);
  for(int i=0; i< nsplines; i++) { 
    symmetry(i)=symmetry_lookup(splinefits[i][0]);
    //splinefits[i].erase(splinefits[i].begin());
  }
  assign_indiv_symmetries();

  //cout << "created " << nfunc() << " total functions " << endl;

  double max_support=0.0;
  for(int i=0; i< nsplines; i++) { 
    int n=splinefits[i].size();
    double supp=atof(splinefits[i][n-2].c_str());
    if(supp > max_support) max_support=supp;
  }
  //cout << "maximum support " << max_support <<endl;
  
  
  for(int s=0; s< nsplines; s++) {
    splines(s).readspline(splinefits[s], enforce_cusp, cusp/double(symmetry_lvalue(symmetry(s))+1));
    if(requested_cutoff > 0) { 
      splines(s).enforceCutoff(requested_cutoff);
    }
  }

  findCutoffs();
  return nfunc();
}
//-------------------------------------------------------------------------
void Cubic_spline::raw_input(ifstream & input)
{

  splines.Resize(nsplines);
  for(int s=0; s< nsplines; s++) {
    splines(s).raw_input(input);
  }
}


//-------------------------------------------------------------



int Cubic_spline::nfunc()
{

  int totf=0; //running total of number of functions.
  //return number of functions with symmetry

  for(int i=0; i< nsplines; i++)
  {
    switch(symmetry(i))
    {
    case sym_S:
      totf+= 1;
      break;
    case sym_P:
    case sym_P_siesta:
      totf+= 3;
      break;
    case sym_5D:
    case sym_D_siesta:
      totf+=5;
      break;
    case sym_6D:
      totf+= 6;
      break;
    case sym_7F:
    case sym_7F_crystal://CRYSTAL F orbital
    case sym_F_siesta:
      totf+=7;
      break;
    case sym_10F:
      totf+= 10;
      break;
    case sym_15G:
      totf+=15;
      break;
    case sym_9G:
      totf+=9;
      break;
    default:
      cout << "Bad symmetry in Cubic_spline::nfunc()\n";
      cout << symmetry(i) << endl;
      return 0;
    }
  }
  return totf;
}

//----------------------------------------------------------------------

doublevar Cubic_spline::cutoff(int n){
  assert(n<nfunctions);
  return rcut(n);
}


void Cubic_spline::findCutoffs()
{
  //Find the cutoff lengths for each function.
  rcut.Resize(nfunctions);
  int totfunc=0;
  for(int i=0; i < nsplines; i++) {
    int nrep=0;
    switch(symmetry(i)) {
    case sym_S:
      nrep=1;
      break;
    case sym_P:
    case sym_P_siesta:
      nrep=3;
      break;
    case sym_5D:
    case sym_D_siesta:
      nrep=5;
      break;
    case sym_6D:
      nrep=6;
      break;
    case sym_7F:
    case sym_7F_crystal://CRYSTAL F orbital
    case sym_F_siesta:
      nrep=7;
      break;
    case sym_10F:
      nrep=10;
      break;
    case sym_9G:
      nrep=9;
      break;
    case sym_15G:
      nrep=15;
      break;
    default:
      error("I don't have this symmetry:", symmetry(i));
    }
    doublevar cutoff=splines(i).findCutoff();
    for(int j=0; j< nrep; j++) {
      rcut(totfunc++)=cutoff;
    }
  }

  threshold=0;
  for(int i=0; i< rcut.GetDim(0); i++) {
    if(threshold < rcut(i)) threshold=rcut(i);
  }
  
  for(int i=0; i< nsplines; i++) { 
    splines(i).pad(threshold);
  }
}
//-------------------------------------------------------


/*!
This takes a vector of words
presumably gotten from an input file and parses them into a basis function.
It takes something very similar to GAMESS input, but it does not allow you
to have a normalization constant.(ie, all function headers should be like
S 1, not S 1 1.0)  It will normalize it anyway.



*/
int Cubic_spline::readbasis(vector <string> & words,unsigned int & pos,
                            Array1 <double> & parms){

  doublevar spacing=parms(0);
  int n= (int) parms(1);
  double yp1=parms(2);
  double ypn=parms(3);
  int symmtype=0;

  if(parms.GetDim(0) > 4)
    symmtype=(int) parms(4);

  assert(symmtype==0 || symmtype==1);

  string symm;
  int ngauss;

  vector <symmetry_type> symmetry_temp;

  for(unsigned int p=pos; p<words.size(); p++)  {
    symm=words[p];
    ngauss=atoi(words[++p].c_str());
    symmetry_temp.push_back(symmetry_lookup(symm));

    vector <doublevar> exp(ngauss);
    vector <doublevar> c(ngauss);
    if(p+ngauss*3 >= words.size()) {
      error("Unexpected end of GAMESS section.  Count is wrong.");
    }
    for(int i=0; i<ngauss; i++) {
      p++;
      exp[i]=atof(words[++p].c_str());
      c[i]=atof(words[++p].c_str());
    }
    exponent.push_back(exp);
    coefficient.push_back(c);
  }


  nsplines=exponent.size();
  Array1 <doublevar> y(n);
  Array1 <doublevar> x(n);
  splines.Resize(nsplines);
  symmetry.Resize(nsplines);
  for(int i=0; i< nsplines; i++) {
    symmetry(i)=symmetry_temp[i];
  }

  //cout << "nfunc " << nfunc() << endl;
  
  const double normtol=1e-5; //tolerance for the normalization of the basis fn
  for(int funcNum=0; funcNum<nsplines; funcNum++)  {
    //calculate interpolation points
    for(int j=0; j<n; j++) {
      x(j)=j*spacing;
    }
    y=0;
    for(unsigned int i=0; i < exponent[funcNum].size(); i++) {
      doublevar norm;
      if(norm_type=="GAMESSNORM") {
        doublevar fac=sqrt(2.*exponent[funcNum][i]/pi);
        doublevar feg=4.*exponent[funcNum][i];
        doublevar feg2=feg*feg;
        
        switch(symmetry(funcNum))   {
          case sym_S:
            norm=sqrt(2.*feg*fac);
            break;
          case sym_P:
            norm=sqrt(2.*feg2*fac/3.);
            break;
          case sym_5D:
          case sym_D_siesta:
          case sym_6D:
            norm=sqrt(2.*feg*feg2*fac/15.);
            break;
          case sym_7F:
          case sym_7F_crystal://CRYSTAL F orbital
          case sym_F_siesta:
          case sym_10F:
            norm=sqrt(2.*feg2*feg2*fac/105.);
            break;
          case sym_9G:
          case sym_15G:
            norm=sqrt(2.*feg2*feg2*feg*fac/945.); //Lubos-done
            break;
            /*  general formula here with n principal quantum numbers
             norm=sqrt[(2 ((4exponent)**n) sqrt(2exponent/pi))/(2n-1)!!]
             */
          default:
            norm=0;
            error("Unknown symmetry in Cubic_spline::readbasis! Shouldn't be here!");
        }
      }
      else if(norm_type=="CRYSTAL") {
        switch(symmetry(funcNum))
        {
          case sym_S:
            norm=1;
            break;
          case sym_P:
            norm=1/sqrt(3.0);
            break;
          case sym_5D:
          case sym_D_siesta:
          case sym_6D:
            norm=1/sqrt(5.0);
            break;
          case sym_7F:
          case sym_7F_crystal://CRYSTAL F orbital
          case sym_F_siesta:
          case sym_10F:
            norm=1/sqrt(7.0);
            break;
          case sym_9G:
          case sym_15G:
            norm=1/sqrt(9.0);
            break;
          default:
            norm=0;
            error("unknown symmetry in cubic_spline: ", symmetry(funcNum));
        }
      }
      else if(norm_type=="NONE") {
        norm=1;
      }
      else {
        norm=0;
        error("Didn't understand NORMTYPE ", norm_type);
      }


      for(int j=0; j<n; j++)
      {
        y(j)+=norm*coefficient[funcNum][i]
              *exp(-exponent[funcNum][i]*x(j)*x(j));
      }
    }

    //Renormalize to one
    if(renormalize) {
      //check normalization
      doublevar sum=0;
      for(int j=0; j<n; j++)  {
        int p=0;
        switch(symmetry(funcNum)) {
          case sym_S:
            p=1;
            break;
          case sym_P:
            p=2;
            break;
          case sym_5D:
          case sym_D_siesta:
          case sym_6D:
            p=3;
            break;
          case sym_7F:
          case sym_7F_crystal://CRYSTAL F orbital
          case sym_F_siesta:
          case sym_10F:
            p=4;
            break;
          case sym_9G:
          case sym_15G:
            p=5;
            break;
          default:
            error("unknown symmetry when checking the normalization");
        }

        //doublevar power=pow(x(j), 2*(symmetry(funcNum)+1) );
        doublevar power=pow(x(j), 2*p );
        sum+=spacing*power*y(j)*y(j);
      }
      //cout << "normalization : " << sqrt(sum) << endl;

      //renormalize if needed

      if(fabs(sum-1) > normtol){
        //debug_write(cout, "Normalizing basis function, sum: ", sum, "\n");
        for(int j=0; j<n; j++){
          y(j)=y(j)/sqrt(sum);
        }
      }
    }

    //Force the function to go to zero if the user requested a cutoff
    if(requested_cutoff > 0) {
      const double smooth=1.2;
      const double cutmax=requested_cutoff-1e-6;
      const double cutmin=cutmax-smooth;
      for(int j=0; j< n; j++) {
        if(x(j) > cutmin) {
          if(x(j) > cutmax) { //if we're beyond the cutoff completely
            y(j)=0;
          }
          else {  //if we're in the smooth cutoff region
            double zz=(x(j)-cutmin)/smooth;
            y(j) *= (1-zz*zz*zz*(6*zz*zz-15*zz+10));
          }
        }
      }
    }


    //--Estimate the derivatives on the boundaries
    if ( zero_derivative ) {
      yp1=0.0;
    } else {
      yp1=(y(1)-y(0))/spacing;
    }
    ypn=(y(n-1) - y(n-2) ) /spacing;

    if(enforce_cusp) { 
      double der=cusp/double(symmetry_lvalue(symmetry(funcNum))+1);      
      yp1=der;

      // KMR:  using an STO near the atom should be a separate decision from 
      // passing the first derivative to the spliner
      if(match_sto) { //inspired from J. Chem. Phys. 130, 114107 (2009)
        //Here, we match the value, first derivative, and second derivative to the 
        //given smooth function at some correction radius rc.  We can then safely
        //replace the function from [0:rc] with the Slater function.
        //cout << "enforcing cusp" << endl;
        double rc=cusp_matching;
        int closest=rc/spacing;
        rc=x(closest);
        //cout << "rc " << rc << endl;
        double curve=(y(closest+1)+y(closest-1)-2*y(closest))/(spacing*spacing);
        double deriv=(y(closest+1)-y(closest-1))/(2*spacing);
        double f=y(closest);
        double b=(deriv*der*der-curve*der)/(curve*der*rc*rc+curve*rc*2-deriv*der*der*rc*rc-deriv*4*der*rc-2*deriv);
        double a=curve/(exp(der*rc)*(der*der*(1+b*rc*rc)+4*der*b*rc+2*b));
        double c=f-a*exp(der*rc)*(1+b*rc*rc);
        //cout << "a " << a << " b " << b << " c " << c << endl;
        for(int j=0; j <= closest; j++) {
          y(j)=a*exp(der*x(j))*(1+b*x(j)*x(j))+c;
        }
      } 
    }
    /*
    for(int j=0; j< n; j++) { 
      cout << "jjjjjjjj " << x(j) << "  " << y(j) << endl;
    }
    cout << "jjjjjjj " << endl;
    */
    
    splines(funcNum).splinefit(x,y,yp1, ypn);
  }
  nfunctions=nfunc();


  assign_indiv_symmetries();
  findCutoffs();
  return nfunctions;
}


//-------------------------------------------------------



void Cubic_spline::assign_indiv_symmetries() {
  nfunctions=nfunc();
  nsplines=symmetry.GetDim(0);
  //---------------------------------------------
  //Assign the individual symmetries for the various codes
  indiv_symmetry.Resize(nfunctions);
  nfuncspline.Resize(nsplines);
  int totfunc=0;
  for(int funcNum=0; funcNum < nsplines; funcNum++){
    //cout << "funcNum " << funcNum << endl;
    //cout << "totfunc " << totfunc << endl;
    switch(symmetry(funcNum))
    {
      case sym_S:
        nfuncspline(funcNum)=1;
        indiv_symmetry(totfunc++)=isym_S;
        break;
      case sym_P:
        nfuncspline(funcNum)=3;
        indiv_symmetry(totfunc++)=isym_Px;
        indiv_symmetry(totfunc++)=isym_Py;
        indiv_symmetry(totfunc++)=isym_Pz;
        break;
      case sym_5D:
        nfuncspline(funcNum)=5;  //Crystal ordering
        indiv_symmetry(totfunc++)=isym_Dz2r2;
        indiv_symmetry(totfunc++)=isym_Dxz;
        indiv_symmetry(totfunc++)=isym_Dyz;
        indiv_symmetry(totfunc++)=isym_Dx2y2;
        indiv_symmetry(totfunc++)=isym_Dxy;
        break;
      case sym_6D:
        nfuncspline(funcNum)=6;
        indiv_symmetry(totfunc++)=isym_Dxx;
        indiv_symmetry(totfunc++)=isym_Dyy;
        indiv_symmetry(totfunc++)=isym_Dzz;
        indiv_symmetry(totfunc++)=isym_Dxy;
        indiv_symmetry(totfunc++)=isym_Dxz;
        indiv_symmetry(totfunc++)=isym_Dyz;
        break;
      case sym_7F_crystal://CRYSTAL F orbital
        nfuncspline(funcNum)=7;
        indiv_symmetry(totfunc++)=isym_F0;
        indiv_symmetry(totfunc++)=isym_Fp1;
        indiv_symmetry(totfunc++)=isym_Fm1;
        indiv_symmetry(totfunc++)=isym_Fp2;
        indiv_symmetry(totfunc++)=isym_Fxyz;
        indiv_symmetry(totfunc++)=isym_Fp3mod;
        indiv_symmetry(totfunc++)=isym_Fm3;
      break;
      case sym_7F:
        nfuncspline(funcNum)=7;
        indiv_symmetry(totfunc++)=isym_F0;
        indiv_symmetry(totfunc++)=isym_Fm3;
        indiv_symmetry(totfunc++)=isym_Fp3mod;
        indiv_symmetry(totfunc++)=isym_Fp2;
        indiv_symmetry(totfunc++)=isym_Fxyz;
        indiv_symmetry(totfunc++)=isym_Fm1;
        indiv_symmetry(totfunc++)=isym_Fp1;
        break;
      case sym_10F:
        nfuncspline(funcNum)=10;
        indiv_symmetry(totfunc++)=isym_Fxxx;
        indiv_symmetry(totfunc++)=isym_Fyyy;
        indiv_symmetry(totfunc++)=isym_Fzzz;
        indiv_symmetry(totfunc++)=isym_Fxxy;
        indiv_symmetry(totfunc++)=isym_Fxxz;
        indiv_symmetry(totfunc++)=isym_Fyyx;
        indiv_symmetry(totfunc++)=isym_Fyyz;
        indiv_symmetry(totfunc++)=isym_Fzzx;
        indiv_symmetry(totfunc++)=isym_Fzzy;
        indiv_symmetry(totfunc++)=isym_Fxyz;
        break;
      case sym_9G:
        nfuncspline(funcNum)=9;
        indiv_symmetry(totfunc++)=isym_G0;
        indiv_symmetry(totfunc++)=isym_G1;
        indiv_symmetry(totfunc++)=isym_G2;
        indiv_symmetry(totfunc++)=isym_G3;
        indiv_symmetry(totfunc++)=isym_G4;
        indiv_symmetry(totfunc++)=isym_G5;
        indiv_symmetry(totfunc++)=isym_G6;
        indiv_symmetry(totfunc++)=isym_G7;
        indiv_symmetry(totfunc++)=isym_G8;
        break;
      case sym_15G:
        nfuncspline(funcNum)=15;
        indiv_symmetry(totfunc++)=isym_Gxxxx;
        indiv_symmetry(totfunc++)=isym_Gyyyy;
        indiv_symmetry(totfunc++)=isym_Gzzzz;
        indiv_symmetry(totfunc++)=isym_Gxxxy;
        indiv_symmetry(totfunc++)=isym_Gxxxz;
        indiv_symmetry(totfunc++)=isym_Gyyyx;
        indiv_symmetry(totfunc++)=isym_Gyyyz;
        indiv_symmetry(totfunc++)=isym_Gzzzx;
        indiv_symmetry(totfunc++)=isym_Gzzzy;
        indiv_symmetry(totfunc++)=isym_Gxxyy;
        indiv_symmetry(totfunc++)=isym_Gxxzz;
        indiv_symmetry(totfunc++)=isym_Gyyzz;
        indiv_symmetry(totfunc++)=isym_Gxxyz;
        indiv_symmetry(totfunc++)=isym_Gyyxz;
        indiv_symmetry(totfunc++)=isym_Gzzxy;
        break;
      case sym_P_siesta:
        nfuncspline(funcNum)=3;
        indiv_symmetry(totfunc++)=isym_Py;
        indiv_symmetry(totfunc++)=isym_Pz;
        indiv_symmetry(totfunc++)=isym_Px;
        break;
      case sym_D_siesta:
        nfuncspline(funcNum)=5;
        indiv_symmetry(totfunc++)=isym_Dxy;
        indiv_symmetry(totfunc++)=isym_Dyz;
        indiv_symmetry(totfunc++)=isym_Dz2r2;
        indiv_symmetry(totfunc++)=isym_Dxz;
        indiv_symmetry(totfunc++)=isym_Dx2y2;
        break;
      case sym_F_siesta:
        nfuncspline(funcNum)=7;
        indiv_symmetry(totfunc++)=isym_Fm3;
        indiv_symmetry(totfunc++)=isym_Fxyz;
        indiv_symmetry(totfunc++)=isym_Fm1;
        indiv_symmetry(totfunc++)=isym_F0;
        indiv_symmetry(totfunc++)=isym_Fp1;
        indiv_symmetry(totfunc++)=isym_Fp2;
        indiv_symmetry(totfunc++)=isym_Fp3mod;
        break;
      default:
        error("unknown symmetry type: ", symmetry(funcNum));
        
    }
  }
}


//----------------------------------------------------------------------

int Cubic_spline::showinfo(string & indent, ostream & os)
{
  os << indent << "Cubic spline for " << atomname << endl;
  if(zero_derivative)
    os << indent << "Zero derivative enforced at the orgin\n";
  if(customspacing!=0.02)
    os << indent << "Using custom spacing of "<<customspacing<<endl;
  os << indent << nsplines << "  radial functions\n";
  os << indent << setw(10) << "function" << setw(10) << "symmetry"
     << setw(10) << "cutoff" << endl;
  
  int totfunc=0;
  for(int i=0; i< nsplines; i++) {
    os << indent << setw(10) << i << setw(10) << symmetry_lookup(symmetry(i))
       << setw(10) << rcut(totfunc) << endl;
    totfunc+=nfuncspline(i);
  }
  


  return 1;
}

//----------------------------------------------------------------------

int Cubic_spline::writeinput(string & indent, ostream & os)
{
  os << indent << atomname << endl;
  os << indent << "AOSPLINE \n";
  os << indent << "NORMTYPE  " << norm_type << endl;
  if(renormalize==false) {
    os << indent << "NORENORMALIZE\n";
  }
  if(zero_derivative) {
    os << indent << "ZERO_DERIVATIVE\n";
  }
  if(customspacing!=0.02)
    os << indent << "SPACING "<<customspacing<<endl;
  if(requested_cutoff > 0 ) {
    os << indent << "CUTOFF " << requested_cutoff << endl;
  }
  if(enforce_cusp) { 
    os << indent << "CUSP " << cusp << endl;
    os << indent << "CUSP_MATCHING " << cusp_matching << endl;
  }
  if(exponent.size() > 0) { 
    os << indent << "GAMESS { \n";
    
    for(int funcNum=0; funcNum<nsplines; funcNum++) {

      os << indent << symmetry_lookup(symmetry(funcNum)) << endl;
      os << indent << exponent[funcNum].size() << endl;
      for(unsigned int i=0; i < exponent[funcNum].size(); i++) {
        os << indent << "  " << i << "  " << exponent[funcNum][i]
        << "   " << coefficient[funcNum][i] << endl;
      }
    }
    os << indent << " } \n";
  }
  else { 
    for(int funcNum=0; funcNum < nsplines; funcNum++) { 
      os << indent << "SPLINE { \n";
      os << indent << symmetry_lookup(symmetry(funcNum)) << endl;
      splines(funcNum).writeinput(indent, os);
      os << indent << " }\n";
    }
  }

  return 1;
}

//------------------------------------------------------------------------
