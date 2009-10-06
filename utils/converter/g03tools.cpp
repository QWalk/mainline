/*
 
 Copyright (C) 2007 Hiori Kino
 
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

#include <fstream>
#include <iostream>
#include <iomanip>
#include <map>
#include <string>
#include <vector>
#include <cmath>

#include <string.h>
#include <stdlib.h>

#include "converter.h"
#include "Pseudo_writer.h"
#include "basis_writer.h"
#include "g03tools.h"



using namespace std;
/************************************************************************************
 *
 * **********************************************************************************/


/* length of char *, but ignore ' ' at the end of char * */
int strlen_trim(/* input */const char *c)
{
  int n=strlen(c);
  if (n==0) { return 0;}
  char const *endc;
  endc=&c[n-1];
  for ( ; endc!=c;endc--) 
  {
     if (*endc!= ' ') break;
  }
  if (endc==c) { return 0; }
  else {
  return (int) (endc-c+1);
  }
}



// atof can not read '0.1D+1' format
// atod can read '0.1D+1' 
double atod(/* input */const char *s)
{
    char *str;
    str=strdup(s);
    char *p;
    for (p=str;*p!='\0';p++) {
       if (*p =='D' || *p=='d' ) *p='E';
    }
    double v= atof(str);
    free(str);
    return v;
}



// do not use atof, use atod
#define atof(x) atod(x)





/* L(integer) -> L(char *) */
const char *L2Lstr(int L)
{
        const char *t="L2Lstr" ;

   if ( L<0 || L>=nlabel_L ) {
        cout << t << " error L out of range "<< L <<endl;
        exit(ERR_CODE);
   }
   return label_Lstr[L];
}



/* "S", "P", ...  -> 0,1,... */
static int g03_angularMomentum(/* input */ string & L)
{
  const char *t="g03_angularMomentum: ";
  int found=0;
  int i;
  for (i=0;i<nlabel_L;i++) {
    if ( L == label_Lstr[i] ) { found=1; break; }
  }
  if (found) {
    return i; 
  }
  else {
    cout << t<<"error in augular momentum "<< L <<endl;
    exit(ERR_CODE);
  }
}

/* S,P,SP,D,... -> 0,1,1,2,... */
int g03_maxangularMomentum(/* input */ string & L)
{
  const char *t="g03_maxangularMomentum: ";
  int found=0;
  int i;

  if (L=="SP" || L=="L" ) { return 1; }

  for (i=0;i<nlabel_L;i++) {
    if ( L == label_Lstr[i] ) { found=1; break; }
  }
  if (found) {
    return i;
  }
  else {
    cout << t<<"error in augular momentum "<< L <<endl;
    exit(ERR_CODE);
  }
}

/* S,P,SP,6D,.. -> 1 3 4 5, 
 * 5D->5,, 6D->6
 */
int g03_basissizeL(/*input*/ string & L )
{
   static const char *t="g03_basissizeL: ";
   int ret;
   if      (L=="S") { ret=1;}
   else if (L=="P") { ret=3;}
   else if (L=="SP" || L=="L" ) {ret=4;}
   else if (L=="D" ||  L=="6D") { ret=6;}
   else if (L=="5D") { ret=5;}
   else if (L=="F" ||  L=="10F") {ret=10;}
   else if (L=="7F") { ret=7;}
   else {
      cout << t << "L= \""<<L<<"\" not supported"<<endl;
   }
   return ret;
}


/* 'P and up' or 'S - P' to set of '1 1' or '0 1' */
 int label_to_L(/* input */string & label, 
		/* output */ int & refL)
{
  const char *t = "label_to_L: ";
  string spc=" ";
  vector<string> words; 

  split(label, spc, words);

  int nwords=words.size();
  int L; // return value

  if (nwords==1) {
     L= g03_angularMomentum(words[0]);
    refL = L;
  }
  else if (nwords==3) {
  // P and up,  type
  if ( words[1]=="and" && words[2] =="up" ) { // local part
    L= g03_angularMomentum(words[0]);
    refL = L;
  }
  //  S - P, type
  else if (words[1]=="-") {
    L = g03_angularMomentum(words[0]);
    refL =  g03_angularMomentum(words[2]); 
  }
  else {
    cout << t<< " unknown label " << label << endl;
    exit(ERR_CODE);
  }
  }
  return L;
}
//
//  from crystal2qmc 
//
double cutoff_divider(vector<Gaussian_basis_set> &basis,  vector< vector<double> > & latvec) 
{
  const char *t="cutoff_divider: ";
  if ( basis.empty() ) {
     cout << t<< "error, basis is empty" <<endl;
     exit(ERR_CODE); 
  }
  if ( latvec.empty() ) {
     cout << t<< "error, latvec is empty" <<endl;
     exit(ERR_CODE);
  }
  double min=1e8;

  for(vector<Gaussian_basis_set>::iterator bas=basis.begin();
      bas != basis.end(); bas++) {
    for(vector <vector<double> >::iterator i=bas->exponents.begin();
        i!=bas->exponents.end(); i++) {
      for(vector<double>::iterator j=i->begin(); j!= i->end(); j++) {
        if(*j < min) min=*j;
      }
    }
  }

  cout << "minimum exponent " << min << endl;
  double cutoff_length= sqrt(-log(1e-8)/min);
  cout << "cutoff length " << cutoff_length << endl;

    double min_latsize=1e8;
    //sysout << "LATTICEVEC { \n";
    for(int i=0; i< 3; i++) {
      double length=0;
      for(int j=0; j< 3; j++) {
        //sysout << latvec[i][j] << "   ";
        length+=latvec[i][j]*latvec[i][j];
      }
      //sysout << endl;
      if(min_latsize > length) min_latsize=length;
   }
  
  double ret = sqrt(min_latsize)/cutoff_length;  

  return ret;
}

/*************************************************************************************
 *
 * ***********************************************************************************/



static int gamessvec_normalize(/* input */
		vector<Atom> & atoms, 
		vector <Gaussian_basis_set> & basis,
		/* input output */
		vector < vector < double> > & moCoeff )
{
  const char *t="gamessvec_normalize:";
  int natoms=atoms.size();

  //Normalize
  const double pi=3.1415926535;
  double snorm=1./sqrt(4.*pi);
  double pnorm=sqrt(3.)*snorm;
  vector <double> dnorm;
  dnorm.push_back(sqrt(5./(4*pi)));
  dnorm.push_back(sqrt(5./(4*pi)));
  dnorm.push_back(sqrt(5./(4*pi)));
  dnorm.push_back(sqrt(15./(4*pi)));
  dnorm.push_back(sqrt(15./(4*pi)));
  dnorm.push_back(sqrt(15./(4*pi)));
  //
  vector <double> fnorm;
  for(int i=0; i< 3; i++) {
    fnorm.push_back(snorm*sqrt(7.)); //xxx, yyy, zzz
  }
  for(int i=3; i< 9; i++) {
    fnorm.push_back(snorm*sqrt(35.0)); //xxy, xxz, yyx, yyz, zzx, zzy
  }
  fnorm.push_back(snorm*sqrt(105.0)); //xyz
  //
  vector <double> gnorm; // I believe g norms are now correct-Lubos
  for(int i=0; i< 3; i++) {
    gnorm.push_back(snorm*sqrt(9.)); //xxxx, yyyy, zzzz
  }
  for(int i=3; i< 9; i++) {
    gnorm.push_back(snorm*sqrt(63.0)); //xxxy, xxxz, yyyx, yyyz, zzzx, zzzy
  }
  for(int i=9; i< 12; i++) {
    gnorm.push_back(snorm*sqrt(105.)); //xxyy, xxzz, yyzz
  }
  for(int i=12; i< 15; i++) {
    gnorm.push_back(snorm*sqrt(315.0)); //xxyz, yyxz, zzxy
  }
  //-----------------------after here
  //
  int totmo2=moCoeff.size();
  for(int mo=0; mo < totmo2; mo++) {
    int func=0;
    for(int at=0; at < natoms; at++) {
      int bas=atoms[at].basis;
      int nbasis=basis[bas].types.size();
      for(int i=0; i< nbasis; i++) {
        if(basis[bas].types[i] == "S") {
          moCoeff[mo][func]*=snorm;
          func++;
        }
        else if(basis[bas].types[i] == "P") {
          for(int j=0; j< 3; j++) {
            moCoeff[mo][func]*=pnorm;
            func++;
          }
        }
        else if(basis[bas].types[i] == "6D") {
          for(int j=0; j< 6; j++) {
            moCoeff[mo][func]*=dnorm[j];
            func++;
          }
        }
        else if(basis[bas].types[i] == "10F") {
          for(int j=0; j< 10; j++) {
            moCoeff[mo][func]*=fnorm[j];
            func++;
          }
        }
        else if(basis[bas].types[i] == "15G") {
          for(int j=0; j< 15; j++) {
            moCoeff[mo][func]*=gnorm[j];
            func++;
          }
        }
        else if(basis[bas].types[i]== "L") {
          moCoeff[mo][func] *=snorm;
          func++;
          for(int j=0; j< 3; j++) {
            moCoeff[mo][func]*=pnorm;
            func++;
          }
        }
        else {
          cout << t << "unknown type " << basis[bas].types[i] << endl;
          exit(1);
        }
      }
    }
    //cout <<" func " << func << endl;
  }

  return 0;
}

int g03_fix_form(
		/* input */
		int nup, int ndown, 
		string & scftype, 
		int dtype, int ftype, 
		vector< vector< double> > & alphamo,
		vector< vector< double> > & betamo,
		vector<PseudoValence> &pseudovalence, 
		vector <Atom>  & atoms,
		/*input and output */
	        vector <Gaussian_basis_set> & basis,
		vector <Gaussian_pseudo_writer> & pseudo,
		/* output */
                Slat_wf_writer & slwriter,
                vector < vector < double> > & moCoeff,
		vector <Center> & centers,
		vector <int> & basisidx 
		)
{
  const char *t = "g03_change_form";
  const int vorb=10;
  vector<double> emptyvector;

  int basissize = alphamo[0].size();

  // slwriter 
  slwriter.nup   = nup;
  slwriter.ndown = ndown;
  int norb = nup + vorb;
  if ( norb > basissize )  norb = basissize; 
  slwriter.spin_dwn_start= norb;
  slwriter.calctype = scftype;

#if 0
  // done when atoms are read 
  //  ang -> au
 for (vector<Atom>::iterator it=atoms.begin(); it!=atoms.end(); it++) {
	 for (int k=0;k<3;k++) {
             it->pos[k] /= ANG;
	 }
	 string_upper(it->name);
 }
#endif

 // upper case
#if 0
 for (vector<Gaussian_basis_set>::iterator it = basis.begin();  it!=basis.end(); it++) {
	 string_upper(it->label);
 }
 for (vector<Gaussian_pseudo_writer>::iterator it = pseudo.begin(); it!=pseudo.end(); it++) {
	 string_upper(it->label);
 }
#endif

  if ( !pseudo.empty() ) {
  for (vector<Atom>::iterator it=atoms.begin(); it!=atoms.end(); it++) {
      string name = it->name; 
      vector<PseudoValence>::iterator itps;
      int found=0;
      // find the same name 
      for (itps=pseudovalence.begin(); itps!=pseudovalence.end(); itps++) {
	      if (name==itps->label) {
		      found=1;
		      break;
	      }
      }
      if (found ) {
          it->charge = itps->valencenum; 
      }
  }
  }

  // "D"  -> 5D or 6D
  // "F" ->  7F or 10F 
  for (vector<Gaussian_basis_set>::iterator it = basis.begin();  it!=basis.end(); it++) {
     for (vector<string>::iterator itlabel=it->types.begin(); itlabel!=it->types.end(); itlabel++ ) {
	     if (*itlabel=="D") *itlabel=dtypestr[dtype];
	     if (*itlabel=="F") *itlabel=ftypestr[ftype];
     }
  }


  //QMC likes the local part at the end of the list, and GAMESS outputs
  //  //it at the beginning, so let's fix that..
  //
  int npseud=pseudo.size();
  for(int ps=0; ps < npseud; ps++) {
    if ( pseudo[ps].exponents.empty() ) continue;
    pseudo[ps].exponents.push_back(pseudo[ps].exponents[0]);
    pseudo[ps].nvalue.push_back(pseudo[ps].nvalue[0]);
    pseudo[ps].coefficients.push_back(pseudo[ps].coefficients[0]);
    pseudo[ps].exponents.erase(pseudo[ps].exponents.begin());
    pseudo[ps].coefficients.erase(pseudo[ps].coefficients.begin());
    pseudo[ps].nvalue.erase(pseudo[ps].nvalue.begin());
  }

#if 0
   cout << "---------- pseudo -----------" <<endl;
   for(int ps=0; ps < npseud; ps++) {
	   pseudo[ps].print_pseudo(cout);
   } 
#endif

   // moCoeff
   // moCoeff = alphamo;  copy only norb of alphamo
   cout << "norb="<< norb<<endl;
   int iorb=0;
   for ( vector< vector<double> >::iterator vec=alphamo.begin();
                   vec != alphamo.end(); vec++ ) {
         moCoeff.push_back( *vec ); iorb++;
         if ( iorb >= norb ) break; 
   }
   cout << "spin up end size="<< moCoeff.size() <<endl;
   iorb=0;
   for ( vector< vector<double> >::iterator vec=betamo.begin(); 
		   vec != betamo.end(); vec++ ) {
         moCoeff.push_back( *vec ); iorb++;
         if ( iorb>= norb ) break;
   }
   cout << "moCoeff total size="<< moCoeff.size()<<endl;

   {
     int factor;
   if ( scftype =="UHF" ) {
      factor=2;
   }
   else {
      factor=1;
   }
      if ( moCoeff.size() != norb*factor ) {
          cout << t << " error size mismatch " << moCoeff.size() <<endl;
	   exit(ERR_CODE);
      }
   }

   gamessvec_normalize(/* input */
                   atoms,
                   basis,
                   /* input output */
                    moCoeff );

  int natoms=atoms.size();; 
  centers.resize(natoms);
  basisidx.resize(natoms);
  for(int at=0; at < natoms; at++) {
    for(int i=0; i< 3; i++) centers[at].pos[i]=atoms[at].pos[i];
    centers[at].equiv_atom=at;
    centers[at].name=atoms[at].name;
    basisidx[at]=basis[atoms[at].basis].nfunc();
  }

#if 0
   cout << "----------vec --------------"<<endl;
   print_orbitals(cout , centers, basisidx, moCoeff);
#endif

   return 0;
}



int exist_file(string & filename)
{
   ifstream is(filename.c_str());
   if (is) { is.close(); return 1; }
   is.close(); return 0;
}


int string_upper(string &str)
{
    int i0; 
    for (int i=0;i<str.size();i++) {
	    i0=i;
	    // use toupper of C, not C++
	   if ( islower(str[i0]) ) str[i0] = (char)toupper(str[i0]); 
    }
    return 0;
}

#if 0
int find_atomname(vector<g03Atom> &atomlist, string& key)
{
   for (int i=0;i<atomlist.size(); i++) {
        if (atomlist[i].name==key) { return i;}
   }
   return -1;
}
#endif

int find_gatomname(vector<Atom> &atomlist, string& key)
{
   for (int i=0;i<atomlist.size(); i++) {
        if (strcasecmp(atomlist[i].name.c_str(),key.c_str())==0 ) { return i;}
   }
   return -1;
}



int atomname2number(string &atom)
{
   static char *t="atomname2number:";
   char const **s = element_lookup_caps_null; s++;
   int id=1;
   for ( ; *s; s++ ) {
       if ( strcasecmp(atom.c_str(), *s)==0 ) {
          return id;
       }
       id++;
   }
   cout <<t << " not found " << atom <<endl;
   exit(ERR_CODE);
   return 0;
}




int g03_process_sharpheader( 
        /* input */
          vector<string> &lines,
        /* output */
          int &found_gen, int &found_pseudo)
{
    const static char *t="g03_com_processheader:";

    int basis_6d=0, basis_10f=0;
    int iop3_24=0, iop3_18=0;
    int gvb_wf=0;

    int error_exit=0; 
    string comment;

    string all;
    for ( vector<string>::iterator it=lines.begin(); it!=lines.end(); it++) {
	    all += it->substr(1,it->size());
    }
    {
           string line_upper=all;
           string_upper(line_upper);
	   //cout << line_upper<<endl;
           if ( strstr(line_upper.c_str(),"/GEN")!=NULL ) {
                   found_gen=1;
                   //cout << "match /gen"<<endl;
		   comment += " /gen";
           }
           if ( strstr(line_upper.c_str(),"PSEUDO=READ")!=NULL ) {
                   found_pseudo=1;
                   //cout << "match pseudo" <<endl;
		   comment += " pseudo=read";
           }
	   if (strstr(line_upper.c_str(),"IOP(3/24=10)")!=NULL) {
               iop3_24=1;
	       comment += " IOP(3/24=10)";
	   }
	   if (strstr(line_upper.c_str(),"IOP(3/18=1)")!=NULL) {
	       iop3_18=1;
	       comment += " IOP(3/18=1)";
	   }
	   if (strstr(line_upper.c_str(),"GVB")!=NULL) {
		   gvb_wf = 1;
	   }
    }


    if (iop3_24==0 || iop3_18==0 ) {
        cout << " add 'IOP(3/24=10) IOP(3/18=1)' in the # section of the com file"<<endl;
	error_exit=1;
    }
    if (gvb_wf) {
        cout << " GVB is not supported" <<endl;
	error_exit=1;
    }
 
    if (error_exit) {
        exit(ERR_CODE);
    }
    return 0;
}

/* reorder orbital, gaussian -> gamess */
int g03_orbital_reorder(/* input */
//		vector<g03basis> & basisset,
		vector<Atom> &atoms,
		vector<Gaussian_basis_set> &gbasisset, 
		int dtype, int ftype, 
                   /* input & output */
		vector< vector<double> > & alphamo,
		vector< vector<double> > & betamo
		)
{
   const char *t="g03_orbital_reorder: ";
 /*  gamess f order
   43  MN 1 XXX   0.000000   3.702658   0.098126   0.000000   0.000000
   44  MN 1 YYY   0.000000  -0.098126   3.702658   0.000000   0.000000
   45  MN 1 ZZZ  -0.016783   0.000000   0.000000   3.843891   0.636826
   46  MN 1 XXY   0.000000  -0.043883   1.655879   0.000000   0.000000
   47  MN 1 XXZ  -0.011991   0.000000   0.000000   1.723675   0.287293
   48  MN 1 YYX   0.000000   1.655879   0.043883   0.000000   0.000000
   49  MN 1 YYZ  -0.011991   0.000000   0.000000   1.723675   0.287293
   50  MN 1 ZZX   0.000000   1.656851   0.043909   0.000000   0.000000
   51  MN 1 ZZY   0.000000  -0.043909   1.656851   0.000000   0.000000
   52  MN 1 XYZ   0.000000   0.000000   0.000000   0.000000   0.000000
*/
 /* gaussian  f order
  43       18XXX         .00000    .00000    .00000    .00213    .00000
  44       18YYY         .00000    .00000    .00213    .00000    .00000
  45       18ZZZ         .00000   -.00010    .00000    .00000    .00255
  46       18XYY         .00000    .00000    .00000    .00095    .00000
  47       18XXY         .00000    .00000    .00095    .00000    .00000
  48       18XXZ         .00000   -.00004    .00000    .00000    .00122
  49       18XZZ         .00000    .00000    .00000    .00097    .00000
  50       18YZZ         .00000    .00000    .00097    .00000    .00000
  51       18YYZ         .00000   -.00004    .00000    .00000    .00122
  52       18XYZ         .00000    .00000    .00000    .00000    .00000
*/
/*
 * gaussian xxx yyy zzz xyy xxy xxz xzz yzz yyz xyz
 *           1   2   3   4   5   6   7   8   9   10
 * gamess   xxx yyy zzz xxy xxz xyy yyz xzz yzz xyz
 *           1   2   3   5   6   4   9   7   8   10
 */
   const int n_f10order=10;
   const static int f10order[n_f10order]={ 1,  2,   3,   5,   6,   4,   9,   7,   8,   10 };

   // now support up to F */
   const int Lmaxorb=3; // max of L=F=3


   // find max of L
   // if L<="D" , do nothing
   int Lmax=0;
   int Lmax2=0;
   for (vector<Gaussian_basis_set>::iterator it=gbasisset.begin(); it != gbasisset.end(); it++) {
         for ( vector<string>::iterator itL = it->types.begin(); itL != it->types.end(); itL++ ) {
		 string Lstr = *itL;
		 int L= g03_maxangularMomentum(Lstr);
		 Lmax2 = max(Lmax2,L);
                 if ( L>Lmaxorb) {
			 cout << t << "L= " <<  L << " not supported"<<endl;
			 exit(ERR_CODE);
		 }
	 }
   }
   // error check, 
   Lmax = Lmax2;
   if ( Lmax==2 ) { // up to D
	   if (dtype==0) {
		   cout << "Error: now accepts only 6D" <<endl;
		   exit(ERR_CODE);
	   }
   } else if (Lmax==3) { // up to F
	   if (dtype==0 || ftype==0) {
		   cout << "Error: now accepts only 6D and 10F"<<endl;
		   exit(ERR_CODE);
	   }
   }

   if ( Lmax<3 ) { // S,P,D ,   no F
	   // reordering is uncecessary
	   // reordering is for F 
	   return 1;
   }

   // calc. vecsize
   int vecsize=0;
   int vecsize2=0;
   for ( vector<Atom>::iterator it=atoms.begin(); it!=atoms.end(); it++) {
	   int id=it->basis; 
	   for (vector<string>:: iterator itL = gbasisset[id].types.begin(); itL!= gbasisset[id].types.end(); itL++) {
		   string Lstr= *itL;
		   vecsize2+= g03_basissizeL(Lstr);
	   }
   }
   vecsize = vecsize2;


   // vecsize == alphamo[0].size()  ? 
   if ( vecsize != alphamo[0].size() ) {
      cout << t << "size of eigenvector != size of basisset" <<endl;
      cout << vecsize << " " << alphamo[0].size()<<endl;
      exit(ERR_CODE);
   }


   cout << t << "There exist F orbials"<<endl;

   // make order 
   vector<int> order;
   int itot=0;
   // use Gaussian_basis_set 
   for (vector<Atom>::iterator it=atoms.begin(); it!=atoms.end(); it++) {
	   int id=it->basis;
     for (vector<string>::iterator itL=gbasisset[id].types.begin(); itL!=gbasisset[id].types.end(); itL++) {
       string Lstr= *itL;
       if (Lstr=="F" || Lstr=="F10" ) {
	 int nbasisL=g03_basissizeL(Lstr);
	 for (int i=0; i < nbasisL  ; i++) {
	   order.push_back(itot+f10order[i]-1);
	 }
	 itot += nbasisL;
       }
       else {
	 int nbasisL=g03_basissizeL(Lstr);
	 for (int i=0; i < nbasisL  ; i++) {
	   order.push_back(itot+i);
	 }
	 itot += nbasisL;
       }
     }
   }


   // check size again
   if ( order.size() != vecsize) {
	   cout << t << "internal error, size differs"<<endl;
	   cout << order.size() << " " << vecsize <<endl;
	   exit(ERR_CODE);
   }

#if 0
    // print order
    cout << t << "reorder eigenvector, order="<<endl;
    for (vector<int>::iterator it = order.begin(); it!=order.end(); it++ ) {
	cout << *it << " " ;
    }
    cout << endl;
#endif

   // reorder vector
   //
   // alphamo

   for ( vector< vector<double> >::iterator it=alphamo.begin(); it!= alphamo.end(); it++ ) {
       vector<double> v;
       v = *it; 
       for (int i=0;i<vecsize;i++) {
	  if (i!=order[i]) {
             (*it)[i] = v[ order[i] ];
          }
       }
   }

   // betamo
   if ( ! betamo.empty() ) {

   for ( vector< vector<double> >::iterator it=betamo.begin(); it!= betamo.end(); it++ ) {
       vector<double> v;
       v = *it;
       for (int i=0;i<vecsize;i++) {
          if (i!=order[i]) {
             (*it)[i] = v[ order[i] ];
          }
       }
   }

   }

   return 0;
}



int basis_uniq_id(/* input */
                int bas,  vector<Atom> atoms)
{
    const static char *t="is_uniq_basis: ";

    for (int id=0; id<bas; id++) {
       if ( atoms[id].name == atoms[bas].name ) { return id; }
    }
    return bas; /* uniq */
}


 int lastkey_is_spc(/* input */ string &str)
{
   if ( str.find_first_not_of(" ")==string::npos )  return 1;
   return 0;
}

 int lastkey_is_spc_substr20(/* input */ string &str)
{
   string sub=str.substr(0,20);
//   cout << "<"<< sub<<">"<<endl;
   if ( sub.find_first_not_of(" ")==string::npos )  return 1;
   return 0;
}

/* return found=1, not found=0 */
int filekeyinput_get_string(
        /* input */
       string & filename,
       string & keyin,
       /* output */
       string & str)
{
   const char *t="filekeyinput_get_string: ";
   ifstream is(filename.c_str());
  if(!is) {
    cout <<t<< "Couldn't open " << filename << endl;
    exit(ERR_CODE);
  }

   string line;
   while ( getline(is,line) ) {
       if ( line.find(keyin)!=string::npos ) {
          str=line;
          return 1;
       }
   }
   return 0;
} 

/* return found=1, not found=0 */
int filekeyinput_get_stringv(
       /* input */
       string & filename, 
       string & keyin, int nskip, 
       int nread, string & keyout, int  (*func)(string &), 
       /* output */
       vector<string> & str )
{
   const char *t="filekeyinput_get_stringv: ";
   ifstream is(filename.c_str());

  if(!is) {
    cout <<t<< "Couldn't open " << filename << endl;
    exit(ERR_CODE);
  }

   string line;

   str.clear();
   while ( getline(is,line) ) {
       if ( line.find(keyin)!=string::npos ) {
 //         cout << "found "<< keyin <<endl;
          for (int i=0;i<nskip;i++) getline(is,line); 
	  if (nskip==-1) str.push_back(line); 
          if (nread>0) {
              for (int i=0;i<nread;i++) {
                  getline(is,line);
                  str.push_back(line);
              }
              goto foundandlast;
          }
          else if ( keyout!="" ) {
              while ( getline(is,line) ) {
                  if ( line.find(keyout) != string::npos ) goto foundandlast;
                  str.push_back(line);
              }
          }
          else {
              while ( getline(is,line) ) {
                  if ( func(line) ) goto foundandlast;
                  str.push_back(line);
              }
          }
       }
   }

   return 0;
 foundandlast: 
   /* found */
   return 1;
}


/* return found=1, not found=0 */
int filekeyinput_get_stringv_last(
       /* input */
       string & filename, 
       string & keyin, int nskip, 
       int nread, string & keyout, int  (*func)(string &), 
       /* output */
       vector<string> & str )
{
   int found_once=0;
   const char *t="filekeyinput_get_stringv_last: ";
   ifstream is(filename.c_str());

  if(!is) {
    cout <<t<< "Couldn't open " << filename << endl;
    exit(ERR_CODE);
  }

   string line;

   str.clear();
   while ( getline(is,line) ) {
       if ( line.find(keyin)!=string::npos ) {
 //         cout << "found "<< keyin <<endl;
          for (int i=0;i<nskip;i++) getline(is,line); 
          if (nread>0) {
              str.clear();
              for (int i=0;i<nread;i++) {
                  getline(is,line);
                  str.push_back(line);
              }
	      found_once=1;
              goto foundandlast;
          }
          else if ( keyout!="" ) {
              str.clear();
              while ( getline(is,line) ) {
                  if ( line.find(keyout) != string::npos ) { 
		        found_once=1;
		        goto foundandlast;
		  }
                  str.push_back(line);
              }
          }
          else {
              str.clear();
              while ( getline(is,line) ) {
                  if ( func(line) ) {
                        found_once=1;
                        goto foundandlast;
                  }
                  str.push_back(line);
              }
          }
       }
foundandlast: ;
   }

   if (found_once) { 
       return 1;
   }
   else { 
       return 0;
   }
}


#if 0
main()
{
   string filename="Mn2O2A_S11_gen.log";

   vector<string> strv;
   string keyin;
   string keyout;
   int nskip, nread;

   keyin="Standard orientation:";
   nskip=4; nread=4; keyout="";
   file_get_stringv(filename, keyin, 
      nskip, 
      nread, keyout, NULL, 
      strv);
   
   for (vector<string>::iterator itstrv=strv.begin();
            itstrv != strv.end();  itstrv++ ) {
      cout << *itstrv <<endl; 
   }


   keyin="AO basis set in the form of";
   nskip=0; nread=0; keyout=""; 
   file_get_stringv(filename, keyin,
      nskip,
      nread, keyout, lastkey_is_spc,
      strv);

   for (vector<string>::iterator itstrv=strv.begin();
            itstrv != strv.end();  itstrv++ ) {
      cout << *itstrv <<endl;
   }


   keyin="Pseudopotential Parameters";
   nskip=4; nread=0; keyout="=================";
   file_get_stringv(filename, keyin,
      nskip,
      nread, keyout, NULL,
      strv);

   for (vector<string>::iterator itstrv=strv.begin();
            itstrv != strv.end();  itstrv++ ) {
      cout << *itstrv <<endl;
   }


   keyin="EIGENVALUES --";
   nskip=0; nread=0; keyout="";
   file_get_stringv(filename, keyin,
      nskip,
      nread, keyout, lastkey_is_spc_substr20,
      strv);

   for (vector<string>::iterator itstrv=strv.begin();
            itstrv != strv.end();  itstrv++ ) {
      cout << *itstrv <<endl;
   }


   keyin="Standard basis:";
   string str;
   file_get_string(
        /* input */
        filename,
        keyin,
       /* output */
        str); 
   cout << str<<endl;

}

#endif


int g03_com_processatomgeometry(/* input */
		vector<string> lines,
		/* output */
	//	vector<g03Atom> &g03atoms,
		vector<Atom> &atomscom)
{
   const static char *t="g03_com_processatomgeometry:";
   string spc=" ";
   vector<string> words;
//   g03atoms.clear();
   atomscom.clear();
   int id=0;


   for ( vector<string>::iterator it=lines.begin(); it!=lines.end(); it++) {
	   if ((*it)[0]==' ') break;
	   words.clear(); split(*it,spc,words);
	   Atom a_atom;
	   if (words[0]=="X") continue;
	   a_atom.name=words[0];
	   atomscom.push_back(a_atom);
   }
   //cout << t << "number of atoms = "<< g03atoms.size() << endl;
   return 0;
}


int g03_com_processPS( /* input */
          vector<string> lines,
	 // vector<g03Atom> &g03atoms,
          vector<Atom> &atoms ,
          /* output */
          // vector<g03PSP> &g03ps,
	  vector<Gaussian_pseudo_writer> &pseudos
	  )
{
   const static char *t="g03_com_processPS:";


   vector<string> words;
   string spc=" ";
   vector<string>::iterator it=lines.begin(); 
   while (1) {
//    g03PSP atom_pseudo;

    Gaussian_pseudo_writer a_pseudo; 

    // next must be atom and 0
    words.clear(); split(*it,spc,words);
    string atomname=words[0];
 //   atom_pseudo.atomnum = atomname2number(atomname);
    string_upper(atomname);
    a_pseudo.label = atomname; 

    int id = find_gatomname(atoms,atomname)+1;
    if (id>0) { 
      a_pseudo.atomnum=id;
    }
    else {
      cout << t << "atomname '" << atomname << "' not found" <<endl;
      exit(ERR_CODE);
    }


    
    if ( words[1] != "0" ) {
      cout << t << "I don't understand "<< *it <<endl;
      exit(ERR_CODE);
    }
    // ecpname, lmax, core
    it++;
    words.clear(); split(*it,spc,words);
    int Lmax = atoi(words[1].c_str());
    //cout << "Lmax= " << Lmax <<endl;
    for (int L0=0;L0<=Lmax; L0++) {
//      g03PSP_L pseudoL;
      vector<int> powlist;
      vector<double> explist;
      vector<double> coeflist;
      it++; // "D", "S - D",..., angular momentum
      int L, refL;
      //cout << *it <<endl;
      L = label_to_L( *it, refL);
      //cout << "L,refL= "<< L << " " << refL <<endl;
//      pseudoL.L=L;
//      pseudoL.refL=refL;
//      pseudoL.power.clear(); pseudoL.expn.clear(); pseudoL.coeff.clear();

      it++; // number of gaussian
      words.clear(); split(*it,spc,words);
      int ng = atoi(words[0].c_str());
      //cout << "L,ng= " << L << " " << ng <<endl;

      for (int ig=0;ig<ng;ig++) {
	it++;
	words.clear(); split(*it,spc,words);
	int pow= atoi(words[0].c_str());
	double expn=atod(words[1].c_str());
	double coeff=atod(words[2].c_str());
//	pseudoL.power.push_back(pow);
//	pseudoL.expn.push_back(expn);
//	pseudoL.coeff.push_back(coeff);

	powlist.push_back( pow -2 );
	explist.push_back( expn );
	coeflist.push_back( coeff );
      }
 //     atom_pseudo.psL.push_back(pseudoL);

      a_pseudo.nvalue.push_back(powlist);
      a_pseudo.exponents.push_back(explist);
      a_pseudo.coefficients.push_back(coeflist);
    } // for

 //   g03ps.push_back(atom_pseudo);

    pseudos.push_back(a_pseudo);
    a_pseudo.nvalue.clear();
    a_pseudo.exponents.clear();
    a_pseudo.coefficients.clear();

    *it++;
    //cout << " next="<<line <<endl;
    if ( it==lines.end() ) {
      break;
    }

  } // while

#if 0
   cout << "pseudo "<<endl;
   for (vector<Gaussian_pseudo_writer>::iterator it=pseudos.begin(); it!=pseudos.end(); it++) {
        it->print_pseudo(cout);
   }
#endif

   return 0;
}





/* read fchk header format */
static int g03_fchk_read_header(/* input */string &line,string & key,
		        /* output */
			      string &type, string& val, 
		      /* input */int format)
{
  const char *t ="g03_fchk_read_header: ";
  string space=" ";
  vector<string> words;
  
  words.clear();
  split(key,space,words); 
  int nkeywords=  words.size();
  words.clear();
  split(line,space,words);
  if (format==1) {
    if ( nkeywords+2 != words.size() ) {
      cout << t << "keywords length error,"
	   << nkeywords << " " << words.size() <<endl;
      exit(ERR_CODE);
    }
    type=words[nkeywords];           
    val =words[nkeywords+1];
  }
  else {

    if ( nkeywords+3 != words.size() ) {
      cout << t<< "keywords length error,"
	   << nkeywords << " " << words.size() <<endl;
      exit(ERR_CODE);
    }
    if (words[nkeywords+1]!="N=") {
      cout << t<< "format error "<<endl;
      exit(ERR_CODE);
    }
    type=words[nkeywords];
    val =words[nkeywords+2];


  }
  return 0;
}


// -------------------------------------------------------
// read g03 fchk file
// -------------------------------------------------------
//
// using interface might be better
// 
/* read fchk,  integer format */
 int g03_fchk_read_i(/* input */ string & filename,string & key, 
		 /* output */ int &i)
{
   const char *t =  "g03_fchk_read_i: ";
  ifstream  is(filename.c_str());

  if(!is) {
    cout << t<< "Couldn't open " << filename << endl;
    exit(ERR_CODE);
  }

  i=0;

  string line;
  string space=" ";
  vector<string> words;
  
  while (getline(is,line)) {
    if ( line.find(key)!= string::npos ) {
      string type,val;
      g03_fchk_read_header(line,key,type,val,1);
      if (type != "I" ) {
	cout << t<< "header error,type mismatch"<<endl;
	exit(ERR_CODE);
      }
      i = atoi( val.c_str() );
      is.close();
      return 0;
    }
  }
  is.close();
  return 1;
}

int g03_fchk_read_r(/* input */ string & filename,string & key,
                 /* output */ double &r)
{
   const char *t =  "g03_fchk_read_r: ";
  ifstream  is(filename.c_str());

  if(!is) {
    cout << t<< "Couldn't open " << filename << endl;
    exit(ERR_CODE);
  }

  r=0.0;

  string line;
  string space=" ";
  vector<string> words;

  while (getline(is,line)) {
    if ( line.find(key)!= string::npos ) {
      string type,val;
      g03_fchk_read_header(line,key,type,val,1);
      if (type != "R" ) {
        cout << t<< "header error,type mismatch"<<endl;
        exit(ERR_CODE);
      }
      r = atod( val.c_str() );
      is.close();
      return 0;
    }
  }
  is.close();
  return 1;
}


/* read fchk, vector of integer format */
 int g03_fchk_read_iv(/* input */ string & filename,string & key, 
			  int n, 
			  /* output */int & nout,vector<int> &iv)
{
  const char *t = "g03_fchk_read_iv: ";
  ifstream  is(filename.c_str());

  if(!is) {
    cout <<t << "Couldn't open " << filename << endl;
    exit(ERR_CODE);
  }

  iv.clear();
  nout=0;

  string line;
  string space=" ";
  vector<string> words;

  while (getline(is,line)) {
    if ( line.find(key)!= string::npos ) {
      string type,val;
      g03_fchk_read_header(line,key,type,val,2);
      if (type != "I" ) {
	cout <<t<< "header error,type mismatch"<<endl;
	exit(ERR_CODE);
      }
      nout = atoi( val.c_str() );
      if ( n!=0 && (n != nout) ) {
	cout <<t<< "size mismatch " << n << " " <<nout << endl;
      }
      iv.clear();
      for (int i=0;i<nout;i++) { 
	int val;
	is >>  val;
	iv.push_back(val);
      }
      is.close();
      return 0;
    }
  }
  is.close();
  return 1;
}

/* read fchk, vector of double format */
 int g03_fchk_read_rv(/* input */string & filename,string & key,
			  int n, 
			  /* output */int & nout,vector<double> &dv)
{
   const char *t = "g03_fchk_read_rv: ";
  ifstream  is(filename.c_str());

  dv.clear();
  nout=0;

  if(!is) {
    cout <<t<< "Couldn't open " << filename << endl;
    exit(ERR_CODE);
  }

  string line;
  string space=" ";
  vector<string> words;

  while (getline(is,line)) {
    if ( line.find(key)!= string::npos ) {
      string type,val;
      g03_fchk_read_header(line,key,type,val,2);
      if (type != "R" ) {
	cout <<t<< "header error,type mismatch"<<endl;
	exit(ERR_CODE);
      }
      nout = atoi( val.c_str() );
      if ( n!=0 && (n != nout) ) {
	cout <<t<< "size mismatch " << n << " " <<nout << endl;
      }
      dv.clear();
      for (int i=0;i<nout;i++) {
	double val;
	is >>  val;
	dv.push_back(val);
      }
      is.close();
      return 0;
    }
  }
  is.close();
  return 1;
}



/* g03fchk shell type (integer) -> shell type (char *) */
static const char *g03fchk_shelltype2orbtype(int i, int dtype, int ftype)
{
   const char *t="shelltype2orbtype";

   const static char *SPstr="SP";

   if (i==-1) return SPstr;
   if (i<0 || i>=nlabel_L ) {
      cout << t << " error in shell id " << i << endl;
      exit(ERR_CODE);
   }
   char const *ret= label_Lstr[i];
   // 0<=dtype<=1
   // 0<=ftype<=1
   // error ?
   if ( dtype<0 || dtype>1  || ftype<0 || ftype>1 ) {
       cout << t << " erro in dtype or ftype" <<endl
	       << dtype << " " << ftype <<endl;
       exit(ERR_CODE);
   }
   if (i==2 && dtype==0) {// 5D
	   ret=dtypestr[dtype];
   }
   if (i==3 && ftype==0) { // 7F
           ret=ftypestr[ftype];
   }
   return ret; 
}

int g03_basisset_delete_dup(
		/*input and output*/
		vector<Atom> &atoms,
		vector<Gaussian_basis_set> &basisset)
{
	const char *t="g03_basisset_delete_dup:";
	vector<string> namesdone;

     vector<Gaussian_basis_set> basissetnew; 

     for (int iatom=0;iatom<basisset.size(); iatom++) {
	     int found=0;
	     string name = basisset[iatom].label;
	     for (int inamesdone=0;inamesdone<namesdone.size(); inamesdone++) {
                   if (name==namesdone[inamesdone]) { 
			   found=1;
			   break;
		   }
	     }
	     if (found==0) {
                   basissetnew.push_back(basisset[iatom]);
		   namesdone.push_back(name);
	     }
     }
     basisset.clear();
     basisset=basissetnew; 
  
  
     // fix atoms.basis
     for (int iatom=0;iatom<atoms.size(); iatom++) {
	     string name=atoms[iatom].name;
	     int found=0;
         for (int ibasis=0;ibasis<basisset.size(); ibasis++) {
		 if (  name== basisset[ibasis].label ) {
			 atoms[iatom].basis=ibasis;
			 found=1;
			 break; 
		 }
	 }
	 if (found==0) {
             cout << t <<  "can not find atom name="<<name<<endl;
	     exit(ERR_CODE);
	 }
     }
  

     return 0;
}

/* read fchkfileame,  set basisset */
int g03_fchk_readbasisset(/* input */ string & fchkfilename,
		vector<Atom> atoms, 
		                /* output */
	//	        vector<g03basis> & basisset ,
			vector<Gaussian_basis_set> &gbasisset
			  )
{
   const char *t="g03_fchk_readps";
   int idummy=0;
   int nout;
   string key;
  
//   basisset.clear();

   int dtype, ftype; 
   key="Pure/Cartesian d shells";
   g03_fchk_read_i(fchkfilename,key,dtype);
   key="Pure/Cartesian f shells";
   g03_fchk_read_i(fchkfilename,key,ftype);

   int n_shelltypes;
   vector<int> shelltypes;
    key="Shell types";
   g03_fchk_read_iv(fchkfilename,key, idummy, n_shelltypes, shelltypes);

   int n_noprimitives;
   vector<int> noprimitives;
    key="Number of primitives per shell";
   g03_fchk_read_iv(fchkfilename,key, idummy, n_noprimitives, noprimitives);

   int n_shellmap;
   vector<int> shellmap;
    key="Shell to atom map";
   g03_fchk_read_iv(fchkfilename,key, idummy, n_shellmap, shellmap);

   int n_primitive_exp;
   vector<double> primitive_exp;
    key="Primitive exponents";
   g03_fchk_read_rv(fchkfilename,key, idummy, n_primitive_exp, primitive_exp);

   int n_contract_coeff;
   vector<double> contract_coeff;
    key="Contraction coefficients";
   g03_fchk_read_rv(fchkfilename,key, idummy, n_contract_coeff, contract_coeff );

   int n_p_contract_coeff=0;
   vector<double> p_contract_coeff;
    key="P(S=P) Contraction coefficients";
   g03_fchk_read_rv(fchkfilename,key, idummy, n_p_contract_coeff, p_contract_coeff );


  // size check
  if ( n_shelltypes != n_noprimitives || n_noprimitives != n_shellmap ) {
     cout << t<< " error in shell number " << n_shelltypes << " " << n_noprimitives << " " 
         << n_shellmap <<endl;
     exit(ERR_CODE);
  }
#if 0
  if ( n_primitive_exp != n_contract_coeff || n_contract_coeff != n_p_contract_coeff ) {
     cout << t<< " error in primitive exp " <<  n_primitive_exp << " " << n_contract_coeff  << " "
          << n_p_contract_coeff <<endl; 
     exit(ERR_CODE);
  }
#else
  if ( n_primitive_exp != n_contract_coeff ) {
     cout << t<< " error in primitive exp " <<  n_primitive_exp << " " << n_contract_coeff <<endl;
     exit(ERR_CODE);
  }
  if (n_p_contract_coeff>0) {
    if (n_contract_coeff != n_p_contract_coeff ) {
	    cout << t<< " error in P(S=P) Contraction coefficients " << n_contract_coeff << " " 
		    << n_p_contract_coeff <<endl;
    }
  }

#endif

  // make basisset
  Gaussian_basis_set a_gbasis; 
  int natom=0;
  for ( int ishell=0; ishell< n_shellmap; ishell++ ) {
	  natom = max(natom, shellmap[ishell]);
  }
  gbasisset.clear();
  gbasisset.resize(natom);
  int iprimitive=0;
  for ( int ishell=0; ishell< n_shellmap; ishell++ ) {
	  int iatom=shellmap[ishell]-1; 
	  string type= g03fchk_shelltype2orbtype(shelltypes[ishell],dtype,ftype );
	  vector<double> explist;
	  vector<double> coeflist;
	  int ncoeff= noprimitives[ishell]; 
          for ( int ip_shell =0 ; ip_shell < ncoeff; ip_shell++ ) {
               explist.push_back(primitive_exp[iprimitive+ip_shell]);
	       coeflist.push_back( contract_coeff[iprimitive+ip_shell]);
	  }
	  if (shelltypes[ishell]>=0) {
	     gbasisset[iatom].types.push_back(type);
	  }
	  else {
	     gbasisset[iatom].types.push_back("S");
	  }
	  gbasisset[iatom].exponents.push_back(explist);
	  gbasisset[iatom].coefficients.push_back(coeflist); 
	  if ( shelltypes[ishell] == -1 && n_p_contract_coeff>0  ) {
		  explist.clear();
		  coeflist.clear();
            for ( int ip_shell =0 ; ip_shell < ncoeff; ip_shell++ ) {
               explist.push_back(primitive_exp[iprimitive+ip_shell]);
               coeflist.push_back( p_contract_coeff[iprimitive+ip_shell]);
	    }
	    if (shelltypes[ishell]<0) {
		    gbasisset[iatom].types.push_back("P");
	    }
	    gbasisset[iatom].exponents.push_back(explist);
	    gbasisset[iatom].coefficients.push_back(coeflist);
	  }
	  iprimitive += ncoeff; 
  }
  // fix label
  for (int iatom=0;iatom<gbasisset.size(); iatom++) {
      gbasisset[iatom].label= atoms[iatom].name;
      string_upper(gbasisset[iatom].label);
  }
#if 0
  for ( vector<Gaussian_basis_set>::iterator it=gbasisset.begin(); 
		  it!=gbasisset.end(); it++) {
	  cout << "basis "<<endl;
	  it->print_basis(cout);
  }
#endif


  return 0;
}



/* read logfile, get orbital names,
 * orbital names are S,X,Y,Z, XX,YY, or  D-2,D-1,... 
 */
int g03_log_analyze_orbname(
               /* input */
		string & logfilename , 
	       /* output */
               vector<string> orbname)
{
  const char *t="g03_log_analyze_orbname";

  string keyin,keyout;
  int nskip, nread;
  vector<string> strv;

   int ret;

   keyin="EIGENVALUES --";
   nskip=0; nread=0; keyout="";/* use lastkey_is_spc_substr20 */
   ret=filekeyinput_get_stringv(logfilename, keyin, nskip,
      nread, keyout, lastkey_is_spc_substr20,
      strv);

   string name;
   orbname.clear();
   for (vector<string>::iterator itstrv=strv.begin();
            itstrv != strv.end();  itstrv++ ) {
//          cout << *itstrv <<endl;
          name = itstrv->substr(11,6);
          orbname.push_back(name);
   }

last:

  int nbasis = orbname.size();

  // error check 
  { 
	  int flag=0; 
	  for (int i=0;i<nbasis;i++) {
	    if ( orbname[i].find("D 0") != string::npos ) { 
	      flag=1; 
	      cout << "It seems that this is 5d orbital." <<endl;
	    }
	    if (orbname[i].find("F 0") != string::npos) {
	      flag=1;
	      cout << "It seems that this is 7f orbital." <<endl;
	    }
	  }
	  if (flag) {
	    cout <<t <<" error orbital\n" <<endl;
	    exit(ERR_CODE);
	  }
  }
  return 0;
}




static int basisset_find_atomlast(vector<string> & basissetstring,
		int startline )
{
   for (int i=startline; i<basissetstring.size(); i++) {
	   string str=basissetstring[i].substr(0,5); 
	   if ( str==" ****" ) { return i; }
   }
   return basissetstring.size()-1; 
}

int g03log_analyze_basisset(/*input*/
		vector<string> & basissetstring,
		vector<Atom> &atoms, 
		/*output*/
		vector<Gaussian_basis_set> &basisset) 
{

  const char *t="g03log_analyze_basisset:";
  Gaussian_basis_set a_basisset; 
  int iline=0;
  int nline=basissetstring.size(); 
  //cout << "nline: " << nline << endl;
  string sep=" ";
  vector<string> words;
  while (iline < nline) {
      /* 1 0  <- atomid */ 
      words.clear(); split(basissetstring[iline], sep,words);
      int atomnum=atoi(words[0].c_str()); 
     // cout << "atomnum=" << atomnum <<endl;
      iline++;
      vector<double> explist;
      vector<double> coeflist;
      vector<double> coeflist2;
      int lastatom=basisset_find_atomlast(basissetstring,iline);
      for ( ; iline <= lastatom; ) {
          words.clear(); split(basissetstring[iline], sep,words);
	  int spflag;
	  if ( words[0]=="SP" || words[0]=="L") {
		  spflag=1;
		  a_basisset.types.push_back("S");
	  }
	  else {
		  spflag=0;
		  a_basisset.types.push_back(words[0]);
	  }
	  int n=atoi(words[1].c_str());
	  //cout << "type,n=" << words[0] << " " << n <<endl;
          iline++;
	  explist.clear(); coeflist.clear(); coeflist2.clear();
	  for (int i=0;i<n;i++) {
              words.clear(); split(basissetstring[iline], sep,words);
	      double val=atod(words[0].c_str());
	      explist.push_back(val);
	      val = atod(words[1].c_str());
	      coeflist.push_back(val);
	      if ( spflag ) {
	      val = atod(words[2].c_str());
	      coeflist2.push_back(val);
	      }
	      iline++;
	  }
	  a_basisset.exponents.push_back(explist); 
	  a_basisset.coefficients.push_back(coeflist); 
          if ( spflag ) {
		a_basisset.types.push_back("P");
		a_basisset.exponents.push_back(explist);
		a_basisset.coefficients.push_back(coeflist2);
	  }

	  string str=basissetstring[iline].substr(0,5);
	  if (str==" ****") iline++; 
      }
      basisset.push_back(a_basisset); 
      a_basisset.types.clear();
      a_basisset.exponents.clear();
      a_basisset.coefficients.clear();
  }
   // fix label
  if ( basisset.size() != atoms.size() ) {
     cout << t << "basisset.size() != atoms.size()"<<endl;
     exit(ERR_CODE);
  }
  for (int iatom=0;iatom<basisset.size(); iatom++) {
    //cout << "basis " << iatom << " name " << atoms[iatom].name << endl ;
     basisset[iatom].label= atoms[iatom].name; 
     string_upper(basisset[iatom].label);
  }
#if 0
  /*printout*/
  for (vector<Gaussian_basis_set>::iterator it=basisset.begin(); 
		  it != basisset.end(); it++) {
	  it->print_basis(cout);
  }
#endif
}


/*
          1         2         3        4          5         6         7         8         9
0123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789

  Center     Atomic      Valence      Angular      Power                                                       Coordinates
    1         25           15                                                                      -1.534430 -2.555963  -.017
0123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789
                                      D and up
                                                     2      138.2855530     6414.98369000
*/

static int psstr_find_atomrange(/*input*/
     vector<string> psstr,
     int start) 
{
   for (int i = start+1; i< psstr.size(); i++) { 
     string str =  psstr[i].substr( 0 ,10); /* find Center */
     if (strlen_trim(str.c_str())>0) { return i-1; }
   }
   return psstr.size()-1; 
}

static int psstr_find_Lrange(/*input*/
     vector<string> psstr,
     int start)
{
   for (int i = start+1; i< psstr.size(); i++) {
     string str =  psstr[i].substr(38,10);  /* find Angular */
     if (strlen_trim(str.c_str())>0) { return i-1; }
     str =  psstr[i].substr( 0 ,10);   /* find Center */
     if (strlen_trim(str.c_str())>0) { return i-1; }
   }
   return psstr.size()-1;
}

int g03log_analyze_ps(/*input*/
                vector<string> & psstr,
                /* output */
        vector<Gaussian_pseudo_writer> & pseudos   ,
	vector<PseudoValence> &pseudovalence )
{
  pseudos.clear();

  string spc=" ";
  vector<string> words;  

  int centernum, atomnum, valencenum;
  Gaussian_pseudo_writer a_pseudo; 
  vector<int> nlist;
  vector<double> explist;
  vector<double> coeflist;
  int iline=0; 
  while ( iline<psstr.size() ) {
      string str =  psstr[iline].substr( 0 ,10); 
      centernum=atoi(str.c_str()); 
      str =  psstr[iline].substr( 10 ,10);
      atomnum=atoi(str.c_str());
      str =  psstr[iline].substr( 20 ,10);
      valencenum=atoi(str.c_str());
      iline++;
      int lastlineatom = psstr_find_atomrange(psstr, iline);
      a_pseudo.nvalue.clear();
      a_pseudo.exponents.clear();
      a_pseudo.coefficients.clear();
      vector<int> Lnumlist;
      vector<int> refLlist; 
      int addit;
      while ( iline <= lastlineatom ) {
         addit=0;
	 if ( psstr[iline].find("No pseudopotential") != string::npos  ) {
		 iline=lastlineatom+1;
		 break;
	 }
	 if ( psstr[iline].find("Pseudopotential same as on center") != string::npos  ) {
		 iline=lastlineatom+1;
		 break;
	 }
	 addit=1;
         str =  psstr[iline].substr(38,10);  /* find Angular */

         int refL; 
         int Lnum = label_to_L(str,refL); 
	 Lnumlist.push_back(Lnum);
	 refLlist.push_back(refL);

         nlist.clear(); explist.clear(); coeflist.clear();
         iline++;
         int lastlineL = psstr_find_Lrange(psstr,iline);
	 //cout << "from " << iline << " to " << lastlineL<<endl;
         for ( ; iline<=lastlineL; iline++) {
            words.clear(); split(psstr[iline], spc, words);            
            int num=atoi(words[0].c_str());
            nlist.push_back(  num -2 );
            double expval=atod(words[1].c_str());
            explist.push_back( expval );
            double coefval; 
            if (words[2]=="**************") { coefval=0.0; }
            else { coefval=atod(words[2].c_str()); }
            coeflist.push_back( coefval );
         }
         a_pseudo.nvalue.push_back(nlist);
         a_pseudo.exponents.push_back(explist);
         a_pseudo.coefficients.push_back(coeflist);
         
      }  /* iline< lastlineatom  */
      if (addit) {
      a_pseudo.label = element_lookup_caps_null[atomnum]; 
      a_pseudo.atomnum = atomnum; 
      pseudos.push_back(a_pseudo); 

      PseudoValence pv;
      pv.label = element_lookup_caps_null[atomnum];
      pv.atomnum=atomnum;
      pv.valencenum=valencenum;
      pv.Lnum=Lnumlist;
      pv.refL=refLlist;
      pseudovalence.push_back(pv);
      }
nextloop: ;
  } /* iline < psstring.size() */

#if 0
  /* printout */
  int npsp=pseudos.size();
  cout << "npsp="<<npsp<<endl;
  for(int psp=0; psp < npsp; psp++) {
    cout << "---------------" << psp <<"-------------" << endl; 
    pseudos[psp].print_pseudo(cout);
  }
#endif

}

int g03_log_readatoms( /*input*/ vector<string> & strv,
		vector<Atom> &atomslog )
{
        const char *t="main_log2atoms: ";
       // atoms.clear();
	atomslog.clear();
  string space=" ";
  vector<string> words;
  int itot=0;
  for (vector<string>::iterator it = strv.begin(); it!=strv.end(); it++ ) {
          itot++;
       words.clear();
       split(*it,space,words);
   //    g03Atom a_atom;
       Atom atom;
       /* 0 id,  1 atomicnumber, 2 atomictype=0, 3-5 positon */
       if ( words.size()!= 6) { cout << t << " error, size != 6 "<<endl; exit(ERR_CODE); }
       int id=atoi(words[0].c_str());
     //  a_atom.id=id;
       if (itot!=id) {  cout << t << " error, wrong id?  "<<endl; exit(ERR_CODE); }
       int atomnum = atoi(words[1].c_str());
       if (atomnum==-2) { continue; } // -2 is translational vector 
      // a_atom.name= element_lookup_caps_null[atomnum];
      // a_atom.atomnum=atomnum;
       atom.name=element_lookup_caps_null[atomnum];
       string_upper(atom.name); 
       atom.charge=atomnum; 
      // a_atom.pos.clear();
       for (int i=3;i<=5;i++) {
          double r = atod(words[i].c_str()); // unit is  Ang 
       //   a_atom.pos.push_back(r);
	  atom.pos[i-3]=r/ANG ;    // ANG -> AU 
       }
   //    atoms.push_back(a_atom);
       atomslog.push_back(atom); 
   }
  return 0;
}



/* write VEC section of gamess */

static int  g03gamess_vectorsection( /* input */
                int nbasis, 
		vector< vector<double> > & alphamo, 
		vector< vector<double> > & betamo ,
                ofstream & os 
		/* no output */
		)
{
   const char *t="g03gamess_vectorsection";

   char buf[30];

   cout << t << " size of vector=" << nbasis <<  endl;

  os << " $VEC" <<endl;

  for (int i=0; i< nbasis ; i++ ) {
     int k=0;
     int row = (nbasis+1)/5;
     if ( (nbasis+1)%5 != 0 )  { row++; }
     for (int j=0; j< row ; j++ ) {
         if (k>=nbasis) break;
         sprintf(buf,"%2d %2d",(i+1)%100, (j+1)%100);
         os << buf;
         for (int m=0;m<5;m++) {
	    if (k>=nbasis) break;
             sprintf(buf,"%15.8E",alphamo[i][k]);
             os << buf ;
             k++;
         }
         os <<endl;
     }
  }


  if ( betamo.size()>0) {

  for (int i=0; i< nbasis ; i++ ) {
     int k=0;	   
     int row = (nbasis+1)/5;
     if ( (nbasis+1)%5 != 0 )  { row++; }
     for (int j=0; j< row ; j++ ) {
         if (k>=nbasis) break;
	 sprintf(buf,"%2d %2d",(i+1)%100, (j+1)%100);
         os << buf;
         for (int m=0;m<5;m++) {
	     if (k>=nbasis) break;
             sprintf(buf,"%15.8E",betamo[i][k]);
             os << buf ; 
             k++;
         }
	 os <<endl;
     }
  }
  }


  os << " $END" <<endl;

  return 0;
}


/* write header section of gamess */
static int  g03gamess_headersection(/* input */
		int multiplicity, string & scftype,int nbasis,
                vector<Gaussian_pseudo_writer> & ecp, 
                ofstream & os )
{
   const char *t ="g03gamess_datasection";

  string ecpread;
  if ( ecp.size() >0 ) {  ecpread="ECP=READ" ; }
  else { ecpread=""; }

  os << " $CONTRL SCFTYP=" << scftype << " MULT=" << multiplicity 
     << " ICHARG=0  RUNTYP=ENERGY MAXIT=70 " << ecpread << " $END" <<endl
     << " $SYSTEM KDIAG=3 $END" <<endl
     << " $DFT DFTTYP=B3LYP NRAD0=48 NTHE0=12 NPHI0=24"<< endl
     << "    SWITCH=1.0E-4 $END" << endl
     << " $SCF DIRSCF=.T. RSTRCT=.F. SOSCF=.F.  NPUNCH=2 $END" <<endl
     << " $GUESS GUESS=MOREAD NORB=" << nbasis << " $END" <<endl<<endl;

  cout << "I add '$DFT DFTTYPE=B3LYP ...' section, check the input file"<<endl;

  return 0;
}

/* write data section of gamess */
static int g03gamess_datasection(   /* input */
                vector<Atom> & atoms,
                vector<Gaussian_basis_set>  & basisset,
                ofstream & os
                )
{
        const char *t="g03gamess_datasection";
        char buf1[30], buf2[30], buf3[30];
        const char *form="%16.9E";

 os << " $DATA" <<endl
    << endl
    << "C1" <<endl;

  for ( vector<Atom>::iterator it = atoms.begin();
        it != atoms.end(); it++) {
//    cout <<"atom="<< vbas->id << " " << vbas->range  << " "
//         <<     vbas->bas_L.size() <<endl;
    int basisid=it->basis; 
    int atomnum =atomname2number(it->name); 
    // pos[0:2] is in AU, AU-> ANG 
    os <<  it->name << " " << atomnum << " "
         <<  it->pos[0]*ANG << " " << it->pos[1]*ANG << " "
         <<  it->pos[2]*ANG << endl;
    for ( int ibasis=0; ibasis< basisset[basisid].types.size(); ibasis++ ){
         string orbtype=  basisset[basisid].types[ibasis];  
//      cout << vbas_L->orbtype << " " << vbas_L->ncoeff << " "
//           << vbas_L->weight << " " << vbas_L->unknown << endl;
#if 0
       if ( orbtype=="SP" ) {
               os << "  L " ;
       } else {
#endif
               os << "  " << orbtype<< " " ;
#if 0
       }
#endif
       int nexp=basisset[basisid].exponents[ibasis].size();
       os << nexp << " "<< 1 <<endl;
#if 0
      if ( vbas_L->orbtype=="SP" ) {
        for (int i=0;i<nexp;i++) {
                sprintf(buf1,form, basisset[id].exponents[ibasis][i]);
                sprintf(buf2,form, basisset[id].coefficients[ibasis][i]);
                sprintf(buf3,form, vbas_L->coeff2[i]);

          os << "    " << i+1 << "   " <<buf1 << " " << buf2 << " "
               << buf3 <<endl;
        }
      } else {
#endif
        for (int i=0;i<basisset[basisid].exponents[ibasis].size();i++) {
                sprintf(buf1,form, basisset[basisid].exponents[ibasis][i]);
                sprintf(buf2,form, basisset[basisid].coefficients[ibasis][i]);
          os << "    " << i+1 << "   " <<  buf1 << " " << buf2 << endl;
        }
#if 0
      }
#endif
    }
    os <<endl;
  }

  os<< " $END"<<endl;
  return 0;
}

/* write ECP section of gamess */
static int g03gamess_ecpsection( /* input */
                vector<Atom> & atoms,
                vector<Gaussian_pseudo_writer> & pseudos,
		vector<PseudoValence> &pseudovalence, 
                ofstream & os
                /* no output */
                )
{
   const char *t="g03gamess_ecpsection";
  os << " $ECP" <<endl;

  for ( vector< Atom  > ::iterator itatom=atoms.begin();
		  itatom!=atoms.end(); itatom++ ) {
	  // find ecp
	  string name = itatom->name; 
	  int found=0;
	  vector<Gaussian_pseudo_writer>::iterator itps;
	  for (itps=pseudos.begin(); itps !=pseudos.end(); itps++) {
               if ( strcasecmp(name.c_str(), itps->label.c_str())==0 ) { found=1; break;}
	  }
	  int found2=0;
	  vector<PseudoValence>::iterator itp;
	  for (itp=pseudovalence.begin(); itp!=pseudovalence.end(); itp++) {
               if ( strcasecmp(name.c_str(), itp->label.c_str())==0 ) { found2=1; break;}
	  }
	  if (( found && !found2)  ||  (!found && found2) ) {
               cout << t << " error "<<endl;
	       exit(ERR_CODE);
	  }
	  if (found && found2) {
		  int maxl = itps->nvalue.size()-1;
		  int corenum = (int) (itp->atomnum - itp->valencenum ); 
            os << name <<"-ECP GEN "  << corenum << " " << maxl <<endl;
	    for (int ibas=0;ibas<itps->nvalue.size(); ibas++) {
		    string Lstr;
		    if (itp->Lnum[ibas]== itp->refL[ibas] ) {
			    Lstr = L2Lstr(itp->Lnum[ibas]);
		    }
		    else {
                           Lstr = string(L2Lstr(itp->Lnum[ibas]))+string( "-") +string( L2Lstr(itp->refL[ibas]));
		    }
		    os << " " <<itps->nvalue[ibas].size() << " " << Lstr <<endl;
		    for (int in=0;in<itps->nvalue[ibas].size(); in++) {
                          os << "   "<<itps->coefficients[ibas][in] << " " 
				  << itps->nvalue[ibas][in] +2 << " " 
				  << itps->exponents[ibas][in] << endl;
		    }
	    }
          }
	  else {
             os << name << "-ECP" <<endl;
	  }
#if 0 
     int id=vvpsp->id-1;
     vector<g03PSP_L> psL = vvpsp->psL;
     string generate;
     if (psL.size()==0) { generate =""; }
     else { generate="GEN "; }
     if ( alreadythesameps(pseudos,id) ) {
        os << atoms[id].name << "-ECP "<< generate << endl;
        cout << "id="<< id+1 << " have the same ps" <<endl;
        continue;
     }
     os << atoms[id].name << "-ECP " <<generate << vvpsp->atomnum-vvpsp->valence  << " "
             << findlmax( vvpsp->psL ) <<endl;

  //    cout << "num of PS " << psL.size()<<endl;;
    for ( vector<g03PSP_L>::iterator vpsp=psL.begin();
          vpsp!=psL.end(); vpsp++ ) {
 //     cout << "L=" << vpsp->L << " refL=" << vpsp->refL <<endl;
      if (vpsp->L == vpsp->refL) {
        os << vpsp->power.size() << "  " <<  L2Lstr(vpsp->L) <<endl;
      }
      else {
        os << vpsp->power.size() << "  " <<  L2Lstr(vpsp->L) << "-" << L2Lstr(vpsp->refL) << endl;
      }
      for (int i=0;i< vpsp->power.size(); i++) {
        os <<"   " << vpsp->coeff[i] << " " << vpsp->power[i] <<  " " << vpsp->expn[i] <<endl;
      }
    }
#endif
  }

  os << " $END" <<endl;

  return 0;
}


/* make gamess input */
int g03gamess_make_inp(/* input */
		int multiplicity,
		string &scftype,
                vector<Atom> & atoms,
                vector<Gaussian_basis_set>  &basisset,
		vector<Gaussian_pseudo_writer> &ecp, 
		vector<PseudoValence> &pseudovalence, 
		vector< vector<double> > & alphamo,
		vector< vector<double> > & betamo,
		string & filename
		/* no output */)
{
   const char *t="g03gamess_convert";
   ofstream os(filename.c_str());

  if(!os) {
    cout <<t<< "Couldn't open " << filename << endl;
    exit(ERR_CODE);
  }

  os.setf(ios_base::scientific | ios_base::uppercase);
  //os.precision(8); 
  os.precision(9);
   
  int nbasis = alphamo[0].size();
  g03gamess_headersection(multiplicity, scftype, nbasis, ecp, os);
  g03gamess_datasection( atoms, basisset, os );
  g03gamess_ecpsection( atoms, ecp, pseudovalence, os);
  g03gamess_vectorsection( nbasis, alphamo, betamo , os);

  cout << "gamess inputfile to \'" << filename <<"\'"<< endl;
  return 0;
}

#if 0
int g03gamess_make_out(
                /* input */
               string & outfilename,
               vector<g03Atom> & atoms,
               string & scftype,
               double eref,
               vector <g03basis> & basis,
               vector <g03PSP>  & psp,
               int nup, int ndown
               /* output */
               )
{
   const char *t = "g03gamess_make_out; ";
   ofstream os(outfilename.c_str());
   if(!os) {
      cout <<t<< "Couldn't open " << outfilename << endl;
      exit(ERR_CODE);
   }

   os.setf(ios_base::scientific | ios_base::uppercase );
   os.precision(9);


         // atoms = Ang unit 
	 
   os << " ATOM      ATOMIC                      COORDINATES (BOHR)" <<endl;
   os<< "           CHARGE         X                   Y                   Z" <<endl;
   for (vector<g03Atom>::iterator itatoms = atoms.begin();
                   itatoms!= atoms.end(); itatoms++) {
       os <<" " << itatoms->name    << " " << itatoms->atomnum << " "
          << itatoms->pos[0]/ANG << " " << itatoms->pos[1]/ANG << " " <<  itatoms->pos[2]/ANG <<endl;
             // ANG -> AU 
   }
   os << endl;

   int zcore;
   if ( psp.empty() ){   /* no PS */
       zcore=0;
   }
   else {
   /* zcore */
   if ( atoms.size() != basis.size() || atoms.size() != psp.size() ) {
           cout << t << "sizes inconsistent " << atoms.size() << " "
                   << basis.size() << " " << psp.size() <<endl;
           exit(ERR_CODE);
   }
   zcore=0;
   for ( vector<g03PSP>::iterator itpsp=psp.begin();
                   itpsp!= psp.end(); itpsp++ ) {
      zcore += itpsp->atomnum-itpsp->valence;
   }
   if (zcore%2==1 ) {
           cout << t << "zcore is odd "<<endl;
           exit(ERR_CODE);
   }
   }

   /* this is all electron */
   os << " NUMBER OF OCCUPIED ORBITALS (ALPHA)          =  " << nup+ zcore/2 <<endl;
   os << " NUMBER OF OCCUPIED ORBITALS (BETA )          =  " << ndown+zcore/2 <<endl;

   os << " FINAL UHF ENERGY IS " << eref <<endl;
   os << " SCFTYP="<< scftype <<endl;

   os << "            --------------" << endl
      << "            ECP POTENTIALS" << endl
      << "            --------------" << endl 
      <<endl;

   for (vector<g03PSP>::iterator ipsp = psp.begin(); ipsp!=psp.end(); ipsp++) {
       int id = ipsp->id-1;
       int zcore = ipsp->atomnum - ipsp->valence; 
       int lmax = findlmax(ipsp->psL);
       if (ipsp->refid==0 ) continue;
       if (ipsp->refid==ipsp->id )  {
       os << " PARAMETERS FOR \""<< atoms[id].name << "-ECP \" ON ATOM  " << id+1 
	       << "  WITH ZCORE  " << zcore << " AND LMAX " << lmax << " ARE "<<endl;
       for ( vector<g03PSP_L>::iterator ipspL = ipsp->psL.begin(); 
		       ipspL != ipsp->psL.end(); ipspL++ ) {
            os << " FOR L= " << ipspL->L << "   COEFF    N           ZETA" <<endl;
            for ( int i=0;i< ipspL->expn.size(); i++) {
		 os << "   " << i+1 << " " << ipspL->coeff[i]<< " " << ipspL->power[i]<< " " << ipspL->expn[i] <<endl;
	    }
       }
       }
       else {
          os << " PARAMETERS FOR \""<< atoms[id].name << "-ECP \" ON ATOM  " << id+1
		  << " ARE THE SAME AS ATOM  " << ipsp->refid <<endl;
       }
       os << endl;

   }

   os <<endl;
   cout << "gamess output to: " << outfilename 
      << ". This is a file for gamess2qmc" << endl;
   return 0;
}

int g03gamess_make_pun(
	/* input */
		string & punfilename,
		vector <g03Atom> & atoms, 
		vector <g03basis> & basis,
		vector< vector<double> > & alphamo,
		vector< vector<double> > & betamo 
		)
{
	const char *t = "g03gamess_make_pun: ";
    ofstream os (punfilename.c_str());
    if (!os) {
         cout <<t<< "Couldn't open " << punfilename << endl;
         exit(ERR_CODE);
    }
    if ( atoms.size()!= basis.size() ) {
	    cout << t <<  " error size inconsistent "<< atoms.size() << " " << basis.size() <<endl;
	    exit(ERR_CODE);
    }

    os.setf(ios_base::scientific  | ios_base::uppercase);
    os.precision(8);

    os << " $DATA"<<endl
     << endl
     << "C1  0 " <<endl;

     vector<g03basis>::iterator itbasis=basis.begin();
     for ( vector<g03Atom>::iterator itatoms=atoms.begin(); 
          itatoms != atoms.end();  
	  itatoms++, itbasis++ ) {
		  os << itatoms->name << " " << itatoms->atomnum << " " << 
			  itatoms->pos[0] << " " << itatoms->pos[1] << " " << itatoms->pos[2] <<endl;
		  for ( vector<g03basis_L>::iterator itbasL = itbasis->bas_L.begin();
				  itbasL != itbasis->bas_L.end(); itbasL++ ) {
                        if ( itbasL->orbtype=="SP" ||  itbasL->orbtype=="L" ) {
				os << "   L  " << itbasL->ncoeff <<endl;
				for (int i=0;i<itbasL->ncoeff ;i++) {
                                    os <<"     " << i+1 << " " << itbasL->expn[i] <<  " "
					  <<  itbasL->coeff1[i] << " " << itbasL->coeff2[i] <<endl;
				}
			}
			else {
                                os << "   "<< itbasL->orbtype << " "  << itbasL->ncoeff <<endl;
                                for (int i=0;i<itbasL->ncoeff ;i++) {
                                    os << "     " << i+1 << " " << itbasL->expn[i] <<  " "
                                        <<    itbasL->coeff1[i] <<endl;
                                }
                        }
		  }
		  os <<endl;
     }
     os << " $END" << endl<<endl;

     char buf1[10];
     char buf2[5][20];
     int nbasis = alphamo.size();
     const char *form="%15.8E";
     int i1;

     os << " $VEC"<<endl;
     i1=0;
     for (vector< vector<double> >::iterator itvec1 = alphamo.begin();
		     itvec1 != alphamo.end(); itvec1++, i1++ ) {
         vector<double>::iterator itvec2 =  itvec1->begin() ; 
         int i2=0; 
	 while ( itvec2 != itvec1->end() ) {
		 int itot=0;
		 sprintf(buf1,"%2d %2d", (i1+1)%100, (i2+1)%100 );
		 for (int i=0; i<5;i++ ) buf2[i][0]='\0';
		 for (int i=0; i<5;i++ ) {
                     sprintf(buf2[i],form, *itvec2);
   	             itvec2++;
                     if ( itvec2 == itvec1->end() ) break;
		 }
                 os << buf1 ;
		for (int i=0;i<5;i++) {
		       os << buf2[i] ; 	
		}
		os <<endl; i2++;
         }
     }

     i1=0;
     for (vector< vector<double> >::iterator itvec1 = betamo.begin();
                     itvec1 != betamo.end(); itvec1++, i1++ ) {
         vector<double>::iterator itvec2 =  itvec1->begin() ; 
         int i2=0;
         while ( itvec2 != itvec1->end() ) {
                 int itot=0;
                 sprintf(buf1,"%2d %2d", (i1+1)%100, (i2+1)%100 );
                 for (int i=0; i<5;i++ ) buf2[i][0]='\0';
                 for (int i=0; i<5;i++ ) {
                     sprintf(buf2[i],form, *itvec2);
                     itvec2++;
                     if ( itvec2 == itvec1->end() ) break;
                 }
                 os << buf1 ;
                for (int i=0;i<5;i++) {
                       os << buf2[i] ;
                }
                os <<endl; i2++;
         }
     }
    os << " $END" <<endl;
    cout << "gamess punch to: " <<punfilename
     << ". This is a file for gamess2qmc." << endl;
    return 0;
}
#endif



// delete dupilicated data 
int g03_ps_delete_dup( /*input and output */
              vector<Gaussian_pseudo_writer> &pseudolog )

{
	if (pseudolog.empty() ) { return 0; }

  int ndata1=pseudolog.size();
  vector<Gaussian_pseudo_writer> pseudolognew;
  vector<string> namesdone;

  for (vector<Gaussian_pseudo_writer>::iterator it=pseudolog.begin(); it!=pseudolog.end(); it++) {
      string name = it->label;
      if (it->nvalue.empty()) { continue; } // dont add it
      int found=0;
      for ( vector<string>::iterator itname=namesdone.begin(); itname!=namesdone.end(); itname++) {
          if (name == *itname) { found=1; break; }
      }
      if (found==0) {
              pseudolognew.push_back(*it);
              namesdone.push_back(name);
      }
  }

  pseudolog.clear();
  pseudolog = pseudolognew;

  int ndata2=pseudolog.size();

  return ndata2-ndata1; 
}

// valid data: return 1
// invalid data, something is broken : return 0
//
int g03_ps_check_validity( /*input*/
              vector<Gaussian_pseudo_writer> &pseudolog )
{
  int flag=0;
  for ( vector<Gaussian_pseudo_writer>::iterator it = pseudolog.begin(); it!=pseudolog.end(); it++) {
          for (int icset=0;icset<it->coefficients.size(); icset++) {
                  for (int icset2 =0; icset2 < it->coefficients[icset].size(); icset2++) {
                          if ( it->coefficients[icset][icset2] ==0.0 ||
                                          it->exponents[icset][icset2] ==0 ) {
                                  flag++;
                                  goto foundzero;
                          }
                  }
          }
  }

foundzero:
   if (flag==0) {
       // do nothing
       return 1;
   }
   return 0;
}


