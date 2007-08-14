/*
 
 Copyright (C) 2007 Hiori Kino (kino dot hiori atmark nims dot go dot jp )
 
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

#include <string.h>
#include <stdlib.h>

#include "converter.h"
#include "g03tools.h"

#include <cmath>

using namespace std;


void g03Atom::print_atom(ostream & os) 
{
       os << setw(5) << name << " "
          << setw(5) <<  id << " "
          << setw(5) << atomnum << " " ;
       for (int i=0;i<pos.size(); i++)  os << setw(15) << pos[i] <<" " ;
       os <<std::endl;
};


int g03basis::print(ostream &os)
{
    g03basis *vbas;
    vbas=this;
    os <<"atom="<< vbas->id << " " << vbas->range  << " "
         <<     vbas->bas_L.size() <<endl;
    for ( vector<g03basis_L>::iterator vbas_L = vbas->bas_L.begin();
          vbas_L != vbas->bas_L.end(); vbas_L++) {
      os << vbas_L->orbtype << " " << vbas_L->ncoeff << " "
           << vbas_L->weight << " " << vbas_L->unknown << endl;
      if ( vbas_L->orbtype=="SP" ) {
        for (int i=0;i<vbas_L->expn.size();i++) {
          os << "   " <<vbas_L->expn[i] << " " << vbas_L->coeff1[i] << " "
               << vbas_L->coeff2[i] <<endl;
        }
      } else {
        for (int i=0;i<vbas_L->expn.size();i++) {
          os << "   " <<  setprecision(10) << vbas_L->expn[i] << " " 
		  << setprecision(10) << vbas_L->coeff1[i] << endl;
        }
      }
    }
    return 0;
}

void g03_print_basisset(vector<g03basis> &bas)
{
    for (vector<g03basis>::iterator it = bas.begin(); it!= bas.end(); it++) {
       it->print(cout);
    }
}



int g03basiskind::print(ostream &os)
{
    g03basiskind *vbas;
    vbas=this;

    os << "atom=" ;
    for (vector<int>::iterator ids = vbas->ids.begin();  ids!=  vbas->ids.end(); ids++) {
       os << " " << *ids ;
    }
    os <<endl;
#if 0
    os <<"atom="<< vbas->id << " " << vbas->range  << " "
         <<     vbas->bas_L.size() <<endl;
#endif
    for ( vector<g03basis_L>::iterator vbas_L = vbas->bas_L.begin();
          vbas_L != vbas->bas_L.end(); vbas_L++) {
      os << vbas_L->orbtype << " " << vbas_L->ncoeff << " "
           << vbas_L->weight << " " << vbas_L->unknown << endl;
      if ( vbas_L->orbtype=="SP" ) {
        for (int i=0;i<vbas_L->expn.size();i++) {
          os << "   " <<vbas_L->expn[i] << " " << vbas_L->coeff1[i] << " "
               << vbas_L->coeff2[i] <<endl;
        }
      } else {
        for (int i=0;i<vbas_L->expn.size();i++) {
          os << "   " << setprecision(10) <<  vbas_L->expn[i] << " " 
		  << setprecision(10) << vbas_L->coeff1[i] << endl;
        }
      }
    }
    return 0;
}

void g03_print_basisetkind(vector<g03basiskind> &baskind) 
{
    for (vector<g03basiskind>::iterator it = baskind.begin(); it!= baskind.end(); it++) {
       it->print(cout);
    }
}

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




/* find max L of ps */
int findlmax( vector<g03PSP_L> & psL )
{
   int lmax=0;
   for ( vector<g03PSP_L>::iterator it = psL.begin(); it!= psL.end();   it++) {
      if (lmax < it->L ) { lmax = it->L ; }
   }
   return lmax;
}



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



#if 0
/* print basis set */
static int g03_print_basiset( /* input */ vector<g03basis> & basisset )
{
  cout << "basisset" <<endl;
  for ( vector<g03basis>::iterator vbas = basisset.begin();
        vbas != basisset.end(); vbas++) {
    cout <<"atom="<< vbas->id << " " << vbas->range  << " "
         <<     vbas->bas_L.size() <<endl;
    for ( vector<g03basis_L>::iterator vbas_L = vbas->bas_L.begin();
          vbas_L != vbas->bas_L.end(); vbas_L++) {
      cout << vbas_L->orbtype << " " << vbas_L->ncoeff << " "
           << vbas_L->weight << " " << vbas_L->unknown << endl;
      if ( vbas_L->orbtype=="SP" ) {
        for (int i=0;i<vbas_L->expn.size();i++) {
          cout << "   " <<vbas_L->expn[i] << " " << vbas_L->coeff1[i] << " "
               << vbas_L->coeff2[i] <<endl;
        }
      } else {
        for (int i=0;i<vbas_L->expn.size();i++) {
          cout << "   " <<  vbas_L->expn[i] << " " << vbas_L->coeff1[i] << endl;
        }
      }
    }
  }
  return 0;
}
#endif

#if 0
/* print basissetkind */
static int g03_print_basisetkind( /* input */ vector<g03basiskind> & basisset )
{
  cout << "basisset" <<endl;
  for ( vector<g03basis>::iterator vbas = basisset.begin();
        vbas != basisset.end(); vbas++) {
        cout << "atom=" ;
        for (vector<int>::iterator ids=vbas->ids; ids!=vbas->ids->end(); ids++) {
           cout << " " << *ids  << endl;
        }
#if 0
    cout <<"atom="<< vbas->id << " " << vbas->range  << " "
         <<     vbas->bas_L.size() <<endl;
#endif
    for ( vector<g03basis_L>::iterator vbas_L = vbas->bas_L.begin();
          vbas_L != vbas->bas_L.end(); vbas_L++) {
      cout << vbas_L->orbtype << " " << vbas_L->ncoeff << " "
           << vbas_L->weight << " " << vbas_L->unknown << endl;
      if ( vbas_L->orbtype=="SP" ) {
        for (int i=0;i<vbas_L->expn.size();i++) {
          cout << "   " <<vbas_L->expn[i] << " " << vbas_L->coeff1[i] << " "
               << vbas_L->coeff2[i] <<endl;
        }
      } else {
        for (int i=0;i<vbas_L->expn.size();i++) {
          cout << "   " <<  vbas_L->expn[i] << " " << vbas_L->coeff1[i] << endl;
        }
      }
    }
  }
  return 0;
}
#endif


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

/* print vector< g03PSP > */
 int g03_print_ps(/* input */
          vector< g03PSP  > pseudos)
{
  cout <<"pseudopotential" <<endl;
  for ( vector< g03PSP  > ::iterator vvpsp= pseudos.begin();
        vvpsp != pseudos.end(); vvpsp++ ) {
    cout <<"atom "<< vvpsp->id << " " << vvpsp->atomnum << " " << vvpsp->valence <<endl;
    vector<g03PSP_L> psL = vvpsp->psL;
    cout << "num of PS " << psL.size()<<endl;;
    for ( vector<g03PSP_L>::iterator vpsp=psL.begin();
          vpsp!=psL.end(); vpsp++ ) {
      cout << "L=" << vpsp->L << " refL=" << vpsp->refL <<endl;
      for (int i=0;i< vpsp->power.size(); i++) {
        cout <<"   " << vpsp->power[i] << " " 
             << setw(15) << setprecision(10) << vpsp->expn[i] << " " 
             << setw(15) << setprecision(10) << vpsp->coeff[i] <<endl;
      }
    }
  }
  return 0;
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
          cout << "unknown type " << basis[bas].types[i] << endl;
          exit(1);
        }
      }
    }
    //cout <<" func " << func << endl;
  }

  return 0;
}

int g03_change_form(
		/* input */
		int nup, int ndown, 
		string & scftype, 
                vector <g03Atom> & g03atoms,
                vector<g03basis> & g03basisset,
		vector<g03PSP> & g03psp, 
		vector< vector< double> > & alphamo,
		vector< vector< double> > & betamo,
		/* output */
                Slat_wf_writer & slwriter,
                vector <Atom>  & atoms,
                vector <Gaussian_basis_set> & basis,
		vector <Gaussian_pseudo_writer> & pseudo,
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

  atoms.clear();
  if ( g03psp.empty() ) {
  for (vector<g03Atom>::iterator it = g03atoms.begin();
          it != g03atoms.end(); it++ ) {
      Atom a_atom;
      a_atom.name=it->name;
      a_atom.charge = it->atomnum;
      a_atom.pos=it->pos; 
#if 1
      { for (int i=0; i<a_atom.pos.size();i++) {
           a_atom.pos[i] /= ANG;   // ANG -> AU 
					       }
      }
#endif
      a_atom.basis=it->id-1;
      atoms.push_back(a_atom);
  }
  }
  else {
  vector<g03PSP>::iterator itps = g03psp.begin();
  for (vector<g03Atom>::iterator it = g03atoms.begin();
          it != g03atoms.end(); it++ ) {
      Atom a_atom;
      a_atom.name=it->name;
      a_atom.charge = itps->valence; 
      a_atom.pos=it->pos; 
#if 1
      { for (int i=0; i<a_atom.pos.size();i++) {
           a_atom.pos[i] /= ANG;   // ANG -> AU
                                               }
      }
#endif
      a_atom.basis=it->id-1;
      atoms.push_back(a_atom);
      itps++;
  }
  }
  for (vector<Atom> :: iterator it= atoms.begin(); it!= atoms.end(); it++ ) {
      it->print_atom(cout);
  }

  // basis 
  for ( vector<g03basis>:: iterator g03atom = g03basisset.begin(); 
		  g03atom != g03basisset.end(); g03atom++) { 
	Gaussian_basis_set tmpbasis;
	int id = g03atom->id-1;
        tmpbasis.label= atoms[id].name; 
	for ( vector < g03basis_L >::iterator basL = g03atom->bas_L.begin(); 
			basL != g03atom->bas_L.end();  basL ++ ) {
		string orbtype = basL->orbtype;
		string orb=orbtype;
		if ( orbtype=="SP" || orbtype=="L" ) {orb="S";} // later set "P"; 
		else if ( orbtype=="D" ) {orb = "6D";}
		else if ( orbtype=="F" ) {orb = "10F";}
		else if ( orbtype=="G" ) {orb = "15G";}
		tmpbasis.types.push_back(orb);
		tmpbasis.exponents.push_back(emptyvector);
		tmpbasis.coefficients.push_back(emptyvector);
		int place = tmpbasis.types.size()-1;
		for (int iex=0;iex< basL->ncoeff; iex++ ) {
			tmpbasis.exponents[place].push_back(basL->expn[iex]);
			tmpbasis.coefficients[place].push_back(basL->coeff1[iex]);
		}
		if ( orbtype=="SP") { 
			orb="P";
                	tmpbasis.types.push_back(orb);
			tmpbasis.exponents.push_back(emptyvector);
			tmpbasis.coefficients.push_back(emptyvector);
                	int place = tmpbasis.types.size()-1;
                	for (int iex=0;iex< basL->ncoeff; iex++ ) {
                        	tmpbasis.exponents[place].push_back(basL->expn[iex]);
                        	tmpbasis.coefficients[place].push_back(basL->coeff2[iex]);
                	}
                }
        }
        basis.push_back(tmpbasis);
  }

#if 0
  cout << "-------basis set--------" <<endl;
  for ( vector<Gaussian_basis_set> ::iterator bas = basis.begin();
		   bas != basis.end(); bas++ ) {
          bas->print_basis(cout);
  }
#endif


  vector<string> addedname;  addedname.clear();

  for ( vector<g03PSP>::iterator atomps =  g03psp.begin() ; 
		  atomps != g03psp.end(); atomps++ ) {
	  cout << "g03psp" <<endl;
        int id=atomps->id-1;
	Gaussian_pseudo_writer pseudotmp;
	pseudotmp.label = atoms[id].name;

        int found=0;
        for (vector<string>::iterator itname=addedname.begin();
              itname!=addedname.end(); itname++) {
            if (pseudotmp.label == *itname )  {
               found=1; break;
            }
        }
        if (found) continue;

        pseudotmp.atomnum = atomps->atomnum;
        if (atomps->psL.empty()) {
           vector<double> emptyr; emptyr.push_back(0.0);
           vector<int> emptyi; emptyi.push_back(-2);
           pseudotmp.exponents.push_back(emptyr);
           pseudotmp.coefficients.push_back(emptyr);
           pseudotmp.nvalue.push_back(emptyi);
        }
        else {
	for (vector<g03PSP_L>::iterator psL = atomps->psL.begin();
			psL != atomps->psL.end(); psL++ ) {
	   pseudotmp.exponents.push_back(psL->expn);
	   pseudotmp.coefficients.push_back(psL-> coeff);
	   vector<int> pw = psL->power; 
	   for ( int i=0; i<pw.size(); i++ ) {
                pw[i] = pw[i]-2;   //  in the print_pseudo, this value +2 again
	   }
	   pseudotmp.nvalue.push_back(pw);

	}
        }
	pseudo.push_back(pseudotmp);
        addedname.push_back(pseudotmp.label);
  }
  cout << "# of pseudo="<< pseudo.size()<<endl;

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

int find_atomname(vector<g03Atom> &atomlist, string& key)
{
   for (int i=0;i<atomlist.size(); i++) {
        if (atomlist[i].name==key) { return i;}
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


int g03_baskindset2basset(/* input */
		vector<g03basiskind> &baskindset,
		/* output */
		vector<g03basis> &basisset)
{
   const char *t = "g03_baskindset2basset: ";
   // find max_of_id in baskindset
   int maxid=0;
   for (vector<g03basiskind>::iterator baskind = baskindset.begin();  baskind!=baskindset.end(); baskind++) {
       for (vector<int>::iterator id = baskind->ids.begin(); id!=baskind->ids.end(); id++) {
           if  (maxid < *id) maxid = *id; 
       }
   }

   //cout << t << "max id ="<<maxid<<endl;

   g03basis bas;
   // make null arrays
   for (int i=0;i<maxid;i++) {
	   basisset.push_back(bas); 
   }

   // copy the contents
   for (vector<g03basiskind>::iterator baskind = baskindset.begin();  baskind!=baskindset.end(); baskind++) {
       for (vector<int>::iterator id = baskind->ids.begin(); id!=baskind->ids.end(); id++) {
	    //   cout << t << "copy id="<< *id <<endl;
            basisset[*id-1].id = *id; basisset[*id-1].range=0;
            basisset[*id-1].bas_L = baskind->bas_L;
       }
   }    

   // check all the basisset is here
   int err_flag=0;
   for (vector<g03basis>::iterator bas= basisset.begin(); bas!=basisset.end(); bas++) {
       if (bas->id==0) {
	       err_flag=1;
	       break;
       }
       if (bas->bas_L.size()==0) {
	       err_flag=1;
	       break;
       }
   }
   if (err_flag) {
	   cout << t << "some basis is missing"<<endl;
	   exit(ERR_CODE);
   }

   // print
   //g03_print_basisset(basisset);

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
#if 0
           if ( strstr(line_upper.c_str()," 5D ")!=NULL ) {
               basis_6d=1;
           }
           if ( strstr(line_upper.c_str()," 6D ")!=NULL ) {
               basis_6d=2;
           }
           if ( strstr(line_upper.c_str()," 7F ")!=NULL ) {
               basis_10f=1;
           }
           if ( strstr(line_upper.c_str()," 10F ")!=NULL ) {
               basis_10f=2;
           }
#endif
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

#if 0
    basis_dftype="";
    if (basis_6d==1) {
	    basis_dftype= " 5D";
    }
    else if (basis_6d==2) {
            basis_dftype= " 6D";
    }

    if (basis_10f==1) {
            basis_dftype+= " 7F ";
    }
    else if (basis_10f==2) {
            basis_dftype+= " 10F ";
    }
    cout << "matches '" << comment << basis_dftype <<"'" << endl;
#endif

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
		vector<g03basis> & basisset,
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

#if 0
   // now accept only 6D and  10F
   if ( dtype==0 || ftype==0 ) {
      cout << "Error: now accepts only 6D and 10F"<<endl;
      exit(ERR_CODE);
   }
#endif

#if 0
   // if one accept 5D, 7F..., one must remake this subroutine 
   //
   const int N_num_Lorb=4;
   const int num_Lorb[N_num_Lorb]={1,3,6,10}; // number of S P 6D and 10F 
#endif

   // find max of L
   // if L<="D" , do nothing
   int Lmax=0;
   for ( vector<g03basis>::iterator itbasis=basisset.begin(); itbasis!= basisset.end(); itbasis++) {
      for ( vector<g03basis_L>::iterator itbasisL=itbasis->bas_L.begin(); itbasisL!=itbasis->bas_L.end(); itbasisL++) {
          string Lstr = itbasisL->orbtype;
          int L= g03_maxangularMomentum(Lstr);
          if (Lmax<L) Lmax=L;
          if ( L > Lmaxorb ) {
             cout << t << "L= " <<  L << " not supported"<<endl;
             exit(ERR_CODE);
          }
      }
   }

   // error check, 
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
   for ( vector<g03basis>::iterator itbasis=basisset.begin(); itbasis!= basisset.end(); itbasis++) {
      for ( vector<g03basis_L>::iterator itbasisL=itbasis->bas_L.begin(); itbasisL!=itbasis->bas_L.end(); itbasisL++) {
          string Lstr = itbasisL->orbtype;
          vecsize += g03_basissizeL(Lstr);
      }
   }

#if 0
   // count size of eigenvector
   int vecsize=0;
   for ( vector<g03basis>::iterator itbasis=basisset.begin(); itbasis!= basisset.end(); itbasis++) {
      for ( vector<g03basis_L>::iterator itbasisL=itbasis->bas_L.begin(); itbasisL!=itbasis->bas_L.end(); itbasisL++) {
          string Lstr = itbasisL->orbtype;
          int L= g03_angularMomentum(Lstr);
	  //cout << Lstr << " " << L <<" " << num_Lorb[L] << endl ;
          if ( L > Lmaxorb ) {
             cout << t << "L= " <<  L << " not supported"<<endl;
             exit(ERR_CODE);
          }
          vecsize+= num_Lorb[L]  ;
      }
   }
#endif

   // vecsize == alphamo[0].size()  ? 
   if ( vecsize != alphamo[0].size() ) {
      cout << t << "size of eigenvector != size of basisset" <<endl;
      cout << vecsize << " " << alphamo[0].size()<<endl;
      exit(ERR_CODE);
   }

#if 0
   // make order of L
   vector< vector<int> >  Lorder;
   for (int iL=0;iL<N_num_Lorb; iL++) {
      vector<int> iorder;
      switch(iL) {
      case 0:
      case 1:
      case 2:
      for (int io=0;io<num_Lorb[iL]; io++) {
          iorder.push_back(io);
      }
      break;
      case 3:
      for (int io=0;io<num_Lorb[iL]; io++) {
          iorder.push_back(forder[io]-1);
      }
      break;
      default:
        // this can not happen, but 
        cout << t << "not supported"<<endl;
        exit(ERR_CODE);
      break;
      }
      Lorder.push_back(iorder);
   }
#endif

   cout << t << "There exist F orbials"<<endl;

   // make order 
   vector<int> order;
   int itot=0;
   for ( vector<g03basis>::iterator itbasis=basisset.begin(); itbasis!= basisset.end(); itbasis++) {
      for ( vector<g03basis_L>::iterator itbasisL=itbasis->bas_L.begin(); itbasisL!=itbasis->bas_L.end(); itbasisL++) {
          string Lstr = itbasisL->orbtype;
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
#if 0
          int L= g03_angularMomentum(Lstr);
	  for (int i=0; i < Lorder[L].size() ; i++) {
              order.push_back(itot+Lorder[L][i]);
	  }
	  itot += Lorder[L].size() ; 
#endif
          
      }
   }

   // check size again
   if ( order.size() != vecsize) {
	   cout << t << "internal error, size differs"<<endl;
	   cout << order.size() << " " << vecsize <<endl;
	   exit(ERR_CODE);
   }

#if 1
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


int ps_uniq_id(/* input */
		int ips ,  vector<g03PSP> pseudo)
{
  const static char *t="ps_uniq_id: ";
    for (int id=0; id<ips; id++) {
       if ( pseudo[id].atomnum == pseudo[ips].atomnum ) { return id; }
    }
    return ips; /* uniq */
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
		vector<g03Atom> &g03atoms)
{
   const static char *t="g03_com_processatomgeometry:";
   string spc=" ";
   vector<string> words;
   g03atoms.clear();
   int id=0;

   for ( vector<string>::iterator it=lines.begin(); it!=lines.end(); it++) {
      if ((*it)[0]==' ') break;
      words.clear(); split(*it,spc,words);
      g03Atom atom;
      if (words[0]=="X") continue;
      atom.name = words[0];
      atom.id = ++id;
      atom.atomnum = atomname2number(atom.name);
      g03atoms.push_back(atom);
   } 
   //cout << t << "number of atoms = "<< g03atoms.size() << endl;

#if 0
   { // count Xatomnum;
   int Xatomnum = 0;
   for (vector<g03Atom>::iterator it= g03atoms.begin();  it!= g03atoms.end() ; it++) {
      if ( it->name=="X" || it->name == "x" ) Xatomnum++;
   }
   if (Xatomnum>0) {
      cout << t << " X or x atom(s) found, not supported" <<endl;
      exit(ERR_CODE);
   }
   }
#endif
   return 0;
}

int g03_com_processPS( /* input */
          vector<string> lines,
	  vector<g03Atom> &g03atoms,
          /* output */
          vector<g03PSP> &g03ps)
{
   const static char *t="g03_com_processPS:";

   vector<string> words;
   string spc=" ";
   vector<string>::iterator it=lines.begin(); 
   while (1) {
    g03PSP atom_pseudo;

    // next must be atom and 0
    words.clear(); split(*it,spc,words);
    string atomname=words[0];
    atom_pseudo.atomnum = atomname2number(atomname);

    int id = find_atomname(g03atoms,atomname)+1;
    if (id>0) {
      atom_pseudo.id=id; // the first id found
    }
    else {
      cout << t << "atomname '" << atomname << "' not found" <<endl;
      exit(ERR_CODE);
    }
    //cout << t << " ps, read " << atomname <<" " << id <<endl;

    
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
      g03PSP_L pseudoL;
      it++; // "D", "S - D",..., angular momentum
      int L, refL;
      //cout << *it <<endl;
      L = label_to_L( *it, refL);
      //cout << "L,refL= "<< L << " " << refL <<endl;
      pseudoL.L=L;
      pseudoL.refL=refL;
      pseudoL.power.clear(); pseudoL.expn.clear(); pseudoL.coeff.clear();

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
	pseudoL.power.push_back(pow);
	pseudoL.expn.push_back(expn);
	pseudoL.coeff.push_back(coeff);
      }
      atom_pseudo.psL.push_back(pseudoL);
    } // for

    g03ps.push_back(atom_pseudo);

    *it++;
    //cout << " next="<<line <<endl;
    if ( it==lines.end() ) {
      break;
    }

  } // while

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

/* read fchkfileame,  set basisset */
int g03_fchk_readbasisset(/* input */ string & fchkfilename,
		                /* output */
		        vector<g03basis> & basisset )
{
   const char *t="g03_fchk_readps";
   int idummy=0;
   int nout;
   string key;
  
   basisset.clear();

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
  g03basis_L basisL;
  g03basis basis;
  int iprimitive=0;
  for ( int ishell=0; ishell< n_shellmap; ishell++ ) {
     if ( basis.id==0 ) {
          basis.id=shellmap[ishell];
          basis.range=0;
     }
     else if (basis.id != shellmap[ishell] ) {
          basisset.push_back(basis);
          basis.clear();
          basis.id=shellmap[ishell];
          basis.range=0;
     }

     basisL.clear();
     basisL.orbtype = g03fchk_shelltype2orbtype( shelltypes[ishell],dtype,ftype );
     basisL.ncoeff = noprimitives[ishell];
     basisL.weight=1.0;
     basisL.unknown=0.0;
     for ( int ip_shell =0 ; ip_shell < noprimitives[ishell]; ip_shell++ ) {
          basisL.expn.push_back(primitive_exp[iprimitive]);
          basisL.coeff1.push_back( contract_coeff[iprimitive]);
          if ( shelltypes[ishell] == -1 && n_p_contract_coeff>0  ) {
//		 cout << "coeff2 " << p_contract_coeff[iprimitive] << endl;
             basisL.coeff2.push_back( p_contract_coeff[iprimitive]);
          }
          iprimitive++;
     }
     basis.bas_L.push_back(basisL);

  }
  if ( basis.id!=0 ) {
        basisset.push_back(basis);
  }

//  g03_print_basiset(basisset);

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

#if 0
  ifstream is(logfilename.c_str());

  if(!is) {
    cout <<t<< "Couldn't open " << logfilename << endl;
    exit(ERR_CODE);
  }

  string line;
  string space=" ";
  vector<string> words;
  string name;
  string sub;

  orbname.clear();

  while (getline(is,line) ) {
     if ( strncmp(line.c_str(),"     EIGENVALUES -- ",20)==0 ) {
	//     cout << line <<endl;

        while ( getline(is,line) ) {
//   1 1   O  1S       
          sub = line.substr(0,20);
	//  cout << "sub="<<sub <<endl;
          if ( sub.find_first_not_of(" ") == string::npos ) { goto last; }
          name = line.substr(11,6);
	//  cout << name <<endl;
          orbname.push_back(name); 
	}
     }
  }
#endif

last:

#if 0
  if ( nbasis != orbname.size() ) {
	  cout << t << " error in size of orbital names " << endl;
	  cout << nbasis << " " << orbname.size() <<endl;
	  exit(ERR_CODE);
  }
#endif

  int nbasis = orbname.size();
#if 0
  cout << "orbital names" <<endl;
  for (int i=0;i<nbasis;i++) {
	 cout << orbname[i] ; 
  }
  cout <<endl;
#endif

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





/* read basissetstring of gaussian format,
 * set basisset 
 */
int g03_log_analyze_basissetfromcards(/* input */
		vector<string> & basissetstring, 
		/* output */
		vector<g03basis> & basisset )
{
  const char *t = "g03_log_analyze_basissetfromcards: ";
  string g03seperator="****";
  vector<string> words;
  string sep=" ";


  // first set basis set to baskindset.
  // then copy baskindset to basisset
 
  g03basiskind baskind;
  vector<g03basiskind> baskindset;

  g03basis_L bas_L;
//  g03basis bas;

  vector<string>::iterator thisline=basissetstring.begin();

  while ( thisline!= basissetstring.end() ) {

    //cout << "first=" << *thisline <<endl;
    //*thisline =  "Centers:  1 3 ..."
    words.clear(); split(*thisline,sep,words);
    if (words[0]!="Centers:") {
        cout << t << "error, unknown format"<<endl;
	cout << *thisline <<endl;
	exit(ERR_CODE);
    }
    // read it
    vector<int> ids;
    for (int i=1;i<words.size();i++) {
	int id = atoi(words[i].c_str());
	ids.push_back(id);
    }

    // read each bas
    baskind.clear();
    baskind.ids=ids;
     
    thisline++; 
    while ( thisline!= basissetstring.end() ) {


      if (thisline->find(g03seperator) != string::npos ) {
//	cout << "found ****"<<endl;	      
	baskindset.push_back(baskind);
	baskind.clear();
	break;
      }

      //  SP   1 1.00        .000000000000
      words.clear(); split(*thisline,sep,words);
      if (words.size()!=4) {
        cout <<t<< "error in basis set, basis set title size!=4 "
	     <<words.size()<<endl;
	cout << *thisline <<endl;
        exit(ERR_CODE);
      }
      string orbtype=words[0];
      int ncoeff = atoi(words[1].c_str());
      double weight=atod(words[2].c_str());
      double unknown=atod(words[3].c_str());
      if (unknown!=0.0) {
        cout <<t<< "error in basis set, unknown field!=0.0 "
	     <<unknown<<endl;
	cout << *thisline<<endl;
        exit(ERR_CODE);
      }

      bas_L.clear();
      bas_L.orbtype=orbtype; bas_L.ncoeff=ncoeff; bas_L.weight=weight; bas_L.unknown=unknown;

      vector<double> expn; expn.clear();
      vector<double> coeff1,coeff2;  coeff1.clear();  coeff2.clear();

      thisline++;
      //       Exponent=  3.7267155000D+01 Coefficients=  3.9037801000D-02    for non-SP
      //      or
      //      Exponent=  3.7267155000D+01 Coefficients=  3.9037801000D-02   3.9037801000D-02   for SP
      for (int i=0;i<ncoeff;i++) {
	words.clear();
	split(*thisline,sep,words);
	if (orbtype=="SP") {
	  if (words.size()!=5) {
	    cout <<t<< "error in basis set, unknown L definition"
		 <<unknown<<endl;
	    cout << *thisline <<endl;
	    exit(ERR_CODE);
	  }
	  double val1,val2,val3;
	  val1 = atod(words[1].c_str());
	  val2 = atod(words[3].c_str());
	  val3 = atod(words[4].c_str());
	  expn.push_back(val1);
	  coeff1.push_back(val2);
	  coeff2.push_back(val3);
	}
	else {
	  if (words.size()!=4) {
	    cout <<t<< "error in basis set, unknown L definition"
		 <<unknown<<endl;
	    cout << *thisline <<endl;
	    exit(ERR_CODE);
	  }
	  double val1,val2;
	  val1 = atod(words[1].c_str());
	  val2 = atod(words[3].c_str());
	  expn.push_back(val1);
	  coeff1.push_back(val2);
	}
	thisline++;
      } // for i=

//      cout << id <<" " << orbtype << endl;
//      for (int i=0;i<expn.size();i++) {
//	cout << expn[i] <<" " << coeff1[i] <<endl;
//      }

      bas_L.expn= expn;
      bas_L.coeff1=coeff1;
      bas_L.coeff2=coeff2;
      baskind.bas_L.push_back(bas_L);
//      cout << "push bas_L " << bas_L.orbtype << " " << bas_L.ncoeff <<endl;


    } //while 

    thisline++;

  } //while  

#if 0
  g03_print_basisetkind(baskindset);
#endif

  // baskindset -> basisset 
  g03_baskindset2basset(baskindset,basisset);  

  return 0;
}


/* read basissetstring of gaussian format,
 * set basisset 
 */
int g03_log_analyze_basisset(/* input */
		vector<string> & basissetstring, 
		/* output */
		vector<g03basis> & basisset )
{
  const char *t = "g03_log_analyze_basisset: ";
  string g03seperator="****";
  vector<string> words;
  string sep=" ";

  int iatom=0;

  g03basis_L bas_L;
  g03basis bas;

  vector<string>::iterator thisline=basissetstring.begin();
  while ( thisline!= basissetstring.end() ) {

    iatom++;

    //   1 0
    words.clear();
    split(*thisline,sep,words);
    if (words.size()!=2) {
      cout <<t<< "error in basis set, size()!=2,"<< words.size() << endl;
      cout << *thisline <<endl;
      exit(ERR_CODE);
    }
    int id=atoi(words[0].c_str());
    int range=atoi(words[1].c_str());  
    if ( id!=iatom ) {
      cout <<t<< "error in basis set, id!=iatom "<< id <<endl;
      cout << *thisline <<endl;
      exit(ERR_CODE);
    }
    if (range!=0) {
      cout <<t<< "error in basis set, range !=0,"<<range <<endl;
      exit(ERR_CODE);
    }

    bas.clear();
    bas.id=id;
    bas.range=range;
	
     
    thisline++; 
    while ( thisline!= basissetstring.end() ) {

      if (thisline->find(g03seperator) != string::npos ) {
//	cout << "found ****"<<endl;	      
	basisset.push_back(bas);
	bas.clear();
	break;
      }

      //  SP   1 1.00        .000000000000
      words.clear();
      split(*thisline,sep,words);
      if (words.size()!=4) {
        cout <<t<< "error in basis set, basis set title size!=4 "
	     <<words.size()<<endl;
	cout << *thisline <<endl;
        exit(ERR_CODE);
      }
      string orbtype=words[0];
      int ncoeff = atoi(words[1].c_str());
      double weight=atod(words[2].c_str());
      double unknown=atod(words[3].c_str());
      if (unknown!=0.0) {
        cout <<t<< "error in basis set, unknown field!=0.0"
	     <<unknown<<endl;
        exit(ERR_CODE);
      }

      bas_L.clear();
      bas_L.orbtype=orbtype;
      bas_L.ncoeff=ncoeff;
      bas_L.weight=weight;
      bas_L.unknown=unknown;

      vector<double> expn; expn.clear();
      vector<double> coeff1,coeff2;  coeff1.clear();  coeff2.clear();

      thisline++;
      //        .2073000000D+01   .1022977051D+00   .4596910139D+00
      //        .6471000000D+00   .9176052296D+00   .6322148170D+00
      for (int i=0;i<ncoeff;i++) {
	words.clear();
	split(*thisline,sep,words);
	if (orbtype=="SP") {
	  if (words.size()!=3) {
	    cout <<t<< "error in basis set, unknown field!=0.0"
		 <<unknown<<endl;
	    exit(ERR_CODE);
	  }
	  double val1,val2,val3;
	  val1 = atod(words[0].c_str());
	  val2 = atod(words[1].c_str());
	  val3 = atod(words[2].c_str());
	  expn.push_back(val1);
	  coeff1.push_back(val2);
	  coeff2.push_back(val3);
	}
	else {
	  if (words.size()!=2) {
	    cout <<t<< "error in basis set, unknown field!=0.0"
		 <<unknown<<endl;
	    exit(ERR_CODE);
	  }
	  double val1,val2;
	  val1 = atod(words[0].c_str());
	  val2 = atod(words[1].c_str());
	  expn.push_back(val1);
	  coeff1.push_back(val2);
	}
	thisline++;
      } // for i=

//      cout << id <<" " << orbtype << endl;
//      for (int i=0;i<expn.size();i++) {
//	cout << expn[i] <<" " << coeff1[i] <<endl;
//      }

      bas_L.expn= expn;
      bas_L.coeff1=coeff1;
      bas_L.coeff2=coeff2;
      bas.bas_L.push_back(bas_L);
//      cout << "push bas_L " << bas_L.orbtype << " " << bas_L.ncoeff <<endl;


    } //while 

    thisline++;

  } //while  

#if 0
  g03_print_basiset(basisset);
#endif

  return 0;
}



/* read psstring, set pseudos */
int g03_log_analyze_ps(/* input */ 
		vector<string> & psstring,
		/* output */
	vector<g03PSP> & pseudos	)
{
  const char *t = "g03_log_analyze_ps: ";
  string sep=" ";
  vector<string> words;
  int iatom;
  vector<int> nval;
  vector<double> expn;
  vector<double> coeff;


  g03PSP_L  pseudoL;
  g03PSP atom_pseudo;

  int id,atomid,valencenum;

  pseudoL.L=-1; // initialize 
  iatom=0;

  for ( vector<string>::iterator thisline=psstring.begin();
	thisline != psstring.end(); thisline++) {

    words.clear();
    split(*thisline,sep,words);
    string centernumberstr= thisline->substr(0,10);
    int centernumber = atoi(centernumberstr.c_str());

    if ( thisline->find("No pseudopotential") !=  string::npos ) {
 //         cout << "found No Pseudopotential "<<endl;
           atom_pseudo.refid=0; 
    }
    else if  ( thisline->find("Pseudopotential same as on center") != string::npos ) {
	    atom_pseudo.refid=atoi(words[5].c_str());
    }
    else if ( centernumber >0 ) {
      //   1        8            6    .000000   .000000   .243053
      //  ---> new ps
  //    cout << "found new ps " << endl;

      if ( pseudoL.L>=0 ) {
	atom_pseudo.psL.push_back(pseudoL);
      }
      if ( atom_pseudo.id>0 ) {
	pseudos.push_back(atom_pseudo);
      }

      iatom++;

      id = atoi(words[0].c_str());
      atomid=atoi(words[1].c_str());
      valencenum=atomid;
      if ( words.size()==6) {
	valencenum=atoi(words[2].c_str());
      }
      if (id!=iatom) {
	cout <<t<< "error in ps, id!=iatom " <<  id <<" " << iatom << endl;
	exit(ERR_CODE);
      }

      // clear pseudo
      atom_pseudo.clear();
      atom_pseudo.id=id;
      atom_pseudo.atomnum=atomid;
      atom_pseudo.valence=valencenum;
      atom_pseudo.refid=id;

      pseudoL.clear(); 
    }
    else if ( isalpha(*(words[0].c_str())) ) {
//	    cout << "found new L " <<endl;
      if ( pseudoL.L >=0 ) {
	// add it 
	atom_pseudo.psL.push_back(pseudoL);
      }
      // new angular momentum
      int refL;
      int L = label_to_L(*thisline,refL);
      pseudoL.L = L;
      pseudoL.refL = refL;
      pseudoL.power.clear();
      pseudoL.expn.clear();
      pseudoL.coeff.clear();
    }
    else if ( words.size()==3 ) {
      // ps
   //   cout << "found ps" << endl;
      int pow = atoi(words[0].c_str());
      double expn = atod(words[1].c_str());
      double coeff = atod(words[2].c_str());
      pseudoL.power.push_back(pow);
      pseudoL.expn.push_back(expn);
      pseudoL.coeff.push_back(coeff);
    }

  }

  if ( pseudoL.L>=0 ) {
    atom_pseudo.psL.push_back(pseudoL);
  }
  if ( atom_pseudo.id>0 ) {
    // add
    pseudos.push_back(atom_pseudo);
  }

  // fix ps
  for ( int id = 0; id < pseudos.size(); id++ ) {
     int uniqid = ps_uniq_id(id, pseudos); 
     cout << " ps id=" << id << ", uniqid= " << uniqid <<endl;
     if ( id != uniqid )  pseudos[id].refid= uniqid+1; 
  }

#if 0
  g03_print_ps(pseudos);
#endif
#if 0
  cout <<"pseudopotential" <<endl;
  for ( vector< g03PSP  > ::iterator vvpsp= pseudos.begin();
	vvpsp != pseudos.end(); vvpsp++ ) {
    cout <<"atom "<< vvpsp->id << " " << vvpsp->atomnum << " " << vvpsp->valence <<endl;
    vector<g03PSP_L> psL = vvpsp->psL;
    cout << "num of PS " << psL.size()<<endl;;
    for ( vector<g03PSP_L>::iterator vpsp=psL.begin();
	  vpsp!=psL.end(); vpsp++ ) {
      cout << "L=" << vpsp->L << " refL=" << vpsp->refL <<endl;
      for (int i=0;i< vpsp->power.size(); i++) {
	cout <<"   " << vpsp->power[i] << " " << vpsp->expn[i] << " " << vpsp->coeff[i] <<endl;
      }
    }
  }
#endif
  return 0;
}


int g03_log_readatoms( /*input*/ vector<string> & strv,
                /*output*/ vector<g03Atom> & atoms )
{
        const char *t="main_log2atoms: ";
        atoms.clear();
  string space=" ";
  vector<string> words;
  int itot=0;
  for (vector<string>::iterator it = strv.begin(); it!=strv.end(); it++ ) {
          itot++;
       words.clear();
       split(*it,space,words);
       g03Atom a_atom;
       /* 0 id,  1 atomicnumber, 2 atomictype=0, 3-5 positon */
       if ( words.size()!= 6) { cout << t << " error, size != 6 "<<endl; exit(ERR_CODE); }
       int id=atoi(words[0].c_str());
       a_atom.id=id;
       if (itot!=id) {  cout << t << " error, wrong id?  "<<endl; exit(ERR_CODE); }
       int atomnum = atoi(words[1].c_str());
       if (atomnum==-2) { continue; } // -2 is translational vector 
       a_atom.name= element_lookup_caps_null[atomnum];
       a_atom.atomnum=atomnum;
       a_atom.pos.clear();
       for (int i=3;i<=5;i++) {
          double r = atod(words[i].c_str()); // unit is  Ang 
          a_atom.pos.push_back(r);
       }
       atoms.push_back(a_atom);
   }
  return 0;
}



/* write data section of gamess */
static int g03gamess_datasection(   /* input */
		vector<g03Atom> & atoms,
		vector<g03basis>  & basisset,
                ofstream & os
		)
{
	const char *t="g03gamess_datasection";
	char buf1[30], buf2[30], buf3[30];
	const char *form="%16.9E"; 

 os << " $DATA" <<endl
    << endl
    << "C1" <<endl;

  for ( vector<g03basis>::iterator vbas = basisset.begin();
        vbas != basisset.end(); vbas++) {
//    cout <<"atom="<< vbas->id << " " << vbas->range  << " "
//         <<     vbas->bas_L.size() <<endl;
    int id= vbas->id-1;
    os <<  atoms[id].name << " " <<atoms[id].atomnum << " " 
	 <<  atoms[id].pos[0] << " " << atoms[id].pos[1] << " "
         <<  atoms[id].pos[2] << endl;  
    for ( vector<g03basis_L>::iterator vbas_L = vbas->bas_L.begin();
          vbas_L != vbas->bas_L.end(); vbas_L++) {
//      cout << vbas_L->orbtype << " " << vbas_L->ncoeff << " "
//           << vbas_L->weight << " " << vbas_L->unknown << endl;
       if ( vbas_L->orbtype=="SP" ) {
	       os << "  L " ;
       } else {
               os << "  " << vbas_L->orbtype<< " " ;
       }
       os << vbas_L->ncoeff << " "<< vbas_L->weight <<endl;
      if ( vbas_L->orbtype=="SP" ) {
        for (int i=0;i<vbas_L->expn.size();i++) {
		sprintf(buf1,form, vbas_L->expn[i]);
		sprintf(buf2,form, vbas_L->coeff1[i]);
		sprintf(buf3,form, vbas_L->coeff2[i]);

          os << "    " << i+1 << "   " <<buf1 << " " << buf2 << " "
               << buf3 <<endl;
        }
      } else {
        for (int i=0;i<vbas_L->expn.size();i++) {
		sprintf(buf1,form, vbas_L->expn[i]);
                sprintf(buf2,form, vbas_L->coeff1[i]);
          os << "    " << i+1 << "   " <<  buf1 << " " << buf2 << endl;
        }
      }
    }
    os <<endl;
  }

  os<< " $END"<<endl;

  return 0;
}






/* check whether psudos[psid] already appeared or not  */

static int alreadythesameps(/* input */
                vector<g03PSP> & pseudos,
                int psid)
{
    const char *t="alreadythesameps: ";
    double psatomnum = pseudos[psid].atomnum; 
    for ( int i=0; i< psid; i++ ) {
#if 0
       cout << "psatomnum=" << psatomnum << " pseudos[i].atomnum=" << pseudos[i].atomnum <<endl;
#endif
       if ( psatomnum == pseudos[i].atomnum ) { return 1; } 
    }
    return 0;
}

/* write ECP section of gamess */
static int g03gamess_ecpsection( /* input */
		vector<g03Atom> & atoms,
		vector<g03PSP> & pseudos,
                ofstream & os
		/* no output */
		)
{
   const char *t="g03gamess_ecpsection";
  os << " $ECP" <<endl;
  
  for ( vector< g03PSP  > ::iterator vvpsp= pseudos.begin();
        vvpsp != pseudos.end(); vvpsp++ ) {
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

  }

  os << " $END" <<endl;

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
                vector<g03PSP> & ecp, 
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

/* make gamess input */
int g03gamess_make_inp(/* input */
		int multiplicity,
		string &scftype,
                vector<g03Atom> & atoms,
                vector<g03basis>  & basisset,
		vector<g03PSP> &ecp, 
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
  g03gamess_ecpsection( atoms, ecp, os);
  g03gamess_vectorsection( nbasis, alphamo, betamo , os);

  cout << "gamess inputfile to \'" << filename <<"\'"<< endl;
  return 0;
}

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

