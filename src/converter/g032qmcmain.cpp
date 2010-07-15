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


//
// bug report to kino dot hiori atmark nims dot go dot jp
// atmark -> @
// dot -> . and delete space
//
// This code is awkward. I know it.
//
//  first version: Nov 6, 2006.
//
//  last version:  Aug. 22, 2007.
//
//
//
// How to run
//
// g03 part
//
// iop(3/24=10)
// Printing of gaussian function table.
// 10      Print as GenBas input.
//
// iop(3/18=1)
// PRINTING OF PSEUDOPOTENTIALS
// 1        PRINT
//
// Pseudo=read
// Read pseudo-potential data from the input stream. Input is described in the next subsection below. Cards is a synonym for Read.
//
// Gen
// The Gen keyword allows a user-specified basis set to be used in a Gaussian calculation. 
//
// better using 
// SCF(tight) 
// and 
// Integral(grid=ultrafine) or Integral(grid=finegrid) 
// if you use diffuse basis functions 
//
//
// Prepare g03 input
// Add   "6D 10F IOP(3/24=10) IOP(3/18=1)" options in the # file
// Make chk file by adding '%chk=foo.chk' line
//
// E.g.,  com file (foo.com) looks like
//    %chk=foo.chk
//    #P B3LYP/gen 6d 10f scf opt test iop(3/24=10) iop(3/18=1) pseudo=read
//    ...
//
// Run g03
// > g03 foo.com
// You will get foo.log. Of course, normally exitted files are necessary.
//
// Then make fchk file after running g03
// >  formchk foo.chk foo.fchk 
//
// log and fchk is necessary. com is optionally necessary in the case when the PS is broken in the log file.
//
// g032qmc
// usage:
// > g032gmc -log foo.log -fchk foo.fchk 
// maybe
// > g032qmc -log foo.log -fchk foo.fchk -com foo.com  -w -gamess 
// 
// You can see options.
// > g032qmc -h 
//
//
//
//  atoms and their geometries 
//  		read fchk, also read log for consistency,  optionally read com for consistency 
//  		comment:	geometries are in the 'standard' orientation (not 'input' orientation)
//  alphamo, betamo 
//  		read fchk
//  basisset 
//  		read fchk
//  		comment: log has more accuracy than fchk 
//  pseudopotential (PS)
//  		read log, also optionally read com if it is broken. 
//  		comment: no PS in fchk.
//  totalenergy 
//  		read fchk
//  nup,ndown,multiplicity
//  		read fchk
//
// The main reason to read log is to check consistency of the system calculated. 
//
//
//
//  PBC extension
//
//  iop(5/33=1) 
//  1        THE EIGENVALUES AND THE M. O. COEFFICIENTS ARE PRINTED AT THE END OF THE SCF.
//  comment: all the eigenvalues, but MO coeffcients are for k=Gamma. 
//
//  iop(5/98=1) 
//  Whether to save eigenvalues and orbitals at all k-points. 
//  1        Yes.
//  comment: formchk writes eigevector of only k=Gamma.
//
//  iop(5/103=N)
//  Number of occupied and virtual orbitals to print for each k-point. 
//  N       N occupieds and N virtuals.
//
//  PBC(NKpoint=N)
//  Do approximately N k-points.
//
//
#if 0

CI calculation
#P hf/6-31G(d,p)  SCF Test cisd  iop(9/28=-1)


 Normalization: A(0)=1
...
 Final wavefunction coefficients:
 Dominant configurations:
 ***********************
 Spin Case        I    J    A    B          Value
    AA            3        10             .141106D-02
    AA            3        13             .180006D-03
...
    AA            8        40            -.840573D-03
   AAAA           3    5   10   13       -.485889D-02
   AAAA           3    5   10   14        .118347D-02
   AAAA           3    5   10   16        .162436D-02
   AAAA           3    5   10   20        .173122D-02
   AAAA           3    5   10   24       -.839045D-03
   AAAA           3    5   10   28       -.216302D-02
...
   AAAA           6    8   34   31       -.724529D-04
   ABAB           3    3   10   10       -.488000D-02
   ABAB           3    3   10   13        .847809D-03
   ABAB           3    3   10   14        .275269D-02
   ABAB           3    3   10   16       -.828540D-03
...
   ABAB           8    6   40   37       -.512640D-03
   ABAB           8    6   40   39        .181167D-02
 Largest amplitude= 1.24D-01

#endif

#include <algorithm>
#include <iostream>
#include <iomanip>
#include <fstream>

#include <string>
#include <vector>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <stdlib.h>

#include "converter.h"
#include "Pseudo_writer.h"
#include "basis_writer.h"
#include "wf_writer.h"

#include "g03tools.h"


using namespace std;





class option_type {
  public:
  char *option_str;
  char *defstr;
  int  (*func)(char*, char**);
  char **ptr;
  int incl;
  char *opton_explanation;
} ;

int main_analyze_option_com(char *opt, char **comfilename)
{
//   cout << "log "<<opt<< endl;
      *comfilename = strdup(opt);
      return 0;
}


int main_analyze_option_log(char *opt, char **logfilename)
{
//   cout << "log "<<opt<< endl;
   *logfilename = strdup(opt);
   return 0;
}

int main_analyze_option_fchk(char *opt,  char **fchkfilename)
{
//  cout << "fchk "<<opt << endl;
  *fchkfilename=strdup(opt);
  return 0;
}

int main_analyze_option_out(char *opt,  char **outfilename)
{
//  cout << "fchk "<<opt << endl;
  if (opt) {
  *outfilename=strdup(opt);
  }
  return 0;
}

int main_print_options(int nopt,option_type *opts)
{
  cout << "options"<<endl;
  for (int i=0;i<nopt;i++) {
    cout << TAB <<setw(5)<<left<< opts[i].option_str
         << " " << opts[i].opton_explanation <<endl;
  }
  return 0;
}

int main_analyze_option_help(char *opt, char **dummy)
{
  return 1;
}

int main_analyze_option_w(char *opt, char **overwrite)
{
  *overwrite="y";
  return 0;
}


int main_analyze_option_gamess(char *opt, char **print_gamess)
{
   *print_gamess = "y";
   return 0;
}


int main_analyze_option(/*input*/
         int argc,char ** argv,
         /*output*/
	 char** comfilename, 
         char** logfilename,char** fchkfilename,
         char** outfilename, char **overwritestr,
	 char ** print_gamess
	 )
{
   static const int nopt=7;
   int flag;
   /* old fashoned ? */
   option_type opts[nopt]=
  {
   {"-log",  NULL,   &main_analyze_option_log,  logfilename, 1,
    "filename \"g03 log filename (necessary)\""},

   {"-fchk", NULL,   &main_analyze_option_fchk, fchkfilename,1,
    "filename \"g03 fchk filename (necessary)\""},

   {"-com",  NULL,   &main_analyze_option_com,  comfilename, 1,
    "filename \"g03 com filename (optional)\""},

   {"-out", "sample",&main_analyze_option_out, outfilename,1,
    "filename \"output filename without extension (default=sample)\""},

   {"-w",    "n",   &main_analyze_option_w, overwritestr, 0,
    "\t\t\"overwrite existing files (default=off)\""},

   {"-gamess","n",   &main_analyze_option_gamess, print_gamess, 0,
    "\t\"print gamess files (default=off)\""},

   {"-h",    NULL,   &main_analyze_option_help, NULL,        0,
    "\t\t\"print help\""}
  };

  for (int io=0;io<nopt;io++) {
    if ( opts[io].defstr ) *(opts[io].ptr) = opts[io].defstr;
  }

  for (int iargc=1;iargc<argc; iargc++) {
    for (int io=0;io<nopt;io++) {
       if ( strcmp(argv[iargc], opts[io].option_str)==0 ) {
           int flag;
           if ( opts[io].incl==0 ) {
            flag =opts[io].func(NULL, opts[io].ptr);
           }
           else if ( iargc+1 < argc ) {
            flag =opts[io].func(argv[iargc+1],opts[io].ptr);
           }
           else {  
            flag=1; /* error */
           }
           if (flag) {
               goto print_help;
           }
           for (int i=0;i<opts[io].incl;i++) { iargc++; }
           goto nextarg;
       }
    }
    cout << "error in arg, unknown option: "<< argv[iargc] <<endl;
    goto print_help;
 nextarg: ;
  }

  if (*logfilename!=NULL && *fchkfilename!=NULL ) { return 0; }
 print_help:
    main_print_options(nopt,opts);
    exit(ERR_CODE);
    return 0;
}


/************************************************************************************/
/************************************************************************************/
/************************************************************************************/




/* 
 * read comfilename 
 *   read   only 
 *   1. atomnames, not geometry 
 *   2. pseudopotential
 */
int main_read_g03comfile( /* input */
		string & comfilename,
		vector<Atom> atoms, 
		/* output */
		vector<Atom> &atomscom, 
		vector<Gaussian_pseudo_writer> &pseudos
		)
{
   const char *t="main_read_g03comfile: ";

   string basis_dftype;

   if (comfilename.size()<=0 ) { return 0; }

   cout << "--- read com file ---" <<endl;

   ifstream is(comfilename.c_str());

   if ( !is) {
      cout << t<< "Couldn't open " << comfilename << endl;
       exit(ERR_CODE);
   }

   string line;
   string spc=" ";
   vector<string> words;
   int found_gen=0;
   int found_pseudo=0;
   // search # line
   words.clear();
   while (getline(is,line)) {
      if ( line[0]=='#' ) {
	   //cout << "# line="<<endl;
	   words.push_back( line );
	   while ( getline(is,line) ) {
               if (  line.find_first_not_of(" ")==string::npos )  goto exit_searchsharp;
	       words.push_back( line );
	   }
      }   
   }
   cout << t << "# line, not found." <<endl;
   exit(ERR_CODE);

exit_searchsharp: 
   for (vector<string>::iterator it=words.begin(); it!=words.end(); it++) {
      cout << *it <<endl;
   }

   g03_process_sharpheader(words, found_gen, found_pseudo);

   getline(is,line); // title
   getline(is,line); // blank
   getline(is,line); // charge multiplicity 

   //
   //  atom and geometry
   //
   words.clear();
   while (getline(is,line)) { 
      if ( line.find_first_not_of(" ")==string::npos ) {
          break;
      }
      words.push_back(line);
   }
   g03_com_processatomgeometry( words, atomscom);
   // may have variables section
   while ( line.find_first_not_of(" ")!=string::npos ) {
        getline(is,line);
   } 

   // next is basis set, if found_gen!=0, 
   // but skip them
   //cout << line <<endl;
   if (found_gen) {
     while (getline(is,line)) { // atom and geometry
        if ( line.find_first_not_of(" ")==string::npos ) {
           break;
        }
	//cout <<"basis:"<< line <<endl;
     } 
   }

   // next is ps
   if (found_pseudo) {
//     g03ps.clear();
     words.clear();

     while (getline(is,line)) {
       if ( line[0]=='!') { continue; }
       else { break; }
     }

     while (1) {
        if ( line.find_first_not_of(" ")==string::npos ) {
           break;
        }
	//cout << "ps:" << line<<endl;
        words.push_back(line);
	getline(is,line);
     }
     g03_com_processPS(words,atoms, pseudos);  
     int ndeleted = g03_ps_delete_dup(/*input and output*/
              pseudos);
     if (ndeleted>0) {
              cout<< ndeleted << " PS(com) data delted"<<endl;
     }

   }
   // check it.
   // g03_print_ps(g03ps);

   return 0;
}

/* read logfilename and get parameters */
int main_read_g03logfile( /* input */
	    	string & logfilename, 
		int pbcsystem, 
		/* output */
			      vector<Atom> &atomslog,
			      vector<Gaussian_basis_set> &gbasisset, 
		              vector<Gaussian_pseudo_writer> &pseudos, 
			      vector<PseudoValence> &pseudovalence 
		)

{
	cout << "--- read log file ---" <<endl;

  const char *t="main_read_g03logfile: ";

  int ret;

  string keyin;
  string keyout;
  int nread, nskip;
  vector<string> strv; 
  string str;

  string space=" ";
  vector<string> words;

  keyin="General basis"; 
  ret=filekeyinput_get_string(
        /* input */
        logfilename,
        keyin,
       /* output */
        str); 
  if ( ret==0 ) {
     keyin="Standard basis:";
     ret=filekeyinput_get_string(
        /* input */
        logfilename,
        keyin,
       /* output */
        str);
  }
  //cout << "logfile: basis df type? " << str<<endl;
  //basis_dftype = str;

  strv.clear();

  keyin="#";
  nskip=-1;  // include # line 
  nread=0; keyout="-------------------";
  ret=filekeyinput_get_stringv(
       /* input */
        logfilename,
        keyin,  nskip,
        nread,  keyout, NULL,
       /* output */
        strv );

  { 
   int found_gen, found_pseudo; 
   string basis_dftype; 
   g03_process_sharpheader(strv, found_gen, found_pseudo);
  }
  



  keyin="Pseudopotential Parameters";
  nskip=4; nread=0; keyout="==================";
  ret=filekeyinput_get_stringv(
       /* input */
        logfilename,
        keyin,  nskip,
        nread,  keyout, NULL,
       /* output */
        strv );

  g03log_analyze_ps(/*input*/ strv,
		  /*output*/ pseudos, pseudovalence); 
  int ndeleted =  g03_ps_delete_dup(/*input and output */
                    pseudos);
  if (ndeleted>0) {
      cout<< ndeleted << " PS(log) data delted"<<endl;
  }



  if (pbcsystem) {
  keyin="Input orientation:";
  nskip=4; nread=0; keyout="-------------------";
  ret=filekeyinput_get_stringv_last(
       /* input */
        logfilename,
        keyin,  nskip,
        nread,  keyout, NULL,
       /* output */
        strv );
  }
  else {
  keyin="Standard orientation:";
  nskip=4; nread=0; keyout="-------------------";
  ret=filekeyinput_get_stringv_last(
       /* input */
        logfilename,
        keyin,  nskip,
        nread,  keyout, NULL,
       /* output */
        strv );
   }

   g03_log_readatoms(strv, // atoms,
		   atomslog);



   if (pbcsystem) {
  keyin="Threshold for overlap eigenvalues";
  nskip=-1;  // include # line
  nread=1; keyout="";
  ret=filekeyinput_get_stringv(
		       /* input */
		        logfilename,
		              keyin,  nskip,
		              nread,  keyout, NULL,
		       /* output */
		        strv );
    cout<< "PBC: "<< strv[0] << ", but discard this information" <<endl;

   }


  // basis set 
  strv.clear();
    keyin="AO basis set in the form of general basis input:";
    nskip= 0; nread=0; keyout="";  /* use lastkey_is_spc */
    ret=filekeyinput_get_stringv(
       /* input */
        logfilename,
        keyin,  nskip,
        nread,  keyout, lastkey_is_spc,
       /* output */
        strv );

      if (strv.size()>0) {
              cout << "basis set is from 'AO basis set' section"<<endl;
        /* string aobasisset -> basisset */
       //  g03_log_analyze_basisset(/* input */ strv,
        //        /* output */ basisset);
         g03log_analyze_basisset(/*input*/ strv,
			atomslog, 
                         /*output*/ gbasisset);
	 g03_basisset_delete_dup(atomslog,gbasisset);

      }


  return 0; 
}




/* read fchkfilename and set parameters */
int main_read_g03fchkfile(/* input */
		string & fchkfilename, 
		/* output */
			       vector<Atom> &gatoms, 
			       vector< vector<double> > &pbcvector, 
			       vector<Gaussian_basis_set> &gbasisset, 
			       vector<vector<double> > &alphamo, 
			       vector<vector<double> > &betamo,
			       int & nup, int & ndown, int & multiplicity ,
			       double &totalenergy ,
			       int &dtype,int &ftype
			       )
{	
	cout << "--- read fchk file ---"<<endl;

  const char *t = "main_read_g03fchkfile: ";
  string key;
  int n,nout;

  int natom;
  vector<int> iv;
  vector<double> rv;
   
  key="Number of atoms";
  g03_fchk_read_i(fchkfilename,key,natom);
  cout << key << ": " << natom <<endl;

  int nbasis; 
  key="Number of basis functions";
  g03_fchk_read_i(fchkfilename,key,nbasis);
  cout << key << ": " << nbasis <<endl;


  key="Multiplicity";
  g03_fchk_read_i(fchkfilename,key,multiplicity);
  cout << key << ": " << multiplicity <<endl;
  key="Number of alpha electrons";
  g03_fchk_read_i(fchkfilename,key,nup);
  cout << key << ": " <<nup <<endl;
  key="Number of beta electrons";
  g03_fchk_read_i(fchkfilename,key,ndown);
  cout << key << ": " <<ndown <<endl;
  key="Total Energy";
  g03_fchk_read_r(fchkfilename,key,totalenergy);
  cout << key << ": " << setprecision(15) << totalenergy <<endl;

  
  key="Pure/Cartesian d shells";
  g03_fchk_read_i(fchkfilename,key,dtype);
     cout << key << ": "<< dtypestr[dtype]<<endl;
  key="Pure/Cartesian f shells";
  g03_fchk_read_i(fchkfilename,key,ftype);
     cout << key << ": "<< ftypestr[ftype] <<endl;
  
  key="Atomic numbers";
  g03_fchk_read_iv(fchkfilename,key,natom,nout,iv);
  for (int i=0;i<natom;i++) {
    Atom a_atom;
    a_atom.charge = iv[i]; /* atomic number , not valence charge */
    a_atom.name = element_lookup_caps_null[iv[i]];
    string_upper(a_atom.name); 
    a_atom.basis= i;
    gatoms.push_back(a_atom);
  }


  g03_fchk_readbasisset(fchkfilename, //basisset, 
		  gatoms, 
		  gbasisset);
           // read dtype and ftype in the subroutine again

  // fix dup
  g03_basisset_delete_dup(gatoms,gbasisset);



  key="Current cartesian coordinates";
  g03_fchk_read_rv(fchkfilename,key,natom*3,nout,rv);
          // unit is AU
  {
    int j=0;
    for (int i=0;i<natom;i++) {
      for (int k=0;k<3;k++) {
//	atoms[i].pos.push_back(rv[j]*ANG); //  AU -> ANG 
	gatoms[i].pos[k]=rv[j];   //*ANG; 
	j++;
      }
    }
  }

  key="Translation vectors";
  g03_fchk_read_rv(fchkfilename,key,9,nout,rv);
  if ( nout>0 ) {
     if (nout==9 ) {
       pbcvector.resize(3);
       for (int i=0;i<3;i++) pbcvector[i].resize(3);
       int itot=0;
       for (int i=0;i<3;i++) 
          for (int j=0;j<3;j++) {
             pbcvector[i][j] = rv[itot++]; //*ANG;
	  }
       cout << "System is PBC!" << endl;
       cout << "Translation vectors (ANG)=" <<endl;
       for (int i=0;i<3;i++) {
         cout << pbcvector[i][0]*ANG << " " << pbcvector[i][1]*ANG  << " "
		 << pbcvector[i][2]*ANG <<endl;
       }
     }
     else {
       cout << t << "found Translation vectors, but their components are not 9"
     	       << endl; 
       exit(ERR_CODE);
     } 
  }


  key="Alpha MO coefficients";
  g03_fchk_read_rv(fchkfilename,key,nbasis*nbasis,nout,rv);
  if (rv.size()!=nbasis*nbasis) {
     cout << t << " error in alpha MO"<<endl;
     exit(ERR_CODE);
  }
  alphamo.clear();
  {
    int j=0;
    for (int i=0;i<nbasis;i++) {
      vector<double> rv1;
      rv1.clear();
      for (int k=0;k<nbasis;k++) {
	rv1.push_back(  rv[j++] );
      }
      alphamo.push_back(rv1);
    }
  }

  key="Beta MO coefficients";
  g03_fchk_read_rv(fchkfilename,key,nbasis*nbasis,nout,rv);
  if (rv.size()!=nbasis*nbasis && rv.size()!=0 ) {
     cout << t << " error in beta MO"<<endl;
     exit(ERR_CODE);
  }
  
  betamo.clear();
  if ( rv.size() == nbasis*nbasis ) 
  {
    int j=0;
    for (int i=0;i<nbasis;i++) {
      vector<double> rv1;
      rv1.clear();
      for (int k=0;k<nbasis;k++) {
	rv1.push_back(  rv[j++] );
      }
      betamo.push_back(rv1);
    }
  }
  return 0;
}



int main_print_qwalk(
   /* input */
     string & outputname, 
     vector<Center> & centers, 
     vector<int> & basisidx,
     vector< vector<double> > & moCoeff,
     vector< Gaussian_basis_set > & basis, 
     Slat_wf_writer &  slwriter , 
     vector< Atom> & atoms,
     vector< vector<double> > & pbcvector, 
     vector<double> &electric_field,
     vector<Gaussian_pseudo_writer> &pseudo,
     vector<double> &origin,
     vector<double> &kpoint, 
     double eref, 
     int write_boilerplate,
     int opt_overwrite
   /* no output */
   )
{
  cout << "writing qwalk input files ... "<<endl;

  int natoms=atoms.size();

  string orboutname=outputname+".orb";
  if ( !opt_overwrite ) {
     if (exist_file(orboutname) ) { 
       cout << orboutname << " exists. stop"<<endl;
       cout << "use -w to overwrite it/them."<<endl;
       exit(ERR_CODE); 
     }
  }
  slwriter.orbname=orboutname;
  string basisoutname=outputname+".basis";
  if ( !opt_overwrite ) {
     if (exist_file(orboutname) ) { 
        cout << orboutname << " exists. stop"<<endl;
        exit(ERR_CODE); 
     }
  }
  slwriter.basisname=basisoutname;

  ofstream orbout(orboutname.c_str());
  print_orbitals(orbout, centers, basisidx, moCoeff);
  orbout.close();

  cout << "output: " << orboutname <<endl;

  ofstream basisout(basisoutname.c_str());
  int nbas=basis.size();
  for(int bas=0; bas < nbas; bas++) {
    //if ( basis_uniq_id(bas,atoms)==bas ) {
    basisout << "BASIS { \n";
    basis[bas].print_basis(basisout);
    basisout << "}\n\n\n";
    //}
  }
  basisout.close();

  cout << "output: " << basisoutname << endl;


  string slateroutname=outputname+".slater";
  if ( !opt_overwrite ) {
     if (exist_file(slateroutname) ) { 
        cout << slateroutname << " exists. stop" <<endl;
        exit(ERR_CODE); 
     }
  }
  ofstream slaterout(slateroutname.c_str());
  slwriter.print_wavefunction(slaterout);
  slaterout.close();
  cout << "output: " << slateroutname <<endl;

  //---------------------------------Jastrow2 output
  string jast2outname=outputname+".jast2";
  if ( !opt_overwrite ) {
     if (exist_file(jast2outname) ) { 
       cout << jast2outname << " exists. stop" <<endl;
       exit(ERR_CODE); 
     }
  }

  double basis_cutoff;
  if (pbcvector.empty()) {  // molecule
     basis_cutoff=7.5; //arbitrary cutoff
  }
  else {   // crystal 
    basis_cutoff=find_basis_cutoff(pbcvector);
  }
  Jastrow2_wf_writer jast2writer;
  jast2writer.set_atoms(atoms);

  ofstream jast2out(jast2outname.c_str());
  print_std_jastrow2(jast2writer, jast2out, basis_cutoff);
  jast2out.close();
  cout << "output: " << jast2outname <<endl;

  //--------------------------------------System output

  string sysoutname=outputname+".sys";
  if ( !opt_overwrite ) {
     if (exist_file(sysoutname) ) { 
       cout << sysoutname << " exists. stop" <<endl;
       exit(ERR_CODE); 
     }
  }
  ofstream sysout(sysoutname.c_str());
  if ( pbcvector.empty() ) { // molecule 
    sysout << "SYSTEM { MOLECULE "<<endl;
    sysout << "  NSPIN { " << slwriter.nup << "  "
	   << slwriter.ndown << " } \n";
  }
  else { // PERIODIC
    sysout << "SYSTEM { PERIODIC" <<endl;
    sysout << "  NSPIN { " << slwriter.nup << "  "
	   << slwriter.ndown << " }" <<endl;
    sysout << " LATTICEVEC {" <<endl;
    for (int i=0;i<3;i++) {
      sysout << pbcvector[i][0] << " " << pbcvector[i][1] << " "
	     << pbcvector[i][2] << endl;
    }
    sysout <<  " }" <<endl;
    sysout << " origin { " << origin[0] << " " << origin[1] << " "
	   <<  origin[2] << " }" << endl;
    sysout << " cutoff_divider " << cutoff_divider(basis, pbcvector) <<endl;
    sysout << " kpoint { " << kpoint[0] << " " << kpoint[1] << " " 
	   << kpoint[2] << " }" <<endl;

  }
  for(int at=0; at <natoms; at++) {
    atoms[at].print_atom(sysout);
  }

  if(electric_field.size() ==3) {
    sysout << "  electric_field { ";
    for(int d=0;d < 3; d++)
      sysout << electric_field[d] << "   ";
    sysout << " } \n";
  }
  sysout << "}\n\n\n";
  // end of SYSTEM section 

  int npsp=pseudo.size();
  for(int psp=0; psp < npsp; psp++) {
    pseudo[psp].print_pseudo(sysout);
  }
  sysout.close();
  cout << "output: " <<  sysoutname <<endl;

  if(write_boilerplate) {
    string hfoutname=outputname+".hf";
  if ( !opt_overwrite ) {
     if (exist_file(hfoutname) ) { 
        cout << hfoutname << " exists. stop"<<endl;
        exit(ERR_CODE); }
  }
    ofstream hfout(hfoutname.c_str());
    print_vmc_section(hfout, outputname, eref);
    hfout << "\n\n";
    hfout << "INCLUDE " << sysoutname << "  \n";
    hfout << "TRIALFUNC { INCLUDE " << slateroutname << "}\n\n";
    hfout.close();
    cout << "output: " << hfoutname <<endl;

    string optoutname=outputname+".opt";
  if ( !opt_overwrite ) {
     if (exist_file(optoutname) ) { 
        cout << optoutname << " exists. stop"<<endl;
        exit(ERR_CODE); }
  }
    ofstream optout(optoutname.c_str());
    print_opt_section(optout, outputname, eref);
    optout << "\n\n";
    optout << "INCLUDE " << sysoutname << " \n";
    optout << "TRIALFUNC { \n  SLATER-JASTROW \n"
           << "  WF1 { INCLUDE " << slateroutname << " } \n"
           << "  WF2 { INCLUDE " << jast2outname   << " } \n"
           << "}\n\n";
    optout.close();
    cout << "output: " << optoutname << endl;
  }

  //cout << "end " << endl;

  return 0;
}


/***************************************************/
int main_check_consistency_g03atoms(
		/*input*/
		vector<Atom> &gatomsfchk,
		vector<Atom> &gatomslog,
		vector<Atom> &gatomcom
		/* no output */)
{
  const char *t="main_check_consistency_g03atoms: ";
  const double eps=1.0e-5;


  // size
  if ( gatomsfchk.size() != gatomslog.size() ) {
	  cout << t << " size differs fchk,log=" << gatomsfchk.size() <<" " << gatomslog.size() <<endl;
	  exit(ERR_CODE);
  }
  for (int iatom=0;iatom< gatomsfchk.size(); iatom++) {
	  if ( gatomsfchk[iatom].name != gatomslog[iatom].name ) {
               cout << t << "name differs at " << iatom <<endl;
	       cout <<  gatomsfchk[iatom].name <<" " <<  gatomslog[iatom].name  <<endl;
	       exit(ERR_CODE);
	  }
	  int flag=0;
	  for (int k=0;k<3;k++) {
		  if ( fabs(gatomsfchk[iatom].pos[k]- gatomslog[iatom].pos[k]) > eps )  {
			  flag++;
		  }
	  }
	  if (flag>0) {
		  cout << t << "position differs at iatom=" << iatom <<endl;
		  cout<< "fchk=" <<endl;
		  gatomsfchk[iatom].print_atom(cout);
		  cout<<"log="<<endl;
                  gatomslog[iatom].print_atom(cout);
		  exit(0);
	  }
  }

  return 0;
}


/************************************************/

int  main_check_consistency_basisset(
				     /* input */
	      vector<Gaussian_basis_set> &gbasissetlog,
	      vector<Gaussian_basis_set> &gbasissetfchk 
				     /* no output */
	      )
{
  const static char *t="main_check_consistency_basisset: ";
  const static double eps=1.0e-5;
  

  if (gbasissetlog.size() != gbasissetfchk.size() ) {
     cout << t << "size differs, log=" <<gbasissetlog.size() << 
	     " fchk=" << gbasissetfchk.size() <<endl;
     exit(ERR_CODE);
  }
  for (int iatom=0;iatom< gbasissetlog.size() ; iatom++) {
	  if ( gbasissetlog[iatom].types.size() !=
			  gbasissetfchk[iatom].types.size() ) {
	    cout << t << "size differs at iatom=" << iatom <<endl;
	    cout << gbasissetlog[iatom].types.size() << " " 
		    << gbasissetfchk[iatom].types.size() <<endl;
	    exit(ERR_CODE);
	  }
          for (int ibasis=0;ibasis< gbasissetlog[iatom].types.size() ; ibasis++) {
              if ( gbasissetlog[iatom].types[ibasis] != 
			      gbasissetfchk[iatom].types[ibasis] ) {
                   cout << t << "type differs at iatom, ibasis="<< iatom 
			   << " " << ibasis <<endl;
		   exit(ERR_CODE);
	      }
	  }
  }
  return 0;
}

int main(int argc, char ** argv)
{
	const char *t="main: ";

   int opt_overwrite=0;
   int opt_printgamess=0;
   string logfilename;
   string fchkfilename;
   string outfilebasename ;
   string comfilename ;
   {
   char *comfile=NULL, *logfile=NULL, *fchkfile=NULL;
   char *outfile=NULL;
   char *overwritestr=NULL;
   char *print_gamess=NULL;
   main_analyze_option(/* input */argc,argv,
         /* output */
         &comfile, &logfile,&fchkfile, &outfile,&overwritestr,  &print_gamess);

   if ( strcmp(overwritestr,"y")==0 ) opt_overwrite=1;
   //cout << "overwrite="<< opt_overwrite <<endl;
   if ( strcmp(print_gamess,"y")==0 ) opt_printgamess=1;


   logfilename=logfile; 
   fchkfilename=fchkfile;
   outfilebasename = outfile;
   if (comfile) comfilename = comfile;
   }

  //
  //---------- end of argv ------------------


//  vector<string>  orbname;
  string scftype;
  double eref;
  int nup, ndown, multiplicity, dtype, ftype;
  int pbcsystem=0; // periodic boundary condition or not

  vector<vector<double> > alphamo, betamo;
  vector<Atom> atomsfchk; 

  // fchk
  //
  // alphamo, beamo, nup, ndown, multiplicity
  // atom and geometry (use this!)
  // basisset

  vector<Gaussian_basis_set> gbasissetfchk;
  vector< vector<double> > pbcvector;
  vector<double> origin(3);
  vector<double> kpoint(3); 
  for (int i=0;i<3;i++) { origin[i]=0.0; kpoint[i]=0.0; }

  main_read_g03fchkfile( /* input */ 
		  fchkfilename, 
		  /* output */ 
	   atomsfchk, pbcvector,
	   gbasissetfchk, alphamo,betamo,nup,ndown,multiplicity, eref,
	   dtype,ftype);
  //cout << "size of vector ="<< alphamo.size() << " " << betamo.size() <<endl;
	  
  if (! pbcvector.empty()) pbcsystem=1;
  if (betamo.empty()) {
      scftype="RHF";
      if (ndown!=nup) {
              cout << t<<" error , ndown must be equal to nup"<<endl;
              exit(ERR_CODE);
      }
  }
  else {
      scftype="UHF";
      cout << "found beta component, ";
  }
  cout << "scftpe="<<scftype <<endl;


  // log 
  vector<Gaussian_basis_set> gbasissetlog; 
  vector<Gaussian_pseudo_writer> pseudoslog;
  vector<PseudoValence> pseudovalence;   // compensate infomation of Gaussian_pseudo_writer
  vector<Atom> atomslog; 

  //
  // atom and geometry
  // basisset  
  // PS
  main_read_g03logfile(/* input */logfilename,  pbcsystem, 
		  /* output */ 
	       atomslog,
	       gbasissetlog,  
	       pseudoslog, 
	       pseudovalence
	       );


  vector<Atom> atomscom; 
  vector<Gaussian_pseudo_writer> pseudoscom;

  // com
  //
  // atoms (only symbols)
  // PS 

  if (! comfilename.empty()) {
	  // atoms are only to check consistency 
	  // use this ps if the ps in the log is broken
      main_read_g03comfile(/* input */ comfilename,  atomsfchk, 
               /* output */   atomscom, pseudoscom ); 
  }

  cout << "--- have read files ---" <<endl;


  // check consistency
  //
  main_check_consistency_g03atoms(
		  atomsfchk,atomslog,atomscom); 

  main_check_consistency_basisset(
		  gbasissetlog, gbasissetfchk );

  
  // PS  may be broken in the log file
  // fix broken PS using pseudoscom
  //
  if (g03_ps_check_validity(/*input*/ pseudoslog)==0 ) { // NG 
	  if (!pseudoscom.empty()){
	  pseudoslog.clear();
	  pseudoslog = pseudoscom;
	  cout << " PS in the log file is broken, use PS in the com file"<<endl;
	  }
	  else {
          cout<< "error: PS is broken, but no additional data  to fix it"<<endl;
	  cout << "supply -com file or fix PS in the log file "<<endl;
	  exit(ERR_CODE);
	  }
  }

  // reorder f orbital, if there are
  // The order of the f orbitals differs. 
  //
  g03_orbital_reorder(/* input */ // basisset,
		  atomsfchk, gbasissetfchk, dtype, ftype, 
		  /* input and output */ alphamo, betamo);


  // ------------------------------------------------------------------
  //  output
  // -------------------------------------------------------------------

  if (opt_printgamess) {
  //output

  {  
    string filename;

    filename=outfilebasename +".inp";
    if ( !opt_overwrite ) {
      if (exist_file(filename) ) { 
	cout << filename << " exists. stop."<<endl;
	exit(ERR_CODE);
      }
    }
    g03gamess_make_inp(
		       /* input */
		       multiplicity,
		       scftype,
		       atomsfchk,
		       gbasissetfchk,
		       pseudoslog,
		       pseudovalence, 
		       alphamo,
		       betamo,
		       filename
		       /* no output */);
    cout << "output: " << filename <<endl;
  }
#if 0
  {
    string outfilename=outfilebasename+".out";
    string punfilename=outfilebasename+".pun";
    if ( !opt_overwrite ) {
      if (exist_file(outfilename) ) { 
        cout << outfilename << " exists. stop"<<endl;
        exit(ERR_CODE);
      }
    }
    g03gamess_make_out(
               /* input */
		       outfilename,
		       g03atoms,
		       scftype,
		       eref,
		       basisset,
		       g03ps,
		       nup,  ndown
		       /* no output */
		       );
    cout << "output: " << outfilename <<endl;
    if ( !opt_overwrite ) {
      if (exist_file(punfilename) ) { 
        cout << punfilename << " exists. stop"<<endl;
        exit(ERR_CODE);
      }
    }

    g03gamess_make_pun(
        /* input */
		       punfilename,
		       g03atoms,
		       basisset,
		       alphamo,
		       betamo
         /* no output */
		       );
    cout << "output: " << punfilename <<endl;

  }
#endif
  }  // opt_printgamess 


#if 1

   Slat_wf_writer  slwriter;
   if (pbcsystem) {
  slwriter.write_centers=false;
  slwriter.use_global_centers=true;
   }
   else {
  slwriter.write_centers=false;
  slwriter.use_global_centers=false;
   }
  slwriter.mo_matrix_type="CUTOFF_MO";

   vector <Gaussian_basis_set> & gamessbasis=gbasissetfchk;
   vector < vector < double> >  moCoeff;
   vector <Gaussian_pseudo_writer> & gamesspseudo=pseudoslog;
   vector <Center>  centers;
   vector <int>  gamessnbasis;
   vector <Atom> & gamessatoms=atomsfchk; 

   // Some parts aren't suitable to use library subroutines 
   // Some parts are not set. 
   // fix them 
   g03_fix_form(
                /* input */
                 nup,  ndown,
                 scftype,
		 dtype,ftype,
                 alphamo,
                 betamo,
		 pseudovalence, 
		 /*input output*/
		 gamessatoms,
		 gamessbasis,
		 gamesspseudo, 
                /* output */
                 slwriter,
                 moCoeff,
		 centers,
		 gamessnbasis
                );

   string outputname=outfilebasename;
   const int write_boilerplate=1; 
   vector<double> electric_field; electric_field.clear();

   main_print_qwalk(
   /* input */
      outputname,
      centers,
      gamessnbasis,
      moCoeff,
      gamessbasis,
       slwriter ,
      gamessatoms,
      pbcvector, 
     electric_field,
     gamesspseudo,
     origin,
     kpoint,
      eref,
      write_boilerplate,
      opt_overwrite
    /* no output */
    );

#endif

  cout << "Done"<<endl;


  cout << endl<< "Warning: Now the basisset must be the same for the same specie of atoms"<<endl<<endl;

  return 0;
}


