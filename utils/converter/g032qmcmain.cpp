//
// bug report to kino dot hiori atmark nims dot go dot jp
// atmark -> @
// dot -> . and delete space
//
// This code is awkward. I know it.
//
//  first version: Nov 6, 2006.
//
//  last version:  Jan. 16, 2007.
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
// log and fchk is necessary. com is optionally necessary.
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
//  		read log, also read fchk for consistency, 
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


/******************************************************************/


#if 0
/* analyze dftypestr and set dtype and ftype if they are 
 *  5d, 6d = 1 or 2  dtype
 *  7f, 10f = 1 or 2  ftype
 *  Now only accepts 6d and 10f basis set 
 *  Maybe no d or f funtions, so this is not necessarily error.
 */
int main_analyze_dftype(/* input */
          string& sub, 
           /* output */
          int &dtype, int &ftype)
{
  const char *t="main_analyze_dftype: ";
  dtype=0; ftype=0;
  // format is  'General basis read from cards:  (5D, 7F)'
  if (  sub.find("5D") != string::npos  ) dtype=1;
  if (  sub.find("6D") != string::npos  ) dtype+=2;
  if (  sub.find("7F") != string::npos  ) ftype=1;
  if (  sub.find("10F") != string::npos  )ftype+=2;

  cout << "dtype="<<dtype<< " ftype="<<ftype<<endl;
  return 0;
}
#endif

#if 0
//
// now qwalk accepts only 6d and 10f basis.
// This subroutine check it.
//
int main_check_basissettype(/* input */
           vector<g03basis> & basisset,
           int dtype, int ftype, 
           string &logfilename
             /* output none */ )
{
  int maxL=0;
  for ( vector<g03basis>::iterator vbas = basisset.begin();
        vbas != basisset.end(); vbas++) {
    for ( vector<g03basis_L>::iterator vbas_L = vbas->bas_L.begin();
          vbas_L != vbas->bas_L.end(); vbas_L++) {
          int tmpL = g03_maxangularMomentum(vbas_L->orbtype);
          if ( tmpL > maxL ) maxL = tmpL;
    }
  }

  if ( maxL <=1 ) { return 0; } // no problem

  int dtype, ftype;
  int flag=1;
#if 0
  main_analyze_dftype(dftypestr, dtype, ftype);
#endif
  if ( maxL<=2)  {
     if (dtype==2)  flag=0;
  }
  else if ( maxL <= 3 ) {
     if ( dtype==2 && ftype==2 ) flag=0;
  }

  if (flag==0) { return 0;}

  cout << "no 6D or 10F keywords, read orbital_names part of the log file again" << endl;

  vector<string> orbname;
  g03_log_analyze_orbname( logfilename,  orbname);
  // if error occurs, stop in g03_log_analyze_orbname
  
  return 0;
}
#endif


/* 
 * read comfilename 
 *   read   only 
 *   1. atomnames, not geometry 
 *   2. pseudopotential
 */
int main_read_g03comfile( /* input */
		string & comfilename,
		/* output */
                vector<g03Atom> &g03atoms, 
		vector<g03PSP> &g03ps
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
#if 0
   if ( strcmp(accept_df,"y")!=0 ) {
     if (basis_dftype.find_first_not_of("6D")==string::npos ||
       basis_dftype.find_first_not_of("10F")==string::npos ) {
       cout << "cann't find 6D or 10F, specify it or confirm that the basis set uses 6D and 10F and use the option."
                <<endl;
       exit(ERR_CODE);
     }
   }
#endif

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
   g03_com_processatomgeometry( words, g03atoms);
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
     g03ps.clear();
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
     g03_com_processPS(words,g03atoms,g03ps);  
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
			      //string &basis_dftype, 
			      vector<g03Atom> &atoms,
                              vector<g03basis> & basisset,
                              vector<g03PSP> &g03ps)
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

  //
  //   BASISSET
  //

  strv.clear();
#if 0
  keyin="General basis read from cards:";
  nskip= 0; nread=0; keyout="AO basis set in the form of general basis input:";
  ret=filekeyinput_get_stringv(
       /* input */
        logfilename,
        keyin,  nskip,
        nread,  keyout, lastkey_is_spc,
       /* output */
        strv );
  if (strv.size()>0) {
	  cout << "basis set is from 'General basis read from cards:' section"<<endl;
     /* string aobasisset -> basisset */
     g03_log_analyze_basissetfromcards(/* input */ strv,
                  /* output */ basisset);
  }
  else {
#endif
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
         g03_log_analyze_basisset(/* input */ strv, 
		  /* output */ basisset);
      }
#if 0
  }
#endif

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
#if 0
   // request 6d and 10f  in the input   
   if ( strcmp(accept_dftype,"y")!=0 ) {
     if (basis_dftype.find_first_not_of("6D")==string::npos || 
       basis_dftype.find_first_not_of("10F")==string::npos ) {
       cout << "cann't find 6D or 10F, specify it or confirm that the basis set uses 6D and 10F and use the option." 
	        <<endl;
       exit(ERR_CODE);
     }
   }
#endif
  }
  

#if 0
  /* check whether the basis set is acceptable to qwalk or not */
  main_check_basissettype(
		  /* input */
              basisset, basis_dftype , logfilename 
	      /* no output */ );
#endif


  keyin="Pseudopotential Parameters";
  nskip=4; nread=0; keyout="==================";
  ret=filekeyinput_get_stringv(
       /* input */
        logfilename,
        keyin,  nskip,
        nread,  keyout, NULL,
       /* output */
        strv );


  /* string psset -> class g03ps */
  g03_log_analyze_ps(/* input */ strv, 
		 /* output */ g03ps);

#if 0
   keyin="SCF Done:";
   ret=filekeyinput_get_string(
        /* input */
        logfilename,
        keyin,
       /* output */
        str);
 
	   words.clear();
	   split(str,space,words);
	   eref=atod(words[4].c_str());
	  // cout << scftype <<" E="<< eref <<endl;
	   cout << "E= "<< setprecision(15) << eref <<endl;
#endif

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

   g03_log_readatoms(strv, atoms);
 


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

  return 0; 
}




/* read fchkfilename and set parameters */
int main_read_g03fchkfile(/* input */
		string & fchkfilename, 
		/* output */
			       vector<g03Atom> &atoms, 
			       vector< vector<double> > &pbcvector, 
			       vector<g03basis> & basisset, 
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
  

  g03_fchk_readbasisset(fchkfilename, basisset);
           // read dtype and ftype in the subroutine again


  key="Atomic numbers";
  g03_fchk_read_iv(fchkfilename,key,natom,nout,iv);
  for (int i=0;i<natom;i++) {
    g03Atom a_atom; a_atom.clear();
    a_atom.id=i+1;
    a_atom.name = element_lookup_caps_null[iv[i]];  
    a_atom.atomnum=iv[i];  /* all electron */
    atoms.push_back(a_atom); 
  }

  key="Current cartesian coordinates";
  g03_fchk_read_rv(fchkfilename,key,natom*3,nout,rv);
          // unit is AU
  {
    int j=0;
    for (int i=0;i<natom;i++) {
      for (int k=0;k<3;k++) {
	atoms[i].pos.push_back(rv[j++]*ANG); //  AU -> ANG 
      }
    }
    for (int i=0;i<natom;i++) {
      atoms[i].print_atom(cout);
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
             pbcvector[i][j] = rv[itot++]*ANG;
	  }
       cout << "System is PBC!" << endl;
       cout << "Translation vectors (ANG)=" <<endl;
       for (int i=0;i<3;i++) {
         cout << pbcvector[i][0] << " " << pbcvector[i][1]  << " "
		 << pbcvector[i][2] <<endl;
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
       cout << orboutname << " exits. stop"<<endl;
       exit(ERR_CODE); 
     }
  }
  slwriter.orbname=orboutname;
  string basisoutname=outputname+".basis";
  if ( !opt_overwrite ) {
     if (exist_file(orboutname) ) { 
        cout << orboutname << " exits. stop"<<endl;
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
    if ( basis_uniq_id(bas,atoms)==bas ) {
    basisout << "BASIS { \n";
    basis[bas].print_basis(basisout);
    basisout << "}\n\n\n";
    }
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
		vector<g03Atom> & atoms,
		vector<g03Atom> & atomsfromlog,
		vector<g03Atom> & atomsfromcom
		/* no output */)
{
  const char *t="main_check_consistency_g03atoms: ";
  const double eps=1.0e-5;

  int check_com = ! atomsfromcom.empty(); 

  if ( atoms.size() != atomsfromlog.size() ) {
    cout<< t<<"number of atoms differs "  
	<<  atoms.size() <<  " " << atomsfromlog.size()
	<<endl;
    exit(ERR_CODE);
  }

  if ( check_com ) {
    if ( atoms.size() != atomsfromcom.size() ) {
      cout<< t<<"number of atoms differs "
	  <<  atoms.size() <<  " " << atomsfromlog.size()
	  << atomsfromcom.size() <<endl;
      exit(ERR_CODE);
    }
  }

  // atomnumber and geometry
  int natom=atoms.size();
  for (int iatom=0;iatom<natom;iatom++) {
    if ( atoms[iatom].atomnum != atomsfromlog[iatom].atomnum ) {
      cout<<t<< "atomname or charge differs (vs log) "
	  <<  atoms[iatom].atomnum << " " << atomsfromlog[iatom].atomnum
	  <<endl;
      cout <<t<< " at atom "<< iatom+1 <<endl;
      exit(ERR_CODE);
    }
    int flag=0;
    int ip;
    for (int ip=0;ip<3;ip++) {
      if ( fabs(atoms[iatom].pos[ip]-atomsfromlog[iatom].pos[ip]) > eps ) {
	flag =1;
	cout <<  fabs(atoms[iatom].pos[ip]-atomsfromlog[iatom].pos[ip]) <<endl;
      }
    }
    if (flag) {
      cout << t<< " posision differs at " << iatom+1 << " "
	   << endl;
      
      atoms[iatom].print_atom(cout);
      atomsfromlog[iatom].print_atom(cout);
      exit(ERR_CODE);
    }
  }

  // data of the com file
  // check only atomnumber
  //
  if ( check_com ) {
    for (int iatom=0;iatom<natom;iatom++) {
      if ( atoms[iatom].atomnum != atomsfromcom[iatom].atomnum ) {
	cout<<t<< "atomname or charge differs (vs com) "
	    <<  atoms[iatom].atomnum << " " << atomsfromcom[iatom].atomnum
	    <<endl;
	cout <<t<< " at atom "<< iatom+1 <<endl;
	exit(ERR_CODE);
      }
    }
  }
  
  return 0;
}


/************************************************/

int  main_check_consistency_basisset(
				     /* input */
              vector<g03basis> & basisset, 
	      vector<g03basis> & basissetfromfchk
				     /* no output */
	      )
{
  const static char *t="main_check_consistency_basisset: ";
  const static double eps=1.0e-5;
  
#if 0
  cout << "-----------------------------"<<endl;
  cout << "basis set" <<endl;
#endif
  vector<g03basis>::iterator vbas2=basissetfromfchk.begin();
  for ( vector<g03basis>::iterator vbas = basisset.begin();
        vbas != basisset.end(); vbas++) {
#if 0
    cout <<"atom="<< vbas->id << " " << vbas->range  << " "
         <<     vbas->bas_L.size() <<endl;
#endif
    if ( vbas->id != vbas2->id ||  vbas->bas_L.size() != vbas2->bas_L.size() ) {
        cout << t << " basisset error, basis differs at " <<  vbas->id  <<endl;
	cout <<"-----------------------------"<<endl;
	vbas->print(cout);
	cout <<"-----------------------------"<<endl;
	vbas2->print(cout);
        exit(ERR_CODE);
    }
    vector<g03basis_L>::iterator vbas2_L = vbas2->bas_L.begin();
    for ( vector<g03basis_L>::iterator vbas_L = vbas->bas_L.begin();
          vbas_L != vbas->bas_L.end(); vbas_L++) {
#if 0
      cout << vbas_L->orbtype << " " << vbas_L->ncoeff << " "
           << vbas_L->weight << " " << vbas_L->unknown << endl;
#endif
      if ( vbas_L->orbtype != vbas2_L->orbtype || 
               vbas_L->ncoeff !=  vbas2_L->ncoeff ||
                vbas_L->weight != vbas2_L->weight )  {
           cout << t << " basisset error, basis differs at " <<  vbas->id  <<endl;
	   cout <<"-----------------------------"<<endl;
	   vbas->print(cout);
	   cout <<"-----------------------------"<<endl;
	   vbas2->print(cout);
           exit(ERR_CODE);
      }
      if ( vbas_L->orbtype=="SP" ) {
        for (int i=0;i<vbas_L->expn.size();i++) {
#if 0
          cout << "   " <<vbas_L->expn[i] << " " << vbas_L->coeff1[i] << " "
               << vbas_L->coeff2[i] <<endl;
#endif
          //if ( vbas_L->expn[i] != vbas2_L->expn[i] ) { 
          if ( fabs( (vbas_L->expn[i] - vbas2_L->expn[i])/vbas2_L->expn[i] ) > eps ) { 
             cout << t << " basisset error, expn differs " <<endl;
	     cout << vbas_L->expn[i]  << " " << vbas2_L->expn[i]  << endl;
	     cout <<"-----------------------------"<<endl;
	     vbas->print(cout);
	     cout <<"-----------------------------"<<endl;
	     vbas2->print(cout);
             exit(ERR_CODE);
          }
        }
      } else {
        for (int i=0;i<vbas_L->expn.size();i++) {
#if 0
          cout << "   " <<  vbas_L->expn[i] << " " << vbas_L->coeff1[i] << endl;
#endif
          //if ( vbas_L->expn[i] != vbas2_L->expn[i] ) {
	  if ( fabs( (vbas_L->expn[i] - vbas2_L->expn[i] )/vbas2_L->expn[i] ) >eps ) {
             cout << t << " basisset error, expn differs " <<endl;
	     cout << vbas_L->expn[i] << " " << vbas2_L->expn[i] <<endl;
	     cout <<"-----------------------------"<<endl;
	     vbas->print(cout);
	     cout <<"-----------------------------"<<endl;
	     vbas2->print(cout);
             exit(ERR_CODE);
          }

        }
      }
      vbas2_L++;
    }
    vbas2++;
  }
#if 0
  cout << "---------------------------"<<endl;
#endif
  return 0;


}

int main_g03psL_comparesize_and_copy(
           /* input */ vector<g03PSP_L> & psLcom,
           /* output */ vector<g03PSP_L> & psL )
{
   const char *t = "main_g03psL_comparesize_and_copy:"; 
   if ( psLcom.size()!= psL.size() ) {
      cout << t<< "size is diferent " << psLcom.size() << " " << psL.size() <<endl;
      return 1;
   }
   vector<g03PSP_L>::iterator itcom=psLcom.begin(); 
   for ( vector<g03PSP_L>::iterator it= psL.begin(); it!=psL.end(); it++) {
      if ( it->L != itcom->L  || it->refL != itcom->refL ) {
          cout << t<< "L or refL is different." << it->L << " " << it->refL <<endl;
          return 1;
      }   
      //cout << t << "copy, L,refL="<< it->L << " " << it->refL <<endl;
      it->power=itcom->power;
      it->expn =itcom->expn;
      it->coeff=itcom->coeff;

      itcom++; 
   }

   return 0; // normal exit 
}

int  main_check_and_make_ps(/* input */ 
              vector<g03PSP>  & g03psfromcom,
              /* input output */ 
              vector<g03PSP>  & g03ps )
{
  static char *t = "main_check_make_ps:";
  int flag=0;// error flag
  //
  // find 0.0 in the coefficient and exponent
  //
  for ( vector< g03PSP  > ::iterator vvpsp= g03ps.begin();
        vvpsp != g03ps.end(); vvpsp++ ) {
    //cout <<"atom "<< vvpsp->id << " " << vvpsp->atomnum << " " << vvpsp->valence <<endl;
    vector<g03PSP_L> psL = vvpsp->psL;
    //cout << "num of PS " << psL.size()<<endl;;
    for ( vector<g03PSP_L>::iterator vpsp=psL.begin();
          vpsp!=psL.end(); vpsp++ ) {
      //cout << "L=" << vpsp->L << " refL=" << vpsp->refL <<endl;
      for (int i=0;i< vpsp->power.size(); i++) {
         // !!! job !!!, look for 0.0 
        if (vpsp->expn[i]==0.0 || vpsp->coeff[i]==0 ) { 
            flag=1;
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

   cout << "-> some parts of the PS in the logfile is broken. "; 
   if ( g03psfromcom.empty() ) {
        cout << "But no com file, or insufficient data in the com file. aborted"<<endl;
        exit(ERR_CODE);
   }

      cout << "I compensate them from the com file." <<endl;


   // g03ps
   //cout << "original"<<endl;
   //g03_print_ps( g03ps );

   //cout << "from com"<<endl;
   // g03psfromcom
   //g03_print_ps( g03psfromcom );

   // g03ps 
   for ( vector<g03PSP>::iterator it=g03ps.begin(); it!=g03ps.end(); it++) {
       int id = it->id;
       int atomnum = it->atomnum;
       int refid=it->refid;
       //cout << "g03ps " << it->id << " " <<it->atomnum<< " " << it->refid <<endl;
       if (refid!=id)  continue; 
       // find atomonum in g03psfromcom 
       int found=0;
       vector<g03PSP>::iterator itcom;
       for ( itcom=g03psfromcom.begin(); itcom!=g03psfromcom.end(); itcom++ ) {
           if ( atomnum == itcom->atomnum ) { found=1; break; }
       }
       if (found) {
          // copy
          if ( main_g03psL_comparesize_and_copy(itcom->psL,it->psL) ) {
              // =1 , if error occured 
              cout << t << " error in id,atomnum="<<it->id<< " " << it->atomnum <<endl;
              exit(ERR_CODE);
          }
       }
       else {
          cout << t << " atomnum="<<atomnum<<" not found"<<endl;
          exit(ERR_CODE);
       }
   }

   //cout << "compensated"<<endl;
   //g03_print_ps( g03ps );
   //
   return 0;
}

int main(int argc, char ** argv)
{
	const char *t="main: ";

   char *comfile=NULL, *logfile=NULL, *fchkfile=NULL;
   char *outfile=NULL;
   char *overwritestr=NULL;
   char *print_gamess=NULL;
   main_analyze_option(/* input */argc,argv,
         /* output */
         &comfile, &logfile,&fchkfile, &outfile,&overwritestr,  &print_gamess);

   int opt_overwrite=0;
   if ( strcmp(overwritestr,"y")==0 ) opt_overwrite=1;
   //cout << "overwrite="<< opt_overwrite <<endl;
   int opt_printgamess=0;
   if ( strcmp(print_gamess,"y")==0 ) opt_printgamess=1;


  string logfilename=logfile; 
  string fchkfilename=fchkfile;
  string outfilebasename = outfile;
  string comfilename ;
  if (comfile) comfilename = comfile;


  //---------- end of argv ------------------


  vector<string>  orbname;
  string scftype;
  double eref;
  int nup, ndown, multiplicity, dtype, ftype;
  int pbcsystem=0; // periodic boundary condition or not

  vector<vector<double> > alphamo, betamo;
  vector<g03Atom> g03atoms;

  // fchk
  //
  // alphamo, beamo, nup, ndown, multiplicity
  // atom and geometry (use this!)
  // basisset

  vector<g03basis> basissetfromfchk;
  vector< vector<double> > pbcvector;
  vector<double> origin(3);
  vector<double> kpoint(3); 
  for (int i=0;i<3;i++) { origin[i]=0.0; kpoint[i]=0.0; }

  main_read_g03fchkfile( /* input */ fchkfilename, 
		  /* output */ 
           g03atoms, pbcvector,
	   basissetfromfchk, alphamo,betamo,nup,ndown,multiplicity, eref,
	   dtype,ftype);
  //cout << "size of vector ="<< alphamo.size() << " " << betamo.size() <<endl;
  // error check
#if 0
  if ( dtype==1 || ftype==1 ) {
	  cout << "Error: now this program accepts only 6D and 10F" <<endl
		  << "But d/f shells are " << dtypestr[dtype] << " " << ftypestr[ftype] <<endl;
	  exit(ERR_CODE);
  }
#endif
	  
  if (! pbcvector.empty()) pbcsystem=1;

  vector<g03basis> basisset;
  vector<g03PSP> g03ps;
  vector<g03Atom> g03atomsfromlog ;

  // log
  //
  // basis type ( 6d, 10f ... )
  // atom and geometry
  // basisset  (use this!) (this basis set has move accuracy than those in the fchk file)
  // PS
  main_read_g03logfile(/* input */logfilename,  pbcsystem, 
		  /* output */ 
	       g03atomsfromlog, basisset, g03ps);


  vector<g03Atom> g03atomsfromcom; 
  vector<g03PSP>  g03psfromcom;

  // com
  //
  // atoms (only symbols)
  // PS 

  if (! comfilename.empty()) {
	  // atoms are only to check consistency 
	  // use this ps if the ps in the log is broken
      main_read_g03comfile(/* input */ comfilename, 
               /* output */ g03atomsfromcom, g03psfromcom ); 
  }

  cout << "--- have read files ---" <<endl;


  // check consistency
  //
  // g03atoms(from fchk), g03atomsfromlog, g03atomsfromcom
  main_check_consistency_g03atoms(g03atoms,g03atomsfromlog,g03atomsfromcom); 

  // basissetfromfchk, basisset(from log)
  main_check_consistency_basisset(basisset, basissetfromfchk);
  

  // PS  may be broken in the log file
  // fix broken PS
  //
  main_check_and_make_ps(/* input */  g03psfromcom,
               /* output */ g03ps );

#if 0
  /* check whether the basis set is acceptable to qwalk or not 
   *  6D or 10F
   */
  main_check_basissettype(
        /* input */
           basisset, basis_dftype , logfilename
        /* no output */ );
#endif

  // reorder f orbital, if there are
  g03_orbital_reorder(/* input */ basisset, dtype, ftype, 
		  /* input and output */ alphamo, betamo);


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


#if 1
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
		       g03atoms,
		       basisset,
		       g03ps,
		       alphamo,
		       betamo,
		       filename
		       /* no output */);
    cout << "output: " << filename <<endl;
  }

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

  }  // opt_printgamess 

#endif

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

   vector <Gaussian_basis_set>  gamessbasis;
   vector < vector < double> >  moCoeff;
   vector <Gaussian_pseudo_writer>  gamesspseudo;
   vector <Center>  centers;
   vector <int>  gamessnbasis;
   vector <Atom>  gamessatoms; 

   g03_change_form(
                /* input */
                 nup,  ndown,
                 scftype,
                 g03atoms,
                 basisset,
                 g03ps,
                 alphamo,
                 betamo,
                /* output */
                 slwriter,
                 gamessatoms,
                 gamessbasis,
                 gamesspseudo,
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


