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
#ifndef G03TOOLS_H_INCLUDED

#define G03TOOLS_H_INCLUDED

#include <string>
#include <vector>
#include <iostream>

#include "converter.h"
#include "wf_writer.h"
#include "basis_writer.h"
#include "Pseudo_writer.h"



#define ERR_CODE 10

#define TAB '\t'

#define ANG 0.5291772108

const static int nlabel_L=7;
const static char *label_Lstr[nlabel_L]={"S", "P", "D", "F",   "G", "I", "J" };
const static char *dtypestr[]={"5D","6D"};
const static char *ftypestr[]={"7F","10F"};



static const char * element_lookup_caps_null[]= {"NULL", "H", "HE",
"LI", "BE", "B", "C", "N", "O", "F", "NE",
"NA", "MG", "AL", "SI", "P", "S", "CL", "AR",
"K", "CA", "SC", "TI", "V", "CR", "MN", "FE",
     "CO", "NI", "CU", "ZN", "GA", "GE", "AS", "SE", "BR", "KR",
"RB", "SR", "Y", "ZR", "NB", "MO", "TC", "RU", "RH", "PD",
      "AG", "CD", "IN", "SN", "SB", "TE", "I", "XE",
"CS", "BA", "LA", "CE", "PR", "ND", "PM", "SM", "EU", "GD", "TB",
      "DY", "HO", "ER", "TM", "YB", "LU", "HF", "TA", "W", "RE",
      "OS", "IR", "PT", "AU", "HG", "TL", "PB", "BI", "PO", "AT", "RN",NULL  };


int strlen_trim(const char *c);

double atod(const char *s);

class g03Atom {
  public: 
  std::string name;
  int id;
  int atomnum;  /* all electron,   Atoms.charge=valence */
  std::vector<double> pos;

  void clear() { name=""; id=0; atomnum=0; pos.clear(); }
  void print_atom(std::ostream & os);
};


class g03basis_L {
 public:
  std::string orbtype;
  int ncoeff;
  double weight,unknown;
  std::vector<double> expn,coeff1,coeff2;
  g03basis_L() { clear(); }
  void clear() {
    orbtype=""; ncoeff=0; weight=0.0; unknown=0.0;
    expn.clear(); coeff1.clear(); coeff2.clear();
  }

};
class g03basis {
 public:
  int id,range;
  std::vector<g03basis_L>  bas_L;
  g03basis() { clear(); }
  void clear() { id=0; range=0; bas_L.clear(); }
    int print(std::ostream &os);

};

class g03basiskind {
  public:
	  std::vector<int>  ids;
	  std::vector<g03basis_L>  bas_L;
	  g03basiskind() { clear(); }
	  void clear() { ids.clear(); bas_L.clear(); }
          int print(std::ostream &os);

};

class g03PSP_L {
 public:
  int L, refL;
  std::vector<int>  power;
  std::vector<double> expn;
  std::vector<double> coeff;
  g03PSP_L() { clear(); }
  void clear() {L=-1;refL=-1; power.clear(); expn.clear(); coeff.clear(); }
};

class g03PSP {
 public:
  int id,atomnum,valence,refid;
  std::vector<g03PSP_L> psL;
  g03PSP() { clear(); }
  void clear() { id=0; atomnum=0; valence=0; psL.clear(); }
};

int strlen_trim(/* input */const char *c);

double atod(/* input */const char *s);

int findlmax( std::vector<g03PSP_L> & psL );

const char *L2Lstr(int L);

int g03_maxangularMomentum(/* input */ std::string & L);

int g03_basissizeL(/*input*/ std::string & L);

int label_to_L(/* input */ std::string & label,
		                 /* output */ int & refL);

double cutoff_divider(std::vector<Gaussian_basis_set> &basis,  std::vector< std::vector<double> > & latvec);


int g03_print_ps(/* input */
		           std::vector< g03PSP  > pseudos);


int g03_change_form(
                /* input */
                int nup, int ndown,
                std::string & scftype,
                std::vector <g03Atom> & g03atoms,
                std::vector<g03basis> & g03basisset,
                std::vector<g03PSP> & g03psp,
                std::vector< std::vector< double> > & alphamo,
                std::vector< std::vector< double> > & betamo,
                /* output */
                Slat_wf_writer & slwriter,
                std::vector <Atom>  & atoms,
                std::vector <Gaussian_basis_set> & basis,
                std::vector <Gaussian_pseudo_writer> & pseudo,
                std::vector < std::vector < double> > & moCoeff,
                std::vector <Center> & centers,
                std::vector <int> & basisidx
                ); 

int exist_file(std::string & filename);

int string_upper(std::string &str);

int find_atomname(std::vector<g03Atom> &atomlist, std::string& key);


int atomname2number(std::string &atom);

void g03_print_basisetkind(std::vector<g03basiskind> &baskind);

void g03_print_basiset(std::vector<g03basis> &baskind);

int g03_baskindset2basset(/* input */
	                std::vector<g03basiskind> &baskindset,
	                /* output */
	                std::vector<g03basis> &basisset)
;

int g03_process_sharpheader( 
		        /* input */
          std::vector<std::string> &lines,
          /* output */
            int &found_gen, int &found_pseudo);

int g03_orbital_reorder(/* input */
                std::vector<g03basis> & basisset, int dtype, int ftype, 
                   /* input & output */
                std::vector< std::vector<double> > & alphamo,
                std::vector< std::vector<double> > & betamo
                );


int basis_uniq_id(/* input */
		  int bas,  std::vector<Atom> atoms);


int ps_uniq_id(/* input */
	          int ips ,  std::vector<g03PSP> pseudo); 




int lastkey_is_spc(/* input */ std::string &str);

int lastkey_is_spc_substr20(/* input */ std::string &str);

int filekeyinput_get_string(
        /* input */
       std::string & filename,
       std::string & keyin,
       /* output */
       std::string & str);

int filekeyinput_get_stringv(
       /* input */
       std::string & filename, 
       std::string & keyin, int nskip, 
       int nread, std::string & keyout, int  (*func)(std::string &), 
       /* output */
       std::vector<std::string> & str );


int filekeyinput_get_stringv_last(
       /* input */
       std::string & filename,
       std::string & keyin, int nskip,
       int nread, std::string & keyout, int  (*func)(std::string &),
       /* output */
       std::vector<std::string> & str );




int g03_com_processheader(
        /* input */
          std::vector<std::string> &lines,
        /* output */
          int &found_gen, int &found_pseudo, std::string &dftype)
;

int g03_com_processatomgeometry(/* input */
                std::vector<std::string> lines,
                /* output */
                std::vector<g03Atom> &g03atoms)
;

int g03_com_processPS( /* input */
          std::vector<std::string> lines,
          std::vector<g03Atom> &g03atoms,
          /* output */
          std::vector<g03PSP> &g03ps)
;





 int g03_fchk_read_i(/* input */ std::string & filename,std::string & key,
                 /* output */ int &i)
;
int g03_fchk_read_r(/* input */ std::string & filename,std::string & key,
                 /* output */ double &r)
;

 int g03_fchk_read_iv(/* input */ std::string & filename,std::string & key,
                          int n,
                          /* output */int & nout,std::vector<int> &iv)
;

 int g03_fchk_read_rv(/* input */std::string & filename,std::string & key,
                          int n,
                          /* output */int & nout,std::vector<double> &dv)
;

int g03_fchk_readbasisset(/* input */ std::string & fchkfilename,
                                /* output */
                        std::vector<g03basis> & basisset )
;







int g03_log_analyze_orbname(
               /* input */
                std::string & logfilename ,
               /* output */
               std::vector<std::string> orbname)
;

int g03_log_analyze_basissetfromcards(/* input */
                std::vector<std::string> & basissetstring,
                /* output */
                std::vector<g03basis> & basisset )
;


int g03_log_analyze_basisset(/* input */
                std::vector<std::string> & basissetstring,
                /* output */
                std::vector<g03basis> & basisset )
;

int g03_log_analyze_ps(/* input */
                std::vector<std::string> & psstring,
                /* output */
        std::vector<g03PSP> & pseudos        )
;

int g03_log_readatoms( /*input*/ std::vector<std::string> & strv,
		                /*output*/ std::vector<g03Atom> & atoms )
;



int g03gamess_make_inp(/* input */
                int multiplicity,
                std::string &scftype,
                std::vector<g03Atom> & atoms,
                std::vector<g03basis>  & basisset,
                std::vector<g03PSP> &ecp,
                std::vector< std::vector<double> > & alphamo,
                std::vector< std::vector<double> > & betamo,
                std::string & filename
                /* no output */)
;

int g03gamess_make_out(
                /* input */
               std::string & outfilename,
               std::vector<g03Atom> & atoms,
               std::string & scftype,
               double eref,
               std::vector <g03basis> & basis,
               std::vector <g03PSP>  & psp,
               int nup, int ndown
               /* output */
               )
;

int g03gamess_make_pun(
        /* input */
                std::string & punfilename,
                std::vector <g03Atom> & atoms,
                std::vector <g03basis> & basis,
                std::vector< std::vector<double> > & alphamo,
                std::vector< std::vector<double> > & betamo
                )
;



#endif

