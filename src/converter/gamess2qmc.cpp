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
#include "converter.h"
#include <iostream>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>
#include <fstream>
#include <cmath>
#include "basis_writer.h"
#include "Pseudo_writer.h"
#include "wf_writer.h"
using namespace std;
void read_gamess_punch(int & add_norb, string & outputfilename,
		       string & punchfilename,
                       vector <Atom> & atoms,
                       Slat_wf_writer & slwriter,
                       vector <Gaussian_basis_set> & basis,
                       vector < vector < double> > & moCoeff);

void read_gamess_psp(istream & is,
                     vector <Atom> & atoms,
                     vector <Gaussian_pseudo_writer> & pseudo,
                     vector <double> & pspremoved);

void read_gamess_output(string & filename,
                        vector <Atom> & atoms,
                        Slat_wf_writer & slwriter, double & eref,
			vector <double> & electric_field,
			int & vorb);

void mo_analysis(vector <Atom> & atoms,
                 vector <Gaussian_basis_set> & basis,
                 vector < vector < double> > & moCoeff);



void usage(const char * name) {
  cout << "usage: " << name <<   " <options> <output> " << endl;
  cout << "Where options can be: \n";
  cout << "-compare_punch \n";
  cout << "-o          Base name for your run case\n";
  cout << "-virtual    Number of virtual orbitals to be read (Default 3)\n";
  cout << "-bp      Do not write boilerplate input files \n";
  cout << "-mo_analysis   Print largest MO coefficients into file mo_analysis \n";
  exit(1);
}

//######################################################################

int main(int argc, char ** argv) {
  string infilename;
  string outputname;
  string old_punchfile;
  int vorb=10; //default number of virtual orbs 
  int write_boilerplate=0;
  int do_mo_analysis=0;
  
  for(int i=1; i< argc-1; i++) {
    if(!strcmp(argv[i], "-o") && argc > i+1) {
      outputname=argv[++i];
    }
    else if(!strcmp(argv[i], "-compare_punch") && argc > i+1) {
            old_punchfile=argv[++i];
    }
    else if(!strcmp(argv[i], "-virtual") && argc > i+1) {
            vorb=int(atof(argv[++i]));
    }
    else if(!strcmp(argv[i], "-bp"))
      write_boilerplate=1;
    else if(!strcmp(argv[i], "-mo_analysis"))
      do_mo_analysis=1;
    else {
    cout << "Didn't understand option " << argv[i]
         << endl;
         usage(argv[0]);
    }
  }
  

  if(argc >= 2) {
    infilename=argv[argc-1];
  }
  else { usage(argv[0]); }

  if(outputname == "") {
    outputname=infilename;
  }


  ///////////////////////////////////////////////////////

  //Assume infilename is the root, and we have .inp, .pun, .out
  //Get atomic coordinates, alpha/beta orbitals, reference energy
  //from output file
  string outputfilename=infilename+".out";
  //also try .log if .out doesn't exist
  ifstream test_in(outputfilename.c_str());
  if(!test_in) {
    outputfilename=infilename+".log";
    test_in.close();
    test_in.clear();
    test_in.open(outputfilename.c_str());
    if(!test_in) { 
      cerr << "Couldn't find a .log or .out file" << endl; 
      exit(1); 
    }
  }
  test_in.close(); test_in.clear();
  
  vector <Atom> atoms;
  vector <Gaussian_pseudo_writer > pseudo;

  vector <double> electric_field;

  Slat_wf_writer slwriter;
  slwriter.write_centers=false;
  slwriter.mo_matrix_type="CUTOFF_MO";

  int nelectrons;
  double eref;
  read_gamess_output(outputfilename, atoms, slwriter, eref,
		     electric_field,vorb);

  ifstream is;
  is.open(outputfilename.c_str());
  vector <double> pspremoved;
  read_gamess_psp(is, atoms, pseudo, pspremoved);
  is.close();


  //Adjust the atomic charges for the pseudopotentials
  int natoms=atoms.size();
  int npseud=pseudo.size();
  for(int at=0; at < natoms; at++) {
    for(int psp=0; psp < npseud; psp++) {
      if(atoms[at].name==pseudo[psp].label) {
        atoms[at].charge-=pspremoved[psp];
        slwriter.nup-=int(pspremoved[psp]/2);
        slwriter.ndown-=int(pspremoved[psp]/2);
      }
    }
  }
  nelectrons=slwriter.nup+slwriter.ndown;



  vector < vector < double> > moCoeff;
  vector < Gaussian_basis_set > basis;

  
  vector < vector <double> >oldMOCoeff;
  if(old_punchfile !="") {
    read_gamess_punch(vorb, outputfilename, 
		      old_punchfile, atoms, slwriter, basis, oldMOCoeff);
  }
  string punchfilename=infilename+".pun";
  
  test_in.open(punchfilename.c_str());
  if(!test_in) {
    punchfilename=infilename+".dat";
    test_in.close();
    test_in.clear();
    test_in.open(punchfilename.c_str());
    if(!test_in) { 
      cerr << "Couldn't find a .dat or .pun file" << endl; 
      exit(1); 
    }
  }
  read_gamess_punch(vorb, outputfilename, 
		    punchfilename, atoms, slwriter, basis, moCoeff);
  if(do_mo_analysis)
    mo_analysis(atoms, basis,moCoeff);

 

  vector < Center> centers;
  vector <int> nbasis;
  centers.resize(atoms.size());
  nbasis.resize(natoms);
  for(int at=0; at < natoms; at++) {
    for(int i=0; i< 3; i++) centers[at].pos[i]=atoms[at].pos[i];
    centers[at].equiv_atom=at;
    centers[at].name=atoms[at].name;
    nbasis[at]=basis[atoms[at].basis].nfunc();
  }




  if(old_punchfile != "") {
    int nspin=1;

    if(slwriter.calctype=="UHF") nspin=2;
    for(int s=0; s< nspin; s++) {
      vector <int> compare_list;
      if(slwriter.calctype=="RHF" || slwriter.calctype=="ROHF") {
        int nmo=max(slwriter.nup, slwriter.ndown);
        for(int i=0; i< nmo; i++) {
          compare_list.push_back(i);
        }
      }
     else {
      if(s==0) {
        for(int i=0; i< slwriter.nup; i++) {
          compare_list.push_back(i);
        }
      }
      if(s==1) {
        for(int i=slwriter.spin_dwn_start;
            i < slwriter.spin_dwn_start+slwriter.ndown; i++) {
          compare_list.push_back(i);
        }
       }
     }

     compare_mo(oldMOCoeff, moCoeff, compare_list);
    }
  }



  //-------------------------------
  //print out the qmc input file
  
  string orboutname=outputname+".orb";
  slwriter.orbname=orboutname;
  string basisoutname=outputname+".basis";
  slwriter.basisname=basisoutname;

  ofstream orbout(orboutname.c_str());
  print_orbitals(orbout, centers, nbasis, moCoeff);
  orbout.close();

  ofstream basisout(basisoutname.c_str());
  int nbas=basis.size();
  for(int bas=0; bas < nbas; bas++) {
    basisout << "BASIS { \n";
    basis[bas].print_basis(basisout);
    basisout << "}\n\n\n";
  }
  basisout.close();

  string slateroutname=outputname+".slater";
  ofstream slaterout(slateroutname.c_str());
  slwriter.print_wavefunction(slaterout);
  slaterout.close();


  //---------------------------------Jastrow2 output
  string jast2outname=outputname+".jast2";
  
  double basis_cutoff=7.5; //arbitrary cutoff
  Jastrow2_wf_writer jast2writer;
  jast2writer.set_atoms(atoms);
  
  ofstream jast2out(jast2outname.c_str());
  print_std_jastrow2(jast2writer, jast2out, basis_cutoff);
  jast2out.close();
  

  string jast3outname=outputname+".jast3";
  ofstream jast3out(jast3outname.c_str());
  vector<string> unique_atoms;
  find_unique_atoms(atoms, unique_atoms);
  print_3b_jastrow2(jast3out,unique_atoms,basis_cutoff);
  jast3out.close();
  
  //--------------------------------------System output
  
  string sysoutname=outputname+".sys";
  ofstream sysout(sysoutname.c_str());
  sysout << "SYSTEM { MOLECULE \n";
  sysout << "  NSPIN { " << slwriter.nup << "  "
         << slwriter.ndown << " } \n";
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

  int npsp=pseudo.size();
  for(int psp=0; psp < npsp; psp++) {
    pseudo[psp].print_pseudo(sysout);
  }
  sysout.close();
  
  if(write_boilerplate) {
    string hfoutname=outputname+".hf";
    ofstream hfout(hfoutname.c_str());
    print_vmc_section(hfout, outputname, eref);
    hfout << "\n\n";
    hfout << "INCLUDE " << sysoutname << "  \n";
    hfout << "TRIALFUNC { INCLUDE " << slateroutname << "}\n\n";
    hfout.close();
    
    string optoutname=outputname+".opt";
    ofstream optout(optoutname.c_str());
    print_opt_section(optout, outputname, eref);
    optout << "\n\n";
    optout << "INCLUDE " << sysoutname << " \n";
    optout << "TRIALFUNC { \n  SLATER-JASTROW \n"
           << "  WF1 { INCLUDE " << slateroutname << " } \n"
           << "  WF2 { INCLUDE " << jast2outname   << " } \n"
           << "}\n\n";
    optout.close();
  }
  
  //cout << "end " << endl;
  return 0;
}

//######################################################################

void read_gamess_output(string & outputfilename,
                        vector <Atom> & atoms,
                        Slat_wf_writer & slwriter,
                        double & eref,
			vector <double> & electric_field,
			int & vorb) {

  ifstream is(outputfilename.c_str());
  if(!is) {
    cout << "Couldn't open " << outputfilename << endl;
    exit(1);
  }

  string line;
  string space=" ";

  vector <string> words;

  slwriter.nup=-1;
  slwriter.ndown=-1;

  while(getline(is, line)) {
    split(line, space, words);
    //cout << line;
    //if(words.size() > 0) cout << " firstword " << words[0] << endl;
       
    //Atoms
    if(words.size() > 2 && words[0]=="CHARGE" && words[1]=="X" && words[2] == "Y") {
      Atom tempatom;
      while(getline(is, line)) {
        //if(line=="") break;
        words.clear();
        split(line, space, words);
        //cout << "atom: " << line ;
        if(words.size()==0) break;
        tempatom.name=words[0];
        tempatom.charge=atof(words[1].c_str());
        for(int i=0; i< 3; i++) {
          tempatom.pos[i]=atof(words[i+2].c_str());
        }
        atoms.push_back(tempatom);
      }
    }
    //If we did a geometry optimization, pick up the last one 
    //(just read them in until there aren't any more)
    if(words.size() > 3 && words[0]=="ATOM" && words[1]=="CHARGE" && words[2]=="X") { 
      atoms.clear();
      is.ignore(180,'\n');  //ignore the line of ------'s
      Atom tempatom;
      while(getline(is, line)) {
        words.clear();
        split(line, space, words);
        if(words.size()==0) break;
	if(words[0]=="COORDINATES") break; //Useful for symmetry reduced coordinates
        tempatom.name=words[0];
	tempatom.charge=atof(words[1].c_str());
        //GAMESS prints out this section in angstroms for some reason
        for(int i=0; i< 3; i++) { 
          tempatom.pos[i]=atof(words[i+2].c_str())/0.529177249;
        }
        atoms.push_back(tempatom);
      }
      
    }      
    //---done atoms
    else if(words.size() > 4 && words[2]=="OCCUPIED" && words[3] == "ORBITALS") {
      if(words[4] == "(ALPHA)" && words[5] == "=") {
        //cout << line << endl;
        slwriter.nup=atoi(words[6].c_str());
      }
      if(words[4] == "(BETA" && words[6] == "=") {
        //cout << line << endl;
        slwriter.ndown=atoi(words[7].c_str());
      }

    }


    else if(words.size() > 4 && words[0] == "FINAL" && words[2] == "ENERGY") {
      eref=atof(words[4].c_str());
    }
    else if(words.size() > 0 && words[0]=="SCFTYP=RHF") {
      slwriter.calctype="RHF";
    }
    else if(words.size() > 0 && words[0]=="SCFTYP=ROHF") {
      slwriter.calctype="ROHF";
    }
    else if(words.size() > 0 && words[0]=="SCFTYP=UHF") {
      slwriter.calctype="UHF";
    }
    else if(words.size() > 0 && words[0]=="SCFTYP=GVB") {
      slwriter.calctype="GVB";
    }
    else if(words.size() > 0 && words[0]=="$EFIELD") {
      electric_field.resize(3);
      getline(is,line); //break line
      getline(is,line);
      words.clear();
      split(line, space, words);
      assert(words.size() > 4);
      for(int d=0; d< 3; d++) 
        electric_field[d]=atof(words[d+1].c_str());
    }
      
    words.clear();
  }

  is.close();
  is.clear();


  if(atoms.size() == 0) {
    cout << "******WARNING*******  Couldn't find any atoms " << endl;
  }
  if(slwriter.calctype=="") {
    cout << "Couldn't find SCFTYP" << endl;
    exit(1);
  }
  if(slwriter.nup <0 ) {
    cout << "Couldn't find the number of alpha orbitals" << endl;
    exit(1);
  }
  if(slwriter.ndown < 0) {
    cout << "Couldn't find the number of beta orbitals" << endl;
    exit(1);
  }


  if(slwriter.calctype=="GVB") {
    is.open(outputfilename.c_str());
    if(!is) {
      cout << "Couldn't open " << outputfilename << endl;
      exit(1);
    }
    while(getline(is,line)) {
      split(line, space, words);
      if(words[0]=="PAIR" && words[1]=="INFORMATION") {
        //cout << "found " << endl;
        is.ignore(180, '\n'); //-----
        is.ignore(180, '\n'); //ORBITAL   CI COEFFICIENTS ...
        is.ignore(180, '\n'); //PAIR   1 2 ...

        vector < vector <int> > orb(2);
        vector < vector <double> > coeff(2);
        int minorb=50000;
        int npair=0;
        while(getline(is, line)) {
          if(line=="") break;
          //cout << line << endl;
          words.clear();
          split(line, space, words);
          orb[0].push_back(atoi(words[1].c_str()));
          orb[1].push_back(atoi(words[2].c_str()));

          for(int s=0; s< 2; s++) {
            if(minorb > orb[s][npair])
              minorb=orb[s][npair];
          }

          coeff[0].push_back(atof(words[3].c_str()));
          coeff[1].push_back(atof(words[4].c_str()));
          npair++;
        }
        //cout << "minorb " << minorb << endl;
        vector <int> base_occ(minorb-1);
        for(int i=0; i< minorb-1; i++) {
          base_occ[i]=i+1;
        }

        vector <int> pow2(npair+1);
        pow2[0]=1;
        for(int i=1; i< npair+1; i++) {
          pow2[i]=2*pow2[i-1];
        }
        int ndet=pow2[npair];
        
        //cout << ndet << " determinants" << endl;
        vector <int> prom(npair);

        for(int d=0; d< ndet; d++) {
          //cout << "determinant " << d << "  : ";
          for(int i=0; i< npair; i++) {
            prom[i]=(d/pow2[i])%2;
            //cout << prom[i] << " ";
          }
          //cout << endl;
          slwriter.occ_up.push_back(base_occ);
          slwriter.occ_down.push_back(base_occ);
          int start=base_occ.size();
          double wt=1;
          for(int i=0; i< npair; i++) {
            slwriter.occ_up[d].push_back(orb[prom[i]][i]);
            slwriter.occ_down[d].push_back(orb[prom[i]][i]);
            wt*=coeff[prom[i]][i];
          }
          slwriter.detwt.push_back(wt);
        }
        vorb=max(vorb,npair);
	
        break;
      }
      words.clear();
    }

    is.close();
  }


}


//#####################################################################

void read_gamess_psp(istream & is,
                     vector <Atom> & atoms,
                     vector <Gaussian_pseudo_writer> & pseudo,
                     vector <double> & pspremoved)
{

  vector <double> emptyvector;
  vector <int> emptyintvector;
  vector <string> words;
  string line;
  string space=" ";
  while(getline(is, line)) {
    split(line, space, words);
    if(words[0]=="ECP" && words[1]=="POTENTIALS") {
      is.ignore(200, '\n'); //line of ---'s
      is.ignore(200, '\n'); //empty line
      while(getline(is, line)) {
        words.clear();
        split(line, space, words);
        if(words[0]=="THE" && words[2] == "RUN") break;

        //-------Read in a pseudopotential
        if(words[0]=="PARAMETERS" && words.size() > 7 ) {
          int atomplace=0;  //placement of the word ATOM on the line
          //Sometimes it's at 4, sometimes it's at 5
          if (words[4]=="ATOM") atomplace=5;
          else if(words[5]=="ATOM") atomplace=6;
          else {
            cout << "couldn't find ATOM in pseudo reading" << endl;
            exit(1);
          }
          int atomnum=atoi(words[atomplace].c_str())-1;
          //cout << "atom number " << atomnum << endl;


          //----equivalent atom
          if(words[atomplace+3]=="SAME") {
            int equivatom=atoi(words[atomplace+6].c_str())-1;
            if(atoms[atomnum].name != atoms[equivatom].name) {
              cout << "atoms " << atomnum << " and " << equivatom
              << " don't match, even though they should." << endl;
              exit(1);
            }
          }


          //---nonequivalent atom
          else {
            Gaussian_pseudo_writer pseudotmp;

            pseudotmp.label=atoms[atomnum].name;
            pspremoved.push_back(atof(words[atomplace+3].c_str()));
            //cout << "pspremoved " << words[atomplace+3] << endl;
            int nblocks=atoi(words[atomplace+6].c_str())+1;
            is.ignore(120, '\n'); //FOR L= COEFF etc..
            for(int block=0; block < nblocks; block++) {
              pseudotmp.exponents.push_back(emptyvector);
              pseudotmp.coefficients.push_back(emptyvector);
              pseudotmp.nvalue.push_back(emptyintvector);
              while(getline(is, line)) {
                words.clear();
                split(line, space, words);
                if(words.size()==0 || words[0]=="FOR") break;
            pseudotmp.coefficients[block].push_back(atof(words[1].c_str()));
                pseudotmp.nvalue[block].push_back(atoi(words[2].c_str())-2);
                pseudotmp.exponents[block].push_back(atof(words[3].c_str()));
              }
            }
            pseudo.push_back(pseudotmp);
          }


        }
        //--------------------------------

        //cout << line << endl;
      }

    }
    words.clear();
  }


  //QMC likes the local part at the end of the list, and GAMESS outputs
  //it at the beginning, so let's fix that..
  int npseud=pseudo.size();
  for(int ps=0; ps < npseud; ps++) {
    pseudo[ps].exponents.push_back(pseudo[ps].exponents[0]);
    pseudo[ps].nvalue.push_back(pseudo[ps].nvalue[0]);
    pseudo[ps].coefficients.push_back(pseudo[ps].coefficients[0]);
    pseudo[ps].exponents.erase(pseudo[ps].exponents.begin());
    pseudo[ps].coefficients.erase(pseudo[ps].coefficients.begin());
    pseudo[ps].nvalue.erase(pseudo[ps].nvalue.begin());
  }


}



//######################################################################

void read_gamess_punch(int & vorb,
		       string & outputfilename,
		       string & punchfilename,
                       vector <Atom> & atoms,
                       Slat_wf_writer & slwriter,
                       vector <Gaussian_basis_set> & basis,
                       vector < vector < double> > & moCoeff) {
  ifstream is(punchfilename.c_str());
  if(!is) {
    cout << "Couldn't open " << punchfilename << endl;
    exit(1);
  }
  string line;
  vector <string> words;
  string space=" ";
  vector <double> emptyvector;
  int natoms=atoms.size();

  moCoeff.clear();
  basis.clear();

  //Get the basis from the punch file
  while(getline(is, line)) {
    words.clear();
    split(line, space, words);

    if(!words.empty() && words[0] == "$DATA") {
      getline(is, line); //line after $DATA
      getline(is, line); //symmetry line
      words.clear();
      split(line, space, words);
      if(words[0] != "C1") {
        getline(is, line); //For non-C1 symmetry, there's an
        //extra line after the symmetry line
      }

      //Basis header
      while(getline(is, line) ) {
        words.clear(); split(line, space, words);
        if(!words.empty() ) {
          if(words[0]=="$END") break;
          //cout << "Label " << words[0] << endl;
          int unique=1;
          for(unsigned int i=0; i< basis.size(); i++) {
            if(basis[i].label == words[0]) {
              unique=0;
              break;
            }
          }
          if(unique) { //if unique, read the basis set
            Gaussian_basis_set tmpbasis;
            tmpbasis.label=words[0];

            while(getline(is, line)) {
              words.clear(); split(line, space, words);
              if(words.empty()) break;
              int haveLtype=0;
              if(words[0]=="D")
                tmpbasis.types.push_back("6D");
              else if(words[0]=="F")
                tmpbasis.types.push_back("10F");
              else if(words[0]=="G")
                tmpbasis.types.push_back("15G");
              else if(words[0]=="L") {
                tmpbasis.types.push_back("S");
                tmpbasis.types.push_back("P");
                haveLtype=1;
              }
              else
                tmpbasis.types.push_back(words[0]);
              
              int nexpansion=atoi(words[1].c_str());
              tmpbasis.exponents.push_back(emptyvector);
              tmpbasis.coefficients.push_back(emptyvector);

              
              int place=tmpbasis.types.size()-1;
              if (haveLtype) { 
                tmpbasis.exponents.push_back(emptyvector); 
                tmpbasis.coefficients.push_back(emptyvector);
                place--; //added two basis functions for the L-type AO
              } 
              for(int i=0; i< nexpansion; i++) {
                getline(is, line);
                words.clear(); split(line, space, words);
                tmpbasis.exponents[place].push_back(atof(words[1].c_str()));
                tmpbasis.coefficients[place].push_back(atof(words[2].c_str()));
                if (haveLtype) { 
                  assert(4==words.size());
                  tmpbasis.exponents[place+1].push_back(atof(words[1].c_str()));
                  tmpbasis.coefficients[place+1].push_back(atof(words[3].c_str()));
                }
              }
              //cout << line << endl;
            }
            basis.push_back(tmpbasis);
          }
          else {//if we've already read this one, skip to the next set
            while(getline(is, line)) {
              //cout << "skipping " << line << endl;
              words.clear(); split(line, space, words);
              if(words.empty()) break;
            }
          }
        }
      }

    }
  }
  is.close();
  is.clear();
  int totbasis=0;
  for(unsigned int at=0; at < atoms.size(); at++) {
    for(unsigned int bas=0; bas< basis.size() ; bas++) {
      if(basis[bas].label==atoms[at].name) {
        atoms[at].basis=bas;
        totbasis+=basis[bas].nfunc();
      }
    }
  }


  //Now get the MO coefficients
  const int norb=slwriter.nup+vorb;
  //Read in vorb more orbitals than the number of up
  //electrons, per spin channel if UHF (default is 3)
  is.open(punchfilename.c_str());

  string temp;
  int lastmo=-1;
  int totmospin=0;
  int moindex=-1;

  int totmos=totbasis;
  //cout << "totbasis " << totbasis << endl;

  ifstream inout(outputfilename.c_str());
  string search="TOTAL NUMBER OF MOS IN VARIATION SPACE";
  while(getline(inout,line)) { 
    if(line.find(search) < line.size()) { 
      vector <string> wrds;
      split(line,space, wrds);
      totmos=atoi(wrds[7].c_str());
      //cout << "set totmos to " << totmos << endl;
    }	  
  }

  inout.close();


  while(getline(is, line)) {
    if(line.find("$VEC") != string::npos) {
      //cout << "found vec " << endl;
      moCoeff.clear(); //we'll end up taking the last $VEC section
      lastmo=-1;
      totmospin=0;
      moindex=-1; 

      while(getline(is, line)) {
        words.clear(); split(line, space, words);
        if(words[0]=="$END") break;
        temp.assign(line.begin(), line.begin()+2);
        int currmo=atoi(temp.c_str());
        if(currmo != lastmo) {
          lastmo=currmo;
          totmospin++;

          //For RHF and ROHF, there shouldn't be more MO's than
          //basis functions
          if(totmospin > totbasis &&
              (slwriter.calctype=="RHF" || slwriter.calctype=="ROHF") ) {
            cout << "too many mo's: there are " << totbasis << " basis functions and "
                 << totmospin << " mo's so far"  << endl;
            cout << "line : " << line << endl;
          }

	  //cout << "totmospin " << totmospin 
	  //     << " currmo " << currmo << endl;
          //For UHF, the down channel MO's are listed after all the
          //up channel ones, so we start over when we find a '1' in
          //the MO field
          if(totmospin >=totmos && currmo==1
             && slwriter.calctype=="UHF") {
            cout << "starting spin down " << line <<  endl;
            slwriter.spin_dwn_start=moindex+1;
            totmospin=0;
          }

          //If we should read this molecular orbital, then push
          //an empty vector onto our coefficient stack
          if(totmospin <=norb) {
            moCoeff.push_back(emptyvector);
            moindex++;
            //cout << "moindex " << moindex << endl;
          }
        }
        //Push the MO coefficients onto moCoeff.
        if(totmospin <= norb) {
          for(int i=5; i < 66; i+=15) {
            if(line.size() >= i+15) {
              temp.assign(line.begin()+i, line.begin()+i+15);
              moCoeff[moindex].push_back(atof(temp.c_str()));
              //cout << "mo " << temp << endl;
            }
          }

        }

      }
    }
  }
  is.close();



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
        else {
          cout << "unknown type " << basis[bas].types[i] << endl;
          exit(1);
        }
        
      }
    }
    //cout <<" func " << func << endl;
  }

}


//######################################################################
void mo_analysis(vector <Atom> & atoms,
                 vector <Gaussian_basis_set> & basis,
                 vector < vector < double> > & moCoeff){

  
  ofstream an_out("mo_analysis");
  const double print_thresh=1e-1;
  int natoms=atoms.size();
  vector <string> pnames(3);
  pnames[0]="x   ";
  pnames[1]="y   ";
  pnames[2]="z   ";
  vector <string> dnames(6);
  dnames[0]="xx  ";
  dnames[1]="yy  ";
  dnames[2]="zz  ";
  dnames[3]="xy  ";
  dnames[4]="xz  ";
  dnames[5]="yz  ";
  vector <string> fnames(10);
  fnames[0]="xxx   ";
  fnames[1]="yyy   ";
  fnames[2]="zzz   ";
  fnames[3]="xxy   ";
  fnames[4]="xxz   ";
  fnames[5]="yyx   ";
  fnames[6]="yyz   ";
  fnames[7]="zzx   ";
  fnames[8]="zzy   ";
  fnames[9]="xyz   ";
  vector <string> gnames(15);
  gnames[0]="xxxx   ";
  gnames[1]="yyyy   ";
  gnames[2]="zzzz   ";
  gnames[3]="xxxy   ";
  gnames[4]="xxxz   ";
  gnames[5]="yyyx   ";
  gnames[6]="yyyz   ";
  gnames[7]="zzzx   ";
  gnames[8]="zzzy   ";
  gnames[9]="xxyy   ";
  gnames[10]="xxzz   ";
  gnames[11]="yyzz   ";
  gnames[12]="xxyz   ";
  gnames[13]="yyxz   ";
  gnames[14]="zzxy   ";
 

  int totmo2=moCoeff.size();
  for(int mo=0; mo < totmo2; mo++) {
    an_out << "\n----------------\n";
    an_out << "MO " << mo+1 << endl;
    int func=0;
    for(int at=0; at < natoms; at++) {
      int bas=atoms[at].basis;
      int nbasis=basis[bas].types.size();
      for(int i=0; i< nbasis; i++) {
        if(basis[bas].types[i] == "S") {
          if(fabs(moCoeff[mo][func]) > print_thresh) {
            an_out << atoms[at].name<< at  << "  S     " << moCoeff[mo][func]<< endl;
          }
          func++;
        }
        else if(basis[bas].types[i] == "P") {
          for(int j=0; j< 3; j++) {
            if(fabs(moCoeff[mo][func]) > print_thresh) {
              an_out << atoms[at].name << at << "  "  << "P" 
		     << pnames[j] << " " << moCoeff[mo][func]
		     << endl;
            }
            func++;
          }
        }
        else if(basis[bas].types[i] == "6D") {
          for(int j=0; j< 6; j++) {
            if(fabs(moCoeff[mo][func]) > print_thresh) {
              an_out << atoms[at].name << at << "  "   << "D" 
		     << dnames[j] << " " <<  moCoeff[mo][func]
		     << endl;
            }
            func++;
          }
        }
        else if(basis[bas].types[i] == "10F") {
          for(int j=0; j< 10; j++) {
            if(fabs(moCoeff[mo][func]) > print_thresh) {
              an_out << atoms[at].name << at << "  "   << "F" 
		     << fnames[j] << " " <<  moCoeff[mo][func]
		     << endl;
            }
            func++;
          }
        }
        else if(basis[bas].types[i] == "15G") {
          for(int j=0; j< 15; j++) {
            if(fabs(moCoeff[mo][func]) > print_thresh) {
              an_out << atoms[at].name << at << "  "   << "G" 
		     << gnames[j] << " " <<  moCoeff[mo][func]
		     << endl;
            }
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
  an_out.close();
}
