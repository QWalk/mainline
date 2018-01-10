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
#include "wf_writer.h"
#include "basis_writer.h"
#include "Pseudo_writer.h"
#include <algorithm>
#include <fstream>
#include <cstring>
#include "vecmath.h"
#include <string>
#include "elements.h"
#include <stdlib.h>
#include <stdio.h>
using namespace std;

void get_crystal_pseudo(istream & infile,
    vector <Gaussian_pseudo_writer> & pseudo);

void get_crystal_latvec(istream & infile,
    vector< vector < double > > & latvec);

void get_crystal_atoms(istream & infile,
    vector < Atom > & atoms);

void get_crystal_basis(istream & infile,
    vector <Gaussian_basis_set> & basis);

// for orbitals with real coefficients
void read_crystal_orbital(istream & is,
    string & fort10file,
    vector <Atom> & atoms,
    Slat_wf_writer & slwriter,
    vector <Gaussian_basis_set> & basis,
    vector <double> & origin,
    vector < vector <double> >& latvec,
    vector < vector < double > > & moCoeff,
    vector <double> & shift //amount to shift the atoms
    );
void read_kpt_eigenpos(istream & is,
    vector <string> & rkpoints,
    vector <long int> & reigen_start, 
    vector <string> & ckpoints,
    vector <long int> & ceigen_start); 
// for orbitals with complex coefficients
void read_crystal_orbital(istream & is,
    vector <Atom> & atoms,
    Slat_wf_writer & slwriter,
    vector <Gaussian_basis_set> & basis,
    vector <double> & origin,
    vector < vector <double> >& latvec,
    vector < vector < dcomplex > > & moCoeff,
    vector <double> & shift //amount to shift the atoms
    );

void read_crystal_orbital_all(istream & is,
    string & fort10file,
    vector <Atom> &atoms,
    Slat_wf_writer & slwriter,
    vector <Gaussian_basis_set> & basis,
    vector <double> & origin,
    vector < vector <double> >& latvec,
    vector < vector < vector < double > > > & moCoeff,
    vector <string> &kpoints, 
    vector < long int > & eigen_start, 
    vector <vector <double> > &kptCoord,
    vector <double> & shift,  //amount to shift the atoms
    vector < int > shifted, 
    vector < vector <int> > nshift
    );

// for orbitals with complex coefficients
void read_crystal_orbital_all(istream & is,
    vector <Atom> &atoms,
    Slat_wf_writer & slwriter,
    vector <Gaussian_basis_set> & basis,
    vector <double> & origin,
    vector < vector <double> >& latvec,
    vector < vector < vector < dcomplex > > > & moCoeff,
    vector <string> &kpoints,
    vector < long int > & eigen_start,
    vector < vector <double> > &kptCoord,
    vector <double> & shift, //amount to shift the atoms
    vector < int > shifted, 
    vector < vector <int> > nshift
    );



void usage(const char * name) {
  cout << "usage: " << name <<   " <options> <crystal output> " << endl;
  cout << "Where options can be: \n\n";
  cout << "-c          Work with complex orbitals. Only real are searched by default.\n\n";
  cout << "-fort10file Formatted fortran unit 10 from read10 utility.\n";
  cout << "            try to match the vectors in fort.10 to the ones in the output.\n";
  cout << "            It may or may not work correctly..\n\n";
  cout << "-o          Base name for your run case\n";
  exit(1);
}

int main(int argc, char ** argv) {

  bool had_error=false;
  bool cmplx=false;
  string infilename;
  string outputname, fort10file;
  vector <double> shift(3);
  for(int d=0; d< 3; d++) shift[d]=.125;

  for(int i=1; i< argc-1; i++) {
    if(!strcmp(argv[i], "-o") && argc > i+1) {
      outputname=argv[++i];
    }
    else if(!strcmp(argv[i],"-fort10file")) {
      if (i+1 < argc) {
        fort10file=argv[i+1];
        i++;
      }
      else {
        cout << "-fort10file needs an argument" << endl;
        exit(1);
      }
    }
    else if(!strcmp(argv[i],"-shift")){
      if(i+1 < argc) {
        for(int d=0; d< 3; d++)  { shift[d]=atof(argv[i+1]);
          i++;
        }
      }
      else { 
        cout << "-shiftz needs an argument" << endl;
        exit(1);
      }
    }
    else if(!strcmp(argv[i], "-noshift")) { 
      for(vector<double>::iterator i=shift.begin();
          i!= shift.end(); i++) *i=0.0;
    }
    else if(!strcmp(argv[i], "-c")) { 
      cmplx=true;
    }

    else {
      cout << "Didn't understand option " << argv[i] << endl;

      had_error=true;
    }
  }

  if(argc >= 2) {
    infilename=argv[argc-1];
  }
  else { had_error=true; }

  if(outputname == "") {
    outputname="qwalk";
    //outputname=infilename;
  }

  if(had_error) usage(argv[0]);

  if ( cmplx && ( fort10file != "" ) ) {
    cout << "-fort10file is not supported for complex orbitals (-c)." << endl;
    exit(1);
  }


  ///////////////////////////////////////////////////////////////

  vector <Atom> atoms;
  vector <Gaussian_pseudo_writer > pseudo;
  vector <Gaussian_basis_set > basis;

  Slat_wf_writer slwriter;
  slwriter.use_global_centers=true;
  slwriter.write_centers=false;
  if ( cmplx ) {
    slwriter.mo_matrix_type="CUTOFF_MO";
    slwriter.orbtype="CORBITALS";
  } else {
    slwriter.mo_matrix_type="CUTOFF_MO";
  }

  vector < vector < double > > latvec;


  ifstream infile(infilename.c_str());
  if(!infile) {
    cout << "couldn't open " << infilename << endl;
    exit(1);
  }
  get_crystal_latvec(infile, latvec);
  infile.close();
  infile.clear();

  infile.open(infilename.c_str());
  get_crystal_atoms(infile, atoms);
  infile.close();
  infile.clear();
  if(latvec.size() > 0) { 
    vector <double> atomshift;
    for(int i=0; i< 3; i++) atomshift.push_back(0);
    for(int i=0; i< 3; i++) {
      for(int j=0; j< 3; j++) {
        atomshift[j]+=shift[i]*latvec[i][j];
      }
    }
    for(int i=0; i < atoms.size(); i++) {
      atoms[i].pos=atoms[i].pos+atomshift;
    }
  }
  Shifter shiftobj; 
  vector < vector <int> > nshift;
  nshift.resize(atoms.size());
  vector <int> shifted(atoms.size()); 
  if(latvec.size() > 0) { 
    for (int at = 0; at<atoms.size(); at++) { 
      shifted[at]=shiftobj.enforcepbc(atoms[at].pos, latvec, nshift[at]);
      cout << "shifted " << at << "  " << shifted[at] << endl;
    }
  }

  infile.open(infilename.c_str());
  get_crystal_basis(infile, basis);
  infile.close();
  infile.clear();
  infile.open(infilename.c_str());
  get_crystal_pseudo(infile, pseudo);
  infile.close();
  infile.clear();

  int nelectrons=-1;
  double totspin=-1e8;
  infile.open(infilename.c_str());
  string calctype, testword;
  int nelectrons_up=-1, nelectrons_down=-1;
  string line, space=" ";
  vector <string> spl;
  while(true) {
    spl.clear();
    getline(infile,line);
    split(line,space,spl);
    if(spl.size() > 2 and spl[1]=="SCF" and spl[2]=="ENDED")  break;
    if(spl.size() > 2 and spl[0]=="TOTAL" and spl[1]=="ATOMIC" and spl[2]=="SPINS") { 
      spl.clear();
      while(true) { 
        getline(infile,line);
        if(line[3]=='T') break;
        split(line,space,spl);
      }
      totspin=0.0; 
      for(vector<string>::iterator i=spl.begin(); i!=spl.end(); i++) { 
        totspin+=atof(i->c_str());
      }
    }
    if(spl.size() > 4 and spl[0]=="N." and spl[2]=="ELECTRONS" and spl[4]=="CELL") { 
      nelectrons=atof(spl[5].c_str());
    }
    if(spl.size() > 4 and spl[0]=="TYPE" and spl[2]=="CALCULATION") {
      if(spl[4]=="RESTRICTED") {
        if(spl[5]=="CLOSED") calctype="RHF";
        else if(spl[5]=="OPEN") calctype="ROHF";
      }
      else if(spl[4]=="UNRESTRICTED") calctype="UHF";
      else { 
        cout << "Didn't understand calctype " << spl[4] << endl;
        exit(1);
      }
    }
  }
  slwriter.calctype=calctype;
  infile.close();

  if(atoms.size() == 0) {
    cout << "I couldn't find any atoms!!" << endl;
    exit(1);
  }

  //----------------------------------------------------------------------
  //Done parsing, now organize and link.

  //Figure out which atoms go with which basis sets.
  int natoms=atoms.size();
  int nbasis=basis.size();
  for(int at=0; at<natoms; at++) {
    int found_basis=0;
    for(int bas=0; bas< nbasis; bas++) {

      if(atoms[at].name == basis[bas].label) {
        atoms[at].basis=bas;
        found_basis=1;
        break;
      }
    }
    if(!found_basis) {
      cout << "Couldn't find basis for atom " << atoms[at].name
        << endl;
      exit(1);
    }
  }


  //We need to get the effective charges from the user, since
  //Crystal isn't nice enough to output them in any consistent way.
  int testelectrons=0;

  testelectrons=0;
  for(int at=0; at < natoms; at++) {
    int isUnique=1;
    for(unsigned int i=0; i< at; i++) {
      if(atoms[at].name == atoms[i].name) {
        isUnique=0;
        atoms[at].charge=atoms[i].charge;
        break;
      }
    }
    if(isUnique) {
      int npseudo=pseudo.size();
      bool found_in_pseudo=false;
      for(int j=0; j< npseudo; j++) { 
        if(atoms[at].name==pseudo[j].label) { 
          atoms[at].charge=pseudo[j].effcharge;
          found_in_pseudo=true;
          break;
        }
      }
      if(!found_in_pseudo) { 
        cout << "Please enter the effective charge of " << atoms[at].name
          << endl;
        cin >> atoms[at].charge;
        cout << "Thanks" << endl;
      }
    }
    testelectrons+=(int) atoms[at].charge;
  }
  if(testelectrons != nelectrons) {
    cout << "The total number of electrons should be " << nelectrons
      << ", \nbut the sum due to the effective charge on the ion is "
      << testelectrons << ".  This may be what you intended, but be careful!\n";
  }



  //Determine how many electrons are spin up or down if it's not RHF
  if(calctype=="RHF" ) {
    if(nelectrons % 2 != 0) {
      cout << "It seems like you're doing RHF, but there is an odd number of"
        << " electrons!  I'm confused." << endl;
      exit(1);
    }
    slwriter.nup=nelectrons/2;
    slwriter.ndown=nelectrons/2;
  }
  else if(calctype=="ROHF" || calctype=="UHF") {
    cout << "Detected a ROHF or UHF calculation." << endl;
    //    "  What is the spin state?(1=singlet, 2=doublet...)" << endl; 
    cout << "N_up-N_down is found to be " << totspin << endl;
    int spin=totspin+0.1;
    if( abs(totspin-spin) > 1e-4) { 
      cout << "WARNING: spin is not close to an integer!" << endl;
    }
    slwriter.nup=(nelectrons-spin)/2 + spin;
    slwriter.ndown=(nelectrons-spin)/2;
    if(slwriter.nup+slwriter.ndown != nelectrons) {
      cout << "problem..  nup and ndown are calculated to be "
        << slwriter.nup << "   " << slwriter.ndown
        << " but they don't add up to be " << nelectrons << endl;
      exit(1);
    }
  }
  else {
    cout << "Don't understand calculation type of "
      << calctype << endl;
    exit(1);
  }

  vector <double> origin(3);
  for(int i=0; i< 3; i++) origin[i]=0;

  vector <vector < vector < double > > > moCoeff;
  vector < vector < vector < dcomplex > > > CmoCoeff;
  infile.clear();
  infile.open(infilename.c_str());
  vector <string> rkpts;
  vector <string> ckpts; 
  vector < vector <double> > ckptCoord; 
  vector < vector <double> > rkptCoord; 
  vector <long int > ceigen_start, reigen_start;  
  read_kpt_eigenpos(infile, rkpts, reigen_start, ckpts, ceigen_start); 
  if(cmplx) 
    read_crystal_orbital_all(infile, atoms, slwriter, basis,
        origin, latvec, CmoCoeff, ckpts, ceigen_start, ckptCoord, shift, shifted, nshift); 
  infile.close();
  infile.clear();
  infile.open(infilename.c_str());
  read_crystal_orbital_all(infile, fort10file, atoms, slwriter, basis,
      origin, latvec, moCoeff, rkpts, reigen_start, rkptCoord, shift, shifted, nshift);
  infile.close();
  natoms=atoms.size();


  vector <Center> centers;
  centers.resize(atoms.size());
  for(int at=0; at < natoms; at++) {
    for(int i=0; i< 3; i++) centers[at].pos[i]=atoms[at].pos[i];
    centers[at].equiv_atom=at;
    centers[at].name=atoms[at].name;
  }

  //-------------------------------
  //print out the qmc input file

  string orboutname=outputname+".orb";
  slwriter.orbname=orboutname;
  string basisoutname=outputname+".basis";
  slwriter.basisname=basisoutname;

  vector <int> nbasis_list;
  nbasis_list.resize(natoms);
  for(int at=0; at < natoms; at++) {
    nbasis_list[at]=basis[atoms[at].basis].nfunc();
  }

  ofstream basisout(basisoutname.c_str());
  int nbas=basis.size();
  for(int bas=0; bas < nbas; bas++) {
    basisout << "BASIS { \n";
    basis[bas].print_basis(basisout);
    basisout << "}\n\n\n";
  }
  basisout.close();



  //--------------------------Jastrow 2 output 

  double basis_cutoff;
  if(latvec.size() >0) 
    basis_cutoff=find_basis_cutoff(latvec);
  else basis_cutoff=7.5;

  Jastrow2_wf_writer jast2writer;
  jast2writer.set_atoms(atoms);
  string jast2outname=outputname+".jast2";
  ofstream jast2out(jast2outname.c_str());
  print_std_jastrow2(jast2writer, jast2out, basis_cutoff);
  jast2out.close();

  // Jastrow 3 output
  string jast3outname=outputname+".jast3";
  ofstream jast3out(jast3outname.c_str());
  vector<string> unique_atoms;
  find_unique_atoms(atoms, unique_atoms);
  print_3b_jastrow2(jast3out,unique_atoms,basis_cutoff);
  jast3out.close();

  //------------------------------------------System output

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

  cout << "Writing real orbitals ...  " << endl; 
  for (int kpt=0; kpt<rkpts.size(); kpt++) {
    slwriter.orbtype="ORBITALS";

    vector <string> words; 
    string space = " "; 
    split(rkpts[kpt], space, words);
    stringstream lk;
    lk << atoi(words[0].c_str())-1;

    string lst = lk.str();
    string forb=outputname+"_"+lst+".orb"; 
    cout << "  " << forb + ": "; 
    ofstream orbout(forb.c_str());
    print_orbitals(orbout, centers, nbasis_list, moCoeff[kpt]);
    orbout.close();
    string slateroutname=outputname+"_"+lst+".slater";
    slwriter.orbname=forb; 
    ofstream slaterout(slateroutname.c_str());
    slwriter.print_wavefunction(slaterout);
    slaterout.close();

    string sysoutname=outputname+"_"+lst+".sys";
    ofstream sysout(sysoutname.c_str());
    sysout.precision(12);
    sysout << "SYSTEM { ";
    if(latvec.size() > 0) sysout << " PERIODIC \n";
    else sysout << " MOLECULE \n";
    sysout << "  NSPIN { " << slwriter.nup << "  "
      << slwriter.ndown << " } \n";


    if(latvec.size() > 0) { 
      sysout << "LATTICEVEC { \n";
      for(int i=0; i< 3; i++) {
        for(int j=0; j< 3; j++) {
          sysout << latvec[i][j] << "   ";
        }
        sysout << endl;
      }
      sysout << " } " << endl;
      sysout << " origin { " << origin[0] << "   " 
        << origin[1] << "   " << origin[2] << "  } " << endl;
      sysout << "  cutoff_divider " 
        << basis_cutoff*2.0/cutoff_length << endl;
      sysout << "  kpoint { " << rkptCoord[kpt][0] 
        << "   " << rkptCoord[kpt][1] 
        << "   " << rkptCoord[kpt][2] << " } " << endl;
    }

    for(int at=0; at <natoms; at++) {
      atoms[at].print_atom(sysout);
    }
    sysout << "}\n\n\n";



    int npsp=pseudo.size();
    for(int psp=0; psp < npsp; psp++) {
      pseudo[psp].print_pseudo(sysout);
    }
    sysout.close();

  }
  //  cout << "Now writing complex orbitals" << endl; 
  if(cmplx)  { 
    cout << "Writing complex orbitals ...  " << endl; 
    for (int kpt=0; kpt<ckpts.size(); kpt++) {
      slwriter.orbtype="CORBITALS";
      vector <string> words; 
      string space = " "; 
      split(ckpts[kpt], space, words);
      //    string lst=string(atoi(words[0].c_str())-1);
      stringstream lk;
      lk << atoi(words[0].c_str())-1;
      string lst = lk.str();
      string forb=outputname+"_"+lst+".orb"; 
      cout << "  " << forb << ": "; 
      ofstream orbout(forb.c_str());
      print_orbitals(orbout, centers, nbasis_list, CmoCoeff[kpt]);
      orbout.close(); 
      string slateroutname=outputname+"_"+lst+".slater";
      slwriter.orbname=forb; 
      ofstream slaterout(slateroutname.c_str());
      slwriter.print_wavefunction(slaterout);
      slaterout.close();

      string sysoutname=outputname+"_"+lst+".sys";
      ofstream sysout(sysoutname.c_str());
      sysout.precision(12);
      sysout << "SYSTEM { ";
      if(latvec.size() > 0) sysout << " PERIODIC \n";
      else sysout << " MOLECULE \n";
      sysout << "  NSPIN { " << slwriter.nup << "  "
        << slwriter.ndown << " } \n";


      if(latvec.size() > 0) { 
        sysout << "LATTICEVEC { \n";
        for(int i=0; i< 3; i++) {
          for(int j=0; j< 3; j++) {
            sysout << latvec[i][j] << "   ";
          }
          sysout << endl;
        }
        sysout << " } " << endl;
        sysout << " origin { " << origin[0] << "   " 
          << origin[1] << "   " << origin[2] << "  } " << endl;
        sysout << "  cutoff_divider " 
          << basis_cutoff*2.0/cutoff_length << endl;
        sysout << "  kpoint { " << ckptCoord[kpt][0] 
          << "   " << ckptCoord[kpt][1] 
          << "   " << ckptCoord[kpt][2] << " } " << endl;
      }

      for(int at=0; at <natoms; at++) {
        atoms[at].print_atom(sysout);
      }
      sysout << "}\n\n\n";



      int npsp=pseudo.size();
      for(int psp=0; psp < npsp; psp++) {
        pseudo[psp].print_pseudo(sysout);
      }
      sysout.close();
    }
  }


}

//######################################################################

void get_crystal_latvec(istream & infile,
    vector< vector < double > > & latvec) {
  string testword;
  while(infile >> testword) {
    //cheesy way to get the converter to work for a molecular case
    //I don't know if it would work in general.
    // if(testword == "MOLECULAR") {

    //  vector <double> tmp;
    //  for(int i=0; i< 3; i++) {
    //tmp.push_back(100);
    //  }
    //  for(int i=0; i< 3; i++) latvec.push_back(tmp);
    //  break;
    //}
    if(testword == "DIRECT") {
      infile >> testword;  //LATTICE
      infile >> testword;  //VECTORS
      infile >> testword;  //COMPON.
      infile >> testword;  //(A.U.)
      if(testword=="(A.U.)") {
        //cout << "found lattice parms" << endl;
        infile.ignore(150, '\n');
        infile >> testword;
        //cout << testword << endl;
        infile >> testword;
        //cout << testword << endl;
        infile.ignore(150,'\n');
        latvec.reserve(3);
        for(int i=0; i< 3; i++) {
          latvec[i].reserve(3);
          vector <double> tmp;
          for(int j=0; j< 3; j++) {
            double dummy;
            infile >> dummy;
            //here we try to avoid having issues with the Ewald summation
            if(dummy > 100.) { 
              cout << "WARNING: rescaling lattice vector to " << 100 << " from " << dummy 
                << ". This will probably not change your results, but if it does, change EWALD_GMAX in the system input to something higher." << endl;
              dummy=100.;
            }
            tmp.push_back(dummy);
          }
          latvec.push_back(tmp);
          //cout << endl;
          infile.ignore(150, '\n');
        }
        break;
      }
    }
  }

  if(latvec.size() == 0) {

    //cout << "Couldn't find lattice vector in output file." << endl;
    //exit(1);
  }
}

//######################################################################

string erasetail(string s) {
  char ms[] = "1"; 
  string ns=""; 
  int i=0; 
  while (s[i]!='1' and i < s.length()) {
    i++; 
  } 
  return s.substr(0, i); 
}


void get_crystal_atoms(istream & infile,
    vector < Atom > & atoms) {
  string line;
  string space=" ";
  Atom temp_atom;
  double bohr=0.5291772083;

  while(getline(infile, line)) {
    vector <string> words; 
    split(line, space, words);

    //2003 detection
    if( words.size() >= 5 && 
        words[1]=="ATOM" && words[2]=="X(ANGSTROM)") {



      //cout << "found atom section " << line << endl;
      infile.ignore(150, '\n'); //ignore line of *'s
      getline(infile, line);
      const char endmatch='*';
      while(true) {
        vector < string> atomwords;
        split(line, space, atomwords);
        if(atomwords.size() > 4) {
          //cout << "line " << line << endl;
          temp_atom.name=erasetail(atomwords[2]);
          temp_atom.pos[0]=atof(atomwords[3].c_str())/bohr;
          temp_atom.pos[1]=atof(atomwords[4].c_str())/bohr;
          temp_atom.pos[2]=atof(atomwords[5].c_str())/bohr;
          atoms.push_back(temp_atom);
          getline(infile, line);
        }
        else break;
      }
      break;
    }
    //crystal2003 molecules
    else if(words.size() >=4 && words[0]=="ATOM"  
        && words[1]=="X(ANGSTROM)") { 
      //cout << "found atom section " << line << endl;
      infile.ignore(150, '\n'); //ignore line of *'s
      getline(infile, line);
      const char endmatch='*';
      while(true) {
        vector < string> atomwords;
        split(line, space, atomwords);
        if(atomwords.size() > 4) {
          //cout << "line " << line << endl;
          temp_atom.name=atomwords[3];
          temp_atom.pos[0]=atof(atomwords[4].c_str())/bohr;
          temp_atom.pos[1]=atof(atomwords[5].c_str())/bohr;
          temp_atom.pos[2]=atof(atomwords[6].c_str())/bohr;
          atoms.push_back(temp_atom);
          getline(infile, line);
        }
        else break;
      }
      break; 
    }
    //crystal98 with COORPRT detection
    else if(words.size() >= 5 && 
        words[0]=="ATOM" && words[1]=="X" ) { //1998

      //cout << "found atom section " << line << endl;
      infile.ignore(150, '\n'); //ignore line of *'s
      getline(infile, line);
      const char endmatch='*';
      while(true) {
        vector < string> atomwords;
        split(line, space, atomwords);
        if(atomwords[0]!="INFORMATION") {
          //cout << "line " << line << endl;
          temp_atom.name=atomwords[1];
          temp_atom.pos[0]=atof(atomwords[2].c_str())/bohr;
          temp_atom.pos[1]=atof(atomwords[3].c_str())/bohr;
          temp_atom.pos[2]=atof(atomwords[4].c_str())/bohr;
          atoms.push_back(temp_atom);
          getline(infile, line);
        }
        else break;
      }
      break;
    }


  }

  //We truncate all atom names to two letters, since 
  //crystal does that wrt the basis set & psp

  for(vector<Atom>::iterator at=atoms.begin();
      at != atoms.end(); at++) {
    if(at->name.size() > 2) 
      at->name.erase(at->name.begin()+2, at->name.end());
  }


}

//##################################################################

void get_crystal_basis(istream & infile,
    vector <Gaussian_basis_set> & basis) {
  string line;
  string space=" ";
  vector <string> words;
  vector < vector < string > > basis_sections;
  vector <string> basis_labels;
  vector <string> blank_strvec;
  while(getline(infile, line)) {
    words.clear();
    split(line, space, words);
    if(words.size() > 4 && words[0]=="ATOM" && words[1] == "X(AU)" && 
        (words[4]=="NO." or words[4]=="N.")) {
      infile.ignore(150,'\n');

      cout << "found basis " << line << endl;
      getline(infile, line);
      cout << "first line " << line << endl;

      const char endmatch='*';
      string currname;
      currname.resize(2);
      string currtype;
      currtype.resize(2);
      int nbasis=-1;
      while(search_n(line.begin(), line.end(), 5, endmatch) == line.end()
          && line.size() != 0) {

        if(line[3] != ' ') {
          currname[0]=line[5];
          currname[1]=line[6];
          getline(infile, line);
          vector <string> this_sec;
          while(line.size() >= 4 && line[3] == ' ') {
            this_sec.push_back(line);              
            getline(infile, line);
          }
          if(this_sec.size() > 0) {
            basis_sections.push_back(this_sec);
            basis_labels.push_back(currname);
          }
        }
      }
      break;
    }
  }

  assert(basis_sections.size() == basis_labels.size());

  //for(unsigned int i=0; i< basis_sections.size(); i++) {
  //  cout << "label " << basis_labels[i] << endl;
  //  for(unsigned int j=0; j< basis_sections[i].size(); j++) {
  //    cout << basis_sections[i][j] << endl;
  //  }
  //}

  basis.resize(basis_sections.size());
  for(unsigned int bas=0; bas < basis_sections.size(); bas++) {
    basis[bas].label=basis_labels[bas];
    vector < string > indiv_types;
    vector < vector < string> > indiv_funcs;

    int fnum=-1;
    string currtype;
    currtype.resize(2);

    //cout << "parsing the basis" << endl;

    for(unsigned int j=0; j< basis_sections[bas].size(); ) {
      if(basis_sections[bas][j][36] != ' ') {
        currtype[0]=basis_sections[bas][j][36];
        currtype[1]=basis_sections[bas][j][37];
        indiv_types.push_back(currtype);
        j++;
        vector <string> thisfunc;
        while(j < basis_sections[bas].size() && basis_sections[bas][j][36] == ' ') {
          thisfunc.push_back(basis_sections[bas][j]);
          j++;
        }
        indiv_funcs.push_back(thisfunc);
      }
    }

    assert(indiv_funcs.size()==indiv_types.size());

    //for(unsigned int i=0; i< indiv_funcs.size(); i++) {
    //  cout << "type " << indiv_types[i] << endl;
    //  for(unsigned int j=0; j< indiv_funcs[i].size(); j++) {
    //    cout << "line " << indiv_funcs[i][j] << endl;
    //  }
    //}

    int nfuncs=0;
    for(vector<string>::iterator btype=indiv_types.begin(); 
        btype != indiv_types.end(); btype++) {
      if((*btype) == "SP") nfuncs+=2;
      else nfuncs++;
    }
    //cout << nfuncs << " total functions " << endl;

    basis[bas].exponents.resize(nfuncs);
    basis[bas].coefficients.resize(nfuncs);
    int currf=0;
    for(unsigned int f=0; f< indiv_funcs.size(); f++) {
      if(indiv_types[f]=="SP") {
        basis[bas].types.push_back("S ");
        basis[bas].types.push_back("P ");
      }
      else if(indiv_types[f]=="D ")
        basis[bas].types.push_back("5D");
      else if(indiv_types[f]=="F ")
        basis[bas].types.push_back("7F_crystal"); 
      else basis[bas].types.push_back(indiv_types[f]);
      //vector <double> tmpexp;
      //vector < double> tmpcoeff;
      //vector < double> tmpcoeff2;//For SP functions
      for(vector <string>::iterator line=indiv_funcs[f].begin(); 
          line != indiv_funcs[f].end(); line++) {
        string exptxt;
        string::iterator lineb=line->begin();
        exptxt.assign(lineb+40, lineb+50);
        double exp=atof(exptxt.c_str());
        basis[bas].exponents[currf].push_back(exp);
        //cout << indiv_types[f] << ": "<< exp; 
        string coefftxt, coefftxt2;
        if(indiv_types[f]=="S ") {
          coefftxt.assign(lineb+50, lineb+60);
        }
        else if(indiv_types[f]=="P ") {
          coefftxt.assign(lineb+60, lineb+70);
        }
        else if(indiv_types[f]=="D ") {
          coefftxt.assign(lineb+70, lineb+80);
        }
        else if(indiv_types[f]=="SP") {
          coefftxt.assign(lineb+50, lineb+60);
          coefftxt2.assign(lineb+60, lineb+70);
        }
        else if(indiv_types[f]=="F ") {//the position of D/F/G
          coefftxt.assign(lineb+70, lineb+80); 
        }
        else {
          cout << "WARNING!!  Don't know type " << indiv_types[f] << endl;
        }
        basis[bas].coefficients[currf].push_back(atof(coefftxt.c_str()));
        //	cout << coefftxt.c_str() << endl; 
        //cout << "index: " << currf << endl; 
        //cout << "ORB:" << indiv_types[f] << " "<< atof(coefftxt.c_str())<< endl; 
        if(indiv_types[f]=="SP") {
          basis[bas].exponents[currf+1].push_back(exp);
          basis[bas].coefficients[currf+1].push_back(atof(coefftxt2.c_str()));        
        }
      }

      currf++;
      if(indiv_types[f]=="SP") currf++;

    }
  }      

  int nbasis=basis.size();
  //Strip whitespace from basis names and set options
  for(int bas=0; bas < nbasis; bas++) {
    if(basis[bas].label[1]==' ') {
      basis[bas].label.erase(basis[bas].label.end()-1, basis[bas].label.end());
    }
    //This isn't necessary in Crystal2009 it seems
    //basis[bas].options=" NORMTYPE CRYSTAL \n NORENORMALIZE \n";

    //Also strip whitespace from the basis types
    for(vector <string>::iterator i=basis[bas].types.begin(); 
        i!= basis[bas].types.end(); i++) {
      if( (*i)[1] == ' ') {
        i->erase(i->begin()+1, i->end());
      }
    }

  }

}


//##################################################################
void get_crystal_pseudo(istream & infile,
    vector <Gaussian_pseudo_writer> & pseudo) {
  string testword;
  while (infile >> testword) {
    //------------------------------
    // Pseudopotentials
    if(testword == "PSEUDOPOTENTIAL") {
      infile >> testword;
      if(testword == "INFORMATION") {
        infile.ignore(125, '\n'); //rest of line
        infile.ignore(125, '\n'); //Line of **'s
        infile.ignore(125, '\n'); //empty line
        string line;
        getline(infile, line);
        const char endmatch='*';
        vector <double> double_blank;
        vector <int> n_blank;
        Gaussian_pseudo_writer pseudo_blank;
        string temp;
        int currpsp=-1;
        int currl=-1;
        //search loop
        vector <string> words;
        string space=" ";
        //cout << "searching " << line << endl;
        while(search_n(line.begin(), line.end(), 20, endmatch) == line.end()) { 
          words.clear();
          split(line,space,words);

          //cout << "pre " << line << endl;
          if(words[0]=="INFORMATION" or words[0]=="NUCLEAR") break;
          if(line.size() > 2 && line[1]=='A' && line[2]=='T') {
            //new atom
            pseudo.push_back(pseudo_blank);
            currpsp++;
            temp.assign(line.begin()+15, line.begin()+18);
            pseudo[currpsp].atomnum=atoi(temp.c_str())%100;
            pseudo[currpsp].label=element_lookup_caps[pseudo[currpsp].atomnum];
            temp.assign(line.begin()+34, line.begin()+41);
            //pseudo[currpsp].effcharge=int(atof(words[5].c_str()));
            //cout << "found " << pseudo[currpsp].label << endl;
            cout << "Found " << pseudo[currpsp].label << " : effective charge " << atoi(temp.c_str()) << endl; 
            pseudo[currpsp].effcharge=int(atoi(temp.c_str())); 
            currl=-1;
          }
          else if(words.size() > 1 && words[0]!="TYPE") {
            //cout << "LINE: "<<  line << " "<< line.size() << line[5] << endl;
            if(words[0][0]=='W' || words[0][0] == 'P') {
              pseudo[currpsp].exponents.push_back(double_blank);
              pseudo[currpsp].coefficients.push_back(double_blank);
              pseudo[currpsp].nvalue.push_back(n_blank);
              currl++;
              words.erase(words.begin());
              words.erase(words.begin());
            }
            if(words.size() > 2) {
              
              string exps, coeffs, ns; 
              exps.assign(line.begin()+8, line.begin()+22);
              coeffs.assign(line.begin()+22, line.begin()+35); 
              ns.assign(line.begin()+35, line.begin()+39); 
              if (coeffs.find("*", 0) == 0) {
                cout << "****WARNING: Fail to read coefficient (****** occurs) for gaussian function with exponent: "<< exps << endl; 
              }

              //              pseudo[currpsp].exponents[currl].push_back(atof(words[0].c_str()));
              pseudo[currpsp].exponents[currl].push_back(atof(exps.c_str())); 
              //              pseudo[currpsp].coefficients[currl].push_back(atof(words[1].c_str()));
              pseudo[currpsp].coefficients[currl].push_back(atof(coeffs.c_str()));
              //              pseudo[currpsp].nvalue[currl].push_back(atoi(words[2].c_str()));
              pseudo[currpsp].nvalue[currl].push_back(atoi(ns.c_str()));
              //cout << "adding " << exps << " " << coeffs << " " << ns << " done " <<  endl;
              if(words.size() > 3) { 
                exps.assign(line.begin()+39, line.begin()+53);
                coeffs.assign(line.begin()+53, line.begin()+66); 
                ns.assign(line.begin()+66, line.begin()+70); 
                if (coeffs.find("*", 0) == 0) {
                  cout << "****WARNING: Fail to read coefficient (****** occurs) for gaussian function with exponent: "<< exps << endl; 
                }
                //                pseudo[currpsp].exponents[currl].push_back(atof(words[3].c_str()));
                //                pseudo[currpsp].coefficients[currl].push_back(atof(words[4].c_str()));
                //                pseudo[currpsp].nvalue[currl].push_back(atoi(words[5].c_str()));
                pseudo[currpsp].exponents[currl].push_back(atof(exps.c_str())); 
                pseudo[currpsp].coefficients[currl].push_back(atof(coeffs.c_str())); 
                pseudo[currpsp].nvalue[currl].push_back(atof(ns.c_str())); 
                //cout << "adding second " << exps << " " << coeffs << " " << ns << endl;
                
              }
            }
            else { 
              cout << "WARNING: Pseudopotential output seems corrupted!" << endl;
            }

          }
          getline(infile, line);
        }
        //----end search loop
      }
      break;
    }
  }

  // cout << "done " << endl;
  int npseud=pseudo.size();
  for(int ps=0; ps < npseud; ps++) {
    vector <double> tmp=pseudo[ps].exponents[0];
    pseudo[ps].exponents.erase(pseudo[ps].exponents.begin());
    pseudo[ps].exponents.push_back(tmp);


    tmp=pseudo[ps].coefficients[0];
    pseudo[ps].coefficients.erase(pseudo[ps].coefficients.begin());
    pseudo[ps].coefficients.push_back(tmp);
    vector<int> tmp2=pseudo[ps].nvalue[0];
    pseudo[ps].nvalue.erase(pseudo[ps].nvalue.begin());
    pseudo[ps].nvalue.push_back(tmp2);

    // pseudo[ps].exponents.push_back(pseudo[ps].exponents[0]);
    // pseudo[ps].nvalue.push_back(pseudo[ps].nvalue[0]);
    // pseudo[ps].coefficients.push_back(pseudo[ps].coefficients[0]);
    // pseudo[ps].exponents.erase(pseudo[ps].exponents.begin());
    // pseudo[ps].coefficients.erase(pseudo[ps].coefficients.begin());
    // pseudo[ps].nvalue.erase(pseudo[ps].nvalue.begin());
    //cout << "size " << pseudo[ps].exponents.size() << endl;
    //pseudo[ps].print_pseudo(cout);
  }

}

//######################################################################
/*!
  Given a starting position and an input stream, append the MO
  coefficients from Crystal to the vector of vectors given.
  moCoeff will be in form [mo][coeff].
  */
int readMO(istream & is, long int start, vector < vector <double> > & moCoeff) {
  is.clear();
  string line;
  is.seekg(start);
  //  cout << "position " << is.tellg() << endl;
  int totmo=moCoeff.size();
  string space=" ";
  vector <double> emptyVector;
  while(getline(is, line) ) {
    if(line.size() > 15 && line[5]=='(' && line[15]==')') break;
    if(line.size() > 45 && line[35]=='(' && line[45]==')') break;

    //check for various words that crystal can end with
    if(line[1]=='E' || line[3] == 'T' || line[6]=='B') break;
    vector <string> words;
    //cout << "start line " << line << endl;
    split(line, space, words);
    int nmo_this=words.size();
    if (words[1]=="NEWK") break; 
    for(int i=0; i< nmo_this; i++)
      moCoeff.push_back(emptyVector);
    //cout << "nmo_this " << nmo_this << endl;
    //is.ignore(150, '\n'); //empty line
    getline(is, line);  //sometimes crystal outputs a newline, and sometimes not
    if(line.size() < 15) getline(is, line);
    while(line.size() > 15 && line[1] !='E') {
      words.clear();
      split(line, space, words);
      //cout << "line " << line << " size " << words.size() << endl;
      if(words[0]=="BETA") break;
      if(words[1]=="NEWK") break; 
      if(words.size()!= nmo_this+1) {
        cerr << "Problem reading in MOs: expected " << nmo_this << " orbitals, got " << words.size()-1 << endl;
        cerr << "Line:" << line << endl;
        exit(1);
      }
      assert(words.size() == nmo_this+1);
      for(int i=0; i< nmo_this; i++) {
        moCoeff[totmo+i].push_back(atof(words[i+1].c_str()));
      }
      getline(is, line);
    }
    totmo+=nmo_this;
  }
  return totmo;
}

// version for complex orbitals
int readMO(istream & is, long int start,
    vector < vector <dcomplex> > & moCoeff) {
  is.clear();
  string line;
  is.seekg(start);
  //  cout << "position " << is.tellg() << endl;
  int totmo=moCoeff.size();
  string space=" ";
  vector <dcomplex> emptyVector;
  while(getline(is, line) ) {
    if(line.size() > 15 && line[5]=='(' && line[15]==')') break;
    if(line.size() > 45 && line[35]=='(' && line[45]==')') break;
    //check for various words that crystal can end with
    if(line[1]=='E' || line[3] == 'T' || line[6]=='B') break;
    vector <string> words;
    //cout << "start line " << line << endl;
    split(line, space, words);
    if (words[1]=="NEWK") break; 
    int nmo_this=words.size()/2;
    for(int i=0; i< nmo_this; i++)
      moCoeff.push_back(emptyVector);
    //cout << "nmo_this " << nmo_this << endl;
    //is.ignore(150, '\n'); //empty line
    getline(is, line);  //sometimes crystal outputs a newline, and sometimes not
    if(line.size() < 15) getline(is, line);
    while(line.size() > 15 && line[1] !='E') {
      words.clear();
      split(line, space, words);
      //cout << "line " << line << " size " << words.size() << endl;
      if(words[0]=="BETA") break;
      if (words[1]=="NEWK") break; 

      if(words.size()!= 2*nmo_this+1) {
        cerr << "Problem reading in MOs: expected " << nmo_this << " orbitals, got " << words.size()-1 << endl;
        cerr << "Line:" << line << endl;
        exit(1);
      }
      
      for(int i=0; i< nmo_this; i++) {
        moCoeff[totmo+i].push_back(
            dcomplex( atof(words[2*i+1].c_str()), atof(words[2*i+2].c_str()) )
            );
        //cout << "imag. part " << atof(words[2*i+2].c_str()) << endl;
      }
      getline(is, line);
    }
    totmo+=nmo_this;
  }
  return totmo;
}

// ---------------------------------------------------------------------------

void fort10input(istream & is,
    string & fort10file,
    vector <Atom> & atoms,
    Slat_wf_writer & slwriter,
    vector <Gaussian_basis_set> & basis,
    vector <double> & origin,
    vector < vector <double> > & latvec,
    vector < vector <double> > & moCoeff,
    vector <double> & shift) {

  string calctype=slwriter.calctype;

  //reading from the output of the readcrys10.f program to get better CO's.
  if(fort10file != "") {
    cout << "Reading from the fort10 formatted file " << endl;
    ifstream fort10(fort10file.c_str());
    if(!fort10) {
      cout << "Couldn't open " << fort10file << endl;
      exit(1);
    }

    int nfunctions=moCoeff[0].size();

    int nspin=1;
    if(calctype=="UHF") {
      cout << "fort10file doesn't work for UHF yet.." << endl;
      exit(1);
    }

    string dummy;
    int monum;
    int nmatch=0;
    cout << "Matching the .o orbitals with the .10 orbitals..this may "
      << "take some time" << endl;
    for(int f=0; f< nfunctions; f++) { //total number of MO's
      vector <double> currmo;
      currmo.reserve(nfunctions);
      fort10.ignore(120, '\n'); //line of =====
      fort10.ignore(120, '\n'); //
      fort10 >> dummy >> monum;
      //cout << "dummy " << dummy << endl;
      if(monum != f+1) {
        cout << "monum in fort10file doesn't match what it should be: "
          << "monum " << monum << " calculated function " << f+1 << endl;
        exit(1);
      }
      int tempfunc;
      double tempval;
      for(int f2=0; f2 < nfunctions; f2++) {
        fort10 >> tempfunc >> tempval;
        if(tempfunc != f2+1) {
          cout << "functions don't match in fort10file " << tempfunc
            << " calculated " << f2+1 << endl;
          exit(1);
        }
        currmo.push_back(tempval);
        //now loop through the mo's we've found and find the equivalent.
      }
      //cout << "currmo size " << currmo.size() << endl;

      for(vector < vector <double > >:: iterator mo=moCoeff.begin();
          mo != moCoeff.end(); mo++) {
        double dot=0, mag_out=0, mag_fort10=0;
        vector <double>::iterator curr=currmo.begin();
        for(vector<double>::iterator moc=mo->begin(); 
            moc != mo->end(); moc++) {
          assert(curr != currmo.end());
          dot+=(*moc)*(*curr);
          mag_out+=(*moc)*(*moc);
          mag_fort10+=(*curr)*(*curr);
          curr++;
        }
        dot /= sqrt(mag_out*mag_fort10);
        if(fabs(dot-1) < .01) {
          //cout << "match: " << dot << endl;
          *mo= currmo;
          nmatch++;
        }
      }
      /*
         for(int mo=0; mo < totmo; mo++) {
         double dot=0;
         double mag_out=0, mag_fort10=0;
         for(int f=0; f< nfunctions; f++) {
         dot+=moCoeff[mo][f]*currmo[f];
         mag_out+=moCoeff[mo][f]*moCoeff[mo][f];
         mag_fort10+=currmo[f]*currmo[f];
         }
         dot /= sqrt(mag_out*mag_fort10);

      //If they're the same molecular orbital, replace the output file
      //one with the read in one.
      if(fabs(dot-1) < .01) {
      cout << monum << " -> " << mo+1 << " dot " << dot <<  endl;
      for(int f=0; f< nfunctions; f++) {
      moCoeff[mo][f]=currmo[f];
      }
      }

      }
      */
    }
    cout << "Matched " << nmatch << " orbitals " << endl;
  }
}

//-----------------------------------------------------------------------------

void MO_analysis(istream & is,
    string & fort10file,
    vector <Atom> & atoms,
    Slat_wf_writer & slwriter,
    vector <Gaussian_basis_set> & basis,
    vector <double> & origin,
    vector < vector <double> > & latvec,
    vector < vector <double> > & moCoeff,
    vector <double> & shift, 
    int totmo, string mo_filename="mo_analysis") {

  int natoms=atoms.size();
  ofstream an_out(mo_filename.c_str());

  const double print_thresh=1e-3;
  //const double pi=3.1415926535897932385;
  const double pi=3.1415926535;
  double snorm=1./sqrt(4.*pi);
  double pnorm=sqrt(3.)*snorm; // sqrt(3/4/pi)
  vector <double> dnorm;
  dnorm.push_back(.5*sqrt(5./(4*pi)));//sqrt(5, 16)
  dnorm.push_back(sqrt(15./(4*pi)));
  dnorm.push_back(sqrt(15./(4*pi)));
  dnorm.push_back(.5*sqrt(15./(4.*pi)));
  dnorm.push_back(sqrt(15./(4*pi)));
  vector <double> fnorm; 
  //f orbital normalizations are from <http://winter.group.shef.ac.uk/orbitron/AOs/4f/equations.html>
  fnorm.push_back( sqrt( 7./(16.*pi)) ); 
  fnorm.push_back( sqrt(21./(32.*pi)) ); 
  fnorm.push_back( sqrt(21./(32.*pi)) ); 
  fnorm.push_back( sqrt(105./(16.*pi)) ); 
  fnorm.push_back( sqrt(105./(4.*pi))  ); //xyz 
  fnorm.push_back( sqrt(35./(32.*pi))  ); 
  fnorm.push_back( sqrt(35./(32.*pi))  ); 

  vector <string> dnames(5);
  dnames[0]="z2r2";
  dnames[1]="xz  ";
  dnames[2]="yz  ";
  dnames[3]="x2y2";
  dnames[4]="xy  ";
  vector <string> pnames(3);
  pnames[0]="x   ";
  pnames[1]="y   ";
  pnames[2]="z   ";

  vector <string> fnames(7);
  fnames[0]="F0   ";
  fnames[1]="Fp1  ";
  fnames[2]="Fm1  ";
  fnames[3]="Fp2   ";
  fnames[4]="Fxyz  ";
  fnames[5]="Fp3   ";
  fnames[6]="Fm3   ";
  for(int mo=0; mo < totmo; mo++) {
    int func=0;
    an_out << "\n----------------\n";
    an_out << "MO " << mo << endl;
    for(int at=0; at < natoms; at++) {
      int bas=atoms[at].basis;
      int nbasis=basis[bas].types.size();
      for(int i=0; i< nbasis; i++) {
        if(basis[bas].types[i] == "S") {
          moCoeff[mo][func]*=snorm;
          if(fabs(moCoeff[mo][func]) > print_thresh) {
            an_out << atoms[at].name<< at  << "  S     " << moCoeff[mo][func]
              << endl;
          }
          func++;
        }
        else if(basis[bas].types[i] == "P") {
          for(int j=0; j< 3; j++) {
            moCoeff[mo][func]*=pnorm;
            if(fabs(moCoeff[mo][func]) > print_thresh) {
              an_out << atoms[at].name << at << "  "  << "P" 
                << pnames[j] << " " << moCoeff[mo][func]
                << endl;
            }
            func++;
          }
        }
        else if(basis[bas].types[i] == "5D") {
          for(int j=0; j< 5; j++) {
            moCoeff[mo][func]*=dnorm[j];
            if(fabs(moCoeff[mo][func]) > print_thresh) {
              an_out << atoms[at].name << at << "  "   << "D" 
                << dnames[j] << " " <<  moCoeff[mo][func]
                << endl;
            }
            func++;
          }
        }
        else if(basis[bas].types[i] == "7F_crystal") {
          for(int j=0; j<7; j++) {
            moCoeff[mo][func]*=fnorm[j];
            if(fabs(moCoeff[mo][func]) > print_thresh) {
              an_out << atoms[at].name << at << "  "   << "F" 
                << fnames[j] << " " <<  moCoeff[mo][func]
                << endl;
            }
            func++;
          }
        }
        else {
          cout << "Error: unknown basis type in read_crystal_orbital"
            << endl;
          exit(1);
        }
      }
    }
  }

  an_out.close();

}

void MO_analysis(istream & is,
    vector <Atom> & atoms,
    Slat_wf_writer & slwriter,
    vector <Gaussian_basis_set> & basis,
    vector <double> & origin,
    vector < vector <double> > & latvec,
    vector < vector <dcomplex> > & moCoeff,
    vector <double> & shift, 
    int totmo, string mo_filename="mo_analysis") {

  int natoms=atoms.size();
  ofstream an_out(mo_filename.c_str());

  const double print_thresh=1e-3;
  //const double pi=3.1415926535897932385;
  const double pi=3.1415926535;
  double snorm=1./sqrt(4.*pi);
  double pnorm=sqrt(3.)*snorm;
  vector <double> dnorm;
  dnorm.push_back(.5*sqrt(5./(4*pi)));
  dnorm.push_back(sqrt(15./(4*pi)));
  dnorm.push_back(sqrt(15./(4*pi)));
  dnorm.push_back(.5*sqrt(15./(4.*pi)));
  dnorm.push_back(sqrt(15./(4*pi)));
  vector <double> fnorm; 
  //f orbital normalizations are from <http://winter.group.shef.ac.uk/orbitron/AOs/4f/equations.html>
  fnorm.push_back( sqrt( 7./(16.*pi)) ); 
  fnorm.push_back( sqrt(21./(32.*pi)) ); 
  fnorm.push_back( sqrt(21./(32.*pi)) ); 
  fnorm.push_back( sqrt(105./(16.*pi)) ); 
  fnorm.push_back( sqrt(105./(4.*pi))  ); //xyz 
  fnorm.push_back( sqrt(35./(32.*pi))  ); 
  fnorm.push_back( sqrt(35./(32.*pi))  ); 
  vector <string> fnames(7);
  fnames[0]="F0   ";
  fnames[1]="Fp1  ";
  fnames[2]="Fm1  ";
  fnames[3]="Fp2   ";
  fnames[4]="Fxyz  ";
  fnames[5]="Fp3   ";
  fnames[6]="Fm3   ";

  vector <string> dnames(5);
  dnames[0]="z2r2";
  dnames[1]="xz  ";
  dnames[2]="yz  ";
  dnames[3]="x2y2";
  dnames[4]="xy  ";
  vector <string> pnames(3);
  pnames[0]="x   ";
  pnames[1]="y   ";
  pnames[2]="z   ";
  for(int mo=0; mo < totmo; mo++) {
    int func=0;
    an_out << "\n----------------\n";
    an_out << "MO " << mo << endl;
    for(int at=0; at < natoms; at++) {
      int bas=atoms[at].basis;
      int nbasis=basis[bas].types.size();
      for(int i=0; i< nbasis; i++) {
        if(basis[bas].types[i] == "S") {
          moCoeff[mo][func]*=snorm;
          if(abs(moCoeff[mo][func]) > print_thresh) {
            an_out << atoms[at].name<< at  << "  S     " << moCoeff[mo][func]
              << endl;
          }
          func++;
        }
        else if(basis[bas].types[i] == "P") {
          for(int j=0; j< 3; j++) {
            moCoeff[mo][func]*=pnorm;
            if(abs(moCoeff[mo][func]) > print_thresh) {
              an_out << atoms[at].name << at << "  "  << "P" 
                << pnames[j] << " " << moCoeff[mo][func]
                << endl;
            }
            func++;
          }
        }
        else if(basis[bas].types[i] == "5D") {
          for(int j=0; j< 5; j++) {
            moCoeff[mo][func]*=dnorm[j];
            if(abs(moCoeff[mo][func]) > print_thresh) {
              an_out << atoms[at].name << at << "  "   << "D" 
                << dnames[j] << " " <<  moCoeff[mo][func]
                << endl;
            }
            func++;
          }
        }
        else if(basis[bas].types[i] == "7F_crystal") {
          for(int j=0; j< 7; j++) {
            moCoeff[mo][func]*=fnorm[j];
            if(abs(moCoeff[mo][func]) > print_thresh) {
              an_out << atoms[at].name << at << "  "   << "F" 
                << fnames[j] << " " <<  moCoeff[mo][func]
                << endl;
            }
            func++;
          }
        }
        else {
          cout << "Error: unknown basis type in read_crystal_orbital"
            << endl;
          exit(1);
        }
      }
    }
  }

  an_out.close();

}


// ----------------------------------------------------------------------------

//Reads in all the molecular orbitals that were outputed by Crystal
void read_crystal_orbital(istream & is,
    string & fort10file,
    vector <Atom> & atoms,
    Slat_wf_writer & slwriter,
    vector <Gaussian_basis_set> & basis,
    vector <double> & origin,
    vector < vector <double> > & latvec,
    vector < vector <double> > & moCoeff,
    vector <double> & shift) {

  assert(atoms.size() > 0);
  assert(basis.size() > 0);
  string dummy;
  string calctype=slwriter.calctype;
  //vector <int> nbasis;
  //vector <double > emptyVector; //to push onto moCoeff to add a mo.

  //Find out how many functions should be there.
  int natoms=atoms.size();
  int totfunctions=0;
  for(int at=0; at < natoms; at++ ) {
    int bas=atoms[at].basis;
    int nfunc=basis[bas].nfunc();
    totfunctions+=nfunc;
  }
  //cout << "Should be " << totfunctions << " functions " << endl;

  string space=" ";
  int totmo=0;
  vector <long int> eigen_start;
  vector <string> kpoints;
  while(is >> dummy) {
    //cout << "dummy " << dummy << endl;
    if(dummy == "FINAL") {
      is >> dummy;
      if(dummy == "EIGENVECTORS" ) {
        is.ignore(125, '\n'); //clear the line with FINAL EIG..
        string line;
        getline(is, line);
        // cout << "line " << line << endl;
        while(getline( is,line)) {
          if(line.size() > 15 && line[5]=='(' && line[15]==')') {
            //cout << line[5] << "  " << line[15] << endl;
            is.ignore(150, '\n');//two blank lines
            is.ignore(150, '\n');
            long int pos=is.tellg();
            string line2;
            getline(is, line2);
            vector <string> words;
            split(line2, space, words);
            if(words[0]!=words[1]) {
              //cout << "real k-point line " << line << endl;
              kpoints.push_back(line);
              eigen_start.push_back(pos);
            }
          }
        }
      }
    }
  }


  int nkpts=eigen_start.size();
  if(calctype=="UHF") {
    assert(nkpts%2==0);
    nkpts/=2;
  }
  cout << "Found " << nkpts << " k-points with real eigenvectors " << endl;
  for(int i=0; i< nkpts; i++) {
    cout << i << " : " << kpoints[i] <<  "   position " << eigen_start[i] 
      << endl;
  }
  int kpt;

  if(nkpts > 1) {
    cout << "Please choose a point[0-" << nkpts-1 << "]:";
    cout.flush();

    cin >> kpt;
    while(kpt < 0 || kpt >= nkpts) {
      cout << "\nout of range..please re-enter:";
      cout.flush();
      if(!(cin >> kpt)) { 
        cout << "\nError reading from STDIN" << endl; 
        exit(1);
      }
    }
    cout << endl;
  }
  else kpt=0;



  is.clear();
  totmo=readMO(is, eigen_start[kpt], moCoeff);
  if(calctype=="UHF") {
    slwriter.spin_dwn_start=moCoeff.size();
    totmo=readMO(is, eigen_start[kpt+nkpts], moCoeff);
  }
  // cout << "nmo's " << moCoeff.size() << endl;

  slwriter.kpoint.resize(3);
  vector<string> kwords;
  kpoints[kpt].erase(kpoints[kpt].find(')'));
  split(kpoints[kpt], space, kwords);
  assert(kwords.size()>=5);
  double max=0;
  for(int i=0; i< 3; i++) {
    slwriter.kpoint[i]=atoi(kwords[i+2].c_str());
    if(slwriter.kpoint[i] > max) max=slwriter.kpoint[i];
  }
  for(int i=0; i< 3; i++) 
    if(abs(max) > 1e-5) slwriter.kpoint[i]/=max;


  cout << "chosen k-point " << slwriter.kpoint[0] <<"   " 
    << slwriter.kpoint[1] << "   " << slwriter.kpoint[2] << endl;


  //cout << "Found " << moCoeff.size() << " MO's. and "
  //<< moCoeff[0].size() << " functions " << endl;
  if(totfunctions != (int) moCoeff[0].size()) {
    cout << "The number of basis functions doesn't match between what was"
      << "read(" << moCoeff[0].size() << ") and calculated("
      << totfunctions << ")" << endl;
    exit(1);
  }
  if(totmo!=moCoeff.size()) {
    cout << "in make_orb, totmo is " << totmo
      << " and moCoeff.size() is " << moCoeff.size()
      << ".  They should be equal." << endl;
  }


  fort10input(is, fort10file, atoms, slwriter, basis, origin,
      latvec, moCoeff, shift);

  // analysis of band character and NORMALIZATION(!) of coefficients
  string mo_filename="xxmo_analysis";
  for(int d=0; d< 3; d++) append_number(mo_filename,slwriter.kpoint[d]);
  MO_analysis(is, fort10file, atoms, slwriter, basis, origin,
      latvec, moCoeff, shift, totmo,mo_filename);


  //if our k-point isn't zero, we need to fix any shifts
  int f=0;
  Shifter shiftobj;
  shiftobj.origin=origin;

  //------------------------------------
  //Shift the atoms away from the edges.  Should
  //reduce the number of centers outside the cell most of the time.

  //vector <double> shift;
  //double shift_amount=0;
  //shift.push_back(shift_amount); shift.push_back(shift_amount); 
  //shift.push_back(shift_amount);
  if(latvec.size() > 0) { 
    //Now enforce the pbc's.. 

    for(unsigned int at=0; at< atoms.size(); at++) {
      vector <int> nshift;
      int bas=atoms[at].basis;
      int nfunc=basis[bas].nfunc();
      if(shiftobj.enforcepbc(atoms[at].pos, latvec, nshift)) {
        //cout << "at " << at << "  shifted " << nshift[0] << "  " << nshift[1] 
        //    << "   " << nshift[2] << endl;
        double kdots=0;
        for(int d=0; d< 3; d++) 
          kdots+=slwriter.kpoint[d]*nshift[d];
        kdots=cos(pi*kdots);
        //cout << "kdots " << kdots << endl;
        for(int i=0; i< nfunc; i++) {
          for(int mo=0; mo < totmo; mo++) {
            moCoeff[mo][f]*=kdots;
          }
          f++;
        }
      }
      else f+=nfunc;
    }
  }


}


//Reads in all the molecular orbitals that were outputed by Crystal
//version for complex coefficients
void read_crystal_orbital(istream & is,
    vector <Atom> & atoms,
    Slat_wf_writer & slwriter,
    vector <Gaussian_basis_set> & basis,
    vector <double> & origin,
    vector < vector <double> > & latvec,
    vector < vector <dcomplex> > & moCoeff,
    vector <double> & shift) {

  assert(atoms.size() > 0);
  assert(basis.size() > 0);
  string dummy;
  string calctype=slwriter.calctype;
  //vector <int> nbasis;

  //Find out how many functions should be there.
  int natoms=atoms.size();
  int totfunctions=0;
  for(int at=0; at < natoms; at++ ) {
    int bas=atoms[at].basis;
    int nfunc=basis[bas].nfunc();
    totfunctions+=nfunc;
  }
  //cout << "Should be " << totfunctions << " functions " << endl;

  string space=" ";
  int totmo=0;
  vector <long int> eigen_start;
  vector <string> kpoints;
  while(is >> dummy) {
    //cout << "dummy " << dummy << endl;
    if(dummy == "FINAL") {
      is >> dummy;
      if(dummy == "EIGENVECTORS" ) {
        is.ignore(125, '\n'); //clear the line with FINAL EIG..
        string line;
        getline(is, line);
        // cout << "line " << line << endl;
        while(getline( is,line)) {
          if(line.size() > 15 && line[5]=='(' && line[15]==')') {
            //cout << line[5] << "  " << line[15] << endl;
            is.ignore(150, '\n');//two blank lines
            is.ignore(150, '\n');
            long int pos=is.tellg();
            string line2;
            getline(is, line2);
            vector <string> words;
            split(line2, space, words);
            if(words[0]==words[1]) {
              //cout << "complex  k-point line " << line << endl;
              kpoints.push_back(line);
              eigen_start.push_back(pos);
            }
          }
        }
      }
    }
  }


  int nkpts=eigen_start.size();
  if(calctype=="UHF") {
    assert(nkpts%2==0);
    nkpts/=2;
  }
  cout << "Found " << nkpts << " k-points with complex eigenvectors " << endl;
  for(int i=0; i< nkpts; i++) {
    cout << i << " : " << kpoints[i] <<  "   position " << eigen_start[i] 
      << endl;
  }
  int kpt;
  if(nkpts > 1) {
    cout << "Please choose a point[0-" << nkpts-1 << "]:";
    cout.flush();

    cin >> kpt;
    while(kpt < 0 || kpt >= nkpts) {
      cout << "\nout of range..please re-enter:";
      cout.flush();
      if(!(cin >> kpt)) { 
        cout << "\nError reading from STDIN" << endl; 
        exit(1);
      }
    }
    cout << endl;
  }
  else kpt=0;

  is.clear();
  is.seekg(1);
  string line;
  int shrink_fact[3];
  while ( getline(is,line) ) {
    vector<string> words;
    split(line, space, words);
    if ( ( words.size()>5 )
        && ( words[0]=="SHRINK." )
        && ( words[1]=="FACT.(MONKH.)" ) ) {
      for (int is = 0; is < 3; is++ ) {
        shrink_fact[is]=atoi(words[2+is].c_str());
      }
      break;
    }
  }
  //cout << "shrinking factor " << shrink_fact << endl; 

  is.clear();
  totmo=readMO(is, eigen_start[kpt], moCoeff);
  if(calctype=="UHF") {
    slwriter.spin_dwn_start=moCoeff.size();
    totmo=readMO(is, eigen_start[kpt+nkpts], moCoeff);
  }
  cout << "nmo's " << moCoeff.size() << endl;

  slwriter.kpoint.resize(3);
  vector<string> kwords;
  kpoints[kpt].erase(kpoints[kpt].find(')'));
  split(kpoints[kpt], space, kwords);
  assert(kwords.size()>=5);
  for(int i=0; i< 3; i++) {
    slwriter.kpoint[i]=atoi(kwords[i+2].c_str());
  }
  for(int i=0; i< 3; i++) 
    slwriter.kpoint[i]/=shrink_fact[i]/2.0;


  cout << "chosen k-point " << slwriter.kpoint[0] <<"   " 
    << slwriter.kpoint[1] << "   " << slwriter.kpoint[2] << endl;


  //cout << "Found " << moCoeff.size() << " MO's. and "
  //<< moCoeff[0].size() << " functions " << endl;
  if(totfunctions != (int) moCoeff[0].size()) {
    cout << "The number of basis functions doesn't match between what was"
      << "read(" << moCoeff[0].size() << ") and calculated("
      << totfunctions << ")" << endl;
    exit(1);
  }
  if(totmo!=moCoeff.size()) {
    cout << "in make_orb, totmo is " << totmo
      << " and moCoeff.size() is " << moCoeff.size()
      << ".  They should be equal." << endl;
  }

  // analysis of band character and NORMALIZATION(!) of coefficients
  MO_analysis(is, atoms, slwriter, basis, origin,
      latvec, moCoeff, shift, totmo);

  //if our k-point isn't zero, we need to fix any shifts
  int f=0;
  Shifter shiftobj;
  shiftobj.origin=origin;

  //------------------------------------
  //Shift the atoms away from the edges.  Should
  //reduce the number of centers outside the cell most of the time.

  //vector <double> shift;
  //double shift_amount=0;
  //shift.push_back(shift_amount); shift.push_back(shift_amount); 
  //shift.push_back(shift_amount);
  if(latvec.size() > 0) { 
    vector <double> atomshift;
    for(int i=0; i< 3; i++) atomshift.push_back(0);
    for(int i=0; i< 3; i++) {
      for(int j=0; j< 3; j++) {
        atomshift[j]+=shift[i]*latvec[i][j];
      }
    }

    for(int i=0; i < natoms; i++) {
      atoms[i].pos=atoms[i].pos+atomshift;
    }
    //Now enforce the pbc's.. 

    for(unsigned int at=0; at< atoms.size(); at++) {
      vector <int> nshift;
      int bas=atoms[at].basis;
      int nfunc=basis[bas].nfunc();
      if(shiftobj.enforcepbc(atoms[at].pos, latvec, nshift)) {
        //cout << "at " << at << "  shifted " << nshift[0] << "  " << nshift[1] 
        //    << "   " << nshift[2] << endl;
        dcomplex kdots=0;
        for(int d=0; d< 3; d++) 
          kdots+=slwriter.kpoint[d]*nshift[d];
        //kdots=cos(pi*kdots);
        kdots=exp(pi*kdots*dcomplex(0.0,1.0));
        //cout << "kdots " << kdots << endl;
        for(int i=0; i< nfunc; i++) {
          for(int mo=0; mo < totmo; mo++) {
            moCoeff[mo][f]*=kdots;
          }
          f++;
        }
      }
      else f+=nfunc;
    }
  }


}

void read_crystal_orbital_all(istream & is,
    string & fort10file,
    vector <Atom> &atoms,
    Slat_wf_writer & slwriter,
    vector <Gaussian_basis_set> & basis,
    vector <double> & origin,
    vector < vector <double> > & latvec,
    vector < vector < vector <double> > > & moCoeff, 
    vector <string> & kpoints, 
    vector <long int> & eigen_start, 
    vector < vector < double > > &kptCoord, 
    vector <double> & shift, 
    vector <int> shifted, 
    vector < vector <int> > nshift) {

  assert(atoms.size() > 0);
  assert(basis.size() > 0);
  string dummy;
  string calctype=slwriter.calctype;
  //vector <int> nbasis;
  //vector <double > emptyVector; //to push onto moCoeff to add a mo.

  //Find out how many functions should be there.
  int natoms=atoms.size();
  int totfunctions=0;
  for(int at=0; at < natoms; at++ ) {
    int bas=atoms[at].basis;
    int nfunc=basis[bas].nfunc();
    totfunctions+=nfunc;
  }
  //cout << "Should be " << totfunctions << " functions " << endl;

  string space=" ";
  int totmo=0;

  int nkpts=eigen_start.size();
  if(calctype=="UHF") {
    assert(nkpts%2==0);
    nkpts/=2;
  }
  moCoeff.resize(nkpts);
  kpoints.resize(nkpts);
  kptCoord.resize(nkpts);
  cout << "Found " << nkpts << " k-points with real eigenvectors " << endl;
  for(int i=0; i< nkpts; i++) {
    cout << i << " : " << kpoints[i] << endl;
  }
  int kpt; 

  is.clear();
  is.seekg(1);
  string line;
  int shrink_fact[3];
  for(int i=0; i<3; i++) shrink_fact[i]=1.0;
  while ( getline(is,line) ) {
    vector<string> words;
    split(line, space, words);
    if ( ( words.size()>5 )
        && ( words[0]=="SHRINK." )
        && ( words[1]=="FACT.(MONKH.)" ) ) {
      for (int is=0; is<3;is++) 
        shrink_fact[is]=atoi(words[2+is].c_str());
      break;
    }
  }
  for (kpt=0; kpt<nkpts; kpt++) {
    is.clear();
    cout << "Reading orbital for " << kpoints[kpt] << "..." << endl; 
    totmo=readMO(is, eigen_start[kpt], moCoeff[kpt]);
    if(calctype=="UHF") {
      slwriter.spin_dwn_start=moCoeff[kpt].size();
      totmo=readMO(is, eigen_start[kpt+nkpts], moCoeff[kpt]);
    }
    // cout << "nmo's " << moCoeff[kpt].size() << endl;

    slwriter.kpoint.resize(3);
    for(int i=0; i<3; i++) slwriter.kpoint[i]=0.0;
    vector<string> kwords;
    kpoints[kpt].erase(kpoints[kpt].find(')'));
    split(kpoints[kpt], space, kwords);
    if(kwords.size()<5) {
      cout << "Not enough words in " << kpoints[kpt] << endl;
      exit(1);
    }
    
    double max=0;
    for(int i=0; i< 3; i++) {
      slwriter.kpoint[i]=atoi(kwords[i+2].c_str());
      //if(slwriter.kpoint[i] > max) max=slwriter.kpoint[i];
      slwriter.kpoint[i]/=shrink_fact[i]/2.;
    }
    kptCoord[kpt].resize(3);
    for(int i=0; i< 3; i++) {
      kptCoord[kpt][i] = slwriter.kpoint[i]; 
    }


    //cout << "Found " << moCoeff[kpt].size() << " MO's. and "
    //<< moCoeff[kpt][0].size() << " functions " << endl;
    if(totfunctions != (int) moCoeff[kpt][0].size()) {
      cout << "The number of basis functions doesn't match between what was"
        << "read(" << moCoeff[kpt][0].size() << ") and calculated("
        << totfunctions << ")" << endl;
      exit(1);
    }
    if(totmo!=moCoeff[kpt].size()) {
      cout << "in make_orb, totmo is " << totmo
        << " and moCoeff[kpt].size() is " << moCoeff[kpt].size()
        << ".  They should be equal." << endl;
    }


    fort10input(is, fort10file, atoms, slwriter, basis, origin,
        latvec, moCoeff[kpt], shift);

    // analysis of band character and NORMALIZATION(!) of coefficients
    string mo_filename="mo_analysis";
    for(vector<double>::iterator d=slwriter.kpoint.begin(); d!=slwriter.kpoint.end(); d++)
      append_number(mo_filename,*d);
    MO_analysis(is, fort10file, atoms, slwriter, basis, origin,
        latvec, moCoeff[kpt], shift, totmo,mo_filename);


    //if our k-point isn't zero, we need to fix any shifts
    int f=0;
    Shifter shiftobj;
    shiftobj.origin=origin;

    //------------------------------------
    //Shift the atoms away from the edges.  Should
    //reduce the number of centers outside the cell most of the time.

    //vector <double> shift;
    //double shift_amount=0;
    //shift.push_back(shift_amount); shift.push_back(shift_amount); 
    //shift.push_back(shift_amount);
    if(latvec.size() > 0) { 

      //Now enforce the pbc's.. 

      for(unsigned int at=0; at< atoms.size(); at++) {
        int bas=atoms[at].basis;
        int nfunc=basis[bas].nfunc();
        if(shifted[at]) {
          //cout << "at " << at << "  shifted " << nshift[0] << "  " << nshift[1] 
          //    << "   " << nshift[2] << endl;
          double kdots=0;
          for(int d=0; d< 3; d++) 
            kdots+=slwriter.kpoint[d]*nshift[at][d];
          kdots=cos(pi*kdots);
          //cout << "kdots " << kdots << endl;
          for(int i=0; i< nfunc; i++) {
            for(int mo=0; mo < totmo; mo++) {
              moCoeff[kpt][mo][f]*=kdots;
            }
            f++;
          }
        }
        else f+=nfunc;
      }
    }
  } 
}



void read_crystal_orbital_all(istream & is,
    vector <Atom> &atoms,
    Slat_wf_writer & slwriter,
    vector <Gaussian_basis_set> & basis,
    vector <double> & origin,
    vector < vector <double> > & latvec,
    vector <vector < vector <dcomplex> > > & moCoeff,
    vector <string> & kpoints, 
    vector <long int> & eigen_start, 
    vector < vector < double > > &kptCoord, 
    vector <double> & shift, 
    vector <int> shifted, 
    vector < vector <int> > nshift) {

  assert(atoms.size() > 0);
  assert(basis.size() > 0);
  string dummy;
  string calctype=slwriter.calctype;
  //vector <int> nbasis;
  //Find out how many functions should be there.
  int natoms=atoms.size();
  int totfunctions=0;
  for(int at=0; at < natoms; at++ ) {
    int bas=atoms[at].basis;
    int nfunc=basis[bas].nfunc();
    totfunctions+=nfunc;
  }
  //cout << "Should be " << totfunctions << " functions " << endl;

  string space=" ";
  int totmo=0;
  int nkpts=eigen_start.size();
  if(calctype=="UHF") {
    assert(nkpts%2==0);
    nkpts/=2;
  }
  kpoints.resize(nkpts);
  cout << "Found " << nkpts << " k-points with complex eigenvectors " << endl;
  for(int i=0; i< nkpts; i++) {
    cout << i << " : " << kpoints[i] << endl;
  }
  int kpt;
  moCoeff.resize(nkpts);
  kptCoord.resize(nkpts);

  is.clear();
  is.seekg(1);
  string line;
  int shrink_fact[3];
  for(int i=0; i<3; i++) shrink_fact[i]=1.0;
  while ( getline(is,line) ) {
    vector<string> words;
    split(line, space, words);
    if ( ( words.size()>5 )
        && ( words[0]=="SHRINK." )
        && ( words[1]=="FACT.(MONKH.)" ) ) {
      for (int is = 0; is<3; is++) {
        shrink_fact[is]=atoi(words[2+is].c_str());
      }
      break;
    }
  }
  //  cout << "shrinking factor: " << shrink_fact[0] << "  " << shrink_fact[1] << "  " << shrink_fact[2] << endl; 
  for (kpt = 0; kpt < nkpts; kpt++) {

    //cout << "shrinking factor " << shrink_fact << endl; 

    is.clear();
    cout << "Reading orbital for " << kpoints[kpt] << "..." << endl; 
    totmo=readMO(is, eigen_start[kpt], moCoeff[kpt]);
    if(calctype=="UHF") {
      slwriter.spin_dwn_start=moCoeff[kpt].size();
      totmo=readMO(is, eigen_start[kpt+nkpts], moCoeff[kpt]);
    }
    //cout << "nmo's " << moCoeff[kpt].size() << endl;

    slwriter.kpoint.resize(3);
    for(int i=0;i < 3; i++) 
      slwriter.kpoint[i]=0.0;

    vector<string> kwords;
    kpoints[kpt].erase(kpoints[kpt].find(')'));
    split(kpoints[kpt], space, kwords);
    if(kwords.size()<5) {
      cout << "Not enough words in " << kpoints[kpt] << endl;
    }
    for(int i=0; i< 3; i++) {
      slwriter.kpoint[i]=atoi(kwords[i+2].c_str());
    }
    kptCoord[kpt].resize(3);
    for(int i=0; i< 3; i++) {
      slwriter.kpoint[i]/=shrink_fact[i]/2.0;
      kptCoord[kpt][i] =  slwriter.kpoint[i];
    }


    //cout << "Found " << moCoeff[kpt].size() << " MO's. and "
    //<< moCoeff[kpt][0].size() << " functions " << endl;
    if(totfunctions != (int) moCoeff[kpt][0].size()) {
      cout << "The number of basis functions doesn't match between what was"
        << "read(" << moCoeff[kpt][0].size() << ") and calculated("
        << totfunctions << ")" << endl;
      exit(1);
    }
    if(totmo!=moCoeff[kpt].size()) {
      cout << "in make_orb, totmo is " << totmo
        << " and moCoeff[kpt].size() is " << moCoeff[kpt].size()
        << ".  They should be equal." << endl;
    }

    // analysis of band character and NORMALIZATION(!) of coefficients
    string mo_filename="mo_analysis";
    for(vector<double>::iterator d=slwriter.kpoint.begin(); d!=slwriter.kpoint.end(); d++)
      append_number(mo_filename,*d);
    MO_analysis(is, atoms, slwriter, basis, origin,
        latvec, moCoeff[kpt], shift, totmo,mo_filename);

    //if our k-point isn't zero, we need to fix any shifts
    int f=0;
    Shifter shiftobj;
    shiftobj.origin=origin;

    //------------------------------------
    //Shift the atoms away from the edges.  Should
    //reduce the number of centers outside the cell most of the time.

    //vector <double> shift;
    //double shift_amount=0;
    //shift.push_back(shift_amount); shift.push_back(shift_amount); 
    //shift.push_back(shift_amount);
    if(latvec.size() > 0) { 
      //Now enforce the pbc's.. 

      for(unsigned int at=0; at< atoms.size(); at++) {

        int bas=atoms[at].basis;
        int nfunc=basis[bas].nfunc();
        if(shifted[at]) {
          //cout << "at " << at << "  shifted " << nshift[0] << "  " << nshift[1] 
          //    << "   " << nshift[2] << endl;
          dcomplex kdots=0;
          for(int d=0; d< 3; d++) 
            kdots+=slwriter.kpoint[d]*nshift[at][d];
          //kdots=cos(pi*kdots);
          kdots=exp(pi*kdots*dcomplex(0.0,1.0));
          //cout << "kdots " << kdots << endl;
          for(int i=0; i< nfunc; i++) {
            for(int mo=0; mo < totmo; mo++) {
              moCoeff[kpt][mo][f]*=kdots;
            }
            f++;
          }
        }
        else f+=nfunc;
      }
    }
  }
}

void read_kpt_eigenpos(istream & is,
    vector <string> & rkpoints,
    vector <long int> & reigen_start, 
    vector <string> & ckpoints,
    vector <long int> & ceigen_start) {

  string dummy; 
  string space = " "; 
  while(is >> dummy) {
    //cout << "dummy " << dummy << endl;
    if(dummy == "FINAL") {
      is >> dummy;
      if(dummy == "EIGENVECTORS" ) {
        is.ignore(125, '\n'); //clear the line with FINAL EIG..
        string line;
        getline(is, line);
        // cout << "line " << line << endl;
        while(getline( is,line)) {
          if(line.size() > 15 && line[5]=='(' && line[15]==')') {
            //cout << line[5] << "  " << line[15] << endl;
            is.ignore(150, '\n');//two blank lines
            is.ignore(150, '\n');
            long int pos=is.tellg();
            string line2;
            getline(is, line2);
            vector <string> words;
            split(line2, space, words);
            if(words[0]==words[1]) {
              //cout << "complex  k-point line " << line << endl;
              ckpoints.push_back(line);
              ceigen_start.push_back(pos);
            } else { 
              rkpoints.push_back(line);
              reigen_start.push_back(pos);
            }
          }
        }
      }
    }
    else if(dummy == "NEWK") {
      is >> dummy;
      if(dummy == "EIGENVECTORS" ) {
        is.ignore(125, '\n'); //clear the line with FINAL EIG..
        string line;
        //getline(is, line);
        // cout << "line " << line << endl;
        while(getline( is,line)) {
          if(line.size() > 45 && line[28]=='K' && line[29]=='='&& line[35]=='(' && line[45]==')') {
            //cout << line[5] << "  " << line[15] << endl;
            is.ignore(150, '\n');//one blank lines
            long int pos=is.tellg();
            string line2;
            getline(is, line2);
            vector <string> words1, words;
            split(line, space, words1);
            string kstr = words1[4] + " ( " + words1[6] + "  " + words1[7] + " " + words1[8]; 
            split(line2, space, words);
            //	    cout << words[0] << " " << words[1] << endl; 
            if(words[0]==words[1]) {
              //cout << "complex  k-point line " << line << endl;
              ckpoints.push_back(kstr);
              ceigen_start.push_back(pos);
            } else { 
              rkpoints.push_back(kstr);
              reigen_start.push_back(pos);
            }
          }
        }
      }
    }
  }
}

// ===========================================================================
