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
#include <string>
#include <fstream>
#include <cmath>

#include "basis_writer.h"
#include "Pseudo_writer.h"
#include "wf_writer.h"
using namespace std;

void read_lattice_vector(istream & is, vector <string> & currline, 
                         vector <vector <double> > & latvec);
void read_atoms(istream & is, vector <string> & currline, 
                         vector <Atom> & atoms);
void read_mo_coefficients(istream & is, vector <string> & currline,
                          Slat_wf_writer & slwriter, vector <vector <double> > & moCoeff);
//reads in the basis for each type of atom
void read_basis(vector <Atom> & atoms, vector<Spline_basis_writer> & basis );
void read_psp(vector <Atom> & atoms, vector <Spline_pseudo_writer> & pseudo);
//###########################################################################


int main(int argc, char ** argv) { 
  if(argc < 2) {
    cout << "usage: siesta2qmc [siesta stdout]\n";
    exit(1);
  }
  string infilename=argv[1];
  string outputname="qwalk";
  
  //------------------------------read in the variables from the output file
  vector <Atom> atoms;
  vector <vector < double> > latvec;
  Slat_wf_writer slwriter;
  vector <vector <double> >  moCoeff;
  double spin_pol=0;
  vector <Spline_pseudo_writer> pseudo;
  
  
  ifstream is(infilename.c_str());
  if(!is) { 
    cout << "Couldn't open " << infilename << endl;
    exit(1);
  }
  string line;
  string space=" ";
  vector <string> currline;
  while(getline(is,line)) { 
    currline.clear();
    split(line, space, currline);
    read_atoms(is,currline,atoms);
    read_lattice_vector(is,currline,latvec);
    read_mo_coefficients(is, currline, slwriter, moCoeff);
    if(currline.size()> 6 && currline[2]=="spin" && currline[3]=="polarization") { 
      spin_pol=atof(currline[6].c_str());
      cout << "spin polarization " << spin_pol << endl;
    }
  }
  is.close();
  slwriter.orbname=outputname+".orb";
  slwriter.basisname=outputname+".basis";
  slwriter.mo_matrix_type="CUTOFF_MO";
  //TEMPORARY!!!
  slwriter.nup=3; slwriter.ndown=1;
  
  
  vector <Spline_basis_writer> basis;
  read_basis(atoms, basis);
  read_psp(atoms, pseudo);
  //---------------------------------------------
  //--Write out all the collected data
  //---------------------------------------------
  
  
  //---------------------Basis file
  ofstream os (slwriter.basisname.c_str());
  os.precision(15);
  for(vector<Spline_basis_writer>::iterator i=basis.begin();
      i!= basis.end(); i++) {
    os << "BASIS { \n";
    i->print_basis(os);
    os << "}\n\n";
  }
  os.close(); os.clear();
  
  //-----------------ORB file
  
  os.open(slwriter.orbname.c_str());
  //this is really redundant with gamess2qmc.cpp..should probably be refactored somehow
  vector < Center> centers;
  vector <int> nbasis;
  centers.resize(atoms.size());
  int natoms=atoms.size();
  nbasis.resize(natoms);
  for(int at=0; at < natoms; at++) {
    for(int i=0; i< 3; i++) centers[at].pos[i]=atoms[at].pos[i];
    centers[at].equiv_atom=at;
    centers[at].name=atoms[at].name;
    nbasis[at]=basis[atoms[at].basis].nfunc();
  }
  print_orbitals(os,centers,nbasis, moCoeff);
  os.close(); os.clear();
  
  //--------------------Slater determinant file
  
  string slateroutname=outputname+".slater";
  os.open(slateroutname.c_str());
  slwriter.print_wavefunction(os);
  os.close(); os.clear();
  
  
  //---------------------sys file
  //May want to add a command-line switch to print out a molecule
  //also, it should choose the cutoff divider properly
  string sysoutname=outputname+".sys";
  ofstream sysout(sysoutname.c_str());
  sysout << "SYSTEM { PERIODIC \n";
  sysout << "  NSPIN { " << slwriter.nup << "  "
    << slwriter.ndown << " } \n";
  for(int at=0; at <natoms; at++) {
    atoms[at].print_atom(sysout);
  }
  cout << "latvec " << latvec.size() << endl; 
  sysout << "LATTICEVEC { \n";
  for(int i=0; i< 3; i++) { 
    cout << "size " << latvec[i].size() << endl;
    for(int j=0; j< 3; j++) sysout << latvec[i][j] << " ";
    sysout << endl;
  }
  sysout << "  }\n";

  sysout << "}\n\n\n";
  
  int npsp=pseudo.size();
  for(int psp=0; psp < npsp; psp++) {
    pseudo[psp].print_pseudo(sysout);
  }
  sysout.close();
  
  //--------------------------Jastrow 2 output

  double basis_cutoff=find_basis_cutoff(latvec);
  string jast2outname=outputname+".jast2";
  Jastrow2_wf_writer jast2writer;
  jast2writer.set_atoms(atoms);
  
  
  ofstream jast2out(jast2outname.c_str());
  print_std_jastrow2(jast2writer, jast2out, basis_cutoff);
  jast2out.close();
  
}

//###########################################################################

void read_atoms(istream & is, vector <string> & currline, 
                vector <Atom> & atoms) {
  string space=" ";
  if(currline.size() > 3 && currline[0]=="outcoor:" && currline[2]=="coordinates") { 
    atoms.clear(); //so we always take the last set of atoms.
    double fac=1;
    if(currline[3]=="(Ang):") { 
      cout << "converting from angstrom" << endl;
      fac=1/.529177249;
    }
    else if(currline[3]!="(Bohr):") { 
      cout << "Don't understand format of the coordinates: " << currline[3] << endl;
      exit(1);
    }
    cout << "reading in coordinates" << endl;
    string line;
    getline(is,line);
    currline.clear(); split(line, space, currline);
    while(currline.size() > 0) { 
      Atom tmpatom;
      tmpatom.name=currline[4];
      //This is the *label*, not the charge, so that will be changed later(in the psp reading)
      tmpatom.charge=atoi(currline[3].c_str()); 
      for(int i=0; i< 3; i++) { 
        tmpatom.pos[i]=atof(currline[i].c_str())*fac;
      }
      atoms.push_back(tmpatom);
      tmpatom.print_atom(cout);
      getline(is,line);
      currline.clear(); split(line,space,currline);
      
    }
      
  }
}

//###########################################################################
void read_lattice_vector(istream & is, vector <string> & currline, 
                         vector <vector <double> > & latvec) {
  string space=" ";
  if(currline.size() > 3 && currline[0]=="outcell:" && currline[3]=="vectors") { 
    latvec.clear();
    double fac=1;
    if(currline[4]=="(Ang):") { 
      fac=1/.529177249;
    }
    else if(currline[4]!="(Bohr):") { 
      cout << "Don't understand units: " << currline[4] << endl;
    }
    cout << "reading in lattice vectors " << endl;
    string line;
    latvec.resize(3);
    for(int i=0; i< 3; i++) { 
      latvec[i].resize(3);
      currline.clear();
      getline(is,line); split(line, space,currline);
      for(int j=0; j< 3; j++) { 
        latvec[i][j]=atof(currline[j].c_str())*fac;
      }
    }
    cout << "lattice vectors " << endl;
    for(int i=0; i < 3; i++) { 
      for(int j=0; j< 3; j++) cout << latvec[i][j] << "  ";
      cout << endl;
    }
  }
}

//###########################################################################
void read_mo_coefficients(istream & is, vector <string> & currline,
                          Slat_wf_writer & slwriter, vector <vector <double> > & moCoeff) { 
  slwriter.use_global_centers=true;
  slwriter.write_centers=false;
  string space=" ";
  if(currline.size() > 3 && currline[0]=="writewave:" && currline[2]=="Functions") { 
    string line;
    getline(is, line); //empty line;
    getline(is, line); //k-points
    getline(is, line); //nspin
    currline.clear(); split(line, space, currline);
    cout << line << endl;
    int nspin=atoi(currline[4].c_str());
    if(nspin==1) { 
      slwriter.calctype="RHF";
    }
    else if(nspin==2) { 
      slwriter.calctype="UHF";
    }
    else { cout << "error reading nspin " << endl; exit(1); }
    
    getline(is,line); //'Number of basis orbs
    currline.clear(); split(line,space, currline);
    int nbasis=atoi(currline[5].c_str());
    cout << "calc type " << slwriter.calctype << "  nbasis " << nbasis << endl;

    getline(is,line); 
    
    for(int spin=0; spin < nspin; spin++) { 
      getline(is, line); getline(is,line); //one blank line and *'s
      getline(is,line); //k-point line, ignore for now
      getline(is,line); //spin component
      getline(is,line); //Num. wavefunctions
      currline.clear(); split(line, space, currline);
      int norb=atoi(currline[3].c_str());
      cout << "norb " << norb << endl;
      //Now start with reading the ****orbitals*** (KS wave functions)
      for(int orb=0; orb < norb; orb++) { 
        for(int i=0; i< 5; i++) getline(is,line); //blank, wf#, energy, ---'s, header
        vector <double> orb(nbasis);
        for(int b=0; b< nbasis; b++) { 
          getline(is, line);
          currline.clear(); split(line, space,currline);
          orb[b]=atof(currline[5].c_str());
        }
        moCoeff.push_back(orb);
        getline(is, line); //a line of --'s
      }
      if(spin==0 && nspin==2) { 
        slwriter.spin_dwn_start=moCoeff.size();
      }
    } //spin loop
    cout << "total norbs " << moCoeff.size() << endl;

  }
  
  
  //does not fill nup, ndown, mo_matrix_type, orbname, or basisname in slwriter
}

//###########################################################################
void find_unique_atoms(const vector<Atom> & atoms, vector<string> & unique_atoms) { 
  unique_atoms.clear();
  for(vector<Atom>::const_iterator at=atoms.begin(); at != atoms.end();
      at++) { 
    if(find(unique_atoms.begin(), unique_atoms.end(),at->name)==unique_atoms.end()){
      unique_atoms.push_back(at->name);
      cout << "unique atom " << at->name << endl;
    }
  }
}

//###########################################################################

//Make a uniform grid out of a 
void make_uniform(vector <double> & r, 
                  vector <double> & vals,
                  vector <double> & ur,
                  vector <double> & uvals) { 
  const double spacing=0.05; //bohr..
  double cutoff=10.0; //mostly for convienence..we work in terms of npts below
                      //this can be set quite high, since QWalk will find where to cut it off anyway
                      //(and fairly agressively)
  int npts=int(cutoff/spacing);
  
  int interval=0;
  for(int i=0; i< npts; i++) { 
    double rad=i*spacing;
    double val=0.0;
    if(i==0)  { 
      val=vals[0];
    }
    else { 
      while(r[interval+1] < rad && interval < r.size()-1) interval++;
      double x=(rad-r[interval])/(r[interval+1]-r[interval]);
      cout << "rad " << rad << " x " << x << " r1 " << r[interval] 
        << " r2 " << r[interval+1] << endl;
      val=(1-x)*vals[interval]+x*vals[interval+1];
    }
    ur.push_back(rad);
    uvals.push_back(val);
  }
  
}

//###########################################################################


void read_basis(vector <Atom> & atoms, vector<Spline_basis_writer> & basis ) {
  vector <string> unique_atoms;
  find_unique_atoms(atoms,unique_atoms);
  
  for(vector<string>::iterator nm=unique_atoms.begin();
      nm != unique_atoms.end(); nm++) { 
    Spline_basis_writer tmp_basis;
    tmp_basis.label=*nm;
    //do S orbitals..
    for(int el=1; el < 5; el++) { //search over possible l-values
      string prefix="ORB.S";
      append_number(prefix, el);
      prefix+=".";
      string postfix="."+*nm;
      int n=1;
      while(true) { 
        string file=prefix;
        append_number(file,n);
        file+=postfix;
        ifstream is(file.c_str());
        if(!is) break;
        cout << "found " << file << endl;
        switch(el) { 
          case 1:
            tmp_basis.types.push_back("S");
            break;
          case 2:
            tmp_basis.types.push_back("P");
            break;
          case 3:
            tmp_basis.types.push_back("D_siesta");
            break;
          default:
            cout << "Don't support this l-value.  Bug the maintainer." << endl;
            exit(1);
        }
        //Read in the file..
        is.ignore(180,'\n'); is.ignore(180,'\n');
        double rad, val;
        vector <double> rads, vals;
        while(is >> rad && is >> val) {
          rads.push_back(rad); vals.push_back(val);
        }
        vector <double> urad, uval;
        make_uniform(rads,vals,urad,uval);
        tmp_basis.rad.push_back(urad);
        tmp_basis.vals.push_back(uval);
        n++;
      }
    }
    basis.push_back(tmp_basis);
  }
  
  //We only have one basis object per atom name type, so 
  //assigning the basis number is pretty easy.
  for(vector <Atom>::iterator at=atoms.begin(); at!= atoms.end(); 
      at++) { 
    int nbasis=basis.size();
    for(int i=0; i< nbasis; i++) { 
      if(basis[i].label==at->name) { 
        at->basis=i;
      }
    }
  }
}





//###########################################################################


void read_psp(vector <Atom> & atoms, vector <Spline_pseudo_writer> & pseudo) {
  vector <string> unique_atoms;
  cout << "****Pseudopotential " << endl;
  string space=" ";
  find_unique_atoms(atoms,unique_atoms);
  for(vector<string>::iterator at=unique_atoms.begin();
      at != unique_atoms.end(); at++) { 
    Spline_pseudo_writer tmp_pseudo;
    string filename=*at+".psf";
    ifstream is(filename.c_str());
    if(!is) { 
      cout << "Couldn't open " << filename << endl;
      exit(1);
    }
    for(int i=0; i< 3; i++) is.ignore(180,'\n'); //first three lines are random bits..
    string line;
    getline(is, line);
    vector <string> spl;
    split(line, space, spl);
    cout << *at << " : effective charge : " << spl[5] <<  " integerized " << atoi(spl[5].c_str()) << endl;
    int zeff=atoi(spl[5].c_str());
    for(vector<Atom>::iterator i=atoms.begin(); i!= atoms.end();
        i++) { 
      if(i->name==*at) i->charge=zeff;
    }
    int npoints=atoi(spl[2].c_str());
    int nl=atoi(spl[0].c_str());
    int nl_up=atoi(spl[1].c_str());
    if(nl_up > 0) { 
      cout << "Can't deal with spin-polarized pseudopotentials..sorry\n";
      exit(1);
    }
    cout << "npoints " << npoints <<  " nl " << nl << endl;
    is.ignore(180,'\n'); //Radial grid..
    vector <double>  rad;
    for(int i=0; i< npoints; i++) { 
      double dum;
      is >> dum;
      rad.push_back(dum);
    }

    for(int l=0; l< nl; l++) { 
      vector <double> val;
      is.ignore(180,'\n');
      getline(is, line);
      cout << line << endl;
      int currl; is >> currl;
      cout << "l " << l << " "  << currl << endl;
      assert(l==currl);
      

      for(int i=0; i< npoints; i++) { 
        double dum; is >> dum; 
        val.push_back(dum);
        val[i]*= 0.5/rad[i];
      }
      vector <double> urad, uval;
      make_uniform(rad,val,urad,uval);
      
      string outname="psp";
      append_number(outname,l);
      ofstream out(outname.c_str());
      int n=urad.size();
      for(int i=0; i< n; i++) {
        out << urad[i] << "  " << uval[i] << endl;
      }
      
      tmp_pseudo.psp_pos.push_back(urad);
      tmp_pseudo.psp_val.push_back(uval);
      
    }
    tmp_pseudo.label=*at;
    int npts=tmp_pseudo.psp_val[0].size();
    for(int i=0; i< npts; i++) { 
      for(int l=0; l< nl-1; l++) { 
        tmp_pseudo.psp_val[l][i]-=tmp_pseudo.psp_val[nl-1][i];
      }
    }
    pseudo.push_back(tmp_pseudo);
  }
}


//###########################################################################

