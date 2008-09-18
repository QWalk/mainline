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
#include <cstdlib>
#include "wf_writer.h"
using namespace std;
void Slat_wf_writer::print_wavefunction(ostream & inputfile ) {
  int nmo=0;
  if(calctype == "RHF" || calctype == "ROHF") {
    nmo=max(nup, ndown);
  }
  else if(calctype == "UHF") {
    nmo=spin_dwn_start+ndown;
  }
  else if(calctype=="GVB") {
    assert(detwt.size() > 0);
    nmo=0;
    for(unsigned int i=0; i< occ_up.size(); i++) {
      for(unsigned int j=0; j < occ_up[i].size(); j++)
        if(nmo < occ_up[i][j]) nmo=occ_up[i][j];
    }

    for(unsigned int i=0; i< occ_down.size(); i++) {
      for(unsigned int j=0; j < occ_down[i].size(); j++)
        if(nmo < occ_down[i][j]) nmo=occ_down[i][j];
    }
  }
  else {
    cout << "strange calctype " << calctype << endl;
    exit(1);
  }

  int nelectrons=nup+ndown;
  if(magnification < 0) 
    magnification=max(nelectrons/20.0, 1.0);

  inputfile << "  SLATER \n"
            << "  ORBITALS {\n"
            << "  " << mo_matrix_type << endl
            << "    MAGNIFY " << magnification << endl
            << "    INCLUDE " << basisname << endl
            << "    NMO " << nmo << endl
            << "    ORBFILE " << orbname << endl;

  if(write_centers) {
    inputfile << "    CENTERS { READ " << centername << " } " << endl;
  }
  else if(use_global_centers) {
    inputfile << "    CENTERS { USEGLOBAL } " << endl;
  }
  else {
    inputfile << "    CENTERS { USEATOMS } " << endl;
  }
  inputfile << "  }\n\n";

  if(detwt.size()==0) {
    //see below about the two-determinant UHF wave function
    //if(calctype=="UHF" && ndown==nup) 
    //  inputfile << "  DETWT { .707 .707 } \n";
    //else 
    inputfile << "  DETWT { 1.0 } \n";
    
    
    inputfile << "  STATES {\n    ";
    
    
    inputfile << " #Spin up orbitals\n   ";
    for(int i=1; i< nup+1; i++) {
      inputfile << i << " ";
      if(i%20==0) inputfile << "\n    ";
    }
    int downstart=1;
    if(calctype=="UHF") {
      downstart=spin_dwn_start+1;
    }
    inputfile << "\n    ";
    
    inputfile << " #Spin down orbitals \n   ";
    for(int i=downstart; i < downstart+ndown; i++) {
      inputfile << i << " ";
      if((i-downstart+1)%20==0) inputfile << "\n    ";
    }

    //This doesn't necessarily help, and can sometimes hurt, so
    //we don't put it in.
    //if(calctype=="UHF" && ndown==nup) {
    //  inputfile <<"\n #for UHF, we put in two determinants \n";
    //  inputfile << "  #Spin up orbitals\n";
    //  for(int i=downstart; i < downstart+ndown; i++) {
    //    inputfile << i << " ";
    //    if((i-downstart+1)%20==0) inputfile << "\n    ";
    //  }
    //  inputfile << endl;
    //  inputfile << " #Spin down orbitals\n   ";
    //  for(int i=1; i< nup+1; i++) {
    //    inputfile << i << " ";
    //    if(i%20==0) inputfile << "\n    ";
    //  }
    //}

    inputfile << "  }\n";
  }
  else {
    int ndet=detwt.size();
    assert(occ_up.size()==detwt.size());
    assert(occ_down.size()==detwt.size());
    inputfile << "DETWT { ";
    for(int d=0; d< ndet; d++) {
      inputfile << detwt[d] << "  ";
    }
    inputfile << " } " << endl;
    
    inputfile << " STATES { \n";
    for(int d=0; d< ndet; d++) {
      assert(occ_up[d].size()==nup);
      inputfile << "# Up spin \n";
      for(int i=0; i < nup; i++) {
        inputfile << occ_up[d][i] << "  ";
      }
      inputfile << endl;

      assert(occ_down[d].size()==ndown);
      inputfile << "# Down spin\n";
      for(int i=0; i< ndown; i++) {
        inputfile <<  occ_down[d][i] << "  ";
      }
      inputfile << endl;
    }
    inputfile << "} " << endl;
  }
  
}


//######################################################################

#include "converter.h"
#include "basis_writer.h"

void Jastrow_wf_writer::set_atoms(vector <Atom> & atoms) {
  int natoms=atoms.size();
  for(int at=0; at < natoms; at++) {
    bool unique=true;
    int nunique=atomnames.size();
    for(int i=0; i< nunique; i++) {
      if(atoms[at].name==atomnames[i]) {
        unique=false;
        break;
      }
    }
    if(unique) {
      atomnames.push_back(atoms[at].name);
    }
  }
}


void Jastrow_wf_writer::add_basis(Basis_writer & b) {
  basis.push_back(&b);
}

void Jastrow_wf_writer::print_wavefunction(ostream & os) {
  os << "JASTROW \n\n";
  for(vector <Basis_writer *>:: iterator i=basis.begin();
      i!=basis.end(); i++) {
      (*i)->print_basis(os);
  }

  for(vector <string>::iterator i=atomnames.begin();
      i!=atomnames.end(); i++) {
      os << "EICORRELATION { " << *i << " 0.0 0.0 0.0 }\n\n";
      os << "EEICORRELATION { " << *i << endl
         << "  1 1 0   0.\n"
         << "  1 0 1   0.\n"
         << "  1 1 1   0.\n"
         << "  2 2 0   0.\n"
         << "  2 0 1   0.\n"
         << "  2 0 2   0.\n"
         << "  2 2 2   0.\n"
         << "  3 3 0   0.\n"
         << "  3 0 2   0.\n"
         << "  3 3 2   0.\n"
         << "  1 2 2   0.\n"
         << "  2 3 2   0.\n"
         << "}\n\n";
  }

}


//#######################################################################
void Jastrow2_wf_writer::set_atoms(vector <Atom> & atoms) {
  int natoms=atoms.size();
  for(int at=0; at < natoms; at++) {
    bool unique=true;
    int nunique=atomnames.size();
    for(int i=0; i< nunique; i++) {
      if(atoms[at].name==atomnames[i]) {
        unique=false;
        break;
      }
    }
    if(unique) {
      atomnames.push_back(atoms[at].name);
    }
  }
}


void Jastrow2_wf_writer::add_ee_basis(Basis_writer & b, int group) {
  assert(group < ngroups);
  ee_basis[group].push_back(&b);
}

void Jastrow2_wf_writer::add_ei_basis(Basis_writer & b, int group) {
  assert(group < ngroups);
  ei_basis[group].push_back(&b);
}



void Jastrow2_wf_writer::print_wavefunction(ostream & os) {
  os << "JASTROW2 \n\n";
  assert(ee_basis.size()==ei_basis.size());
  assert(ee_basis.size()==ngroups);
  //cout << "ngroups " << ngroups << endl;
  string indent= "   ";

  for(int g=0; g< ngroups; g++) {
    os << "GROUP { \n";
    os << "  OPTIMIZEBASIS " << endl;
    for(vector <Basis_writer *>:: iterator i=ei_basis[g].begin();
        i!=ei_basis[g].end(); i++) {
        os << "EIBASIS { " << endl;
        (*i)->print_basis(os);
        os << "}" << endl;
    }

    if(ei_basis[g].size() > 0) {
    os << "ONEBODY { " << endl;

    for(vector <string>::iterator i=atomnames.begin();
        i!=atomnames.end(); i++) {
        int nfunc=0;
        for(unsigned int bas=0;
          bas< ei_basis[g].size(); bas++) {
          if(*i == ei_basis[g][bas]->label) {
            nfunc+=ei_basis[g][bas]->nfunc();
          }
        }

        if(nfunc > 0) {
          os << "  COEFFICIENTS { " << *i << "  ";
          for(int i=0; i< nfunc; i++)
            os << " 0.0  ";
          os << " } " << endl;
        }
    }
    os << " } " << endl;
    }

    for(vector <Basis_writer*> :: iterator i=ee_basis[g].begin();
        i!=ee_basis[g].end(); i++) {
        os << "EEBASIS { " << endl;
        (*i)->print_basis(os);
        os << "}" << endl;
    }

    if(spindiff[g]) {
      assert(ee_basis[g].size()==2);

      os << "TWOBODY_SPIN { \n";
      os << "  FREEZE \n";
      os << "  LIKE_COEFFICIENTS { 0.25  0.0 } \n";
      os << "  UNLIKE_COEFFICIENTS { 0.0 0.5 } \n";
      os << " }\n";
    }
    else {
      int nfunc=0;
      for(vector <Basis_writer *>:: iterator i=ee_basis[g].begin();
        i!=ee_basis[g].end(); i++) {
        nfunc+=(*i)->nfunc();
      }
      if(nfunc > 0) {
        os << "TWOBODY { " << endl;
        os << "  COEFFICIENTS { " ;
        for(int j=0; j< nfunc; j++)
          os << " 0.0  ";
        os << " }" << endl;
        os << " } " << endl;
      }
    }
    os << "} " << endl;

  }

}


//----------------------------------------------------------------------

void print_std_jastrow2(Jastrow2_wf_writer & jast2writer, ostream & os,
                        double basis_cutoff) {
  jast2writer.set_groups(2);
  jast2writer.set_spin_diff(0);

  vector <Cutoff_cusp_basis> eecusp; //Group 0 is the cusp
  eecusp.resize(2);
  eecusp[0].cutoff=eecusp[1].cutoff=basis_cutoff;
  jast2writer.add_ee_basis(eecusp[0], 0);
  jast2writer.add_ee_basis(eecusp[1], 0);

  //ugly hack..  The jast writer saves pointers to the basis functions,
  //but when we modify a vector, those pointers can change, so
  //we have to add all the basis functions, then assign them.
  //This is pretty unstable..there has to be a better way of doing
  //it without losing polymorphism.

  vector <Poly_pade_basis> eibas; //group 2 is the regular basis
  for(unsigned int i=0; i< jast2writer.atomnames.size(); i++) {
    Poly_pade_basis tbas;
    tbas.label=jast2writer.atomnames[i];
    tbas.cutoff=basis_cutoff;
    tbas.beta0=.2;
    tbas.nfunc_=2;
    eibas.push_back(tbas);

  }
  for(unsigned int i=0; i< eibas.size(); i++) {
   jast2writer.add_ei_basis(eibas[i], 1);
  }
  Poly_pade_basis eebas;
  eebas.cutoff=basis_cutoff;
  eebas.beta0=.5;
  eebas.nfunc_=1;
  jast2writer.add_ee_basis(eebas,1);
  jast2writer.print_wavefunction(os);
}


//----------------------------------------------------------------------
//this is done in a much more straightforward way than the Jastrow2_writer
//I don't think that we need all the polymorphism, etc, that we had before
void print_3b_jastrow2(ostream & os, std::vector<std::string> & unique_atoms, double cutoff) { 
  os << "JASTROW2" << endl;
  os << "GROUP { \n"
    << " OPTIMIZEBASIS\n"
    << " EEBASIS { EE CUTOFF_CUSP GAMMA 24.0 CUSP 1.0 CUTOFF " << cutoff << " } \n"
    << " EEBASIS { EE CUTOFF_CUSP GAMMA 24.0 CUSP 1.0 CUTOFF " << cutoff << " } \n"
    << " TWOBODY_SPIN {  FREEZE \n"
    << "   LIKE_COEFFICIENTS { 0.25  0.0 } \n"
    << "   UNLIKE_COEFFICIENTS { 0.0 0.5 } \n }\n"
    << "}\n";
  os << "GROUP { \n OPTIMIZEBASIS\n"
    <<" EEBASIS { EE POLYPADE BETA0 0.5 NFUNC 3 RCUT " << cutoff << " }\n";
  for(vector<string>::iterator i=unique_atoms.begin(); i!= unique_atoms.end();
      i++) { 
    os << " EIBASIS { " << *i << " POLYPADE BETA0 0.2 NFUNC 4 RCUT " << cutoff << " } \n";
  }
  os << " ONEBODY { \n";
  for(vector<string>::iterator i=unique_atoms.begin(); i!= unique_atoms.end();
      i++) { 
    os << "  COEFFICIENTS { " << *i << " 0 0 0 0  } \n";
  }
  os << " }\n";
  os << " TWOBODY { \n";
  os << "  COEFFICIENTS { 0 0 0 } \n";
  os << " }\n";
  os << " THREEBODY { \n";
  for(vector<string>::iterator i=unique_atoms.begin(); i!= unique_atoms.end();
      i++) { 
    os << "  COEFFICIENTS { " << *i << " 0 0 0 0 0 0 0 0  0 0 0 0 } \n";
  }
  os << " }\n";
  os << "}\n";
  
}
//----------------------------------------------------------------------
void fold_kpoint(Slat_wf_writer & slwriter, 
                 std::vector <std::vector <double> > & latvec,
                 int dir,
                 std::vector <std::vector < double> > & moCoeff,
                 std::vector <Atom> & atoms) { 
  //Note: this won't work for UHF wavefunctions..
  vector <Atom> natoms;
  for(vector<Atom>::iterator i=atoms.begin(); i!=atoms.end(); i++) { 
    Atom tmp=*i;
    for(int d=0;d < 3; d++) { 
      tmp.pos[d]+=latvec[dir][d];
    }
    natoms.push_back(tmp);
  }
  atoms.insert(atoms.end(), natoms.begin(), natoms.end());

  vector <vector <double> > nmocoeff;
  vector<double> motmp;
  int count=1;
  for(vector<vector<double> >::iterator i=moCoeff.begin(); i!=moCoeff.end();
      i++) { 
    motmp.clear();
    //gamma point
    for(vector<double>::iterator j=i->begin(); j!= i->end(); j++) { 
      motmp.push_back(*j);
    }
    for(vector<double>::iterator j=i->begin(); j!= i->end(); j++) { 
      motmp.push_back(*j);
    }
    nmocoeff.push_back(motmp);
    motmp.clear();
    cout << "count " << count << " "  << nmocoeff.size() << " ";
    //X point
    for(vector<double>::iterator j=i->begin(); j!= i->end(); j++) { 
      motmp.push_back(*j);
    }
    for(vector<double>::iterator j=i->begin(); j!= i->end(); j++) { 
      motmp.push_back(-*j);
    }
    nmocoeff.push_back(motmp);
    cout << " " << nmocoeff.size() << endl;
    count++;
  }
  moCoeff=nmocoeff;
  
  for(int d=0; d< 3; d++) {
    latvec[dir][d]*=2;
  }
  
}

//----------------------------------------------------------------------
