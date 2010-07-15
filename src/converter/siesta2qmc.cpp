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
#include <sstream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include "basis_writer.h"
#include "Pseudo_writer.h"
#include "wf_writer.h"
#include "vecmath.h"

using namespace std;
enum coord_t { bohr, ang, scaled, fractional };

void read_lattice_vector(istream & is, vector <string> & currline, 
                         vector <vector <double> > & latvec);
void read_atoms(istream & is, vector <string> & currline, 
                         vector <Atom> & atoms,coord_t & coord_type);

void read_mo_coefficients(istream & is, vector <string> & currline,
                          Slat_wf_writer & slwriter, vector < vector <vector <dcomplex> > > & moCoeff,
			  vector < vector <double > > &kpoints);

//reads in the basis for each type of atom
void read_basis(vector <Atom> & atoms, vector<Spline_basis_writer> & basis );
void read_psp(vector <Atom> & atoms, vector <Spline_pseudo_writer> & pseudo);

void fix_basis_norm(vector <Atom> & atoms, vector <Spline_basis_writer> & basis,
                    vector < vector <dcomplex> > & moCoeff);

void orbital_split(string line, vector <string> & currline);
void reciprocal_lattice_vectors(vector < vector <double > > &latvec,
				vector < vector <double > > &recip_latvec);
void kpoint_to_frac_coords(vector < vector <double> > &recip_latvec,
			   vector < vector <double> > &kpoints,
			   vector < vector <double> > &kpoints_frac);
//###########################################################################


int main(int argc, char ** argv) { 
  if(argc < 2) {
    cout << "usage: siesta2qmc < -o basename > < -fold fold_dir nfold > [siesta stdout]\n";
    cout << "  -o          use basename as the prefix for QWalk.  Defaults to qwalk\n";
    exit(1);
  }
  
  int result = 0;

  string outputname="qwalk";
  int fold_dir=-1;
  int nfold=0;
  
  for(int i=1; i< argc-1; i++) { 
    if(!strcmp(argv[i],"-o") && i+1 < argc) { 
      outputname=argv[++i];
    }
    if(!strcmp(argv[i],"-fold") && i+2 < argc) { 
      fold_dir=atoi(argv[++i]);
      nfold=atoi(argv[++i]);
    }
  }
  string infilename=argv[argc-1];
  
  //------------------------------read in the variables from the output file
  vector <Atom> atoms;
  vector <vector < double> > latvec;
  vector <vector < double > > recip_latvec;
  coord_t lattice_const_type;
  double lattice_constant;

  Slat_wf_writer slwriter;
  
  //Default Magnification
  slwriter.magnification=50.0;
 

  //Assume complex, if imaginary part happens to be 0 later, copy to a real
  //set of orbitals and proceed
  vector < vector <vector <dcomplex> > > moCoeff;

  int spin_pol=0;
  vector <Spline_pseudo_writer> pseudo;
  
  
  ifstream is(infilename.c_str());
  if(!is) { 
    cout << "Couldn't open " << infilename << endl;
    exit(1);
  }
  string line;
  string space=" ";
  vector <string> currline;
  int netCharge=0;
  vector < vector < double > > kpoints;
  vector < vector < double > > kpoints_frac;
  coord_t coord_type;
  while(getline(is,line)) { 
    currline.clear();
    split(line, space, currline);
    read_atoms(is,currline,atoms, coord_type);
    read_lattice_vector(is,currline,latvec);
    read_mo_coefficients(is, currline, slwriter, moCoeff, kpoints);

    if(currline.size()> 6 && currline[2]=="spin" && currline[3]=="polarization") { 
      spin_pol=atoi(currline[6].c_str());
      double tmp=atof(currline[6].c_str());
      if(fabs(spin_pol-tmp) > 0.5) spin_pol+=1;
      assert(fabs(spin_pol-tmp) < 0.1);
      cout << "spin polarization " << spin_pol << endl;
    }
    if(currline.size() > 8 && currline[0]=="redata:" && currline[1]=="Net" && currline[2]=="charge") {
      cout << "getting charge on " << line << endl;
      netCharge=atoi(currline[7].c_str());
      double tmp=atof(currline[7].c_str());
      if(fabs(netCharge-tmp) > 0.5) spin_pol+=1;
      assert(fabs(netCharge-tmp) < 0.1);
      cout << "net charge " << netCharge << endl;
    }
  }
  is.close();

  if(slwriter.calctype=="") { 
    cout << "Couldn't find wave function coefficients.  You may want to make sure that \n"
    << "LongOutput and WaveFuncKPoints are set.\n";
    exit(1);
  }
  
  int num_kpoints = kpoints.size();
  reciprocal_lattice_vectors(latvec, recip_latvec);
  kpoint_to_frac_coords(recip_latvec, kpoints, kpoints_frac);

  cout << " coord_type " << coord_type << endl;
  if(coord_type==fractional) { 
    cout << "rescaling " << endl;
    vector <Atom> tmpatoms=atoms;
    int natoms=atoms.size();
    for(int at=0; at < natoms; at++) { 
      for(int i=0; i< 3; i++) atoms[at].pos[i]=0.0;
      for(int i=0; i< 3; i++) { 
        for(int j=0; j< 3; j++) { 
          atoms[at].pos[i]+=tmpatoms[at].pos[j]*latvec[j][i];
        }
      }
    }
  } else if (coord_type==scaled) {
    cout << "BROKEN" << endl;
    for (int n=0; n<atoms.size(); ++n)
      for (int i=0; i<3; ++i)
	atoms[n].pos[i] *= lattice_constant;
  }

 
  vector <Spline_basis_writer> basis;
  read_basis(atoms, basis);
  read_psp(atoms, pseudo);
  
  for (int i = 0; i < num_kpoints; ++i)
    {
      fix_basis_norm(atoms, basis, moCoeff[i]);
      for(int j=0; j< nfold; j++) 
	fold_kpoint(slwriter, latvec,fold_dir,moCoeff[i],atoms);
    }
  
  int nelectrons=0;
  for(vector<Atom>::iterator at=atoms.begin(); 
      at != atoms.end(); at++) { 
    nelectrons+=at->charge;
  }
  nelectrons-=netCharge;
  cout << "nelectrons " << nelectrons << endl;
  assert((nelectrons-spin_pol)%2==0);
  slwriter.nup=int(nelectrons-spin_pol)/2+spin_pol;
  slwriter.ndown=int(nelectrons-spin_pol)/2;




  Shifter shiftobj;
  vector <double> origin(3);
  for(int i=0; i< 3; i++) origin[i]=0;
  shiftobj.origin=origin;

  //Correct phase shift for atoms not at the origin in complex wfns
  //Easiest to use k-points in siesta format which are reciprocal of
  //normal bohr coordinates
  int f = 0; //index for the correct basis functions to shift
  
  for(unsigned int at=0; at< atoms.size(); at++) {
    int bas=atoms[at].basis;
    int nfunc=basis[bas].nfunc();
    int temp_f = f;
    for (int n = 0; n < num_kpoints; ++n) {
      f = temp_f;
      dcomplex kdotr=0;
      for(int d=0; d< 3; d++) 
	kdotr += kpoints[n][d] * (atoms[at].pos[d] - origin[d]);
      kdotr = exp(kdotr*dcomplex(0.0,1.0));
      for(int i=0; i< nfunc; i++) {
	for(int mo=0; mo < moCoeff[n].size(); mo++) {
	  moCoeff[n][mo][f] *= kdotr;
	}
	f++;
      }
    }
  }
  
  //Make sure all the atoms are inside the simulation cell
  //Taking care of phase factors for general k-points
  //This is problem if left to qwalk for complex k-points
  
  f = 0;
  
  for(unsigned int at=0; at< atoms.size(); at++) {
    vector <int> nshift;
    int bas=atoms[at].basis;
    int nfunc=basis[bas].nfunc();
    if(shiftobj.enforcepbc(atoms[at].pos, latvec, nshift)) {
      cout << "at " << at << "  shifted " << nshift[0] << "  " << nshift[1] 
	   << "   " << nshift[2] << endl;
      int temp_f = f;
      for (int n = 0; n < num_kpoints; ++n) {
	f = temp_f;
	dcomplex kdots=0;
	for(int d=0; d< 3; d++) 
	  kdots+=kpoints_frac[n][d]*nshift[d];
	//kdots=cos(pi*kdots);
	kdots=exp(pi*kdots*dcomplex(0.0,1.0));
	//cout << "kdots " << kdots << endl;
	for(int i=0; i< nfunc; i++) {
	  for(int mo=0; mo < moCoeff[n].size(); mo++) {
	    moCoeff[n][mo][f]*=kdots;
	  }
	  f++;
	}
      }
    }
    else f+=nfunc;
  }
  
  
  
  
  //---------------------------------------------
  //--Write out all the collected data
  //---------------------------------------------
  

  //---------------------Basis file
  slwriter.basisname=outputname+".basis";
  ofstream os (slwriter.basisname.c_str());
  os.precision(15);
  for(vector<Spline_basis_writer>::iterator i=basis.begin();
      i!= basis.end(); i++) {
    os << "BASIS { \n";
    i->print_basis(os);
    os << "}\n\n";
  }
  os.close(); os.clear();


  //--------------------------Jastrow 2 output
  string jast2outname=outputname+".jast2";
  double basis_cutoff=find_basis_cutoff(latvec);
  Jastrow2_wf_writer jast2writer;
  jast2writer.set_atoms(atoms);
  
  
  ofstream jast2out(jast2outname.c_str());
  print_std_jastrow2(jast2writer, jast2out, basis_cutoff);
  jast2out.close();




  //Loop over all kpoints and create appropriate input files for each
  for (int n = 0; n < num_kpoints; ++n)
    {
      stringstream buffer;
      vector < vector < double > > realMo;
      buffer << n;
      string tag("_" + buffer.str());
      if (num_kpoints == 1) tag = "";

      slwriter.orbname=outputname + tag + ".orb";
      string slateroutname=outputname + tag +".slater";
      string sysoutname=outputname + tag +".sys";

      //Check if the kpoint is complex
      int is_complex = 0;
      for (int i = 0; i < moCoeff[n].size(); ++i)
	{
	  for (int j = 0; j < moCoeff[n][i].size(); ++j)
	    {
	      if (fabs(moCoeff[n][i][j].imag()) > 0)
		{
		  is_complex = 1;
		  break;
		}
	    }
	}

      //Set the labels, and if it is real, copy to a real equivalent
      if (is_complex)
	{
	    slwriter.orbtype = "CORBITALS";
	    slwriter.mo_matrix_type = "CCUTOFF_MO";
	}
      else
	{
	  slwriter.orbtype = "ORBITALS";
	  slwriter.mo_matrix_type = "CUTOFF_MO";

	  //Copy real part of orbital, to real orbitals
	  realMo.resize(moCoeff[n].size());
	  for (int i = 0; i < moCoeff[n].size(); ++i)
	    {
	      realMo[i].resize(moCoeff[n][i].size());
	      for (int j = 0; j < moCoeff[n][i].size(); ++j)
		realMo[i][j] = moCoeff[n][i][j].real();
	    }
			       
	}


  
      //-----------------ORB file
      
      os.open(slwriter.orbname.c_str());
      //this is really redundant with gamess2qmc.cpp..should probably be refactored 
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
      
      if (is_complex)
	print_orbitals(os, centers, nbasis, moCoeff[n]);
      else
	print_orbitals(os, centers, nbasis, realMo);
      
      os.close(); os.clear();
      
      //--------------------Slater determinant file
      
      
      os.open(slateroutname.c_str());
      slwriter.print_wavefunction(os);
      os.close(); os.clear();
      
      
      //---------------------sys file
      //May want to add a command-line switch to print out a molecule
      
      ofstream sysout(sysoutname.c_str());
      sysout << "SYSTEM { PERIODIC \n";
      sysout << "  NSPIN { " << slwriter.nup << "  "
	     << slwriter.ndown << " } \n";
      for(int at=0; at <natoms; at++) {
	atoms[at].print_atom(sysout);
      }
      
      sysout << "LATTICEVEC { \n";
      double min_latsize=1e8;
      for(int i=0; i< 3; i++) { 
	double length=0;
	for(int j=0; j< 3; j++) {
	  sysout << latvec[i][j] << " ";
	  length+=latvec[i][j]*latvec[i][j];
	}
	sysout << endl;
	if(min_latsize > length) min_latsize=length;
      }
      sysout << "  }\n";
      //Here we assume that nothing has a larger cutoff than 12.0 bohr..is that a good guess?
      sysout << "  cutoff_divider " << sqrt(min_latsize)/12.0 << endl;
      
      //Print out the k-point to the sys file
      sysout << "  kpoint { " << kpoints_frac[n][0] 
	     << "   " << kpoints_frac[n][1] 
	     << "   " << kpoints_frac[n][2] << " } " << endl;
      
      sysout << "}\n\n\n";
      
      int npsp=pseudo.size();
      for(int psp=0; psp < npsp; psp++) {
	pseudo[psp].print_pseudo(sysout);
      }
      sysout.close();
      
    }
}
//###########################################################################

void read_atoms(istream & is, vector <string> & currline, 
                vector <Atom> & atoms, coord_t & coord_type) {
  string space=" ";
  if(currline.size() > 3 && currline[0]=="outcoor:" && currline[2]=="coordinates") { 
    atoms.clear(); //so we always take the last set of atoms.
    double fac=1;
    if(currline[3]=="(Bohr):") {
      fac=1;
      coord_type=bohr;
    }
    else if(currline[3]=="(Ang):") { 
      cout << "converting from angstrom" << endl;
      fac=1/.529177249;
      coord_type=ang;
    }
    else if(currline[3]=="(scaled):") { 
      cout << "scaled coordinates" << endl;
      fac=1.0;
      coord_type=scaled;
    }
    else if(currline[3]=="(fractional):") { 
      cout << "fractional coordinates: WARNING: if you have non-orthogonal unit cells, the conversion may be incorrect.  Check atomic coordinates!!" << endl;
      fac=1.0;
      coord_type=fractional;
    }
    else { 
      cout << "Don't understand format of the coordinates: " << currline[3] << endl;
      exit(1);
    }
    //cout << "reading in coordinates" << endl;
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
      //tmpatom.print_atom(cout);
      getline(is,line);
      currline.clear(); split(line,space,currline);
      
    }
      
  }
}

//##########################################################################
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
    //cout << "reading in lattice vectors " << endl;
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
  }
}
//##########################################################################
//Need reciprocal lattice vectors to properly normalize k-points under certain
//input conditions
void reciprocal_lattice_vectors(vector < vector <double > > &latvec,
				vector < vector <double > > &recip_latvec)
{
  //Allocate vectors
  vector <double> temp_vec;
  temp_vec.resize(3);

  recip_latvec.resize(3);
  for (int i = 0; i < 3; ++i)
    recip_latvec[i].resize(3);
  
  //Compute vectors
  recip_latvec[0] = cross(latvec[1], latvec[2]);
  temp_vec = cross(latvec[1], latvec[2]);
  recip_latvec[0] /= dot(latvec[0], temp_vec);

  recip_latvec[1] = cross(latvec[2], latvec[0]);
  temp_vec = cross(latvec[2], latvec[0]);
  recip_latvec[1] /= dot(latvec[1],temp_vec);

  recip_latvec[2] = cross(latvec[0], latvec[1]);
  temp_vec = cross(latvec[0], latvec[1]);
  recip_latvec[2] /= dot(latvec[2],temp_vec);

  for (int i = 0; i < 3; ++i)
    recip_latvec[i] *= 2.0 * pi;
}
//#################################################################
//Compute to qwalk-friendly fractional coordinates in the reciprocal
//lattice basis
void kpoint_to_frac_coords(vector < vector <double> > &recip_latvec,
			   vector < vector <double> > &kpoints,
			   vector < vector <double> > &kpoints_frac)
{
  //K-points are stored in the standard basis in reciprocal lattice space
  //so to convert, simply need to invert the matrix of reciprocal lattice
  //vectors and apply the resulting matrix to the k-point vectors...

  vector <double> temp_vec;
  temp_vec.resize(3);

  kpoints_frac.resize(kpoints.size());
  for (int n = 0; n < kpoints.size(); ++n)
    kpoints_frac[n].resize(3);

  vector < vector <double> > change_matrix;
  matrix_inverse(recip_latvec, change_matrix);
  

  for (int n = 0; n < kpoints.size(); ++n)
    {
      temp_vec[0] = change_matrix[0][0] * kpoints[n][0] +
	change_matrix[0][1] * kpoints[n][1] +
	change_matrix[0][2] * kpoints[n][2];
      temp_vec[1] = change_matrix[1][0] * kpoints[n][0] +
	change_matrix[1][1] * kpoints[n][1] +
	change_matrix[1][2] * kpoints[n][2];
      temp_vec[2] = change_matrix[2][0] * kpoints[n][0] +
	change_matrix[2][1] * kpoints[n][1] +
	change_matrix[2][2] * kpoints[n][2];


      //This factor of two is a little bit mysterious to me,
      //but it matches the qwalk convention and everything works out...
      for (int i = 0; i < 3; ++i)
	kpoints_frac[n][i] = 2.0 * temp_vec[i];
    }
     
}

//###########################################################################
//It turns out, that when there are imaginary coeff's that are both double,
//digit and negative, that the real and imaginary coeff's are not bridged by
//a space, but rather by the - sign.  Thus the more general split() function
//needs a modification
void orbital_split(string line, vector <string> & currline)
{
  string space = " ";
  split(line, space, currline);
  if (currline.size() < 7) //If last two were joined by a -...
    {
      int loc = currline[5].rfind("-", currline[5].size());
      string temp_1 = currline[5].substr(0, loc - 1);
      string temp_2 = currline[5].substr(loc, currline[5].size());
      currline[5] = temp_1;
      currline.push_back(temp_2);
    }
}


//###########################################################################
void read_mo_coefficients(istream & is, vector <string> & currline,
                          Slat_wf_writer & slwriter, vector <vector <vector <dcomplex> > > & moCoeff,
			  vector < vector <double> > &kpoints) { 
  slwriter.use_global_centers=true;
  slwriter.write_centers=false;
  string space=" ";
  int num_kpoints = 0;

  //Temp value for storing kpoints
  vector <double> temp_kpoint;
  temp_kpoint.resize(3);

  if(currline.size() > 3 && currline[0]=="writewave:" && currline[2]=="Functions") { 
    string line;
    getline(is, line); //empty line;

    //Read number of kpoints
    getline(is, line); //number of k-points
    currline.clear(); split(line, space, currline);
    num_kpoints = atof(currline[4].c_str());
    moCoeff.resize(num_kpoints);
    kpoints.resize(num_kpoints);
    cout << "Detected " << num_kpoints << " k-points in Siesta output file." << endl;

    getline(is, line); //nspin
    currline.clear(); split(line, space, currline);
    //cout << line << endl;
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
    
    for (int i=0; i < num_kpoints; ++i)
      {
	for(int spin=0; spin < nspin; spin++) { 
	  getline(is, line); getline(is,line); //one blank line and *'s
	  
	  //Load k-point value into vector for printing in sys file
	  getline(is,line); //k-point line
	  currline.clear(); split(line, space, currline);
	  for (int j = 0; j < 3; ++j)
	    temp_kpoint[j] = atof(currline[j+3].c_str());
	  kpoints[i] = temp_kpoint;
	  
	  getline(is,line); //spin component
	  getline(is,line); //Num. wavefunctions
	  currline.clear(); split(line, space, currline);
	  int norb=atoi(currline[3].c_str());
	  //cout << "norb " << norb << endl;
	  //Now start with reading the ****orbitals*** (KS wave functions)
	  for(int orb_num=0; orb_num < norb; orb_num++) { 
	    for(int j=0; j< 5; j++) getline(is,line); //blank, wf#, energy, ---'s, header
	    vector <dcomplex> orb(nbasis);
	    for(int b=0; b< nbasis; b++) { 
	      getline(is, line);
	      currline.clear(); 
	      orbital_split(line, currline);
	      orb[b] = complex<double>(atof(currline[5].c_str()), atof(currline[6].c_str()));
	    }
	    moCoeff[i].push_back(orb);
	    getline(is, line); //a line of --'s
	  }
	  if(spin==0 && nspin==2) { 
	    slwriter.spin_dwn_start=moCoeff[i].size();
	  }
	} //spin loop
	//cout << "total norbs " << moCoeff.size() << endl;
      }
  }
  
  
  //does not fill nup, ndown, mo_matrix_type, orbname, or basisname in slwriter
}



//###########################################################################

void smooth_grid(vector <double> & r,
                 vector <double> & vals) { 
  vector <double> newvals;
  newvals.resize(vals.size());
  int n=vals.size();
  assert(r.size()==vals.size());
  for(int i=2; i< n-2; i++) { 
    newvals[i]=(vals[i-2]+2*vals[i-1]+3*vals[i]+2*vals[i+1]+vals[i+2])/9.0;
  }
  newvals[1]=(vals[0]+vals[1]+vals[2])/3.0;
  newvals[0]=vals[0];
  newvals[n-2]=(vals[n-3]+vals[n-2]+vals[n-1])/3.0;
  newvals[n-1]=vals[n-1];
  vals=newvals;
}


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
      //cout << "rad " << rad << " x " << x << " r1 " << r[interval] 
      //  << " r2 " << r[interval+1] << endl;
      val=(1-x)*vals[interval]+x*vals[interval+1];
    }
    ur.push_back(rad);
    uvals.push_back(val);
  }
  
}


//###########################################################################

void pad_spline(vector <double> & r, vector <double> & vals) { 
  double spacing=r[1]-r[0];
  double extent=r[r.size()-1];
  while(extent < 14.0) {  
    extent+=spacing;
    r.push_back(extent);
    vals.push_back(0.0);
  }
}


//Here we smoothly cut off the function instead of just chopping it off, 
//as Siesta seems happy to do.
void smooth_cutoff(vector <double> & r, 
                  vector <double> & vals,
                  double smooth=1.2) { 
  int n=r.size();
  double threshold=1e-8;
  double cut=r[n-1];
  for(int i=n-1; i> 0; i--) { 
    if(fabs(vals[i]) > threshold) { 
      cut=r[i];
      cout << "Found cutoff radius of " << cut << " at " << i <<  endl;
      break;
    }
  }
  
  double cutmin=cut-smooth;
  vector <double> newvals;
  newvals.resize(n);
  double a,b,c,d; //extrapolation parameters
  for(int j=0; j< n; j++) {
    if(r[j] > cutmin) {
      if(r[j] > cut) { //if we're beyond the cutoff completely
        newvals[j]=0;
      }
      else {  //if we're in the smooth cutoff region
        double zz=(r[j]-cutmin)/smooth;
        newvals[j] = vals[j]*(1-zz*zz*zz*(6*zz*zz-15*zz+10));
      }
      
        //better, but makes a small jump right at the cutoff
      //if(j+1 < n)
      //  newvals[j]=0.25*vals[j-1]+0.5*vals[j]+0.25*vals[j+1];
      if(j+2 < n) 
        newvals[j]=(vals[j-2]+2*vals[j-1]+3*vals[j]+2*vals[j+1]+vals[j+2])/9.0;
      else newvals[j]=0.0;
      
    }
    else newvals[j]=vals[j];
    
  }
/*
if(r[j] >= cutmin && r[j-1] < cutmin) { 
  double v, der, der2;
  v=vals[j];
  double diff=r[j+1]-r[j-1];
  der=(vals[j+1]-vals[j-1])/diff;
  der2=(vals[j+1]+vals[j-1])/(diff*diff);
  double cutmin3=cutmin*cutmin*cutmin;
  double cut3=cut*cut*cut;
  smooth=-smooth;
  cout << "der " << der << " der2 " << der2 << endl;
  a=(der2-2/smooth*(der-v/smooth))/(6.*cutmin-2/smooth*(3*cutmin*cutmin-(cutmin3-cut3)/smooth));
  b=1/smooth*(der-v/smooth-a*(3*cutmin*cutmin-(cutmin3-cut3)/smooth));
  c=(v-a*(cutmin3-cut3))/smooth-b*(cutmin+cut);
  d=-a*cut3-b*cut*cut-c*cut;
  cout << "a,b,c,d " << a << " " << b << " " << c <<"  "<< d << endl;
  newvals[j]=vals[j];
  cout <<"val " << vals[j] << " function " << a*r[j]*r[j]*r[j]+b*r[j]*r[j]+c*r[j]+d << endl;
}
else if(r[j] > cutmin && r[j] < cut) { 
  
  newvals[j]=a*r[j]*r[j]*r[j]+b*r[j]*r[j]+c*r[j]+d;
}
else if(r[j] > cut) newvals[j]=0.0;
    else newvals[j]=vals[j];
  }
  vals=newvals;
*/
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
    for(int st=1; st < 20; st++) { //search over possible states
      string prefix="ORB.S";
      append_number(prefix, st);
      prefix+=".";
      string postfix="."+*nm;
      int n=1;
      while(true) { 
        string file=prefix;
        append_number(file,n);
        file+=postfix;
        ifstream is(file.c_str());
        if(!is) break;
        string dummy;
        is >> dummy >> dummy; 
        if(dummy != *nm) { 
          cout << "file " << file << " doesn't match name in file " << dummy << endl;
          exit(1);
        }
        int el;
        is >> el;
        is.ignore(180,'\n');
        cout << "found " << file << " l-value " << el <<  endl;
        switch(el) { 
          case 0:
            tmp_basis.types.push_back("S");
            break;
          case 1:
            tmp_basis.types.push_back("P_siesta");
            break;
          case 2:
            tmp_basis.types.push_back("5D_siesta");
            break;
          case 3:
            tmp_basis.types.push_back("7F_siesta");
            break;
          default:
            cout << "Don't support this l-value.  Bug the maintainer." << endl;
            exit(1);
        }
        //Read in the file..
        is.ignore(180,'\n');
        double rad, val;
        vector <double> rads, vals;
        while(is >> rad && is >> val) {
          rads.push_back(rad); vals.push_back(val);
        }
        vector <double> urad, uval;
        smooth_grid(rads,vals);
        pad_spline(rads,vals);
        urad=rads; uval=vals;
        //make_uniform(rads,vals,urad,uval);
        //smooth_cutoff(urad,uval);
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
  //cout << "****Pseudopotential " << endl;
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
    //cout << *at << " : effective charge : " << spl[5] <<  " integerized " << atoi(spl[5].c_str()) << endl;
    int zeff=atoi(spl[5].c_str());
    for(vector<Atom>::iterator i=atoms.begin(); i!= atoms.end();
        i++) { 
      if(i->name==*at) i->charge=zeff;
    }
    int npoints=atoi(spl[2].c_str());
    int nl=atoi(spl[0].c_str());
    int nl_up=atoi(spl[1].c_str());
    int nl_down=nl;
    if(nl_up > 0) { 
      tmp_pseudo.spin_dep=true;
      cout << "Spin-dependent pseudopotential found for " << *at << endl;
      nl+=nl_up;
      //cout << "Can't deal with spin-polarized pseudopotentials..sorry\n";
      //exit(1);
    }
    //cout << "npoints " << npoints <<  " nl " << nl << endl;
    is.ignore(180,'\n'); //Radial grid..
    vector <double>  rad;
    for(int i=0; i< npoints; i++) { 
      double dum;
      is >> dum;
      rad.push_back(dum);
    }

    vector <int> lvalue;
    vector <vector <double> > up_rad(nl_down), up_val(nl_down);
    vector <vector <double> > down_rad(nl_down), down_val(nl_down);
    for(int l=0; l< nl; l++) { 
      vector <double> val;
      is.ignore(180,'\n');
      getline(is, line);
      //cout << line << endl;
      int currl; is >> currl;
      cout << "l " << l << " "  << currl << endl;
      lvalue.push_back(currl);
     // assert(l==currl || l==currl-nl_down);
      for(int i=0; i< npoints; i++) { 
        double dum; is >> dum; 
        val.push_back(dum);
        val[i]*= 0.5/rad[i];
      }
      vector <double> urad, uval;
      make_uniform(rad,val,urad,uval);
      
      if(l < nl_down) { 
        down_rad[currl]=urad;
        down_val[currl]=uval;
      }
      else { 
        up_rad[currl]=urad;
        up_val[currl]=uval;
      }
      //tmp_pseudo.psp_pos.push_back(urad);
      //tmp_pseudo.psp_val.push_back(uval);
    }
    tmp_pseudo.label=*at;
    int npts=down_val[0].size();

    for(int i=0; i < nl_down; i++) { 
      assert(down_rad[i].size() > 0);
      cout << i << " : " << up_rad[i].size() << "  " << down_rad[i].size() << endl;
      if(up_rad[i].size()==0) {
        cout << "setting up channel " << i << " to down channel "<< endl;
        up_rad[i]=down_rad[i];
        up_val[i]=down_val[i];
      }
      else { 
        for(int j=0; j< npts; j++) { 
          assert(fabs(up_rad[i][j]-down_rad[i][j]) < 1e-9);
          up_val[i][j]+=down_val[i][j];
        }
      }
    }

    
    for(int l=0; l < nl_down-1; l++) { 
      for(int i=0; i< npts; i++) { 
        up_val[l][i]-=up_val[nl_down-1][i];
        down_val[l][i]-=down_val[nl_down-1][i];
      }
    }
     
    
    for(int l=0; l < nl_down; l++) { 
      tmp_pseudo.psp_pos.push_back(down_rad[l]);
      tmp_pseudo.psp_val.push_back(down_val[l]);
    }
    if(tmp_pseudo.spin_dep) { 
      for(int l=0; l < nl_down; l++) { 
        tmp_pseudo.psp_pos.push_back(up_rad[l]);
        tmp_pseudo.psp_val.push_back(up_val[l]);
      }
    }
    
    /*
    int npts=tmp_pseudo.psp_val[0].size();
        
    if(nl_up > 0) { 
      //Add the local part onto the up spin channel, too..
      tmp_pseudo.psp_pos.push_back(tmp_pseudo.psp_pos[nl_down-1]);
      tmp_pseudo.psp_val.push_back(tmp_pseudo.psp_val[nl_down-1]);
    }
        
    for(int i=0; i< npts; i++) { 
      for(int l=0; l< nl_down-1; l++) { 
        tmp_pseudo.psp_val[l][i]-=tmp_pseudo.psp_val[nl_down-1][i];
      }
      for(int l=nl_down; l< nl-1; l++) { 
        tmp_pseudo.psp_val[l][i]-=tmp_pseudo.psp_val[nl-1][i];
      }
    }
         */
    pseudo.push_back(tmp_pseudo);
  }
}


//###########################################################################
void fix_basis_norm(vector <Atom> & atoms, 
                    vector <Spline_basis_writer> & basis,
                    vector < vector <dcomplex> > & moCoeff) {
  
  int nmo=moCoeff.size();
  int funcnum=0;
  
  //Looping over each symmetry in the basis attached to each atom
  for(vector<Atom>::iterator at=atoms.begin();
      at != atoms.end(); at++) { 
    vector<Spline_basis_writer>::iterator bas=basis.begin();
    //Assuming that there is just one basis per atom,
    //which should be correct..
    while(bas->label != at->name) bas++;
    
    for(vector<string>::iterator type=bas->types.begin();
        type!=bas->types.end(); type++) { 
      int L=0;
      if(*type=="S") L=0;
      else if(*type=="P_siesta") L=1;
      else if(*type=="5D_siesta") L=2;
      else if(*type=="7F_siesta") L=3;
      else { cout << "unknown type! " << *type << endl; exit(3); }
      double norm=(2.0*L+1)/(4*pi);
      for(int m=-L; m< L+1; m++) { 
        int absm=abs(m);
        double i=1;
        for(int j=L-absm+1; j < L+absm+1; j++) {
          i*=j;
        }
        double mnorm=norm/i;
        if(m!=0) mnorm*=2.;
        mnorm=sqrt(mnorm);
        if(L==1) { 
          switch(m) { 
            case -1:
              mnorm*=-1; break;
            case 0:
              mnorm*=1; break;
            case 1:
              mnorm*=-1; break;
          }
        }
        else if(L==2) { 
          switch(m) { 
            case -2:
              mnorm*=6;
              break;
            case -1:
              mnorm*=-3; break;
            case 0:
              mnorm*=0.5;break;
            case 1:
              mnorm*=-3; break;
            case 2:
              mnorm*=3; break;
            default:
              cout << "error in mnorm assignment" << endl; exit(3);
          }
        }
        else if(L==3) { 
          switch(m) { 
            case -2:
              mnorm*=2;
              break;
          }
        }
        
              
                
        //cout << "i " << i << "  mnorm " << mnorm << endl;

        for(int mo=0; mo < nmo; mo++) { 
          //cout << moCoeff[mo].size() << " " << funcnum << endl;
          if(moCoeff[mo].size() <= funcnum) { 
            cout << "You have less MO coefficients than atomic orbitals in the ORB files.  This could be from a remnant\n"
            "of a previous run.  The easiest solution is to rm ORB* and rerun the siesta calculation or try to guess\n"
            "which ORB files are extraneous\n";
          }
          assert(moCoeff[mo].size() > funcnum);
          moCoeff[mo][funcnum]*=mnorm;
        }
        funcnum++;
      }
    }
    
  }
  
  
}

