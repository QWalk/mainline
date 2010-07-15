/*
 
Copyright (C) 2009 Michal Bajdich

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
#include "elements.h"  
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include "vecmath.h"
#include "elements.h"
#include "CBasis_function.h"
#ifdef USE_MPI
#include "mpi.h"
#endif

using namespace std;

class mpi_info_struct
{
public:
  mpi_info_struct()
  {
    node=0;
    nprocs=1;
  }
  int node;
  int nprocs;
};
extern mpi_info_struct mpi_info;
mpi_info_struct mpi_info;


int roundoff(double x){
  double y;
  y=floor(x);
  if (x>=y+0.5) return int(y+1);
  else return int(y);
}

string givesymmetry(int & i){
  if(i==0)
    return "S";
  else if(i==1)
    return "P";
  else if(i==2)
    return "D";
  else if(i==3)
    return "F";
  else if(i==4)
    return "G";
  else
    return "what was that?";
}

double linear_extrapolate(vector <double> & x, vector <double> & y, double xpos){
  int i=0;
  while(i<x.size()){
    if(x[i]>xpos)
      break;
    i++;
  }
  //return y[i-1];
  if(i!=0)
    return ((y[i]-y[i-1])/(x[i]-x[i-1]))*(xpos-x[i-1])+y[i-1];
  else
    return y[0];
}


double Determinant(vector < vector <double> > & A){
  if(A.size()!=3){
    cout <<"ERROR in Determinant"<<endl;
    exit(1);
  }
  return A[0][0]*A[1][1]*A[2][2]+A[0][1]*A[1][0]*A[2][0]+A[0][2]*A[1][0]*A[2][1]
        -A[0][2]*A[1][1]*A[2][0]-A[0][1]*A[1][0]*A[2][2]-A[0][0]*A[1][2]*A[2][1];

}


void read_abinit_sys(std::istream & is, vector < Atom > & atoms, vector < Atom > & primatoms,
		     vector <Spline_pseudo_writer> & abinit_pseudo,
		     vector <double> & origin,
		     vector < vector < double> > & latvec,
		     vector < vector < double> > & primlatvec,
		     vector <int> & factor,
		     int & nelectrons, int & spin, double & eref, int & debug);

void read_abinit_wf(string & filename,  CBasis_function * abinitbasis, int & debug, int & norbs, int & node);

void write_abinit_basis(ostream & os, CBasis_function * abinitbasis, int & norbs, vector < vector <double > > & primlatvec, string & outname) {
  string indent="  ";
  if(mpi_info.node==0)
    abinitbasis->writeinput(indent, os, norbs, primlatvec, outname);
}

void plot_orbitals(string & outputname,  CBasis_function * abinitbasis, int & norbs, vector < Atom > & atoms,
                   vector < vector < double> > & primlatvec, vector < double> & origin, double & resolution, int & i);

void usage(const char * name) {
  if(mpi_info.node==0){
    cout << "usage: " << name <<   " <options> <pwfn.data file> " << endl;
    cout << "Where options can be: \n";
    cout << "-o <string> Base name for your run case\n";
    cout << "-debug to have more informative printout\n";
    cout << "-res <double default=0.1> resolution for plots \n";
    cout << "-norbs <int> number of orbitals \n";
    cout << "-spin <int default=1> 1=singlet, 2=doublet...\n";
    cout << "-factor <int i, default=1> <int j, default=1> <int k, default=1> simulation cell is ixjxk of primitive cell\n";
    cout << "-kpoint <double> <double> <double> for the simulation cell (in qwalk units of cell/Pi), e.g.: for L-point: 1 1 1 \n";
  }
  exit(1);
}


int main(int argc, char ** argv) {
#ifdef USE_MPI
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &(mpi_info.nprocs));
  MPI_Comm_rank(MPI_COMM_WORLD, &(mpi_info.node));
  cout << "processor " << mpi_info.node << " alive " << endl;
#endif


  int norbs=-1;
  int debug=0;
  int spin=1;
  vector <int> factor(3);
  vector <double> kpoint;
  factor[0]=factor[1]=factor[2]=1;
  double resolution=0.1;
  string infilename, outputname;
  
  if(argc >= 2) {
    infilename=argv[argc-1];
  }
  else { usage(argv[0]); }

  for(int i=1; i< argc-1; i++) {
    if(!strcmp(argv[i], "-o") && argc > i+1) {
      outputname=argv[++i];
    }
    else if(!strcmp(argv[i], "-debug")) {
     debug=1;
    }
    else if(!strcmp(argv[i], "-res") && argc > i+1) {
      resolution=atof(argv[++i]);
    }
    else if(!strcmp(argv[i], "-norbs") && argc > i+1) {
      norbs=atoi(argv[++i]);
    }
    else if(!strcmp(argv[i], "-spin") && argc > i+1) {
      spin=atoi(argv[++i]);
    }
    else if(!strcmp(argv[i], "-factor") && argc > i+3) {
      for(int d=0;d<3;d++)
	factor[d]=atoi(argv[++i]);
    }
    else if(!strcmp(argv[i], "-kpoint") && argc > i+3) {
      for(int d=0;d<3;d++)
	kpoint.push_back(atof(argv[++i]));
    }
    else {
      usage(argv[0]);
    }
  }
 

  /*
  if(kpoint.size()!=3){
    if(mpi_info.node==0)
      cout << "Please supply the kpoint in the simulation cell by: -kpoint for the simulation cell (in units of cell/Pi), e.g.: for gamma: 1 1 1 \n" << endl;
    exit(1);
  }
  */

  if(outputname == "") {
    outputname=infilename;
  }

  //objects to read into
  //vector <Atom> atoms; 
  //Blochwave_function * abinitbasis; 
  //vector <Center> centers;
  //vector < vector < double> > moCoeff;
  //vector <Spline_pseudo_writer> pseudo;
  //vector <double> origin(3);
  //Slat_wf_writer slwriter;
  //Jastrow2_wf_writer jast2writer;

  ///////////////////////////////////////////////////////
  vector <Spline_pseudo_writer> pseudo;
  vector < vector <double> > latvec, primlatvec;
  vector <double> origin(3);
  vector <Atom> atoms;
  vector <Atom> primatoms;
  Slat_wf_writer slwriter;
  vector <double> occupation; //electronic occupation

  CBasis_function * abinitbasis;
  abinitbasis=new Blochwave_function;

  slwriter.mo_matrix_type="BSPLINE_MO";
  double eref=0;
  string sysfilename=infilename;
  string wffilename=infilename;
  ifstream sysfile(sysfilename.c_str());
  if(!sysfile) {
    cout << "couldn't open " << sysfilename << endl;
    exit(1);
  }

  int nelectrons=0;
  int nspins=0;
  //here read the system stuff
  read_abinit_sys(sysfile, atoms, primatoms, pseudo, origin, latvec, primlatvec, factor, nelectrons, nspins, eref, debug);
  sysfile.close();

  slwriter.magnification=1.0/sqrt((double)factor[0]*factor[1]*factor[2]);


  //------------reciprocal cell
  vector < vector <double> >  crossProduct(3);

  crossProduct[0].push_back(primlatvec[1][1]*primlatvec[2][2]-primlatvec[1][2]*primlatvec[2][1]);
  crossProduct[0].push_back(primlatvec[1][2]*primlatvec[2][0]-primlatvec[1][0]*primlatvec[2][2]);
  crossProduct[0].push_back(primlatvec[1][0]*primlatvec[2][1]-primlatvec[1][1]*primlatvec[2][0]);

  crossProduct[1].push_back(primlatvec[2][1]*primlatvec[0][2]-primlatvec[2][2]*primlatvec[0][1]);
  crossProduct[1].push_back(primlatvec[2][2]*primlatvec[0][0]-primlatvec[2][0]*primlatvec[0][2]);
  crossProduct[1].push_back(primlatvec[2][0]*primlatvec[0][1]-primlatvec[2][1]*primlatvec[0][0]);

  crossProduct[2].push_back(primlatvec[0][1]*primlatvec[1][2]-primlatvec[0][2]*primlatvec[1][1]);
  crossProduct[2].push_back(primlatvec[0][2]*primlatvec[1][0]-primlatvec[0][0]*primlatvec[1][2]);
  crossProduct[2].push_back(primlatvec[0][0]*primlatvec[1][1]-primlatvec[0][1]*primlatvec[1][0]);
  
  vector < vector <double> >  reciplatvec(3);
  double det=Determinant(primlatvec);
  if(mpi_info.node==0 && debug)
    cout<<"Primitive cell volume: "<<det<<endl;
 

  if(mpi_info.node==0 && debug)
    cout<<"Reciprocal Lattice Vectors to Primitive Lattice: "<<endl;
  for(int i=0; i< 3; i++) {
    for(int j=0; j< 3; j++) {
      reciplatvec[i].push_back(crossProduct[i][j]/det);
      if(mpi_info.node==0 && debug ){
	cout <<reciplatvec[i][j]<<" ";
      }
    }
    if(mpi_info.node==0 && debug)
      cout <<endl;
  }
  

  if(nspins==1) {
    if(nelectrons % 2 != 0) {
      cout << "It seems like you're doing RHF, but there is an odd number of"
      << " electrons!  I'm confused." << endl;
      exit(1);
    }
    slwriter.nup=nelectrons/2;
    slwriter.ndown=nelectrons/2;
    slwriter.calctype="RHF";
    if(norbs<nelectrons/2)
      norbs=nelectrons/2;
  }
  else if(nspins==2) {
    if(spin==1)
      cout << "Detected a ROHF or UHF calculation."
	"  Supply the the spin state with -spin <int> (where 1=singlet, 2=doublet...)" << endl;
   
    spin-=1;
    slwriter.calctype="UHF";
    slwriter.nup=(nelectrons-spin)/2 + spin;
    slwriter.ndown=(nelectrons-spin)/2;
    
    //cout <<slwriter.nup<<"  "<<slwriter.ndown<<endl;
    
    if(slwriter.nup+slwriter.ndown != nelectrons) {
      cout << "problem..  nup and ndown are calculated to be "
           << slwriter.nup << "   " << slwriter.ndown
           << " but they don't add up to be " << nelectrons << endl;
      exit(1);
    }

    if(norbs<nelectrons){
      norbs=nelectrons;
    }
    slwriter.spin_dwn_start=norbs-slwriter.ndown;
  }
  else {
    cout << "Failed to read: Spin polarized: .true./.false."<< endl;
    exit(1);
  }


  //reading abinit_wf
  read_abinit_wf(wffilename, abinitbasis, debug, norbs, mpi_info.node);

  //try to determine the K-point
  if(kpoint.size()!=3){
    vector < vector <double> > allkvectors;
    abinitbasis->getKvectors(primlatvec,allkvectors);
    //taking the first one
    if(allkvectors.size()==1)
      for(int d=0;d<3;d++)
	kpoint.push_back(allkvectors[0][d]);
    else{
      if(mpi_info.node==0)
	cout << "Please supply the kpoint in the simulation cell by: -kpoint for the simulation cell (in qwalk units of cell/Pi), e.g.: for L-point: 1 1 1 \n" << endl;
      exit(1);
    }
  }    

  //plot abinit_wf
  plot_orbitals(outputname, abinitbasis, norbs, primatoms, primlatvec, origin, resolution, debug);

  //write all the other files
  if(mpi_info.node==0){
    string orboutname=outputname+".orb";
    slwriter.orbname=orboutname;
    //one center at the origin only
    vector <Center> centers;
    Center tempcenter;
    tempcenter.name="origin";
    tempcenter.equiv_atom=0;
    tempcenter.basis=0;
    centers.push_back(tempcenter);
    
    
    bool status;
    slwriter.write_centers=true;
    string centeroutname;
    centeroutname=outputname+".centers";
    slwriter.centername=centeroutname;
    string basisoutname=outputname+".basis";
    slwriter.basisname=basisoutname;
    
    ofstream basisout(basisoutname.c_str());
    write_abinit_basis(basisout, abinitbasis, norbs, primlatvec, outputname);
    basisout.close();
    
    ofstream centerout(centeroutname.c_str());
    //centerout << centers.size() << endl;
    for(vector < Center>:: iterator i=centers.begin();
	i != centers.end(); i++) {
      i->print_center(centerout);
    }
    centerout.close();
  
    ofstream orbout(orboutname.c_str());
    int cycles=int(norbs/2)+norbs%2;
    for(int orbs=0;orbs<cycles;orbs++)
      orbout<<orbs+1<<"  1"<<endl;
    orbout.close();
  
    string slateroutname=outputname+".slater";
    
    ofstream slaterout(slateroutname.c_str());
    slwriter.print_wavefunction(slaterout);
    slaterout.close();
    
    
    //--------------------------Jastrow 2 output
    
    
    double basis_cutoff=find_basis_cutoff(latvec);
    string jast2outname=outputname+".jast2";
    Jastrow2_wf_writer jast2writer;
    jast2writer.set_atoms(atoms);
    
    
    ofstream jast2out(jast2outname.c_str());
    print_std_jastrow2(jast2writer, jast2out, basis_cutoff);
    jast2out.close();
    
    
    //-----------------------------------System output
    
    string sysoutname=outputname+".sys";
    ofstream sysout(sysoutname.c_str());
    sysout << "SYSTEM { PERIODIC \n";
    sysout << "  NSPIN { " << slwriter.nup << "  "
	   << slwriter.ndown << " } \n";
    
    sysout << "  CUTOFF_DIVIDER  5.00001" << endl;
    
    string intend="  ";

    sysout << "LATTICEVEC { \n";
    for(int i=0; i< 3; i++) {
      sysout <<intend;
      for(int j=0; j< 3; j++) {
	sysout << latvec[i][j] << "   ";
      }
      sysout << endl;
    }
    sysout << " }\n\n";

    

    sysout << "PRIMLATTICEVEC { \n";
    for(int i=0; i< 3; i++) {
      sysout <<intend;
      for(int j=0; j< 3; j++) {
	sysout << primlatvec[i][j] << "   ";
      }
      sysout << endl;
    }
    sysout << " }\n\n";
    
    sysout << "  ORIGIN { ";
    for(int i=0; i< 3; i++) sysout << origin[i] << "  ";
    sysout << "}\n";

    sysout << "  KPOINT { ";
    for(int i=0; i< 3; i++) sysout << kpoint[i] << "  ";
    sysout << "}\n";
    
    for(vector <Atom>::iterator at=atoms.begin(); at != atoms.end(); at++) {
      at->print_atom(sysout);
    }
    sysout << "}\n\n\n";
    
    
    int npsp=pseudo.size();
    for(int psp=0; psp < npsp; psp++) {
      pseudo[psp].print_pseudo(sysout);
    }
    sysout.close();
    
    //---------------------------HF and OPT outputs
    
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
  }//mpi-info.node==0

#ifdef USE_MPI
  MPI_Finalize();
#endif
  
}//main end;


//-------------------------------------------------------------------


void read_abinit_sys(std::istream & is, vector < Atom > & atoms, vector < Atom > & primatoms,
		     vector <Spline_pseudo_writer> & abinit_pseudo,
		     vector <double> & origin,
		     vector < vector < double> > & latvec,
		     vector < vector < double> > & primlatvec,
		     vector < int > & factor,
		     int & nelectrons,
		     int & spin,
		     double & eref,
		     int & debug) {
  // cout <<"Reading abinits system"<<endl;
  vector <string> uniquenames;
  string line;
  string space=" ";
  int natoms=0;
  origin[0]=0.0;origin[1]=0.0;origin[2]=0.0;
   
  vector <string> words;
  //read in atoms and latvec
  while(getline(is, line)) {
    words.clear();
    split(line, space, words);
    //if(words[0]=="Number")
    // cout <<words[0]<<" "<<words[1]<<" "<<words[2]<<" "<<words[3]<<" "<<words[4]<<" "<<words[4]<<endl;
    //reading "Number of electrons per primitive cell
   
    if(words[0]=="Number" && words[1]=="of" && words[2]=="electrons"){
      getline(is,line);
      words.clear();
      split(line, space, words);
      nelectrons=atoi(words[0].c_str());
      if(mpi_info.node==0)
	cout <<"Number of electrons per primitive cell "<<nelectrons<<endl;
    }

    //reading "Number of atoms per primitive cell"
    if(words[0]=="Number" && words[1]=="of" && words[2]=="atoms"){
      getline(is,line);
      words.clear();
      split(line, space, words);
      natoms=atoi(words[0].c_str());
      if(mpi_info.node==0)
	cout <<"Number of atoms per primitive cell "<<natoms<<endl;
    }

    //reading "Spin polarized:"
    if(words[0]=="Spin" && words[1]=="polarized:"){
      getline(is,line);
      words.clear();
      split(line, space, words);
      if(words[0]==".true.")
	spin=2;
      else if(words[0]==".false.")
	spin=1;
      else{
	cout <<"Did not understand Spin polarized: "<<words[0]<<endl;
	exit(1);
      }
      if(mpi_info.node==0)
	cout <<"Spin polarized: "<<words[0]<<endl;
    }

    //reading "Total energy (au per primitive cell)"
    if(words[0]=="Total" && words[1]=="energy" && words[2]=="(au"){
      getline(is,line);
      words.clear();
      split(line, space, words);
      eref=atof(words[0].c_str());
      if(mpi_info.node==0)
	cout <<"Total energy (au per primitive cell) "<<eref<<endl;
    }
    
    //reading "Atomic numbers and positions of atoms (au)"
    Atom tempatom;
    if(words[0]=="Atomic" && words[1]=="numbers" && words[3]=="positions")
      for(int i=0;i<natoms;i++){
	getline(is,line);
	words.clear();
	split(line, space, words);
	tempatom.charge=atoi(words[0].c_str());;
        tempatom.pos[0]=atof(words[1].c_str());
	tempatom.pos[1]=atof(words[2].c_str());
	tempatom.pos[2]=atof(words[3].c_str());
	tempatom.name=element_lookup_caps[int(tempatom.charge)];
	primatoms.push_back(tempatom);
	int nunique=uniquenames.size();
	bool unique=true;
	for(int i=0; i< nunique; i++) {
	  if(tempatom.name==uniquenames[i]) {
	    unique=false;
	    break;
	  }
	}
	if(unique) {
	  uniquenames.push_back(tempatom.name);
	}
      }
    


    //reading "Primitive lattice vectors (au)"
    vector < double> latv(3);
    if(words[0]=="Primitive" && words[1]=="lattice" && words[2]=="vectors")
      for(int i=0;i<3;i++){
	getline(is,line);
	words.clear();
	split(line, space, words);
        latv[0]=atof(words[0].c_str());
	latv[1]=atof(words[1].c_str());
	latv[2]=atof(words[2].c_str());
	primlatvec.push_back(latv);
      }
    if(words[0]=="G" && words[1]=="VECTORS"){
      //cout <<" done reading system part"<<endl;
      break;
    }
  }

  if(mpi_info.node==0)
    cout <<"Atoms in the primitive cell"<<endl; 
  for(unsigned int at=0; at < primatoms.size(); at++) {
    if(mpi_info.node==0)
      cout << primatoms[at].name << "  " << primatoms[at].pos[0] << "  " << primatoms[at].pos[1]
	   << "   " << primatoms[at].pos[2] << endl;
  }
  
  if(mpi_info.node==0 && debug){
    cout <<"converting the PPs"<<endl;
  }
  
  //There should be pseudopotential files for each of the atom names..
  const string comment="#";
  for(unsigned int file=0; file < uniquenames.size(); file++) {
    string pspfilename=uniquenames[file];
    ifstream pspfile(pspfilename.c_str());
    if(!pspfile) {
      cout << "Couldn't open Pseudopotential file: " << pspfilename << endl;
      exit(1);
    }
    string line;
    vector <string> words;
    Spline_pseudo_writer temppsp;
    temppsp.label=uniquenames[file];
    //reading OPIUM generated 
    getline(pspfile, line);
    words.clear();
    split(line, space, words);
    if( (words[0]=="OPIUM" && words[1]=="generated" )||(words[1]=="fhi98PP")){
      //reading second line of FHI code (something like this): 7.000  5.000    021003              zatom,zion,pspdat
      getline(pspfile, line);
      words.clear();
      split(line, space, words);
      double zatom=atof(words[0].c_str());
      if(element_lookup_caps[int(zatom)]!=uniquenames[file]){
	cout <<"Charge of atoms and its PP file do not match\n"<<endl;
	exit(1);
      }
      double effcharge=atof(words[1].c_str()); //effective charge
      //Go ahead and assign the effective charge to the atoms.
      for(unsigned int at=0; at< primatoms.size(); at++) {
	if(primatoms[at].name==temppsp.label)
	  primatoms[at].charge=effcharge;
      }
      if(mpi_info.node==0 && debug)
	cout <<"for "<<uniquenames[file]<<" found effective charge "<<effcharge<<endl;
      //reading the 3-rd line of FHI code (something like this):  6   7    3   2    467     0       pspcod,pspxc,lmax,lloc,mmax,r2well
      getline(pspfile, line);
      words.clear();
      split(line, space, words);
      int npoints=atoi(words[4].c_str());  //Number of points in this expansion
      int maxL=atoi(words[2].c_str());
      int locL=atoi(words[3].c_str());
      


      if(mpi_info.node==0 && debug)
	cout <<"Lmax: "<<givesymmetry(maxL)<<", Lloc: "<<givesymmetry(locL)<<", points on the grid: "<<npoints<<endl;
      vector <double> positions; //for this l-value
      vector <double> values;
      positions.resize(npoints);
      values.resize(npoints);
      vector < vector <double> > psp_pos;
      vector < vector <double> > psp_val;
      for(int i=0; i< maxL+1; i++) {
	while(getline(pspfile, line)){
	  //cout <<line<<endl;
	  words.clear();
	  split(line, space, words);
	  if(atoi(words[0].c_str())==npoints)
	    break;
	}
	for(int point=0; point < npoints; point++) {
	  getline(pspfile, line);
	  words.clear();
	  split(line, space, words);
	  positions[point]=atof(words[1].c_str());
	  values[point]=atof(words[3].c_str());;
	  //cout <<positions[point]<<" "<<values[point]<<endl;
	}
	psp_pos.push_back(positions);
	psp_val.push_back(values);
      }//i
      if(maxL!=locL){
	for(int p=0; p < npoints; p++) {
	  values[p]=psp_val[locL][p];
	  positions[p]=psp_pos[locL][p];
	}
	if(mpi_info.node==0 && debug)
	  cout <<"adding extra channel "<<endl;
	psp_pos.push_back(positions);
	psp_val.push_back(values);
      }
      
      
      //To get the psp's in the normal form, we:
      //Subtract the last l-value from the first ones
      //Add Z/r to the last one--this is done in the core code via a switch
      
      int newmaxL=psp_val.size();
      
      for(int i=0; i < newmaxL-1; i++) {
	for(int p=0; p < npoints; p++) {
	  psp_val[i][p]-= psp_val[newmaxL-1][p];
	}
      }

      //needs to be on the linear grid starting at x=0; this bellow works reasonably well
      double xmax=10.00;
      //cout <<"#xmin: "<<xmin<<" xmax: "<<xmax<<" number of points:  "<<x.size()<<endl;
      int npts=10000;
      positions.resize(npts);
      values.resize(npts);
      double spacing=xmax/npts;
      for(int l=0; l < newmaxL; l++){
	double xpos=0;
	for(int i=0;i<npts;i++){
	  positions[i]=xpos;
	  values[i]=linear_extrapolate(psp_pos[l],psp_val[l],xpos);
	  //cout <<positions[i]<<"  "<<values[i]<<endl;
	  xpos+=spacing;
	}
	temppsp.psp_pos.push_back(positions);
	temppsp.psp_val.push_back(values);
      }//l
      if(mpi_info.node==0 && debug)
	cout <<"Final number of channels "<<temppsp.psp_val.size()<<endl;
      abinit_pseudo.push_back(temppsp);
    }
    else{
      if(mpi_info.node==0)
	cout <<"This is not fhi98PP/OPIUM PP format (*.FHI), skipping the pseudopotential file: " << pspfilename<<endl;
    }
  }//for each unique atom

  //resizing everything to simulation cell
  nelectrons*=factor[0]*factor[1]*factor[2];
  eref*=factor[0]*factor[1]*factor[2];

  for(int i=0;i<3;i++){
    vector <double> vec(3);
    for(int j=0;j<3;j++)
      vec[j]=primlatvec[i][j]*factor[i];
    latvec.push_back(vec);
  }


  int primatomssize=atoms.size();
  for(int i=0;i<factor[0];i++){
    for(int j=0;j<factor[1];j++){
      for(int k=0;k<factor[2];k++){
	//if(!(i==0 && j==0 && k==0)){
       	  for(unsigned int at=0; at < primatoms.size(); at++){
	    Atom newatom;
	    newatom.name=primatoms[at].name;
	    newatom.charge=primatoms[at].charge;
	    for(int d=0;d<3;d++)
	      newatom.pos[d]=primatoms[at].pos[d]+i*primlatvec[0][d]+j*primlatvec[1][d]+k*primlatvec[2][d];
	    //if(at==0){
	      //cout <<primlatvec[][0]<<" "<<primlatvec[k][1]<<" "<<primlatvec[k][2]<<endl;
	    // cout <<newatom.name<<": "<<newatom.pos[0]<<" "<<newatom.pos[1]<<" "<<newatom.pos[2]<<endl;
	    //}
	    atoms.push_back(newatom);
	  }//at
	  //}
      }//k
    }//j
  }//i

  if(mpi_info.node==0 && debug){
    cout <<"Atoms in the simulation cell"<<endl; 
    for(unsigned int at=0; at < atoms.size(); at++) {
      if(mpi_info.node==0)
	cout << atoms[at].name << "  " << atoms[at].pos[0] << "  " << atoms[at].pos[1]
	     << "   " << atoms[at].pos[2] << endl;
    }
  }


}

//----------------------------------------------------------------------

void read_abinit_wf(string & filename, CBasis_function * abinitbasis, int & debug, int & norbs, int & node) {
  abinitbasis->read(filename, debug, norbs, node);
}

//----------------------------------------------------------------------

void plot_orbitals(string & outputname,  CBasis_function * abinitbasis, int & norbs, vector < Atom > & atoms,
                   vector < vector < double> > & latvec, vector < double> & origin, double & resolution, int & debug ){

  
  if(mpi_info.node==0)
    cout<<"#################### Ploting Bloch wave orbitals ###################"<<endl;
  vector < vector <double> > resolution_array(3); 
  vector <int> D_array1(3); //dummy array1
  for(int i=0;i<3;i++){
      double lenght=0;
      for(int j=0;j<3;j++){
	if(mpi_info.node==0 && debug)
	  cout <<latvec[i][j]<<" ";
	lenght+=latvec[i][j]*latvec[i][j];
      }
      if(mpi_info.node==0 && debug)
	cout <<endl;
      lenght=sqrt(lenght);
      if(mpi_info.node==0 && debug)
	cout <<"lenght "<<lenght<<" resolution "<<resolution<<endl;//"roundoff "<<roundoff(lenght/resolution)<<endl;
      D_array1[i]= roundoff(lenght/resolution); 
      for(int j=0;j<3;j++){
	resolution_array[i].push_back(latvec[i][j]/D_array1[i]);
      }
  }

  int npts=D_array1[0]*D_array1[1]*D_array1[2];
  double MBs=npts*8/(1024*1024);

  vector < complex <double> > grid;
  //calculate value of each molecular orbital at each grid point and store in an Array1
  // grid values with x=fastest running variable, and z=slowest
  if(mpi_info.node==0){
    cout<<"calculating "<<D_array1[0]<<"x"<<D_array1[1]<<"x"<<D_array1[2]<<" = "<<npts<<" grid points\n";
    cout<<"for "<<norbs<<" real orbitals\n";
    cout<<"using "<<MBs<<" Mb per orbital in memory or "<<MBs*norbs<<" Mb in total memory\n";
  }

  int cycles=int(norbs/2)+norbs%2;

#ifdef USE_MPI
  MPI_Barrier(MPI_COMM_WORLD);
#endif
    
  for(int orb=mpi_info.node;orb<cycles;orb+=mpi_info.nprocs){  
    vector < complex <double> > grid;
    vector < double > xyz(3);
    for(int xx=0;xx<D_array1[0];xx++){
      if(mpi_info.node==0){
	 cout<<100*xx/(D_array1[0]-1)<<" % of plot \n";
	 flush(cout);
      }
      for(int yy=0; yy<D_array1[1];yy++){
	for(int zz=0; zz<D_array1[2];zz++){
	  for(int d=0;d<3;d++){
	    xyz[d]=origin[d]+xx*resolution_array[0][d]+yy*resolution_array[1][d]+zz*resolution_array[2][d];
	    //cout <<xyz[0]<<" "<<xyz[1]<<" "<<xyz[2]<<endl;
	  }
	  complex <double> value;
	  abinitbasis->calcVal(xyz,value,orb);
	  grid.push_back(value);
       }//zz
      }//yy
    }//xx

    
    string basename,basename2;
    char strbuff[40];
    sprintf(strbuff, "%d", orb+1);
    basename2 = outputname;
    basename2 += ".orb";
    basename2 += strbuff;
    ofstream os;
    

    for(int part=0;part<2;part++){
      if(part==0)
	basename=basename2+".real";
      else
	basename=basename2+".imag";
      string cubename=basename+".cube";
      if(debug)
	cout <<" node: "<<mpi_info.node<<" is storing orbital: "<<2*orb+1+part<<" to file: "<<cubename<<endl;

      os.open(cubename.c_str());
      int natoms=atoms.size();
      os << "GOS plot output\n";
      os << "Molecular orbital " << 2*orb+part << endl;
      for(int d=0;d<3;d++)
	os << D_array1[d]<<"   "<< resolution_array[d][0] <<"  "<<resolution_array[d][1]<<"  "<<resolution_array[d][2]<< endl;
      for(int k=0;k<3;k++) //need 3 empty lines to mimic jeep like file
	os << endl;
      if(natoms){
	os << "  " << natoms<<endl;
	for(int at=0; at< natoms; at++) {
	  os << "   " << atoms[at].charge << "   0.0000    " << atoms[at].pos[0]
	     <<"    " << atoms[at].pos[1] << "   " << atoms[at].pos[2] << endl;
	}
      }
      else //if natoms=0 like in HEG
	os << "  1" << endl<<endl;
      os <<"  1" << endl;
      os <<"  "<< origin[0]<<"  "<< origin[2]<<"  "<< origin[4]<<endl;
      os <<"  "<<resolution_array[0][0]+resolution_array[1][0]+resolution_array[2][0]
	 <<"  "<<resolution_array[0][1]+resolution_array[1][1]+resolution_array[2][1] 
	 <<"  "<<resolution_array[0][2]+resolution_array[1][2]+resolution_array[2][2] <<endl;
      os <<endl<<endl;
      os << "  "<<D_array1[0]<<"  "<<D_array1[1]<<"  "<<D_array1[2]<<endl;
      os.setf(ios::scientific);
      for(int j=0; j< grid.size(); j++) {
	if(part==0)
	  os <<setw(20)<<setprecision(10)<<grid[j].real();
	else
	  os <<setw(20)<<setprecision(10)<<grid[j].imag();
	if(j%6 ==5) os << endl;
      }
      os << endl;
      os.unsetf(ios::scientific);
      os<<setprecision(6);
      os.close();
    }//part
  }//orb
  if(mpi_info.node==0)
    cout<<"#################### done ploting all orbitals ###################"<<endl;
#ifdef USE_MPI
  MPI_Barrier(MPI_COMM_WORLD);
#endif

}
