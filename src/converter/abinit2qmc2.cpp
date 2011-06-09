/*
 
Copyright (C) 2011 Lucas Wagner

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
#include "Pseudo_writer.h"
#include <fstream>
#include <cmath>
#include "elements.h"
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include "vecmath.h"
#include <iomanip>
using namespace std;

/*
 * -o          # Output base for qwalk files
 * -ao     # output file from abinit
 * */


//--------------------------------------------------------------
/*!
 * Return a list of k-points with which we can construct a wave function for the supercell given.  
 * We're not able to do this for all supercell matrices, so we error out if this
 * algorithm doesn't work.
 * */
void supercell_kpoints(const vector <vector <double> > & supercell, 
    vector <vector <double> > & kpoints) { 
  assert(supercell.size()==3);
  assert(supercell[0].size()==3);
  kpoints.clear();
  vector <vector <double> >  inversesupercell(3);
  for(int d=0; d< 3; d++) inversesupercell.resize(3);
  matrix_inverse(supercell,inversesupercell);
  vector <double> kpoint(3);

  int nkpts=1;
  for(int d=0; d< 3; d++) { 
    for(int d1=0; d1 < 3; d1++) {
      cout << inversesupercell[d][d1] << " ";
      nkpts*=max(fabs(supercell[d][d1]),1.0);
    }
    cout << endl;
  }

  int nsearch=10;
  for(int ii=-nsearch; ii < nsearch; ii++) 
  for(int jj=-nsearch; jj < nsearch; jj++) 
  for(int kk=-nsearch; kk < nsearch; kk++) { 
    for(int d=0; d< 3; d++) {
      kpoint[d]=ii*inversesupercell[0][d]
        +jj*inversesupercell[1][d]
        +kk*inversesupercell[2][d];
    }
    if(kpoint[0] >= 0 && kpoint[0] < 1  &&
        kpoint[1] >= 0 && kpoint[1] < 1 &&
        kpoint[2] >= 0 && kpoint[2] < 1) { 
      cout << "kpoint " << kpoint[0] << " " << kpoint[1] << " " << kpoint[2] 
        << endl;
      kpoints.push_back(kpoint);
    }
  }

  if(nkpts!=kpoints.size()) { 
    cout << "Did not find a number of kpoints equal to the number required by the supercell." << endl;
    cout << "nreq " << nkpts << " nfound " << kpoints.size() << endl;
    exit(1);
  }

}
//----------------------------------------------------------------------
void extend_supercell(const vector <vector <double> > & supercell,
    const vector <Atom> & primatoms, const vector < vector <double> > & primlat,
    vector <Atom> & superatoms, vector <vector <double> > & superlat) { 
 assert(supercell.size()==3);
 assert(supercell[0].size()==3);
 superlat=primlat;
 for(int d=0; d< 3; d++) { 
   for(int d1=0; d1 < 3; d1++) {
     superlat[d][d1]=0.0;
     for(int d2=0; d2 < 3; d2++) { 
       superlat[d][d1]+=supercell[d][d2]*primlat[d2][d1];
     }
     cout << superlat[d][d1] << " ";
   }
   cout << endl;
 }

  vector <vector <double> >  inversesuperlat(3);
  for(int d=0; d< 3; d++) inversesuperlat.resize(3);
  matrix_inverse(superlat,inversesuperlat);


 vector < vector <double> > deltas;

 int nsearch=10;
 for(int ii=-nsearch; ii< nsearch; ii++) 
 for(int jj=-nsearch; jj< nsearch; jj++)
 for(int kk=-nsearch; kk< nsearch; kk++) { 
   vector<double> delta(3);
   for(int d=0; d< 3; d++) { 
     delta[d]=ii*primlat[0][d]+jj*primlat[1][d]+kk*primlat[2][d];
   }
   vector <double> u(3);
   for(int d=0; d< 3; d++) u[d]=0.0;
   for(int d=0; d< 3; d++) {
     for(int d1=0; d1 < 3; d1++) { 
       u[d]+=inversesuperlat[d1][d]*delta[d];
     }
   }
   if(u[0]>=0 && u[0] < 1 && u[1]>=0 && u[1]<1 && u[2]>=0 && u[2]<1) {
     deltas.push_back(delta);
     cout << "delta " << delta[0] << "  " << delta[1] << " " << delta[2] << endl;
   }
 }
 superatoms.clear();

 for(vector<Atom>::const_iterator i=primatoms.begin(); i!=primatoms.end(); i++) { 
   for(vector <vector <double> >::iterator d=deltas.begin(); d!=deltas.end(); d++) { 
     Atom tmp_at=*i;
     for(int d1=0; d1< 3; d1++) { 
       tmp_at.pos[d1]+=(*d)[d1];
     }
     superatoms.push_back(tmp_at);
   }
 }

}


//----------------------------------------------------------------------

class WF_kpoint  { //all the information necessary to calculate the periodic part of the wave function at a given k-point.
private:
  vector <vector <double> > gvec;
  vector <vector <complex <double> > > coeff;
  //vector <double> kpoint;
  vector <double> occupation;
  vector <vector <double> > latvec;
public:
  int norbitals() { 
    return coeff.size();
  }
  int npw() { 
    return gvec.size();
  }
  void occ(vector <double> & occ_) { 
    occ_=occupation;
  }
  void read_wf_from_wfk(FILE * wffile,vector <vector <double> > & latvec);
  void evaluate_orbital(vector <double> & r, vector <complex <double> > & orbitals);

};




//----------------------------------------------------------------------
class Abinit_converter { 
  public:
    void readabinitout(string filename,vector <vector<double> > & selected_kpt);
    void add_psp(string filename);

    void write_files(string basename);
    void write_orbitals(string  filename);
    void write_sys(string filename);
    void write_slater(string filename);
  private:
    void read_wfk(string filename, vector <vector <double> > & selected_kpt);
    //void read_wf_from_wfk(FILE * wffile,vector<double> & occupation);
    void skip_wf(FILE * wffile);
    void assign_occupations(vector<double> &occupations,bool spin_polarized); 
    void prune_coefficients();
    //void evaluate_orbital(vector<double> & r, vector <complex <double> > & orbitals);
    void reassign_z(); 
  
    //Information from the conversion.
    Slat_wf_writer slater;
    Jastrow_wf_writer jast;
    vector <Atom> atoms;
    vector <vector<double> > latvec;
    vector <Spline_pseudo_writer> pspspline; //one for each atom type
/*
    vector < vector <double> > gvec;
    vector < vector <complex < double> > > coeff;
    vector <vector <double> > kpoint_orbital; //k-point of a given orbital
    vector <double> occupation; 
*/    
    vector <WF_kpoint> wavefunctions;
    vector <vector <double> > kpoint_orbital;
    vector<double> kpoint; //overall k-point
    bool complex_wavefunction;
    double grid_resolution;
};
//--------------------------------------------------------------

int main(int argc, char ** argv) { 
  string outbase="qwalk",abinit_out="abinit.out";

  vector < vector<double> > selected_kpt(1);
  selected_kpt[0].resize(3);
  selected_kpt[0][0]=selected_kpt[0][1]=selected_kpt[0][2]=0.0;

  for(int i=1; i< argc; i++) { 
    if(!strcmp(argv[i],"-o") && i<argc-1)
      outbase=argv[++i];
    else if(!strcmp(argv[i],"-ao") && i < argc-1) 
      abinit_out=argv[++i];
    else if(!strcmp(argv[i],"-kpoint") && i < argc-3) { 
      selected_kpt[0][0]=atof(argv[++i]);
      selected_kpt[0][1]=atof(argv[++i]);
      selected_kpt[0][2]=atof(argv[++i]);
    }
    else { 
      cout << "Error parsing command line " << endl;
      exit(1);
    }
  }

  Abinit_converter abconverter;
  abconverter.readabinitout(abinit_out,selected_kpt);
  abconverter.write_files(outbase);
  return 0;
}

//--------------------------------------------------------------


void Abinit_converter::readabinitout(string  filename,vector <vector<double> > & selected_kpt) { 
  ifstream is(filename.c_str());
  if(!is) { 
    cout << "couldn't open " << filename << endl;
    exit(2);
  }
  string line;
  vector<string> words;
  string sep=" ";
 // vector <int> typat;
 // vector<double> znucl;
  string outputroot;

  while(getline(is,line)) { 
    if(line.find("- pspini: atom type")!=string::npos) { 
      words.clear(); split(line,sep,words);
      assert(words.size() > 8);
      add_psp(words[8]);
    }
    else if(line.find("root for output files")!=string::npos) { 
      words.clear(); split(line,sep,words);
      assert(words.size()>6);
      outputroot=words[6];
    }
  }

  read_wfk(outputroot+"_DS1_WFK",selected_kpt);

}
//--------------------------------------------------------------
int read_int(FILE * file) {
  int i;
  fread(&i, sizeof(int), 1, file);
  return i;
}

double read_double(FILE * file) {
  double f;
  fread(&f, sizeof(double), 1, file);
  return f;
}
void uread(int * p,int ncount,FILE * file) { 
  fread(p,sizeof(int),ncount,file);
}
void uread(double * p,int ncount,FILE * file) { 
  fread(p,sizeof(double),ncount,file);
}

void uread(char * p,int ncount,FILE * file) { 
  fread(p,sizeof(char),ncount,file);
}

/*!
 * Fortran for some reason puts in a 4-byte header in front and back of every write statement,
 * so we have to clear it since we know when the writes come from the excellent
 * abinit documentation.
 * */
void clear_header(FILE * file) { 
  char buff[4];
  fread(buff,sizeof(char),4,file);
}
//----------------------------------------------------------------------
void Abinit_converter::read_wfk(string filename, 
    vector <vector <double> >  & selected_kpts) { 
  char codvsn[6];
  int headform,fform,bantot,date,intxc,ixc,natom,ngfft[3],nkpt,npsp,nspden;
  int nspinor,nsppol,nsym,ntypat,occopt,pertcase,usepaw;
  int usewvl,cplex;
  double acell[3],ecut,ecutdg,ecutsm,ecut_eff,qptn[3],rprimd[9],stmbias;
  double tphysel,tsmear;
  FILE * wffile=fopen(filename.c_str(), "r");
  if(wffile==NULL) { 
    cout << "Couldn't open " << filename << endl;
    exit(1);
  }
  //First block
  clear_header(wffile); 
  fread(codvsn,sizeof(char),6,wffile);
  headform=read_int(wffile);
  fform=read_int(wffile);
  clear_header(wffile);
  codvsn[6]='\0';
  cout << "codvsn " << codvsn << " headform " << headform << " fform " << fform << endl;

  clear_header(wffile);
  bantot=read_int(wffile);
  date=read_int(wffile);
  intxc=read_int(wffile);
  ixc=read_int(wffile);
  natom=read_int(wffile);
  fread(ngfft,sizeof(int),3,wffile);
  nkpt=read_int(wffile);
  nspden=read_int(wffile);
  nspinor=read_int(wffile);
  nsppol=read_int(wffile);
  nsym=read_int(wffile);
  npsp=read_int(wffile);
  ntypat=read_int(wffile);
  occopt=read_int(wffile);
  pertcase=read_int(wffile);
  usepaw=read_int(wffile);
  ecut=read_double(wffile);
  ecutdg=read_double(wffile);
  ecutsm=read_double(wffile);
  ecut_eff=read_double(wffile);
  fread(qptn,sizeof(double),3,wffile);
  fread(rprimd,sizeof(double),9,wffile);
  stmbias=read_double(wffile);
  tphysel=read_double(wffile);
  tsmear=read_double(wffile);
  usewvl=read_int(wffile);
  clear_header(wffile);
  cout << "nkpt " << nkpt << " bantot " << bantot << endl;
  cout << "date " << date << " natom " << natom << " nspinor " << nspinor 
    << " ecut " << ecut << " usepaw " << usepaw << " npsp " << npsp << endl;
  int * istwfk=new int[nkpt],*nband=new int[nkpt*nsppol],
         *npwarr=new int[nkpt], *so_psp=new int[npsp],
         *symafm=new int[nsym], *symrel=new int[9*nsym],
         *typat=new int[natom];

  double * kpt=new double[3*nkpt],*occ=new double[bantot],*tnons=new double[3*nsym],
        *znucltypat=new double[ntypat],*wtk=new double[nkpt];
  clear_header(wffile);
  uread(istwfk,nkpt,wffile);
  cout << "istwfk ";
  for(int i=0; i< nkpt; i++) cout << istwfk[i] << " ";
  cout << endl;
  uread(nband,nkpt*nsppol,wffile);
  uread(npwarr,nkpt,wffile);
  uread(so_psp,npsp,wffile);
  uread(symafm,nsym,wffile);
  uread(symrel,9*nsym,wffile);
  uread(typat,natom,wffile);
  uread(kpt,3*nkpt,wffile);
  uread(occ,bantot,wffile);
  uread(tnons,3*nsym,wffile);
  uread(znucltypat,ntypat,wffile);
  uread(wtk,nkpt,wffile);
  clear_header(wffile);

  cout << "typatom ";
  for(int i=0; i< natom; i++) cout << typat[i] << " ";
  cout << endl;

  for(int i=0; i< npsp; i++) { 
    clear_header(wffile);
    char title[132];
    uread(title,132,wffile);
    read_double(wffile); read_double(wffile); 
    for(int j=0; j< 5; j++) read_int(wffile);
    clear_header(wffile);
  }

  double residm, *xred=new double[3*natom],etotal,fermie;
  clear_header(wffile);
  residm=read_double(wffile);
  uread(xred,3*natom,wffile);
  etotal=read_double(wffile);
  fermie=read_double(wffile);
  clear_header(wffile);
  cout << "etot " << etotal << endl;
// At this point, we've read in the header part and gotten all the billions 
// of variables that abinit puts out.  We now print out what we need from the
// file, just as a check.
  cout << "Number of spin channels " << nsppol << endl;

  grid_resolution=0.4;

  cout << "Lattice vectors " << endl;
  latvec.resize(3);
  for(int i=0; i< 3; i++) { 
    latvec[i].resize(3);
    for(int j=0; j< 3; j++) { 
      cout << rprimd[i*3+j] << " ";
      latvec[i][j]=rprimd[i*3+j];
    }
    cout << endl;
  }
  cout << "Atom positions (reduced)" << endl;
  for(int at=0; at < natom; at++) { 
    Atom tmpat;
    tmpat.charge=znucltypat[typat[at]-1];
    tmpat.name=element_lookup_caps[int(tmpat.charge)];
    for(int i=0; i< 3; i++) { 
      tmpat.pos[i]=0.0;
      for(int j=0; j< 3; j++) {
        tmpat.pos[i]+=latvec[i][j]*xred[at*3+j];
      }
    }
    tmpat.print_atom(cout);
    atoms.push_back(tmpat);
  }
  cout << "K-points " << endl;
  kpoint.resize(3);
  wavefunctions.resize(nsppol*nkpt);
  for(int isppol=0; isppol < nsppol; isppol++) { 
    vector <double> kpoint_tmp(3);
    for(int k=0; k< nkpt; k++) { 
      for(int i=0; i< 3; i++) kpoint_tmp[i]=2*kpt[k*3+i];
      kpoint_orbital.push_back(kpoint_tmp);
      wavefunctions[isppol*nkpt+k].read_wf_from_wfk(wffile,latvec);
    }
  }
  
  complex_wavefunction=true;

  //assign_occupations(occupation,nsppol==2);
  //cout << "nup " << slater.nup << " ndown " << slater.ndown << endl;
  

  delete [] istwfk,nband,npwarr,so_psp,symafm,symrel,typat;
  delete [] kpt,occ,tnons,znucltypat,wtk;
}
//----------------------------------------------------------------------
//

void Abinit_converter::skip_wf(FILE * wffile) { 
  clear_header(wffile);
  int npw=read_int(wffile);
  int nspinor=read_int(wffile);
  int nband=read_int(wffile);
  //cout << "skipping npw " << npw << " nspinor " <<  nspinor << " nband " << nband << endl;
  clear_header(wffile);
  fseek(wffile,ftell(wffile)
      +4*sizeof(char)*2*(2+nband) //the nband+2 write statements
      +sizeof(int)*3*npw //The plane wave coordinates
      +sizeof(double)*2*nband //eigenvalues and occupation numbers
      +sizeof(double)*2*npw*nband, //The wf coefficients
      SEEK_SET);
}
//----------------------------------------------------------------------

void Abinit_converter::assign_occupations(vector<double>& occupations,bool spin_polarized) { 
 slater.detwt.resize(1);
 slater.detwt[0]=1.0;
 slater.occ_up.resize(1);
 slater.occ_down.resize(1);
 int nstates=occupations.size();
 //------spin polarized
 if(spin_polarized) {
   slater.calctype="UHF";
   assert(nstates%2==0);
   nstates/=2;
   for(int i=0; i< nstates; i++) { 
     if(fabs(occupations[i]-1.0) < 1e-3)  { 
       slater.occ_up[0].push_back(i+1);
     }
     else if(fabs(occupations[i]) > 1e-3) { 
       cout << "occupation !=0 or 1 for spin-polarized case" << endl;
       exit(1);
     }
     if(fabs(occupations[i+nstates]-1.0) < 1e-3) 
       slater.occ_down[0].push_back(i+nstates+1);
     else if(fabs(occupations[i+nstates]) > 1e-3) { 
       cout << "occupation != 0 or 1 for spin-polarized case" << endl;
       exit(1);
     }
   }
 }
 else  {  //spin unpolarized
   slater.calctype="RHF";
   for(int i=0; i< nstates; i++) { 
     if(fabs(occupations[i]-2.0) < 1e-3) { 
       slater.occ_up[0].push_back(i+1);
       slater.occ_down[0].push_back(i+1);
     }
     else if(fabs(occupations[i]) > 1e-3) { 
       cout << "occupation !=0 or 2 for spin unpolarized case" << endl;
       exit(1);
     }
   }
 }
 slater.nup=slater.occ_up[0].size();
 slater.ndown=slater.occ_down[0].size();
}
//--------------------------------------------------------------

void Abinit_converter::write_files(string basename) { 

  write_orbitals(basename+".orb");

  vector <vector <double> > supercell(3);
  for(int d=0; d< 3; d++) supercell[d].resize(3);
  supercell[0][0]=2.0; supercell[0][1]=0.0; supercell[0][2]=0.0;
  supercell[1][0]=0.0; supercell[1][1]=1.0; supercell[1][2]=0.0;
  supercell[2][0]=0.0; supercell[2][1]=0.0; supercell[2][2]=1.0;


  vector <Atom> oldatoms=atoms;
  vector <vector <double> > oldlatvec=latvec;
  extend_supercell(supercell,oldatoms,oldlatvec,atoms,latvec);
  vector <vector <double> > supercell_kpts;
  supercell_kpoints(supercell,supercell_kpts);
  int nkpts=wavefunctions.size();
  if(nkpts < supercell_kpts.size()) { 
    cout << "Not enough k-points for this supercell. Need " 
      << nkpts << " and only found " << atoms.size()/oldatoms.size() << endl;
    exit(2);
  }
  vector <double> tot_occupation;
  for(int k=0; k< nkpts; k++) { 
    bool matching_kpt=false;
    for(vector <vector <double> >::iterator skpt=supercell_kpts.begin(); 
        skpt!=supercell_kpts.end(); skpt++) { 
      bool match=true;
      for(int d=0;d < 3; d++) {
        if(fabs(2*(*skpt)[d]-kpoint_orbital[k][d]) > 1e-3)
          match=false;
      }
      if(match) {
        matching_kpt=true;
        break;
      }
    }
    if(matching_kpt) { 
      vector <double> occ;
      wavefunctions[k].occ(occ);
      tot_occupation.insert(tot_occupation.end(),occ.begin(),occ.end());
    }
    else { 
      vector <double> occ(wavefunctions[k].norbitals());
      for(vector<double>::iterator i=occ.begin(); i!=occ.end(); i++) *i=0.0;
      tot_occupation.insert(tot_occupation.end(),occ.begin(),occ.end());
    }
  }

  assign_occupations(tot_occupation,false);  

  reassign_z();
  write_sys(basename+".sys");
  write_slater(basename+".slater");

  string jast2outname=basename+".jast2";
  double basis_cutoff=find_basis_cutoff(latvec);
  Jastrow2_wf_writer jast2writer;
  jast2writer.set_atoms(atoms);
  
  
  ofstream jast2out(jast2outname.c_str());
  print_std_jastrow2(jast2writer, jast2out, basis_cutoff);
  jast2out.close();
  
}


//--------------------------------------------------------------

void Abinit_converter::prune_coefficients() { 
   //In principle, we can save a bit of time by taking out the G,-G pairs.
   //For the moment, I won't worry about it.
   //We may be able to save a lot of time by making the wave functions real
   //instead of detecting it post facto, as well.
}

//--------------------------------------------------------------

//Here we just dumbly write all the orbitals in wavefunctions to a file.
void Abinit_converter::write_orbitals(string  filename) { 
  slater.orbname=filename;
  //vector <complex < double> > orbitals(coeff.size());
  int totorbitals=0,maxorbitals=0;
  for(vector<WF_kpoint>::iterator w=wavefunctions.begin(); w!=wavefunctions.end();w++) { 
    totorbitals+=w->norbitals();
    maxorbitals=max(w->norbitals(),maxorbitals);
  }
  vector<complex <double> > orbitals(maxorbitals);
  vector <double> res(3);
  vector <int> npts(3);
  for(int d=0;d < 3; d++) {
    double len=0;
    for(int d1=0; d1 < 3; d1++) len+=latvec[d][d1]*latvec[d][d1];
    len=sqrt(len);
    npts[d]=int(len/grid_resolution)+1;
    res[d]=1.0/npts[d];
  }
  cout << "point resolution " << grid_resolution << endl;
  vector<double> r(3),prop(3);
  int count=0;
  long int ntotalpts=npts[0]*npts[1]*npts[2]*totorbitals;
  cout << "size of all orbitals: " << 2*ntotalpts*8.0/1024.0/1024.0 << " megabytes ("<< ntotalpts << " total values)" << endl;
  vector <complex<double> > allorbitals(ntotalpts);
  vector <complex<double> > orbnorm(totorbitals);
  for(vector<complex<double> >::iterator i=orbnorm.begin(); i!=orbnorm.end(); i++) 
    *i=0.0;

  for(int ii=0; ii < npts[0]; ii++) { 
    cout << "converting: " << 100*double(ii)/double(npts[0]) << "\% done" << endl;
    for(int jj=0; jj< npts[1]; jj++) { 
      for(int kk=0; kk< npts[2]; kk++) { 
        for(int i=0; i< 3; i++) r[i]=0.0;
        prop[0]=ii*res[0]; prop[1]=jj*res[1]; prop[2]=kk*res[2];
        for(int i=0; i< 3; i++) { 
          for(int j=0; j< 3; j++) { 
            r[i]+=prop[j]*latvec[i][j];
          }
        }
        vector< complex <double> >::iterator 
          ptr=allorbitals.begin()+kk+jj*npts[2]+ii*npts[2]*npts[1],
          orbn=orbnorm.begin();
        
        for(vector<WF_kpoint>::iterator w=wavefunctions.begin(); 
            w!=wavefunctions.end(); w++) { 
          w->evaluate_orbital(r,orbitals);
          vector< complex <double> >::iterator i=orbitals.begin();
          for(; i!= orbitals.end(); i++,orbn++) { 
            *ptr=*i;
            *orbn+=complex<double>(i->real()*i->real(),i->imag()*i->imag());
            ptr+=npts[1]*npts[2]*npts[0];
          }
          count++;
        }
      }
    }
  }

  for(vector<complex<double> >::iterator orbn=orbnorm.begin(); orbn!=orbnorm.end();
      orbn++) cout <<"norm " << *orbn << endl;

  ofstream os(filename.c_str());
  os.precision(15);
  os << "Orbitals file" << endl;
  os << "norbitals " << totorbitals << endl;
  os << "K-point of each orbital " << endl;
  vector<WF_kpoint>::iterator w=wavefunctions.begin();
  for(vector<vector <double> >::iterator k=kpoint_orbital.begin() ;
      k!=kpoint_orbital.end() && w!=wavefunctions.end(); k++,w++) {
    for(int i=0; i< w->norbitals(); i++) 
      os << (*k)[0] << " " << (*k)[1] << " " << (*k)[2] << endl;
  }
  os << "Lattice vectors " << endl;
  for(int i=0; i< 3; i++) { 
    for(int j=0; j< 3; j++) { 
      os << latvec[i][j] << " ";
    }
    os << endl;
  }
  os << "resolution ";
  for(int i=0; i< 3; i++) os << res[i] << " ";
  os << endl;
  os << "npoints ";
  for(int i=0; i< 3; i++) os << npts[i] << " ";
  os << endl;
  os << "orbitals follow (orbital,x,y,z indices) " << endl;
  vector<complex<double> >::iterator ptr=allorbitals.begin();
  for(int orb=0; orb < orbnorm.size(); orb++) { 
    int nxyz=npts[0]*npts[1]*npts[2];
    if(complex_wavefunction) { 
      for(int i=0; i< nxyz; i++) { 
        os << *ptr << "\n";
        ptr++;
      }
    }
    else { 
      if(orbnorm[orb].real() > orbnorm[orb].imag())  { 
        cout << "orb " <<orb << " real " << endl;
        for(int i=0; i< nxyz; i++) {
          os << ptr->real() << "\n";
          ptr++;
        }
      }
      else for(int i=0; i< nxyz; i++) { 
        os << ptr->imag() << "\n";
        ptr++;
      }
    }
  }
  
}
//--------------------------------------------------------------
void Abinit_converter::write_sys(string filename) { 
  ofstream os(filename.c_str());
  os << "System { periodic \n";
  os << "  nspin { " << slater.nup << " " << slater.ndown << " } \n";
  for(vector<Atom>::iterator at=atoms.begin(); at!=atoms.end(); at++) { 
    at->print_atom(os);
  }
  os << "LatticeVec { \n";
  for(int i=0; i< 3; i++) { 
    for(int j=0; j< 3; j++) { 
      os << latvec[i][j] << " ";
    }
    os << endl;
  }
  os <<  "}\n";
  os << "   kpoint { " << kpoint[0] << " " << kpoint[1] << " " << kpoint[2] << " } \n"; 
  
  os << "} ";


  for(vector<Spline_pseudo_writer>::iterator i=pspspline.begin(); i!=pspspline.end();i++) { 
    i->print_pseudo(os);
  }
}


//--------------------------------------------------------------
void Abinit_converter::write_slater(string filename) { 
  ofstream os(filename.c_str());
  if(!complex_wavefunction) 
    slater.orbtype="ORBITALS";
  else 
    slater.orbtype="CORBITALS";
  slater.mo_matrix_type="EINSPLINE_MO";
  slater.print_wavefunction(os);

}

//--------------------------------------------------------------

void Abinit_converter::add_psp(string filename) { 
  if(filename.find(".fhi")==string::npos) { 
    cout << "Seems like " << filename << " isn't a psp file \n";
    return;
  }
  ifstream is(filename.c_str());
  string line;
  getline(is,line);
  getline(is,line);
  vector<string> words; string sep=" ";
  split(line,sep,words);
  int zatom=atoi(words[0].c_str());
  int zion=atoi(words[1].c_str());
  getline(is,line); words.clear(); split(line,sep,words);
  int npoints=atoi(words[4].c_str());
  int lmax=atoi(words[2].c_str());
  int lloc=atoi(words[3].c_str());
  cout << "lmax " << lmax << " lloc " << lloc << " zatom " << zatom << " zion " << zion << endl;
  assert(lloc == lmax); //for now, we'll assume the highest is  the local one
                        //I think that QWalk also assumes this at the time of this writing, so 
                        //we'd have to change it there too.
  Spline_pseudo_writer psp;
  for(int l=0; l< lmax; l++) { 
    //find the start of the values
    double rcut=0;
    while(getline(is,line)) { 
      words.clear(); split(line,sep,words);
      if(atoi(words[0].c_str())==npoints) { 
        rcut=atof(words[1].c_str());
        break;
      }
    }
    //
    vector<double> pos(npoints),val(npoints);
    vector<double>::iterator p=pos.begin(),v=val.begin();
    for(;p!=pos.end() && v!=val.end(); p++,v++) { 
      words.clear(); getline(is,line); split(line,sep,words);
      *p=atof(words[1].c_str());
      *v=atof(words[3].c_str());
    }
    vector <double> upos,uval;
    make_uniform(pos,val,upos,uval,0.01);
    psp.psp_pos.push_back(upos);
    psp.psp_val.push_back(uval);
  }
  psp.label=element_lookup_caps[int(zatom)];
  psp.effcharge=zion;
  int newnpoints=psp.psp_pos[0].size();
  for(int i=0; i< lloc-1; i++) { 
    assert(psp.psp_pos[i].size()==psp.psp_pos[lloc-1].size());
    for(int j=0; j< newnpoints; j++) {
      assert(fabs(psp.psp_pos[i][j]-psp.psp_pos[lloc-1][j]) < 1e-5);
      psp.psp_val[i][j]-=psp.psp_val[lloc-1][j];
    }
  }

  psp.print_pseudo(cout);
  pspspline.push_back(psp);

}
//--------------------------------------------------------------
void Abinit_converter::reassign_z() { 
  for(vector<Spline_pseudo_writer>::iterator p=pspspline.begin(); 
      p!=pspspline.end(); p++) { 
    for(vector<Atom>::iterator a=atoms.begin(); 
        a!=atoms.end(); a++) {
      if(a->name==p->label) { 
        a->charge=p->effcharge;
        a->print_atom(cout);
      }
    }
  }
}


//----------------------------------------------------------------------

void WF_kpoint::read_wf_from_wfk(FILE * wffile,vector <vector <double> > & latvec_) { 
  latvec=latvec_;
  vector <vector <double> > gprim;
  matrix_inverse(latvec,gprim);
  
  clear_header(wffile);
  int npw=read_int(wffile);
  int nspinor=read_int(wffile);
  int nband=read_int(wffile);
  clear_header(wffile);
  cout << "npw " << npw << " nspinor " << nspinor << " nband " << nband << endl;
  if(gvec.size() > 0 && gvec.size()!=npw) { 
    cout << "Something wrong--number of g-vectors is changing" << endl;
    exit(1);
  }
  gvec.resize(npw);
  clear_header(wffile);
  //these are the reduced g-vectors.  They can be transformed by multiplying 
  //by gprim, which is the matrix inverse of rprimd obtained above.
  int * tmpgvecs=new int[3*npw];
  uread(tmpgvecs,3*npw,wffile);
  for(vector<vector <double> >::iterator g=gvec.begin(); g!=gvec.end(); g++) 
    g->resize(3);
  for(int g=0; g< npw; g++) { 
    for(int i=0; i< 3; i++) { 
      gvec[g][i]=0.0;
      for(int j=0; j< 3; j++) { 
        gvec[g][i]+=2*pi*gprim[j][i]*tmpgvecs[g*3+j];
      }
    }
  }
  delete [] tmpgvecs;
  //cout << "First few g-vectors : \n";
  //for(int g=0; g< 10; g++) { 
  //  cout << gvec[g][0] << " " << gvec[g][1] << " " << gvec[g][2] << endl;
  //}
  clear_header(wffile);
//coeff.resize(nband);

  double *eigen=new double[nband],*occ=new double[nband];
  clear_header(wffile);
  uread(eigen,nband,wffile);
  uread(occ,nband,wffile);
  clear_header(wffile);
  //cout << "occupations " << endl;
  for(int i=0; i< nband; i++) { 
  //  cout << occ[i] << " ";
    cout << "eigenvals " << eigen[i] << endl;
    occupation.push_back(occ[i]);
  }
  //cout << endl;
  delete[] eigen,occ;

  double * tmpcoeff=new double[2*npw*nspinor];
  vector <complex <double> > tmpvec(npw);
  
  assert(nspinor==1);
  for(int i=0; i< nband; i++) { 
    clear_header(wffile);
    uread(tmpcoeff,2*npw*nspinor,wffile);
    for(int ipw=0; ipw < npw; ipw++) { 
      tmpvec[ipw]=complex<double>(tmpcoeff[ipw*2],tmpcoeff[ipw*2+1]);
    }
    coeff.push_back(tmpvec);
    clear_header(wffile);
  }
  delete[] tmpcoeff;
}

//----------------------------------------------------------------------
void WF_kpoint::evaluate_orbital(vector<double> & r, 
    vector <complex< double> > & orbitals) { 
  assert(orbitals.size()==coeff.size());
  int ngvec=gvec.size();
  vector < complex<double> > basis(ngvec);
  vector <complex <double> >::iterator b=basis.begin();
  for(vector<vector<double> >::iterator g=gvec.begin(); g!=gvec.end(); g++,b++) { 
    double dot=0.0;
    for(int d=0; d< 3; d++) { 
      dot+=(*g)[d]*r[d];
    }
    *b=complex<double>(cos(dot),sin(dot));
  }
  
  vector< complex < double> >::iterator o=orbitals.begin();
  vector<vector<complex <double> > >::iterator co=coeff.begin();
  for(;o!=orbitals.end() && co!=coeff.end(); co++,*o++) { 
    *o=0.0;
    vector <complex<double> >::iterator co_j=co->begin();
    b=basis.begin();
    for(;b!=basis.end() && co_j!=co->end(); b++,co_j++) { 
      *o+=(*co_j)*(*b);
    }
  }
}

//----------------------------------------------------------------------
