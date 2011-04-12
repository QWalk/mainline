/*
 
Copyright (C) 2009 Lucas Wagner

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
#include "vecmath.h"
using namespace std;

/*
 * -files      # files file for abinit
 * -o          # Output base for qwalk files
 * -wfdata     # wf output from abinit
 * */


//--------------------------------------------------------------

class Abinit_converter { 
  public:
    void readfile(string filename);
    void readabinitout(string filename);
    void add_psp(string filename);
    void write_files(string basename);
    void write_orbitals(string  filename);
    void write_sys(string filename);
    void write_slater(string filename);
  private:
    void prune_coefficients();  //Cancel out the G and -G degeneracies
    void read_wfk(string filename);
    void read_wf_from_wfk(FILE * wffile);
    
    void evaluate_orbital(vector<double> & r, vector <complex <double> > & orbitals);
    void reassign_z(); 
    //Information from the conversion.
    Slat_wf_writer slater;
    Jastrow_wf_writer jast;
    vector <Atom> atoms;
    vector <vector<double> > latvec;
    vector <Spline_pseudo_writer> pspspline; //one for each atom type
    vector < vector <double> > gvec;
    vector < vector <complex < double> > > coeff;
    vector<double> kpoint;
    bool spin_polarized;
    double grid_resolution;
};
//--------------------------------------------------------------

int main(int argc, char ** argv) { 

  string files="",outbase="qwalk",wfdata="qwalk.in";

  for(int i=1; i< argc; i++) { 
    if(!strcmp(argv[i],"-o") && i<argc-1)
      outbase=argv[++i];
    else if(!strcmp(argv[i],"-wfdata") && i < argc-1) 
      wfdata=argv[++i];
    else if(!strcmp(argv[i],"-files") && i < argc-1)
      files=argv[++i];
    else { 
      cout << "Error parsing command line " << endl;
      exit(1);
    }
  }


  Abinit_converter abconverter;
  abconverter.readabinitout("abinit.out");
  abconverter.write_files("nwabinit");
  exit(0);
  abconverter.readfile(wfdata);
  ifstream is(files.c_str());
  if(!is) { cout << "Couldn't open files file " << files << endl; exit(1); }
  string line;
  while(getline(is,line)) { 
    abconverter.add_psp(line);
  }
  abconverter.write_files(outbase);
  return 0;
}

//--------------------------------------------------------------

//Only read in the information from the converter output file.  All other processing goes in other
//functions.
void Abinit_converter::readfile(string  filename) { 
  ifstream is(filename.c_str());
  if(!is) { cout << "Couldn't open " << filename << endl; exit(1); }
  string line;
  vector <string> words;
  int nelectrons=0;
  string sep=" ";
  while(getline(is,line)) { 
    words.clear();
    if(line=="Number of electrons per primitive cell") {
      getline(is,line);
      nelectrons=atoi(line.c_str());
      cout <<  "number of electrons " <<  nelectrons << endl;
      slater.nup=nelectrons/2;
      slater.ndown=nelectrons/2;
    }
    else if(line=="GEOMETRY") { 
      getline(is,line); getline(is,line); getline(is,line);
      int nat=atoi(line.c_str());
      cout << "Number of atoms " << nat << endl;
      getline(is,line);
      for(int i=0; i< nat; i++) { 
        words.clear();
        getline(is,line);
        cout << "line " << line; 
        split(line,sep,words);
        Atom at;
        at.charge=atoi(words[0].c_str());
        for(int d=0;d < 3; d++) at.pos[d]=atof(words[d+1].c_str());
        at.name=element_lookup_caps[int(at.charge)];
        atoms.push_back(at);
      }
    }
    else if(line=="Primitive lattice vectors (au)") { 
      latvec.resize(3);
      for(int d1=0; d1< 3; d1++) { 
        latvec[d1].resize(3);
        words.clear();
        getline(is,line);
        split(line,sep,words);
        for(int d2=0; d2 < 3; d2++) latvec[d1][d2]=atof(words[d2].c_str());
      }
    }
    else if(line=="G VECTORS") { 
      getline(is,line); getline(is,line); getline(is,line);
      int ngvec=atoi(line.c_str());
      getline(is,line); cout << line << endl; assert(line=="Gx Gy Gz (au)");
      gvec.resize(ngvec);
      for(vector< vector<double> >::iterator i=gvec.begin(); i!=gvec.end(); i++) { 
        getline(is,line);
        words.clear(); split(line,sep,words);
        for(vector<string>::iterator j=words.begin(); j!=words.end(); j++) 
          i->push_back(atof(j->c_str()));
      }
    }
    else if(line=="WAVE FUNCTION") { 
      getline(is,line); getline(is,line); getline(is,line);
      int nkpts=atoi(line.c_str());
      assert(nkpts==1);
      getline(is,line); getline(is,line);
      words.clear(); split(line,sep,words);
      int kptnum=atoi(words[0].c_str());
      int nbandsup=atoi(words[1].c_str());
      int nbandsdown=atoi(words[2].c_str());
      for(int d=0; d< 3; d++) kpoint.push_back(atof(words[3+d].c_str()));
      assert(nbandsdown==0);
      coeff.resize(nbandsup);
      for(vector<vector < complex < double> > >::iterator band=coeff.begin(); band!=coeff.end(); band++) {
        getline(is,line); getline(is,line); getline(is,line);
        cout << line << endl;
        assert(line=="Eigenvector coefficients");
        band->resize(gvec.size());
        for(vector< complex <double> >::iterator c=band->begin(); c!=band->end(); c++) { 
          is >> *c;
        }
        is.ignore(180,'\n');
      }
    }
    else if(line=="Plane wave cutoff (au)") { 
      getline(is,line);
       //seems safe (magic number, I know!)
      grid_resolution=3.0/atof(line.c_str());    
    }
    else if(line=="Spin polarized:") { 
      getline(is,line);
      if(line.find(".true.")!=string::npos) 
        spin_polarized=true;
      else spin_polarized=false;
    }
  }

  prune_coefficients();
  //add_psp("pp.Si.burkatzki.fhi");
}

//--------------------------------------------------------------


void Abinit_converter::readabinitout(string  filename) { 
  ifstream is(filename.c_str());
  if(!is) { 
    cout << "couldn't open " << filename << endl;
    exit(2);
  }
  string line;
  vector<string> words;
  string sep=" ";
  vector <int> typat;
  vector<double> znucl;
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

  read_wfk(outputroot+"_DS1_WFK");

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
 * Fortran for some reason puts in a 4-byte header in front of every write statement,
 * so we have to clear it since we know when the writes come from the excellent
 * abinit documentation.
 * */
void clear_header(FILE * file) { 
  char buff[4];
  fread(buff,sizeof(char),4,file);
}

void Abinit_converter::read_wfk(string filename) { 
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
  assert(nsppol==1);
  if(nsppol==1) spin_polarized=false;
  int nelectrons=0;
  for(int i=0; i< nband[0]; i++) { 
    nelectrons+=int(occ[i]);
  }
  slater.nup=nelectrons/2;
  slater.ndown=nelectrons/2;
  cout << "nelectrons " << nelectrons << endl;

  grid_resolution=0.2;

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
        tmpat.pos[i]+=latvec[i][j]*xred[at*3+i];
      }
    }
    tmpat.print_atom(cout);
    atoms.push_back(tmpat);
  }

  cout << "K-points " << endl;
  for(int k=0; k< nkpt; k++) { 
    for(int i=0; i<3; i++) { 
      cout << kpt[k*3+i] << " ";
    }
    cout << endl;
  }
  
  delete [] istwfk,nband,npwarr,so_psp,symafm,symrel,typat;
  delete [] kpt,occ,tnons,znucltypat,wtk;
  read_wf_from_wfk(wffile);
}
//----------------------------------------------------------------------
void Abinit_converter::read_wf_from_wfk(FILE * wffile) { 
  vector <vector <double> > gprim;
  matrix_inverse(latvec,gprim);
  cout << "Inverse Matrix " << endl;
  for(int i=0; i< 3; i++) { 
    for(int j=0; j< 3; j++) cout << gprim[i][j] << " ";
    cout << endl; 
  }
  
  clear_header(wffile);
  int npw=read_int(wffile);
  int nspinor=read_int(wffile);
  int nband=read_int(wffile);
  clear_header(wffile);
  gvec.resize(npw);
  clear_header(wffile);
  //these are the reduced g-vectors.  They can be transformed by multiplying 
  //by gprim, which is the matrix inverse of rprimd obtained above.
  int * tmpgvecs=new int[3*npw];
  uread(tmpgvecs,3*npw,wffile);
  for(vector<vector <double> >::iterator g=gvec.begin(); g!=gvec.end(); g++) g->resize(3);
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
  coeff.resize(nband);

  double *eigen=new double[nband],*occ=new double[nband];
  clear_header(wffile);
  uread(eigen,nband,wffile);
  uread(occ,nband,wffile);
  clear_header(wffile);
  cout << "occupations " << endl;
  for(int i=0; i< nband; i++) { 
    cout << occ[i] << " ";
  }
  cout << endl;
  delete[] eigen,occ;

  double * tmpcoeff=new double[2*npw*nspinor];
  assert(nspinor==1);
  for(int i=0; i< nband; i++) { 
    clear_header(wffile);
    uread(tmpcoeff,2*npw*nspinor,wffile);
    coeff[i].resize(npw);
    for(int ipw=0; ipw < npw; ipw++) { 
      coeff[i][ipw]=complex<double>(tmpcoeff[ipw*2],tmpcoeff[ipw*2+1]);
    }
    clear_header(wffile);
  }
  delete[] tmpcoeff;
  //cout << "first few coefficients from each band \n";
  //for(int i=0; i< nband; i++) { 
  //  cout << "band " <<i << endl;
  //  for(int ipw=0; ipw < 10; ipw++) { 
  //    cout << coeff[i][ipw] << endl;
  //  }
  //   cout << endl;
  // }

}


//--------------------------------------------------------------

void Abinit_converter::write_files(string basename) { 
  reassign_z();
  write_orbitals(basename+".orb");
  write_sys(basename+".sys");
  write_slater(basename+".slater");
}


//--------------------------------------------------------------

void Abinit_converter::prune_coefficients() { 
   //In principle, we can save a bit of time by taking out the G,-G pairs.
   //For the moment, I won't worry about it.
   //We may be able to save a lot of time by making the wave functions real
   //instead of detecting it post facto, as well.
}

//--------------------------------------------------------------

void Abinit_converter::evaluate_orbital(vector<double> & r, 
    vector <complex< double> > & orbitals) { 
  assert(orbitals.size()==coeff.size());
  int ngvec=gvec.size();
  vector < complex<double> > basis(ngvec);
  vector <complex <double> >::iterator b=basis.begin();
  for(vector<vector<double> >::iterator g=gvec.begin(); g!=gvec.end(); g++,b++) { 
    double dot=(*g)[0]*r[0]+(*g)[1]*r[1]+(*g)[2]*r[2];
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

//--------------------------------------------------------------

void Abinit_converter::write_orbitals(string  filename) { 
  slater.orbname=filename;
  vector <complex < double> > orbitals(coeff.size());
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
  long int ntotalpts=npts[0]*npts[1]*npts[2]*orbitals.size();
  cout << "size of all orbitals: " << ntotalpts*8.0/1024.0/1024.0 << " megabytes ("<< ntotalpts << ") total points" << endl;
  complex<double> * allorbitals=new complex<double>[ntotalpts];
  vector <complex<double> > orbnorm(coeff.size());
  for(vector<complex<double> >::iterator i=orbnorm.begin(); i!=orbnorm.end(); i++) 
    *i=0.0;

  for(int ii=0; ii < npts[0]; ii++) { 
    cout << "ii " << ii << endl;
    for(int jj=0; jj< npts[1]; jj++) { 
      for(int kk=0; kk< npts[2]; kk++) { 
        for(int i=0; i< 3; i++) r[i]=0.0;
        prop[0]=ii*res[0]; prop[1]=jj*res[1]; prop[2]=kk*res[2];
        for(int i=0; i< 3; i++) { 
          for(int j=0; j< 3; j++) { 
            r[i]+=prop[j]*latvec[i][j];
            //os << r[j] << " ";
          }
        }
        evaluate_orbital(r,orbitals);
        complex<double> * ptr=allorbitals+kk+jj*npts[2]+ii*npts[2]*npts[1];
        vector< complex <double> >::iterator i=orbitals.begin(),orbn=orbnorm.begin();
        for(; i!= orbitals.end(); i++,orbn++) { 
          *ptr=*i;
          *orbn+=complex<double>(i->real()*i->real(),i->imag()*i->imag());
          ptr+=npts[1]*npts[2]*npts[0];
        }
        //cout << r[0] <<  " " << r[1] << " " << r[2] << " "
        //  << orbitals[15] << endl;
        //cout << orbitals[0] << " " ;
        //if(count%6==5) cout << endl;
        count++;
      }
    }
  }

  for(vector<complex<double> >::iterator orbn=orbnorm.begin(); orbn!=orbnorm.end();
      orbn++) cout <<"norm " << *orbn << endl;

  ofstream os(filename.c_str());
  os << "Orbitals file" << endl;
  os << "norbitals " << orbitals.size() << endl;
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
  complex<double> * ptr=allorbitals;
  for(int orb=0; orb < orbnorm.size(); orb++) { 
    int nxyz=npts[0]*npts[1]*npts[2];
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
  delete [] allorbitals;
  
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
  os << "   kpoint { 0 0 0 } \n";
  
  os << "} ";


  for(vector<Spline_pseudo_writer>::iterator i=pspspline.begin(); i!=pspspline.end();i++) { 
    i->print_pseudo(os);
  }
}


//--------------------------------------------------------------
void Abinit_converter::write_slater(string filename) { 
  ofstream os(filename.c_str());
  slater.orbtype="ORBITALS";
  slater.mo_matrix_type="EINSPLINE_MO";
  slater.calctype="RHF";
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
//--------------------------------------------------------------
