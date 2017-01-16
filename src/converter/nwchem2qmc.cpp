#include "converter.h"
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>
#include <fstream>
#include <cmath>
#include "basis_writer.h"
#include "Pseudo_writer.h"
#include "wf_writer.h"
/*
   Unlike gamess, nwchem prints out orbitals in a binary file called
   [input].movecs. So there is an additional step to first convert
   the .movecs to ASCII file. Go to nwchem_directory/contrib/mov2asc/,
   then make it. We need executable file called "mov2asc".

   Then, run
   mov2asc [large number]  [input].movecs [input].nw.vecs.
   Be sure to name the nwchem output file [input].nw.out. So
   [input].nw.vecs and [input].nw.out are the two files we need. Now run
   /qwalk_directory/src/converter/nwchem2qmc [input],
   this generates all the files needed by qwalk.

   I tested lots of cases, eg. atoms molecules with different kinds of
   ecps and different methods(RHF,ROHF,DFT), and the results matches
   well. Currently it only supports 6D,10F....basis. 
   */
using namespace std;


void usage(const char * name) {
  cout << "usage: " << name <<   " <options> <output> " << endl;
  cout << "Where options can be: \n";
  cout << "-o          Base name for your run case\n";
  exit(1);
}


void read_nwchem_vecs(string & outputfilename,
    string & vecfilename,
    vector <Atom> & atoms,
    Slat_wf_writer & slwriter,
    vector <Gaussian_basis_set> & basis,
    vector <string> & basis_function_type,
    vector < vector <double> > & moCoeff
    );

void read_nwchem_ecp(string & outputfilename,
    vector <Atom> & atoms,
    vector <Gaussian_pseudo_writer> & ecp,
    vector <double> & ecpremoved);

void read_nwchem_basis( string & outputfilename,
    vector <Atom> & atoms,
    vector <Gaussian_basis_set> & basis,
    vector <string> & basis_function_type);

void read_nwchem_out(string & filename, vector<Atom> & atoms,
    Slat_wf_writer & slwriter);

void switch_moCoeff_order(string & outputfilename, int num_basis_func,
    vector <string> basis_function_type,  
    vector < vector <double> > & tmpCoeff,
    vector < vector <double> > & moCoeff);



int main(int argc, char ** argv)
{
  string infilename;
  infilename=argv[1];
  string outputname;

  for(int i=1; i< argc-1; i++) {
    if(!strcmp(argv[i], "-o") && argc > i+1) {
      outputname=argv[++i];
    }
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
  if (outputname=="")  { outputname=infilename; } 

  string outputfilename=infilename+".nw.out";
  string vecfilename=infilename+".nw.vecs";
  


  vector <Atom> atoms;
  vector <Gaussian_basis_set>  basis;
  vector <string> basis_function_type;
  vector <Gaussian_pseudo_writer> ecp;
  Slat_wf_writer slwriter;
  slwriter.mo_matrix_type="CUTOFF_MO"; 
  int nelectrons;

  // read system and wf type    
  read_nwchem_out(outputfilename,atoms,slwriter);

  // read basis 
  read_nwchem_basis( outputfilename, atoms, basis, basis_function_type);
  // read vecs
  vector < vector< double > > tmpCoeff;
  vector < vector< double > > moCoeff;
  read_nwchem_vecs(outputfilename, vecfilename,
      atoms, slwriter, basis, basis_function_type, moCoeff);
  //  cout<<"hello "<<moCoeff[0][0]<<" "<<moCoeff[0][1]<<endl;   

  // read ecp
  vector <double> ecpremoved;
  read_nwchem_ecp(outputfilename, atoms, ecp, ecpremoved);
  //   cout<<ecp[0].coefficients[0][0] <<" "<<ecp[0].coefficients[0][1]<<" "<< ecp[0].coefficients[0][2] << endl ;
  int natoms=atoms.size();
  //   int necp=ecp.size();
  //   for (int at=0; at < natoms; at++)
  //   { for (int psp=0; psp < necp; psp++ )
  //     {  if (atoms[at].name==ecp[psp].label)
  //       {  atoms[at].charge-=ecpremoved[psp];
  //          slwriter.nup-=int(ecpremoved[psp]/2);
  //          slwriter.ndown-=int(ecpremoved[psp]/2);


  //       }
  //     }   
  //   }
  nelectrons=slwriter.nup+slwriter.ndown;


  vector < Center > centers;
  vector <int> nbasis;
  //  cout<<"hello "<<basis[0].nfunc()<<endl;
  centers.resize(atoms.size());
  nbasis.resize(natoms);
  for(int at=0; at < natoms; at++) {
    for(int i=0; i< 3; i++) centers[at].pos[i]=atoms[at].pos[i];
    centers[at].equiv_atom=at;
    centers[at].name=atoms[at].name;
    nbasis[at]=basis[atoms[at].basis].nfunc();
  }
  //print out the qmc input file  
  cout<<"Wring QMC input files...";
  string orboutname=outputname+".orb";
  slwriter.orbname=orboutname;
  string basisoutname=outputname+".basis";
  slwriter.basisname=basisoutname; 

  //orbital output
  ofstream orbout(orboutname.c_str());
  print_orbitals(orbout, centers, nbasis, moCoeff);
  orbout.close();

  //basis output
  ofstream basisout(basisoutname.c_str());
  int nbas=basis.size();
  for(int bas=0; bas < nbas; bas++) {
    basisout << "BASIS { \n";
    basis[bas].print_basis(basisout);
    basisout << "}\n\n\n";
  }
  basisout.close();


  //slater output
  string slateroutname=outputname+".slater";
  ofstream slaterout(slateroutname.c_str());
  slwriter.print_wavefunction(slaterout);
  slaterout.close();

  // jastrow 2 ouput
  string jast2outname=outputname+".jast2";

  double basis_cutoff=15; //arbitrary cutoff
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

  //system output
  string sysoutname=outputname+".sys";
  ofstream sysout(sysoutname.c_str());
  sysout << "SYSTEM { MOLECULE \n";
  sysout << "  NSPIN { " << slwriter.nup << "  "
    << slwriter.ndown << " } \n";
  for(int at=0; at <natoms; at++) {
    atoms[at].print_atom(sysout);
  }

  sysout << "}\n\n\n";

  int necp=ecp.size();
  for(int ecp_index=0; ecp_index < necp; ecp_index++) {
    ecp[ecp_index].print_pseudo(sysout);
  }

  //  sysout.close();

  cout<< "Done with preparing QMC input files" <<endl;   
  return 0;

}   

void read_nwchem_out(string & outputfilename,
    vector <Atom> & atoms,
    Slat_wf_writer & slwriter)
{  
  ifstream inFile;
  inFile.open(outputfilename.c_str());
  if (!inFile) {
    cout << "Unable to open output file" << endl;
    exit(1);
  }
  cout<<"Reading system geometry and wavefunction type..." ;
  string line;
  string space=" ";

  vector <string> words;
  vector <double> emptyvector;
  slwriter.nup=-1;
  slwriter.ndown=-1;
  int read_atoms=0;
  int read_units=0;
  double length_scale;    
  while(getline(inFile,line))
  {  words.clear();
    split(line,space,words);
    //read length unit
    if (words.size() > 4 && words[0]=="Output" && words[1]=="coordinates" && read_units==0)
    {  if (words[3]=="a.u.")
      { length_scale=1.000000; read_units=1; }
      else if(words[3]=="angstroms")
      { length_scale=1.889725989; read_units=1; }
    }
    //read atoms
    if(words.size() == 6 && words[2]=="Charge" && words[3]=="X" && words[4]=="Y" && read_atoms==0)
    {  getline(inFile,line);


      Atom tempatom;
      while(getline(inFile,line))
      {  words.clear();
        split(line, space, words);
        if(words.size()==0) { read_atoms=1;  break;}
        tempatom.name=words[1];
        //     cout<<"hello "<< tempatom.name <<endl;
        tempatom.charge=atof(words[2].c_str());
        for (int i=0; i<3; i++) 
        {
          tempatom.pos[i]=atof(words[i+3].c_str())*length_scale;
        }  
        atoms.push_back(tempatom);  
      }

      //  atoms.push_back(tempatom);
    }

    //atoms done
    //     cout<<"hello "<< atoms[0].name <<"hello "<<atoms[1].name<<endl;
    //   else if(words.size() > 0 && words[0]=="functions" )
    //   {  basis.size()=words[2]; }

    //read electron configuration

    if(words.size()==4 && words[0]=="closed" && words[1]== "shells" )
    {  slwriter.nup=atoi(words[3].c_str());
      slwriter.ndown=atoi(words[3].c_str());
    }

    if(words.size()==4 && words[0]=="open" && words[1]=="shells" )
    {  slwriter.nup=slwriter.nup+atoi(words[3].c_str()); }
    if(words.size()==4 && words[0]=="Alpha" && words[1]=="electrons")
    {  slwriter.nup=atoi(words[3].c_str()); }
    if(words.size()==4 && words[0]=="Beta" && words[1]=="electrons") 
    {  slwriter.ndown=atoi(words[3].c_str()); } 
    //electron configuation done 

    //read wavefunction type
    if(words.size() == 3  && words[0]=="wavefunction" && words[2]=="RHF")
    {  slwriter.calctype="RHF"; }
    if(words.size() == 3 && words[0]=="wavefunction" && words[2]=="ROHF")
    {  slwriter.calctype="ROHF"; }
    if(words.size() == 4 && words[0]=="SCF" && words[3]== "DFT")
    { getline(inFile,line);
      words.clear();
      split(line,space,words);
      if (words.size()==4 && words[2]=="closed" && words[3]=="shell.")
      { slwriter.calctype="RHF"; }
      else if (words.size()==4 && words[2]=="spin"  && words[3]=="polarized.")
      { slwriter.calctype="UHF";  }
      else { cout<< "Couldn't understand DFT wavefunction type"; }
    }



    words.clear();
  }
  inFile.close();
  inFile.clear();

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

  cout<<"Done. "<<endl;
}



void read_nwchem_basis( string & outputfilename,
    vector <Atom> & atoms,
    vector <Gaussian_basis_set> & basis,
    vector <string> & basis_function_type)
{  
  ifstream inFile;
  inFile.open(outputfilename.c_str());
  if (!inFile) {
    cout << "Unable to open output file";
    exit(1);
  }
  cout<<"Reading basis set...";
  cout.flush();
  string line;
  string space=" ";

  vector <string> words;
  vector <double> emptyvector;
  int read_basis=0;
  int read_basis_label=0;
  int basis_index=0; 

  while (getline(inFile,line))
  {  words.clear();
    split(line,space,words);
    if (words.size()>5 && words[0]=="Basis" && read_basis==0)
    {  
      int tmpbasis_index=0; //index for each ao basis for an atom
      int basis_index=-1;  //index for each atom
      Gaussian_basis_set tmpbasis;
      while(getline(inFile,line))
      {  

        words.clear();
        split(line,space,words);
        if(words.size()>0 && words[0]=="Summary") 
        {
          read_basis=1;
          break;
        }
        if(words.size()==2 && words[0]=="Exponent" && words[1]=="Coefficients")
        { basis.push_back(tmpbasis);   
          tmpbasis_index=0;
          basis_index=basis_index+1;
        }
        if(words.size()==4 && tmpbasis_index == atoi(words[0].c_str()))  //continue to read existing basis set
        {      basis[basis_index].exponents[atoi(words[0].c_str())-1].push_back(atof(words[2].c_str()));
          basis[basis_index].coefficients[atoi(words[0].c_str())-1].push_back(atof(words[3].c_str()));
        }

        if (words.size()==4 && tmpbasis_index == (atoi(words[0].c_str())-1 ))  // read a new ao basis set
        {  basis[basis_index].exponents.push_back(emptyvector);
          basis[basis_index].coefficients.push_back(emptyvector);

          if (words[1]=="D")  
          {  basis[basis_index].types.push_back("6D"); 
            basis[basis_index].exponents[atoi(words[0].c_str())-1].push_back(atof(words[2].c_str()));
            //            cout<<"hello"<<endl;  
            basis[basis_index].coefficients[atoi(words[0].c_str())-1].push_back(atof(words[3].c_str()));

          }
          else if (words[1]=="F")
          {  basis[basis_index].types.push_back("10F");
            basis[basis_index].exponents[atoi(words[0].c_str())-1].push_back(atof(words[2].c_str()));
            basis[basis_index].coefficients[atoi(words[0].c_str())-1].push_back(atof(words[3].c_str()));

          }
          else if (words[1]=="G")  
          {  basis[basis_index].types.push_back("15G");
            basis[basis_index].exponents[atoi(words[0].c_str())-1].push_back(atof(words[2].c_str()));
            basis[basis_index].coefficients[atoi(words[0].c_str())-1].push_back(atof(words[3].c_str()));
          }
          else  
          {  basis[basis_index].types.push_back(words[1]);
            basis[basis_index].exponents[atoi(words[0].c_str())-1].push_back(atof(words[2].c_str()));
            basis[basis_index].coefficients[atoi(words[0].c_str())-1].push_back(atof(words[3].c_str()));
          }
          tmpbasis_index=tmpbasis_index+1;

        }
      }
    }
    if (words.size()==6 && words[0]=="Tag" && words[1]=="Description" && read_basis_label==0)
    { getline(inFile,line);
      for (int i=0; i< basis.size(); i++)
      {  getline(inFile,line);
        words.clear();
        split(line,space,words);
        basis[i].label=words[0];
      }

    }
  }
  // for (int i=0; i<13; i++)
  //   { for (int j=0; j<basis[0].exponents[i].size(); j++)
  //     { cout<< "hello basis "<< i<< " " <<basis[0].exponents[i][j] << endl;} 

  //    } 


  inFile.close();
  inFile.clear();


  //correspond each index of basis of atoms to the index of basis
  for (int at=0; at< atoms.size(); at++) { 
    atoms[at].basis=-1000;
    for (int bas=0; bas<basis.size(); bas++)
      if (basis[bas].label==atoms[at].name)
      {  atoms[at].basis=bas;
      } 
  }


  // setting basis spherical harmonic function labels in nwchem fashion  ie. s px py pz  dxx dxy...
  for(int at=0; at < atoms.size(); at++)
  {  
    //cout << "basis for " << at << " : " << atoms[at].basis <<endl;
    if(atoms[at].basis < 0) {
      cerr << "Didn't find basis for " << atoms[at].name << endl;
      exit(1);
    }

    for (int j=0; j< basis[atoms[at].basis].types.size(); j++)
    {
      if(basis[atoms[at].basis].types[j]=="S") basis_function_type.push_back("s");
      else if (basis[atoms[at].basis].types[j]=="P") 
      { basis_function_type.push_back("px");
        basis_function_type.push_back("py");
        basis_function_type.push_back("pz");
      }
      else if(basis[atoms[at].basis].types[j]=="6D")
      { basis_function_type.push_back("dxx");
        basis_function_type.push_back("dxy");
        basis_function_type.push_back("dxz");
        basis_function_type.push_back("dyy");
        basis_function_type.push_back("dyz");
        basis_function_type.push_back("dzz");
      }
      else if(basis[atoms[at].basis].types[j]=="10F")
      { basis_function_type.push_back("fxxx");
        basis_function_type.push_back("fxxy");
        basis_function_type.push_back("fxxz");
        basis_function_type.push_back("fxyy");
        basis_function_type.push_back("fxyz");
        basis_function_type.push_back("fxzz");
        basis_function_type.push_back("fyyy");
        basis_function_type.push_back("fyyz");
        basis_function_type.push_back("fyzz");
        basis_function_type.push_back("fzzz");

      }
      else if(basis[atoms[at].basis].types[j]=="15G")
      { basis_function_type.push_back("gxxxx");
        basis_function_type.push_back("gxxxy");
        basis_function_type.push_back("gxxxz");
        basis_function_type.push_back("gxxyy");
        basis_function_type.push_back("gxxyz");
        basis_function_type.push_back("gxxzz");
        basis_function_type.push_back("gxyyy");
        basis_function_type.push_back("gxyyz");
        basis_function_type.push_back("gxyzz");
        basis_function_type.push_back("gxzzz");
        basis_function_type.push_back("gyyyy");
        basis_function_type.push_back("gyyyz");
        basis_function_type.push_back("gyyzz");
        basis_function_type.push_back("gyzzz");
        basis_function_type.push_back("gzzzz");

      }
      else {
        cerr << "unknown basis type "<< basis[atoms[at].basis].types[j] << endl;
        exit(1);
      }
    }
  }

  //  for (int i=0; i<basis_function_type.size(); i++)
  //  { cout<<i <<"  "<< basis_function_type[i] <<endl;
  //  }    
  cout<<"Done."<<endl;
}

void read_nwchem_ecp( string & outputfilename,
    vector <Atom> & atoms,
    vector <Gaussian_pseudo_writer> & ecp,
    vector <double> & ecpremoved
    )

{ 
  ifstream inFile;
  inFile.open(outputfilename.c_str());
  if (!inFile)
  { cout<<"Couldn't open"<< outputfilename<< endl; 
    exit(1);
  }
  cout<<"Reading ECP...";
  vector <double> emptyvector;
  vector <int> emptyintvector;
  vector <string> words;
  string line;
  string space=" ";
  int ecp_read=0;
  int index=0;
  vector <string> ecplabels;

  while(getline(inFile,line))
  {  //cout<<"hello1 ecp"<<endl; 
    if (ecp_read==1) break;
    words.clear();
    // cout<<"hello2"<<endl;
    split(line,space,words);
    //cout <<"hello3" <<endl;
    //cout << line << endl;
    //cout << words[2] << endl;
    if(words.size()==5 && words[2]=="Replaces")
    {  
      // cout<<"hello ecp"<<endl;
      ecpremoved.push_back(atof(words[3].c_str()));  
      ecplabels.push_back(words[0]);
    }
    //     cout<<"hello ecp"<<endl;
    if(words.size()==3 && words[0]=="R-exponent" && words[1]=="Exponent")
    {  
      Gaussian_pseudo_writer ecptmp;
      ecptmp.label=atoms[index].name;
      int block=-1;
      int have_stu_ecp=0;
      while(getline(inFile,line))
      {  words.clear();
        split(line,space,words);
        if (words.size()>2 && words[0]=="NWChem" && words[2]=="Module") 
        { ecp.push_back(ecptmp);
          ecp_read=1;
          break;
        }
        else if (words.size()==5 && words[2] == "Replaces") 
        {  ecp.push_back(ecptmp);
          ecpremoved.push_back(atof(words[3].c_str()));
          ecplabels.push_back(words[0]);
          index=index+1;
          break;
        }

        // if using Stu ecp, then the local term is ignored, so reorder it
        else if ( words.size()==5 && words[0]=="1" && words[1]=="U-s" )
        {  have_stu_ecp=1; //cout<<"hello "<<have_stu_ecp<<endl;
          if (block==-1)
          { // cout<< block <<endl;
            ecptmp.nvalue.push_back(emptyintvector);
            ecptmp.exponents.push_back(emptyvector);
            ecptmp.coefficients.push_back(emptyvector);
            //cout << " Here " << endl;
            ecptmp.nvalue[0].push_back(0);
            ecptmp.exponents[0].push_back(1.0000);
            ecptmp.coefficients[0].push_back(0.0000);
            //cout<< ecptmp.nvalue[0].size() <<endl;
            //cout<< ecptmp.nvalue[0][0] << endl;
            ecptmp.nvalue.push_back(emptyintvector);
            ecptmp.exponents.push_back(emptyvector);
            ecptmp.coefficients.push_back(emptyvector);
            block=1;
            ecptmp.nvalue[block].push_back(atoi(words[2].c_str())-2);
            ecptmp.exponents[block].push_back(atof(words[3].c_str()));
            ecptmp.coefficients[block].push_back(atof(words[4].c_str()));
          }
          else if (block==1)
          {  ecptmp.nvalue[block].push_back(atoi(words[2].c_str())-2);
            ecptmp.exponents[block].push_back(atof(words[3].c_str()));
            ecptmp.coefficients[block].push_back(atof(words[4].c_str()));
          } 

        }

        else  if (words.size()==6 && words[0]=="1" && words[2]=="L")
        {    
          if (block==-1)
          { ecptmp.nvalue.push_back(emptyintvector);
            ecptmp.exponents.push_back(emptyvector);
            ecptmp.coefficients.push_back(emptyvector);
            block=0; //cout<<"nvaluesize"<<ecptmp.nvalue.size()<<endl;
            ecptmp.nvalue[block].push_back(atoi(words[3].c_str())-2);
            ecptmp.exponents[block].push_back(atof(words[4].c_str()));
            ecptmp.coefficients[block].push_back(atof(words[5].c_str()));
          }
          else if (block==0)
          { ecptmp.nvalue[block].push_back(atoi(words[3].c_str())-2);
            ecptmp.exponents[block].push_back(atof(words[4].c_str()));
            ecptmp.coefficients[block].push_back(atof(words[5].c_str()));

          }    

        } 
        else if (words.size()==5)
        {  if ( block==atoi(words[0].c_str())-2+have_stu_ecp )
          { ecptmp.nvalue.push_back(emptyintvector);
            ecptmp.exponents.push_back(emptyvector);
            ecptmp.coefficients.push_back(emptyvector);
            block=block+1;
            ecptmp.nvalue[block].push_back(atoi(words[2].c_str())-2);
            ecptmp.exponents[block].push_back(atof(words[3].c_str()));
            ecptmp.coefficients[block].push_back(atof(words[4].c_str()));

          }
          else if ( block==atoi(words[0].c_str())-1+have_stu_ecp )
          {  ecptmp.nvalue[block].push_back(atoi(words[2].c_str())-2);
            ecptmp.exponents[block].push_back(atof(words[3].c_str()));
            ecptmp.coefficients[block].push_back(atof(words[4].c_str()));

          }

        }



      }


    }


  }
  inFile.close();
  inFile.clear();

  for(int i=0; i< ecp.size(); i++)
  { ecp[i].label=ecplabels[i];}

  for(int at=0;at<atoms.size();at++)  
  { for(int i=0; i< ecp.size(); i++)
    {
      if(ecp[i].label==atoms[at].name)
      { 
        atoms[at].charge=atoms[at].charge-ecpremoved[i];
      }
    }
  }

  for(int ecp_index=0; ecp_index< ecp.size(); ecp_index++)
  { 
    ecp[ecp_index].exponents.push_back(ecp[ecp_index].exponents[0]);
    ecp[ecp_index].nvalue.push_back(ecp[ecp_index].nvalue[0]);
    ecp[ecp_index].coefficients.push_back(ecp[ecp_index].coefficients[0]);
    ecp[ecp_index].exponents.erase(ecp[ecp_index].exponents.begin());
    ecp[ecp_index].coefficients.erase(ecp[ecp_index].coefficients.begin());
    ecp[ecp_index].nvalue.erase(ecp[ecp_index].nvalue.begin());

  } 
  cout<<"Done."<<endl;
}
//###############################
void read_nwchem_vecs(string & outputfilename,
    string & vecfilename,
    vector <Atom> & atoms,
    Slat_wf_writer & slwriter,
    vector <Gaussian_basis_set> & basis,
    vector <string> & basis_function_type,
    vector < vector <double> > & moCoeff) 
{ 
  ifstream inFile;
  inFile.open(vecfilename.c_str());
  if(!inFile) 
  { cout<<"Couldn't open"<< vecfilename <<endl;
    exit(1);
  }
  cout<<"Reading orbitals...";
  string line;
  vector <string> words;
  string space=" ";
  vector < vector <double > >  tmpCoeff;
  vector <double> emptyvector;
  int natoms=atoms.size();
  moCoeff.clear();

  int num_basis_func;
  int num_orbital_set;
  //  cout<<"hello "<<endl; 
  while(getline(inFile,line)) 
  { words.clear();
    split(line,space,words);
    if (words.size() > 0 && words[0]=="ao" && words[1]=="basis")
    { 
      moCoeff.clear();
      //lastmo=-1;
      //totmospin=0;
      //moindex=-1;

      getline(inFile,line);
      words.clear();
      split(line,space,words);
      num_orbital_set=atoi(words[0].c_str());
      getline(inFile,line);
      words.clear(); 
      split(line,space,words);
      num_basis_func=atoi(words[0].c_str());
      getline(inFile,line);
      words.clear();
      split(line,space,words);
      int norb[num_orbital_set];
      for (int l=0; l<num_orbital_set; l++)
      {  norb[l]=atoi(words[l].c_str()); }
      int m=num_basis_func/3; 
      int n=num_basis_func%3;
      if(n!=0) { m=m+1; }
      if (num_orbital_set==1)
      {  for (int i=0; i<m; i++)
        { getline(inFile,line); }
        for (int i=0; i<m; i++)
        { getline(inFile,line); }
        for (int i=0; i<norb[0]; i++)
        { tmpCoeff.push_back(emptyvector);
          for (int j=0; j<m; j++) 
          { getline(inFile,line);
            words.clear();
            split(line,space,words);
            for (int k=0; k<words.size(); k++)
              tmpCoeff[i].push_back(atof(words[k].c_str()));   
          }
        }
      }
      if (num_orbital_set==2)
      {  slwriter.spin_dwn_start=norb[0];
        for (int i=0; i<m; i++)
        { getline(inFile,line); }
        for (int i=0; i<m; i++)
        { getline(inFile,line); }
        for (int i=0; i<norb[0]; i++)
        { tmpCoeff.push_back(emptyvector);
          for (int j=0; j<m; j++)
          { getline(inFile,line);
            words.clear();
            split(line,space,words);
            for (int k=0; k<words.size(); k++)
              tmpCoeff[i].push_back(atof(words[k].c_str()));
          }
        }

        for (int i=0; i<m; i++)
        { getline(inFile,line); }
        for (int i=0; i<m; i++)
        { getline(inFile,line); }
        for (int i=norb[0]; i<norb[0]+norb[1]; i++)
        { tmpCoeff.push_back(emptyvector);
          for (int j=0; j<m; j++)
          { getline(inFile,line);
            words.clear();
            split(line,space,words);
            for (int k=0; k<words.size(); k++)
              tmpCoeff[i].push_back(atof(words[k].c_str()));
          }
        }

      }


    }  

  }
  inFile.close();
  inFile.clear();
  //  for (int i=0; i<10; i++)
  //  {  cout <<i+1<<" "<< tmpCoeff[0][i]<<endl; 
  //  } 
  switch_moCoeff_order(outputfilename, num_basis_func, basis_function_type,tmpCoeff, moCoeff);
  cout<<"Done."<<endl;
}   

// since nwchem and gamess print orbitals corresponding to a different order of basis, so we need switch the order
void switch_moCoeff_order(string & outputfilename, int num_basis_func,
    vector <string> basis_function_type,
    vector < vector <double> > & tmpCoeff,
    vector < vector <double> > & moCoeff
    )
{  ifstream inFile;
  inFile.open(outputfilename.c_str());
  if(!inFile)
  { cout<<"Couldn't open"<< outputfilename <<endl;
    exit(1);
  }
  //   vector<string> basis_function_type;
  //   string line;
  //   vector <string> words;
  //  string space=" ";
  vector<double> emptyvector;
  //   int have_basis_function_type;
  //   while(getline(inFile,line))
  //  { words.clear();
  //   split(line,space,words);
  //  if(words.size()>0 && words[1]=="function" && words[2]=="labels")
  //     {  have_basis_function_type=1;
  //        getline(inFile,line);
  //       getline(inFile,line);
  //      getline(inFile,line);
  //     while(getline(inFile,line))
  //    { words.clear();
  //          split(line,space,words);
  //         if(words.size()==0)  break;
  //        basis_function_type.push_back(words[3]);
  //     }
  //    }
  // }     



  //  if (have_basis_function_type==0) { cout <<"Can't find any basis function type. " <<endl;
  //                                     cout <<"Did you set print option to high in the input? "<< endl; 
  //                                     exit(1);
  //                                   }
  for (int i=0; i<tmpCoeff.size(); i++)
  { moCoeff.push_back(emptyvector);
    for (int j=0; j< num_basis_func; j++)
    {  moCoeff[i].push_back(0);  }

  }
  const double pi=3.1415926535;
  for (int i=0; i<tmpCoeff.size(); i++)
  {  for (int j=0; j<num_basis_func; j++)
    {  if(basis_function_type[j]=="s")
      {  moCoeff[i][j]=tmpCoeff[i][j]*1./sqrt(4.*pi); }
      else if(basis_function_type[j]=="px" || basis_function_type[j]=="py" || basis_function_type[j]=="pz")
      {  moCoeff[i][j]=tmpCoeff[i][j]*sqrt(3.)*1./sqrt(4.*pi); }
      else if(basis_function_type[j]=="dxx")
      {  moCoeff[i][j]=tmpCoeff[i][j]*sqrt(5./(4*pi)); }
      else if(basis_function_type[j]=="dxy")
      {  moCoeff[i][j+2]=tmpCoeff[i][j]*sqrt(15./(4*pi)); }
      else if(basis_function_type[j]=="dxz") 
      {  moCoeff[i][j+2]=tmpCoeff[i][j]*sqrt(15./(4*pi)); }
      else if(basis_function_type[j]=="dyy")
      {  moCoeff[i][j-2]=tmpCoeff[i][j]*sqrt(5./(4*pi)); }
      else if(basis_function_type[j]=="dyz")
      {  moCoeff[i][j+1]=tmpCoeff[i][j]*sqrt(15./(4*pi)); }
      else if(basis_function_type[j]=="dzz")
      {  moCoeff[i][j-3]=tmpCoeff[i][j]*sqrt(5./(4*pi)); }
      else if(basis_function_type[j]=="fxxx")
      {  moCoeff[i][j]=tmpCoeff[i][j]*1./sqrt(4.*pi)*sqrt(7.);}
      else if(basis_function_type[j]=="fxxy")
      {  moCoeff[i][j+2]=tmpCoeff[i][j]*1./sqrt(4.*pi)*sqrt(35.); }
      else if(basis_function_type[j]=="fxxz")
      {  moCoeff[i][j+2]=tmpCoeff[i][j]*1./sqrt(4.*pi)*sqrt(35.); }
      else if(basis_function_type[j]=="fxyy")
      {  moCoeff[i][j+2]=tmpCoeff[i][j]*1./sqrt(4.*pi)*sqrt(35.); }
      else if(basis_function_type[j]=="fxyz")
      {  moCoeff[i][j+5]=tmpCoeff[i][j]*1./sqrt(4.*pi)*sqrt(105.); }
      else if(basis_function_type[j]=="fxzz")
      {  moCoeff[i][j+2]=tmpCoeff[i][j]*1./sqrt(4.*pi)*sqrt(35.); }
      else if(basis_function_type[j]=="fyyy")
      {  moCoeff[i][j-5]=tmpCoeff[i][j]*1./sqrt(4.*pi)*sqrt(7.); }
      else if(basis_function_type[j]=="fyyz")
      {  moCoeff[i][j-1]=tmpCoeff[i][j]*1./sqrt(4.*pi)*sqrt(35.); }
      else if(basis_function_type[j]=="fzzz")
      {  moCoeff[i][j-7]=tmpCoeff[i][j]*1./sqrt(4.*pi)*sqrt(7.); }
      else if(basis_function_type[j]=="gxxxx")
      {  moCoeff[i][j]=tmpCoeff[i][j]*1./sqrt(4.*pi)*sqrt(9.); }
      else if(basis_function_type[j]=="gxxxy")
      {  moCoeff[i][j+2]=tmpCoeff[i][j]*1./sqrt(4.*pi)*sqrt(63.); }
      else if(basis_function_type[j]=="gxxxz")
      {  moCoeff[i][j+2]=tmpCoeff[i][j]*1./sqrt(4.*pi)*sqrt(63.); }
      else if(basis_function_type[j]=="gxxyy")
      {  moCoeff[i][j+6]=tmpCoeff[i][j]*1./sqrt(4.*pi)*sqrt(105.); }
      else if(basis_function_type[j]=="gxxyz")
      {  moCoeff[i][j+8]=tmpCoeff[i][j]*1./sqrt(4.*pi)*sqrt(315.); }
      else if(basis_function_type[j]=="gxxzz")
      {  moCoeff[i][j+5]=tmpCoeff[i][j]*1./sqrt(4.*pi)*sqrt(105.); }
      else if(basis_function_type[j]=="gxyyy")
      {  moCoeff[i][j-1]=tmpCoeff[i][j]*1./sqrt(4.*pi)*sqrt(63.); }
      else if(basis_function_type[j]=="gxyyz")
      {  moCoeff[i][j+6]=tmpCoeff[i][j]*1./sqrt(4.*pi)*sqrt(315.); }
      else if(basis_function_type[j]=="gxyzz")
      {  moCoeff[i][j+6]=tmpCoeff[i][j]*1./sqrt(4.*pi)*sqrt(315.); }
      else if(basis_function_type[j]=="gxzzz")
      {  moCoeff[i][j-2]=tmpCoeff[i][j]*1./sqrt(4.*pi)*sqrt(63.); }
      else if(basis_function_type[j]=="gyyyy")
      {  moCoeff[i][j-9]=tmpCoeff[i][j]*1./sqrt(4.*pi)*sqrt(9.); }
      else if(basis_function_type[j]=="gyyyz")
      {  moCoeff[i][j-5]=tmpCoeff[i][j]*1./sqrt(4.*pi)*sqrt(63.); }
      else if(basis_function_type[j]=="gyyzz")
      {  moCoeff[i][j-1]=tmpCoeff[i][j]*1./sqrt(4.*pi)*sqrt(105.); }
      else if(basis_function_type[j]=="gyzzz")
      {  moCoeff[i][j-5]=tmpCoeff[i][j]*1./sqrt(4.*pi)*sqrt(63.); }
      else if(basis_function_type[j]=="gzzzz")
      {  moCoeff[i][j-12]=tmpCoeff[i][j]*1./sqrt(4.*pi)*sqrt(9); }
    }
  }
  //   cout << moCoeff[0][14] << endl; 

}




