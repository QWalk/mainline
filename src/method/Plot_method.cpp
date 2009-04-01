/*
 
Copyright (C) 2007 Zachary Helms
 with further modifications by Lucas K. Wagner

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

#include "Program_options.h"
#include "Plot_method.h"
#include "qmc_io.h"
#include "System.h"
#include "MatrixAlgebra.h"
/*!
Read the "words" from the method section in the input file
via doinput() parsing, and store section information in private
variables orbs, resolution, and minmax.
Set up MO_matrix and Sample_point objects for wavefunction from input
*/
void Plot_method::read(vector <string> words,
                       unsigned int & pos,
                       Program_options & options)
{
  pos=0;	//always start from first word
  vector <string> Torbs;
  vector <string> Tminmax;
  vector <string> orbtext;


  use_complex=0;
  if(!readsection(words, pos=0, orbtext, "ORBITALS")){
    use_complex=1;
    if(!readsection(words, pos=0, orbtext, "CORBITALS"))
      error("Need ORBITALS/CORBITALS in PLOT section");
  }

  if(! readvalue(words,pos=0,resolution,"RESOLUTION"))
    resolution=.2;

 

  if(haskeyword(words, pos=0, "JEEP_CUBE") )
    jeep_like_cube_file=1;
  else  jeep_like_cube_file=0; 

  if(haskeyword(words, pos=0, "PERIODIC") )
    periodic=1;
  else periodic=0; 

  sysprop=NULL;
  allocate(options.systemtext[0],  sysprop);

  if(!use_complex)
    allocate(orbtext, sysprop, mymomat);
  else
    allocate(orbtext, sysprop, cmymomat);

  if(!use_complex)
    totmo=mymomat->getNmo();
  else
    totmo=cmymomat->getNmo();
  int maxorb=0;
  if(readsection(words,pos=0,Torbs,"PLOTORBITALS")) {
    orbs.Resize(Torbs.size());
    for(unsigned int i=0; i<Torbs.size(); i++) {
      orbs(i)=atoi(Torbs[i].c_str());
      if(orbs(i)> maxorb) maxorb=orbs(i);
    }
  }
  else {
    orbs.Resize(totmo);
    for(int i=0; i< totmo; i++) {
      orbs(i)=i+1;
    }
  }

  
  if(maxorb > totmo) {
    error("Too high orbital requested in PLOTORBITALS.  Try increasing " 
	  "NORB in the ORBITALS/CORBITALS section.");
  }

  Array1 <Array1 <int> > orblist(1);
  orblist_pernode.Resize(1);
  int norbs_pernode;
  norbs_pernode=int(orbs.GetDim(0)/mpi_info.nprocs);
  if(norbs_pernode==0)
    norbs_pernode=1;
  

  orblist(0).Resize(orbs.GetDim(0));
  int counter=0;
  for(int i=0; i< orbs.GetDim(0); i++) {
    //cout << "i " << i << endl;
    //cout << "orbs " << orbs(i) << endl;
    orblist(0)(i)=orbs(i)-1;
    if((i>=norbs_pernode*mpi_info.node && i<norbs_pernode*(mpi_info.node+1)) || 
       (i==norbs_pernode*mpi_info.nprocs+mpi_info.node)  )
      counter++;
  }
  orblist_pernode(0).Resize(counter);
  counter=0;
  for(int i=0; i< orbs.GetDim(0); i++) {
    if((i>=norbs_pernode*mpi_info.node && i<norbs_pernode*(mpi_info.node+1)) ||
       (i==norbs_pernode*mpi_info.nprocs+mpi_info.node)  )
      orblist_pernode(0)(counter++)=orblist(0)(i);
  }
  cout <<" on node: "<<mpi_info.node<<" plotting these: ";
  for(int i=0;i<orblist_pernode(0).GetSize();i++)
    cout <<orblist_pernode(0)(i)+1<<" ";
  cout<<endl;


  if(!use_complex)
    mymomat->buildLists(orblist_pernode);
  else
    cmymomat->buildLists(orblist_pernode);

  mywalker=NULL;
  sysprop->generateSample(mywalker);

  if(!use_complex)
    mymovals.Resize(totmo,5);
  else
    cmymovals.Resize(totmo,5);

  periodic=0; 
  origin.Resize(3);
  origin=0;
  LatticeVec.Resize(3,3);
   

  minmax.Resize(6);


  if(readsection(words,pos=0,Tminmax,"MINMAX")) { 
    if(Tminmax.size() != 6)
      error("MINMAX needs 6 values");
    for(unsigned int i=0; i<Tminmax.size(); i++) {
      minmax(i)=atof(Tminmax[i].c_str());
    }
  }
  else if(sysprop->getBounds(LatticeVec)){
    periodic=1;
    single_write(cout,"Using periodic boundaries\n");
    sysprop->kpoint(kpoint);
    sysprop->getorigin(origin);
    //these is periodic MINMAX
    for(int d=0;d<3;d++){
      minmax(2*d)=origin(d);
      minmax(2*d+1)=LatticeVec(0,d)+LatticeVec(1,d)+LatticeVec(2,d)+origin(d);
    }
  }
  else  { error("Need MINMAX in METHOD section");} 

}

/*!
 
*/
void Plot_method::run(Program_options & options, ostream & output) {
  ofstream os; //for writing to *.plt files
  unsigned int electron=0; //# of the electron that will move through grid
  string pltfile; //name of plotfile being written
  Array1 <doublevar> xyz(3);//position of electron "in" MO
    Array2 <doublevar> resolution_array(3,3); 
    Array1 <int> D_array1(3); //dummy array1
    D_array1=0; //sets all 3 components to 0. use as counter for gridpoints
    resolution_array=0;
    if(periodic){
      for(int i=0;i<3;i++){
        doublevar lenght=0;
        for(int j=0;j<3;j++){
          if(mpi_info.node==0)
            cout <<LatticeVec(i,j)<<" ";
          lenght+=LatticeVec(i,j)*LatticeVec(i,j);
        }
        if(mpi_info.node==0)
          cout <<endl;
        lenght=sqrt(lenght);
        if(mpi_info.node==0)
          cout <<"lenght "<<lenght<<endl;
        D_array1(i)= roundoff(lenght/resolution); 
        for(int j=0;j<3;j++){
          resolution_array(i,j)=LatticeVec(i,j)/D_array1(i);
        }
      }
      
    }  
  else{
    for(int d=0;d<3;d++)
      D_array1(d)=roundoff((minmax(2*d+1)-minmax(2*d))/resolution);
    for(int d=0;d<3;d++)
      resolution_array(d,d)=(minmax(2*d+1)-minmax(2*d))/(D_array1(d)-1);
  }

  int npts=D_array1(0)*D_array1(1)*D_array1(2);
  
  Array3 <doublevar> grid;
  if(!use_complex)
    grid.Resize(1,orblist_pernode(0).GetSize(),npts);
  else
    grid.Resize(2,orblist_pernode(0).GetSize(),npts);

  Array1 <doublevar> density(npts);
  
  //generate .xyz file for gOpenMol to view coordinates
  if(mpi_info.node==0){
    pltfile=options.runid + ".xyz";
    os.open(pltfile.c_str());
    write_xyz(sysprop,os);
    os.close();
  }


  //calculate value of each molecular orbital at each grid point and store in an Array1
  // grid values with x=fastest running variable, and z=slowest
  double MBs;
  if(!use_complex)
    MBs=npts*8/(1024*1024);
  else
    MBs=2.0*npts*8/(1024*1024); 
  if(mpi_info.node==0){
    cout<<"calculating "<<D_array1(0)<<"x"<<D_array1(1)<<"x"<<D_array1(2)<<" = "<<npts<<" grid points\n";
    if(!use_complex)
      cout<<"for "<<orbs.GetDim(0)<<" real orbitals\n";
    else
      cout<<"for "<<orbs.GetDim(0)<<" complex orbitals\n";
    cout<<"using "<<MBs<<" Mb per orbital in memory or "<<MBs*orbs.GetDim(0)<<" Mb in total memory\n";
  }
  
  int count=0;
  xyz=0;
  for(int xx=0;xx<D_array1(0);xx++){
    doublevar maxatborder=0;
    if(!periodic){
      xyz(0)=minmax(0)+xx*resolution_array(0,0); //move forward on z axis one resolution unit
      if(mpi_info.node==0){
        cout << "x " << xyz(0) << endl;
        flush(cout);
      }
    }
    if(mpi_info.node==0){
      cout << 100*xx/(D_array1(0)-1) <<" % of plot"<<endl;
      flush(cout);
    }
    for(int yy=0; yy<D_array1(1);yy++){
      if(!periodic)
        xyz(1)=minmax(2)+yy*resolution_array(1,1);  
      for(int zz=0; zz<D_array1(2);zz++){
        if(periodic){
          for(int i=0;i<3;i++)
            xyz(i)=origin(i)+xx*resolution_array(0,i)+yy*resolution_array(1,i)+zz*resolution_array(2,i);
          //cout <<xyz(0)<<" "<<xyz(1)<<" "<<xyz(2)<<endl;
        }
        else
          xyz(2)=minmax(4)+zz*resolution_array(2,2);
        mywalker->setElectronPos(electron,xyz); //move elec#1 to point specified by xyz
        if(!use_complex)
          mymomat->updateVal(mywalker,electron,0,mymovals); //recalculate MO value for elec#1
        else
          cmymomat->updateVal(mywalker,electron,0,cmymovals); //recalculate MO value for elec#1
        density(count)=0;
        for(int i=0; i<orblist_pernode(0).GetSize(); i++) {
          if(!use_complex){
            grid(0,i,count)=mymovals(i,0);
            density(count)+=mymovals(i,0)*mymovals(i,0);
          }
          else{
            grid(0,i,count)=cmymovals(i,0).real();
            grid(1,i,count)=cmymovals(i,0).imag();
            density(count)+=cmymovals(i,0).real()*cmymovals(i,0).real()+
              cmymovals(i,0).imag()*cmymovals(i,0).imag();
          }
        }
        if(!periodic){
          if(zz==D_array1(2)-1 || yy==D_array1(1)-1 ||  xx==D_array1(1)-1){
            for(int i=0; i<orblist_pernode(0).GetSize(); i++) {
              if(!use_complex){
                if(fabs(mymovals(i,0))>maxatborder )
                  maxatborder=fabs(mymovals(i,0));
              }
              else{
                if(cabs(cmymovals(i,0))>maxatborder )
                  maxatborder=cabs(cmymovals(i,0));
              }
              
            }
          }
        }
        count++; //index for cycling through grid points
      }
    }
    if(!periodic){
      cout <<"Max value at the boundary "<<maxatborder<<endl;
    }
  }//xx
  
  //Loop through and generate plot files for each orbital requested
  if(orblist_pernode(0).GetSize()<=0)
    error("number of orbitals requested is not a positive number");
  cout<<"saving data for "<<orbs.GetSize()<<" molecular orbitals"<<endl;
  
  
  for(int i=0; i<orblist_pernode(0).GetSize(); i++) {
    //output to file with orbital number in it
    string basename,basename2;
    char strbuff[40];
    sprintf(strbuff, "%d", orblist_pernode(0)(i)+1);
    basename2 = options.runid;
    basename2 += ".orb";
    basename2 += strbuff;
    
    for(int part=0;part<grid.GetDim(0);part++){
      if(use_complex){
        if(part==0)
          basename=basename2+".real";
        else
          basename=basename2+".imag";
      }
      else{
        basename=basename2;
      }
      
      int dopltfile=0;
      int cubefile=1;
      //Note that the pltfile will be rotated, since it requires z,y,x, 
      //but cube requires x,y,z.
      if(dopltfile) {
        pltfile = basename+".plt"; 
        os.open(pltfile.c_str());
        cout<<"writing to "<<pltfile<<endl;
        // http://www.csc.fi/gopenmol/developers/plt_format.phtml
        os<<"3 "; //rank=3 always
        os<<"2\n"; //dummy variable => "Orbital/density surface"
                   //number of grid points for x, y, & z direction
        os <<D_array1(2)<<" "<<D_array1(1)<<" "<<D_array1(0)<<endl;
        os <<minmax(4)<<" "<<minmax(5)<<" "<<minmax(2)<<" "<<minmax(3)
          <<" "<<minmax(0)<<" "<<minmax(1)<<endl;
        
        os.setf(ios::scientific);
        for(int j=0; j<(D_array1(0)*D_array1(1)*D_array1(2)); j++) {
          os<<setw(16)<<setprecision(8)<<grid(part,i,j)<<endl;
        }
        os.unsetf(ios::scientific);
        os<<setprecision(6);
        os.close();
      }
      if(cubefile) {
        string cubename=basename+".cube";
        cout <<" node: "<<mpi_info.node<<" is storing orbital: "<<orblist_pernode(0)(i)+1<<" to file: "<<cubename<<endl;
        os.open(cubename.c_str());
        int natoms=sysprop->nIons();
        if(!jeep_like_cube_file){
          os << "QWalk plot output\n";
          os << "Molecular orbital " << orblist_pernode(0)(i)+1 << endl;	
          os << "  " << natoms << "   " << minmax(0) << "   "
            << minmax(2) << "   " << minmax(4) << endl;
          for(int i=0;i<3;i++)
            os << D_array1(i) << "   " << resolution_array(i,0) <<"  "<<resolution_array(i,1)<<"  "<<resolution_array(i,2)<< endl;
          Array1 <doublevar> pos(3);
          for(int at=0; at< natoms; at++) {
            mywalker->getIonPos(at,pos);
            os << "   " << mywalker->getIonCharge(at) << "   0.0000    " << pos(0) 
              <<"    " << pos(1) << "   " << pos(2) << endl;
          }
        }
        //optional output for mesh orbital
        if(jeep_like_cube_file){
          os << "QWalk plot output\n";
          os << "Molecular orbital " << orblist_pernode(0)(i)+1 << endl;
          for(int i=0;i<3;i++)
            os << D_array1(i) << "   " << resolution_array(i,0) <<"  "<<resolution_array(i,1)<<"  "<<resolution_array(i,2)<< endl;
          for(int k=0;k<3;k++) //need 3 empty lines to mimic jeep like file
            os << endl;
          if(natoms){
            os << "  " << natoms<<endl;
            Array1 <doublevar> pos(3);
            for(int at=0; at< natoms; at++) {
              mywalker->getIonPos(at,pos);
              os << "   " << mywalker->getIonCharge(at) << "   0.0000    " << pos(0) 
                <<"    " << pos(1) << "   " << pos(2) << endl;
            }
          }
          else //if natoms=0 like in HEG
            os << "  1" << endl<<endl;
          os << "  1" << endl;
          os << "  "<< minmax(0)<<"  "<< minmax(2)<<"  "<< minmax(4)<<endl;
          os << "  "<< minmax(1)-minmax(0)<<"  "<< minmax(3)-minmax(2)<<"  "<< minmax(5)-minmax(4)<<endl;
          os <<endl<<endl;
          os << "  "<<D_array1(0)<<"  "<<D_array1(1)<<"  "<<D_array1(2)<<endl;
        }
        os.setf(ios::scientific);
        for(int j=0; j< npts; j++) {
          os <<setw(20)<<setprecision(10)<<grid(part,i,j);
          if(j%6 ==5) os << endl;
        }
        os << endl;
        os.unsetf(ios::scientific);
        os<<setprecision(6);
        os.close();
        
      }
    }//part
  }//orbital
  
  
  //add density from each node
  parallel_sum(density);
  if(mpi_info.node==0){
    string cubename=options.runid+".dens.cube";
    os.open(cubename.c_str());
    os << "QWalk plot output\n";
    os << "Electron density" << endl;
    int natoms=sysprop->nIons();
    os << "  " << natoms << "   " << minmax(0) << "   "
      << minmax(2) << "   " << minmax(4) << endl;
    for(int i=0;i<3;i++)
      os << D_array1(i) << "   " << resolution_array(i,0) <<"  "<<resolution_array(i,1)<<"  "<<resolution_array(i,2)<< endl;
    Array1 <doublevar> pos(3);
    for(int at=0; at< natoms; at++) {
      mywalker->getIonPos(at,pos);
      os << "   " << mywalker->getIonCharge(at) << "   0.0000    " << pos(0) 
        <<"    " << pos(1) << "   " << pos(2) << endl;
    }
    
    os.setf(ios::scientific);
    for(int j=0; j< npts; j++) {
      os <<setw(20)<<setprecision(10)<<density(j);
      if(j%6 ==5) os << endl;
    }
    os << endl;
    os.unsetf(ios::scientific);
    os<<setprecision(6);
    os.close();    
  }

#ifdef USE_MPI
  MPI_Barrier(MPI_Comm_grp);
#endif


}


/*!
Print information about private variables {orbs,resolution,minmax}
*/
int Plot_method::showinfo(ostream & os)
{
  os<<"#############Plot_method#################\n";
  os<<"orbs="<<orbs(0);
  for(int i=1;i<orbs.GetSize();i++)
    os<<", "<<orbs(i);
  os<<endl;
  os<<"resolution="<<resolution<<endl;
  os<<"xmin="<<minmax(0)<<" xmax="<<minmax(1)<<endl;
  os<<"ymin="<<minmax(2)<<" ymax="<<minmax(3)<<endl;
  os<<"zmin="<<minmax(4)<<" zmax="<<minmax(5)<<endl;
  os<<"done"<<endl;
  return 1;
}
