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
//------------------------------------------------------------------------
//src/Localize_method.cpp
#include "Localize_method.h"
#include "conjugategrad.h"
#include "MatrixAlgebra.h"
#include "Array45.h"
#include "qmc_io.h"
#include "ulec.h"
#include "System.h"
#include "Program_options.h"
#include "Qmc_std.h"
#include "Cubic_spline3D.h"
using namespace std;

//using namespace NR;
/*!
Read the "words" from the method section in the input file
via doinput() parsing, and store section information in private
variables orbs, resolution, and minmax.
Set up MO_matrix and Sample_point objects for wavefunction from input
*/
void Localize_method::read(vector <string> words,
                       unsigned int & pos,
                       Program_options & options)
{
  pos=0;	//always start from first word
  doublevar Tres;
  doublevar Tres2;
  doublevar Tres3;
  int dummy;
  vector <string> Torbs;
  vector <string> Tminmax;
  vector <string> Lcenters;
  vector <string> orbtext;

  if(! readsection(words,pos,Torbs,"PLOTORBITALS"))
    error("Need PLOTORBITALS in METHOD section");
  orbs.Resize(Torbs.size());
  for(unsigned int i=0; i<Torbs.size(); i++)
  {
    orbs(i)=atoi(Torbs[i].c_str());
  }
  
  pos=0;
  if(! readsection(words,pos,Lcenters,"LCENTERS"))
    error("Need LCENTERS in METHOD section");
  dummy=Lcenters.size();

  cout<<"I found "<<orbs.GetSize()<< " orbitals"<<endl;
  cout<<"and "<<dummy<<" components"<<endl;

  if((dummy%3)!=0)
    error("Need more/less values in LCENTERS section");
  if(dummy>3*orbs.GetSize())
    error("More center vectors in LCENTERS then orbitals involved");
  cout <<"Reading localization centers"<<endl; 
  r_max.Resize(dummy/3,3);
  for(int i=0; i<r_max.GetDim(0); i++)
  {
    for(int j=0; j<3; j++){
      r_max(i,j)=atof(Lcenters[i*3+j].c_str());
      cout <<r_max(i,j)<<" ";
      }
    cout <<endl;
  }

  pos=0;
  if(!readsection(words, pos, orbtext, "ORBITALS")) 
    error("Need ORBITALS in PLOT section");

  pos=0;
  if(! readvalue(words,pos,Tres,"RESOLUTION"))
    error("Need RESOLUTION in METHOD section");
  resolution=Tres;
  
  pos=0;
  if(! readvalue(words,pos,Tres2,"RADIUS"))
    error("Need RADIUS in METHOD section");
  radius=Tres2;

  pos=0;
  if(! readvalue(words,pos,Tres3,"NORM"))
    error("Need NORM in METHOD section");
  delta_norm=Tres3;
  
  pos=0;
  if(! readsection(words,pos,Tminmax,"MINMAX"))
    error("Need MINMAX in METHOD section");
  if(Tminmax.size() != 6)
    error("MINMAX needs 6 values");
  minmax.Resize(6);
  for(unsigned int i=0; i<Tminmax.size(); i++)
  {
    minmax(i)=atof(Tminmax[i].c_str());
  }

  sysprop=NULL;
  allocate(options.systemtext[0], sysprop);
  Array1 <Array1 <int> > orblist(1);
  orblist(0).Resize(orbs.GetDim(0));
  
  for(int i=0; i< orbs.GetDim(0); i++)
  {
    cout << "i " << i << endl;
    cout << "orbs " << orbs(i) << endl;
    orblist(0)(i)=orbs(i)-1;
  }
  //  unsigned int startpos=0;
  //mymomat.read(orbtext,startpos,sysprop);
  allocate(orbtext, sysprop, mymomat);
  mymomat->buildLists(orblist);
  mywalker=NULL;
  sysprop->generateSample(mywalker);
  mymovals.Resize(mymomat->getNmo(),5);
 }



bool sphere_check(Array1 <doublevar> & pos,
                      doublevar r_max_x, doublevar r_max_y,  doublevar r_max_z, 
                      doublevar radius)
  //checks whether vector pos is inside the sphere 
{
  if (((pos(0)-r_max_x)*(pos(0)-r_max_x)+(pos(1)-r_max_y)*(pos(1)-r_max_y)
	  +(pos(2)-r_max_z)*(pos(2)-r_max_z))<radius*radius)
    return true;
  else 
    return false;
}

bool pbc_check(Array1 <doublevar> & pos,Array1 <doublevar> & prim_trans,
                      doublevar r_max_x, doublevar r_max_y,  doublevar r_max_z, 
                      doublevar radius,Array1 <doublevar> & new_pos)
  //constructs the periodic images for neighbouring cells
  //end checks whether are inside the sphere
{
  for(int i=0;i<=2;i++){
    if (i<2) new_pos(0)=pos(0)+i*prim_trans(0);
    else new_pos(0)=pos(0)-prim_trans(0);
    for(int j=0;j<=2;j++){
      if (j<2) new_pos(1)=pos(1)+j*prim_trans(1);
      else new_pos(1)=pos(1)-prim_trans(1);
      for(int k=0;k<=2;k++){
        if (k<2) new_pos(2)=pos(2)+k*prim_trans(2);
        else new_pos(2)=pos(2)-prim_trans(2);
        if (sphere_check(new_pos,r_max_x,r_max_y,r_max_z,radius)) return true;
      }
    }
  }
  return false;
}


doublevar norm_i(Array1 <doublevar> & a)
{
  int m=a.GetSize();
  double norm=0.0;
  for (int k=0;k<m;k++)
    norm+=a(k)*a(k);
  return norm;
}

void normalize(Array1 <doublevar> & a)
{
  int m=a.GetSize();
  double norm;
  norm=norm_i(a);
  for (int k=0;k<m;k++)
    a(k)=a(k)/sqrt(norm);
}

class Minimization_function {

private: 
  Array2<doublevar> _overlap_ma_R;

public:  
  void init_overlap_ma_R(Array2<doublevar> &overlap_ma_R);  
  doublevar func(Array1 <doublevar> &a);
  void dfunc(Array1 <doublevar> &a, Vec_O_doublevar &xi);
  Minimization_function();
  ~Minimization_function(); 
};

Minimization_function::Minimization_function(){
}

Minimization_function::~Minimization_function(){
}

void Minimization_function::init_overlap_ma_R(Array2<doublevar> &overlap_ma_R){
  assert(overlap_ma_R.GetDim(0)==overlap_ma_R.GetDim(1));
  _overlap_ma_R.Resize(overlap_ma_R.GetDim(0),overlap_ma_R.GetDim(0));
  for(int k=0;k<overlap_ma_R.GetDim(0);k++)
    for(int l=0;l<=k;l++)
      _overlap_ma_R(l,k)=_overlap_ma_R(k,l)=overlap_ma_R(k,l);
}

doublevar Minimization_function::func(Array1 <doublevar>  &a){
  // \begin{equation}
  //   {\cal F}=1-\frac{\sum_{i,j}^{N}a_i a_j S_{ij}}{\sum_{i}^{N}a_{i}^2}
  //  \end{equation}
  int m=a.GetSize();
  double tmp=0.0;
  double norm=norm_i(a);
  for (int k=0;k<m;k++)
    for (int l=0;l<m;l++){
      tmp+=a(k)*a(l)*_overlap_ma_R(k,l);
    }
  //   cout << 1.0 -tmp/norm<<endl;
  return 1.0 -tmp/norm;
}


void Minimization_function::dfunc(Array1 <doublevar> &a, Array1 <doublevar> &xi){
  // \begin{equation}
  //    \frac{\partial {\cal F}}{\partial a_k}=
  //         -2\frac{\sum_{i}^{N}a_iS_{ik}+2a_k({\cal F}-1)}{\sum_{i}^{N}a_{i}^2}
  //  \end{equation} 
  int m=a.GetSize();
  double fa_i=func(a);
  double norm=norm_i(a);
  double temp=0.0;
  for (int j=0;j<m;j++){
    temp=0.0;
    for (int k=0;k<m;k++)
      temp+=a(k)*_overlap_ma_R(j,k);
    xi(j)=-2.0*(temp+a(j)*(fa_i-1.0))/norm;
    //    cout <<xi(j)<<"  ";
  }
  //  cout <<endl;
} 

/*
obsolete
double local_orbital(int j, Array1 <doublevar> &a, Array2<doublevar> & grid)
{
  double tmp=0.0;
  
  for (int k=0;k<a.GetSize();k++)
    tmp+=a(k)*grid(k,j);
      
  return tmp;
}
*/


void Localize_method::run(Program_options & options, ostream & output)
{
  ofstream os; //for writing to *.plt files
  Array1 <ofstream> os_array;
  unsigned int electron=0; //# of the electron that will move through grid
  string pltfile; //name of plotfile being written
  string orbfile; //name of orbfile being written
  Array1 <string> orbfile_array;
  Array1 <doublevar> xyz(3),resolution_array(3); //position of electron "in" MO
  Array1 <int> D_array1(3); //dummy array1
  D_array1=0; //sets all 3 components to 0. use as counter for gridpoints

  const double ftol=0.00001;//precision in minimalization routine
  double fret; //values of minima
  int iter=0; //number of iterations
  
  Array1 <doublevar> p0(orbs.GetSize()); //starting vector for minimization
  Array1 <doublevar> xi(orbs.GetSize()); //gradient of minimization
  Array2 <doublevar> RM(orbs.GetSize(),orbs.GetSize()); //rotation_matrix
  Array3 <doublevar> overlap_ma_R_all(r_max.GetDim(0),orbs.GetSize(),orbs.GetSize());
  //all overlaps
  Array1 <int> temp_orbs(orbs.GetSize()); //same like orbs-1
  Array1 <doublevar> orb_norm(orbs.GetSize()); //normalization of orbitals
  Array1 <doublevar> orb_norm_new(r_max.GetDim(0)); //normalization of new orbitals
  doublevar determinant_RM; //determinat of rotation_matrix
  doublevar backup,dummy; //some temporary variables
  Array1 <doublevar> sphere_norm(r_max.GetDim(0)); //
  Array1 <doublevar> new_pos(3),prim_trans(3);
  //zeroing
  overlap_ma_R_all=0.0;
  orb_norm=0.0;
  sphere_norm=0.0;
  RM=0.0;
  

  D_array1(0)=roundoff((minmax(1)-minmax(0))/resolution);
  D_array1(1)=roundoff((minmax(3)-minmax(2))/resolution);
  D_array1(2)=roundoff((minmax(5)-minmax(4))/resolution);

  resolution_array(0)=(minmax(1)-minmax(0))/(D_array1(0)-1);  
  resolution_array(1)=(minmax(3)-minmax(2))/(D_array1(1)-1);
  resolution_array(2)=(minmax(5)-minmax(4))/(D_array1(2)-1);

  prim_trans(0)=minmax(1)-minmax(0);
  prim_trans(1)=minmax(3)-minmax(2);
  prim_trans(2)=minmax(5)-minmax(4);


   if(orbs.GetSize()<=0)
    error("number of orbitals requested is not a positive number");
  
  //calculate value of each molecular orbital at each grid point 
  // grid values with x=fastest running variable, and z=slowest

  cout<<"calculating "<<D_array1(0)*D_array1(1)*D_array1(2)
  <<" grid points"<<endl;
  cout<<"for "<< orbs.GetDim(0) <<" old molecular orbitals"<<endl;
  

  cout <<"Should go to "<<D_array1(2)<<" of #"<<endl;
  for(int zz=0;zz<D_array1(2);zz++)
  {
    xyz(2)=minmax(4)+zz*resolution_array(2); //move forward on z axis one resolution unit
    cout <<"#";
    for(int yy=0; yy<D_array1(1);yy++)
    {
      xyz(1)=minmax(2)+yy*resolution_array(1); //move forward on y axis one resolution unit
      for(int xx=0; xx<D_array1(0);xx++)
      {
        xyz(0)=minmax(0)+xx*resolution_array(0); //move forward on x axis one resolution unit
        mywalker->setElectronPos(electron,xyz); //move elec#1 to point specified by xyz
        mymovals=-1e99; //MO_matrix::updateVal() ADDS to mymovals instead to replace
        //mymomat->updateLap(mywalker,electron,0,mymovals); //recalculate MO value for elec#1
	mymomat->updateVal(mywalker,electron,0,mymovals); //recalculate MO value for elec#1
        for(int i=0; i<orbs.GetSize(); i++){
          orb_norm(i)+=(mymovals(i,0)*mymovals(i,0))*resolution_array(0)*resolution_array(1)
            *resolution_array(2);
        }
        for(int i=0; i<r_max.GetDim(0); i++){
          if(pbc_check(xyz,prim_trans,r_max(i,0),r_max(i,1),r_max(i,2),radius,new_pos)){
            sphere_norm(i)+=0.75/(pi*radius*radius*radius);
            for(int k=0;k<orbs.GetSize();k++)
              for(int l=0;l<=k;l++){
                overlap_ma_R_all(i,k,l)+=mymovals(k,0)*mymovals(l,0);
                if (k==l) RM(k,l)=1.0; //Init of Rot. Mat. as identity
              }
          }
        }
      }
    }
  }
  cout <<endl;

  //correction for different normalization
  for(int i=0; i<r_max.GetDim(0); i++){
     for(int k=0; k<orbs.GetSize(); k++)
       for(int l=0;l<=k;l++){
         overlap_ma_R_all(i,k,l)*=resolution_array(0)*resolution_array(1)*resolution_array(2)/
           sqrt(orb_norm(k)*orb_norm(l));
         overlap_ma_R_all(i,l,k)=overlap_ma_R_all(i,k,l);
       }
     cout <<"Estimated error on sphere integral at "<<i+1<<"-th center "<<1.0-sphere_norm(i)
       *resolution_array(0)*resolution_array(1)*resolution_array(2)<<endl;
   }
  

  //printing  renormalization   
  for(int i=0; i<orbs.GetSize(); i++){
    int ii=orbs(i);
    cout <<"total overlap of "<<ii<<"-th orbital "<< orb_norm(i) <<endl;
  }
  cout <<endl;

 
  //starting localization
  cout<<"Starting the localization section with "<<orbs.GetSize()<<" orbitals"<<endl;
  cout <<"number of grid points for z, y, & x direction"<<endl;
  cout <<D_array1(2)<<" "<<D_array1(1)<<" "<<D_array1(0)<<endl;
  cout <<minmax(4)<<" "<<minmax(5)<<" "<<minmax(2)<<" "<<minmax(3)
       <<" "<<minmax(0)<<" "<<minmax(1)<<endl;
  cout <<endl;
    
  for(int i=0; i<r_max.GetDim(0); i++) {
    Minimization_function minfunc; //initialization of object
    Array2 <doublevar> overlap_ma_R(overlap_ma_R_all(i)); 
    minfunc.init_overlap_ma_R(overlap_ma_R);
    fret=1.0;
    int tries=0;
    while (fret > delta_norm){

      //initialization of p0
      for (int ii=0;ii<orbs.GetSize();ii++){
       	p0(ii)=1.0-2.0*unif();
      }
      normalize(p0);

      cout <<"Performing maximization of sum of overlaps at "<<i+1<<"-th center"<<endl;
      iter=0;
      NR_conjugate_gradient <Minimization_function> minimizer(& minfunc);
      minimizer.frprmn(p0, ftol, iter, fret);
      if (fret > delta_norm ) //testing whether better than required limit on normalization
	cout  <<"iter:"<<iter<< "  value=" << fret <<" retrying ...\n";
      if (tries>200) 
        error("Couldn't find value small enough, increase R_cut");
      tries++;
    }
    backup=norm_i(p0);
    cout <<"Current norm is "<<backup<<endl;
    normalize(p0); //normalization of result
    //printout of values of p0
    for (int j=0; j<orbs.GetSize() ; j++){
      RM(i,j)=p0(j)/sqrt(orb_norm(j));
      cout <<"p0["<<j<<"]= "<< p0(j)<<endl;
    }
    cout  <<"@ For "<<i+1<<" center, the overlap outside " << fret <<" after "
          <<iter<<" iterations\n";

  }

  cout <<"Writing it to plot files" <<endl;
  os_array.Resize(r_max.GetDim(0));
  orbfile_array.Resize(r_max.GetDim(0));

  for(int i=0; i<r_max.GetDim(0); i++){
    char strbuff[40];
    sprintf(strbuff, "%d", i+1);
    orbfile_array(i) = options.runid;
    orbfile_array(i) += ".orb";
    orbfile_array(i) += strbuff;
    orbfile_array(i) += ".plt"; /*FIGURE OUT HOW TO CONVERT INT TO STRING*/

    os_array(i).open(orbfile_array(i).c_str());
    cout<<"writing to "<<orbfile_array(i)<<endl;
    os_array(i)<<"3 "; //rank=3 always
    os_array(i)<<"2\n"; //dummy variable => "Orbital/density surface"
    //number of grid points for x, y, & z direction
    os_array(i) <<D_array1(2)<<" "<<D_array1(1)<<" "<<D_array1(0)<<endl;
    os_array(i) <<minmax(4)<<" "<<minmax(5)<<" "<<minmax(2)<<" "<<minmax(3)
                <<" "<<minmax(0)<<" "<<minmax(1)<<endl;
  }

  orb_norm_new=0.0;
  cout <<"Should go to "<<D_array1(2)<<" of #"<<endl;
  for(int zz=0;zz<D_array1(2);zz++){
    xyz(2)=minmax(4)+zz*resolution_array(2); //move forward on z axis one resolution unit
    cout <<"#";
    for(int yy=0; yy<D_array1(1);yy++){
      xyz(1)=minmax(2)+yy*resolution_array(1); //move forward on y axis one resolution unit
      for(int xx=0; xx<D_array1(0);xx++){
        xyz(0)=minmax(0)+xx*resolution_array(0); //move forward on x axis one resolution unit
        mywalker->setElectronPos(electron,xyz); //move elec#1 to point specified by xyz
        mymovals=-1e99; //MO_matrix::updateVal() ADDS to mymovals instead to replace
        //mymomat->updateLap(mywalker,electron,0,mymovals); //recalculate MO value for elec#1
	mymomat->updateVal(mywalker,electron,0,mymovals); //recalculate MO value for elec#1
        for(int i=0; i<r_max.GetDim(0); i++){
          /* http://www.csc.fi/gopenmol/developers/plt_format.phtml */
          dummy=0.0;
          for(int j=0; j<orbs.GetSize(); j++)
              dummy+=RM(i,j)*mymovals(j,0);
          os_array(i)<<dummy<<endl; //writting orbital
          orb_norm_new(i)+=dummy*dummy*resolution_array(0)*resolution_array(1)*resolution_array(2);
        }
      }
    }
  }
  cout <<endl;
  for(int i=0; i<r_max.GetDim(0); i++){
    cout <<"New total overlap of  "<<i+1<<"-th orbital "<< orb_norm_new(i) <<endl;
    os_array(i).close();
  }
  
  determinant_RM=Determinant(RM,orbs.GetSize()); 
  cout << "determinant of RM= "<<determinant_RM<<endl;
  if (abs(determinant_RM)<1.0e-15) error("Determinant of rotation < 1.0e-15");
  
  //print new orbitals to file
  orbfile = options.runid;
  orbfile += "_loc.orb";
  for (int k=0;k<orbs.GetSize();k++)
    temp_orbs(k)=orbs(k)-1;
   
  ofstream testorb(orbfile.c_str());
  mymomat->writeorb(testorb,  RM, temp_orbs);
  testorb.close();
}


/*!
Print information about private variables {orbs,resolution,minmax}
*/
int Localize_method::showinfo(ostream & os)
{
  os << " MO_matrix: ";
  mymomat->showinfo(os);
  os << endl;

  os<<"#############Localize_method#################\n";
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
