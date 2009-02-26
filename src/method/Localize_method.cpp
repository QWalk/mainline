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
#include "Program_options.h"
#include "Localize_method.h"
#include "MatrixAlgebra.h"
#include "Array45.h"
#include "qmc_io.h"
#include "ulec.h"
#include "System.h"
#include "Qmc_std.h"
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
  vector <string> Torbs;
  vector <string> Tminmax;
  vector <string> Lcenters;
  vector <string> orbtext;

  pos=0;
  
  sysprop=NULL;
  allocate(options.systemtext[0], sysprop);

  if(!readsection(words, pos, orbtext, "ORBITALS")) 
    error("Need ORBITALS in PLOT section");
  //allocate orbitals
  allocate(orbtext, sysprop, mymomat);
  nmo=mymomat->getNmo();
  
  pos=0;
  if(readsection(words,pos,Torbs,"ORBITALS_TO_LOCALIZE")){
    orbs.Resize(Torbs.size());
    for(unsigned int i=0; i<Torbs.size(); i++){
      orbs(i)=atoi(Torbs[i].c_str())-1;
    }
    norbs=orbs.GetSize();
    cout<<"I found "<<norbs<< " orbitals"<<endl;
  }
  else{
    orbs.Resize(nmo);
    for(unsigned int i=0; i<nmo; i++)
      orbs(i)=i;
    norbs=orbs.GetSize();
    cout<<"Using all NMO orbitals for localization"<<endl;
  }
  
  pos=0;
  if(! readsection(words,pos,Lcenters,"LCENTERS"))
    error("Need LCENTERS in METHOD section");
  int dummy=Lcenters.size();
  if((dummy%3)!=0)
    error("Need more/less values in LCENTERS section");
  if(dummy>3*norbs)
    error("More center vectors in LCENTERS then orbitals involved");
  
  
  center.Resize(dummy/3);
  ncenters=center.GetDim(0);
  cout <<"Reading "<<ncenters<<" localization centers"<<endl; 
  
  for(int i=0; i<ncenters; i++)
  {
    center(i).Resize(3);
    for(int j=0; j<3; j++){
      center(i)(j)=atof(Lcenters[i*3+j].c_str());
      cout <<center(i)(j)<<" ";
      }
    cout <<endl;
  }

  
  

  pos=0;
  if(! readvalue(words,pos,Tres,"RESOLUTION"))
    error("Need RESOLUTION in METHOD section");
  resolution=Tres;
  
  pos=0;
  if(! readvalue(words,pos,Tres2,"RADIUS"))
    error("Need RADIUS in METHOD section");
  radius=Tres2;

  //read in system an latice vectors  
  int ndim=3;
  latVec.Resize(ndim, ndim);
  origin.Resize(3);
  if(!sysprop->getBounds(latVec)){
    cout <<" Not a periodic system, will need latice vectors"<<endl;
    pos=0;
    vector <string> latvectxt;
    if(!readsection(words, pos, latvectxt, "LATTICEVEC"))
      error("LATTICEVEC is required in PERIODIC");
    if(latvectxt.size() != ndim*ndim)
      error("LATTICEVEC must have exactly ",ndim*ndim, " values");
    for(int i=0; i< ndim; i++) {
      for(int j=0; j< ndim; j++) {
        latVec(i,j)=atof(latvectxt[i*ndim+j].c_str());
      }
    }
    vector <string> origintxt;
    if(readsection(words, pos=0, origintxt, "ORIGIN")) {
      if(origintxt.size() < 3) error("ORIGIN section must have at least 3 elements.");
      for(int i=0; i< 3; i++) origin(i)=atof(origintxt[i].c_str());
    }
    else {
      origin=0;   //defaulting the origin to zero
    }
  }
  else{
    cout <<"Periodic system: Using its latice vectors"<<endl;
    for(int i=0;i<ndim;i++){
      for(int j=0;j<ndim;j++)
        cout<<latVec(i,j)<<"  ";
      cout <<endl;
    }
    cout <<"Warning: Origin is set to (0, 0, 0)"<<endl;
    origin=0;
  }
  
 

  
  if(orbs.GetDim(0)>nmo)
    error("ORBITALS_TO_LOCALIZE > NMO");
  all_orbs.Resize(nmo);
  for(int i=0;i<nmo;i++){
    all_orbs(i)=i;
  }
  mymovals.Resize(nmo,5);
  //make list
  Array1 <Array1 <int> > orblist(1);
  orblist(0).Resize(orbs.GetDim(0));
  for(int i=0; i< orbs.GetDim(0); i++)
    orblist(0)(i)=orbs(i);
  mymomat->buildLists(orblist);

  mywalker=NULL;
  sysprop->generateSample(mywalker);
  //cout <<" done read in"<<endl;

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


doublevar volume3D(Array2 <doublevar> & latVec ){
  doublevar volume=fabs(latVec(2,0)*(latVec(0,1)*latVec(1,2)-latVec(0,2)*latVec(1,1))+
                        latVec(2,1)*(latVec(0,2)*latVec(1,0)-latVec(0,0)*latVec(1,2))+
                        latVec(2,2)*(latVec(0,0)*latVec(1,1)-latVec(0,1)*latVec(1,0)));
  return volume;
}

doublevar distance(Array1 <doublevar> & a){
  doublevar tmp=0.0;
  for(int i=0;i<a.GetSize();i++)
    tmp+=a(i)*a(i);
  return sqrt(tmp);
}

doublevar distance_between_a_b(Array1 <doublevar> & a, Array1 <doublevar> & b){
  assert(a.GetSize()==b.GetSize());
  Array1 <doublevar> c(a.GetSize());
  for(int i=0;i<a.GetSize();i++)
    c(i)=a(i)-b(i);
  return distance(c);
}



/*
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
  //cout <<"norm "<<norm<<endl;
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

*/




void check_the_eignevalue_solver(Array2 <doublevar> & A, 
                            Array2 <doublevar> & B, 
                            Array1 <doublevar> & eigenvals, 
                            Array2 <doublevar> & eigenvecs, 
                            int index){
  int dim=A.GetDim(0);
  doublevar left,right,ratio;
  left=right=0;
  for(int m=0;m<dim;m++)
    for(int n=0;n<dim;n++){
      left+=A(m,n)*eigenvecs(n,index);
      right+=B(m,n)*eigenvecs(n,index);
    }
  if(right!=0.0)
    ratio=left/right;
  else
    error("Right hand side is not well defined");
  
  //cout <<" exact value "<<ratio<<" calculated eigenvalue "<<eigenvals(index)<<endl;
  if(fabs(ratio-eigenvals(index))>1e-5){
    cout <<"WARNING: Did not pass the check for the_eignevalue solver, are you using LAPACK?!"<<endl;
  }
}


void Localize_method::run(Program_options & options, ostream & output)
{
  cout <<" start Localize_method::run"<<endl;
  unsigned int electron=0; //# of the electron that will move through grid
  if(sysprop->nelectrons(0)+sysprop->nelectrons(1)<2)
    error("The Localization methods requires at least two electrons in the system");
  Array1 <doublevar> xyz(3);
  Array2 <doublevar> resolution_array(3,3); 
  Array1 <int> D_array1(3); //dummy array1
  D_array1=0; //sets all 3 components to 0. use as counter for gridpoints

  //const double ftol=1.e-7;//precision in minimalization routine
  double fret; //values of minima
  //int iter=0; //number of iterations
  
  Array1 <doublevar> p0(norbs); //starting vector for minimization
  Array1 <doublevar> xi(norbs); //gradient of minimization
  Array2 <doublevar> RM(nmo,nmo); //rotation_matrix
  Array1 < Array2 <doublevar> > overlap_ma_R_all(ncenters); //overlap matrix within radius R
  Array2 <doublevar> overlap_total(norbs,norbs); //full overlap matrix
  doublevar determinant_RM; //determinat of rotation_matrix
  Array1 <doublevar> sphere_norm(ncenters); //
  Array1 <doublevar> dist(5);
  

  //zeroing arrays
  for(int i=0;i<ncenters;i++){
    overlap_ma_R_all(i).Resize(norbs,norbs);
    overlap_ma_R_all(i)=0.0;
  }
  overlap_total=0.0;
  sphere_norm=0.0;
  RM=0.0;
  for(int i=0;i<nmo;i++)
    RM(i,i)=1.0;


  int ndim=3;
  doublevar Volume= volume3D(latVec);
  //cout <<" Volume of the simulation cell "<<Volume<<endl;
  Array1 < Array1 <doublevar> > latice(ndim);

  for(int i=0;i<ndim;i++){
    latice(i).Resize(ndim);
    for(int j=0;j<ndim;j++)
      latice(i)(j)=latVec(i,j);
    D_array1(i)=roundoff(distance_between_a_b(latice(i),origin)/resolution);
    cout <<D_array1(i)<<"  ";
    for(int j=0;j<ndim;j++){
      resolution_array(i,j)=fabs(latice(i)(j)-origin(j))/(D_array1(i)); 
      cout << resolution_array(i,j)<<"  ";
    }
    cout <<endl;
  }

  doublevar deltaVolume=volume3D(resolution_array);
  //cout <<" deltaVolume "<<deltaVolume<<endl;;
 
  if(deltaVolume<1e-4)
    error("deltaVolume might be too small!");
    
  if(norbs<=0)
    error("number of orbitals requested is not a positive number");
  
  //calculate value of each molecular orbital at each grid point 
  // grid values with x=fastest running variable, and z=slowest

  output<<"calculating "<<D_array1(0)*D_array1(1)*D_array1(2)<<" grid points"<<endl;
  output<<"for "<< orbs.GetDim(0) <<" old molecular orbitals"<<endl;  
  doublevar percentage;
  int lenght=50;

  doublevar check_total_volume=0;
  for(int zz=0;zz<D_array1(2);zz++){
    for(int yy=0; yy<D_array1(1);yy++){
      for(int xx=0; xx<D_array1(0);xx++){
        for(int ii=0;ii<ndim;ii++){
          xyz(ii)=origin(ii)+xx*resolution_array(0)(ii)+yy*resolution_array(1)(ii)+zz*resolution_array(2)(ii);
        }
        //cout <<xyz(0)<<"  "<<xyz(1)<<"  "<<xyz(2)<<endl;
        mywalker->setElectronPos(electron,xyz); //move elec#1 to point specified by xyz
        mymomat->updateVal(mywalker,electron,0,mymovals); //recalculate MO value for elec#1
        for(int k=0;k<norbs;k++){
          for(int l=0;l<=k;l++){
            overlap_total(k,l)+=mymovals(k,0)*mymovals(l,0)*deltaVolume;
          }
        }
        check_total_volume+=deltaVolume;
        for(int i=0; i<ncenters; i++){
          mywalker->setElectronPos(electron+1,center(i));
          mywalker->updateEEDist();
          mywalker->getEEDist(electron,electron+1,dist);
          if(dist(0)<radius){
            sphere_norm(i)+=deltaVolume;
            for(int k=0;k<norbs;k++){
              for(int l=0;l<=k;l++){
                overlap_ma_R_all(i)(k,l)+=mymovals(k,0)*mymovals(l,0)*deltaVolume;
              }//l
            }//k
          }//if inside
        }//center
      }//xx
    }//yy
    percentage=doublevar(zz+1)/doublevar(D_array1(2));
    banner(percentage,lenght,cout);
  }//zz
  //cout << " Integrated total volume - calculated= "<<check_total_volume-Volume<<endl;
  if(fabs(check_total_volume-Volume)>1e-6)
    error("Check the code for dividing the simulation cell");
  
  
  //make overlap matrix symmetric
  for(int k=0; k<norbs; k++){
    for(int l=0;l<=k;l++){
      overlap_total(l,k)=overlap_total(k,l);
    }
  }
  
  doublevar det=Determinant(overlap_total,norbs);
  if(fabs(det)<1e-6)
    error("Full Overlap Matrix has linear dependency vectors");
  //cout <<" Determinant of overlap matrix "<<det<<endl;
  //doublevar multiply=exp(-log(fabs(det))/norbs);
  //cout <<"Renormalization factor for matrices: "<<multiply<<endl;
  
  cout <<"<-------------------------Full Overlap Matrix -------------------------------------->"<<endl;
  for(int k=0; k<norbs; k++){
    for(int l=0;l<=k;l++){
      cout <<overlap_total(k,l)<<"   ";
    }
    cout <<endl;
  }
  cout <<"<------------------------------------------------------------------------------------>"<<endl;
   
  
  doublevar spherevolume=1.33333333333333333333*pi*radius*radius*radius;
  for(int i=0; i<ncenters; i++){
    cout << "<-------------------------"<<i+1<<"-th center's bare overlap matrix ------------------------>" <<endl;
    for(int k=0; k<norbs; k++){
      for(int l=0;l<=k;l++){
        overlap_ma_R_all(i)(l,k)=overlap_ma_R_all(i)(k,l);
        cout << overlap_ma_R_all(i)(k,l)<<"  ";
      }
      cout <<endl; 
    }
    cout <<"<------------------------------------------------------------------------------------>"<<endl;
    output <<"Relative error on sphere integral at "<<i+1<<"-th center "<<(sphere_norm(i)-spherevolume)/spherevolume<<endl;
  }
  
  //maximize the overlap in the sphere on each center
  output<<"Starting the localization section with "<<norbs<<" orbitals"<<endl;
  output <<"number of grid points for each primitive vector direction"<<endl;
  output <<D_array1(2)<<" "<<D_array1(1)<<" "<<D_array1(0)<<endl;
  output <<endl;
  
  doublevar fret_min;
  Array1 <doublevar> pbest(norbs);
  for(int i=0; i<ncenters; i++) {
    cout <<"Center "<<i+1<<"  ";
    fret=1.0;
    fret_min=1.0;
    pbest=0.0;
    //find the bet coeficients using conjugate-grtadient minimization 
    /*
    int trymax=10;
    Minimization_function minfunc; //initialization of object
    minfunc.init_overlap_ma_R(overlap_ma_R_all_renorm(i));
    for (int ttry=0;ttry<trymax;ttry++){
      for (int ii=0;ii<norbs;ii++){
        p0(ii)=1.0-2.0*unif();
      }
      normalize(p0);
      iter=0;
      NR_conjugate_gradient <Minimization_function> minimizer(& minfunc);
      minimizer.frprmn(p0, ftol, iter, fret);
      if(fret<fret_min){
        fret_min=fret;
        pbest=p0;
      }
    }
    */
    //new version, use the largest eigenvalue of the renormalized overlap matrix to find maximum overlap
    Array1 <doublevar> gevals(norbs);
    Array2 <doublevar> gevecs(norbs,norbs);
    GeneralizedEigenSystemSolverRealSymmetricMatrices(overlap_ma_R_all(i),overlap_total,gevals,gevecs);
    //printing eigenvalues
    //for(int index=0; index<norbs; index++){
    // cout << index+1<<":  gevals= "<<gevals(index)<<"  ";
    // for (int j=0; j<norbs ; j++)
    //   cout << gevecs(j,index)<<"  ";
    // cout <<endl;
    //}

    check_the_eignevalue_solver(overlap_ma_R_all(i),overlap_total, gevals, gevecs,0);

    fret_min=1-gevals(0);
    for (int j=0; j<norbs ; j++){
      pbest(j)=gevecs(j,0);
    }
    normalize(pbest);
     
    //if(fabs(minfunc.func(pbest) -fret_min)>1e-5){
    //cout <<" minfunc.func(pbest): "<<minfunc.func(pbest)<<"  fret_min: "<<fret_min<<endl;
    //  error("Eigenvalue solver in Localize method has a bug!");
    //}
    //normalization, nt needed whe using eigen-solver
    //printout of values 
    
    cout <<"overlap= "<<gevals(0)<<endl;
    for (int j=0; j<norbs ; j++){
      RM(orbs(i),orbs(j))=pbest(j);
      cout << pbest(j)<<"  ";
    }
    cout <<endl;
    output  <<"For center "<<i+1<<", norm inside the sphere of R="<<radius<<" is " << 1-fret_min <<endl;
  }
  
  string matrix_file=options.runid;
  matrix_file+=".rotation_matrix";
  ofstream matrix_out(matrix_file.c_str());
 
  matrix_out.precision(15);
  for(int i=0;i<nmo;i++){
    int counter=1;
    for(int j=0;j<nmo;j++){
      matrix_out <<RM(i,j)<<"    ";
      if(counter % 5 ==0) matrix_out << endl;
      counter++;
    }
    matrix_out <<endl;
  }
  matrix_out.close();
    
  determinant_RM=Determinant(RM,nmo); 
  output << "determinant of RM= "<<determinant_RM<<endl;
  if (abs(determinant_RM)<1.0e-3) 
    output<<"WARNING:  Determinant of the rotation matrix = "<<determinant_RM<<" , possible linear dependency"<<endl;
  
  //print new orbitals to file
  
  string orbfile = options.runid;
  orbfile += "_loc.orb";
  // for (int k=0;k<norbs;k++)
  // temp_orbs(k)=orbs(k)-1;

  ofstream testorb(orbfile.c_str());
  mymomat->writeorb(testorb,  RM, all_orbs);
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
  for(int i=1;i<norbs;i++)
    os<<", "<<orbs(i)+1;
  os<<endl;
  os<<"resolution="<<resolution<<endl;
  return 1;
}
