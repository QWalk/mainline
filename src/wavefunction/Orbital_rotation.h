#include <stdio.h>
#include "Qmc_std.h"
#include "Array.h"
#include "qmc_io.h"
#include "MatrixAlgebra.h"
#include "MO_matrix.h"
#include "Wavefunction.h"

/* 
  0. The new nparms() and getnmo() (good)
  1. The initial matrices (good)
  2. The new totoccupation (good)
  3. Make sure that this implementation does the exact same thing as before, but generalized! (good)
    a. It alters totoccupation predictably
    b. It assigns the initial matrices and parameters correctly
  4. Replace 0s with fs, then continue implementing (good)
  5. For Test_method for checking spatial derivatives this involves the following
        i. calcLap (called in updateLap) (done)
        ii. updateInverse(e) called in updateLap(e)) (done)
        iii. updateValNoInverse(e) (called in updateVal(e) if using Clark updates) (done)
        iv. getDetLap(e) (called in getLap) (done)
        v. evalTestPos (called in external programs) (done)
  6. Implemented recently  
     i. 5(iii,v) from above (done) 
     ii. The parameter derivatives, the spatial derivatives work well (done)
     iii. getParms(done),setParms(done),lockIn(done)
     iv. Implement getind (done), getParmDeriv (done) and writeinput
     v. Test Linear method for multiple determinants (done)
       a. Testing: Run Li + C with single determinants (S + S-J) with plots and with 
       multiple determinants (S + S-J) with plots as well (done)
     vi. Check runtime (4*ndet*N^4 for gradient/gradderiv,4*ndet*N^3 for lapDet/Inv), document code (done)
     vii. Code up writeinput (done) and the current_parms printout 
     viii. Template the functions not the class (just let rotation matrices be doublevar)
     ix. Add in the features for not varying some parameters from the other implementation 
      a. Involves changing ACTIVE -> VIRTUAl (done)
      b. Then figure out how to freeze certain parameters (done)
     x. Fetch updates from qwalk-devel and merge, then test for a last time
        a. Check that VMC works with the orbital parameters in place (done)
        b. Test read-in and initialization (what if I don't include the right number of det blocks, etc) (done)
  7. Next steps
     i. Write up slides, clean-up code and REPO, and pull request 
      a. In slides: write up all of the features that are built in (done) 
      b. Clean-up via documentation and removing any superfluous comments (done)
      c. Figure out how to make the REPO smaller (drafts), and then tell Lucas/pull request!
  
  WHILE TESTING: Make sure that you check the 1 and 2 determinant cases!

*/

class Orbital_rotation{
public:
  void read(vector<string> &words, int in_ndet, 
  Array3<Array1<int> > & occupation_orig, 
  Array3<Array1<int> > & occupation, Array1<Array1<int> > & totoccupation);
  int nparms(void);                     //Number of variational parameters
  void lockIn(Array1<doublevar> & parms);       //Lock-in variational parameters
  void getParms(Array1<doublevar> & parmsout);  //get current variational parameters
  void setParms(Array1<doublevar> & parmsin);   //set variational parameters
  
  template <class T> 
  void rotMoVals(int det, int s, Array1<T> & orbvals);      //Rotate the moVals based on current R
  //Get parameter derivatives
  template <class T>
  void getParmDeriv(Array1<log_value<T> > & detwt,
  Array3<T> & moVal, Array3<Array2<T> >& inverse, Array3<log_value<T> > & detVal,  
  Parm_deriv_return & deriv); 
  int getnmo(void){    
    return max(activetotoccupation(0).GetDim(0),activetotoccupation(1).GetDim(0));
  }
  void writeinput(string & indent, ostream & os); //Write input file after each optimization iteration
  Array2<int> Nocc; //Number of orbitals occupied (det,spin)
  Array2<int> Nact; //Number of orbitals in active space (det,spin)
private:
  Array2<Array1<doublevar> > parms;	//Variational parameters (det,spin)(number)
  Array2<Array2<doublevar> > theta;	//Anti-hermitian matrix (det,spin)(i,j)
  Array2<Array2<doublevar> > r;		//Compounded rotation matrix (det,spin)(i,j)
  Array2<Array2<doublevar> > rvar;      //Variable rotation matrix (det,spin)(i,j)
  Array2<Array2<doublevar> > rvarinv;   //Inverse (det,spin)(i,j)
  Array2<Array2<doublevar> > rmult;     //Used to compound the matrices (det,spin)(i,j)

  //For a parameter of spin, gives the indices in theta (det,spin)(parm)(0 or 1 -> i or j)
  Array2<Array1<Array1<int> > > parmsindex;  
  Array3<Array1<int> > activeoccupation;    //occupation but also including Active Space orbitals!
  Array3<Array1<int> >activeoccupation_orig; //occupation_origin but also including Active Space orbitals!
  Array1<Array1<int> >activetotoccupation; //Totoccupation but also including Active Space orbitals!
  Array1<int> isactive; //Which parameters are active and which arent

  int ndet;        //Number of determinants
  int notactive=0; //Total number of inactive parameters 
  void setTheta(void); //Set theta matrix based on parms 
  void setRvar(void); //Rvar=exp(theta)
  void getind(int n,int& det,int& s,int& i,int& j); //Get index of parameter n
};

template <class T>
void Orbital_rotation::rotMoVals(int det, int s, Array1<T> & orbvals) {
  Array1<T> tmp;
  tmp.Resize(Nact(det,s));

  for(int i=0;i<Nact(det,s);i++){
    tmp(i)=0.0;
    for(int c=0;c<Nact(det,s);c++){
      tmp(i)+=r(det,s)(i,c)*orbvals(c);
    }
  }
  orbvals=tmp;

}
template <class T>
void Orbital_rotation::getParmDeriv(Array1<log_value<T> > & detwt,
Array3<T> & moVal, Array3<Array2<T> >& inverse, 
Array3<log_value<T> > & detVal, Parm_deriv_return & deriv){ 
  
  //moVal() dimension, electron, orbital 
  //inverse() function #, det #, spin 
  //detVal() function #, determinant #, spin 
 
  int f=0;
  Array1<T> u; //Vector used for Sherman-Morrison udpate in calculation of gradient
  Array1<T> w; //Vector used for Sherman-Morrison update in calculation of gradderiv
  Array1<T> v; //Vector used for Sherman-Morrison update in calculation of lapInv 
  //Inverse of spatial derivatives of wavefunction matrix, (det,electron,dimension)(i,j)
  Array3<Array2<T> > lapInv; 
  Array3<log_value<T> >lapDet;  //Contains 1/det(lapInv), (det,electron,dimension)
  lapInv.Resize(ndet,Nocc(0,0)+Nocc(0,1),4);
  lapDet.Resize(ndet,Nocc(0,0)+Nocc(0,1),4); 

  //Calculate lapInv and lapDet  
  //lapDet(det,e,d)=det(grad_e^d(D_det^s(e))
  //lapInv(det,e,d)=grad_e^d(D_det^s(e))^-1
  int ts=0;
  for(int det=0;det<ndet;det++){
    for(int e=0;e<Nocc(det,0)+Nocc(det,1);e++){
      if(e<Nocc(det,0)){ts=0;}
      else{ts=1;}

      //Calculate elements of row to update in inverse
      v.Resize(Nocc(det,ts));
      for(int d=0;d<4;d++){
        v=0;
        for(int k=0;k<Nocc(det,ts);k++){
          for(int c=0;c<Nact(det,ts);c++){
            v(k)+=r(det,ts)(k,c)*moVal(d+1,e,activeoccupation(f,det,ts)(c));
          }
        }
          
        lapInv(det,e,d).Resize(Nocc(det,ts),Nocc(det,ts)); 
        for(int x=0;x<Nocc(det,ts);x++){
          for(int y=0;y<Nocc(det,ts);y++){
            lapInv(det,e,d)(x,y)=inverse(f,det,ts)(y,x);
          }
        }
        //Get new inverse and ratio of determinants
        if(ts==0){
          lapDet(det,e,d)=1.0/InverseUpdateRow(lapInv(det,e,d),v,e,Nocc(det,ts));
        }else{
          lapDet(det,e,d)=1.0/InverseUpdateRow(lapInv(det,e,d),v,e-Nocc(det,0),Nocc(det,ts));
        }
        lapDet(det,e,d)*=detVal(f,det,ts);
      }
    }
  }

  //Calculate gradient and gradderiv  
  //gradient(n)=d/d_pn ln(det(psi))
  //gradderiv(n,e,t)=d/d_pn (grad_e^(t) ln(det(psi))), 
  //where grad_e^(1,2,3)=(d/dx_e, d/dy_e, d/dz_e) and grad_e^4=lap_e
  int q=0;
  for(int n=0;n<nparms()+notactive;n++){
    if(isactive(n)){
      //These are the four parameters related to n
      int ni=0;
      int nj=0;
      int nd=0;
      int ns=0;
      int no=0;
     
      getind(n,nd,ns,ni,nj);
      if(ns==0){no=1;}
      
      //Calculate elements of column we want to replace 
      u.Resize(Nocc(nd,ns));
      u=0;
      for(int l=0;l<Nocc(nd,ns);l++){
        for(int c=0;c<Nact(nd,ns);c++){
          u(l)+=r(nd,ns)(nj,c)*moVal(0,l+ns*Nocc(nd,0),activeoccupation(f,nd,ns)(c));
        }
      }

      //Calculate trace of (psi^-1 d/d_pn psi) 
      T tmp=0;
      for(int l=0;l<Nocc(nd,ns);l++){
        tmp+=inverse(f,nd,ns)(l,ni)*u(l);
      }
      
      //Calculate gradient
      deriv.gradient(n-q)=tmp;
      deriv.gradient(n-q)*=detwt(nd).val()*detVal(f,nd,ns).val()*detVal(f,nd,no).val();
      T wval=0;
      for(int d=0;d<ndet;d++){
        wval+=detwt(d).val()*detVal(f,d,0).val()*detVal(f,d,1).val();
      }
      deriv.gradient(n-q)/=wval;
       
      int se=0;
      int oe=0;
      for(int t=0;t<4;t++){
        for(int e=0;e<Nocc(0,0)+Nocc(0,1);e++){
          if(e<Nocc(0,0)){se=0;oe=1;}
          else{se=1;oe=0;}
          
          //grad_e(psi)
          T tmp1=0;
          for(int i=0;i<ndet;i++)
            tmp1+=detwt(i).val()*lapDet(i,e,t).val()*detVal(f,i,oe).val();

          //psi
          T tmp2=0;
          for(int i=0;i<ndet;i++)
            tmp2+=detwt(i).val()*detVal(f,i,se).val()*detVal(f,i,oe).val();
          
          if(se==no){
            //gradderiv when s(e)!=s(n),i.e. p_n and e in different spin spaces
            deriv.gradderiv(n-q,e,t)=deriv.gradient(n-q)*(lapDet(nd,e,t).val()/detVal(f,nd,no).val());
            deriv.gradderiv(n-q,e,t)-=deriv.gradient(n-q)*tmp1/tmp2;
          }else{
            //gradderiv when s(e)=s(n),i.e. p_n and e in same spin spaces
            w.Resize(u.GetDim(0));
            w=u;
            w(e-ns*Nocc(nd,0))=0;

            //Calculate elements of column we want to update
            for(int c=0;c<Nact(nd,ns);c++){
              w(e-ns*Nocc(nd,0))+=r(nd,ns)(nj,c)*moVal(t+1,e,activeoccupation(f,nd,ns)(c));
            }

            //Calculate trace of grad_e(psi)^-1 d/d_pn grad_e(psi) 
            T tmp3=0;
            for(int l=0;l<Nocc(nd,ns);l++){
              tmp3+=lapInv(nd,e,t)(ni,l)*w(l);  
            }

            //Calculate gradderiv
            deriv.gradderiv(n-q,e,t)=detwt(nd).val()*detVal(f,nd,no).val()*lapDet(nd,e,t).val()*tmp3;
            deriv.gradderiv(n-q,e,t)-=deriv.gradient(n-q)*tmp1;
            deriv.gradderiv(n-q,e,t)/=tmp2;
          }
        }//e loop
      }//t loop
    }else{
      q++; 
    }
  }//n loop
}
