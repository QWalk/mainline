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
#include "Wannier_method.h"
#include "qmc_io.h"
#include "System.h"
#include "MatrixAlgebra.h"
#include "ulec.h"
/*!
Read the "words" from the method section in the input file
via doinput() parsing, and store section information in private
variables orbs, resolution, and minmax.
Set up MO_matrix and Sample_point objects for wavefunction from input
*/
void Wannier_method::read(vector <string> words,
                       unsigned int & pos,
                       Program_options & options)
{

  
  sys=NULL;
  allocate(options.systemtext[0],  sys);

  if(! readvalue(words,pos=0,resolution,"RESOLUTION"))
    resolution=.2;
  if(!readvalue(words,pos=0, out_orbs, "OUT_ORB"))
    out_orbs=options.runid+".orb";
  if(!readvalue(words,pos=0, shake,"SHAKE"))
    shake=1./10.;

  vector <vector < string> > orbgroup_txt;
  pos=0;
  vector < string> dummy;
  while( readsection(words,pos,dummy,"ORB_GROUP")) { 
    orbgroup_txt.push_back(dummy);
  }
  orbital_groups.Resize(orbgroup_txt.size());
  for(int i=0; i< orbital_groups.GetDim(0); i++) { 
    int norb_this=orbgroup_txt[i].size();
    orbital_groups[i].Resize(norb_this);
    for(int j=0; j< norb_this; j++) { 
      orbital_groups[i][j]=atoi(orbgroup_txt[i][j].c_str())-1;
    }
  }



  vector <string> orbtext;
  if(!readsection(words, pos=0, orbtext, "ORBITALS")){
      error("Need ORBITALS");
  }
  allocate(orbtext, sys, mymomat);

  mywalker=NULL;
  sys->generateSample(mywalker);

  origin.Resize(3);
  origin=0;
  LatticeVec.Resize(3,3);
   
  if(!sys->getBounds(LatticeVec,origin))
    error("need getBounds");

}

//----------------------------------------------------------------------
int Wannier_method::showinfo(ostream & os) { 
  os << "Wannier localization method " << endl;
  return 1;
}

//----------------------------------------------------------------------

//A is upper triangular and gives the angle to rotate for a pair i,j
void make_rotation_matrix(const Array2 <doublevar> & A, 
    Array2 <doublevar> & B) {
  int n=A.GetDim(0);
  assert(A.GetDim(1)==n);
  B.Resize(n,n);
  B=0.;
  for(int i=0; i< n; i++) B(i,i)=1.0;

  Array2 <doublevar> tmp(n,n),tmp2(n,n);
  Array1 <doublevar> tmpveci(n),tmpvecj(n);
  doublevar s,co;
  for(int i=0; i< n; i++) { 
    for(int j=i+1; j < n; j++) { 
      if(fabs(A(i,j)) > 1e-12) {
        co=cos(A(i,j));
        s=sin(A(i,j));
        /*
        tmp2=0.0;
        tmp=0.0;
        for(int k=0; k< n; k++) tmp2(k,k)=1.0;
        tmp2(i,i)=co;
        tmp2(i,j)=-s;
        tmp2(j,i)=s;
        tmp2(j,j)=co;
        
        //very very inefficient, for testing..
        //MultiplyMatrices(tmp2,B,tmp,n);
        tmp=0;
        for(int a=0; a< n; a++) { 
          for(int c=0; c< n; c++) { 
            tmp(a,c)+=tmp2(a,a)*B(a,c);
          }
        }
        for(int c=0; c< n; c++) { 
          tmp(i,c)+=tmp2(i,j)*B(j,c);
          tmp(j,c)+=tmp2(j,i)*B(i,c);
        }
        B=tmp;
        */
        //------
        for(int c=0; c< n; c++) {
          tmpveci(c)=B(i,c);
          tmpvecj(c)=B(j,c);
        }
        for(int c=0;c < n; c++) {
          B(i,c)*=co;
          B(j,c)*=co;
        }
        for(int c=0; c< n; c++) {
          B(i,c)+= -s*tmpvecj(c);
          B(j,c)+= s*tmpveci(c);
        }
      }
    }
  }

  /*
  cout << "rotation matrix " << endl;
  for(int i=0; i< n; i++) { 
    for(int j=0; j< n; j++) { 
      cout << B(i,j) << " ";
    }
    cout << endl;
  }
  */


}
//A is an antisymmetric matrix and B is the output rotation matrix
void make_rotation_matrix_notworking(const Array2 <doublevar> & A,
      Array2 <doublevar> & B) { 
  int n=A.GetDim(0);
  assert(A.GetDim(1)==n);
  B.Resize(n,n);

  Array2 <dcomplex> skew(n,n),VL(n,n),VR(n,n);
  Array1 <dcomplex> evals(n);
  for(int i=0; i< n; i++) { 
    for(int j=0; j < n; j++) { 
      skew(i,j)=A(i,j);
    }
  }
  GeneralizedEigenSystemSolverComplexGeneralMatrices(skew,evals,VL,VR);

  cout << "evals " << endl;
  for(int i=0; i< n; i++) cout << evals(i) << " ";
  cout << endl;
  cout << "VR " << endl;
  for(int i=0; i< n; i++) {
    for(int j=0; j< n; j++) { 
      cout << VR(i,j) << " ";
    }
    cout << endl;
  }
  cout << "VL " << endl;
  for(int i=0; i< n; i++) { 
    for(int j=0; j< n; j++) { 
      cout << VL(i,j) << " ";
     }
    cout << endl;
  }
  //this is horribly inefficient,most likely

  skew=dcomplex(0.0,0.); //we don't need that any more so we reuse it
  Array2 <dcomplex> work(n,n);
  work=dcomplex(0.0,0.);
  for(int i=0; i< n; i++) { 
    skew(i,i)=exp(evals(i));
  }
  for(int i=0; i< n; i++) { 
    for(int j=0; j<n; j++) { 
      for(int k=0; k< n; k++) { 
        work(i,k)+=skew(i,j)*VR(j,k);
      }
    }
  }
  skew=dcomplex(0.,0.);
  for(int i=0; i< n; i++) { 
    for(int j=0; j<n; j++) { 
      for(int k=0; k< n; k++) { 
//skew(i,k)+=conj(VL(i,j))*work(j,k);
        skew(i,k)=conj(VR(j,i))*work(j,k);
      }
    }
  }


  cout << "rotation " << endl;
  for(int i=0; i< n; i++) { 
    for(int j=0; j< n; j++) { 
      cout << skew(i,j) << " ";
     }
    cout << endl;
  }
  
  
}
//----------------------------------------------------------------------
/*!
 
*/
void Wannier_method::run(Program_options & options, ostream & output) {
  int ngroups=orbital_groups.GetDim(0);
  
  int ntotorbs=0;
  for(int i=0; i< ngroups; i++) { 
    ntotorbs+=orbital_groups(i).GetDim(0);
  }
  Array1 <int> allorbs(ntotorbs);


  int count=0;
  for(int i=0; i< ngroups; i++) { 
    for(int j=0; j< orbital_groups(i).GetDim(0); j++) { 
      allorbs(count++)=orbital_groups(i)(j);
    }
  }

  /*
  Array2 <doublevar> rtmp(allorbs.GetDim(0),allorbs.GetDim(0));
  rtmp=0.0;
  for(int i=0; i< allorbs.GetDim(0); i++) { 
    rtmp(i,i)=1.0;
  }

  ofstream test("test.orb");
  mymomat->writeorb(test, rtmp,allorbs);
  test.close();

  return;
*/

  Array3 <dcomplex> eikr;
  Array2 <doublevar> phi2phi2;
  cout << "calc overlap " << endl;
  Array1 <Array2 <doublevar> > R(ngroups);

  for(int i=0; i< ngroups; i++) { 
    calculate_overlap(orbital_groups(i),eikr,phi2phi2);
    optimize_rotation(eikr,R(i));
  }

  //int nmo_tot=mymomat->getNmo();
  //Array2 <doublevar> Rtot(nmo_tot,nmo_tot);
  //Rtot=0.0;
  //for(int i=0; i< nmo_tot; i++) Rtot(i,i)=1.0;
  
  Array2 <doublevar> Rtot(ntotorbs,ntotorbs);
  Rtot=0.0;

  count=0;
  for(int i=0; i< ngroups; i++) { 
    int norb=orbital_groups(i).GetDim(0);
    for(int j=0; j < norb; j++) { 
      for(int k=0; k< norb; k++) { 
        Rtot(count+j,count+k)=R(i)(j,k);
      }
    }

    count+=norb;
  }
  

  
  ofstream testorb(out_orbs.c_str());
  mymomat->writeorb(testorb, Rtot,allorbs);
  testorb.close();
  


}

//----------------------------------------------------------------------
void Wannier_method::calculate_overlap(Array1 <int> & orb_list, 
    Array3 <dcomplex> & eikr, Array2 <doublevar> & phi2phi2) { 


  Array1 <Array1 <int> > tmp_list(1);
  tmp_list(0)=orb_list;
  mymomat->buildLists(tmp_list);
  int list=0;
  
  Array2 <doublevar>  gvec(3,3);
  if(!sys->getPrimRecipLattice(gvec) ) error("Need getPrimRecipLattice");

  int norb=orb_list.GetDim(0);
  eikr.Resize(3,norb,norb);
  phi2phi2.Resize(norb,norb);
  eikr=0.0;
  phi2phi2=0.0;

  Array1 <doublevar> prop(3);
  Array1 <doublevar> n_lat(3);
  n_lat=0.0;
  for(int d=0;d < 3; d++) { 
    for(int d1=0; d1 < 3; d1++) 
      n_lat(d)+=LatticeVec(d,d1)*LatticeVec(d,d1);
    n_lat(d)=sqrt(n_lat(d));
    prop(d)=resolution/double(n_lat(d));
    //cout << "prop " << prop(d) <<  " res " << resolution << " norm " << n_lat(d) << endl;
  }

  cout << "calculating...." << endl;
  mymovals.Resize(norb,1);
  Array1 <doublevar> r(3),x(3);
  Array1 <doublevar> norm_orb(norb);
  norm_orb=0.0;
  int totpts=0;
  for(r(0)=origin(0); r(0) < 1.0; r(0)+=prop(0)) {
    cout << "percent done: " << r(0)*100 << endl;
    for(r(1)=origin(1); r(1) < 1.0; r(1)+=prop(1)) {
      for(r(2)=origin(2); r(2) < 1.0; r(2)+=prop(2)) {
        
        totpts++;
        for(int d=0; d< 3; d++) { 
          x(d)=0.0;
          for(int d1=0; d1 < 3; d1++) { 
            x(d)+=r(d1)*LatticeVec(d1,d);
          }
        }
        
        mywalker->setElectronPos(0,x);
        mymomat->updateVal(mywalker,0,list,mymovals);
        for(int i=0; i < norb; i++) norm_orb(i)+=mymovals(i,0)*mymovals(i,0);
        for(int d=0; d< 3; d++) { 
          doublevar gdotr=0;
          for(int d1=0; d1 < 3; d1++) {
            gdotr+=gvec(d,d1)*x(d1);
          }
          dcomplex exp_tmp=exp(dcomplex(0,-1)*2.*pi*gdotr);
          for(int i=0; i< norb; i++) { 
            for(int j=0; j< norb; j++) { 
              eikr(d,i,j)+=exp_tmp*mymovals(i,0)*mymovals(j,0);
            }
          }
        }
        for(int i=0; i< norb; i++) { 
          for(int j=0; j< norb; j++) { 
            phi2phi2(i,j)+=mymovals(i,0)*mymovals(i,0)*mymovals(j,0)*mymovals(j,0);
          }
        }
      }
    }
  }

  for(int i=0; i< norb; i++) norm_orb(i)=sqrt(norm_orb(i));

  for(int d=0; d< 3; d++) { 
    for(int i=0; i < norb; i++) { 
      for(int j=0; j< norb; j++) { 
        eikr(d,i,j)/=norm_orb(i)*norm_orb(j);
        //cout << d << " " << i << " " <<  j << " " << eikr(d,i,j) << endl;
      }
    }
  }
  cout << "square overlap " << endl;
  for(int i=0; i < norb; i++) { 
    for(int j=0; j< norb; j++) { 
      phi2phi2(i,j)/=(norm_orb(i)*norm_orb(i)*norm_orb(j)*norm_orb(j));
      cout << sqrt(phi2phi2(i,j)) << " ";
    }
    cout << endl;
  }


}
//----------------------------------------------------------------------
//
//
//
doublevar Wannier_method::evaluate_local(const Array3 <dcomplex> & eikr,
    Array2 <doublevar> & Rgen, Array2 <doublevar> & R) { 
  int norb=Rgen.GetDim(0);
  Array2 <doublevar> gvec(3,3);
  sys->getPrimRecipLattice(gvec);
  Array1 <doublevar> gnorm(3);
  gnorm=0;
  for(int d=0;d < 3; d++) {
    for(int d1=0; d1 < 3; d1++) { 
      gnorm(d)+=gvec(d,d1)*gvec(d,d1);
    }
    gnorm(d)=sqrt(gnorm(d));
  }

  Array2 <dcomplex> tmp(norb,norb),tmp2(norb,norb);
  
  make_rotation_matrix(Rgen,R);
  doublevar func=0;
  for(int d=0; d< 3; d++) {
    tmp=0.0;
    for(int i=0; i< norb; i++) { 
      for(int j=0; j< norb; j++) {
        for(int k=0; k< norb; k++) { 
          tmp(i,k)+=eikr(d,i,j)*R(k,j);
        }
      }
    }
    tmp2=0;
    for(int i=0; i< norb; i++) { 
      for(int j=0; j< norb; j++) { 
        for(int k=0; k< norb; k++) { 
          tmp2(i,k)+=R(i,j)*tmp(j,k);
        }
      }
    }
    //cout << "========for d=" << d << endl;

    for(int i=0; i< norb; i++)  { 
      doublevar f=  -log(norm(tmp2(i,i)))/(gnorm(d)*gnorm(d));
      func+=f;

      // cout << setw(9) << f;
    }
  }
    //cout << endl;
  return func/(3*norb);

}
doublevar Wannier_method::eval_tstep(Array3 <dcomplex> & eikr, Array2 <doublevar> & Rgen,
    Array2 <doublevar> & Rgen_save, Array2 <doublevar> & deriv, doublevar tstep,
    Array2 <doublevar> & R) { 
  int norb=Rgen.GetDim(0);
  for(int ii=0; ii< norb;ii++) { 
    for(int jj=ii+1; jj < norb; jj++) { 
      Rgen(ii,jj)=Rgen_save(ii,jj)-tstep*deriv(ii,jj);
    }
  }
  return evaluate_local(eikr,Rgen,R);
}

//----------------------------------------------------------------------
void Wannier_method::optimize_rotation(Array3 <dcomplex> &  eikr,
    Array2 <doublevar> & R ) { 
    
  int norb=eikr.GetDim(1);
  Array2 <doublevar> gvec(3,3);
  sys->getPrimRecipLattice(gvec);
  Array1 <doublevar> gnorm(3);
  gnorm=0;
  for(int d=0;d < 3; d++) {
    for(int d1=0; d1 < 3; d1++) { 
      gnorm(d)+=gvec(d,d1)*gvec(d,d1);
    }
    gnorm(d)=sqrt(gnorm(d));
  }
  for(int i=0; i< norb; i++) { 
    cout << "rloc2 " << i << " ";
    for(int d=0; d< 3; d++) { 
        cout << -log(norm(eikr(d,i,i)))/(gnorm(d)*gnorm(d)) << " ";
    }
    cout << endl;
  }


  Array2 <doublevar> Rgen(norb,norb),Rgen_save(norb,norb); 
  //R(norb,norb);
  R.Resize(norb,norb);
  //Array2 <dcomplex> tmp(norb,norb),tmp2(norb,norb);
  //Shake up the angles, since often the original orbitals 
  //are at a maximum and derivatives are zero.
  Array2 <doublevar> deriv(norb,norb);
  Rgen=0.0;
  for(int ii=0; ii< norb; ii++) { 
    for(int jj=ii+1; jj< norb; jj++) { 
      Rgen(ii,jj)=rng.gasdev()*pi*shake;
    }
  }
  for(int step=0; step < 800; step++) { 
    doublevar fbase=evaluate_local(eikr,Rgen,R);
    for(int ii=0; ii <norb; ii++) { 
      cout << "deriv ";
      for(int jj=ii+1; jj < norb; jj++) { 
        doublevar save_rgeniijj=Rgen(ii,jj);
        doublevar h=1e-6;
        Rgen(ii,jj)+=h;
        doublevar func=evaluate_local(eikr,Rgen,R);
        deriv(ii,jj)=(func-fbase)/h;
        Rgen(ii,jj)=save_rgeniijj;
        cout << deriv(ii,jj) << " ";
      }
      cout << endl;
    }

    doublevar rloc_thresh=0.0001;
    

    Rgen_save=Rgen;
    doublevar best_func=1e99, best_tstep=0.0;
    doublevar bracket_tstep=0.0;
    doublevar last_func=fbase;
    for(doublevar tstep=0.01; tstep < 20.0; tstep*=2.0) { 
      doublevar func=eval_tstep(eikr,Rgen,Rgen_save,deriv,tstep,R);
      cout << "tstep " << tstep << " func " << func << endl;
      if(func > fbase or func > last_func) {
        bracket_tstep=tstep;
        break;
      }
      else last_func=func;
    }

    cout << "bracket_tstep " << bracket_tstep << endl;
    doublevar resphi=2.-(1.+sqrt(5.))/2.;
    doublevar a=0, b=resphi*bracket_tstep, c=bracket_tstep;
    doublevar af=fbase, bf=eval_tstep(eikr,Rgen,Rgen_save,deriv,b,R), cf=eval_tstep(eikr,Rgen,Rgen_save,deriv,bracket_tstep,R);
    cout << "first step  a,b,c " << a << " " << b << "  " << c 
        << " funcs " << af << " " << bf << " " << cf << endl;
    
    for(int it=0; it < 20; it++) { 
      doublevar d,df;
      if( (c-b) > (b-a))   
        d=b+resphi*(c-b);
      else 
        d=b-resphi*(b-a);
      df=eval_tstep(eikr,Rgen,Rgen_save,deriv,d,R);
      if(df < bf) { 
        if( (c-b) > (b-a) ) {
          a=b;
          af=bf;
          b=d;
          bf=df;
        }
        else { 
          c=b;
          cf=bf;
          b=d;
          bf=df;
        }
      }
      else { 
        if( (c-b) > (b-a) ) { 
          c=d;
          cf=df;
        }
        else { 
          a=d;
          af=df;
        }
      }
      cout << "step " << it << " a,b,c " << a << " " << b << "  " << c 
        << " funcs " << af << " " << bf << " " << cf << endl;
    }
    best_tstep=b;
    /*
    bool made_move=false;
    while (!made_move) { 
      for(doublevar tstep=0.00; tstep < max_tstep; tstep+=0.1*max_tstep) { 
        for(int ii=0; ii< norb;ii++) { 
          for(int jj=ii+1; jj < norb; jj++) { 
            Rgen(ii,jj)=Rgen_save(ii,jj)-tstep*deriv(ii,jj);
          }
        }
        doublevar func=evaluate_local(eikr,Rgen,R);
        if(func < best_func) { 
          best_func=func;
          best_tstep=tstep;
        }
        cout << "    tstep " << tstep << "   " << func << endl;
      }
      if(abs(best_tstep) < 0.2*max_tstep)
        max_tstep*=0.5;
      else if(abs(best_tstep-max_tstep) < 1e-14)
        max_tstep*=2.0;
      else made_move=true;
    }
    */


    for(int ii=0; ii< norb;ii++) { 
      for(int jj=ii+1; jj < norb; jj++) { 
        Rgen(ii,jj)=Rgen_save(ii,jj)-best_tstep*deriv(ii,jj);
      }
    }
    doublevar func2=evaluate_local(eikr,Rgen,R);
    doublevar max_change=0;
    for(int ii=0; ii < norb; ii++) { 
      for(int jj=ii+1; jj< norb; jj++) { 
         doublevar change=abs(Rgen(ii,jj)-Rgen_save(ii,jj));
         if(change > max_change) max_change=change;
      }

    }    
    cout << "tstep " << best_tstep << " rms " << sqrt(func2) <<  " bohr max change " << max_change <<endl;
    doublevar threshold=0.0001;
    if(max_change < threshold) break;
    if(abs(best_func-fbase) < rloc_thresh) break;
    
    

    /*
    bool moved=false;
    
    for(int ii=0; ii< norb; ii++) { 
      for(int jj=ii+1; jj< norb; jj++) { 
        doublevar save_rgeniijj=Rgen(ii,jj);
        doublevar best_del=0;
        doublevar best_f=1e99;
        for(doublevar del=-0.5; del < 0.5; del+=0.05) { 
          cout << "############ for del = " << del << endl;

          Rgen(ii,jj)=save_rgeniijj+del;
          doublevar func=evaluate_local(eikr,Rgen,R);

          if(func < best_f) { 
            best_f=func;
            best_del=del;
          }
          cout << "func " << func << endl;
        }

        Rgen(ii,jj)=save_rgeniijj+best_del;
        if(abs(best_del) > 1e-12) moved=true;
      }
    }
    if(!moved) break;
    */
  }
  make_rotation_matrix(Rgen,R);
  
}
