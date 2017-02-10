#include "Orbital_rotation.h"

void Orbital_rotation::read(vector<string> &words, int in_ndet, 
Array3<Array1<int> > & occupation_orig, 
Array3<Array1<int> > & occupation, Array1<Array1<int> > & totoccupation){
  
  //Multiple determinant read=in
  int f=0;
  int nfunc=1;
  ndet=in_ndet;
  notactive=0;

  //Occupied list sizes
  Nocc.Resize(ndet,2);
  Nact.Resize(ndet,2);
  parms.Resize(ndet,2);

  //Assign Nocc first
  for(int det=0;det<ndet;det++){
    Nocc(det,0)=occupation_orig(f,det,0).GetDim(0);
    Nocc(det,1)=occupation_orig(f,det,1).GetDim(0);
  }
  
  //Read in active spaces, active parameters, and initial parameters 
  unsigned int pos=0;
  Array1<vector<string> > actudetstring;
  actudetstring.Resize(ndet);
  for(int det=0;det<ndet;det++){
    if(!readsection(words,pos,actudetstring(det),"VIRTUAL_SPACE_U")){
      error("Didn't list enough virtual_space_u, expected ",ndet, " got ", det);
    }
    Nact(det,0)=actudetstring(det).size()+Nocc(det,0);
  }

  pos=0;
  Array1<vector<string> > actddetstring;
  actddetstring.Resize(ndet);
  for(int det=0;det<ndet;det++){
    if(!readsection(words,pos,actddetstring(det),"VIRTUAL_SPACE_D")){
      error("Didn't list enough vritual_space_d, expected ",ndet, " got ", det);
    }
    Nact(det,1)=actddetstring(det).size()+Nocc(det,1);
  }

  pos=0;
  vector <string> activestring;
  if(!readsection(words,pos,activestring,"ACTIVE_PARMS")){
    error("Require section ACTIVE_PARMS");
  }
  if(activestring.size()!=nparms()){
    error("Require ", nparms()," terms in ACTIVE_PARMS, but only ", activestring.size()," provided.");
  }
  //Assign active parameters
  isactive.Resize(nparms());
  for(int i=0;i<isactive.GetDim(0);i++){
    isactive(i)=atoi(activestring[i].c_str());
    if(!isactive(i)){
      notactive++;
    }
  }
  
  pos=0;
  vector <string> initparmstring;
  if(!readsection(words,pos,initparmstring,"INIT_PARMS")){
    error("Require section INIT_PARMS");
  }

  //Assign active arrays
  activeoccupation_orig.Resize(nfunc,ndet,2);
  activeoccupation.Resize(nfunc,ndet,2);
  activetotoccupation.Resize(2);
  for(int det=0;det<ndet;det++){
    for(int s=0;s<2;s++){
      activeoccupation_orig(f,det,s).Resize(Nact(det,s));
      activeoccupation(f,det,s).Resize(Nact(det,s));
    }
  }
 
  
  for(int det=0;det<ndet;det++){
    for(int n=0;n<Nocc(det,0);n++){
      activeoccupation_orig(f,det,0)(n)=occupation_orig(f,det,0)(n); 
    }
    for(int n=0;n<Nocc(det,1);n++){
      activeoccupation_orig(f,det,1)(n)=occupation_orig(f,det,1)(n); 
    }
    for(int n=Nocc(det,0);n<Nact(det,0);n++){
      activeoccupation_orig(f,det,0)(n)=atoi(actudetstring(det)[n-Nocc(det,0)].c_str())-1; 
    }
    for(int n=Nocc(det,1);n<Nact(det,1);n++){
      activeoccupation_orig(f,det,1)(n)=atoi(actddetstring(det)[n-Nocc(det,1)].c_str())-1; 
    }
  }
 
  for (int s=0; s<2; s++) {
    vector <int> totocctemp;
    for(int det=0; det<ndet; det++) {
      for(int mo=0; mo < Nact(det,s); mo++) {
        int place=-1;
        int ntot=totocctemp.size();
        for(int i=0; i< ntot; i++) {
          if(activeoccupation_orig(f,det,s)(mo)==totocctemp[i]) {
            place=i;
            break;
          }
        }
        if(place==-1) { //if we didn't find the MO
          activeoccupation(f,det,s)(mo)=totocctemp.size();
          totocctemp.push_back(activeoccupation_orig(f,det,s)(mo));
        }
        else {
          //cout << "found it" << endl;
          activeoccupation(f,det,s)(mo)=place;
        }
      }
    }
    activetotoccupation(s).Resize(totocctemp.size());
    for(int i=0; i<activetotoccupation(s).GetDim(0); i++)
     activetotoccupation(s)(i) = totocctemp[i];
  }
  //Alter matrices
  occupation_orig.Resize(nfunc,ndet,2);
  occupation.Resize(nfunc,ndet,2);
  totoccupation.Resize(2);
  
  for(int s=0;s<2;s++){
    totoccupation(s).Resize(activetotoccupation(s).GetDim(0));
    for(int det=0;det<ndet;det++){
      occupation_orig(f,det,s).Resize(Nact(det,s));
      occupation(f,det,s).Resize(Nact(det,s));
    }
  }
  occupation_orig=activeoccupation_orig;
  occupation=activeoccupation;
  totoccupation=activetotoccupation;
  
  //Assign initial parameters
  for(int det=0;det<ndet;det++){
    parms(det,0).Resize(Nocc(det,0)*(Nact(det,0)-Nocc(det,0)));
    parms(det,1).Resize(Nocc(det,1)*(Nact(det,1)-Nocc(det,1)));
  }
  if(initparmstring.size()==0){
    cout<<"No initial parameters specified, setting all to zero."<<endl;
    for(int det=0;det<ndet;det++){
      for(int s=0;s<2;s++){
        for(int i=0;i<Nocc(det,s)*(Nact(det,s)-Nocc(det,s));i++){
          parms(det,s)(i)=0;
        }
      }
    }
  }else if(initparmstring.size()!=nparms()+notactive){
    error("You have an incorrect number of initial parameters, require ",nparms());
  }else{
    int offset=0;
    int k=0;
    for(int det=0;det<ndet;det++){
      for(int s=0;s<2;s++){
        for(int i=0;i<Nocc(det,s)*(Nact(det,s)-Nocc(det,s));i++){
          if(isactive(i+offset)){
            parms(det,s)(i)=atof(initparmstring[i+offset].c_str());
            //parms(det,s)(i)=atof(initparmstring[i+offset-k].c_str());
           //}else{
           // parms(det,s)(i)=0;
           // k++;
          }
        }
        offset+=Nocc(det,s)*(Nact(det,s)-Nocc(det,s));
      }
    }
  }
  //Set initial matrices, theta = 0, Rvar = I, parms = 0, R = exp(theta(init_parms))
  theta.Resize(ndet,2);
  r.Resize(ndet,2);
  rvar.Resize(ndet,2);
  rvarinv.Resize(ndet,2);
  rmult.Resize(ndet,2);
  parmsindex.Resize(ndet,2);
  for(int det=0;det<ndet;det++){
    for(int s=0;s<2;s++){
      parmsindex(det,s).Resize(Nocc(det,s)*(Nact(det,s)-Nocc(det,s))); 
      theta(det,s).Resize(Nact(det,s),Nact(det,s)); 
      r(det,s).Resize(Nact(det,s),Nact(det,s)); 
      rvar(det,s).Resize(Nact(det,s),Nact(det,s)); 
      rvarinv(det,s).Resize(Nact(det,s),Nact(det,s)); 
      rmult(det,s).Resize(Nact(det,s),Nact(det,s)); 
    }
  }
  setTheta();
  setRvar();

  for(int det=0;det<ndet;det++){
    for(int s=0;s<2;s++){
      for(int i=0;i<Nact(det,s);i++){
        for(int j=0;j<Nact(det,s);j++){
          r(det,s)(i,j)=rvar(det,s)(i,j);
        }
      }
    }
  }
 
  for(int det=0;det<ndet;det++){
    for(int s=0;s<2;s++){
      for(int i=0;i<Nocc(det,s)*(Nact(det,s)-Nocc(det,s));i++){
        parms(det,s)(i)=0;
      }
    }
  }
  
  setTheta();
  setRvar();

}

void Orbital_rotation::writeinput(string & indent, ostream & os){
 
  os<<indent<<"OPTIMIZE_MO"<<endl;
  os<<indent<<"OPTIMIZE_DATA {"<<endl;
  string indent2=indent+"  ";
  int f=0;
  for(int det=0;det<ndet;det++){
    os<<indent2<<"VIRTUAL_SPACE_U { ";
    for(int i=Nocc(det,0);i<Nact(det,0);i++){
      os<<activeoccupation_orig(f,det,0)(i)+1<<" ";
    }
    os<<"}"<<endl;
    os<<indent2<<"VIRTUAL_SPACE_D { ";
    for(int i=Nocc(det,1);i<Nact(det,1);i++){
      os<<activeoccupation_orig(f,det,1)(i)+1<<" ";
    }
    os<<"} "<<endl;
  }
  os<<indent2<<"ACTIVE_PARMS { ";
  for(int i=0;i<nparms()+notactive;i++){
    os<<isactive(i)<<" ";
  }
  os<<"}"<<endl;

  os<<indent2<<"INIT_PARMS { "<<endl;
  os<<indent+indent2;

  //Get final rotation parameters!
  Array2<Array2<dcomplex> > Ain;
  Array2<Array1<dcomplex> > W;
  Array2<Array2<dcomplex> > VL;
  Array2<Array2<dcomplex> > VR;
  Array2<Array2<dcomplex> > VRinvT;
  Array2<Array2<dcomplex> > VRinv;
  Array2<Array2<doublevar> > tmptheta;

  Ain.Resize(ndet,2);
  W.Resize(ndet,2);
  VL.Resize(ndet,2);
  VR.Resize(ndet,2);
  VRinv.Resize(ndet,2);
  VRinvT.Resize(ndet,2);
  tmptheta.Resize(ndet,2);

  int offset=0;
  for(int det=0;det<ndet;det++){
    for(int s=0;s<2;s++){
      Ain(det,s).Resize(Nact(det,s),Nact(det,s));
      tmptheta(det,s).Resize(Nact(det,s),Nact(det,s));
      for(int i=0;i<Nact(det,s);i++){
        for(int j=0;j<Nact(det,s);j++){
          Ain(det,s)(i,j)=r(det,s)(i,j);
        }
      }
      if(Nact(det,s)>0){
        VRinv(det,s).Resize(Nact(det,s),Nact(det,s));
        VRinvT(det,s).Resize(Nact(det,s),Nact(det,s));
        GeneralizedEigenSystemSolverComplexGeneralMatrices(Ain(det,s), W(det,s), VL(det,s), VR(det,s));
        TransposeInverseMatrix(VR(det,s),VRinvT(det,s),Nact(det,s));

        for(int i=0;i<Nact(det,s);i++){
          for(int j=0;j<Nact(det,s);j++){
            VRinv(det,s)(i,j)=VRinvT(det,s)(j,i);
          }
        }
        tmptheta(det,s)=0;
        for(int i=0;i<Nact(det,s);i++){
          for(int j=0;j<Nact(det,s);j++){
            for(int l=0;l<Nact(det,s);l++){
              tmptheta(det,s)(i,j)+=(VRinv(det,s)(i,l)*log(W(det,s)(l))*VR(det,s)(l,j)).real();
            }
          }
        }
      }else{
        tmptheta(det,s)=0;
      }
      int k,l;
      for(int i=0;i<Nocc(det,s)*(Nact(det,s)-Nocc(det,s));i++){
        getind(i+offset,det,s,k,l);
        if(isactive(i+offset))
          os<<tmptheta(det,s)(k,l)<<" ";
        if(i==Nocc(det,s)*(Nact(det,s)-Nocc(det,s))-1){
          offset+=Nocc(det,s)*(Nact(det,s)-Nocc(det,s));
          os<<" "<<endl;
          os<<indent+indent2;
        }
      }
    }
  }

  os<<"}"<<endl;
  os<<indent<<"}"<<endl;

}

void Orbital_rotation::setTheta(void){
  //Set elements of theta and parmsindex based on values of parms
  for(int det=0;det<ndet;det++){
    for(int s=0;s<2;s++){
      int ind=0;
      theta(det,s)=0;
      for(int i=0;i<Nact(det,s);i++){
        for(int j=0;j<=i;j++){
          if(i>(Nocc(det,s)-1) && j<=(Nocc(det,s)-1)){
            if(parms(det,s)(ind)!=0){
              theta(det,s)(j,i)=parms(det,s)(ind);  
              theta(det,s)(i,j)=-parms(det,s)(ind);
            }else{
              theta(det,s)(j,i)=parms(det,s)(ind);
              theta(det,s)(i,j)=parms(det,s)(ind);
            }
            parmsindex(det,s)(ind).Resize(2);
            parmsindex(det,s)(ind)(0)=j;
            parmsindex(det,s)(ind)(1)=i;
            ind++;
          }
        }
      }
    }
  }

}

void Orbital_rotation::setRvar(void){
  Array2<Array2<dcomplex> > Ain;
  Array2<Array1<dcomplex> > W;
  Array2<Array2<dcomplex> > VL;
  Array2<Array2<dcomplex> > VR;
  Array2<Array2<dcomplex> > VRinvT;
  Array2<Array2<dcomplex> > VRinv;

  Ain.Resize(ndet,2);
  W.Resize(ndet,2);
  VL.Resize(ndet,2);
  VR.Resize(ndet,2);
  VRinv.Resize(ndet,2);
  VRinvT.Resize(ndet,2);

  //Set Rvar=exp(theta), and Rvarinv=Rvar.Transpose
  for(int det=0;det<ndet;det++){
    for(int s=0;s<2;s++){
      Ain(det,s).Resize(Nact(det,s),Nact(det,s));
      for(int i=0;i<Nact(det,s);i++){
        for(int j=0;j<Nact(det,s);j++){
          Ain(det,s)(i,j)=theta(det,s)(i,j);
        }
      }
      //SVD
      if(Nact(det,s)>0){
        VRinv(det,s).Resize(Nact(det,s),Nact(det,s));
        VRinvT(det,s).Resize(Nact(det,s),Nact(det,s));
        GeneralizedEigenSystemSolverComplexGeneralMatrices(Ain(det,s), W(det,s), VL(det,s), VR(det,s));
        TransposeInverseMatrix(VR(det,s),VRinvT(det,s),Nact(det,s));

        for(int i=0;i<Nact(det,s);i++){
          for(int j=0;j<Nact(det,s);j++){
            VRinv(det,s)(i,j)=VRinvT(det,s)(j,i);
          }
        }

        rvar(det,s)=0.0;
        //Exponentiate D
        for(int i=0;i<Nact(det,s);i++){
          for(int j=0;j<Nact(det,s);j++){
            for(int l=0;l<Nact(det,s);l++){
              rvar(det,s)(i,j)+=(VRinv(det,s)(i,l)*exp(W(det,s)(l))*VR(det,s)(l,j)).real();
            }
          }
        }
      }else{
        rvar(det,s)=1;
      }
    }

    //Set inverse
    for(int s=0;s<2;s++){
      for(int i=0;i<Nact(det,s);i++){
        for(int j=0;j<Nact(det,s);j++)
          rvarinv(det,s)(i,j)=rvar(det,s)(j,i);
      }
    }
  }

}

void Orbital_rotation::lockIn(Array1<doublevar> & parmsout){
  for(int s=0;s<2;s++){
    for(int det=0;det<ndet;det++){
      parms(det,s)=0;
    }
  }
  setTheta();
  setRvar();
  getParms(parmsout);

}

void Orbital_rotation::getParms(Array1<doublevar> & parmsout){
  parmsout.Resize(nparms());
  parmsout=0;
}

void Orbital_rotation::setParms(Array1<doublevar> & parmsin){
  //Revert previous rotation by Rvar
  rmult=r;
  for(int det=0;det<ndet;det++){
    for(int s=0;s<2;s++){
      MultiplyMatrices(rvarinv(det,s),rmult(det,s),r(det,s),Nact(det,s));  
    }
  }

  //Apply new rotation by Rvar
  int offset=0;
  int k=0;
  //VERY IMPORTANT: the det loop must be oustide, else offset has the incorrect meaning
  for(int det=0;det<ndet;det++){
    for(int s=0;s<2;s++){
      for(int i=0;i<Nocc(det,s)*(Nact(det,s)-Nocc(det,s));i++){
        if(isactive(offset+i)){
          parms(det,s)(i)=parmsin(i+offset-k);
        }else{
          parms(det,s)(i)=0;
          k++;
        }
      }
      offset+=Nocc(det,s)*(Nact(det,s)-Nocc(det,s));
    }
  }

  setTheta();
  setRvar();
  rmult=r;
  for(int det=0;det<ndet;det++){
    for(int s=0;s<2;s++){
      MultiplyMatrices(rvar(det,s),rmult(det,s),r(det,s),Nact(det,s));  
    }
  }
}

int Orbital_rotation::nparms(void){
  int val=0;
  for(int det=0;det<ndet;det++){
    for(int s=0;s<2;s++)
      val+=Nocc(det,s)*(Nact(det,s)-Nocc(det,s));
  }
  return val-notactive;
}

void Orbital_rotation::getind(int n, int& mydet, int& mys, int& i, int& j){
  //For a given parameter n, get the determinant it lives in, which spin, and the indices 
  //inside theta(mydet,mys) that corresponds to parameter n
  mydet=0;
  mys=0;
  int val=0;
  int found=0;
  for(int det=0;det<ndet;det++){
    if(!found){
      for(int s=0;s<2;s++){
        val+=Nocc(det,s)*(Nact(det,s)-Nocc(det,s));
        if(n<val){
          mydet=det;
          mys=s;
          found=1;
          break;
        }
      }
    }else{
      break;
    }
  }
  if(mys==0){
    for(int det=0;det<mydet;det++){
      for(int s=0;s<2;s++){
        n-=Nocc(det,s)*(Nact(det,s)-Nocc(det,s));
      }
    }
  }else{
    for(int det=0;det<mydet;det++){
      for(int s=0;s<2;s++){
        n-=Nocc(det,s)*(Nact(det,s)-Nocc(det,s));
      }
    }
    n-=Nocc(mydet,0)*(Nact(mydet,0)-Nocc(mydet,0));
  }
  
  i=parmsindex(mydet,mys)(n)(0);
  j=parmsindex(mydet,mys)(n)(1);

}

