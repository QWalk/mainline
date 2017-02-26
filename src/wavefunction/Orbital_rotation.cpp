#include "Orbital_rotation.h"

void Orbital_rotation::read(vector<string> &words, int in_ndet, 
Array3<Array1<int> > & occupation_orig, 
Array3<Array1<int> > & occupation, Array1<Array1<int> > & totoccupation){
  //This read-in is somewhat long because I am using a new input file 
  //but the old data structure. Hence blocks of code with "Manipulate data"
  //are just there to take the new input format and map it to the old data structures.

  //New read-in
  int f=0;
  int nfunc=1;
  ndet=in_ndet;
  notactive=0;
  randomparms=0;

  //Occupied list sizes
  Nocc.Resize(ndet,2);
  Nact.Resize(ndet,2);
  parms.Resize(ndet,2);

  //Assign Nocc 
  for(int det=0;det<ndet;det++){
    Nocc(det,0)=occupation_orig(f,det,0).GetDim(0);
    Nocc(det,1)=occupation_orig(f,det,1).GetDim(0);
  }

  //Check if we want random initial parms
  unsigned int pos=0;
  if(haskeyword(words,pos,"RANDOM_PARMS")){
    randomparms=atof(words[pos+1].c_str()); 
    if(randomparms==0)
      error("Please specify a magnitude for random parameters");
  }
 
  //Manipulate data to get actudetstring and actddetstring
  pos=0;
  vector <string> initparmstring;
  Array1<vector<string> > actudetstring;
  Array1<vector<string> > actddetstring;
  actudetstring.Resize(ndet);
  actddetstring.Resize(ndet);
  Array1<Array1<vector<string> > >groupparms;
  groupstrings.Resize(ndet);
  groupparms.Resize(ndet);
  int ngroup=0;
  for(int det=0;det<ndet;det++){
    vector<string> detstring;
    ngroup=0;
    if(!readsection(words,pos,detstring,"DET")){
      error("Did not list enough determinants in OPTIMIZE_DATA"); 
    }
    unsigned int detpos=0;
    vector<string> orbgroupstring;
    while(readsection(detstring,detpos,orbgroupstring,"ORB_GROUP")){
      ngroup++;
    }
    groupstrings(det).Resize(ngroup);
    groupparms(det).Resize(ngroup);
    detpos=0;
    ngroup=0;
    while(readsection(detstring,detpos,orbgroupstring,"ORB_GROUP")){
      //Initial parameters
      unsigned int tmp=detpos;
      if(!strncmp(detstring[detpos].c_str(),"PARAMETERS",10)) {
        readsection(detstring,tmp,groupparms(det)(ngroup),"PARAMETERS");
      }else
        groupparms(det)(ngroup).resize(0);
      
      //Other stuff 
      groupstrings(det)(ngroup)=orbgroupstring;
      ngroup++;
      int is0=0;
      int is1=0;
      for(int i=0;i<orbgroupstring.size();i++){
        int val=atoi(orbgroupstring[i].c_str())-1;
        is0=0;
        is1=0;
        //Check if element is occupied
        for(int k=0;k<occupation_orig(f,det,0).GetDim(0);k++){
          if(val==occupation_orig(f,det,0)(k)){
            is0=1;
          }
        }
        for(int k=0;k<occupation_orig(f,det,1).GetDim(0);k++){
          if(val==occupation_orig(f,det,1)(k)){
            is1=1;
          }
        }
        vector<string> tmp=orbgroupstring;
        tmp.erase(tmp.begin()+i);
        for(int j=0;j<tmp.size();j++){
          if(is0){
            int one=1;
            int two=1;
            //See if element is occupied
            for(int k=0;k<occupation_orig(f,det,0).GetDim(0);k++)
              if(atoi(tmp[j].c_str())-1==occupation_orig(f,det,0)(k))
                one=0;
            //See if element exists in actudetstring
            for(int k=0;k<actudetstring(det).size();k++)
              if(atoi(tmp[j].c_str())==atoi(actudetstring(det)[k].c_str()))
                two=0;
            //Check for two conditions and add
            if(one && two)
              actudetstring(det).push_back(tmp[j]);
          }
          if(is1){
            int one=1;
            int two=1;
            //See if element is occupied
            for(int k=0;k<occupation_orig(f,det,1).GetDim(0);k++)
              if(atoi(tmp[j].c_str())-1==occupation_orig(f,det,1)(k))
                one=0;
            //See if element exists in actddetstring
            for(int k=0;k<actddetstring(det).size();k++)
              if(atoi(tmp[j].c_str())==atoi(actddetstring(det)[k].c_str()))
                two=0;
            //Check for two conditions and add
            if(one && two)
              actddetstring(det).push_back(tmp[j]);
          }
        } 
      }
    }//End read group
  }
  
  //Assign Nact
  for(int det=0;det<ndet;det++){
    Nact(det,0)=actudetstring(det).size()+Nocc(det,0);
    Nact(det,1)=actddetstring(det).size()+Nocc(det,1);
  }

  //Manipulate data to get activestring
  int off=0;
  vector <string> activestring;
  activestring.resize(nparms());
  Array1<Array1<int> >numup;
  Array1<Array1<int> >nump;
  numup.Resize(ndet);
  nump.Resize(ndet);
  for(int det=0;det<ndet;det++){
    numup(det).Resize(groupstrings(det).GetDim(0));
    numup(det)=0;
    nump(det).Resize(groupstrings(det).GetDim(0));
    nump(det)=0;
  }
  for(int det=0;det<ndet;det++){
    for(int s=0;s<2;s++){
      for(int i=0;i<Nact(det,s)-Nocc(det,s);i++){
        for(int j=0;j<Nocc(det,s);j++){
          //Figure out if i and j share an orbital group
          int pass=0;
          int passone=0;
          int passtwo=0;
          int one=0;
          if(s==0){
            one=atoi(actudetstring(det)[i].c_str());
          }else{
            one=atoi(actddetstring(det)[i].c_str());
          } 
          int two=occupation_orig(f,det,s)(j)+1;
          //Loop over groups
          for(int k=0;k<groupstrings(det).GetDim(0);k++){
            passone=0;
            passtwo=0;
            //Loop over elements in a group
            for(int l=0;l<groupstrings(det)(k).size();l++){
               if(one==atoi(groupstrings(det)(k)[l].c_str()))
                 passone=1;
               if(two==atoi(groupstrings(det)(k)[l].c_str()))
                 passtwo=1;
            }
            if(passone && passtwo) {
              if(s==0) numup(det)(k)++; 
              nump(det)(k)++;
              pass=1;
              break;
            }
          }
          if(pass) {
            activestring[i*Nocc(det,s)+j+off]=to_string(1); 
          }else { 
            activestring[i*Nocc(det,s)+j+off]=to_string(0);
          }
        }
      }
      off+=Nocc(det,s)*(Nact(det,s)-Nocc(det,s));
    }
  }

  //Manipulate data to get initparmstring and parminfodet/g
  vector<string> parminfodet;
  vector<string> parminfog;
  for(int det=0;det<ndet;det++){
    vector<string> parmu;
    vector<string> parmd;
    vector<string> parminfodetu;
    vector<string> parminfodetd;
    vector<string> parminfogu;
    vector<string> parminfogd;
    for(int g=0;g<groupstrings(det).GetDim(0);g++){
      if(groupparms(det)(g).size()!=0){
        for(int i=0;i<numup(det)(g);i++){
          parmu.push_back(groupparms(det)(g)[i]);
          parminfodetu.push_back(to_string(det));
          parminfogu.push_back(to_string(g));
        }
        for(int i=0;i<nump(det)(g)-numup(det)(g);i++){
          parmd.push_back(groupparms(det)(g)[i+numup(det)(g)]);
          parminfodetd.push_back(to_string(det));
          parminfogd.push_back(to_string(g));
        }
      }else{
        for(int i=0;i<numup(det)(g);i++){
          parmu.push_back(to_string(0));
          parminfodetu.push_back(to_string(det));
          parminfogu.push_back(to_string(g));
        }
        for(int i=0;i<nump(det)(g)-numup(det)(g);i++){
          parmd.push_back(to_string(0));
          parminfodetd.push_back(to_string(det));
          parminfogd.push_back(to_string(g));
        }
      }
    }
    for(int k=0;k<parmu.size();k++){
      initparmstring.push_back(parmu[k]);
      parminfodet.push_back(parminfodetu[k]);
      parminfog.push_back(parminfogu[k]);
    }
    for(int k=0;k<parmd.size();k++){
      initparmstring.push_back(parmd[k]);
      parminfodet.push_back(parminfodetd[k]);
      parminfog.push_back(parminfogd[k]);
    }
  }

  //Assign parminfo
  parminfo.Resize(parminfodet.size());
  for(int i=0;i<parminfo.GetDim(0);i++){
    parminfo(i).Resize(2);
    parminfo(i)(0)=atoi(parminfodet[i].c_str());
    parminfo(i)(1)=atoi(parminfog[i].c_str());
  }

  //Assign active parameters
  isactive.Resize(nparms());
  for(int i=0;i<isactive.GetDim(0);i++){
    isactive(i)=atoi(activestring[i].c_str());
    if(!isactive(i)){
      notactive++;
    }
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
        if(place==-1) { 
          activeoccupation(f,det,s)(mo)=totocctemp.size();
          totocctemp.push_back(activeoccupation_orig(f,det,s)(mo));
        }
        else {
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
  if(randomparms==0){
    if(initparmstring.size()==0){
      cout<<"No initial parameters specified, setting all to zero."<<endl;
      for(int det=0;det<ndet;det++){
        for(int s=0;s<2;s++){
          for(int i=0;i<Nocc(det,s)*(Nact(det,s)-Nocc(det,s));i++){
            parms(det,s)(i)=0;
          }
        }
      }
    }else if(initparmstring.size()!=nparms()){
      error("You have an incorrect number of initial parameters, require ",nparms());
    }else{
      int offset=0;
      int k=0;
      for(int det=0;det<ndet;det++){
        for(int s=0;s<2;s++){
          for(int i=0;i<Nocc(det,s)*(Nact(det,s)-Nocc(det,s));i++){
            if(isactive(i+offset)){
              parms(det,s)(i)=atof(initparmstring[i+offset-k].c_str());
             }else{
              parms(det,s)(i)=0;
              k++;
            } 
          }
          offset+=Nocc(det,s)*(Nact(det,s)-Nocc(det,s));
        }
      }
    }
  }else{
    cout<<"Randomizing initial parameters"<<endl;
    int offset=0;
    int k=0;
    srand (time(NULL));
    //Warmup
    for(int det=0;det<ndet;det++)
      for(int s=0;s<2;s++)
        for(int i=0;i<Nocc(det,s)*(Nact(det,s)-Nocc(det,s));i++)
          doublevar tmp=(rand()*2.0*randomparms/float(RAND_MAX))-randomparms;
    
    //Randomize initparms
    cout<<"Initial parameters:"<<endl;
    for(int det=0;det<ndet;det++){
      for(int s=0;s<2;s++){
        for(int i=0;i<Nocc(det,s)*(Nact(det,s)-Nocc(det,s));i++){
          if(isactive(i+offset)){
            //Random numbers chosen between [-randomparms, randomparms] uniformly
            parms(det,s)(i)=(rand()*2.0*randomparms/float(RAND_MAX))-randomparms;
            cout<<parms(det,s)(i)<<endl;
           }else{
            parms(det,s)(i)=0;
            k++;
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
  
  //This doesn't work correctly yet
  string indent2=indent+"  ";
  //Get final rotation parameters!
  Array2<Array2<dcomplex> > Ain;
  Array2<Array1<dcomplex> > W;
  Array2<Array2<dcomplex> > VL;
  Array2<Array2<dcomplex> > VR;
  Array2<Array2<dcomplex> > VRinvT;
  Array2<Array2<dcomplex> > VRinv;
  Array2<Array2<doublevar> > tmptheta;
  Array1<doublevar> parms;

  Ain.Resize(ndet,2);
  W.Resize(ndet,2);
  VL.Resize(ndet,2);
  VR.Resize(ndet,2);
  VRinv.Resize(ndet,2);
  VRinvT.Resize(ndet,2);
  tmptheta.Resize(ndet,2);
  parms.Resize(nparms());

  int offset=0;
  int q=0;
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
        if(isactive(i+offset)){
          parms(i+offset-q)=tmptheta(det,s)(k,l);
        }else
          q++;
        if(i==Nocc(det,s)*(Nact(det,s)-Nocc(det,s))-1)
          offset+=Nocc(det,s)*(Nact(det,s)-Nocc(det,s));
      }
    }
  }

  //Write input file
  os<<indent<<"OPTIMIZE_MO"<<endl;
  os<<indent<<"OPTIMIZE_DATA {"<<endl;
  int f=0;
  for(int det=0;det<ndet;det++){
    os<<indent2<<"DET {"<<endl;
    for(int g=0;g<groupstrings(det).GetDim(0);g++){
      os<<indent2<<indent<<"ORB_GROUP { ";
      for(int l=0;l<groupstrings(det)(g).size();l++)
        os<<groupstrings(det)(g)[l].c_str()<<" ";
      os<<"}"<<endl;
      //Print out final parameters 
      os<<indent2<<indent<<"PARAMETERS { ";
      for(int i=0;i<parms.GetDim(0);i++)
        if(parminfo(i)(0)==det && parminfo(i)(1)==g)
          os<<parms(i)<<" ";
      os<<"}"<<endl;
    }
    os<<indent2<<" } "<<endl;
  }
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

