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

  //Check if we want orthogonal 
  pos=0;
  orthog=0;
  if(haskeyword(words,pos,"ORTHOGONAL"))
    orthog=1;

  

  //Manipulate data to get actudetstring and actddetstring
  pos=0;
  vector <string> initparmstring;
  Array1<vector<string> > actudetstring;
  Array1<vector<string> > actddetstring;
  actudetstring.Resize(ndet);
  actddetstring.Resize(ndet);
  Array1<vector<string> >groupparms;
  groupstrings.Resize(ndet);
  groupparms.Resize(ndet);
  int ngroup=0;
  for(int det=0;det<ndet;det++){
    vector<string> detstring;
    int givenparms=0;
    unsigned int parmpos=0;
    ngroup=0;
    if(!orthog){
      if(!readsection(words,pos,detstring,"DET")){
        error("Did not list enough determinants in OPTIMIZE_DATA"); 
      }
    }else{
      //This isn't the most memory efficient becuase now 
      //every determinant stores the same information. We can 
      //change this later if it becomes necessary.
      detstring=words;
    }
    unsigned int detpos=0;
    if(!readsection(detstring,detpos,groupparms(det),"PARAMETERS"))
      groupparms(det).resize(0);

    detpos=0;
    vector<string> orbgroupstring;
    while(readsection(detstring,detpos,orbgroupstring,"ORB_GROUP")){
      ngroup++;
    }
    groupstrings(det).Resize(ngroup);
    detpos=0;
    ngroup=0;
    
    while(readsection(detstring,detpos,orbgroupstring,"ORB_GROUP")){
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

  //Assign nparms_
  if(!orthog) nparms_=nparms();
  else{
    nparms_=0;
    for(int det=0;det<ndet;det++){
      for(int s=0;s<2;s++){
        nparms_+=(int)Nact(det,s)*(Nact(det,s)-1)/2;
      }
    }
  }

  //Useful maps, tells us whether an orbital is globally active
  if(orthog){
    //Actives
    int z=0;
    for(int det=0;det<ndet;det++){
      for(int s=0;s<2;s++){
        for(int l=0;l<occupation_orig(f,det,s).Size();l++){
          if(globalactive.count(occupation_orig(f,det,s)(l))==0) {
            globalactive[occupation_orig(f,det,s)(l)]=true;
            globalindex[occupation_orig(f,det,s)(l)]=z;
            z++;
          }
        }
      }
    }
    
    //Virtuals
    for(int det=0;det<ndet;det++){
      for(int s=0;s<2;s++){
        if(s==0){
          for(int l=0;l<actudetstring(det).size();l++){
            int k=atoi(actudetstring(det)[l].c_str())-1;
            if(globalactive.count(k)==0) { 
              globalactive[k]=false;
              globalindex[k]=z;
              z++;
            }
          }
        }else{
          for(int l=0;l<actddetstring(det).size();l++){
            int k=atoi(actddetstring(det)[l].c_str())-1;
            if(globalactive.count(k)==0) {
              globalactive[k]=false;
              globalindex[k]=z;
              z++;
            }
          }
        }
      }
    }
  }
  
  //DEBUG
  /*for (map<int,bool>::iterator it=globalactive.begin(); it!=globalactive.end(); ++it)
      std::cout << it->first << " => " << it->second << '\n';
  for (map<int,int>::iterator it=globalindex.begin(); it!=globalindex.end(); ++it)
        std::cout << it->first << " => " << it->second << '\n';
  */

  //Manipulate data to get activestring
  int off=0;
  vector <string> activestring;
  activestring.resize(nparms_);
  //Non orthogonal
  if(!orthog){
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
                pass=1;
                break;
              }
            }
            if(pass) {
              activestring[i*Nocc(det,s)+j+off]="1"; 
            }else { 
              activestring[i*Nocc(det,s)+j+off]="0";
            }
          }
        }
        off+=Nocc(det,s)*(Nact(det,s)-Nocc(det,s));
      }
    }
  }else{ //Orthogonal. Could be merged with above loop, but I want to keep them separated for debugging purposes.
    int p=0;
    for(int det=0;det<ndet;det++){
      for(int s=0;s<2;s++){
        for(int i=0;i<Nact(det,s);i++){
          for(int j=i+1;j<Nact(det,s);j++){
            int pass=0;
            int passone=0;
            int passtwo=0;
            int one=0;
            int two=0;
            if(i>=Nocc(det,s)) one=(!s)?atoi(actudetstring(det)[i-Nocc(det,s)].c_str()):atoi(actddetstring(det)[i-Nocc(det,s)].c_str());
            else one=occupation_orig(f,det,s)(i)+1;
            
            if(j>=Nocc(det,s)) two=(!s)?atoi(actudetstring(det)[j-Nocc(det,s)].c_str()):atoi(actddetstring(det)[j-Nocc(det,s)].c_str());
            else two=occupation_orig(f,det,s)(j)+1;
          
            //Check that they share the same group
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
                pass=1;
                break;
              }
            }
           
            //Set active parameters
            if(pass) activestring[p]="1";
            else activestring[p]="0"; 
            
            p++;
          }
        }
      }
    }
  }//[GOOD]
  //Manipulate data to get initparmstring 
  if(!orthog) {
    for(int det=0;det<ndet;det++){
      vector<string> parmu;
      vector<string> parmd;
      int numup=(int)Nact(det,0)*(Nact(det,0)-1)/2;
      int nump=numup+(int)Nact(det,1)*(Nact(det,1)-1)/2;
      if(groupparms(det).size()!=0){
        if(groupparms(det).size()!=nump)
          error("Expected ",nump," parameters, got ",groupparms(det).size());
        for(int i=0;i<numup;i++)
          parmu.push_back(groupparms(det)[i]);
        for(int i=0;i<nump-numup;i++)
          parmd.push_back(groupparms(det)[i+numup]);
      }else{
        for(int i=0;i<numup;i++)
          parmu.push_back("0");
        for(int i=0;i<nump-numup;i++)  
          parmd.push_back("0");
      }
      for(int k=0;k<parmu.size();k++){
        initparmstring.push_back(parmu[k]);
      }
      for(int k=0;k<parmd.size();k++){
        initparmstring.push_back(parmd[k]);
      }
    }
  }else{
    int nump=(int) globalactive.size()*(globalactive.size()-1)/2;
    if(groupparms(0).size()!=0){
      if(groupparms(0).size()!=nump) error("Expected ",nump," parameters, got ",groupparms(0).size());
      for(int i=0;i<nump;i++) initparmstring.push_back(groupparms(0)[i]); //Goes directly into thetaGlobal;
    }else{
      for(int i=0;i<nump;i++) initparmstring.push_back("0");
    }
  }
  
  //Assign active parameters
  isactive.Resize(nparms_);
  for(int i=0;i<isactive.GetDim(0);i++){
    isactive(i)=atoi(activestring[i].c_str());
    if(!isactive(i)){
      notactive++;
    }
  }
  
  //alter nparms_
  if(!orthog) nparms_=nparms();
  else nparms_-=notactive;

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
  int totparms=0;
  if(!orthog){
    for(int det=0;det<ndet;det++){
      parms(det,0).Resize((int)Nact(det,0)*(Nact(det,0)-1)/2);
      parms(det,1).Resize((int)Nact(det,1)*(Nact(det,1)-1)/2);
      totparms+=(int)Nact(det,0)*(Nact(det,0)-1)/2 + (int)Nact(det,1)*(Nact(det,1)-1)/2;
    }
  }else{
    for(int det=0;det<ndet;det++){
      parms(det,0).Resize((int)Nact(det,0)*(Nact(det,0)-1)/2);
      parms(det,1).Resize((int)Nact(det,1)*(Nact(det,1)-1)/2);
    }
    totparms=(int)globalactive.size()*(globalactive.size()-1)/2; 
  }
  if(randomparms==0){
    if(initparmstring.size()==0){
      single_write(cout,"No initial parameters specified, setting all to zero.\n");
      for(int det=0;det<ndet;det++){
        for(int s=0;s<2;s++){
          for(int i=0;i<(int)Nact(det,s)*(Nact(det,s)-1)/2;i++)
            parms(det,s)(i)=0;
        }
      }
    }else if(initparmstring.size()!=totparms){
      error("You have an incorrect number of initial parameters, require ",totparms);
    }else{
      if(!orthog){
        int offset=0;
        int k=0;
        for(int det=0;det<ndet;det++){
          for(int s=0;s<2;s++){
            for(int i=0;i<(int)Nact(det,s)*(Nact(det,s)-1)/2;i++){
              parms(det,s)(i)=atof(initparmstring[i+offset].c_str());
            }
            offset+=(int)Nact(det,s)*(Nact(det,s)-1)/2;
          }
        }
      }else{
        for(int det=0;det<ndet;det++){
          for(int s=0;s<2;s++){
            int p=0;
            for(int i=0;i<Nact(det,s);i++){
              for(int j=i+1;j<Nact(det,s);j++){
                int iloc=activeoccupation_orig(f,det,s)(i);
                int jloc=activeoccupation_orig(f,det,s)(j);
                if(globalactive[iloc] || globalactive[jloc]){
                  int iglob=globalindex[iloc];
                  int jglob=globalindex[jloc];
                  int offset=0;
                  if(iglob<jglob){
                    for(int k=0;k<iglob;k++) offset+=globalactive.size()-1-k;
                    parms(det,s)(p)=atof(initparmstring[offset+jglob-iglob-1].c_str());
                  }else{
                    for(int k=0;k<jglob;k++) offset+=globalactive.size()-1-k;
                    parms(det,s)(p)=-atof(initparmstring[offset+iglob-jglob-1].c_str());
                  }
                }else{
                  parms(det,s)(p)=0;                
                } 
                p++;
              } 
            }
          }             
        }
      }
    }
  }else{
    //cout<<"Randomizing initial parameters"<<endl;
    single_write(cout,"Randomizing initial parameters.\n");
    int offset=0;
    int k=0;
    srand (time(NULL));
    //Warmup
    for(int det=0;det<ndet;det++)
      for(int s=0;s<2;s++)
        for(int i=0;i<Nocc(det,s)*(Nact(det,s)-Nocc(det,s));i++)
          doublevar tmp=(rand()*2.0*randomparms/float(RAND_MAX))-randomparms;
    
    //Randomize initparms
    for(int det=0;det<ndet;det++){
      for(int s=0;s<2;s++){
        for(int i=0;i<Nact(det,s)*(Nact(det,s)-1)/2;i++){
          parms(det,s)(i)=(rand()*2.0*randomparms/float(RAND_MAX))-randomparms;
        }
        offset+=(int)Nact(det,s)*(Nact(det,s)-1)/2;
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
  
  //Set Theta manually
  for(int det=0;det<ndet;det++){
    for(int s=0;s<2;s++){
      int ind=0;
      theta(det,s)=0;
      for(int i=0;i<Nact(det,s);i++){
        for(int j=i;j<Nact(det,s);j++){
          if(i==j) { theta(det,s)(i,j)=0; }
          else{
            theta(det,s)(i,j)=parms(det,s)(ind);
            theta(det,s)(j,i)=-parms(det,s)(ind);
            ind++;
          }
        }
      }
    }
  }
  
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

  if(!orthog){
    //Resize parameters back to shape
    for(int det=0;det<ndet;det++){
      parms(det,0).Resize(Nocc(det,0)*(Nact(det,0)-Nocc(det,0)));
      parms(det,1).Resize(Nocc(det,1)*(Nact(det,1)-Nocc(det,1)));
    }
  }

  //Zero out parameters and reset theta and Rvar
  for(int det=0;det<ndet;det++){
    for(int s=0;s<2;s++){
      for(int i=0;i<parms(det,s).GetDim(0);i++){
        parms(det,s)(i)=0;
      }
    }
  }

  setTheta();
  setRvar();
  //[GOOD]
  //Make restriction matrix
  if(orthog){
    restMat.Resize(nparms_,nparms_);
    int q=0;
    for(int n=0;n<nparms_+notactive;n++){  
      if(isactive(n)){
        int ni,nj,ns,nd; //local indices
        int n_i,n_j; //global indices
        getind(n,nd,ns,ni,nj);
        n_i=activeoccupation_orig(0,nd,ns)(ni);
        n_j=activeoccupation_orig(0,nd,ns)(nj);
        int p=0;
        for(int m=0;m<nparms_+notactive;m++){
          if(isactive(m)){
            int mi,mj,ms,md; //local indices
            int m_i,m_j; //global indices
            getind(m,md,ms,mi,mj);
            m_i=activeoccupation_orig(0,md,ms)(mi);
            m_j=activeoccupation_orig(0,md,ms)(mj);

            if((n_i==m_i) && (n_j==m_j)) { restMat(n-q,m-p)=1; restMat(m-p,n-q)=1; }
            else if((n_i==m_j) && (n_j==m_i)) { restMat(n-q,m-p)=-1; restMat(m-p,n-q)=-1; }
            else { restMat(n-q,m-p)=0; restMat(m-p,n-q)=0; }
            
          }else { p++; }
        }
      }else { q++; }
    }

    //Get number of independent parameters using a BFS
    //http://math.hws.edu/eck/cs327_s04/chapter9.pdf
    bool visited[nparms_];
    for(int i=0;i<nparms_;i++){
      visited[i]=false;
    }
    
    nindep=0;
    for(int n=0;n<nparms_;n++){
      if(!visited[n]){
        nindep++;
        queue<int> q;
        q.push(n);
        visited[n]=true;
        vector<int>tmpcomp;
        tmpcomp.push_back(n);
        while(!q.empty()){
          int w=q.front();
          q.pop();
          for(int m=0;m<nparms_;m++){
            if(restMat(w,m)!=0){
              if(!visited[m]){
                visited[m]=true;
                q.push(m);
                tmpcomp.push_back(m);
              }
            }
          }
        }
        concomp.push_back(tmpcomp);
      }
    }
  }//[GOOD]
}

void Orbital_rotation::writeinput(string & indent, ostream & os){
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
    }
  }
  //Write input file
  os<<indent<<"OPTIMIZE_MO"<<endl;
  os<<indent<<"OPTIMIZE_DATA {"<<endl;
  int f=0;
  if(!orthog){
    for(int det=0;det<ndet;det++){
      os<<indent2<<"DET {"<<endl;
      for(int g=0;g<groupstrings(det).GetDim(0);g++){
        os<<indent2<<indent<<"ORB_GROUP { ";
        for(int l=0;l<groupstrings(det)(g).size();l++)
          os<<groupstrings(det)(g)[l].c_str()<<" ";
        os<<"}"<<endl;
      }
      //Print out final parameters
      os<<indent2<<indent<<"PARAMETERS { ";
      for(int s=0;s<2;s++){
        for(int i=0;i<Nact(det,s);i++){
          os<<indent2<<indent2;
          for(int j=i+1;j<Nact(det,s);j++){
            os<<fixed<<showpoint<<setprecision(7);
            os<<tmptheta(det,s)(i,j)<<" ";
          }     
          os<<"\n";
        }
      }
      os<<indent2<<indent<<"}"<<endl;
      os<<indent2<<"} "<<endl;
    }
  }else{
    os<<indent2<<"ORTHOGONAL"<<endl;
    for(int g=0;g<groupstrings(0).GetDim(0);g++){
      os<<indent2<<"ORB_GROUP { ";
      for(int l=0;l<groupstrings(0)(g).size();l++)
        os<<groupstrings(0)(g)[l].c_str()<<" ";
      os<<"}"<<endl;
    }
   
    //Figure out how to print out parameters
    Array2<doublevar> globaltheta;
    globaltheta.Resize(globalactive.size(),globalactive.size());
    for(int i=0;i<globalactive.size();i++){
      for(int j=0;j<globalactive.size();j++){
        globaltheta(i,j)=0;
      }
    }

    //Build globaltheta
    for(int det=0;det<ndet;det++){
      for(int s=0;s<2;s++){
        for(int i=0;i<Nact(det,s);i++){
          for(int j=i+1;j<Nact(det,s);j++){
            //Get global indices
            int iglob=activeoccupation_orig(f,det,s)(i);
            int jglob=activeoccupation_orig(f,det,s)(j);
            //Make sure it hasn't been filled in
            if(globaltheta(globalindex[iglob],globalindex[jglob])==0){
              //Make sure that it's a-a or a-v
              if(globalactive[iglob] || globalactive[jglob]){
                globaltheta(globalindex[iglob],globalindex[jglob])=tmptheta(det,s)(i,j);
                globaltheta(globalindex[jglob],globalindex[iglob])=tmptheta(det,s)(j,i);
              }
            }
          }
        }
      }
    }
    //Print it out
    os<<indent2<<"PARAMETERS {"<<endl;
    for(int i=0;i<globalactive.size();i++){
      os<<indent2<<indent;
      for(int j=i+1;j<globalactive.size();j++){
        os<<fixed<<showpoint<<setprecision(7);
        os<<globaltheta(i,j)<<" ";
      }
      os<<"\n";
    }
    os<<indent2<<"}"<<endl;
  }
  os<<indent<<"}"<<endl;
}

void Orbital_rotation::setTheta(void){
  //Set elements of theta and parmsindex based on values of parms
  for(int det=0;det<ndet;det++){
    for(int s=0;s<2;s++){
      int ind=0;
      theta(det,s)=0;
      if(!orthog){
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
      }else{
        parmsindex(det,s).Resize((int)Nact(det,s)*(Nact(det,s)-1)/2);
        for(int i=0;i<Nact(det,s);i++){
          for(int j=i+1;j<Nact(det,s);j++){
            theta(det,s)(i,j)=parms(det,s)(ind);
            theta(det,s)(j,i)=-parms(det,s)(ind);
            parmsindex(det,s)(ind).Resize(2);
            parmsindex(det,s)(ind)(0)=i;
            parmsindex(det,s)(ind)(1)=j;
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

  if(!orthog){
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
  }else{
    //Expand the input parmsin to be appropriate size
    Array1<doublevar> intparmsin; //Internal parameter representation
    intparmsin.Resize(nparms_);
    for(int i=0;i<parmsin.Size();i++){
      for(int j=0;j<concomp[i].size();j++){
        intparmsin(concomp[i][j])=restMat(concomp[i][j],concomp[i][0])*parmsin(i);
      }
    }

    //Now assign parms(det,s)
    int offset=0;
    int k=0;
    for(int det=0;det<ndet;det++){
      for(int s=0;s<2;s++){
        for(int i=0;i<(int)Nact(det,s)*(Nact(det,s)-1)/2;i++){
          if(isactive(offset+i)){
            parms(det,s)(i)=intparmsin(i+offset-k);
          }else{
            parms(det,s)(i)=0;
            k++;
          }
        }
        offset+=(int)Nact(det,s)*(Nact(det,s)-1)/2;
      }
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
  if(!orthog){
    int val=0;
    for(int det=0;det<ndet;det++){
      for(int s=0;s<2;s++)
        val+=Nocc(det,s)*(Nact(det,s)-Nocc(det,s));
    }
    return val-notactive;
  }else{
    return nindep;
  }
}

void Orbital_rotation::getind(int n, int& mydet, int& mys, int& i, int& j){
  //For a given parameter n, get the determinant it lives in, which spin, and the indices 
  //inside theta(mydet,mys) that corresponds to parameter n

  //Since I set isactive in a particular order, getind has to satisfy that ordering 
  //as well.
  mydet=0;
  mys=0;
  int val=0;
  int found=0;
  for(int det=0;det<ndet;det++){
    if(!found){
      for(int s=0;s<2;s++){
        if(!orthog) val+=Nocc(det,s)*(Nact(det,s)-Nocc(det,s));
        else val+=(int)Nact(det,s)*(Nact(det,s)-1)/2;
        
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
        if(!orthog) n-=Nocc(det,s)*(Nact(det,s)-Nocc(det,s));
        else n-=(int)Nact(det,s)*(Nact(det,s)-1)/2;
      }
    }
  }else{
    for(int det=0;det<mydet;det++){
      for(int s=0;s<2;s++){
        if(!orthog) n-=Nocc(det,s)*(Nact(det,s)-Nocc(det,s));
        else n-=(int)Nact(det,s)*(Nact(det,s)-1)/2; 
      }
    }
    if(!orthog) n-=Nocc(mydet,0)*(Nact(mydet,0)-Nocc(mydet,0));
    else n-=(int)Nact(mydet,0)*(Nact(mydet,0)-1)/2;
  }

  i=parmsindex(mydet,mys)(n)(0);
  j=parmsindex(mydet,mys)(n)(1);

}

