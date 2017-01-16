//------------------------------------------------------------------------
//src/Pfaff_wf_data.cpp

#include "Qmc_std.h"
#include "qmc_io.h"
#include "Pfaff_wf_data.h"
#include "Wavefunction_data.h"
#include "Pfaff_wf.h"
#include "MatrixAlgebra.h"
#include <algorithm>



/*!
*/
void Pfaff_wf_data::read(vector <string> & words, unsigned int & pos,
                        System * sys)
{
  coef_eps=1e-12; //default value
  //cout <<"Start wf read"<<endl;
  vector <string> strpfwt;
  vector <vector <string> > statevec;
  unsigned int startpos=pos;
  //unsigned int tmp;

  vector <string> npairsstr;
  if( mpi_info.node==0 )
    cout<<"Pfaffian"<<endl;

  readsection(words, pos, npairsstr, "NPAIRS");
  if(npairsstr.size() != 3)
    error("NPAIRS must have 3 elements");
  npairs.Resize(3);
  npairs(0)=atoi(npairsstr[0].c_str());
  npairs(1)=atoi(npairsstr[1].c_str());
  npairs(2)=atoi(npairsstr[2].c_str());
  assert(npairs(0)>=npairs(1)); 
  if(npairs(0)!=sys->nelectrons(0))
    error("Number of upup pairs must be equal to number of spin up electrons  in pfaffian matrix");
  if(npairs(1)!=sys->nelectrons(1))
    error("Number of down-up pairs must be equal to number of spin down electrons in pfaffian matrix");
  if(npairs(2)>sys->nelectrons(0)+sys->nelectrons(1))
    error("Too many upaired orbitals in pfaffian matrix");
  if( (npairs(0)+ npairs(1)+ npairs(2))%2!=0)
    error("need even number of elements in pffafian matrix");
  
  pos=startpos;
  if (readsection(words, pos, strpfwt, "PFWT")){
    npf=strpfwt.size();
    pfwt.Resize(npf);
    for(int pf=0; pf < npf; pf++)
      pfwt(pf)=atof(strpfwt[pf].c_str());
  }
  else {
    npf=1;
    pfwt.Resize(1);
    pfwt(0)=1.0;
  }

  if(haskeyword(words, pos=startpos, "OPTIMIZE_PFWT")) {
    optimize_pfwt=1;
  }
  else optimize_pfwt=0;

  if(haskeyword(words, pos=startpos, "CHECK_PFWT_SIGN")) {
    check_pfwt_sign=1;
  }
  else check_pfwt_sign=0;
  
  if(readvalue(words, pos=0, coef_eps, "THRESHOLD")){
     if( mpi_info.node==0 )
       cout<<"Using treshold = "<<coef_eps<<" for pairing coeficients"<<endl;
  }

  vector <string> str_order;
  pos=startpos;
  int max=0;
  int npairstotal=npairs(0)+ npairs(1)+ npairs(2);
  if(!readsection(words, pos, str_order,"ORBITAL_ORDER")){
    //error("need SINGLET_ORDER in pfaffian section");
    if( mpi_info.node==0 )
      cout << "Assuming the defaul configuration of i-th orbital per i-th pfaffian "<<endl;
    order_in_pfaffian.Resize(npf);
    max=npf-1;
    for(int pf=0;pf<npf;pf++){
      order_in_pfaffian(pf).Resize(npairstotal,npairstotal);
      order_in_pfaffian(pf)=0;
      for(int i=0;i<npairstotal;i++){
	for(int j=i+1;j<npairstotal;j++){
	  order_in_pfaffian(pf)(i,j)=order_in_pfaffian(pf)(j,i)=pf;
	}
      }
    }
  }
  else {
    if (int(str_order.size())!=npairstotal*(npairstotal-1)*npf/2)
      error("supply the right amount of values in the ORBITAL_ORDER section");
    
    order_in_pfaffian.Resize(npf);
    int counter=0;
    for(int pf=0;pf<npf;pf++){
      order_in_pfaffian(pf).Resize(npairstotal,npairstotal);
      order_in_pfaffian(pf)=0;
      for(int i=0;i<npairstotal;i++){
	for(int j=i+1;j<npairstotal;j++){
	  order_in_pfaffian(pf)(i,j)=atoi(str_order[counter].c_str())-1; 
	  order_in_pfaffian(pf)(j,i)=order_in_pfaffian(pf)(i,j);
	  if( order_in_pfaffian(pf)(i,j)>max )
	    max=order_in_pfaffian(pf)(i,j);
	  counter++;
	}
      }
    }
  }
  if( mpi_info.node==0 ){ 
    cout << "ORBITAL_ORDER"<<endl;
    for(int pf=0;pf<npf;pf++){
      for(int i=0;i<npairstotal;i++){
	for(int j=0;j<npairstotal;j++){
	  if((i<npairs(0)+npairs(1) || j<npairs(0)+npairs(1)) && (i!=j))
	    cout << order_in_pfaffian(pf)(i,j)+1<<"  ";
	  else 
	    cout << order_in_pfaffian(pf)(i,j)<<"  "; 
	}
	cout<<endl;
      }
      cout <<endl;
    }

  }

    
  nsfunc=max+1;
  if( mpi_info.node==0 )
    cout <<"Number of functions needed: "<<nsfunc<<endl;

  occupation.Resize(nsfunc);
  occupation_pos.Resize(nsfunc);
  totoccupation.Resize(1);
  tripletorbuu.Resize(nsfunc);
  tripletorbdd.Resize(nsfunc);
  singletorb.Resize(nsfunc);
  unpairedorb.Resize(nsfunc);
  normalization.Resize(nsfunc);
  optimize_string.Resize(nsfunc);
  optimize_total.Resize(nsfunc);
  noccupied.Resize(nsfunc);
  optimize_pf=0;
  ntote_pairs.Resize(nsfunc);
 
  pos=startpos;
  vector <string> mowords;
  if(readsection(words,pos, mowords, "ORBITALS"))
  {
    allocate(mowords, sys, molecorb);
    is_my_mo=1;
    nmo=molecorb->getNmo();
  }
  else
  {
    single_write(cout, "Waiting to recieve the molecular orbital.\n");
  }

  //some important initializations 
  
 
  vector <vector <string> > strpairingorbs;
  vector <string> strpairingorbs_tmp;
  pos=startpos;
  while(readsection(words, pos, strpairingorbs_tmp, "PAIRING_ORBITAL")){
    strpairingorbs.push_back(strpairingorbs_tmp);
  }
  if( int(strpairingorbs.size()) < nsfunc)
    error("Could not find enough PAIRING_ORBITALs");
  for (int pf=0;pf<nsfunc;pf++){
    allocate_pairing_orbital(strpairingorbs[pf],pf); 
  }
    
  // here will start to read the particular pairing functions for each weight
  if( mpi_info.node==0 )
    cout<<"Size of pfaffian matrix is "<<npairs(0)+ npairs(1)+ npairs(2)
	<<"x"<<npairs(0)+ npairs(1)+ npairs(2)<<endl;

  nelectrons.Resize(2);

  nelectrons(0)=sys->nelectrons(0);
  nelectrons(1)=sys->nelectrons(1);

  tote=nelectrons(0)+nelectrons(1);
  ndim=3;
 

  Array1 <int> totocc_temp(nmo);
  totocc_temp=0;
  int counter=0;
  // fill the array for totoccupation
  for (int mo=0;mo<nmo;mo++){
    for (int pf=0;pf<nsfunc;pf++){
      occupation_pos(pf).Resize(occupation(pf).GetSize());
      for(int orb=0;orb<occupation(pf).GetSize();orb++){
        if(occupation(pf)(orb)==mo){
          totocc_temp(mo)=1;
        }
      }
    }
    if(totocc_temp(mo))
      counter++;
  }

  totoccupation(0).Resize(counter);
  counter=0;
  for (int mo=0;mo<nmo;mo++){
    if(totocc_temp(mo)){
      totoccupation(0)(counter)=mo;
      counter++;
    }
  }

  if( mpi_info.node==0 ){ 
    cout << "These are all orbitals involved in the Wavefunction: "<<endl;
  for(int mo=0;mo<totoccupation(0).GetSize();mo++)
    cout << totoccupation(0)(mo)+1 <<" ";
  cout << endl;
  }
  
  for (int pf=0;pf<nsfunc;pf++){
    if( mpi_info.node==0 )
      cout <<" Orbitals positions in totoccupation for Pfaffian "<<pf<<" : "<<endl;
    for(int orb=0;orb<occupation(pf).GetSize();orb++){
      for(int mo=0;mo<totoccupation(0).GetSize();mo++)
        if(totoccupation(0)(mo)==occupation(pf)(orb))
          occupation_pos(pf)(orb)=mo;
      if( mpi_info.node==0 )
	cout << occupation_pos(pf)(orb)<<" ";
    }
    if( mpi_info.node==0 )
      cout <<endl;
  }
  molecorb->buildLists(totoccupation);
}

//----------------------------------------------------------------------
void Pfaff_wf_data::allocate_pairing_orbital(vector <string> & pf_place,
                                           int orbint){

  unsigned int startpos=0;
  

  //ORBITALS_IN_PAIRING
  vector <string> strobs_pairs;
  unsigned int pos=startpos;
  if(!readsection(pf_place, pos, strobs_pairs, "ORBITALS_IN_PAIRING")){
    if( mpi_info.node==0 )
      cout << "Assuming you want to use all orbitals up to NMO in PFAFFIAN section \n";
    occupation(orbint).Resize(nmo);
     for (int i=0;i<nmo;i++){
       occupation(orbint)(i)=i;
     }
  }
  else {
    if(strobs_pairs.size()>nmo)
      error("too many orbitals requested in ORBITALS_IN_PAIRING");
    occupation(orbint).Resize(strobs_pairs.size());
    for (unsigned int i=0;i<strobs_pairs.size();i++){
      occupation(orbint)(i)=atoi(strobs_pairs[i].c_str())-1;
    }
  }

  if( mpi_info.node==0 )
    cout <<occupation(orbint).GetSize()<<" orbitals involved in pairing for "<<
      orbint+1<<" pairing orb"<<endl;

  ntote_pairs(orbint)=occupation(orbint).GetSize();
  if( mpi_info.node==0 )
    cout << "ntote_pairs("<<orbint<<") : "<<ntote_pairs(orbint)<<endl;


  //OPTIMIZE PART
  vector <string> pf_place2;
  
  pos=startpos;
  if(readsection(pf_place, pos, pf_place2, "OPTIMIZE_PF")){
    optimize_pf=1;
    optimize_string(orbint).Resize(4);
    optimize_total(orbint).Resize(4);
    Pfaffian_optimize_read(pf_place2, orbint, pos);
  }
  else {
    optimize_string(orbint).Resize(0);
    optimize_total(orbint).Resize(0);
  }
  
  vector <string> strtripletorb1,strtripletorb2;
  vector <string> strsingletorb;
  vector <string> strunpairedorb;
  vector <string> strnormalization;
  Array1 <Array1 <doublevar> > unpairedorb_temp;

  normalization(orbint).Resize(3+npairs(2));
  normalization(orbint)=1.0;  
  
  
  pos=startpos;
  normalization(orbint).Resize(3+npairs(2));

  if( readsection(pf_place,pos,strnormalization,"ALPHA")){
    if( mpi_info.node==0 )
      cout<<"ALPHA section is no longer used"<<endl;
  }

  
  pos=startpos;
  //default setting
  tripletorbuu(orbint).Resize(ntote_pairs(orbint)*(ntote_pairs(orbint)-1)/2);
  tripletorbuu(orbint)=0; 
  if(npairs(0)-npairs(1)-npairs(2)>0){
    int kounter=0;
    int kk=0;
    for(int i=0;i<ntote_pairs(orbint);i++)
      for(int j=i+1;j<ntote_pairs(orbint);j++){
	if (j==i+1){
	  if( kk>=npairs(1) && kk<(npairs(0)-npairs(2)) && !(kk-npairs(1))%2)
	    tripletorbuu(orbint)(kounter)=1.0;
	  kk++;
	}
	kounter++;
      }
  }
  if(readsection(pf_place, pos, strtripletorb1,"TRIPLET_UU_COEF")){
    if(int(strtripletorb1.size()) <  ntote_pairs(orbint)*(ntote_pairs(orbint)-1)/2){
      if( mpi_info.node==0 )
	cout <<"Will fill up upto "<< ntote_pairs(orbint)*(ntote_pairs(orbint)-1)/2 << " values in TRIPLET_UU_COEF \n";
    }
    else if (int(strtripletorb1.size()) >  ntote_pairs(orbint)*(ntote_pairs(orbint)-1)/2){
      error("Too many TRIPLET_UU_COEF in PFAFFIAN section");
    }
    for(unsigned int i=0; i<strtripletorb1.size(); i++)
      {
	tripletorbuu(orbint)(i)=atof(strtripletorb1[i].c_str());
      }
  }
  
  Renormalization(tripletorbuu(orbint), normalization(orbint)(0), coef_eps);
  if (npairs(1)) {
    pos=startpos;
    //default setting
    tripletorbdd(orbint).Resize(ntote_pairs(orbint)*(ntote_pairs(orbint)-1)/2);
    tripletorbdd(orbint)=0; 
    if(readsection(pf_place,pos,strtripletorb2,"TRIPLET_DD_COEF")){
      if(int(strtripletorb2.size()) <  ntote_pairs(orbint)*(ntote_pairs(orbint)-1)/2){
	if( mpi_info.node==0 )
	  cout <<"Will fill up upto "<< ntote_pairs(orbint)*(ntote_pairs(orbint)-1)/2 << " values in TRIPLET_DD_COEF \n";
    }
      else if (int(strtripletorb2.size()) >  ntote_pairs(orbint)*(ntote_pairs(orbint)-1)/2){
      error("Too many TRIPLET_DD_COEF in PFAFFIAN section");
    }
    for(unsigned int i=0; i<strtripletorb2.size(); i++)
      {
	tripletorbdd(orbint)(i)=atof(strtripletorb2[i].c_str());
      }
    }
    Renormalization(tripletorbdd(orbint), normalization(orbint)(1),coef_eps);
    //default setting
    pos=startpos;
    singletorb(orbint).Resize(ntote_pairs(orbint)*(ntote_pairs(orbint)+1)/2);
    singletorb(orbint)=0;
    int kounter=0;
    int kk=0;
    for(int i=0;i<ntote_pairs(orbint);i++)
      for(int j=i;j<ntote_pairs(orbint);j++){
	if (i==j && kk< npairs(1)){
	  singletorb(orbint)(kounter)=1.0;
	  kk++;
	}
	kounter++;
      }
    if(readsection(pf_place,pos,strsingletorb,"SINGLET_COEF")){
      if(int(strsingletorb.size()) < ntote_pairs(orbint)*(ntote_pairs(orbint)+1)/2)
	if( mpi_info.node==0 )
	  cout <<"Will fill up upto "<< ntote_pairs(orbint)*(ntote_pairs(orbint)+1)/2 << " values in SINGLET_COEF \n";
	else if (int(strsingletorb.size()) >  ntote_pairs(orbint)*(ntote_pairs(orbint)+1)/2){
	  error("Too many SINGLET_COEF in PFAFFIAN section");
        }
      for(unsigned int i=0; i<strsingletorb.size(); i++)
	singletorb(orbint)(i)=atof(strsingletorb[i].c_str());
    }
    Renormalization(singletorb(orbint), normalization(orbint)(2),coef_eps);
  }
  else {
    tripletorbdd(orbint).Resize(0);
    singletorb(orbint).Resize(0);
  }
  if (npairs(2)) {
    pos=startpos;
    unpairedorb(orbint).Resize(ntote_pairs(orbint)*npairs(2));
    unpairedorb(orbint)=0;
    
    if(!readsection(pf_place,pos,strunpairedorb,"UNPAIRED_COEF")){
      error("Need UNPAIRED_COEF in PFAFFIAN section");}
    if(int(strunpairedorb.size()) < ntote_pairs(orbint)*npairs(2)){
      if( mpi_info.node==0 ) 
	cout <<"Will fill up upto "<< ntote_pairs(orbint)*npairs(2) << " values in UNPAIRED_COEF \n";
    }
    else if (int(strunpairedorb.size()) >  ntote_pairs(orbint)*npairs(2)){
      error("Too many UNPAIRED_COEF's in PFAFFIAN section");
    }
    unpairedorb_temp.Resize(npairs(2));
    for(int i=0; i< npairs(2); i++)
      {
        unpairedorb_temp(i).Resize(ntote_pairs(orbint));
        for(int j=0; j<ntote_pairs(orbint); j++)  
              unpairedorb_temp(i)(j)=atof(strunpairedorb[i*ntote_pairs(orbint)+j].c_str());
        Renormalization(unpairedorb_temp(i), normalization(orbint)(3+i),coef_eps);//we have to normalize for every row
        for(int j=0; j<ntote_pairs(orbint); j++){
          unpairedorb(orbint)(i*ntote_pairs(orbint)+j)=unpairedorb_temp(i)(j);
	}
      }
  }
  else {
    unpairedorb(orbint).Resize(0);
  }
  if(optimize_total(orbint).GetSize()){
    if(optimize_string(orbint)(0)!="" &&  normalization(orbint)(0)<coef_eps)
      error("Optimization of TRIPLET_UU_COEF when identical to zero is not supported!");
    if(optimize_string(orbint)(1)!="" &&  normalization(orbint)(1)<coef_eps)
      error("Optimization of TRIPLET_DD_COEF when identical to zero is not supported!");
    if(optimize_string(orbint)(2)!="" &&  normalization(orbint)(2)<coef_eps)
      error("Optimization of SINGLET_COEF when identical to zero is not supported!");
    for(int i=0; i< npairs(2); i++)
      if(optimize_string(orbint)(3)!="" &&  normalization(orbint)(3+i)<coef_eps)
	error("Optimization of UNPAIRED_COEF when identical to zero is not supported!");
  }

}

//----------------------------------------------------------------------
void Pfaff_wf_data::Pfaffian_optimize_read(vector <string> & pf_place,
                                           int pf, unsigned int startpos){
  
  
  unsigned int pos=startpos;


  if(!readvalue(pf_place, pos=0, noccupied(pf), "NOCCUPIED")){
    noccupied(pf)=npairs(0);
  }
  if( mpi_info.node==0 )
    cout << "NOCCUPIED : " <<noccupied(pf)<<endl;
  
  optimize_string(pf)(0)="";
  optimize_string(pf)(1)="";
  optimize_string(pf)(2)="";
  optimize_string(pf)(3)="";

    
  optimize_total(pf)(0).Resize(((ntote_pairs(pf)-1)*ntote_pairs(pf))/2);
  optimize_total(pf)(1).Resize(((ntote_pairs(pf)-1)*ntote_pairs(pf))/2);
  optimize_total(pf)(2).Resize(((ntote_pairs(pf)+1)*ntote_pairs(pf))/2);
  optimize_total(pf)(3).Resize(ntote_pairs(pf)*npairs(2));
  
  
  for (int i=0;i<optimize_total(pf).GetSize();i++)
    for (int j=0;j<optimize_total(pf)(i).GetSize();j++)
      optimize_total(pf)(i)(j)=0;

  if( mpi_info.node==0 )
    cout <<"Pfaffian parameters to optimize: "<<endl;
  int counter=0;
  int fullcounter=0;
  
  //------------Triplets up up------------------------
  vector <string> strtripletuupos;
  if(haskeyword(pf_place, pos=0, "TRIPLET_UU_DIAG")){
      counter=fullcounter=0;
      for (int i=0;i<ntote_pairs(pf);i++){
        for (int j=i+1;j<ntote_pairs(pf);j++){
          if((j==i+1)){//&&(i%2==0)){
            optimize_total(pf)(0)(fullcounter)=1;
            counter++;
          }
          fullcounter++;
        }
      }
      optimize_string(pf)(0)="TRIPLET_UU_DIAG";
      if( mpi_info.node==0 )
	cout <<"   "<< "TRIPLET_UU_DIAG: "<<counter<<endl;
    }
    else if(haskeyword(pf_place, pos=0,"TRIPLET_UU_ALL")){
      counter=fullcounter=0;
      for (int i=0;i<ntote_pairs(pf);i++){
        for (int j=i+1;j<ntote_pairs(pf);j++){
          //if(i%2==0){
	    optimize_total(pf)(0)(fullcounter)=1;
	    counter++;
	  //}
	  fullcounter++;
        }
      }
      optimize_string(pf)(0)="TRIPLET_UU_ALL";
      if( mpi_info.node==0 )
	cout <<"   "<< "TRIPLET_UU_ALL: "<<counter<<endl;
    }
    else if(haskeyword(pf_place, pos=0, "TRIPLET_UU_VIRTUAL")){
      counter=fullcounter=0;
      for (int i=0;i<ntote_pairs(pf);i++){
        for (int j=i+1;j<ntote_pairs(pf);j++){
          if((i>=noccupied(pf))){//&&((i-noccupied(pf))%2==0)){
            optimize_total(pf)(0)(fullcounter)=1;
            counter++;
          }
          fullcounter++;
        }
      }
      optimize_string(pf)(0)="TRIPLET_UU_VIRTUAL";
      if( mpi_info.node==0 )
	cout <<"   "<< "TRIPLET_UU_VIRTUAL: "<<counter<<endl;
    }
    else if(haskeyword(pf_place, pos=0, "TRIPLET_UU_VIRTUAL_DIAG")){
      counter=fullcounter=0;
      for (int i=0;i<ntote_pairs(pf);i++){
        for (int j=i+1;j<ntote_pairs(pf);j++){
          if((i>=noccupied(pf))&&(j==i+1)){//&&((i-noccupied(pf))%2==0)){
            optimize_total(pf)(0)(fullcounter)=1;
            counter++;
          }
          fullcounter++;
        }
      }      
      optimize_string(pf)(0)="TRIPLET_UU_HF2VIRTUALS";
      if( mpi_info.node==0 )
	cout <<"   "<< "TRIPLET_UU_VIRTUAL_DIAG: "<<counter<<endl;
    }
    else if(haskeyword(pf_place, pos=0, "TRIPLET_UU_HF2VIRTUALS")){
      counter=fullcounter=0;
      for (int i=0;i<ntote_pairs(pf);i++){
        for (int j=i+1;j<ntote_pairs(pf);j++){
          if((i>= npairs(1))&&(i< npairs(0))&&(j>=npairs(0))){
            optimize_total(pf)(0)(fullcounter)=1;
            counter++;
          }
          fullcounter++;
        }
      }      
      optimize_string(pf)(0)="TRIPLET_UU_HF2VIRTUALS";
      if( mpi_info.node==0 )
	cout <<"   "<< "TRIPLET_UU_HF2VIRTUALS: "<<counter<<endl;
    }
    else if(readsection(pf_place, pos=0, strtripletuupos, "TRIPLET_UU_FROM_TO")){
      if(strtripletuupos.size()!=2)
        error("expected 2 values in TRIPLET_UU_FROM_TO");
      //cout <<strsingletpos[0]<<"  "<<strsingletpos[1]<<endl;
      Array1 <int> tripletpos(2);
      tripletpos(0)=atoi(strtripletuupos[0].c_str())-1;
      tripletpos(1)=atoi(strtripletuupos[1].c_str())-1;

      if(tripletpos(0)>=tripletpos(1))
        error("tripletpos(0)>=tripletpos(1) in TRIPLET_DD_FROM_TO");
      if(tripletpos(0)<0 || tripletpos(1)>=ntote_pairs(pf))
        error("tripletpos(0)<0 || tripletpos(1)>=ntote_pairs(pf) in TRIPLET_DD_FROM_TO");
      
      //if( mpi_info.node==0 )
      //	cout <<"   "<< "SINGLET_FROM "<<singletpos(0)<<" TO "<<singletpos(1)<<endl;
      
      counter=fullcounter=0;
      for (int i=0;i<ntote_pairs(pf);i++){
        for (int j=i+1;j<ntote_pairs(pf);j++){
          if(i>=tripletpos(0) && i<=tripletpos(1))
            if(j>=tripletpos(0) && j<=tripletpos(1)){
              optimize_total(pf)(0)(fullcounter)=1;
              counter++;
            }
          fullcounter++;
        }
      }

      optimize_string(pf)(0)="TRIPLET_UU_FROM_TO { "+strtripletuupos[0]+"  "+strtripletuupos[1]+" }";
      //cout << optimize_string(pf)(2)<<endl;
      if( mpi_info.node==0 )
        cout <<"   "<< "TRIPLET_UU_FROM_TO: "<<counter<<endl;
    }   

//------------Triplets down down------------------------
    vector <string> strtripletddpos;
    if(haskeyword(pf_place, pos=0, "TRIPLET_DD_DIAG")){
      counter=fullcounter=0;
      for (int i=0;i<ntote_pairs(pf);i++){
        for (int j=i+1;j<ntote_pairs(pf);j++){
          if((j==i+1)){//&&(i%2==0)){
            optimize_total(pf)(1)(fullcounter)=1;
            counter++;
          }
          fullcounter++;
        }
      }
      optimize_string(pf)(1)="TRIPLET_DD_DIAG";
      if( mpi_info.node==0 )
	cout <<"   "<< "TRIPLET_DD_DIAG: "<<counter<<endl;
    }
    else if(haskeyword(pf_place, pos=0,"TRIPLET_DD_ALL")){
      counter=fullcounter=0;
      for (int i=0;i<ntote_pairs(pf);i++){
        for (int j=i+1;j<ntote_pairs(pf);j++){
          //if(i%2==0){
	    optimize_total(pf)(1)(fullcounter)=1;
	    counter++;
	  //}
	  fullcounter++;
        }
      }
      optimize_string(pf)(1)="TRIPLET_DD_ALL";
      if( mpi_info.node==0 )
	cout <<"   "<< "TRIPLET_DD_ALL: "<<counter<<endl;
    }
    else if(haskeyword(pf_place, pos=0, "TRIPLET_DD_VIRTUAL")){
      counter=fullcounter=0;
      for (int i=0;i<ntote_pairs(pf);i++){
        for (int j=i+1;j<ntote_pairs(pf);j++){
          if((i>=noccupied(pf))&&((i-noccupied(pf))%2==0)){
            optimize_total(pf)(1)(fullcounter)=1;
            counter++;
          }
          fullcounter++;
        }
      }
      optimize_string(pf)(1)="TRIPLET_DD_VIRTUAL";
      if( mpi_info.node==0 )
	cout <<"   "<< "TRIPLET_DD_VIRTUAL: "<<counter<<endl;
    }
    else if(haskeyword(pf_place, pos=0, "TRIPLET_DD_VIRTUAL_DIAG")){
      counter=fullcounter=0;
      for (int i=0;i<ntote_pairs(pf);i++){
        for (int j=i+1;j<ntote_pairs(pf);j++){
          if((i>=noccupied(pf))&&(j==i+1)){//&&((i-noccupied(pf))%2==0)){
            optimize_total(pf)(1)(fullcounter)=1;
            counter++; 
          }
          fullcounter++;
        }
      }      
      optimize_string(pf)(1)="TRIPLET_DD_VIRTUAL_DIAG";
      if( mpi_info.node==0 )
	cout <<"   "<< "TRIPLET_DD_VIRTUAL_DIAG: "<<counter<<endl;
    }
    else if(readsection(pf_place, pos=0, strtripletddpos, "TRIPLET_DD_FROM_TO")){
      if(strtripletddpos.size()!=2)
        error("expected 2 values in TRIPLET_DD_FROM_TO");
      //cout <<strsingletpos[0]<<"  "<<strsingletpos[1]<<endl;
      Array1 <int> tripletpos(2);
      tripletpos(0)=atoi(strtripletddpos[0].c_str())-1;
      tripletpos(1)=atoi(strtripletddpos[1].c_str())-1;

      if(tripletpos(0)>=tripletpos(1))
        error("tripletpos(0)>=tripletpos(1) in TRIPLET_DD_FROM_TO");
      if(tripletpos(0)<0 || tripletpos(1)>=ntote_pairs(pf))
        error("tripletpos(0)<0 || tripletpos(1)>=ntote_pairs(pf) in TRIPLET_DD_FROM_TO");
      
      //if( mpi_info.node==0 )
      //	cout <<"   "<< "SINGLET_FROM "<<singletpos(0)<<" TO "<<singletpos(1)<<endl;
      
      counter=fullcounter=0;
      for (int i=0;i<ntote_pairs(pf);i++){
        for (int j=i+1;j<ntote_pairs(pf);j++){
          if(i>=tripletpos(0) && i<=tripletpos(1))
            if(j>=tripletpos(0) && j<=tripletpos(1)){
              optimize_total(pf)(1)(fullcounter)=1;
              counter++;
            }
          fullcounter++;
        }
      }

      optimize_string(pf)(1)="TRIPLET_DD_FROM_TO { "+strtripletddpos[0]+"  "+strtripletddpos[1]+" }";
      //cout << optimize_string(pf)(2)<<endl;
      if( mpi_info.node==0 )
        cout <<"   "<< "TRIPLET_DD_FROM_TO: "<<counter<<endl;
    }
//------------Singlets------------------------

    vector <string> strsingletpos;
    if(haskeyword(pf_place, pos=0,"SINGLET_DIAG")){
      counter=fullcounter=0;
      for (int i=0;i<ntote_pairs(pf);i++){
        for (int j=i;j<ntote_pairs(pf);j++){
          if((j==i)){
            optimize_total(pf)(2)(fullcounter)=1;
            counter++;
          }
          fullcounter++;
        }
      }
      optimize_string(pf)(2)="SINGLET_DIAG";
      if( mpi_info.node==0 )
	cout <<"   "<< "SINGLET_DIAG: "<<counter<<endl;
    }
    else if(haskeyword(pf_place, pos=0,"SINGLET_ALL")){
      counter=0;
      for (int i=0;i<ntote_pairs(pf);i++){
        for (int j=i;j<ntote_pairs(pf);j++){
          optimize_total(pf)(2)(counter)=1;
          counter++;
        }
      }
      optimize_string(pf)(2)="SINGLET_ALL";
      if( mpi_info.node==0 )
      cout <<"   "<< "SINGLET_ALL: "<<counter<<endl;
    }
    else if(readsection(pf_place, pos=0, strsingletpos, "SINGLET_FROM_TO")){
      if(strsingletpos.size()!=2)
        error("expected 2 values in SINGLET_FROM_TO");
      //cout <<strsingletpos[0]<<"  "<<strsingletpos[1]<<endl;
      Array1 <int> singletpos(2);
      singletpos(0)=atoi(strsingletpos[0].c_str())-1;
      singletpos(1)=atoi(strsingletpos[1].c_str())-1;

      if(singletpos(0)>=singletpos(1))
        error("singletpos(0)>=singletpos(1) in SINGLET_FROM_TO");
      if(singletpos(0)<0 || singletpos(1)>=ntote_pairs(pf))
        error("singletpos(0)<0 || singletpos(1)>=ntote_pairs(pf) in SINGLET_FROM_TO");
      
      //if( mpi_info.node==0 )
      //	cout <<"   "<< "SINGLET_FROM "<<singletpos(0)<<" TO "<<singletpos(1)<<endl;
      
      counter=fullcounter=0;
      for (int i=0;i<ntote_pairs(pf);i++){
        for (int j=i;j<ntote_pairs(pf);j++){
          if(i>=singletpos(0) && i<=singletpos(1))
            if(j>=singletpos(0) && j<=singletpos(1)){
              optimize_total(pf)(2)(fullcounter)=1;
              counter++;
            }
          fullcounter++;
        }
      }

      optimize_string(pf)(2)="SINGLET_FROM_TO { "+strsingletpos[0]+"  "+strsingletpos[1]+" }";
      //cout << optimize_string(pf)(2)<<endl;
      if( mpi_info.node==0 )
        cout <<"   "<< "SINGLET_FROM_TO: "<<counter<<endl;
    }
    /*
    else if(haskeyword(pf_place, pos=0,"SINGLET_DIAG_SIMPLE")){
      optimize_string(pf)(2)(2)=npairs(1)+1;
        cout <<"   "<< "SINGLET_DIAG_SIMPLE: "<<optimize_string(pf)(2)(2)<<endl;
    }
    else if(haskeyword(pf_place, pos=0,"SINGLET_DIAG_XYZ")){
      optimize_string(pf)(2)(3)=3;
      cout <<"   "<< "SINGLET_DIAG_XYZ: "<<optimize_string(pf)(2)(3)<<endl;
      }
    else if(haskeyword(pf_place, pos=0,"SINGLET_XYZ")){
      optimize_string(pf)(2)(4)=6;
        cout <<"   "<< "SINGLET_XYZ: "<<optimize_string(pf)(2)(4)<<endl;
    }
    */
    else if(haskeyword(pf_place, pos=0, "SINGLET_VIRTUAL")){
      counter=fullcounter=0;
      for (int i=0;i<ntote_pairs(pf);i++){
        for (int j=i;j<ntote_pairs(pf);j++){
          if((i>=noccupied(pf))){
            optimize_total(pf)(2)(fullcounter)=1;
            counter++;
          }
          fullcounter++;
        }
      }
      optimize_string(pf)(2)="SINGLET_VIRTUAL";
      if( mpi_info.node==0 )
	cout <<"   "<< "SINGLET_VIRTUAL: "<<counter<<endl;
    }
    else if(haskeyword(pf_place, pos=0, "SINGLET_VIRTUAL_DIAG")){
      counter=0;
      for (int i=0;i<ntote_pairs(pf);i++){
        for (int j=i;j<ntote_pairs(pf);j++){
          if((i>=noccupied(pf))&&(j==i)){
            optimize_total(pf)(2)(fullcounter)=1;
            counter++;
          }
          fullcounter++;
        }
      }      
      optimize_string(pf)(2)="SINGLET_VIRTUAL_DIAG";
      if( mpi_info.node==0 )
	cout <<"   "<< "SINGLET_VIRTUAL_DIAG: "<<counter<<endl;
    }
    else if(haskeyword(pf_place, pos=0, "SINGLET_SINGLES_AND_DOUBLES")){
      counter=0;
      for (int i=0;i<ntote_pairs(pf);i++){
        for (int j=i;j<ntote_pairs(pf);j++){
          if((i<noccupied(pf))||(i==j)){
            optimize_total(pf)(2)(fullcounter)=1;
            counter++;
          }
          fullcounter++;
        }
      }      
      optimize_string(pf)(2)="SINGLET_SINGLES_AND_DOUBLES";
      if( mpi_info.node==0 )
	cout <<"   "<< "SINGLET_SINGLES_AND_DOUBLES: "<<counter<<endl;
    }
    else if(haskeyword(pf_place, pos=0, "SINGLET_HF2VIRTUALS")){
      counter=fullcounter=0;
      for (int i=0;i<ntote_pairs(pf);i++){
        for (int j=i;j<ntote_pairs(pf);j++){
          if((i< npairs(1))&&(j>= npairs(1))){
            optimize_total(pf)(2)(fullcounter)=1;
            counter++;
          }
          fullcounter++;
        }
      }      
      optimize_string(pf)(2)="SINGLET_HF2VIRTUALS";
      if( mpi_info.node==0 )
	cout <<"   "<< "SINGLET_HF2VIRTUALS: "<<counter<<endl;
    }

//------------Unpaired------------------------
    if(haskeyword(pf_place, pos=0,"UNPAIRED_ALL")){
      counter=0;
      for(int i=0;i<ntote_pairs(pf);i++)
        for(int j=0;j<npairs(2);j++){
          optimize_total(pf)(3)(counter)=1; 
          counter++;
        }
      optimize_string(pf)(3)="UNPAIRED_ALL";
      if( mpi_info.node==0 )
	cout <<"   "<< "UNPAIRED_ALL: "<<counter<<endl;
    }
//------------Alpha------------------------
    if(haskeyword(pf_place, pos=0,"ALPHA")){
      error("ALPHA section of pairing orbital can no longer be optimized");
    }
  
}


//----------------------------------------------------------------------
int Pfaff_wf_data::supports(wf_support_type support) {
  switch(support) {
  case laplacian_update:
    return 1;
  case density:
    return 1;
  case parameter_derivatives:
    return 1;
  default:
    return 0;
  }
}
//----------------------------------------------------------------------

void Pfaff_wf_data::init_mo() {
  
  molecorb->buildLists(totoccupation);

}

//----------------------------------------------------------------------

int Pfaff_wf_data::useMo(MO_matrix * mo) {
  assert(mo);
  
  if(molecorb) delete molecorb;

  is_my_mo=0;
  molecorb=mo;
  init_mo();
  return 1;
}

//----------------------------------------------------------------------

void Pfaff_wf_data::generateWavefunction(Wavefunction *& wf)
{
  assert(wf==NULL);
  if(!molecorb)
    error("Pfaff_wf_data::Need to allocate molecular orbital before generating any Wavefunction objects");
  wf=new Pfaff_wf;
  Pfaff_wf * pfaffwf;
  recast(wf, pfaffwf);
  pfaffwf->init(this);
  attachObserver(pfaffwf);
}

int Pfaff_wf_data::showinfo(ostream & os)
{
  
  string indent="     ";
  if(!molecorb)
    error("Pfaff_wf_data::showinfo() : Molecular orbital not allocated");

  os << "Pfaffian " << endl;

  os <<"Number of pairs: ";
  for (int i=0;i< npairs.GetSize();i++)
    os <<npairs(i)<<" ";
  //  os << "." << endl;
  os <<endl;

  os << "Pffafian weights:  ";
  for (int pf=0;pf< npf;pf++)
    os <<pfwt(pf)<<"  ";
  os << endl << endl;

  os << "Order of pairing orbital functions per pfaffian"<<endl;
  for(int pf=0;pf<npf;pf++){
    os << "Pfaffian: " <<pf+1<<endl;
    for(int i=0;i<order_in_pfaffian(pf).GetDim(0);i++){
      for(int j=0;j<order_in_pfaffian(pf).GetDim(1);j++){
	char tmp[10];
	char tmp2[10];
	sprintf(tmp, "(%d,%d) ", i+1,j+1);
	sprintf(tmp2, "(%d)(%d) ", i+1,j+1);
	if((i<npairs(0)+npairs(1) || j<npairs(0)+npairs(1)) && (i!=j)){
	  if(i<npairs(0)){
	    if(j<npairs(0))
	      os << " Chi_"<<order_in_pfaffian(pf)(i,j)+1<<tmp;
	    else if(j<npairs(0)+npairs(1))
	      os << " Phi_"<<order_in_pfaffian(pf)(i,j)+1<<tmp;
	    else 
	      os << " varphi_"<<order_in_pfaffian(pf)(i,j)+1<<tmp2;
	  }
	  else if(i<npairs(0)+npairs(1)){
	    if(j<npairs(0))
	      os << "-Phi_"<<order_in_pfaffian(pf)(i,j)+1<<tmp;
	    else if(j<npairs(0)+npairs(1))
	      os << " Chi_"<<order_in_pfaffian(pf)(i,j)+1<<tmp;
	    else
	      os << " varphi_"<<order_in_pfaffian(pf)(i,j)+1<<tmp2;
	  }
	  else{
	    if(j<npairs(0))
	      os << "-varphi_"<<order_in_pfaffian(pf)(i,j)+1<<tmp2;
	    else
	       os << "-varphi_"<<order_in_pfaffian(pf)(i,j)+1<<tmp2;
	  }
	  
	}
	else 
	  os <<"     "<<order_in_pfaffian(pf)(i,j)<<"      "; 
      }
      os <<endl;
    }
    os <<endl;
  }
  
  for (int pf=0;pf<nsfunc;pf++){
    
    os <<"Pairing Orbital "<<pf+1<<endl;
    
    os <<indent<< "Orbitals involved in pairing: \n";
    os <<indent;
    for (int i=0;i< occupation(pf).GetSize();i++)
      os <<occupation(pf)(i)+1<<"  ";
    //os << "." << endl;
    os << endl<<endl;

  
     os <<indent<< "Triplet up up coeficients: " << endl;
     int count=0;
     // os.setf(ios::showpoint|ios::scientific);
     for(int i=0;i<ntote_pairs(pf)-1;i++){
       os <<indent;
       for(int j=0;j<ntote_pairs(pf);j++){
         if(j>i){
           os<<tripletorbuu(pf)(count)<< " ";
           count++;
         }
         //else 
         //  os << "                      ";
       }
       os <<endl;
     }
     
     //os.unsetf(ios::showpoint|ios::scientific);
     os <<endl;
     
     if (npairs(1)){
       os <<indent<< "Triplet down down coeficients: " << endl;
       count=0;
       // os.setf(ios::showpoint|ios::scientific);
       for(int i=0;i<ntote_pairs(pf)-1;i++){
         os <<indent;
         for(int j=0;j<ntote_pairs(pf);j++){
           if(j>i){
             os<<tripletorbdd(pf)(count)<< " ";
             count++;
           }
           //else 
           //  os << "                      ";
         }
         os <<endl;
       }
       os <<endl;

       
       count=0;
       os <<indent<<"Singlet coeficients: " << endl;
       //os.setf(ios::showpoint|ios::scientific);
       for(int i=0;i<ntote_pairs(pf);i++){
         os <<indent;
         for(int j=0;j<ntote_pairs(pf);j++){
           if(j>=i){
             os<<singletorb(pf)(count)<< " ";
             count++;
           }
           //  else 
           //  os << "                      ";
         }
         os <<endl;
       }
       //os.unsetf(ios::showpoint|ios::scientific);
       os <<endl;
     }

     if (npairs(2)){
       count=0;
       os <<indent<<"Unpaired coeficients: " << endl;
       //os.setf(ios::showpoint|ios::scientific);
       for(int i=0;i< npairs(2);i++){
         os <<indent;
         for(int j=0;j<ntote_pairs(pf);j++){
        os<<unpairedorb(pf)(count)<< " ";
        count++;
         }
         os <<endl;
       }
       //os.unsetf(ios::showpoint|ios::scientific);
       os <<endl;
     }
  }
  os << "Molecular Orbital object : ";
  molecorb->showinfo(os);
  
  return 1;
}

//----------------------------------------------------------------------

int Pfaff_wf_data::writeinput(string & indent, ostream & os)
{

  if(!molecorb)
    error("Pfaff_wf_data::writeinput() : Molecular orbital not allocated");
  
  os << indent << "PFAFFIAN" << endl;
  os << indent << "NPAIRS {  ";
  for (int i=0;i< npairs.GetSize();i++)
    os <<npairs(i)<<"  ";
  os << "}" << endl;
  
  os << indent << "PFWT {  ";
  for (int pf=0;pf< npf;pf++)
    os <<pfwt(pf)<<"  ";
  os << "}" << endl;
  if (optimize_pfwt)
    os << indent << "OPTIMIZE_PFWT" <<endl;

  os << indent << "ORBITAL_ORDER {  "<<endl;
  for(int pf=0;pf<npf;pf++){
    for(int i=0;i<order_in_pfaffian(pf).GetDim(0);i++){
      for(int j=i+1;j<order_in_pfaffian(pf).GetDim(1);j++){
	os << indent << order_in_pfaffian(pf)(i,j)+1 <<"  ";
      }
      os << indent << endl;
    }
    os << indent << endl;
  }
  os << indent << "}" << endl;
    
  
  for(int pf=0;pf<nsfunc;pf++){
    os << indent << "PAIRING_ORBITAL {  "<<endl;
    
    os << indent << indent << "ORBITALS_IN_PAIRING {  ";
    for (int i=0;i< occupation(pf).GetSize();i++)
      os <<occupation(pf)(i)+1<<"  ";
    os << "}" << endl;

    if(optimize_pf) {
      if(optimize_string(pf).GetSize()){
        os << indent << indent << "OPTIMIZE_PF {" << endl;
        os << indent << indent << indent << "NOCCUPIED " <<noccupied(pf)<<endl;
        for(int i=0;i<optimize_string(pf).GetSize();i++){
          if (optimize_string(pf)(i)!="")
            os << indent << indent << indent << optimize_string(pf)(i) << endl;
        }
        os << indent << indent << "}" << endl;
      }
    }
    
    
    
    os <<  indent<< indent << "TRIPLET_UU_COEF {" << endl;
    Renormalization(tripletorbuu(pf), normalization(pf)(0),coef_eps);
    int count=0;
    //  os.setf(ios::showpoint|ios::scientific);
    for(int i=0;i<ntote_pairs(pf)-1;i++){
      os << indent << indent;
      for(int j=0;j<ntote_pairs(pf);j++){
        if(j>i){
          os<<tripletorbuu(pf)(count)<< " ";
          count++;
        }
      // else 
        //  os << "                      ";
      }
      os <<endl;
    }

    //os.unsetf(ios::showpoint|ios::scientific);
    os << indent <<indent << "}"<<endl;

    if (npairs(1)){
      os << indent << indent << "TRIPLET_DD_COEF {" << endl;
      Renormalization(tripletorbdd(pf), normalization(pf)(1),coef_eps);
      count=0;
      //  os.setf(ios::showpoint|ios::scientific);
      for(int i=0;i<ntote_pairs(pf)-1;i++){
        os << indent << indent;
        for(int j=0;j<ntote_pairs(pf);j++){
          if(j>i){
            os<<tripletorbdd(pf)(count)<< " ";
            count++;
          }
          // else 
          //  os << "                      ";
        }
        os <<endl;
      }
    
      //os.unsetf(ios::showpoint|ios::scientific);
      os << indent << indent << "}"<<endl;
  
      
      
      count=0;
      os << indent << indent<< "SINGLET_COEF {" << endl;
      Renormalization(singletorb(pf), normalization(pf)(2),coef_eps);
      //os.setf(ios::showpoint|ios::scientific);
      for(int i=0;i<ntote_pairs(pf);i++){
        os<< indent << indent;
        for(int j=0;j<ntote_pairs(pf);j++){
          if(j>=i){
            os<<singletorb(pf)(count)<< " ";
            count++;
          }
          // else 
          // os << "                      ";
        }
        os <<endl;
      }
      //os.unsetf(ios::showpoint|ios::scientific);
      os << indent << indent << "}"<<endl;
    }
    
    if (npairs(2)){
      count=0;
      os << indent << indent << "UNPAIRED_COEF {" << endl;
      Array1 <Array1 <doublevar> > unpairedorb_temp(npairs(2));
      for(int i=0; i<npairs(2); i++){
	unpairedorb_temp(i).Resize(ntote_pairs(pf));
	for(int j=0; j<ntote_pairs(pf); j++)  
	  unpairedorb_temp(i)(j)=unpairedorb(pf)(i*ntote_pairs(pf)+j);
	Renormalization(unpairedorb_temp(i), normalization(pf)(3+i),coef_eps);//we have to normalize for every row
	for(int j=0; j<ntote_pairs(pf); j++)
	  unpairedorb(pf)(i*ntote_pairs(pf)+j)=unpairedorb_temp(i)(j);
      }
      //os.setf(ios::showpoint|ios::scientific);
      for(int i=0;i<npairs(2);i++){
        os<< indent << indent ;
        for(int j=0;j<ntote_pairs(pf);j++){
          os<<unpairedorb(pf)(count)<< " ";
          count++;
        }
        os <<endl;
      }
      //os.unsetf(ios::showpoint|ios::scientific);
      os << indent << indent << "}"<<endl;
    }
    os << indent << "}"<<endl;
  }

  os << indent << "ORBITALS {\n";
  string indent2=indent+"  ";
  molecorb->writeinput(indent2, os);
  os << indent << "}\n";

  return 1;
}

void Renormalization(Array1 <doublevar> & coef, doublevar & norm, doublevar coef_eps)
{
  doublevar tmp=0.0;
  doublevar norm2=0.0;

  //cout <<" coef_eps "<<coef_eps<<endl;
  
  for (int i=0;i< coef.GetSize();i++)
    tmp+=coef(i)*coef(i);
  norm=sqrt(tmp);

  if ( norm > coef_eps ) {
    norm2=1.0/norm;
    norm=1.0;
    for (int i=0;i<coef.GetSize();i++){
      if (abs(coef(i))>coef_eps)
	coef(i)*=norm2;
      else 
	coef(i)=0.0;
    }
  }
  else {
    norm=0.0;
    coef=0.0;
  }
  return; 
}

void Getnormalization(Array1 <doublevar> & coef, doublevar & norm)
{
  doublevar tmp=0.0;
  for (int i=0;i< coef.GetSize();i++)
    tmp+=coef(i)*coef(i);
  norm=sqrt(tmp);
  return; 
}




int Pfaff_wf_data::nparms(){
  if(optimize_pf || optimize_pfwt) {
    
    int counter=0;
    for(int pf=0;pf<nsfunc;pf++)
      for(int k=0;k<optimize_total(pf).GetSize();k++)
        for(int l=0;l<optimize_total(pf)(k).GetSize();l++){
          if(optimize_total(pf)(k)(l)){
            counter++;
          }
        }
    if (optimize_pfwt){
      for(int pf=1; pf< pfwt.GetSize(); pf++) {
        counter++;
        }
    }
    return counter;
  }
  else 
    return 0;
}


void Pfaff_wf_data::getPFCoeff(Array1 <doublevar> & parms){
  //cout <<"Start getPFCoeff"<<endl;
  int counter=0;
  for(int pf=0;pf<nsfunc;pf++)
    for(int k=0;k<optimize_total(pf).GetSize();k++)
      for(int l=0;l<optimize_total(pf)(k).GetSize();l++){
        if(optimize_total(pf)(k)(l)){
          switch(k){
          case 0:
            parms(counter)=tripletorbuu(pf)(l);
            break;
          case 1:
            parms(counter)=tripletorbdd(pf)(l); 
            break;
          case 2:
            parms(counter)=singletorb(pf)(l);
            break;
          case 3: 
            parms(counter)=unpairedorb(pf)(l);
            break;
	  }
          counter++;
        }
      }
  if (optimize_pfwt){
    for(int pf=1; pf< pfwt.GetSize(); pf++) {
      parms(counter)=pfwt(pf);
      counter++;
      }
  }
  
  // if(mpi_info.node==0) {
  //  cout<<"Get parameters values"<<endl;
  // for (int i=0;i<parms.GetSize();i++)
  //    cout << i+1<<"   "<<parms(i)<<endl;
  //}
  
}

void Pfaff_wf_data::setPFCoeff(Array1 <doublevar> & parms){
  // cout <<"Start setPFCoeff"<<endl;
  int counter=0;
  Array1 <Array1 <doublevar> > unpairedorb_temp(npairs(2));

  for(int pf=0;pf<nsfunc;pf++){
    for(int k=0;k<optimize_total(pf).GetSize();k++){
      for(int l=0;l<optimize_total(pf)(k).GetSize();l++){
        if(optimize_total(pf)(k)(l)){
          switch(k){
          case 0:
            tripletorbuu(pf)(l)=parms(counter);
            break;
          case 1:
            tripletorbdd(pf)(l)=parms(counter); 
            break;
          case 2:
            singletorb(pf)(l)=parms(counter);
            break;
          case 3: 
            unpairedorb(pf)(l)=parms(counter);
            break;
	  }
          counter++;
        }
      }
    }
    //cout<<"before"<<endl;
    //for(int i=0;i<normalization(pf).GetSize();i++)
    // cout <<" normalization(pf)("<<i<<"):  "<<normalization(pf)(i)<<"  ";
    //cout<<endl;

    Getnormalization(tripletorbuu(pf), normalization(pf)(0));
    Getnormalization(tripletorbdd(pf), normalization(pf)(1));
    Getnormalization(singletorb(pf), normalization(pf)(2));
        
    //Unpairedorb renormalization
    for(int i=0; i<npairs(2); i++){
      unpairedorb_temp(i).Resize(ntote_pairs(pf));
      for(int j=0; j<ntote_pairs(pf); j++)  
        unpairedorb_temp(i)(j)=unpairedorb(pf)(i*ntote_pairs(pf)+j);
      Getnormalization(unpairedorb_temp(i), normalization(pf)(3+i));//we have to normalize for every row
      for(int j=0; j<ntote_pairs(pf); j++)
        unpairedorb(pf)(i*ntote_pairs(pf)+j)=unpairedorb_temp(i)(j);
    }
    //cout<<"after"<<endl;
    //for(int i=0;i<normalization(pf).GetSize();i++)
    // cout <<" normalization(pf)("<<i<<"):  "<<normalization(pf)(i)<<"  ";
    //cout<<endl;
  }
  if (optimize_pfwt){
    for(int pf=1; pf< pfwt.GetSize(); pf++) {
      pfwt(pf)=parms(counter);
      counter++;
    }
  }
  
    
  // if(mpi_info.node==0) {
  // cout<<"Set parameters values"<<endl;
  // for (int i=0;i<parms.GetSize();i++)
  //  cout << i+1<<"   "<<parms(i)<<endl;
  //}
  
}

void Pfaff_wf_data::getVarParms(Array1 <doublevar> & parms)
{

  if(optimize_pf || optimize_pfwt) {
    parms.Resize(nparms());
    getPFCoeff(parms);
  }
  else {
    parms.Resize(0);
  }
  
}



//----------------------------------------------------------------------
void Pfaff_wf_data::setVarParms(Array1 <doublevar> & parms)
{
  if(optimize_pf || optimize_pfwt) {
    setPFCoeff(parms);
  }
  else {
    parms.Resize(0);
  }
  
  int max=wfObserver.size();
  for(int i=0; i< max; i++) {
    wfObserver[i]->notify(all_wf_parms_change, 0);
  }
}

//----------------------------------------------------------------------

