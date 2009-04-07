/*
 
Copyright (C) 2008 Jindrich Kolorenc

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

#include "CBasis_function.h"
#include "converter.h"
#include <cstdlib>


void sort_lowest_first(vector <double> & vals, vector <int> & list){
  int n=vals.size();
  list.resize(n); 
  for (int i=0; i < n; i++) 
    list[i] = i;
  
  for (int i=1; i < n; i++) {
    const double temp = vals[i];
    //const double abstemp =  fabs(vals[i]);
    int j;
    for (j=i-1; j>=0 && vals[j]>temp; j--) {
      vals[j+1] = vals[j];
      list[j+1] = list[j];
    }
    vals[j+1] = temp;
    list[j+1] = i;
  }
  //for (int i=0; i < n; i++) 
  //cout << list[i]<<endl;
}

int Blochwave_function::read(string & filename, int & debug, int & norbs, int & node){
  if(node==0)
    cout <<"#################### Reading in  Bloch-WF ########################"<<endl;
  ifstream os(filename.c_str());
  //vector < vector < double > > k_vector_tmp;
  //vector <int>  bmax_tmp; // number of bands for all k points
  //vector < vector < vector < complex < double > > > > ckg_tmp;

  vector < vector < vector < complex < double > > > > ckg;
  vector < vector < double > > k_vector;
  vector <int>  bmax;
  int kmax;
  int bn_total; 

  string line;
  string space=" ";
  vector <string> words;
  //read ""Number of G-vectors
  while(getline(os, line)) {
    words.clear();
    split(line, space, words);
    if(words[0]=="Number" && words[1]=="of" && words[2]=="G-vectors"){
      getline(os,line);
      words.clear();
      split(line, space, words);
      nmax=atoi(words[0].c_str());
      if(node==0)
	cout <<"Number of G-vectors "<<nmax<<endl;
      vector <double> gtmp(3);
      getline(os,line);
      for(int i=0;i<nmax;i++){
	getline(os,line);
	//cout <<line<<endl;
	words.clear();
	split(line, space, words);
	gtmp[0]=atof(words[0].c_str());gtmp[1]=atof(words[1].c_str());gtmp[2]=atof(words[2].c_str());
	//cout <<gtmp[0]<<" "<<gtmp[1]<<" "<<gtmp[2]<<endl;
	g_vector.push_back(gtmp);
      }
    }
    if(words[0]=="Number" && words[1]=="of" && words[2]=="k-points"){
      getline(os,line);
      words.clear();
      split(line, space, words);
      kmax=atoi(words[0].c_str());
      bn_total=0;
      //cout <<"Number of  k-points "<<kmax<<endl;
      //read "k-point # ; # of bands (up spin/down spin) ; k-point coords (au)"
      //for each k-point read the bands
      getline(os, line);
      words.clear();
      split(line, space, words);
      for(int k=0;k<kmax;k++){
   	if(!(words[0]=="k-point" && words[1]=="#")){
	  cout <<"expected the k-point information"<<endl;
	  exit(1);
	}
	getline(os,line);
	words.clear();
	split(line, space, words);
	vector <double> vec(3);
	//example: 1    36    36 0.18699956271368 0.18699956271368 0.18699956271368
	vec[0]=atof(words[3].c_str());vec[1]=atof(words[4].c_str());vec[2]=atof(words[5].c_str());
	k_vector.push_back(vec);
	int b_max_up=atoi(words[1].c_str());
	int b_max_down=atoi(words[2].c_str());
	if(debug && node==0){
	  cout <<" k-point: "<<vec[0]<<" "<<vec[1]<<" "<<vec[2]<<endl;
	  cout <<" bands-up: "<<b_max_up<<" bands-down: "<<b_max_down<<endl;
	}
	bmax.push_back(b_max_up+b_max_down);
	vector < vector < complex < double > > > eigenveck;
	vector < double> eigenval_tmp;
	//read "Band, spin, eigenvalue (au)"
	getline(os, line);
	words.clear();
	split(line, space, words);
	for(int bup=0;bup<b_max_up;bup++){
	  if(!(words[0]=="Band," && words[1]=="spin,"&& words[2]=="eigenvalue")){
	    cout <<"expected Band, spin, eigenvalue (au)"<<endl;
	    exit(1);
	  }
	  getline(os, line);
	  words.clear();
	  split(line, space, words);
	  if(debug && node==0)
	    cout <<" Band: "<<words[0]<<" spin: "<<words[1]<<" eigenvalue (au): "<<words[2]<<endl;
	  eigenval_tmp.push_back(atof(words[2].c_str()));
	  getline(os, line);
	  if(line!="Eigenvector coefficients"){
	    cout <<"need Eigenvector coefficients"<<endl;
	    exit(1);
	  }
	  int kounter=0;
	  complex <double> eigenvec;
	  vector < complex <double> > eigenvecband;
	  while(getline(os, line)){
	    words.clear();
	    split(line, space, words);
	    if(words[0]=="Band,"|| words[0]=="k-point"){
	      break;
	    }
	    //cout <<line<<endl;
	    //cout <<line.substr(2,23)<<"  "<<line.substr(26,23)<<endl;
	    //eigenvec.real()=atof(line.substr(2,23).c_str());
	    //eigenvec.imag()=atof(line.substr(26,23).c_str());
            eigenvec=dcomplex(atof(line.substr(2,23).c_str()),
			      atof(line.substr(26,23).c_str()));
	    eigenvecband.push_back(eigenvec);
	    kounter++;
	  }
	  eigenveck.push_back(eigenvecband);
	  //cout <<"read coefs: "<<kounter<<" expected:  "<<nmax<<endl;
	  if(kounter!=nmax){
	    cout <<"did not find all the eigenvectors"<<endl;
	    exit(1);
	  }
	}//for each band-up;
	for(int bdown=0;bdown<b_max_down;bdown++){
	  if(!(words[0]=="Band," && words[1]=="spin,"&& words[2]=="eigenvalue")){
	    cout <<"expected Band, spin, eigenvalue (au)"<<endl;
	    exit(1);
	  }
	  getline(os, line);
	  words.clear();
	  split(line, space, words);
	  if(debug && node==0)
	    cout <<" Band: "<<words[0]<<" spin: "<<words[1]<<" eigenvalue (au): "<<words[2]<<endl;
	  eigenval_tmp.push_back(atof(words[2].c_str()));
	  getline(os, line);
	  if(line!="Eigenvector coefficients"){
	    cout <<"need Eigenvector coefficients"<<endl;
	    exit(1);
	  }
	  int kounter=0;
	  complex <double> eigenvec;
	  vector < complex <double> > eigenvecband;
	  while(getline(os, line)){
	    words.clear();
	    split(line, space, words);
	    if(words[0]=="Band," || words[0]=="k-point"){
	      break;
	    }
	    //cout <<line<<endl;
	    //eigenvec.real()=atof(line.substr(2,23).c_str());
	    //eigenvec.imag()=atof(line.substr(26,23).c_str());
            eigenvec=dcomplex(atof(line.substr(2,23).c_str()),
			      atof(line.substr(26,23).c_str()));
	    eigenvecband.push_back(eigenvec);
	    kounter++;
	  }
	  eigenveck.push_back(eigenvecband);
	  //cout <<"read coefs: "<<kounter<<" expected:  "<<nmax<<endl;
	  if(kounter!=nmax){
	    cout <<"did not find all the eigenvectors"<<endl;
	    exit(1);
	  }
	}//for each band-down;
	ckg.push_back(eigenveck);
	bn_total+=bmax[k];
	eigenvalues.push_back(eigenval_tmp);
      }//for each k point
    }
  }

  /*
  cout <<"gsize "<<g_vector.size()<<endl;
  //cout <<g_vector[0][0]<<" "<<g_vector[1][0]<<endl;

  cout <<"ksize "<<ckg.size()<<endl;
  for(int k=0;k<ckg.size();k++){
    cout <<"bsize "<<ckg[k].size()<<" ";
    for(int b=0;b<ckg[k].size();b++){
      cout <<ckg[k][b].size()<<"  ";
    }
    cout <<endl;
  }
  */
  
  if(debug && node==0){
    cout <<"Total number of k-points found : "<<ckg.size()<<endl;
    cout <<"Total number of bands found : "<<bn_total<<endl;
  }

  if(norbs>bn_total){
    cout <<"ERROR: number of requested orbitals: "<<norbs<<" is too large, found only "<<bn_total<<" orbitals in the file"<<endl;
    exit(1);
  }

  os.close();

  if(debug && node==0 )
    cout <<"#################### K-point analysis ############################"<<endl;

  vector <double> r(5); 
  //vector <int> ireducible_kvec;
  //vector <int> conjugate_kvec;
  //vector < vector <int> > nonzero_part_of_ireducible_kvec;
  vector < vector < complex <double> > > vals;
  //0.602221,1.40376,1.47251
  r[2]=0.602221;// 0.123/0.612098;
  r[3]=1.40376; //0.783/0.612098;	
  r[4]=1.47251; //0.3124/0.612098;

  if(debug && node==0)
    cout <<" r= ("<<r[2]<<","<<r[3]<<","<<r[4]<<")  "<<endl;
  //double glatdotr=0.0;
  //for(int d=0; d< 3; d++){
  //  glatdotr+=(0.612098)*r[d+2];
  //}
  //complex <double> phase=exp(I*glatdotr);
  //cout <<"phase "<<phase<<endl;
    for(int kn=0; kn< kmax; kn++) {
      //cout <<"kvec "<< k_vector(kn,0)<<" "<<k_vector(kn,1)<<" "<<k_vector(kn,2)<<endl;
      vector < complex <double> > valsb;
      for(int bn=0; bn< bmax[kn]; bn++){
	//cout <<" bn: "<<bn<<endl;
	complex <double> tmp=0.0;
	for(int i=0; i< nmax; i++){
	  double gdotr=0.0;
	  for(int d=0; d< 3; d++)
	    gdotr+=( g_vector[i][d]+k_vector[kn][d] )*r[d+2];

	  tmp+=ckg[kn][bn][i]*exp(I*gdotr);
	}//i
	valsb.push_back(tmp);
	//cout <<"val.real "<<tmp.real()<<" val.imag: "<<tmp.imag()<<"ratio "<<tmp.imag()/tmp.real()<<" phase: "<<atan(tmp.imag()/tmp.real())<<endl;
      }//band
      vals.push_back(valsb);
    }//kvec

    
 
    for(int kn1=0; kn1< kmax; kn1++) {
      for(int kn2=kn1+1; kn2< kmax; kn2++){
	for(int bn1=0; bn1< 1; bn1++)
	  for(int bn2=0; bn2< 1; bn2++)
	    if(abs(vals[kn1][bn1].real()-vals[kn2][bn2].real())<1e-5 &&
	       abs(vals[kn1][bn1].imag()+vals[kn2][bn2].imag())<1e-5){
	      //cout <<"kvec "<< k_vector[kn1][0]<<" "<<k_vector[kn1][1]<<" "<<k_vector[kn1][2]
	      //   <<"and kvec "<< k_vector[kn2][0]<<" "<<k_vector[kn2][1]<<" "<<k_vector[kn2][2]
	      //	   <<"are complex conjugate"<<endl;
	      ireducible_kvec.push_back(kn1);
	      conjugate_kvec.push_back(kn2);
	      vector <int> dummy(1);
	      dummy[0]=-1;
	      nonzero_part_of_ireducible_kvec.push_back(dummy);
	    }
      }//kn2
    }//kn1
    
  vector <int> not_paired_k;
  vector < vector <int> > which_part_not_paired;

  for(int kn1=0; kn1< kmax; kn1++){
    int found_one=0;
    for(int i=0;i<ireducible_kvec.size();i++){
      if(ireducible_kvec[i]==kn1 || conjugate_kvec[i]==kn1)
	found_one=1;
    }
    if(!found_one){
      if(debug && node==0)
	cout <<" k-point: "<<kn1+1<<" not paired"<<endl; 
      not_paired_k.push_back(kn1);
      vector <int> part;
      for(int bn=0; bn< bmax[kn1]; bn++){
	//doublevar phase=atan(vals(0,kn1,bn).imag()/vals(0,kn1,bn).real());
	if(debug && node==0)
	  cout <<"band "<<bn+1<<" value(r): "<<vals[kn1][bn]<<endl;
	if(fabs(vals[kn1][bn].imag())<1e-3){
	  part.push_back(0);
	}
	else if (fabs(vals[kn1][bn].real())<1e-3){
	  part.push_back(1);
	}
	else{
	  cout <<"the value of the band for the given K vector has no pure projection"<<endl;
	  exit(1);
	}
      }
      which_part_not_paired.push_back(part);
    }
  }
  //cout <<"Not paired K vectors: "<<endl;
  for(int i=0;i<not_paired_k.size();i++){
    //cout <<not_paired_k[i]<<endl;
    //for(int bn=0; bn< bmax[not_paired_k[i]]; bn++)
    //  cout <<which_part_not_paired[i][bn]<<"  ";
    //cout <<endl;
    ireducible_kvec.push_back(not_paired_k[i]);
    conjugate_kvec.push_back(-1);
    nonzero_part_of_ireducible_kvec.push_back(which_part_not_paired[i]);
  }
  if(debug && node==0){
    cout <<"Pair of  K, -K  points: "<<endl;
    cout <<"note: 0 means no pair "<<endl;
  }
  for(int i=0;i<ireducible_kvec.size();i++){
    if(debug && node==0)
      cout <<ireducible_kvec[i]+1<<" "<<conjugate_kvec[i]+1<<"  ";
    if(conjugate_kvec[i]==-1){
      if(debug && node==0)
	cout <<"("<<k_vector[ireducible_kvec[i]][0]<<", "<<k_vector[ireducible_kvec[i]][1]<<", "<<k_vector[ireducible_kvec[i]][2]<<") pure projections: "; 
      for(int j=0;j<nonzero_part_of_ireducible_kvec[i].size();j++){
	if(nonzero_part_of_ireducible_kvec[i][j]==0){
	  if(debug && node==0)
	    cout <<"real("<<j+1<<"), ";
	}
	else if (nonzero_part_of_ireducible_kvec[i][j]==1){
	  if(debug && node==0)
	      cout <<"imag("<<j+1<<"), ";
	}
	else{
	  cout <<" nonzero_part_of_ireducible_kvec[i][j] should contain only 0 or 1 "<<endl;
	  exit(1);
	}
      }//j
    }
    else{
      if(debug && node==0)
	cout <<"("<<k_vector[ireducible_kvec[i]][0]<<", "<<k_vector[ireducible_kvec[i]][1]<<", "<<k_vector[ireducible_kvec[i]][2]<<") with "<<
	  "("<<k_vector[conjugate_kvec[i]][0]<<", "<<k_vector[conjugate_kvec[i]][1]<<", "<<k_vector[conjugate_kvec[i]][2]<<")";
    }
    if(debug && node==0)
      cout <<endl;
  }
  
  if(debug && node==0)
    cout <<"#################### done K-point analysis ############################"<<endl;

  //reorder the bands by eigenvalues (lowest to highest)

  //vector <double> eigenvalueslinear;
  //vector <int> list;

  if(debug && node==0)
    cout <<"#################### creating ireducible eigenvectors for "<<norbs<<" orbitals ####"<<endl;
  for(int k=0;k<ckg.size();k++){
    for(int b=0;b<ckg[k].size();b++){
      eigenvalueslinear.push_back(eigenvalues[k][b]);
    }
  }
  sort_lowest_first(eigenvalueslinear,list);

  /*
  //vector <int> which_band, which_k ;
  for(int kk=0;kk<list.size();kk++){
    int ll=0;
    for(int k=0;k<ckg.size();k++){
      for(int b=0;b<ckg[k].size();b++){
	if(list[kk]==ll){
	  which_band.push_back(b);
	  which_k.push_back(k);
	}
	ll++;
      }
    }
  }
  
  */
  
  //for(int kk=0; kk<list.size();kk++){
  //  if(kk==norbs)
  //    cout <<"################# here is the cutoff ###########################"<<endl;
  //  cout <<kk+1<<" : "<<eigenvalueslinear[kk]<<endl;//" k: "<<which_k[kk]+1<<" b: "<<which_band[kk]+1<<endl;
  //}
  

  

  
  double eigenmax;
  if(norbs%2==0)
    eigenmax=eigenvalueslinear[norbs-1];
  else if(norbs%2==1 && norbs<eigenvalueslinear.size())
    eigenmax=eigenvalueslinear[norbs];
  else
    eigenmax= eigenmax=eigenvalueslinear[eigenvalueslinear.size()-1];

  if(debug && node==0)
    cout <<"eigenmax= "<<eigenmax<<endl;
  bn_total_final=0;

  int k_kount=0;
  for(int kn=0; kn< ireducible_kvec.size(); kn++) {
    k_vector_final.push_back(k_vector[ireducible_kvec[kn]]);
    vector < vector < complex < double > > > ckg_tmp;
    if(conjugate_kvec[kn]!=-1){ //if conjugate
      int b_kount=0;
      for(int bn=0; bn< bmax[ireducible_kvec[kn]]; bn++) {
	if( eigenvalues[ireducible_kvec[kn]][bn]< eigenmax+1e-5 ){
	  if(debug && node==0)
	    cout <<"kn: "<<ireducible_kvec[kn]<<" bn (both real and imag): "<<bn<<endl;
	  vector < complex < double > > ckg_tmp2;
	  for(int i=0; i< nmax; i++)
	    ckg_tmp2.push_back(ckg[ireducible_kvec[kn]][bn][i]);
	  ckg_tmp.push_back(ckg_tmp2);
	  orb_to_kn.push_back(k_kount);
	  orb_to_bn.push_back(b_kount++);
	}//if eigenvalues
      }//bn
    }//if conjugate
    else{ //not conjugate
      int bn1=0;
      int bn2;
      int b_kount=0;
      while(bn1< bmax[ireducible_kvec[kn]]){
	if(bn1+1==bmax[ireducible_kvec[kn]]) //if there is odd number of bands for the k-vector;
	  bn2=bn1;
	else
	  bn2=bn1+1;
	while(bn2< bmax[ireducible_kvec[kn]]){
	  //cout <<"bn1: "<<bn1<<" bn2: "<<bn2<<endl;
	  if( (eigenvalues[ireducible_kvec[kn]][bn1]< eigenmax+1e-5) && (eigenvalues[ireducible_kvec[kn]][bn2]< eigenmax+1e-5) ){
	    if(debug && node==0)
	      cout <<"bn1: "<<bn1<<" bn2: "<<bn2<<endl;
	    //cout <<"yes"<<endl;
	    vector < complex < double > > ckg_tmp2;
	    for(int i=0; i< nmax; i++){
	      complex < double > ckg_tmp3; 
	      double a1=ckg[ireducible_kvec[kn]][bn1][i].real();
	      double a2=ckg[ireducible_kvec[kn]][bn2][i].real();
	      double b1=ckg[ireducible_kvec[kn]][bn1][i].imag();
	      double b2=ckg[ireducible_kvec[kn]][bn2][i].imag();
	      if((nonzero_part_of_ireducible_kvec[kn][bn1]==0)
		 && (nonzero_part_of_ireducible_kvec[kn][bn2]==0) )
		{
		  //ckg_tmp3.real()=a1-b2; ckg_tmp3.imag()=b1+a2;
		  ckg_tmp3=dcomplex(a1-b2,b1+a2);
		}
	      else if ((nonzero_part_of_ireducible_kvec[kn][bn1]==0)
		       && (nonzero_part_of_ireducible_kvec[kn][bn2]==1) )
		{
		  //ckg_tmp3.real()=a1+a2; ckg_tmp3.imag()=b1+b2;
		  ckg_tmp3=dcomplex(a1+a2,b1+b2);
		}
	      else if ((nonzero_part_of_ireducible_kvec[kn][bn1]==1)
		       && (nonzero_part_of_ireducible_kvec[kn][bn2]==0) )
		{
		  //ckg_tmp3.real()=b1+b2; ckg_tmp3.imag()=a1+a2;
		  ckg_tmp3=dcomplex(b1+b2,a1+a2);
		}
	      else
		{
		  //ckg_tmp3.real()=a1+b2; ckg_tmp3.imag()=b1-a2;
		  ckg_tmp3=dcomplex(a1+b2,b1-a2);
		}
	      ckg_tmp2.push_back(ckg_tmp3);
	    }//i
	    ckg_tmp.push_back(ckg_tmp2);
	    orb_to_kn.push_back(k_kount);
	    orb_to_bn.push_back(b_kount++);
	    break;
	  }//if eigenvalues          
	  bn2++;
	}//end while
	bn1=bn2+1;
      }//end while
    }//not conjugate
    ckg_final.push_back(ckg_tmp);
    bmax_final.push_back(ckg_tmp.size());
    if(debug && node==0){
      cout <<" for  k-vector: "<<ireducible_kvec[kn]+1
	   <<" keeping lowest "<<ckg_tmp.size()<<" bands "<<endl;
    }
    bn_total_final+=ckg_tmp.size();
    k_kount++;
  }//kn
  kmax_final=ckg_final.size();
  ckg.clear();
  
  if(norbs>2*bn_total_final){
    cout <<"ERORR number of orbitals > 2xnbands assigned"<<endl;
    exit(1);
  }


  if(debug && node==0){
    for (int l=0;l<orb_to_kn.size();l++)
      cout <<"kn "<<orb_to_kn[l]+1<<" bn "<<orb_to_bn[l]+1<<endl;
  }

  if(node==0){
    cout <<"Total number of k-points kept : "<<kmax_final<<endl;
    cout <<"Total number of bands kept : "<<bn_total_final<<endl;
  }
    
  if(debug && node==0)
    cout <<"#################### done creating ireducible eigenvectors #####"<<endl;

  if(node==0)
    cout <<"#################### done reading in blochwave wf ################"<<endl;
  return 0;  
}


int Blochwave_function::nfunc()
{
  // cout << "Blochwave::nfunc() " << nmax << endl;
  return bn_total_final;
}

int Blochwave_function::showinfo(string & indent, ostream & os)
{
  os << indent << "Blochwave function\n";
  os << indent << kmax_final << " k points" << endl;
  os << indent << bn_total_final << " Blochwaves"  << endl;
  os << indent << nmax << " plane waves for each Bloch wave" << endl;
  return 1;
}


void Blochwave_function::getKvectors(vector < vector <double > > & primlatvec, vector < vector <double > > & VEC){
  for(int kn=0; kn< kmax_final;kn++){
      vector <double> kvector_to_print(3);
      for(int d1=0;d1<3;d1++){
	kvector_to_print[d1]=0;
	for(int d2=0;d2<3;d2++){
	  kvector_to_print[d1]+=primlatvec[d1][d2]*k_vector_final[kn][d2];
	}
	kvector_to_print[d1]/=pi;
      }
      VEC.push_back(kvector_to_print);
    }//kn
}

int Blochwave_function::writeinput(string & indent, ostream & os, int & norbs, vector < vector <double > > & primlatvec, string & outputname)
{
  cout <<"#################### Writing out  Bloch-WF #######################"<<endl;
  string indent2=indent+indent;
  int width=24;
  int prec=14;
  int cycles=int(norbs/2)+norbs%2;
  //os.setf(ios::scientific);
  int orbcounter=0;
  for(int kn=0; kn< kmax_final;kn++){
    if(kn==orb_to_kn[orbcounter]){
      os << indent << "KVECTOR {\n";
      vector <double> kvector_to_print(3);
      for(int d1=0;d1<3;d1++){
	kvector_to_print[d1]=0;
	for(int d2=0;d2<3;d2++){
	  kvector_to_print[d1]+=primlatvec[d1][d2]*k_vector_final[kn][d2];
	}
	kvector_to_print[d1]/=pi;
      }
      os << indent <<"  "<<kvector_to_print[0]<<"   "<<kvector_to_print[1]<<"  "<<kvector_to_print[2] << endl;

      for(int bn=0; bn< bmax_final[kn]; bn++) {
	if(bn==orb_to_bn[orbcounter]){
	  os << indent2 << "BAND {  \n";
	  string basename,basename2;
	  char strbuff[40];
	  sprintf(strbuff, "%d", orbcounter+1);
	  basename2 = outputname;
	  basename2 += ".orb";
	  basename2 += strbuff;
	  for(int part=0;part<2;part++){
	    if(part==0)
	      basename=basename2+".real";
	    else
	      basename=basename2+".imag";
	    
	    string cubename=basename+".cube";
	    os << indent2 <<cubename<<endl;
	  }
	  os << indent2 << "}" << endl;
	  orbcounter++;
	}
      }//bn
      os << indent << "}" << endl;
    }
  }//kn
  



  /* original bloch wave basis output
  os << "CBASIS {\n";
  os << indent << "origin\n";
  os << indent << "BLOCHWAVE\n";
  os << indent << "GVECTOR { \n";
  int width=24;
  int prec=14;
  os.setf(ios::scientific);
  for(int i=0; i< nmax; i++) {
    os << indent << "  " <<setw(width)<<setprecision(prec)<< g_vector[i][0] 
       << "  " <<setw(width)<<setprecision(prec)<< g_vector[i][1]
       << "  " <<setw(width)<<setprecision(prec)<< g_vector[i][2] << endl;
  }
  os << indent << "}" << endl;

  for(int kn=0; kn< kmax_final;kn++){
     os << indent << "KVECTOR {\n";
     os << indent <<"  "<<setw(width)<<setprecision(prec)<<k_vector_final[kn][0]
	                <<setw(width)<<setprecision(prec)<<k_vector_final[kn][1]
                        <<setw(width)<<setprecision(prec)<<k_vector_final[kn][2] << endl;
     for(int bn=0; bn< bmax_final[kn]; bn++) {
        os << indent2 << "BAND {  \n";
	for(int i=0; i< nmax; i++){
	  os << indent2 <<setw(width)<<setprecision(prec)<<ckg_final[ireducible_kvec[kn]][bn][i].real()
	                <<setw(width)<<setprecision(prec)<<ckg_final[ireducible_kvec[kn]][bn][i].imag() << endl;
	}
	os << indent2 << "}" << endl;
     }
     os << indent << "}" << endl;
  }
  os<<"}\n";
  os.unsetf(ios::scientific);
  os<<setprecision(6);

  */
  cout <<"#################### done writing out Bloch-WF ###################"<<endl;
  return 1;
}



void Blochwave_function::calcVal(const vector <double> & r,
                                 complex <double> & symvals,
                                 int & orb)
{
  //cout << "calcVal " << endl;
  // assert(r.GetDim(0) >= 5);
  //assert(symvals.GetDim(0) >= bn_total+startfill);
  //int index=0;
  

  //for(int kn=0; kn< kmax_final; kn++)
  //for(int bn=0; bn< bmax_final[kn]; bn++) {
  int kn=orb_to_kn[orb];
  int bn=orb_to_bn[orb];
  //cout <<"kn "<<kn<<" bn "<<bn<<endl;
  symvals=0.0;
  for(int i=0; i< nmax; i++) {
    double gdotr=0.0;
    for(int d=0; d< 3; d++)
      gdotr+=( g_vector[i][d]+k_vector_final[kn][d] )*r[d];
    symvals+=ckg_final[kn][bn][i]*exp(I*gdotr);
  }
      //index++;
      //}//each bn 
  //cout << "done" << endl;
}


//------------------------------------------------------------------------
