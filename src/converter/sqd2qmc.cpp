#include "converter.h"
#include "wf_writer.h"
#include "basis_writer.h"
#include "Pseudo_writer.h"
#include "elements.h"  
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include "vecmath.h"
#include "elements.h"
#include "sqd2qmc.h"
#include <libxml/parser.h>
#include <libxml/xmlmemory.h>
#include "hdf5.h"




using namespace std;

void usage(const char * name) {
  cout << "usage: " << name <<   " <options> <.out file> " << endl;
  cout << "Where options can be: \n";
  cout << "-o <string> Base name for your run case\n";
  cout << "-debug to have more informative printout\n";
  exit(1);
}


string givesymmetry(int & i){
  if(i==0)
    return "S";
  else if(i==1)
    return "P";
  else if(i==2)
    return "5D";
  else if(i==3)
    return "7F";
  else if(i==4)
    return "9G";
  else if(i==5)
    return "11H";
  else
    return "what was that?";
}


void read_parameter(xmlDocPtr doc, xmlNodePtr cur, group  & g){
  cur = cur->xmlChildrenNode;
  while (cur != NULL) {
    if ((!xmlStrcmp(cur->name, (const xmlChar *)"parameter"))){
      if((!xmlStrcmp(xmlGetProp(cur, (const xmlChar *)"name"),(const xmlChar *)"charge"))){
	g.charge= atoi((const char*) xmlNodeListGetString(doc, cur->xmlChildrenNode, 1));
	//cout <<" charge "<<g.charge<<endl;
      }
      
    }
    cur = cur->next; 
  }
  
  return;
}




void read_group(xmlDocPtr doc, xmlNodePtr cur, particleset  & p){
  cur = cur->xmlChildrenNode;
  group tmpg;

  while (cur != NULL) {
    if ((!xmlStrcmp(cur->name, (const xmlChar *)"group"))) {
      //cout <<"found group "<<endl;
      tmpg.name=((char *) xmlGetProp(cur, (const xmlChar *)"name"));\

      if(xmlGetProp(cur, (const xmlChar *)"size")!= NULL){
	tmpg.size=atoi((const char*)xmlGetProp(cur, (const xmlChar *)"size"));
      }
      read_parameter(doc,cur,tmpg);
      p.g.push_back(tmpg);
    }
    cur = cur->next;
  }
  
  return;
}





void read_particleset(xmlDocPtr doc, xmlNodePtr cur, xml_system_data  & qmc_xml){
  particleset ptmp;
    if ((!xmlStrcmp(cur->name, (const xmlChar *)"particleset"))){
      ptmp.name=((char *) xmlGetProp(cur, (const xmlChar *)"name"));
      if(xmlGetProp(cur, (const xmlChar *)"size")!=NULL)
	ptmp.size=atoi((const char*)xmlGetProp(cur, (const xmlChar *)"size"));
      read_group(doc,cur,ptmp);
      qmc_xml.p.push_back(ptmp);
    }
}


void read_basisGroup(xmlDocPtr doc, xmlNodePtr cur, vector<basisGroup> & bg){
  cur = cur->xmlChildrenNode;
  while (cur != NULL) {
    if ((!xmlStrcmp(cur->name, (const xmlChar *)"basisGroup"))){
      basisGroup orbtmp;
      //cout <<"found basisGroup "<<" for orbital "<<xmlGetProp(cur, (const xmlChar *)"rid") <<endl;
      orbtmp.rid=(char *)xmlGetProp(cur, (const xmlChar *)"rid");
      orbtmp.ds=(char *)xmlGetProp(cur, (const xmlChar *)"ds");
      orbtmp.n=atoi((const char*)xmlGetProp(cur, (const xmlChar *)"n"));
      orbtmp.m=atoi((const char*)xmlGetProp(cur, (const xmlChar *)"m"));
      orbtmp.l=atoi((const char*)xmlGetProp(cur, (const xmlChar *)"l"));
      orbtmp.s=atoi((const char*)xmlGetProp(cur, (const xmlChar *)"s"));
      bg.push_back(orbtmp);
     }
    cur = cur->next;
  } 
}




void read_atomicBasisSet(xmlDocPtr doc, xmlNodePtr cur, atomicBasisSet & abs){
   cur = cur->xmlChildrenNode;
   while (cur != NULL) {
    if ((!xmlStrcmp(cur->name, (const xmlChar *)"atomicBasisSet"))){
      //cout <<"found atomicBasisSet"<<endl;
      abs.filename=(char *)xmlGetProp(cur, (const xmlChar *)"href");
      abs.elementType=(char *)xmlGetProp(cur, (const xmlChar *)"elementType");
      read_basisGroup(doc,cur,abs.orbs);
    }
    cur = cur->next;
  }
}



void read_basisset(xmlDocPtr doc, xmlNodePtr cur, atomicBasisSet & abs){
  cur = cur->xmlChildrenNode;
  while (cur != NULL) {
    if ((!xmlStrcmp(cur->name, (const xmlChar *)"basisset"))){
      //cout <<"found basisset"<<endl;
      read_atomicBasisSet(doc,cur,abs);
    }
    cur = cur->next;
  }
}


void read_determinantset(xmlDocPtr doc, xmlNodePtr cur, wavefunction  & wf){
  cur = cur->xmlChildrenNode;
  while (cur != NULL) {
    // cout <<"cur->name "<<cur->name<<endl;
     if ((!xmlStrcmp(cur->name, (const xmlChar *)"determinantset"))){
       //cout <<"found determinantset of type "<< xmlGetProp(cur, (const xmlChar *)"type")
       //    <<" for source "<<xmlGetProp(cur, (const xmlChar *)"source")<<endl;
       read_basisset(doc,cur,wf.abs);
     }
     cur = cur->next;
   }
}


void read_wavefunction(xmlDocPtr doc, xmlNodePtr cur, wavefunction  & wf){
  if ((!xmlStrcmp(cur->name, (const xmlChar *)"wavefunction"))){
    //cout <<"found wavefunction of name "<< xmlGetProp(cur, (const xmlChar *)"name")<<endl;
    read_determinantset(doc,cur,wf);
  }
}


static void parse_from_xml(xmlDocPtr doc, xml_system_data  & qmc_xml){
  xmlNode *cur= NULL;
  cur = xmlDocGetRootElement(doc);
  if (xmlStrcmp(cur->name, (const xmlChar *) "qmcsystem")){
    cout <<"This xml file does not contain qmcsystem!"<<endl; exit(-1);
  }
  cur = cur->xmlChildrenNode;
  while (cur != NULL) {
    read_particleset(doc,cur,qmc_xml);
    read_wavefunction(doc,cur,qmc_xml.wf);
    cur = cur->next;
  } 

}

void read_into_array(hid_t grp, const char* name, vector <double> & ref){
  // Turn off error printing
  //H5E_auto_t func;
  //void *client_data;
  //H5Eget_auto (&func, &client_data);
  //H5Eset_auto (NULL, NULL);
  
  hid_t dataset = H5Dopen(grp,name,H5P_DEFAULT);
  if (dataset > -1) {
    hid_t datatype=H5Dget_type(dataset);
    hsize_t dim_out,dims_out;

    hid_t dataspace = H5Dget_space(dataset);
    hid_t status = H5Sget_simple_extent_dims(dataspace, &dim_out, NULL);
    H5Sclose(dataspace);
    dims_out=H5Tget_size(datatype);
    //  cout <<" dim_out "<<dim_out<<" dims_out"<<dims_out<<endl;
    ref.resize(dim_out);
    hid_t ret = H5Dread(dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT,&(ref[0]));
    H5Tclose(datatype);
    H5Dclose(dataset);
  }
  // Turn error printing back on
  //H5Eset_auto (func, client_data);
}

void read_into_array(hid_t grp, const char* name, vector <int> & ref){
  // Turn off error printing
  //H5E_auto_t func;
  //void *client_data;
  //H5Eget_auto (&func, &client_data);
  //H5Eset_auto (NULL, NULL);
  
  hid_t dataset = H5Dopen(grp,name,H5P_DEFAULT);
  if (dataset > -1) {
    hid_t datatype=H5Dget_type(dataset);
    hsize_t dim_out,dims_out;

    hid_t dataspace = H5Dget_space(dataset);
    hid_t status = H5Sget_simple_extent_dims(dataspace, &dim_out, NULL);
    H5Sclose(dataspace);
    dims_out=H5Tget_size(datatype);
    //  cout <<" dim_out "<<dim_out<<" dims_out"<<dims_out<<endl;
    ref.resize(dim_out);
    hid_t ret = H5Dread(dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT,&(ref[0]));
    H5Tclose(datatype);
    H5Dclose(dataset);
  }
  // Turn error printing back on
  //H5Eset_auto (func, client_data);
}


void read_into_array(hid_t grp, const char* name, double & value){
  // Turn off error printing
  //H5E_auto_t func;
  //void *client_data;
  //H5Eget_auto (&func, &client_data);
  //H5Eset_auto (NULL, NULL);
  
  hid_t dataset = H5Dopen(grp,name,H5P_DEFAULT);
  if (dataset > -1) {
    hid_t datatype=H5Dget_type(dataset);
    hsize_t dim_out,dims_out;

    hid_t dataspace = H5Dget_space(dataset);
    hid_t status = H5Sget_simple_extent_dims(dataspace, &dim_out, NULL);
    H5Sclose(dataspace);
    
    //dims_out=H5Tget_size(datatype);
    //  cout <<" dim_out "<<dim_out<<" dims_out"<<dims_out<<endl;
    vector <double> ref;
    ref.resize(dim_out);
    hid_t ret = H5Dread(dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT,&(ref[0]));
    value=ref[0];
    H5Tclose(datatype);
    H5Dclose(dataset);
  }
  // Turn error printing back on
  //H5Eset_auto (func, client_data);
}


void read_into_array(hid_t grp, const char* name, int & value){
  // Turn off error printing
  //H5E_auto_t func;
  //void *client_data;
  //H5Eget_auto (&func, &client_data);
  //H5Eset_auto (NULL, NULL);
  
  hid_t dataset = H5Dopen(grp,name,H5P_DEFAULT);
  if (dataset > -1) {
    hid_t datatype=H5Dget_type(dataset);
    hsize_t dim_out,dims_out;

    hid_t dataspace = H5Dget_space(dataset);
    hid_t status = H5Sget_simple_extent_dims(dataspace, &dim_out, NULL);
    H5Sclose(dataspace);
    
    //dims_out=H5Tget_size(datatype);
    //  cout <<" dim_out "<<dim_out<<" dims_out"<<dims_out<<endl;
    vector <int> ref;
    ref.resize(dim_out);
    hid_t ret = H5Dread(dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT,&(ref[0]));
    value=ref[0];
    H5Tclose(datatype);
    H5Dclose(dataset);
  }
  // Turn error printing back on
  //H5Eset_auto (func, client_data);
}


void splinefit(vector <double>& x, vector <double> & y, double yp1, double ypn, 
	       vector < vector <double> > & coef, vector <double> & pos)
{

  //following stolen from Numerical Recipes, more or less
  int n=x.size();
  pos.resize(n);
  pos=x;

  vector <double> y2(n), u(n);
  double sig, p, qn, un, hi;

  if(yp1 > .99e30) {
    y2[0]=0;
    u[0]=0;
  }
  else {
    y2[0]=-.5;
    u[0]=(3./(x[1]-x[0]))*((y[1]-y[0])/(x[1]-x[0])-yp1);
  }

  for(int i=1; i<n-1; i++)
  {
    sig=(x[i]-x[i-1])/(x[i+1]-x[i-1]);
    p=sig*y2[i-1]+2.;
    y2[i]=(sig-1.)/p;
    u[i]=(6.*((y[i+1]-y[i])/(x[i+1]-x[i])
              -(y[i]-y[i-1])/(x[i]-x[i-1]))
          /(x[i+1]-x[i-1])-sig*u[i-1])/p;
  }

  if(ypn>.99e30)
  {
    qn=0;
    un=0;
  }
  else
  {
    qn=.5;
    un=(3./(x[n-1]-x[n-2]))*(ypn-(y[n-1]-y[n-2])/(x[n-1]-x[n-2]));
  }

  y2[n-1]=(un-qn*u[n-2])/(qn*y2[n-2]+1.);

  for(int k=n-2; k>=0; k--)
  {
    y2[k]=y2[k]*y2[k+1]+u[k];
  }

  for(int i=0; i<n-1; i++)
  {
    vector <double> coeff_tmp(4);
    coeff_tmp[0]=y[i];
    hi=x[i+1]-x[i];
    coeff_tmp[1]=(y[i+1]-y[i])/hi - hi*(2.*y2[i]+y2[i+1])/6.;
    coeff_tmp[2]=y2[i]/2.;
    coeff_tmp[3]=(y2[i+1]-y2[i])/(6.*hi);
    coef.push_back(coeff_tmp);
  }
  vector <double> coeff_tmp(4);
  coeff_tmp[0]=y[n-1];
  coeff_tmp[1]=0;
  coeff_tmp[2]=0;
  coeff_tmp[3]=0;
  coef.push_back(coeff_tmp);

}


int getInterval(double r, vector <double>& pos ) {
  int npts=pos.size();
  int upper=npts;
  int lower=0;
  int guess=int(npts/2);

  //cout << "r= " << r << endl;

  while(true) {
    if(r >= pos[guess]) {
      //cout << "greater " << pos(guess) << endl;
      if(r < pos[guess+1] ) {
        //cout << "found it " << pos(guess+1) << endl;
        return guess;
      }
      else { 
        lower=guess;
        guess=(lower+upper)/2;
        //cout << "didn't find it: new lower " << lower 
        //     << " new guess " << guess << endl;
      }
    }
    if(r < pos[guess] ) {
      upper=guess;
      guess=(lower+upper)/2;
      //cout << "less than " << pos(guess) << endl;
      //cout << "new upper " << upper 
      //     << " new guess " << guess << endl;
    }  
  } 
  
}



double getVal(double r, vector <double>& pos, vector < vector <double> > & coeff) {
  int i=getInterval(r,pos);
  double height=r-pos[i];
  return coeff[i][0]+height*(coeff[i][1]
			     +height*(coeff[i][2]
				      +height*coeff[i][3]));
}



double find_cuttoff(vector <double>& pos, vector < vector <double> > & coeff){
  const double cutoff_threshold=1e-12;
  for(int i=pos.size()-1; i >=0; i--) {

    if(fabs(coeff[i][0]) > cutoff_threshold) {
      //cout << "cutoff " << x+step << endl;
      if(i<pos.size()-1){
	coeff.resize(i+1);
	pos.resize(i+1);
	return pos[i+1];
      }
      else
	return pos[i];
    }
  }
  return 0;
}


void extrapolate_to_linear_grid(vector <double> & x, vector <double> & y , double & cusp, double & dismin2){
  //copy the arrays
  
  vector <double> xtmp;
  vector <double> ytmp;
  
  double xxi=-1;
  double dismin=0.001;
  if(dismin2<dismin){
    cout <<"extrapolate_to_linear_grid: using to fine grid "<<endl;
  }

 
  for(int i=0;i<y.size();i++){
    if(((x[i]-xxi) > dismin) && x[i] > 1.0e-4){
      xtmp.push_back(x[i]);
      ytmp.push_back(y[i]);
      xxi=x[i];
    }
  }
  
  int ne=2;
  double yp1=(ytmp[ne]-ytmp[0])/(xtmp[ne]-xtmp[1]);
  double y0=ytmp[1]+(0.-xtmp[0])*yp1;
  //cout <<" "<<ytmp[0]<<" "<<ytmp[ne]<<" "<<xtmp[0]<<" "<<xtmp[ne]<<endl;;
  cout <<" derivative at the origin "<<yp1<<" yp1/y "<<yp1/ytmp[1]<<endl;

  double ypn=(ytmp[xtmp.size()-1] - ytmp[xtmp.size()-3])/(xtmp[xtmp.size()-1]-xtmp[xtmp.size()-3]);
  //cout <<" derivative at the end " <<ypn<< " ypn/yn "<<ypn/ytmp[xtmp.size()-2]<<endl;
  
  vector < vector <double> > coef;
  vector <double> pos;

  //fix the derivatives and values at the origin
  yp1=cusp*ytmp[0];
  //ypn=0.0
  ytmp[0]=y0;
  xtmp[0]=0.0;

  
  splinefit(xtmp, ytmp, yp1, ypn, coef, pos);
  double cutoff=find_cuttoff(pos,coef);
  if(cutoff<x[x.size()-1])
    cout <<" found cutoff "<<cutoff<<endl;
  
  //double dismin2=0.01;
  double xmin=0.0;
  double xmax=pos[pos.size()-1];
  x.resize(0);
  y.resize(0);
  double r=xmin;
  while(r<=xmax){
    x.push_back(r);
    y.push_back(getVal(r,pos,coef));
    r+=dismin2;
  }
  
  return;
}




void writebasis(ostream & os, xml_system_data  & qmc_xml, string output, string offset){
  /*
  cout <<"found the following atomicBasisSet: "<<endl;
  for (int i=0; i<qmc_xml.wf.abs.orbs.size();i++){
    cout <<qmc_xml.wf.abs.orbs[i].rid<<" "<<qmc_xml.wf.abs.orbs[i].ds<<" "
	 <<qmc_xml.wf.abs.orbs[i].n<<" "
	 <<qmc_xml.wf.abs.orbs[i].l<<" "
	 <<qmc_xml.wf.abs.orbs[i].m<<endl;
  }
  */

  hid_t   file1;
  herr_t ret;
  

  string filename=qmc_xml.wf.abs.filename;
  cout <<" Using orbitals from file "<<filename<<endl;
  file1 = H5Fopen (qmc_xml.wf.abs.filename, H5F_ACC_RDWR, H5P_DEFAULT);

  double rf,ri, log_mesh_spacing;
  int npts;

  char grid_type[256];
  hid_t dataset = H5Dopen(file1, "radial_basis_states/grid/type",H5P_DEFAULT);
  hid_t datatype=H5Dget_type(dataset);
  ret = H5Dread(dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, grid_type);
  H5Tclose(datatype);
  H5Dclose(dataset);

  if(strcmp(grid_type,"log")!=0){
    cout <<" only log mesh is supported for now, exit"; exit(-1);
  }

  read_into_array(file1,"radial_basis_states/grid/ri",ri);
  read_into_array(file1,"radial_basis_states/grid/rf",rf);
  read_into_array(file1,"radial_basis_states/grid/npts",npts);

  //cout <<" ri "<<ri<<" rf "<<rf<<" npts "<<npts<<endl;
  
  log_mesh_spacing=exp(log(rf/ri)/(npts-1));
  /*
  cout.precision(20);                 
  cout.width(25); 
  cout <<" calculated log_mesh_spacing "<<log_mesh_spacing<<endl;
  */
  vector <double> r_log;
  double x=ri;
  for(int j=0;j<npts;j++){
    r_log.push_back(x);
    x*=log_mesh_spacing;
  }
  

  vector <radial_orbital> rad_orb_log_mesh;
  

  char * unique_rid="";
  vector <string> unique_rid_orb_name;

  for (int i=0; i<qmc_xml.wf.abs.orbs.size();i++){  
    if(strcmp(unique_rid,qmc_xml.wf.abs.orbs[i].ds)) { 
      radial_orbital orbtmp;
      os << offset<<offset<<"SPLINE { "<<endl;
      os << offset<<offset<<offset<<givesymmetry(qmc_xml.wf.abs.orbs[i].l)<<endl;
      os << offset<<offset<<offset<<"INCLUDE "<<output<<"."<<qmc_xml.wf.abs.orbs[i].ds<<endl;
      os << offset<<offset<<"       } "<<endl;
      unique_rid=qmc_xml.wf.abs.orbs[i].ds;
      orbtmp.orb_name=unique_rid;
      string tmps=unique_rid;
      string nameofthedataset;

      cout <<" reading data for spline "<<unique_rid<< " of "<< givesymmetry(qmc_xml.wf.abs.orbs[i].l)<<" symmetry"<<endl;
      
      nameofthedataset="radial_basis_states/"+tmps+"/eigenvalue";
      //cout <<" reading in to array dataset "<<nameofthedataset<<endl;    
      read_into_array(file1,nameofthedataset.c_str(),orbtmp.eigenvalue);

      nameofthedataset="radial_basis_states/"+tmps+"/power";
      //cout <<" reading in to array dataset "<<nameofthedataset<<endl;    
      read_into_array(file1,nameofthedataset.c_str(),orbtmp.power);
      
      nameofthedataset="radial_basis_states/"+tmps+"/quantum_numbers";
      //cout <<" reading in to array dataset "<<nameofthedataset<<endl;    
      read_into_array(file1,nameofthedataset.c_str(),orbtmp.quantum_numbers);

      nameofthedataset="radial_basis_states/"+tmps+"/radial_orbital";
      //cout <<" reading in to array dataset "<<nameofthedataset<<endl;    
      read_into_array(file1,nameofthedataset.c_str(),orbtmp.radial_part);

      nameofthedataset="radial_basis_states/"+tmps+"/uofr";
      //cout <<" reading in to array dataset "<<nameofthedataset<<endl;    
      read_into_array(file1,nameofthedataset.c_str(),orbtmp.uofr);

      orbtmp.r_log=r_log;

      rad_orb_log_mesh.push_back(orbtmp);
    }
  }
  //cout <<"closing file"<<endl;
  ret = H5Fclose (file1);

  


  double cusp=0.0;
 
  

  //print the orbitals to files
  for(int i=0;i<rad_orb_log_mesh.size();i++){
    //ofstream orbout(rad_orb_log_mesh[i].orb_name.c_str());
     FILE * orbout;
     string filename=output+"."+rad_orb_log_mesh[i].orb_name;
     orbout = fopen (filename.c_str(),"w");
     fprintf (orbout,"# eigenvalue %g \n",rad_orb_log_mesh[i].eigenvalue);
     fprintf (orbout,"# power %d \n", rad_orb_log_mesh[i].power);
     fprintf (orbout,"# n %d l %d \n",rad_orb_log_mesh[i].quantum_numbers[0], rad_orb_log_mesh[i].quantum_numbers[1]);

     /*
     orbout<<"# eigenvalue "<<rad_orb_log_mesh[i].eigenvalue<<endl;
     orbout<<"# power "<<rad_orb_log_mesh[i].power<<endl;
     orbout<<"# n, l "<<rad_orb_log_mesh[i].quantum_numbers[0]<<" "<<rad_orb_log_mesh[i].quantum_numbers[1]<<endl;
     */

     

     double spacing=0.005;
     extrapolate_to_linear_grid(rad_orb_log_mesh[i].r_log, rad_orb_log_mesh[i].radial_part,cusp, spacing);

     //orbout.precision(16);                 
     //orbout.width(20); 
     for(int j=0;j<rad_orb_log_mesh[i].radial_part.size();j++){
       fprintf (orbout, "%15.12e %15.12e \n",rad_orb_log_mesh[i].r_log[j],rad_orb_log_mesh[i].radial_part[j]);
       //orbout << rad_orb_log_mesh[i].r_log[j] <<" "<<rad_orb_log_mesh[i].radial_part[j]<<endl;
     }
     fclose (orbout);
     //orbout.close();
  }
}


void getoccupations(xml_system_data  & qmc_xml, std::string & calculatetype, vector <int> & occ_up, vector <int> & occ_down){
  char * unique_rid="";
  int n,l,m,s, nold, lold, mold, sold;
  int nup=0;
  int ndown=0;
  vector < vector <int> > unique;
  vector <int>  tmp;
  
  for (int i=0; i<qmc_xml.wf.abs.orbs.size();i++){
    s=qmc_xml.wf.abs.orbs[i].s;
    l=qmc_xml.wf.abs.orbs[i].l;
    if(s>0)
      nup++;
    else
      ndown++;
    
    if(strcmp(unique_rid,qmc_xml.wf.abs.orbs[i].ds)) { 
      unique_rid=qmc_xml.wf.abs.orbs[i].ds;
      tmp.push_back(i);
      for (int j=i+1; j<qmc_xml.wf.abs.orbs.size();j++){
	if(!strcmp(unique_rid,qmc_xml.wf.abs.orbs[j].ds)){
	  tmp.push_back(j);
	}
      }
      unique.push_back(tmp);
      tmp.resize(0);
    }
  }

  // for(int i=0;i<unique.size();i++){
  //  for(int j=0;j<unique[i].size();j++){
  //   int k=unique[i][j];
  //   cout << qmc_xml.wf.abs.orbs[k].ds <<" "<<qmc_xml.wf.abs.orbs[k].n<<" "<<qmc_xml.wf.abs.orbs[k].l<<" "<<qmc_xml.wf.abs.orbs[k].m<<" "<<qmc_xml.wf.abs.orbs[k].s<<endl;
  // }
  // cout <<endl;
  //}

  //vector <int> occ_up;
  //vector <int> occ_down;
  int all_same=1;

  int occ=0;
  for(int i=0;i<unique.size();i++){
    int ss=qmc_xml.wf.abs.orbs[unique[i][0]].s;
    int j=1;
    while (j<unique[i].size() && all_same){
      s=qmc_xml.wf.abs.orbs[unique[i][j]].s;
      if(s!=ss)
	all_same=0;
      j++;
    }
  }
  //string calculatetype;
  if (all_same){ //UHF
    calculatetype="UHF";
    int occ=1;
    for(int i=0;i<unique.size();i++){
      int ll=qmc_xml.wf.abs.orbs[unique[i][0]].l;
      for(int j=0;j<unique[i].size();j++){
	int s=qmc_xml.wf.abs.orbs[unique[i][j]].s;
	if(s>0)
	  occ_up.push_back(occ++);
	else
	  occ_down.push_back(occ++);
      }
      occ+=(2*ll+1-unique[i].size());
    }
  }
  else { //RHF or ROHF
    int occup=1;
    int occdown=1;
    for(int i=0;i<unique.size();i++){
      int ll=qmc_xml.wf.abs.orbs[unique[i][0]].l;
      int tmpup=0;
      int tmpdown=0;
      for(int j=0;j<unique[i].size();j++){
	int s=qmc_xml.wf.abs.orbs[unique[i][j]].s;
	if(s>0){
	  occ_up.push_back(occup++);
	  tmpup++;
	}
	else{
	  occ_down.push_back(occdown++);
	  tmpdown++;
	}
      }
      occup+=(2*ll+1-tmpup);
      occdown+=(2*ll+1-tmpdown);
    }
    if(occ_up.size()==occ_down.size())
      calculatetype="RHF";
    else
      calculatetype="ROHF";
  }      
    

  /*
  for (int i=0;i<occ_up.size();i++)
    cout <<occ_up[i]+1<<" ";
  cout <<endl;

  for (int i=0;i<occ_down.size();i++)
    cout <<occ_down[i]+1<<" ";
  cout <<endl;
  */

}



void read_grid(xmlNodePtr cur, vector <double> & pos){
  int npts= atoi((const char*) xmlGetProp(cur, (const xmlChar *)"npts"));
  pos.resize(npts);
  //cout <<"npts "<< npts<<endl;
  xmlNodePtr cur1=cur->children;
  while(cur1 != NULL)
    {
      string cname((const char*)cur1->name);
      pos.resize(npts);
      if(cname == "data"){
	putContent(pos,cur1);
      }
      cur1 = cur1->next;
    }
}

void read_semilocal(xmlNodePtr cur, vector <double> & pos, double & cutoff){
  if(xmlStrcmp(xmlGetProp(cur, (const xmlChar *)"format"),(const xmlChar *)"r*V")){
    cout <<"expected format in the semilocal should be r*V, exit"<<endl;exit(-1);
  }
  xmlNodePtr cur1=cur->children;
  while(cur1 != NULL)
    {
      if(!xmlStrcmp(cur1->name, (const xmlChar *) "vps")){
	cutoff=atof((const char*) xmlGetProp(cur1, (const xmlChar *)"cutoff"));
	string symmetry=((const char*) xmlGetProp(cur1, (const xmlChar *)"l"));
	cout <<" cutoff "<<cutoff<<" L "<<symmetry<<endl;
	if(symmetry!="s"){
	  cout <<"need L=s in the semilocal section"<<endl; exit(-1);
	}
	xmlNodePtr cur2=cur1->children;
	while(cur2 != NULL)
	  {
	    xmlNodePtr cur3=cur2->children;
	    
	    while(cur3 != NULL)
	      {
		string cname((const char*)cur3->name);
		if(cname == "data"){
		  putContent(pos,cur3);
		}
		cur3 = cur3->next;
	      }
	    cur2 = cur2->next;
	  }
      }
      cur1 = cur1->next;
    }
  
    
  


}




void get_pseudo_from_the_grid(string filename, Spline_pseudo_writer & pseudo){

  cout <<endl;
  cout <<" Getting the external potential from the grid"<<endl;

  pseudo.psp_pos.resize(1);
  pseudo.psp_val.resize(1);
  double cutoff;

  //default values
  /*
  pseudo.psp_pos[0].resize(1);
  pseudo.psp_pos[0][0]=0.0;
  pseudo.psp_val[0].resize(1);
  pseudo.psp_val[0][0]=0.0;
  */

  // build an XML tree from a the file;
  xmlDoc *doc = NULL;
  
  // LIBXML_TEST_VERSION
  const char *psysfilename;
  psysfilename=filename.c_str();
  if ((doc = xmlReadFile(psysfilename, NULL, 0)) == NULL){
    printf("error: could not parse file %s\n", psysfilename);
    exit(-1); 
  }
  else{
    cout <<" Vext filename is "<<filename<<endl;
  }

  
  xmlNode *cur= NULL;
  cur = xmlDocGetRootElement(doc);
  if (xmlStrcmp(cur->name, (const xmlChar *) "pseudo")){
    cout <<"The "<<filename<<" of type xml does not contain pseudo system!"<<endl; exit(-1);
  }
  cur = cur->xmlChildrenNode;
  while (cur != NULL) {
    //cout <<"cur->name "<<cur->name<<endl;
    if(!xmlStrcmp(cur->name, (const xmlChar *) "grid")){
      read_grid(cur, pseudo.psp_pos[0]);
      pseudo.psp_val[0].resize(pseudo.psp_pos[0].size());
    }
    if(!xmlStrcmp(cur->name, (const xmlChar *) "semilocal"))
      read_semilocal(cur, pseudo.psp_val[0], cutoff);
    cur = cur->next;
  } 
  xmlFreeDoc(doc);
  xmlCleanupParser();


  //adjust for V*r format
  for(int i=0;i<pseudo.psp_val[0].size();i++)
    pseudo.psp_val[0][i]/=pseudo.psp_pos[0][i];


  //extrapolate
  double cusp=0.0;
  double spacing=0.02;
  extrapolate_to_linear_grid(pseudo.psp_pos[0], pseudo.psp_val[0],cusp,spacing);
  

  
  //apply cuttoff
  int i=0;
  while(pseudo.psp_pos[0][i]<= cutoff){
    //cout <<" pseudo.psp_pos[0][i] "<<pseudo.psp_pos[0][i]<<endl;
    i++;
  }
  if(i+1<pseudo.psp_pos[0].size()){
    pseudo.psp_pos[0].resize(i+1);
    pseudo.psp_val[0].resize(i+1);  
  }

  //cout <<"i "<<i<<" pos "<<pseudo.psp_pos[0][i-1]<<" val "<<pseudo.psp_val[0][i-1]<<endl;
}




int main(int argc, char ** argv) {
  string infilename, outputname;
  int debug=0;
  
  if(argc >= 2) {
    infilename=argv[argc-1];
  }
  else { usage(argv[0]); }

  for(int i=1; i< argc-1; i++) {
    if(!strcmp(argv[i], "-o") && argc > i+1) {
      outputname=argv[++i];
    }
    else if(!strcmp(argv[i], "-debug")) {
     debug=1;
    }
    else {
      usage(argv[0]);
    }
  }


  if(outputname == "") {
    outputname=infilename;
  }

  vector <Spline_pseudo_writer> pseudo;
  vector <Atom> atoms;
  Slat_wf_writer slwriter;
  vector <double> occupation; //electronic occupation
  
  string sysfilename=infilename+".qmc.xml";
  string vextonthegrid=infilename+".Vext.xml";
  string wffilename=infilename+".orb.dat";
 

  //
  int norbs;
  int nelectrons=0;
  int nup;
  int ndown;
  int basissize=0;
  int natoms;
  int total_charge=0;
  int found_sqd=0;
  double eref;

  
  // build an XML tree from a the file;
  xmlDoc *doc = NULL;
 
  // LIBXML_TEST_VERSION
  const char *psysfilename;
  psysfilename=sysfilename.c_str();
  if ((doc = xmlReadFile(psysfilename, NULL, 0)) == NULL){
    printf("error: could not parse file %s\n", psysfilename);
    exit(-1); 
  }
  else{
    cout <<"######### SQD2QMC by M.B. ########"<<endl;
    cout <<" converting SQD (qmcpack format) file "<<sysfilename<<" to QWalk format"<<endl;
    cout <<endl;

  }
  
  xml_system_data  qmc_xml;
  parse_from_xml(doc,qmc_xml);
  xmlFreeDoc(doc);
  xmlCleanupParser();

  //cout<<" stored charge variable "<<qmc_xml.p[0].g[0].charge<<endl; 
 

  for(int i=0;i<qmc_xml.p.size();i++){
    //cout << " p.name "<< qmc_xml.p[i].name<<" p.size "<< qmc_xml.p[i].size <<endl;
    for(int j=0;j<qmc_xml.p[i].g.size();j++){
      //cout << " g.name "<< qmc_xml.p[i].g[j].name<<" g.size "
      //   << qmc_xml.p[i].g[j].size <<" charge "<<qmc_xml.p[i].g[j].charge<<endl;

      if(strcmp(qmc_xml.p[i].g[j].name,"sqd")==0 && found_sqd==0){
	found_sqd=1;
	if(qmc_xml.p[i].size!=1)
	  {cout <<"sqd should have only one center, exit"<<endl; exit(-1);}
	natoms=1;
	Atom tmpatom;
	tmpatom.pos[0]=tmpatom.pos[1]=tmpatom.pos[2]=0.0;
	tmpatom.charge=qmc_xml.p[i].g[0].charge;
	tmpatom.name=outputname;
	tmpatom.name+="_origin";
	pseudo.resize(1);
	pseudo[0].label=tmpatom.name;
	atoms.push_back(tmpatom);
	slwriter.write_centers=false;
	slwriter.use_global_centers=false;
      }
      if(qmc_xml.p[i].g[j].size>0)
	total_charge+=qmc_xml.p[i].g[j].charge*qmc_xml.p[i].g[j].size;
      else
	total_charge+=qmc_xml.p[i].g[j].charge;
      if(strcmp(qmc_xml.p[i].name,"e")==0)
	if(strcmp(qmc_xml.p[i].g[j].name,"u")==0)
	  nup=qmc_xml.p[i].g[j].size;
	else if (strcmp(qmc_xml.p[i].g[j].name,"d")==0)
	  ndown=qmc_xml.p[i].g[j].size;
    }
  }


  if(!found_sqd) 
    {cout <<"did not find sqd system, exit"<<endl; exit(-1);}

  nelectrons=nup+ndown;
  
  cout <<" Sqd charge "<<atoms[0].charge<<", number of up e- "<<nup<<" and down e- "<<ndown
       <<" so the total charge is "<<total_charge<<endl;
  
  slwriter.mo_matrix_type="BASFUNC_MO";
  slwriter.magnification=1.0;
  slwriter.nup=nup;
  slwriter.ndown=ndown;
  
  

  vector < vector <int> > occ_up(1);
  vector < vector <int> > occ_down(1);
  

  string basisoutname=outputname+".basis";
  slwriter.basisname=basisoutname;

  ofstream basisout(basisoutname.c_str());
  string offset="  ";
  basisout<< "BASIS { "<<endl;
  basisout<< offset<<atoms[0].name <<endl;
  basisout<< offset<<"AOSPLINE" <<endl;
  basisout<< offset<<"NORENORMALIZE"<<endl;
  basisout<< offset<<"CUSP 0.0" <<endl;
  writebasis(basisout, qmc_xml, outputname, offset);
  basisout<< "      } "<<endl;
  basisout.close();
  




  //get the pseudopotential
  get_pseudo_from_the_grid(vextonthegrid,pseudo[0]);
  
  string slateroutname=outputname+".slater";
  getoccupations(qmc_xml,slwriter.calctype, occ_up[0], occ_down[0]);
  

  double weight=1.0;
  slwriter.detwt.push_back(weight);
  slwriter.occ_up.push_back(occ_up[0]);
  slwriter.occ_down.push_back(occ_down[0]);
   
  int largest_orb=0; 
  cout <<" Spin-up occupation "<<endl;
  for (int i=0;i<occ_up[0].size();i++){
    cout <<"  "<<occ_up[0][i]<<" ";
    if(occ_up[0][i]> largest_orb)
      largest_orb=occ_up[0][i];
  }
  cout <<endl;

   cout <<" Spin-down occupation "<<endl;
  for (int i=0;i<occ_down[0].size();i++){
    cout <<"  "<<occ_down[0][i]<<" ";
    if(occ_down[0][i]> largest_orb)
      largest_orb=occ_down[0][i];
  }
  cout <<endl;

  if(slwriter.calctype=="UHF"){
    slwriter.spin_dwn_start=largest_orb-occ_down[0].size();
  }

  //all the outputs will be writen here
  string orboutname=outputname+".orb";
  slwriter.orbname=orboutname;

  ofstream orbout(orboutname.c_str());
  for(int orbs=0;orbs<largest_orb;orbs++)
    orbout<<orbs+1<<"  1"<<endl;
  orbout.close();

  
  ofstream slaterout(slateroutname.c_str());
  slwriter.print_wavefunction(slaterout);
  slaterout.close();
    
    
  //--------------------------Jastrow 2 output
    
  
  
  string jast2outname=outputname+".jast2";
  Jastrow2_wf_writer jast2writer;
  jast2writer.set_atoms(atoms);
    
  double basis_cutoff=7.5; //arbitrary cutoff
  ofstream jast2out(jast2outname.c_str());
  print_std_jastrow2(jast2writer, jast2out, basis_cutoff);
  jast2out.close();
    
    
  //-----------------------------------System output
    
  string sysoutname=outputname+".sys";
  ofstream sysout(sysoutname.c_str());
  sysout << "SYSTEM { MOLECULE \n";
  sysout << "  NSPIN { " << slwriter.nup << "  "
	 << slwriter.ndown << " } \n";
  
      
  for(vector <Atom>::iterator at=atoms.begin(); at != atoms.end(); at++) {
    at->print_atom(sysout);
  }
  sysout << "}\n\n\n";
    
    
  int npsp=pseudo.size();
  for(int psp=0; psp < npsp; psp++) {
    pseudo[psp].print_pseudo(sysout);
  }

  sysout.close();
    
  //---------------------------HF and OPT outputs
    
  string hfoutname=outputname+".hf";
  ofstream hfout(hfoutname.c_str());
  print_vmc_section(hfout, outputname, eref);
  hfout << "\n\n";
  hfout << "INCLUDE " << sysoutname << "  \n";
  hfout << "TRIALFUNC { INCLUDE " << slateroutname << "}\n\n";
  hfout.close();
    
  string optoutname=outputname+".opt";
  ofstream optout(optoutname.c_str());
  optout << "\n\n";
  optout << "INCLUDE " << sysoutname << " \n";
  optout << "TRIALFUNC { \n  SLATER-JASTROW \n"
	 << "  WF1 { INCLUDE " << slateroutname << " } \n"
	 << "  WF2 { INCLUDE " << jast2outname   << " } \n"
	 << "}\n\n";
  print_vmc_section(optout, outputname, eref);
  print_opt_section(optout, outputname, eref);
  print_vmc_section(optout, outputname, eref);
  print_dmc_section(optout, outputname, eref);
  optout.close();

  
}

