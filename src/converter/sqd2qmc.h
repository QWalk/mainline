#ifndef SQD2QMC_H_INCLUDED
#define SQD2QMC_H_INCLUDED
#include <iostream>
#include <sstream>
#include <iosfwd>
#include <string>
#include <vector>
#include <algorithm>
#include <libxml/xmlmemory.h>
#include <libxml/parser.h>
#include <libxml/xpath.h>
using namespace std;

class group {
 public:
  int charge;
  char* name;
  int size;
  group (){
    name="";
    charge=0;
    size=-1;
  }
};


class particleset {
 public:
  vector <group> g;
  char* name;
  int size;
  particleset() {
    name="";
    size=-1;
  }
};


class basisGroup{
 public: 
  char* rid;
  char* ds;
  int n;
  int l;
  int m;
  int s;
  basisGroup(){}
};

class atomicBasisSet{
 public:
  char* type;
  char* elementType;
  char* filename;
  vector <basisGroup> orbs; 
  atomicBasisSet(){}
};


class wavefunction{
 public:
  char* name;
  char *target;
  atomicBasisSet abs;
  wavefunction(){}
};



class xml_system_data {
public:
  vector <particleset> p;
  wavefunction wf;
  xml_system_data() {
  }
};


class radial_orbital {
  public:
  vector <double> radial_part;
  vector <double> uofr;
  vector <double> r_log;
  int power;
  vector <int> quantum_numbers;
  double eigenvalue;
  string orb_name;
  radial_orbital (){
    radial_part.resize(0);
    uofr.resize(0);
    r_log.resize(0);
    quantum_numbers.resize(0);
    eigenvalue=0.0;
    power=0;
    orb_name="";
  }
};



bool putContent(vector <double> & in, xmlNodePtr cur) {
  //istringstream stream((const char*)(xmlNodeListGetString(cur->doc, cur->xmlChildrenNode, 1)));
  //for(int i=0;i<in.size();i++){
  istringstream stream((const char*)(xmlNodeListGetString(cur->doc, cur->xmlChildrenNode, 1)));
  for(int i=0;i<in.size();i++){
    stream >> in[i];
    //cout <<in[i]<<endl;
  }
}
#endif
