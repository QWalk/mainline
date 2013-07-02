#include <iostream>
#include <vector>
#include <cstring>
#include <fstream>
#include <cmath>
#include <cassert>

using namespace std;

class Cube_info { 
public:
  vector <double> origin;
  vector<int> n;
  vector<double> resolution;
  vector<int> atoms;
  vector < vector <double> > atompos; //0, x,y,z
  int nsamples;

  string line1;
  string line2;
  vector <double> density;

  void read_cube(istream & is);
  void write_cube(ostream & os);
  void smooth();
  void normalize();
  void enhance(); //!< make the largest value=1, for plotting
  void add(Cube_info &, double fac );
  void fold(int , bool); //!< fold into a box 1/nth the size.  The bool controls whether then to re-expand it into the original size
  double dens(int x,int y, int z) const { 
    assert(x<n[0]); 
    assert(y<n[1]);
    assert(z<n[2]);
    int p=x*n[1]*n[2]+y*n[2]+z;
    assert(p < density.size());
    return density[p];
  }

  double & dens(int x,int y, int z) { 
    //cout << "dens " << x << " " << y << " " << z <<" " << n[0] << " " << n[1] << " " << n[2] << endl;
    assert(x<n[0]); 
    assert(y<n[1]);
    assert(z<n[2]);
    return density[x*n[1]*n[2]+y*n[2]+z];
  } 

  double volume_voxel() { 
    return resolution[0]*resolution[1]*resolution[2];
  }

  void write_projection(vector <double> & proj, ostream & os);
  void write_plane_proj(int d, ostream & os);
  

};


void usage() { 
  cout << "smooth_cube [options] < input.cube > output.cube\n";
  cout << "------------------------------------------------\n";
  cout << "Where [options] can be any of the following:\n";
  cout << "-smooth                  Smooth the file to remove jaggies\n";
  cout << "-sub  [cubefile]         Subtract the given cube file from the one from the input \n";
  cout << "                         (both are normalized before this is done)\n";
  cout << "-add  [cubefile]         Add the two cube files together (also normalizes here)\n";
  cout << "-fold  [integer]         Fold the density back onto itself.  Useful for supercells of nxnxn\n";
  cout << "                         Currently only works for orthorhombic cells\n";
  cout << "-proj  [x y z]           Project on the line given by [x,y,z].  \n";
  cout << "-noenhance               Don't change the norm to make it easier to visualize. \n";
  cout << "-normalize               Normalize the cube file to have sum 1\n"; 
  cout << "\nIn all cases, the cube file will be 'enhanced' to make the largest value equal to one.\n"
       << "This presumably makes the field easier to deal with in a visualization program\n";
  cout << "\nThe options are interpreted in order, so one can do several smooths by putting \n"
       << "-smooth -smooth, or one can smooth,  then subtract, etc.  Some things are\n"
       << "better ideas than others, of course.\n";
}


int main(int argc, char ** argv) { 
  for(int i=1; i< argc; i++) { 
    if(!strcmp(argv[i],"-h") || !strcmp(argv[i],"--help")) { 
      usage();
      exit(0);
    }
  }
  bool no_write=false;
  Cube_info cube;
  cube.read_cube(cin);
  //cube.normalize();
  int enhance=1;
  for(int i=1; i < argc; i++) {
    if(!strcmp(argv[i],"-noenhance")) enhance=0;
    if(!strcmp(argv[i],"-add") && i<argc+1) {
      cerr << "adding " << argv[i+1] << endl;
      ifstream is(argv[++i]);
      if(!is) {
        cerr << "Couldn't open " << argv[i]
          << endl;
      }
      Cube_info ncube;
      ncube.read_cube(is);
      cerr << "normalizing " << endl;
      cube.normalize();
      ncube.normalize();
      cerr << "adding " << endl;
      cube.add(ncube, 1.0);
    }
    if(!strcmp(argv[i],"-sub") && i<argc+1) {
      cerr << "subtracting " << argv[i+1] << endl;
      ifstream is(argv[++i]);
      if(!is) {
        cerr << "Couldn't open " << argv[i]
          << endl;
      }
      Cube_info ncube;
      ncube.read_cube(is);
      cerr << "normalizing " << endl;
      cube.normalize();
      ncube.normalize();
      cerr << "adding " << endl;
      cube.add(ncube, -1.0);
    }    

    if(!strcmp(argv[i],"-fold") && i < argc+1) {
      cerr << "folding " << argv[i+1] << endl;
      int nfold=atoi(argv[++i]);
      cube.fold(nfold, false);
    }

    if(!strcmp(argv[i],"-foldavg") && i < argc+1) {
      cerr << "fold-averaging " << argv[i+1] << endl;
      int nfold=atoi(argv[++i]);
      cube.fold(nfold, true);
    }

    if(!strcmp(argv[i], "-smooth"))
      cube.smooth();
    if(!strcmp(argv[i], "-normalize"))
      cube.normalize();
    if(!strcmp(argv[i], "-proj")) { 
      if(i >= argc+4) {
        cerr << "proj needs 4 arguments" << endl;
        exit(1);
      }
      vector <double> proj(3);
      proj[0]=atof(argv[++i]);
      proj[1]=atof(argv[++i]);
      proj[2]=atof(argv[++i]);
      ofstream projout(argv[++i]);
      cube.write_projection(proj,projout);
    }
    if(!strcmp(argv[i], "-xy")) { 
      cube.write_plane_proj(2, cout);
      no_write=true;
    }
    

  }
  if(no_write) return 0;
  if(enhance) { 
    cerr << "enhancing " << endl;
    cube.enhance();
  }
  cube.write_cube(cout);
  return 0;

}

//----------------------------------------------------------------------

void Cube_info::read_cube(istream & is) {
  getline(is,line1);
  nsamples=1;
  getline(is,line2);
  int natoms;
  is >> natoms;
  double dumdouble;
  for(int i=0; i< 3; i++) { 
    is >> dumdouble; 
    origin.push_back(dumdouble);
  }

  int dumint;
  for(int i=0; i< 3; i++) {
    is >> dumint;
    n.push_back(dumint);
    for(int j=0; j< 3; j++) {
      is >> dumdouble;
      if(i==j) resolution.push_back(dumdouble);
    }
  }

  vector <double> dumvec(4);
  for(int at=0; at< natoms; at++) {
    is >> dumint;
    atoms.push_back(dumint);
    
    for(int j=0; j< 4; j++) {
      is >> dumvec[j];
    }
    atompos.push_back(dumvec);
  }

  int npoints=n[0]*n[1]*n[2];
  //cout << "natoms " << natoms 
  //     << " npoints " << npoints << endl;

  if(!is) { 
    cerr << "input file too short!" << endl;
    exit(1);
  }

  density.reserve(npoints);
  for(int i=0; i< npoints; i++){
    is >> dumdouble;
    density.push_back(dumdouble);
  }
     
  
}

//----------------------------------------------------------------------

void Cube_info::write_cube(ostream & os) {
  os << line1 << endl;
  os << line2 << endl;
  int natoms=atoms.size();
  os << natoms << "  ";
  double dumdouble;
  for(int i=0; i< 3; i++) {
    os << origin[i] << "  ";
  }
  os << endl;

  int dumint;
  for(int i=0; i< 3; i++) {
    os << n[i] << "  ";
    for(int j=0; j< 3; j++) {
      if(i==j) os << resolution[i] << "  ";
      else os << " 0.0000  ";
    }
    os << endl;
  }

  for(int at=0; at< natoms; at++) {
    os << atoms[at] << "  ";
    for(int j=0; j< 4; j++) {
      os << atompos[at][j] << "  ";
    }
    os << endl;
  }

  int count=0;
  for(vector<double>::iterator i=density.begin();
      i!=density.end(); i++) {
    os << *i << "  ";
    if((++count)%6==0) os << endl;
  }
  os << endl;

}


//----------------------------------------------------------------------

void Cube_info::smooth() {
  int npoints=n[0]*n[1]*n[2];
  vector <double> newdens(npoints);
  int count=0;
  double wt=1/7.0;
  for(int x=0; x< n[0]; x++) { 
    for(int y=0; y < n[1]; y++) { 
      for(int z=0; z < n[2]; z++) {
        double tot=0;
        double pt=dens(x,y,z);
        tot+=wt*pt;
        if(x>0) tot+=wt*dens(x-1,y,z);
        else tot+=wt*pt;
        if(x < n[0]-1) tot+=wt*dens(x+1,y,z);
        else tot+=wt*pt;
        if(y >0) tot+=wt*dens(x,y-1,z);
        else tot+=wt*pt;
        if(y < n[1]-1) tot+=wt*dens(x,y+1,z);
        else tot+=wt*pt;
        if(z > 0) tot+=wt*dens(x,y,z-1);
        else tot+=wt*pt;
        if(z < n[2]-1) tot+=wt*dens(x,y,z+1);
        else tot+=wt*pt;
        newdens[count++]=tot;
      }
    }
  }
  
  count=0;
  for(int x=0; x< n[0]; x++) { 
    for(int y=0; y< n[1] ; y++) { 
      for(int z=0; z< n[2]; z++) { 
        dens(x,y,z)=newdens[count++];
      }
    }
  }

}

//----------------------------------------------------------------------

void Cube_info::normalize() {
  double sum=0.0;
  double vol=volume_voxel();
  for(vector<double>::iterator i=density.begin();
      i!=density.end(); i++) {
    sum+=*i;
  }
  sum*=vol;

  for(vector<double>::iterator i=density.begin();
      i!=density.end(); i++) {
    *i/=sum;
  }  

}

//----------------------------------------------------------------------
void Cube_info::enhance() { 
  double max=0.0;
  for(vector<double>::iterator i=density.begin();
      i!=density.end(); i++) {
    if(*i>max) max=*i;
  }

  for(vector<double>::iterator i=density.begin();
      i!=density.end(); i++) {
    *i/=max;
  }
  
}
  
//----------------------------------------------------------------------
void Cube_info::add(Cube_info & ncube, double fac) {
  for(int d=0; d< 3; d++) {
    if(n[d]!=ncube.n[d]) {
      cerr << "Error: adding cube, but dimensions don't match"
        << ":d= " << d << " n[d] " << n[d] << " ncube.n[d] "
        << ncube.n[d] 
        << endl;
      exit(1);
    }
  }

  if(density.size()!=ncube.density.size()) {
    cerr << "densities aren't the same size in add()"
      << endl;
    exit(1);
  }

  vector <double>::iterator my=density.begin();
  for(vector<double>::iterator add=ncube.density.begin();
      add !=ncube.density.end(); add++) {
    *my+=fac*(*add);
    my++;
  }

}

//----------------------------------------------------------------------

void Cube_info::fold(int fac, bool avg) { 
  vector <int> nn(3);
  for(int d=0; d< 3;d++) nn[d]=n[d]/fac;
  int npoints=nn[0]*nn[1]*nn[2];
  vector <double> ndens(npoints);
  for(vector<double>::iterator i=ndens.begin();
      i!= ndens.end(); i++) *i=0.0;

  vector <double>::iterator my=ndens.begin();
  for(int x=0; x< nn[0]; x++) {
    for(int y=0; y< nn[1]; y++) { 
      for(int z=0; z< nn[2]; z++) {

	for(int ii=0; ii< fac; ii++){
	  for(int jj=0; jj< fac; jj++) {
	    for(int kk=0; kk< fac; kk++) {
	      *my+=dens(x+ii*nn[0],y+jj*nn[1],z+kk*nn[2])/(fac*fac*fac);
	    }
	  }
	}
	my++;
      }
    }
  }

  if(!avg) { 
    n=nn;
    density=ndens;
  }
  else { 
    my=ndens.begin();
    for(int x=0; x< nn[0]; x++) {
      for(int y=0; y< nn[1]; y++) { 
	for(int z=0; z< nn[2]; z++) {
	  
	  for(int ii=0; ii< fac; ii++){
	    for(int jj=0; jj< fac; jj++) {
	      for(int kk=0; kk< fac; kk++) {
		dens(x+ii*nn[0],y+jj*nn[1],z+kk*nn[2])=*my;
	      }
	    }
	  }
	  my++;
	}
      }
    }
  }

}

//----------------------------------------------------------------------

void normalize(vector <double> & x) { 
  assert(x.size()==3);
  double tot=0;
  for(int i=0; i< 3; i++) tot+=x[i]*x[i];
  tot=sqrt(tot);
  for(int i=0; i< 3; i++) x[i]/=tot;
}

double dot(vector <double> & x, vector <double> & y) { 
  assert(x.size()==3);
  assert(y.size()==3);
  double t=0;
  for(int i=0; i< 3; i++) t+=x[i]*y[i];
  return t;
}
//----------------------------------------------------------------------
void Cube_info::write_projection(vector <double> & proj, ostream & os) { 
  
  ::normalize(proj);
  double min=dot(origin,proj);
  vector <double> end(3);
  for(int i=0; i< 3; i++) end[i]=origin[i]+resolution[i]*n[i];
  double max=dot(end,proj);
  double res=resolution[2];
  assert(max > min);
  int npts=int((max-min)/res)+2;
  //don't use max from here on..

  vector <double> projdens(npts);

  for(vector<double>::iterator p=projdens.begin();
      p != projdens.end(); p++) *p=0.0;
  


  vector <double> pt(3);
  for(int x=0; x< n[0]; x++) { 
    for(int y=0; y< n[1]; y++) { 
      for(int z=0; z< n[2]; z++) { 
        pt[0]=x*resolution[0]+origin[0];
        pt[1]=y*resolution[1]+origin[1];
        pt[2]=z*resolution[2]+origin[2];
        double pos=dot(pt,proj);
        //assert(pos < max);
        int p=int((pos-min)/res);
        p=z; //hardcode for 001, since the other stuff is noisy..	
        assert(p<npts && p >=0);

        projdens[p]+=dens(x,y,z);

      }
    }
  }

  double norm=0;
  for(int i=0; i< npts; i++) { 
    norm+=projdens[i]*res;
  }


  for(int i=0; i< npts; i++) { 
    os << i*res+min << "   " << projdens[i]/norm << endl;
  }
  /*

  for(int z=0; z< n[2]; z++) { 
    double sum=0;
    for(int x=0; x< n[0]; x++) {
      for(int y=0; y < n[1]; y++) { 
	int p=x*n[1]*n[2]+y*n[2]+z;
	sum+=dens(x,y,z);
      }
    }
    os << z*resolution[2]+origin[2] << "   " 
	 << sum << endl;
  }
  */
}

//----------------------------------------------------------------------

void Cube_info::write_plane_proj(int d, ostream & os) { 
  vector <double> projdens(n[0]*n[1]);

  for(vector<double>::iterator p=projdens.begin();
      p != projdens.end(); p++) *p=0.0;
  vector <double> pt(3);
  for(int x=0; x< n[0]; x++) { 
    for(int y=0; y< n[1]; y++) { 
      for(int z=0; z< n[2]; z++) { 
//pt[0]=x*resolution[0]+origin[0];
//        pt[1]=y*resolution[1]+origin[1];
//        pt[2]=z*resolution[2]+origin[2];
//        double pos=dot(pt,proj);
        //assert(pos < max);
        
        //assert(p<npts && p >=0);

        projdens[x*n[0]+y]+=dens(x,y,z);

      }
    }
  }

  double norm=0;
  for(int i=0; i< n[0]; i++) { 
    for(int j=0; j< n[1]; j++) { 
      norm+=projdens[i*n[0]+j]*resolution[0]*resolution[1];
    }
  }


  for(int i=0; i< n[0]; i++) { 
    for(int j=0; j< n[1]; j++) { 
      os << i*resolution[0]+origin[0] << " "
        << j*resolution[1]+origin[1] << "   " << projdens[i*n[0]+j]/norm << endl;
    }
  }
}

