
#include <iostream>
#include <cstdio>
using namespace std;


double checksum(double * a, int n) { 
  //Here we just take the average of the array.
  double check=0;
  for(int i=0; i< n; i++) 
    check+=a[i]/n;
  return check;
}

int binary_write_checksum(double * a, int n, FILE * f) { 
  double n_d=n;
  double zero=0.0;
  fwrite(&zero,sizeof(double),1,f);
  fwrite(&zero,sizeof(double),1,f);
  fwrite(&n_d,sizeof(double),1,f);
  double check=checksum(a,n);
  fwrite(a,sizeof(double),n,f);
  fwrite(&check,sizeof(double),1,f);
  if(ferror(f)) { 
    return 0;
  }
  return 1;
}

int main(int argc, char ** argv) { 
  if(argc < 4) { 
    cout << "Swap from the old trace files to the new, checksummed ones" << endl;
    cout << "The old files are binary files with the electron positions and then the weights." <<endl;
    cout << "The new ones have a header for each walker, followed by a checksum. " << endl;
    cout << "Usage : " << argv[0] << " [nelectrons] infile outfile " << endl;
    return 1;
  }
  FILE * fin=fopen(argv[2],"r");
  FILE * fout=fopen(argv[3],"w");
  int n=atoi(argv[1])*3+1;
  double * tmp=new double[n];
  while(fread(tmp, sizeof(double),n,fin)) { 
    binary_write_checksum(tmp,n,fout);
  }
  fclose(fin);
  fclose(fout);
  delete [] tmp;
  return 0;
}

