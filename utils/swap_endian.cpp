
#include <iostream>
#include <cstdio>
using namespace std;
double swap(double d){
   double a;
   unsigned char *dst = (unsigned char *)&a;
   unsigned char *src = (unsigned char *)&d;

   dst[0] = src[7];
   dst[1] = src[6];
   dst[2] = src[5];
   dst[3] = src[4];
   dst[4] = src[3];
   dst[5] = src[2];
   dst[6] = src[1];
   dst[7] = src[0];

   return a;
}

int main(int argc, char ** argv) { 
  if(argc < 3) { 
    cout << "Swap from big-endian to little endian and vice-versa" << endl;
    cout << "Usage : " << argv[0] << " infile outfile " << endl;
    return 1;
  }
  FILE * fin=fopen(argv[1],"r");
  FILE * fout=fopen(argv[2],"w");
  double tmp;
  while(fread(&tmp, sizeof(double),1,fin)) { 
    double tmp2=swap(tmp);
    fwrite(&tmp2,sizeof(double),1,fout);
  }
  fclose(fin);
  fclose(fout);
  return 0;
}

