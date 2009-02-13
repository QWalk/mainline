/*
 
Copyright (C) 2007 Lucas K. Wagner

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

#ifndef CBASIS_FUNCTION_H_INCLUDED
#define CBASIS_FUNCTION_H_INCLUDED

#include <string>
#include <vector>
#include <iostream>
#include <iomanip>
#include <complex>
#include <cstdio>
#include <fstream>
using namespace std;


#ifndef I
const complex <double> I(0.0,1.0);
#endif

class CBasis_function
{
public:


  virtual ~CBasis_function()
  {}

  virtual int read(string & filename, int &, int &, int &)=0;
  virtual int nfunc()=0;
  virtual double cutoff(int )=0;
  virtual string label()=0;
  virtual int showinfo(string & indent, ostream & os)
  {
    os << indent 
       <<  "Showing basis function info..this basis function doesn't"
      " have the showinfo function\n";
    return 1;
  }
  virtual int writeinput(string &, ostream &, int &, vector < vector <double > > &, string &)=0;
  virtual void calcVal(const vector <double> & r,
	       complex <double> & symvals,
	       int & orb)=0;
  virtual void getKvectors(vector < vector <double > > &, vector < vector <double > > &)=0;
};

//int allocate(vector <string> & basistext, CBasis_function * & bptr);
//int deallocate(CBasis_function * & bptr);

class Blochwave_function: public CBasis_function
{
public:

  Blochwave_function()
  {}
  ~Blochwave_function()
  {}

  //-----------------------------------------------------


  virtual int read(string & filename, int &, int &, int &);
  int nfunc();
  double cutoff(int )
  {
    return 1e99;
  }
  virtual string label()
  {
    return centername;
  }
  int showinfo(string & indent, ostream & os);
  int writeinput(string &, ostream &, int &, vector < vector <double > > &, string &);

  void calcVal(const vector <double> & r,
	       complex <double> & symvals,
	       int & orb);

  void getKvectors(vector < vector < double > > & primlatvec, vector < vector < double > > & VEC);

private:

  vector < vector < double > > g_vector;
  vector < vector < double > > k_vector_final;
  vector <int>  bmax_final; // number of bands for all k points
  vector < vector < double > > eigenvalues;
  vector < vector < vector < complex < double > > > > ckg_final; // complex Ck(g) coefficients
  vector <int> ireducible_kvec;
  vector <int> conjugate_kvec;
  vector <double> eigenvalueslinear;
  vector <int> list;
  vector < vector <int> > nonzero_part_of_ireducible_kvec;
  int nmax; // number of G vectors
  int kmax_final; // number of k vectors
  int bn_total_final; // total number of Bloch waves
  string centername;
  vector <int> orb_to_kn;
  vector <int> orb_to_bn;

};






#endif //CBASIS_FUNCTION_H_INCLUDED
//--------------------------------------------------------------------------
