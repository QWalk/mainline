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

#ifndef QMC_IO_H_INCLUDED
#define QMC_IO_H_INCLUDED


#include "Qmc_std.h"
#include <iomanip>

const string startsec="{";
const string endsec="}";
const string comment="#";


int caseless_eq(const string & , const string & );



/*!
This function reads the next word out of a vector of strings, doing error
correction and converting to int, double, or string, depending on what
you give it.  If stringstreams worked, it could be a template, which would
work with everything.  @see readnext
*/
void readnext(vector <string> & s,unsigned int & i, int & t);
void readnext(vector <string> & s,unsigned int & i, doublevar & t);
void readnext(vector <string> & s,unsigned int & i, string & t);

/*!
Starting from pos, finds the next instance of sectionname 
at the same logical
level(will not descend into {}'s) and puts the entire thing
into section

\todo
Perhaps we should have this and the other read functions remove
the keyword and values from the input.  That way, we wouldn't have 
to deal with pos, and we could check to make sure the input was 
completely used.

 */
int readsection(vector <string> & input,
                unsigned int & pos,
                vector <string> & section,
                const char * sectionname);

/*!
Starting from pos, finds the next instance of valuename, and
puts the next word into value.  Uses readnext to do so, so
it supports every type that readnext is overloaded to 
support.
 */
template <class T>
int readvalue(vector <string> & input,
              unsigned int & pos,
              T & value,
              const char * valuename)
{
  int level=0;  //how many brackets we've descended into.

  string name(valuename);
  for(; pos< input.size(); pos++)
  {
    if(input[pos]==startsec)
    {
      level++;
    }
    else if(input[pos]==endsec)
    {
      level--;
    }

    if(level < 0)
    {
      pos++;
      return 0;
    }

    if(level==0 && caseless_eq(input[pos],valuename)) //input[pos]==valuename)
    {
      readnext(input, pos, value);
      
      //input.erase(input.begin()+pos-1,input.begin()+ pos+1);
      //pos-=2;
      return 1;
    }
  }
  //if we get here, we didn't find a value
  return 0;
}

/*!
  Search for a keyword and return whether it's in the level
or not.
 */
int haskeyword(vector <string> & input, unsigned int & pos,
               const char * valuename);


/*!
Separates the given file into words, removing comments and
including any include files.
 */
void parsefile(ifstream & inputfile, vector <string> & words);
/*!
  Count the number of brackets in the words, and make sure that
they add up correctly.
 */
int checkbrackets(vector <string> & words);

/*!
Writes an array of objects to a stream, provided they have
a rawOutput function.
 */
template <class Type>
void write_array(int nfill, Array1 <Type *> & array, ostream & os)
{
  //cout << "write_array " << endl;
  if(!os) error("bad output stream in write_array");
  if(nfill > array.GetDim(0) )
  {
    error("While writing an array in write_array, nfill is ", nfill,
          "and the array size is ", array.GetDim(0));
  }
  os << "Array " << startsec <<" " << nfill << endl;

  for(int i=0; i< nfill; i++)
  { 
    //cout << " i " << i << endl;
    array(i)->rawOutput(os);
  }
  os << " "<< endsec << " \n";
}


/*!
 
 */
template <class Type>
int read_array(Array1 <Type *> & array, istream & is)
{
  string text;
  //cout << "****Entering read_array\n";
  int numobjects;
  is >> text;
  if(text != "Array")
    error("expected Array, got ", text);
  //cout << "first read " << text << endl;
  is >> text >> numobjects;
  if(text != startsec)
    error("expected ", startsec, " got ", text);

  if(numobjects > array.GetDim(0))
  {
    cout << "WARNING: in read_array, the number of objects in the file"
    " is greater than the number of objects passed. In file: " <<
    numobjects << "  passed: " << array.GetDim(0) << endl;
  }

  numobjects=min(numobjects, array.GetDim(0));

  for(int i=0; i< numobjects; i++)
  {
    array(i)->rawInput(is);
  }

  //now find the end of the array.
  int numsec=0;

  while(is >> text)
  {
    if(text == startsec)
    {
      numsec++;
    }
    else if(text==endsec)
    {
      numsec--;
    }

    if(numsec < 0)
    {
      return numobjects;
    }
  }
  error("Reached the end of file without finding closing bracket.");
  return 0;
}




/*!
Modifies the name of a file for parallel runs.  We should
always use this when writing to a file; otherwise, there 
will be collisions.(or use single_write)
 */
void canonical_filename(string & );
void canonical_filename(string &, int num);

void append_number(string &, int num);

/*!
 
 */
template <class T>
void single_write(ostream & os, const T & t)
{
#ifdef USE_MPI
  if(mpi_info.node==0)
  {
#endif
    os << t;
#ifdef USE_MPI
  }
#endif
  os.flush();
}
template <class T, class U>
void single_write(ostream & os, const T & t, const U & u)
{
#ifdef USE_MPI
  if(mpi_info.node==0)
  {
#endif
    os << t << u;
#ifdef USE_MPI
  }
#endif
  os.flush();
}

template <class T, class U, class V>
void single_write(ostream & os, const T & t, const U & u, const V & v)
{
#ifdef USE_MPI
  if(mpi_info.node==0)
  {
#endif
    os << t << u << v;
#ifdef USE_MPI
  }
#endif
  os.flush();
}

template <class T, class U, class V, class W>
void single_write(ostream & os, const T & t,
                  const U & u, const V & v, const W & w)
{
#ifdef USE_MPI
  if(mpi_info.node==0)
  {
#endif
    os << t << u << v << w;
#ifdef USE_MPI
  }
#endif
  os.flush();
}

//----------------------------------------------------------------------

template <class T>
void debug_write(ostream & os, const T & t)
{
#ifdef DEBUG_WRITE
#ifdef USE_MPI
  os << mpi_info.node << ": ";
#endif

  os << t;
  os.flush();
#endif
}
template <class T, class U>
void debug_write(ostream & os, const T & t, const U & u)
{
#ifdef DEBUG_WRITE
#ifdef USE_MPI
  os << mpi_info.node << ": ";
#endif

  os << t << u;
  os.flush();
#endif
}

template <class T, class U, class V>
void debug_write(ostream & os, const T & t, const U & u, const V & v)
{
#ifdef DEBUG_WRITE
#ifdef USE_MPI
  os << mpi_info.node << ": ";
#endif

  os << t << u << v;
  os.flush();
#endif
}

template <class T, class U, class V, class W>
void debug_write(ostream & os, const T & t,
                 const U & u, const V & v, const W & w)
{
#ifdef DEBUG_WRITE
#ifdef USE_MPI
  os << mpi_info.node << ": ";
#endif

  os << t << u << v << w;
  os.flush();
#endif
}

//#include <sstream>
/*  Unfortunately, some gcc libraries don't support this yet,
but it should be very nice when they do..
template <class T>
void readnext(vector <string> & s,unsigned int & i, T & t) {
	stringstream sstemp;

  if( i >= s.size() ) {
  	cout << "Error!  Unexpected end of file at " << s[i] << endl;
    Terminate();
  }
  cout << "readnext:\n";
  cout << s[i+1] << endl;
  sstemp << s[++i];
  sstemp >> t;
  cout << t << endl;
}

template <class T> void read(string & s, T & t) {
	stringstream sstemp;
  sstemp << s;
  sstemp >> t;
}
*/



#endif  //QMC_IO_H_INCLUDED

//--------------------------------------------------------------------------
