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
void append_number(string &, double num);

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
//-----------------------------------------------------------------------------
//A pair of functions for storing and retrieving configurations from all 
//nodes to/from a single file. 
//ConfigType should have the following functions:
//read(istream & is)
//write(ostream & os)
//mpiSend(int node)
//mpiReceive(int node)
//The read() function should be particularly careful not to read past the 
//end of its section; otherwise the retrieval will not go well.
template <class ConfigType> void write_configurations(string & filename, 
                                                      Array1 <ConfigType> & configs) { 
  int nconfigs=configs.GetDim(0);
  time_t starttime;
  time(&starttime);
  single_write(cout,"Writing configurations..\n");  
  string tmpfilename=filename+".qw_tomove";
  wait_turn();
  if(mpi_info.node==0) { remove(tmpfilename.c_str()); } 
  ofstream os(tmpfilename.c_str(),ios::app);
  os.precision(15);
  if(!os) { error("Could not open ", tmpfilename); } 
  for(int i=0; i< nconfigs; i++) {
    os << " walker { \n";
    configs(i).write(os);
    os << "} \n";
  }
  os.close();
  finish_turn();
  /*
  if(mpi_info.node==0) { 
    ofstream os;
    os.precision(15);
    os.open(tmpfilename.c_str());
    if(!os) { error("Could not open ", tmpfilename); } 
    for(int i=0; i< nconfigs; i++) { 
      os << "walker { \n";
      configs(i).write(os);
      os << "} \n";
    }
    time_t midtime;
    time(&midtime);
    single_write(cout, "first write took ", difftime(midtime,starttime), " seconds \n");
#ifdef USE_MPI
    ConfigType tmpconf;
    for(int node=1; node < mpi_info.nprocs; node++) { 
      int nconf_node;
      MPI_Status status;
      int dummy=1;
      MPI_Send(&dummy, 1,MPI_INT, node, 0, MPI_Comm_grp);
      MPI_Recv(&nconf_node, 1, MPI_INT,
               node, 0, MPI_Comm_grp, &status);
      //cout << mpi_info.node << ": receiving " << nconf_node << " from " << node << endl;
      for(int i=0; i< nconf_node; i++) { 
        tmpconf.mpiReceive(node);
        os << "walker { \n";
        tmpconf.write(os);
        os << "}\n";
      }
      //time(&midtime);
      //single_write(cout, "writing to ", node, " took ");
      //single_write(cout, difftime(midtime, starttime), " seconds\n");
    }
#endif
  }
  else { 
#ifdef USE_MPI
    int dummy=1;
    MPI_Status status;
    //Wait to send our configurations so the event queue doesn't get filled up for
    //large numbers of processors.
    MPI_Recv(&dummy,1,MPI_INT, 0, 0, MPI_Comm_grp, &status);
    MPI_Send(&nconfigs, 1, MPI_INT,0, 0, MPI_Comm_grp);
    for(int i=0; i< nconfigs; i++) {
      configs(i).mpiSend(0);
    }
#endif
  }
  */
  if(mpi_info.node==0) 
    rename(tmpfilename.c_str(), filename.c_str());
  time_t endtime;
  time(&endtime);
  single_write(cout, "Write took ", difftime(endtime, starttime), " seconds\n");
}

//Reads configurations from the file and gives an array with the configurations 
//for this 
template <class ConfigType> void read_configurations(string & filename, 
                                                     Array1 <ConfigType> & configs) { 
  cout << "read_configurations" << mpi_info.node << endl;
  vector <ConfigType> allconfigs; 
  time_t starttime;  time(&starttime);
  if(mpi_info.node==0) { 
    ifstream is(filename.c_str());
    if(!is) { error("Could not open ", filename); } 
    string dummy;
    //This may be rough on memory, but gets deleted quickly and is relatively
    //easy to implement
    ConfigType tmpconf;
    while(is >> dummy) { 
      if(caseless_eq(dummy, "walker")) { 
        is >> dummy;
        if(dummy!="{") error("expected { in read_configurations()");
        tmpconf.read(is);
        allconfigs.push_back(tmpconf);
        //is >> dummy;
        //if(dummy!="}") error("expected } in read_configurations()");
      }
    }
    //cout << mpi_info.node << ": Done reading " << endl;
  }
  
  
#ifdef USE_MPI
  if(mpi_info.node==0) { 
    int totconfs=allconfigs.size();
    if(totconfs%mpi_info.nprocs != 0) {
      error("Non-even number of walkers/node.  Change the number of processors to something that divides ",totconfs);
    }
    int n_per_node=totconfs/mpi_info.nprocs;
    //cout << mpi_info.node << ": n_per_node " << n_per_node << endl;
    configs.Resize(n_per_node);
    
    int count=0;
    for(int i=0; i< n_per_node; i++) { 
      configs(i)=allconfigs[count++];
    }
    //cout << mpi_info.node << ": ......" << endl;
    for(int node=1; node < mpi_info.nprocs; node++) { 
      MPI_Send(&n_per_node, 1, MPI_INT,node, 0, MPI_Comm_grp);
      for(int i=0; i< n_per_node; i++) { 
        //cout << mpi_info.node << ": count " << count << endl;
        allconfigs[count++].mpiSend(node);
      }
    }
    //cout << mpi_info.node << ": done sending "<< count << endl;
  }
  else { 
    int nconf_node;
    MPI_Status status;
    MPI_Recv(&nconf_node, 1, MPI_INT,
             0, 0, MPI_Comm_grp, &status);
    //cout << mpi_info.node << ": receiving " << nconf_node << "configurations " << endl;
    configs.Resize(nconf_node);
    for(int i=0; i< nconf_node; i++) { 
      configs(i).mpiReceive(0);
    }
    //cout << mpi_info.node << ": done receiving " << endl;
  }
#else
  int totconfs=allconfigs.size();
  //cout << mpi_info.node << ": totconfs " << totconfs << endl;
  configs.Resize(totconfs);
  for(int i=0; i< totconfs; i++) {
    configs(i)=allconfigs[i];
  }
#endif //USE_MPI

  time_t endtime;
  time(&endtime);
  single_write(cout, "Read took ", difftime(endtime, starttime), " seconds\n");  
  
}


#endif  //QMC_IO_H_INCLUDED

//--------------------------------------------------------------------------
