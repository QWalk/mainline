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


#include "Qmc_std.h"
#include "Program_options.h"
#include "Qmc_method.h"
#include "qmc_io.h"
#include <ctime>
#include "ulec.h"
#include <new>

using namespace std;


/*!
Sorts the text in the ifstream and puts the corresponding
values into options.
*/
void doinput(Program_options & options,
             ifstream & inputfile);

/*!
Decides where we want to go..  We'll assume that the first
argument is the file name of the input file.

\todo
I'd like to have it so that the input errors are accumulated
and then printed out after it's attempted to allocate everything.
Would be a nice feature.  Could be done through exceptions.

*/
int main(int argc, char* argv[])
{
#ifdef USE_MPI
  MPI_Init(&argc, &argv);
#endif

  int inputfilestart=1;
  for(int i=1; i< argc; i++) {
    string arg(argv[i]);
    if(caseless_eq(arg,"-rappture")) {
      global_options::rappture=1;
      inputfilestart=i+1;
    }
  }

  if ( argc <= inputfilestart )
    error("usage: ", argv[0], " [-rappture] filename(s)");

  // No. of indeendent processes in the pack = No. of inputfiles
  int processcount=argc-inputfilestart;
  int group;   // indexes groups of independent processes (needed to assign
               // correct inputfile)

#ifdef USE_MPI
  int nprocs, node, group_size, ingroup;
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &node);
  group_size=nprocs/processcount;
  if ( (nprocs-1)/group_size > processcount-1 )
    error("Number of CPUs not divisible by the number of inputfiles; \n don't know what to do with leftover nodes.");
  group=node/group_size;                    // modulo (selects group)
  ingroup=node%group_size;                  // remainder (gives ordering in the group)
  MPI_Comm_split(MPI_COMM_WORLD,group,ingroup,&MPI_Comm_grp);
  MPI_Comm_size(MPI_Comm_grp, &(mpi_info.nprocs));
  MPI_Comm_rank(MPI_Comm_grp, &(mpi_info.node));
#else
  group=0;
  if ( processcount>1 ) 
    error("More than one input file specified, but this is not a parallel (MPI) version.");
#endif

  vector<string> inputfiles(processcount);
  for(int i=0; i<processcount; i++) inputfiles[i]=argv[i+inputfilestart];

  // allocate per-process variables
  Program_options options;
  ifstream inputfile;
  string outputfile;

  inputfile.open(inputfiles[group].c_str());
  if(!inputfile) {
    error("Couldn't open ",inputfiles[group]);
  }
  options.runid=inputfiles[group];
  outputfile=inputfiles[group];
  outputfile+=".o";
  
  ofstream output;
  if(mpi_info.node==0)
  {
    debug_write(cout, "output to ", outputfile, "\n");
    output.open(outputfile.c_str());
  }

  time_t starttime;
  time(&starttime);
  const size_t buflen=80;
  char hostname[buflen];
  qmcgethostname(hostname, buflen);
//cout << "node " << mpi_info.node << " alive on " << hostname << endl;
  if(output)
  {
    output << "------------------------------------------------\n";
    output << "Quantum Walk development version\n";
    output << "Contributions from: (in no particular order) \n";
    output << "Michal Bajdich, Shuming Hu, Jindrich Kolorenc, Kevin Rasch, \n";
    output << "Pavel Vagner, Rene Derian, Paul Kent, Jarrod McClean, \n";
    output << "Fernando Reboredo, Lucas Wagner, Hiori Kino, ";
    output << "and Lubos Mitas.\n";
    output << "Originated at North Carolina State University \n";
    output << "Please cite J. Comp. Phys. v 228 pp 3390-3404 (2009) \n when publishing results from this program.\n";
    output << "------------------------------------------------\n";
    output << "Running on " << hostname
           << endl << ctime(&starttime) << endl;
  }

  doinput(options, inputfile);
  inputfile.close();

#ifdef USE_MPI
  long int is1, is2;
  rng.getseed(is1, is2);
  rng.seed(123754*mpi_info.node+is1, 234456*mpi_info.node+is2);
  rng.getseed(is1, is2);
#endif
  
  Qmc_method * method=NULL;
  for(unsigned int i=0; i< options.methodtext.size(); i++)
  {

#ifndef NO_EXCEPTIONS
    try {
#endif

    allocate(options.methodtext[i], options, method);
    method->generateVariables(options);
    method->showinfo(output);
    time_t methodstart;
    time(&methodstart);
    //clock_t run_start=clock();
    method->run(options, output);
    //clock_t run_end=clock();
    time_t endtime;
    time(&endtime);
    if(output)
    {
      //On Linux, this method overflows after only about an hour, so it's fairly useless.
      //Better to use difftime so the results are credible.
      //output << "CPU time for run only: " << double(run_end-run_start)/double(CLOCKS_PER_SEC) << endl;
      output << "Wall time for this method: " << difftime(endtime, methodstart) << " seconds." << endl;
      output << "Total wall time so far: " << difftime(endtime, starttime) << " seconds." << endl;
    }
    deallocate(method);
#ifdef USE_MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif

#ifndef NO_EXCEPTIONS
  }
  catch(bad_alloc) {
    error("Had an allocation error and couldn't handle it.  Perhaps the run is too large.");
#ifdef USE_MPI
    MPI_Abort(MPI_COMM_WORLD, 1);
#endif
  }
  catch(Qmc_error err) {
    cout << "Caught qmc_error" << endl;
#ifdef USE_MPI
    MPI_Abort(MPI_COMM_WORLD, 1);
#endif
    break;
  }
  catch(...) {
    cout << "caught a general exception " << endl;
#ifdef USE_MPI
    MPI_Abort(MPI_COMM_WORLD, 1);
#endif
  }

#endif



  }
  if(output)
  {
    output.close();
  }

#ifdef USE_MPI
  MPI_Finalize();
#endif

}

//--------------------------------------------------------------------------

void doinput(Program_options & options,
             ifstream & inputfile
            )
{

  //Atom tempatom;
  //string temp;
  vector <string> words;
  words.reserve(50);

  //cout << "parsefile\n";

  parsefile(inputfile, words);
  //cout << "parse done\n";

  checkbrackets(words);


  unsigned int pos=0;
  vector <string> atomsec;
  if(readsection(words, pos, atomsec, "ATOM") != 0) {
    error("There shouldn't be an ATOM section in the global space any "
          "more.  Please put them in the SYSTEM section.");
  }




  //cout << "atoms done\n";
  //Normal sections that will be parsed by the instantiated
  //class

  vector < string > pseudosec;
  pos=0;
  while( readsection(words, pos, pseudosec, "PSEUDO") != 0) {
    options.pseudotext.push_back(pseudosec);
  }


  vector < string > wfsec;
  pos=0;
  while( readsection(words, pos, wfsec, "TRIALFUNC") != 0)
  {
    options.twftext.insert(options.twftext.end(), wfsec);
  }

  vector <string> syssec;
  pos=0;
  while(readsection(words, pos, syssec, "SYSTEM") != 0)
  {
    options.systemtext.push_back(syssec);
  }
  if(options.systemtext.size()==0)
  {
    error("missing system section.  Try SYSTEM { MOLECULE } or similar");
  }
  //cout << "bah " << endl;

  vector <string> methsec;
  pos=0;
  while( readsection(words, pos, methsec, "METHOD") != 0 )
  {
    options.methodtext.push_back(methsec);
  }

  //cout << "reads done" << endl;;

  vector <string > spinsec;
  pos=0;
  if(readsection(words, pos, spinsec, "NSPIN"))
  {
    error("Please put NSPIN in the SYSTEM section, not in the global space.");
  }



  //Set the random seed if specified
  pos=0;
  vector <string> randnum;
  if(readsection(words, pos, randnum, "RANDOMSEED")) {
    if(randnum.size() != 2) {
      error("RANDOMSEED needs two numbers");
    }
    long int is1=atoi(randnum[0].c_str());
    long int is2=atoi(randnum[1].c_str());
    rng.seed(is1, is2);
  }
  else {
    rng.seed(12345, 98234 );
  }

  pos=0;
  readvalue(words, pos, options.runid, "RUNID");

  //cout << "doinput done\n";

}
//--------------------------------------------------------------------------
