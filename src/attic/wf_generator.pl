#!/usr/bin/perl

print "This is outdated.  If you update it, please remove this line.\n";
exit(1);


if(scalar @ARGV < 1) {
    print "usage:  wf_generator <object name>\n";
    exit(1);
}
$wf_name=$ARGV[0];
$wf_caps=$wf_name;
$wf_caps =~ tr/a-z/A-Z/;

$headerfile="
//--------------------------------------------------------------------------
// $wf_name.h
//  

#ifndef $wf_caps\_H_INCLUDED
#define $wf_caps\_H_INCLUDED

#include \"Wavefunction.h\"
#include \"$wf_name\_data.h\"
#include \"Toolkit.h\"

/*
Don't forget to put allocation information in Wavefunction.cpp and
Wavefunction_data.cpp, and to add the .cpp files in the makefile.
*/
class $wf_name : public Wavefunction {
public:
       virtual void read(vector <string> & words,
                      unsigned int & pos,
                      Program_options & options
                     );

    /*!
    Link the wavefunction with a data object.  From this point
    forward, the wavefunction will use the parameters in the Wavefunction_data
    object to do the calculations.
    */
    virtual void link(Wavefunction_data * wfdata ) {
         recast(wfdata, dataptr);
    }



    virtual void calcVal(System & , Sample_point &);
    virtual void updateVal(System &, Sample_point &, int, doublevar &);
    virtual void rejectVal(int );
    virtual void getVal(int, doublevar &);

    virtual void calcGrad(System &, Sample_point &);
    virtual void updateGrad(System &, Sample_point &, int, Array1 <doublevar> &);
    virtual void rejectGrad(int);
    virtual void getGrad(int, Array1 <doublevar> &);

    virtual void calcLap(System &, Sample_point &);
    virtual void updateLap(System &, Sample_point &, int, Array1 <doublevar> &);
    virtual void rejectLap(int);
    virtual void getLap(int, Array1 <doublevar> &);

    virtual void calcForceBias(System &, Sample_point &);
    virtual void updateForceBias(System &, Sample_point &, int,  Array1 <doublevar> &);
    virtual void rejectForceBias(int);
    virtual void getForceBias(int, Array1 <doublevar> &);

private:
    $wf_name\_data * dataptr;
    //Put the information that needs to be stored for each walker here.

};

#endif //$wf_caps\_H_INCLUDED
//--------------------------------------------------------------------------
";

$dataheaderfile="
//------------------------------------------------------------------------
//$wf_name\_data.h

#ifndef $wf_caps\_DATA_H_INCLUDED
#define $wf_caps\_DATA_H_INCLUDED

#include \"Qmc_std.h\"
class Program_options;
#include \"Wavefunction_data.h\"

class $wf_name\_data : public Wavefunction_data {
public:

   void read(vector <string> & words,
                      unsigned int & pos,
                      Program_options & options
                     );
   void getVarParms(Array1 <doublevar> & parms);
   void setVarParms(Array1 <doublevar> & parms);
   int showinfo(ostream & os);
   int nparms();

private:
   friend class $wf_name;
   //put the static data elements here, that are common across all
   //the walkers.

};

#endif //$wf_caps\_DATA_H_INCLUDED
//------------------------------------------------------------------------
";

$sourcefile="
//--------------------------------------------------------------------------
// $wf_name.cpp
//  
//
#include \"$wf_name.h\"
#include \"qmc_io.h\"
#include \"Qmc_std.h\"
#include \"Sample_point.h\"

/*
Size any arrays from the input section.
*/
void $wf_name\::read(vector <string> & words, unsigned int & pos,
                   Program_options & options) {
}


/*
Calculation methods.  The only absolutely necessary one is the Lap methods;
all the others are by default hardwired into that one.  However, you should
at the very least write a value-only calculation, which drastically 
improves the efficiency of pseudopotential calculations.
*/
//------------------------------------------------------------------------

void $wf_name\::calcVal(System & system, Sample_point & sample) {
    doublevar val;
    for(int e=0; e< sample.size(); e++) {
        updateVal(system, sample, e, val);
    }
}

//------------------------------------------------------------------------

void $wf_name\::updateVal(System & system, Sample_point & sample,int e, doublevar & val) {
    assert(dataptr != NULL);
}

//------------------------------------------------------------------------

void $wf_name\::rejectVal(int e) {
}

//------------------------------------------------------------------------

void $wf_name\::getVal(int e, doublevar & val) {
}

//----------------------------------------------------------------------------

void $wf_name\::calcGrad(System & system, Sample_point & sample) {
    Array1 <doublevar> grad(4);
    for(int e=0; e< sample.size(); e++) {
        updateGrad(system, sample, e, grad);
    }
}

void $wf_name\::updateGrad(System & system, Sample_point & sample,int e,  Array1 <doublevar> & newgrad){
    assert(newgrad.GetDim(0) >= 4);
    Array1 <doublevar> temp(5);
    updateLap(system, sample, e, temp);
    for(int i=0; i< 4; i++) {
        newgrad(i)=temp(i);
    }
}
void $wf_name\::rejectGrad(int e){
   rejectLap(e);
}
void $wf_name\::getGrad(int e, Array1 <doublevar> & grad) {
    assert(grad.GetDim(0) >= 4);
    Array1 <doublevar> temp(5);
    getLap(e, temp);
    for(int i=0; i< 4; i++) {
        grad(i)=temp(i);
    }   
}

void $wf_name\::calcForceBias(System & system, Sample_point & sample) {
    Array1 <doublevar> bias(4);
    for(int e=0; e< sample.size(); e++) {
        updateForceBias(system, sample, e, bias);
    }
}
void $wf_name\::updateForceBias(System & system, Sample_point & sample,int e, Array1 <doublevar> & newbias) {
    updateGrad(system, sample, e, newbias);
}
void $wf_name\::rejectForceBias(int e) {
    rejectGrad(e);
}
void $wf_name\::getForceBias(int e, Array1 <doublevar> & bias) {
    getGrad(e, bias);
}

//------------------------------------------------------------------------

void $wf_name\::calcLap(System & system, Sample_point & sample) {
    Array1 <doublevar> lap(5);
    for(int e=0; e< sample.size(); e++) {
        updateLap(system, sample, e, lap);
    }
}

//------------------------------------------------------------------------

void $wf_name\::rejectLap(int e) {
    assert(dataptr != NULL);
}

//-------------------------------------------------------------------------

void $wf_name\::getLap(int e, Array1 <doublevar> & lap) {
    assert(lap.GetDim(0) >= 5);
    assert(dataptr != NULL);
}

//-------------------------------------------------------------------------


void $wf_name\::updateLap(
    System & system,
    Sample_point & sample,
    int e,
    //!<electron number
    Array1 <doublevar> & newlap
) {
    assert( newlap.GetDim(0) >= 5);
    assert(dataptr != NULL);
}

//-------------------------------------------------------------------------
";

$datasourcefile="
//------------------------------------------------------------------------
//$wf_name\_data.cpp

#include \"Qmc_std.h\"
#include \"qmc_io.h\"
#include \"$wf_name\_data.h\"
#include \"Program_options.h\"

/*!
*/
void $wf_name\_data::read(vector <string> & words, unsigned int & pos,
                   Program_options & options)
{
}

void $wf_name\_data::getVarParms(Array1 <doublevar> & parms) {
}

void $wf_name\_data::setVarParms(Array1 <doublevar> & parms) {
}

int $wf_name\_data::nparms() {
}

int $wf_name\_data::showinfo(ostream & os) {
}

//------------------------------------------------------------------------
";

open(OUTPUT, ">$wf_name.h");
print OUTPUT $headerfile;
close(OUTPUT);

open(OUTPUT, ">$wf_name\_data.h");
print OUTPUT $dataheaderfile;
close(OUTPUT);

open(OUTPUT, ">$wf_name.cpp");
print OUTPUT $sourcefile;
close(OUTPUT);

open(OUTPUT, ">$wf_name\_data.cpp");
print OUTPUT $datasourcefile;
close(OUTPUT);
