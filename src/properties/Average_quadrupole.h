#include "Average_generator.h"

#ifndef AVERAGE_QUADRUPOLE_H_INCLUDED
#define AVERAGE_QUADRUPOLE_H_INCLUDED
#include "Array45.h"
/*!
*/
class Average_quadrupole:public Average_generator { 
public:
  virtual void evaluate(Wavefunction_data * wfdata, Wavefunction * wf,
                        System * sys, Sample_point * sample, Average_return & avg);
  virtual void read(System * sys, Wavefunction_data * wfdata, vector
                    <string> & words);
  virtual void write_init(string & indent, ostream & os);
  virtual void read(vector <string> & words);
  virtual void write_summary(Average_return &,Average_return &, ostream & os);
  virtual void jsonOutput(Average_return &,Average_return &, ostream & );
  
private:
  Array2 <doublevar> gvecs; 
};


class Average_joint_pbc:public Average_generator { 
public:
  virtual void evaluate(Wavefunction_data * wfdata, Wavefunction * wf,
                        System * sys, Sample_point * sample, Average_return & avg);
  virtual void read(System * sys, Wavefunction_data * wfdata, vector
                    <string> & words);
  virtual void write_init(string & indent, ostream & os);
  virtual void read(vector <string> & words);
  virtual void write_summary(Average_return &,Average_return &, ostream & os);
  virtual void jsonOutput(Average_return &,Average_return &, ostream & );
  
private:
  Array2 <doublevar> gvecs; 
};



#endif //AVERAGE_QUADRUPOLE_H_INCLUDED

