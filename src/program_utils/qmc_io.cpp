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

#include "qmc_io.h"

//---------------------------------------------------------------------

int caseless_eq(const string & s1, const string & s2) {
  if(s1.size()!=s2.size()) return 0;

  string::const_iterator i1=s1.begin();
  string::const_iterator i2=s2.begin();
  //cout << s1 << "  "    << s2 << endl;
  while( (i1!=s1.end()) && (i2 != s2.end()) ) {
    //cout << *i1 << "  " << *i2 << endl;
    if( toupper(*i1) != toupper(*i2))
      return 0;
    i1++; i2++;
  }

  return 1;
}


void canonical_filename(string & str)
{
#ifdef USE_MPI
  char strbuff[40];
  sprintf(strbuff, "%d", mpi_info.node);
  str+="_";
  str+=strbuff;
#endif
}
void canonical_filename(string & str, int num)
{
#ifdef USE_MPI
  char strbuff[40];
  sprintf(strbuff, "%d", num);
  str+="_";
  str+=strbuff;
#endif
}

void append_number(string & str, int num)
{
  char strbuff[40];
  sprintf(strbuff, "%d", num);
  str+=strbuff;
}

void append_number(string & str, double num) { 
  char strbuff[180];
  sprintf(strbuff, "%g.15", num);
  str+=strbuff;
}



int haskeyword(vector <string> & input, unsigned int & pos,
               const char * valuename)
{
  int level=0;  //how many brackets we've descended into.

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

    if(level==0 && caseless_eq(input[pos], valuename))//input[pos]==valuename)
    {
      return 1;
    }
  }
  //if we get here, we didn't find a value
  return 0;
}


void readnext(vector <string> & s,unsigned int & i, int & t)
{
  if( i >= s.size() )
  {
    error("Unexpected end of file at ",s[i]);
  }
  //cout << "readnext:\n";
  //cout << s[i+1] << endl;
  t = atoi(s[++i].c_str());
  //cout << t << endl;
}

void readnext(vector <string> & s,unsigned int & i, doublevar & t)
{
  if( i >= s.size() )
  {
    error("Unexpected end of file at ",s[i]);
  }
  //cout << "readnext:\n";
  //cout << s[i+1] << endl;
  t = atof(s[++i].c_str());
  //cout << t << endl;
}

void readnext(vector <string> & s,unsigned int & i, string & t)
{
  if( i+1 >= s.size() )
  {

    error("Unexpected end of section at ",s[i]);
  }
  //cout << "readnext:\n";
  //cout << s[i+1] << endl;
  t = s[++i].c_str();
  //cout << t << endl;
}

//----------------------------------------------------------------------
/*!
Puts the next section after pos with name sectionname
into section.  Stays within the same level that we started;
if the section is ended before finding anything, we return
0.
*/
int readsection(vector <string> & input,
                unsigned int & pos,
                vector <string> & section,
                const char * sectionname)
{
  string temp;
  section.erase(section.begin(), section.end());
  int level=0;
  string name(sectionname);
  //cout << "searching for " << name << endl;
  //for(int i=0; i< input.size(); i++) 
  //  cout << input[i] << endl;
  //cout << "____________________________________" << endl;

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
    //if(input[pos]==sectionname && level == 0)
    //if(level==0) cout << "looking at " << input[pos] 
    //                  << " : " << caseless_eq(input[pos], name)
    //                  << endl;
    if(level==0 && caseless_eq(input[pos], name))
    {
      int nsec=0;
      pos++;  //check for a section start
      if(pos >= input.size()) 
        error("Section starter ", 
              sectionname, " but nothing following it.");
      if(input[pos] != startsec) 
        error("You must follow ", sectionname, 
              "with a ", startsec);
      pos++;
      int begin=pos;
      for(; pos < input.size(); pos++) {
        if(input[pos]==startsec) nsec++;
        if(input[pos]==endsec) nsec--;
        if(nsec < 0) { break;}
        //section.push_back(input[pos]);
      }
      pos++; //step past the end
      int end=pos;
      section.insert(section.end(), 
                     input.begin()+begin, 
                     input.begin()+end-1);
      //input.erase(input.begin()+begin-2, input.begin()+end);
      //for(int i=0; i< section.size(); i++) 
      //  cout << section[i] << "   ";
      //cout << endl;
      //cout << "---------------input left " << endl;
     //for(int i=0; i< input.size(); i++) 
      //  cout << input[i] << "   ";
      //cout << endl;


      /*
      readnext(input, pos, temp);
      if(temp != startsec)
      {
        error("You must follow ", sectionname, "with a ", startsec);
      }
      while(pos < input.size() )
      {
        readnext(input, pos, temp);

        if(temp == startsec)
        {
          nsec++;
        }

        if(temp == endsec)
        {
          nsec--;
        }

        if(nsec<0)
        {
          pos++;
          break;
        }
        section.insert(section.end(), temp);
      
      }
      */
      return 1;
    }
  }
  //if we get here, we didn't find a section.
  return 0;
}

//----------------------------------------------------------------------

/*!
\bug
When you have a section completely in one  block in characters
(ie {1}), it does not realize that the start and end are there, 
and chaos ensues. 
 
 */
void parsefile(ifstream & inputfile, vector <string> & words)
{
  string temp;
  vector <string> tempwords;
  while(inputfile >> temp)
  {
    string temp2, temp3;
    unsigned int pos;

    //Ignore comments
    pos=temp.find(comment);
    if( pos < temp.size() )
    {
      temp.erase(pos,temp.size()-pos);
      inputfile.ignore(800, '\n');
    }

    //Put start section markers in their own word
    pos = temp.find(startsec);
    if( pos < temp.size() )
    {
      temp.erase(pos,startsec.size());
      temp2=startsec;
      temp3.append(temp, pos, temp.size());
      temp.erase(pos, temp.size()-pos);
    }

    //Put any end section markers in their
    //own word
    pos = temp.find(endsec);
    if( pos < temp.size() )
    {
      temp.erase(pos,endsec.size());
      temp2=endsec;
      temp3.append(temp, pos, temp.size());
      temp.erase(pos, temp.size()-pos);
    }

    if( ! temp.empty() )
      tempwords.push_back(temp);
    if( ! temp2.empty() )
      tempwords.push_back(temp2);
    if( ! temp3.empty() )
      tempwords.push_back(temp3);
  }

  // Include any include files(It's recursive here..)
  for(unsigned int i=0; i < tempwords.size(); i++)
  {
    if(caseless_eq(tempwords[i],"INCLUDE"))
    {
      ifstream includefile(tempwords[i+1].c_str());
      vector <string> includewords;
      if(!includefile)
      {
        error("couldn't open ", tempwords[i+1]);
      }
      parsefile(includefile, includewords);
      for(unsigned int j=0; j< includewords.size(); j++)
      {
        words.push_back(includewords[j]);
      }
      includefile.close();
      i++;
    }
    else
    {
      words.push_back(tempwords[i]);
    }
  }


}


//----------------------------------------------------------------------



int checkbrackets(vector <string> & words)
{
  int size=words.size();
  int nopen=0;
  int nclose=0;
  for(int i=0; i< size; i++)
  {
    if(words[i]==startsec)
    {
      nopen++;
    }
    else if(words[i]==endsec)
    {
      nclose++;
    }
  }
  if(nopen!=nclose) {
    single_write(cout, "you have unbalance brackets: ", nopen, " open ");
    single_write(cout, " and ", nclose, " close\n");
    error("The open and close brackets don't match.");
  }
  return 1;
}

//----------------------------------------------------------------------
