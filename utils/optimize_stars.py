#!/usr/bin/env python
import sys
import os
import scipy.optimize
import re
from optparse import OptionParser

HELP="""This program will optimize the energy with respect to any numerical parameter
which is preceded by a key character (typically a *, but this can be changed by 
changing the 'key' variable in the code).

The final optimized file will be in 'staropt.inp'

Supported DFT/quantum chemistry programs are listed below.  As part of the 
installation, you may need to change the executables to the correct location for
your computer.

usage: %s [options] inputfile
"""%sys.argv[0]
key="*"
GAMESS_EXE="gms"
CRYSTAL_EXE="crystal"

def gen_parm_list(lines):
  """Search lines for numbers with *'s next to them and extract
  them into a list of parameters"""
  parms=[]
  for line in lines:
    spl=line.split()
    for word in spl:
      if word[0]==key:
        parms.append(float(word[1:]))

  return parms

def apply_parm_to_list(parms,lines):
  count=0
  newlines=[]
  for line in lines:
    if line.count(key) == 0:
      newlines.append(line)
    else:
      spl=line.split()
      newline=""
      for i,word in enumerate(spl):
        if word[0]==key:
          newline+=str(parms[count])+" "
          count+=1
        else:
          newline+=word+ " "
      newlines.append(newline+"\n")
  return newlines

        
def get_en_gamess(parms,lines):
  nwlines=apply_parm_to_list(parms,lines)
  f=open("staropt.inp",'w')
  for line in nwlines:
    f.write(line)
  f.close()
  os.system(GAMESS_EXE+" staropt.inp")
  f=open("staropt.log",'r')
  en=0.0
  for l in f.readlines():
    if l.count("TOTAL ENERGY =")>0:
      spl=l.split()
      en=float(spl[3])
  print "en ",en, parms
  sys.stdout.flush()
  return en
  

def get_en_crys(parms,lines):
  nwlines=apply_parm_to_list(parms,lines)
  f=open("staropt.inp",'w')
  for line in nwlines:
    f.write(line)
  f.close()
  os.system(CRYSTAL_EXE+" < staropt.inp > staropt.inp.o")
  f=open("staropt.log",'r')
  en=0.0
  for l in f.readlines():
    if l.count("SCF ENDED")>0:
      spl=l.split()
      en=float(spl[8])
  print "en ",en, parms
  sys.stdout.flush()
  return en

if __name__ == "__main__":
  parser = OptionParser(usage=HELP)
  parser.add_option("-g", "--gamess",action="store_true", default=False,
                  help="use GAMESS to evaluate file.")
  parser.add_option("-c", "--crystal",action="store_true", default=False,
                  help="use CRYSTAL to evaluate file.")

  (options, args) = parser.parse_args()
  if len(args) < 1:
    parser.print_usage()
    quit()

  f=open(args[0],'r')
  lines=f.readlines()
  f.close()
  parms=gen_parm_list(lines)
  if options['gamess']:
    parmopt=scipy.optimize.fmin_powell(get_en_gamess,parms,args=[lines])
    get_en_gamess(parmopt,lines)
  elif options['crystal']:
    parmopt=scipy.optimize.fmin_powell(get_en_crystal,parms,args=[lines])
    get_en_crystal(parmopt,lines)
  else:
    parser.print_usage()
    quit()



