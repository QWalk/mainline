from __future__       import print
from subprocess       import call
from os               import getcwd
from mython           import gen_qsub
import json

class Ldict(dict):
  """
  Dictionary that can append data rather than overwriting it when key is
  assigned more than once.
  """
  def __init__(self,**kwargs):
    dict.__init__(self,kwargs)

  def append(self,key,value):
    if key in dict.keys(self):
      if isinstance(dict.__getitem__(self,key),list):
        dict.__getitem__(self,key).append(value)
      elif isinstance(dict.__getitem__(self,key),bool):
        dict.__setitem__(self,key,value)
      else:
        dict.__setitem__(self,key,[dict.__getitem__(self,key),value])
    else:
      dict.__setitem__(self,key,value)

# Reads a qwalk input section.
def read_section(inp,key,pos):
  res = Ldict()
  while inp[pos] != '}':
    if isinstance(inp[pos],float):           # If it's a number.
      while isinstance(inp[pos],float):
        res.append(key,inp[pos])
        pos += 1
    elif inp[pos] == '{':                    # Else if it's demarking a section,
      if isinstance(inp[pos+1],str):         # which could have keywords,
        label = inp[pos+1].lower()
        pos += 2
        val,pos = read_section(inp,key,pos)
        if label != False:
          val['label'] = label
        res.append(key,val)
      else:                                  # or just numbers.
        pos += 1
        val = []
        while isinstance(inp[pos],float):
          val.append(inp[pos])
          pos += 1
        if len(val) == 1: val = val[0]
        res[key] = val
        pos += 1
    else:                                    # Else it's a keyword.
      key = inp[pos].lower()
      if key not in res.keys():
        res[key] = True
      pos += 1
  pos += 1
  return res,pos

# Reads a qwalk input file.
def read_qfile(inpf):
  inpstr = ''
  for line in inpf:
    if '#' in line: # TODO: Needs to be fixed when '#' isn't the first thing in line.
      print("qfiles_io.py: Warning, reading commented lines is incomplete!")
      continue
    inpstr += line
  # Ensure correct splitting. This is inefficient for large files.
  inpstr = inpstr.replace('{',' { ')
  inpstr = inpstr.replace('}',' } ')
  inp    = inpstr.split() + ['}'] # Now I can use "read_section" on the file!
  for i in range(len(inp)): # Make everything numeric into floats.
    try:                inp[i] = float(inp[i])
    except ValueError:  pass
  return read_section(inp,inpf.name,0)[0]
