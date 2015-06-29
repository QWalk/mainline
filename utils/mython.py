from numpy import ndarray,frombuffer
from base64 import b64encode,b64decode
from json import JSONEncoder
from os import getcwd
from datetime import datetime

# Lines should be a list of lists of words.
# Words are separated by spaces, lines by \n.
def lines2str(lines):
  outlines = []
  for line in lines:
    outlines.append(' '.join(map(str,line)))
  outstr = '\n'.join(outlines)
  return outstr

class Ldict(dict):
  """
  Dictionary that can append data rather than overwriting it when append item is
  used more than once.
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

# Source (simplified by me):
# http://stackoverflow.com/questions/3488934/simplejson-and-numpy-array/24375113#24375113
class NumpyEncoder(JSONEncoder):
  def default(self, obj):
    if isinstance(obj, ndarray):
      print obj
      return dict(__ndarray__=obj.tolist(),
                  dtype=str(obj.dtype),
                  shape=obj.shape)
    # Let the base class default method raise the TypeError
    return JSONEncoder(self, obj)
class NumpyToListEncoder(JSONEncoder):
  def default(self, obj):
    if isinstance(obj, ndarray):
      return obj.tolist()
    # Let the base class default method raise the TypeError
    return JSONEncoder(self, obj)
def json_numpy_obj_hook(dct):
  """ Hook for facier numpy encoder

  Example usage:
  expected = np.arange(100, dtype=np.float)
  dumped = json.dumps(expected, cls=NumpyEncoder)
  result = json.loads(dumped, object_hook=json_numpy_obj_hook)
  """

  if isinstance(dct, dict) and '__ndarray__' in dct:
    data = ndarray(dct['__ndarray__'],dtype=dct['dtype'])
    return data.reshape(dct['shape'])
  return dct

def gen_qsub(exe,stdout='',loc='',name='',time='72:00:00',nn=1,np=1,
    prep_commands=[],final_commands=[]):
  """ Generate a qsub file.
  
  Blank strings will generate useful defaults."""

  if stdout=='': stdout='qsub.out'
  if loc=='': loc=getcwd()
  if name=='': name=str(datetime.now()).replace(' ','_')
  header = []
  header.append('#!/bin/bash')
  header.append('#PBS -l nodes=%d:ppn=%d'%(nn,np))
  header.append('#PBS -l walltime=%s'%time)
  header.append('#PBS -j oe')
  header.append('#PBS -m n')
  header.append('#PBS -N %s'%name)
  header.append('#PBS -o qsub.out')
  mpiexe = 'mpirun -n %d %s &> %s'%(nn*np, exe, stdout)
  commands = header + ['cd %s'%loc] + prep_commands + [mpiexe] + final_commands
  outstr = '\n'.join(commands)
  with open(loc+'/qsub.in','w') as qsin:
    qsin.write(outstr)
