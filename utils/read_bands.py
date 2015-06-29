import numpy as np

def read_band_section(inp,line):
  # Read in bands segment.
  spl = line.split()
  nband = int(spl[1])
  nkp   = int(spl[2])
  dkp   = float(line[-24:-13])
  efermi = float(line[-13:])

  line = inp.readline()
  emin = float(line[:12])
  emax = float(line[12:])

  line = inp.readline().split()
  k0 = ''.join(line[:3])
  k1 = ''.join(line[3:])

  dat = []
  for v in range(nband*nkp):
    dat.append(float(inp.read(12)))
    if (v+1)%6 == 0: inp.read(1)
  dat = np.array(dat).reshape(nkp,nband)
  return [dat,k0,k1,dkp,efermi]

def read_dos_section(inp,line):
  spl = line.split()
  negy  = int(spl[2])
  degy  = float(line[-24:-13])
  efermi = float(line[-13:])
  #print line,
  #print nband, negy, degy, efermi

  line = inp.readline()
  estart = float(line[12:])
  #print line,estart

  line = inp.readline().split()
  pi = int(line[0])
  natom = line[1] 
  #print line, pi, oi

  dat = []
  for v in range(negy):
    dat.append(float(inp.read(12)))
    if (v+1)%6 == 0: inp.read(1)
  return [np.array(dat),pi,estart,degy,efermi]

def read_fort25(fn):
  """ Takes fort.25 filename, returns data contained within """
  inp = open(fn,'r')

  alldat = {}
  alldat['bands'] = []
  alldat['dos'] = []
  
  line = 'start'
  while line != '':
    line = inp.readline()
    #print line

    if 'BAND' in line: 
      keys = ['dat','k0','k1','dkp','efermi']
      alldat['bands'].append( dict(zip(keys,read_band_section(inp,line))) )
    if 'DOSS' in line: 
      keys = ['dat','pi','estart','degy','efermi']
      alldat['dos'].append( dict(zip(keys,read_dos_section(inp,line))) )
  return alldat
