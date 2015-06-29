#!/usr/bin/python
from numpy  import array,dot
from mython import lines2str

# Reads a CRYSTAL input file. Can easily be modified to get more information.
def read_cryinp(filename):
  inpstr = ''
  res = {}
  titleline = True
  try:
    with open(filename,'r') as F:
      for line in F:
        if titleline: titleline=False; continue
        inpstr += line
  except IOError:
    return None
  lines = inpstr.split('\n')
  pos = 2
  # Read in geometry section which is primarily position-based.
  res['group']    = int(lines[pos]);                pos += 1
  res['latparms'] = map(float,lines[pos].split());  pos += 1
  res['natoms']   = int(lines[pos]);                pos += 1
  res['atypes'], res['apos'] = [],[]
  for i in range(res['natoms']):
    line = lines[pos].split()
    res['atypes'].append(int(line[0]))
    res['apos'].append(map(float,line[1:]))
    pos += 1

  # Rest of input not position based, and may not be separated correctly by
  # newlines. TODO: this might be true of geometry as well.

  inpl = inpstr.split()

  # Now read in important keyword arguements.
  res['calculation'] = 'hf'
  res['effchg'] = []
  while pos < len(inpl):

    if 'DFT' == inpl[pos]:
      res['calculation'] = 'dft'
      pos += 1
      continue

    if 'CORRELAT' == inpl[pos]:
      res['correlation'] = inpl[pos+1].lower()
      pos += 2
      continue

    if 'EXCHANGE' == inpl[pos]:
      res['exchange'] = inpl[pos+1].lower()
      pos += 2
      continue

    if 'HYBRID' == inpl[pos]:
      res['mixing'] = float(inpl[pos+1])
      pos += 2
      continue

    if 'PBE0' == inpl[pos]:
      res['correlation'] = 'pbe'
      res['exchange'] = 'pbe'
      res['mixing'] = 25.0
      pos += 1
      continue

    # This currently depends on the order of the entry when it doesn't have to.
    if 'INPUT' == inpl[pos]:
      res['effchg'].append(float(inpl[pos+1]))
      pos += 1
      continue

    if 'SHRINK' == inpl[pos]:
      res['kdens'] = int(inpl[pos+1].split()[0])
      pos += 2
      continue

    if 'TOLINTEG' == inpl[pos]:
      res['tolinteg'] = int(inpl[pos+1].split()[0])
      pos += 2
      continue

    if 'TOLDEE' == inpl[pos]:
      res['tole'] = int(inpl[pos+1])
      pos += 2
      continue

    if 'SPINLOCK' == inpl[pos]:
      res['spinlock'] = int(inpl[pos+1].split()[0])
      pos += 2
      continue

    if 'FMIXING' == inpl[pos]:
      res['fmixing'] = int(inpl[pos+1])
      pos += 2
      continue

    if 'BROYDEN' == inpl[pos]:
      res['broyden'] = [float(inpl[pos+1])] + map(int,inpl[pos+2:pos+4])
      pos += 5
      continue

    if 'ATOMSPIN' == inpl[pos]:
      nspin = int(inpl[pos+1])
      spinl = inpl[pos+2:pos+2+2*nspin]
      res['atomspin'] = [(int(spinl[i]),int(spinl[i+1])) for i in range(nspin)]

    pos += 1
  return res

def gen_properties(cryinp,natoms,kpath,denom,projs,
                   get_ef=1,nprint=0,title='blank',
                   above=5,below=8,
                   kresfactor=20,eres=200,
                   npoly=25,):
  """ User-friendly properties input generator.

      Factor "denom" out of kpath to make it integer.
      TODO: negative DOSS projects (onto atoms)
  """
  inp = read_cryinp(cryinp)
  shrink = inp['kdens']
  epera  = array(inp['effchg'])
  
  # This may only make sense when the spin channels are equal.
  filled = int(round(dot(array(natoms),epera))) / 2

  newklines = [['NEWK']]
  newklines.append([shrink,shrink*2])
  newklines.append([get_ef,nprint])

  nline       = len(kpath)
  kres        = kresfactor*nline
  highest     = filled + above
  lowest      = filled - below
  outtofile   = 1
  printevals  = 0
  bandlines = [['BAND']]
  bandlines.append([title])
  bandlines.append([nline,denom,kres,lowest,highest,outtofile,printevals])
  for si in range(len(kpath)-1):
    bandlines.append(kpath[si]+[' ']+kpath[si+1])
  bandlines.append(kpath[-1]+[' ']+kpath[0])

  nproj      =  len(projs)
  printopts  = 0
  doslines = [['DOSS']]
  doslines.append([nproj,eres,lowest,highest,outtofile,npoly,printopts])
  for pi in range(len(projs)):
    doslines.append([len(projs[pi]),' ']+projs[pi])

  outlines = newklines + bandlines + newklines + doslines + [["END"]]
  return lines2str(outlines)


