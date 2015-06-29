#!/usr/bin/python
from numpy            import array
from mython           import Ldict
from subprocess       import call
from dm_tools         import read_dm
from read_numberfluct import read_number_dens,moments
from cryfile_io       import read_cryinp
import sys

# Reads a qwalk input section.
def read_section(inp,key,pos):
  res = Ldict()
  while inp[pos] != '}':
    if isinstance(inp[pos],float):           # If it's a number.
      while isinstance(inp[pos],float):
        #print 'Found data',key,inp[pos]
        res.append(key,inp[pos])
        pos += 1
    elif inp[pos] == '{':                    # Else if it's demarking a section,
      if isinstance(inp[pos+1],str):         # which could have keywords,
        label = inp[pos+1].lower()
        pos += 2
        #print 'Reading section',key,label
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
        #print 'Found data',key,val
        res[key] = val
        pos += 1
    else:                                    # Else it's a keyword.
      key = inp[pos].lower()
      if key not in res.keys():
        #print 'Setting',key
        res[key] = True
      pos += 1
  pos += 1
  return res,pos

# Reads a qwalk input file.
def read_qfile(filename):
  inpstr = ''
  try:
    with open(filename,'r') as F:
      for line in F:
        if '#' in line: # TODO: Needs to be fixed when '#' isn't the first thing in line.
          print 'Warning, reading commented lines is incomplete!'
          continue
        inpstr += line
  except IOError:
    return None
  # Ensure correct splitting. This is inefficient for large files.
  inpstr = inpstr.replace('{',' { ')
  inpstr = inpstr.replace('}',' } ')
  inp    = inpstr.split() + ['}'] # Now I can use "read_section" on the file!
  for i in range(len(inp)): # Make everything numeric into floats.
    try:                inp[i] = float(inp[i])
    except ValueError:  pass
  return read_section(inp,filename,0)[0]

# Temporary function to convert file names to metadata about the calculations.
# In the future, an optional metadata file should contain overall qualitiative
# descriptions of the data, like magnetic ordering, or pressure. These things
# are redundant, but convenient for interpreting sets of numbers together.
def convert_to_metadata(froot):
  mconv = {'che':'checkerboard',
           'str':'collinear',
           'bic':'bicollinear',
           'fst':'collinear, flip 1',
           'dim':'collinar, flip 2',
           'sta':'staggered'}
  basename = froot.split('/')[-1]
  mag = mconv[basename[:3]]
  prs = float(basename[3:].translate(None,'_prs'))
  return {'mag':mag,'prs':prs}

# Use golsing to read out average energy and error.
def read_qenergy(logfile,gosling='./gosling'):
  statfile = logfile.replace('.log','.stat')
  with open(statfile,'w') as out:
    call([gosling, logfile], stdout = out)
  with open(statfile,'r') as F:
    for line in F:
      if 'total_energy0' in line:
        spl = line.split()
        return {'egy':spl[1],'err':spl[3],'var':spl[5]}
  print 'ERROR: cannot find total_energy0 in stat file.'
  return {'egy':None,'err':None,'var':None}

# Read a directory which has multiple files with data into a single dictionary
# with relevant information.
def read_dir_forlucas(froot,gosling='./gosling'):
  # This first section should be edited to reflect naming conventions!
  dftfile   = froot+'.d12'
  bres = {} # Data that is common to all k-points.
  ress = [] # List of all k-point data in directory.

  metad = convert_to_metadata(froot)
  bres['authors']  = 'BL'
  bres['ordering'] = metad['mag']

  print dftfile
  dftdat = read_cryinp(dftfile)
  bres['a'] = dftdat['latparms'][0]
  bres['c'] = dftdat['latparms'][1]
  bres['se_height'] = dftdat['apos'][dftdat['atypes'].index(234)][-1]
  bres['hybrid_mixing'] = dftdat['mixing']

  # Start by taking only the real k-points, since I'm sure these are on solid
  # ground, and have enough sample points.
  realk = array([1,3,8,10,27,29,34,36]) - 1
  for rk in realk:
    kroot = froot + '_' + str(rk)

    print "Working on",kroot+"..." 
    sysdat = read_qfile(kroot+'.sys')
    if sysdat==None:
      print "  (cannot find QMC output, skipping)"
      continue
    else:
      ress.append(bres.copy())
      ress[-1]['kpoint'] = sysdat['system']['kpoint']
    
    print "  energies..." 
    egydat = read_qenergy(kroot+'.dmc.log','/home/brian/bin/gosling')
    if egydat==None:
      print "  (cannot find ground state energy)"
    else:
      ress[-1]['total_energy']     =  egydat['egy']
      ress[-1]['total_energy_err'] =  egydat['err']

    ogpdat = read_qenergy(kroot+'.ogp.log','/home/brian/bin/gosling')
    if egydat==None:
      print "  (cannot find excited state energy)"
    else:
      ress[-1]['excited_energy']     =  ogpdat['egy']
      ress[-1]['excited_energy_err'] =  ogpdat['err']

    print "  fluctuations..." 
    fludat, fluerr = read_number_dens(kroot+'.ppr.o')
    if fludat==None:
      print "  (cannot find number fluctuation)"
    else:
      avg, var, cov, avge, vare, cove = moments(fludat,fluerr)
      ress[-1]['covariance']          =  cov
      ress[-1]['site_charge']         =  [] # TODO

    print "  1-RDM..." 
    odmdat = read_dm(kroot+'.ordm.o')
    if odmdat==None:
      print "(cannot find 1-RDM)"
    else:
      ress[-1]['1rdm'] = odmdat

    print "  done."; 

  return ress

# Currently works on one-body density matrix, can be extended (TODO).
def read_dm(dmfile):
  obdm_section = False
  linecount = 0
  try:
    f=open(dmfile,'r')
  except IOError:
    return None

  for line in f:
    if 'tbdm: nmo' in line:
      nmo = int(line.split()[-1])
    if obdm_section:
      orbdm_data += line
      linecount  += 1
      if linecount > nmo**2:
        obdm_section = False
    else:
      if 'One-body density matrix' in line:
        obdm_section = True
    
  f.close()
