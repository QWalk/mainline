# Lucas Wagner, edited by Brian Busemeyer

import numpy as np
def read_number_dens(filename):
  try:
    f=open(filename,'r')
  except IOError:
    return None,None
  data=np.zeros((0,0,0,0,0,0))
  data_err=np.zeros((0,0,0,0,0,0))
  while True:
    line=f.readline()
    #print line
    if line.find("Region fluctuation")!=-1:
      line=f.readline()
      spl=line.split()
      nspin=int(spl[1])
      maxn=int(spl[3])
      nregion=int(spl[5])
      data=np.zeros((nspin,nspin,nregion,nregion,maxn,maxn))
      #print "data ", data.shape
      data_err=np.zeros((nspin,nspin,nregion,nregion,maxn,maxn))
      
      for s1 in range(0,nspin):
        for s2 in range(0,nspin):
          for r1 in range(0,nregion):
            for r2 in range(0,nregion):
              line=f.readline()
              #print line
              for n1 in range(0,maxn):
                spl=f.readline().split()
                for n2 in range(0,maxn):
                  data[s1,s2,r1,r2,n1,n2]=float(spl[n2])
                for n2 in range(maxn,2*maxn):
                  data_err[s1,s2,r1,r2,n1,n2-maxn]=float(spl[n2])
      break
  return data,data_err

def moments(data,data_err):
  nspin=data.shape[0]
  nregion=data.shape[2]
  maxn=data.shape[4]
  
  avg=np.zeros((nspin,nregion))
  avg_err=np.zeros((nspin,nregion))
  for s in range(0,nspin):
    for r in range(0,nregion):
      for n in range(0,maxn):
        avg[s,r]+=n*data[s,s,r,r,n,n]
        avg_err[s,r]+=n**2*data_err[s,s,r,r,n,n]**2
  avg_err = avg_err**.5

  var=np.zeros((nspin,nregion))
  var_err=np.zeros((nspin,nregion))
  for s in range(0,nspin):
    for r in range(0,nregion):
      for n in range(0,maxn):
        var[s,r] += (n-avg[s,r])**2 * data[s,s,r,r,n,n]
        var_err[s,r]+=data_err[s,s,r,r,n,n]**2*(n-avg[s,r])**4 + 2*data[s,s,r,r,n,n]*avg_err[s,r]**2*(n-avg[s,r])**2
  var_err = var_err**.5

  covar=np.zeros((nspin,nspin,nregion,nregion))
  covar_err=np.zeros((nspin,nspin,nregion,nregion))
  for s1 in range(0,nspin):
    for s2 in range(0,nspin):
      for r1 in range(0,nregion):
        for r2 in range(0,nregion):
          corr=0.0
          corr_err=0.0
          for n1 in range(0,maxn):
            for n2 in range(0,maxn):
              p=data[s1,s2,r1,r2,n1,n2]
              pe=data_err[s1,s2,r1,r2,n1,n2]
              corr += p*(n2-avg[s2,r2])*(n1-avg[s1,r1])
              corr_err+= pe**2*(n1-avg[s1,r1])**2*(n2-avg[s2,r2])**2 + p**2*avg_err[s,r]**2*(n2-avg[s2,r2])**2 + p**2*avg_err[s,r]**2*(n1-avg[s1,r1])**2
          covar[s1,s2,r1,r2]=corr
          covar_err[s1,s2,r1,r2]=corr_err
  covar_err = covar_err**.5

  return avg,var,covar,avg_err,var_err,covar_err
