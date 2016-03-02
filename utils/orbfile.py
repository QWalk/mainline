import numpy as np

class Orbfile:
  ###################
  def __init__(self,f):
    """Only works for a 'flat' .orb file, in which all coefficients, all basis functions, 
    all centers are listed. Needs extension otherwise."""
    in_coeff=False
    self.mo=[]
    self.center=[]
    self.basis=[]
    self.index=[]
    coeffs_linear=[]
    for line in f:
      #line=f.readline()
      if in_coeff:
        coeffs_linear.extend(map(float,line.split()))
      if "COEFFICIENTS" in line:
        in_coeff=True
      elif not in_coeff:
        spl=line.split()
        self.mo.append(int(spl[0])-1)
        self.basis.append(int(spl[1])-1)
        self.center.append(int(spl[2])-1)
        self.index.append(int(spl[3])-1)

    self.nmo=np.max(self.mo)+1
    ncoeff=len(coeffs_linear)
    print("Nmo",self.nmo, "ncoeff",len(coeffs_linear), "nfunc",len(coeffs_linear)/self.nmo)
    if(ncoeff%self.nmo==0):
      self.nfunc=int(ncoeff/self.nmo)
    else:
      print("Can't read this file because it is kind of non-standard.")
      
    #Should be fixed if we support general .orb files
    self.coeff_mat=np.array(coeffs_linear)
    self.coeff_mat=self.coeff_mat.reshape((self.nmo,self.nfunc))
##########################
  def rotate(self,orb_group,rotation):
    tmp=np.copy(self.coeff_mat)
    for ii,i in enumerate(orb_group):
      for ji,j in enumerate(orb_group):
        tmp[i,j]=0
        for ki,k in enumerate(orb_group):
          tmp[i,j]+=rotation[ii,ki]*self.coeff_mat[k,j]
    self.coeff_mat=np.copy(tmp)
    pass

#########################    
  def write(self,f):
    for m,c,b,i in zip(self.mo,self.center,self.basis,self.index):
      f.write("%i %i %i %i\n"%(m+1,b+1,c+1,i+1))
    f.write("COEFFICIENTS\n")
    count=0
    for i in self.coeff_mat.flatten():
      f.write(" %.15f "%i)
      count+=1
      if count%5==0:
        f.write("\n")
  ######################


