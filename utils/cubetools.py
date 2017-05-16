import sys
import numpy as np

def read_cube(filename):
  f=open(filename, 'r')
  cube={}
  cube['comment']=f.readline()
  cube['type']=f.readline()
  spl=f.readline().split()
  cube['natoms']=int(spl[0])
  cube['origin']=list(map(float, spl[1:]))
  cube['npoints']=np.array([0,0,0])
  cube['latvec']=np.zeros((3,3))
  for i in range(0,3):
    spl=f.readline().split()
    cube['npoints'][i]=int(spl[0])
    cube['latvec'][i,:]=list(map(float,spl[1:]))
  natoms=cube['natoms']
  cube['atomname']=[]
  cube['atomxyz']=np.zeros((natoms,3))
  for i in range(0,natoms):
    spl=f.readline().split()
    cube['atomname'].append(spl[0])
    cube['atomxyz'][i,:]=list(map(float,spl[2:]))
  cube['data']=np.zeros(cube['npoints'])
  vector=[]
  while True:
    spl=f.readline().split()
    if len(spl) < 1:
      break
    vector.extend(map(float,spl))
  nread=len(vector)
  count=0
  for x in range(0,cube['npoints'][0]):
    for y in range(0,cube['npoints'][1]):
      for z in range(0,cube['npoints'][2]):
        cube['data'][x,y,z]=vector[count]
        count+=1
        #if count >= nread:
        #  break;
      #if count>=nread:
      #  break
    #if count >= nread:
    #  break

  return cube

def write_cube(cube, filename):
  f=open(filename,'w')
  f.write(cube['comment'])
  f.write(cube['type'])
  natoms=cube['natoms']
  
  f.write(" %i "%natoms)
  f.write("\n")
  for i in range(0,3):
    f.write("%i "%cube['npoints'][i])
    f.write(" %g %g %g \n"%(cube['latvec'][i,0],cube['latvec'][i,1],cube['latvec'][i,2]))
  for i in range(0,natoms):
    f.write("%s 0.0 "%cube['atomname'][i])
    f.write(" %g %g %g \n"%(cube['atomxyz'][i,0],cube['atomxyz'][i,1],cube['atomxyz'][i,2]))
  count=0
  for x in range(0,cube['npoints'][0]):
    for y in range(0,cube['npoints'][1]):
      for z in range(0,cube['npoints'][2]):
        f.write("%g "%cube['data'][x,y,z])
        count+=1
        if count%5==0:
          f.write('\n')
  f.write('\n')
  
  

def write_xsf(cube,filename):
  f=open(filename,'w')
  f.write("CRYSTAL\n")
  f.write("PRIMVEC\n")
  natoms=cube['natoms']
  for i in range(0,3):
    npts=cube['npoints'][i]
    f.write(" %g %g %g \n"%(npts*cube['latvec'][i,0],npts*cube['latvec'][i,1],npts*cube['latvec'][i,2]))
  f.write("PRIMCOORD\n")
  f.write("%i 1\n"%natoms)
  for i in range(0,natoms):
    f.write("%s "%cube['atomname'][i])
    f.write(" %g %g %g \n"%(cube['atomxyz'][i,0],cube['atomxyz'][i,1],cube['atomxyz'][i,2]))
  f.write("BEGIN_BLOCK_DATAGRID_3D\n cube_file_conversion \n")
  f.write("BEGIN_DATAGRID_3D\n")
  f.write("%i %i %i\n"%(cube['npoints'][0],cube['npoints'][1],cube['npoints'][2]))
  f.write("0.0 0.0 0.0\n")
  for i in range(0,3):
    npts=cube['npoints'][i]
    f.write(" %g %g %g \n"%(npts*cube['latvec'][i,0],npts*cube['latvec'][i,1],npts*cube['latvec'][i,2]))
  
  count=0
  for z in range(0,cube['npoints'][2]):
    for y in range(0,cube['npoints'][1]):
      for x in range(0,cube['npoints'][0]):
        f.write("%g "%cube['data'][x,y,z])
        count+=1
        if count%5==0:
          f.write('\n')
  f.write('\n')
  f.write("END_DATAGRID_3D\n")
  f.write("END_BLOCK_DATAGRID_3D\n")
  


def normalize_abs(cube):
  vol=np.abs(np.linalg.det(cube['latvec']))
  norm=np.sum(np.abs(cube['data']))*vol
  cube['data']/=norm
  return cube
