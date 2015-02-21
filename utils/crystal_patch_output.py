import sys
outname=sys.argv[1]
propname=sys.argv[2]
patchname=sys.argv[3]


prop=open(propname,'r')
shrink=[1,1,1]
for line in prop:
  if "SHRINK FACTORS(MONK.)" in line:
    spl=line.split()
    shrink[0]=int(spl[2])
    shrink[1]=int(spl[3])
    shrink[2]=int(spl[4])


patch=open(patchname,'w')

out=open(outname,'r')
for line in out:
  if "SHRINK. FACT.(MONKH.)" in line:
    patch.write("SHRINK. FACT.(MONKH.)    %i  %i  %i  NUMBER OF K POINTS IN THE IBZ    XXX\n"%(shrink[0],shrink[1],shrink[2]))
  else:
    patch.write(line)
    
  if "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT EDFT" in line:
    break
out.close()

prop=open(propname,'r')
patch.write("NEWK EIGENVECTORS\n \n")
found_hamil=False
for line in prop:
  if "HAMILTONIAN EIGENVECTORS" in line:
    found_hamil=True
  if found_hamil:
    patch.write(line)

