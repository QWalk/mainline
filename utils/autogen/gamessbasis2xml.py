from __future__ import print_function
import sys

f=open(sys.argv[1],'r')

for line in f:
  spl=line.split()
  
  if len(spl) > 1 and ("S" in spl[0] or "s" in spl[0] \
          or "P" in spl[0] or "p" in spl[0] \
          or "D" in spl[0] or "d" in spl[0] \
          or "F" in spl[0] or "f" in spl[0]):
    print('</Contraction>')
    print('<Contraction nterms="%s" Angular_momentum="%s">'%(spl[1],spl[0].lower()))
  elif len(spl) > 1:
    print('<Basis-term Exp="%s" Coeff="%s"/>'%(spl[0],spl[1]))

    
print('</Contraction>')

