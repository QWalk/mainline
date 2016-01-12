from __future__ import print_function
from __future__ import absolute_import
from __future__ import division
from __future__ import unicode_literals

from xml.etree.ElementTree import ElementTree
from pymatgen.core.periodic_table import Element
import os
from io import StringIO 
import sys
import shutil
import string
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from cif2crystal import pseudopotential_section
import pymatgen

library_directory="../"

######################################################################
def generate_basis(symbol,xml_name,initial_charges={}):
    """
    Author: "Kittithat (Mick) Krongchon" <pikkienvd3@gmail.com> and Lucas K. Wagner
    Returns a string containing the basis section.  It is modified according to a simple recipe:
    Args:
        symbol (str): The symbol of the element to be specified in the
            D12 file.
        xml_name (str): The name of the XML pseudopotential and basis
            set database.
    Returns:
        str: The pseudopotential and basis section.
    """
    
    maxorb=3
    basis_name="vtz"
    nangular={"s":1,"p":1,"d":1,"f":1,"g":0}
    maxcharge={"s":2,"p":6,"d":10,"f":15}
    basis_index={"s":0,"p":2,"d":3,"f":4}
    transition_metals=["Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn"]
    if symbol in transition_metals:
      maxorb=4
      nangular['s']=2
    
    tree = ElementTree()
    tree.parse(xml_name)
    element = tree.find('./Pseudopotential[@symbol="{}"]'.format(symbol))
    atom_charge = int(element.find('./Effective_core_charge').text)
    if symbol in initial_charges.keys():
      atom_charge-=initial_charges[symbol]
    basis_path = './Basis-set[@name="{}"]/Contraction'.format(basis_name)
    found_orbitals = []
    totcharge=0
    ret=[]
    ncontract=0
    for contraction in element.findall(basis_path):
        angular = contraction.get('Angular_momentum')

        #Figure out which coefficients to print out based on the minimal exponent
        nterms = 0
        basis_part=[]
        for basis_term in contraction.findall('./Basis-term'):
            exp = basis_term.get('Exp')
            coeff = basis_term.get('Coeff')
            basis_part += ['  {} {}'.format(exp, coeff)]
            nterms+=1
        #now write the header 
        if nterms > 0:
          found_orbitals.append(angular)          
          charge=min(atom_charge-totcharge,maxcharge[angular])
          #put in a special case for transition metals:
          #depopulate the 4s if the atom is charged
          if symbol in transition_metals and symbol in initial_charges.keys() \
                  and initial_charges[symbol] > 0 and found_orbitals.count(angular) > 1 \
                  and angular=="s":
            charge=0
          totcharge+=charge
          ret+=["0 %i %i %g 1"%(basis_index[angular],nterms,charge)] + basis_part
          ncontract+=1


    return ["%i %i"%(Element(symbol).number+200,ncontract)] +\
            pseudopotential_section(symbol,xml_name) +\
            ret

######################################################################

def basis_section(struct,params=[0.2,2,3],initial_charges={}):
  sites=struct.as_dict()['sites']
  elements=set()
  for s in sites:
    nm=s['species'][0]['element']
    elements.add(nm)
  basislines=[]
  for e in elements:
    basislines+=generate_basis(e,library_directory+"BFD_Library.xml",
            initial_charges=initial_charges)
  return basislines


######################################################################

def xyz2geom(xyzfile):
  struct=pymatgen.io.xyz.XYZ.from_string(xyzfile).molecule
  geomlines=["MOLECULE","1"]
  sites=struct.as_dict()['sites']
  geomlines+=["%i"%len(sites)]
  for v in sites:
      nm=v['species'][0]['element']
      nm=str(Element(nm).Z+200)
      print(v)
      geomlines+=[nm+" %g %g %g"%(v['xyz'][0],v['xyz'][1],v['xyz'][2])]

  return geomlines,struct

######################################################################

class Xyz2Crystal:
  _name_="Xyz2Crystal"
  # Currently, check_status() and runcrystal requires user not to modify outfn. 
  # In the future, we should probably store name of d12 in record, 
  # so this is not needed.
  def run(self,job_record,outfn="autogen.d12"):
    #TODO: support  kmesh,charge
    if job_record['pseudopotential']!='BFD':
      print("ERROR: only support BFD pseudoptentials for now")
      quit()

    geomlines,primstruct=xyz2geom(job_record['xyz'])
    basislines=basis_section(primstruct,job_record['dft']['basis'],
                              job_record['dft']['initial_charges'])
    modisym=[]
    if len(job_record['dft']['initial_spin']) > 0:
      ninit=len(job_record['dft']['initial_spin'])
      #modisym += ["MODISYMM",str(ninit)]
      count=1
      for i,s in enumerate(job_record['dft']['initial_spin']):
        if s!=0:
          modisym += [str(i+1)+" "+str(count)+" "]
          count+=1
      modisym.insert(0,str(count-1))
      modisym.insert(0,"MODISYMM")

    # List of lines will be joined by '\n'.
    outlines = ["Generated by cif2crystal"] +\
                geomlines +\
                modisym + \
                ["END"] +\
                basislines +\
                ["99 0"] +\
                ["CHARGED"] +\
                ["END"]
    outlines += ["DFT"]
    if job_record['dft']['spin_polarized']:
      outlines += ["SPIN"]
    outlines += [ 
      "EXCHANGE",
      job_record['dft']['functional']['exchange'],
      "CORRELAT",
      job_record['dft']['functional']['correlation'],
      "HYBRID", 
      str(job_record['dft']['functional']['hybrid']),
      "XLGRID",
      "END",
      "SCFDIR",
      "BIPOSIZE",
      "100000000",
      "EXCHSIZE",
      "10000000",
      "TOLDEE",
      str(job_record['dft']['edifftol']),
      "FMIXING",
      str(job_record['dft']['fmixing']),
      "BROYDEN",
      ' '.join(map(str,job_record['dft']['broyden'])),
      "TOLINTEG",
      ' '.join(map(str,job_record['dft']['tolinteg'])),
      "MAXCYCLE",
      str(job_record['dft']['maxcycle'])
    ]
    if job_record['dft']['spin_polarized']:
      outlines += [
        "SPINLOCK",
        str(job_record['total_spin'])+' 200'
      ]
      if len(job_record['dft']['initial_spin']) > 0:
        ninit=len(job_record['dft']['initial_spin'])
        outlines += ["ATOMSPIN",str(ninit)]
        for i,s in enumerate(job_record['dft']['initial_spin']):
          outlines += [str(i+1)+" "+str(s)+" "]
    if job_record['dft']['restart_from'] != None:
      try: 
        shutil.copy(job_record['dft']['restart_from'],'fort.20')
        outlines += ["GUESSP"]
      except IOError:
        print("Error: couldn't find restart file")
        return 'failed'
    outlines += ["END"]
    with open(outfn,'w') as outf:
      outf.write('\n'.join(outlines))
      outf.close()
    return 'ok'

  def check_status(self,job_record):
    if not os.path.isfile("autogen.d12"):
      return 'not_started'
    status = self.run(job_record,outfn="new.autogen.d12")
    new = open("new.autogen.d12",'r').read()
    old = open("autogen.d12",'r').read()
    if new.split() != old.split():
      print("Warning: job record inconsistent with past input")
      return 'ok'
    else:
      return 'ok'

  def retry(self,job_record):
    return self.run(job_record)
  def output(self,job_record):
    return job_record
        


######################################################################

if __name__ == "__main__":
  geomstring,primstruct=cif2geom(sys.argv[1])
  basisstring=basis_section(primstruct)
  print("Generated by cif2crystal\n",geomstring,"END\n",
          basisstring,"99 0\nEND\n",end="",sep="")
  print("""SHRINK
8 16
DFT
EXCHANGE
PBE
CORRELAT
PBE
END
SCFDIR
BIPOSIZE
100000000
EXCHSIZE 
10017422
FMIXING
99
BROYDEN
.01 60 8
END
""")
