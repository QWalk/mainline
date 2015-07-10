from __future__ import print_function
from xml.etree.ElementTree import ElementTree
from collections import OrderedDict
import pymatgen as mg
from pymatgen.io.cifio import CifParser
from pymatgen.core.periodic_table import Element
import os
import cStringIO 
import sys

library_directory="/projects/wagner/apps/qwalk_src/utils/autogen/" #this is relative to the run directorys. You can also put an absolute position

######################################################################
def pseudopotential_section(symbol, xml_name):
    """
    Author: "Kittithat (Mick) Krongchon" <pikkienvd3@gmail.com>
    Returns a string of the pseudopotential section which is to be written
    as part of the basis set section.

    Args:
        symbol (str): The symbol of the element to be specified in the
            D12 file.
        xml_name (str): The name of the XML pseudopotential and basis
            set database.

    Returns:
        str: The pseudopotential section.
    """
    tree = ElementTree()
    tree.parse(xml_name)
    element = tree.find('./Pseudopotential[@symbol="{}"]'.format(symbol))
    eff_core_charge = element.find('./Effective_core_charge').text
    local_path = './Gaussian_expansion/Local_component'
    non_local_path = './Gaussian_expansion/Non-local_component'
    local_list = element.findall(local_path)
    non_local_list = element.findall(non_local_path)
    nlocal = len(local_list)
    m = [0, 0, 0, 0, 0]
    proj_path = './Gaussian_expansion/Non-local_component/Proj'
    proj_list = element.findall(proj_path)
    for projector in proj_list:
        m[int(projector.text)] += 1
    pseudopot_strlist = []
    ecp_keyword = 'INPUT'
    pseudopot_strlist.append('{}\n'.format(ecp_keyword))
    pseudopot_strlist.append('{}    '.format(eff_core_charge))
    pseudopot_strlist.append('{} {} {} '.format(nlocal, m[0], m[1]))
    pseudopot_strlist.append('{} {} {}\n'.format(m[2], m[3], m[4]))
    for local_component in local_list:
        exp_gaus = local_component.find('./Exp').text
        coeff_gaus = local_component.find('./Coeff').text
        r_to_n = local_component.find('./r_to_n').text
        pseudopot_strlist.append('  {} {} '.format(exp_gaus, coeff_gaus))
        pseudopot_strlist.append('{}\n'.format(r_to_n))
    for non_local_component in non_local_list:
        exp_gaus = non_local_component.find('./Exp').text
        coeff_gaus = non_local_component.find('./Coeff').text
        r_to_n = non_local_component.find('./r_to_n').text
        pseudopot_strlist.append('  {} {} '.format(exp_gaus, coeff_gaus))
        pseudopot_strlist.append('{}\n'.format(r_to_n))
    pseudopot_string = ''.join(pseudopot_strlist)
    return pseudopot_string


######################################################################
def modify_basis(symbol,xml_name,cutoff=0.2,params=[0.2,2,3],initial_charges={}):
    """
    Author: "Kittithat (Mick) Krongchon" <pikkienvd3@gmail.com> and Lucas K. Wagner
    Returns a string containing the basis section.  It is modified according to a simple recipe:
    1) The occupied atomic orbitals are kept, with exponents less than 'cutoff' removed.
    2) These atomic orbitals are augmented with uncontracted orbitals according to the formula 
        e_i = params[0]*params[2]**i, where i goes from 0 to params[1]
        These uncontracted orbitals are added for every occupied atomic orbital (s,p for most elements and s,p,d for transition metals)

    Args:
        symbol (str): The symbol of the element to be specified in the
            D12 file.
        xml_name (str): The name of the XML pseudopotential and basis
            set database.
        cutoff: smallest allowed exponent
        params: parameters for generating the augmenting uncontracted orbitals

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
    print(initial_charges)
    if symbol in initial_charges.keys():
      atom_charge-=initial_charges[symbol]
    basis_path = './Basis-set[@name="{}"]/Contraction'.format(basis_name)
    found_orbitals = []
    totcharge=0
    ret=""
    ncontract=0
    for contraction in element.findall(basis_path):
        angular = contraction.get('Angular_momentum')
        if found_orbitals.count(angular) >= nangular[angular]:
            continue
        found_orbitals.append(angular)

        #Figure out which coefficients to print out based on the minimal exponent
        nterms = 0
        basis_part=""
        for basis_term in contraction.findall('./Basis-term'):
            exp = basis_term.get('Exp')
            coeff = basis_term.get('Coeff')
            if float(exp) > cutoff:
                line = '  {} {}\n'.format(exp, coeff)
                basis_part+=line
                nterms+=1
        #now write the header 
        if nterms > 0:
          charge=min(atom_charge-totcharge,maxcharge[angular])
          #put in a special case for transition metals:
          #depopulate the 4s if the atom is charged
          if symbol in transition_metals and symbol in initial_charges.keys() \
                  and initial_charges[symbol] > 0 and found_orbitals.count(angular) > 1 \
                  and angular=="s":
            charge=0
          totcharge+=charge
          ret+="0 %i %i %g 1\n"%(basis_index[angular],nterms,charge) +\
                basis_part
          ncontract+=1

    #Add in the uncontracted basis elements
    angular_uncontracted=['s','p']
    if symbol in transition_metals:
        angular_uncontracted.append('d')

    for angular in angular_uncontracted:
        for i in range(0,params[1]):
            exp=params[0]*params[2]**i
            line='  {} {}\n'.format(exp,1.0)
            ret+="0 %i %i %g 1\n"%(basis_index[angular],1,0.0) +\
                line
            ncontract+=1

    return "%i %i\n"%(mg.Element(symbol).number+200,ncontract) +\
            pseudopotential_section(symbol,xml_name) +\
            ret

######################################################################

def basis_section(struct,params=[0.2,2,3],initial_charges={}):
  sites=struct.as_dict()['sites']
  elements=set()
  for s in sites:
    nm=s['species'][0]['element']
    elements.add(nm)
  basisstring=""
  for e in elements:
    basisstring+=modify_basis(e,library_directory+"BFD_Library.xml",
            params=params,initial_charges=initial_charges)
  return basisstring


######################################################################

def cif2geom(ciffile):
  parser=CifParser(ciffile)
  struct=parser.get_structures()[0]
  primstruct=struct.get_primitive_structure()
  lat=primstruct.as_dict()['lattice']
  sites=primstruct.as_dict()['sites']
  geomstring="CRYSTAL\n0 0 0 \n1\n%g %g %g %g %g %g\n"%\
          (lat['a'],lat['b'],lat['c'],lat['alpha'],lat['beta'],lat['gamma'])

  geomstring+="%i\n"%len(sites)
  for v in sites:
      nm=v['species'][0]['element']
      nm=str(Element(nm).Z+200)
      geomstring+=nm+" %g %g %g\n"%(v['abc'][0],v['abc'][1],v['abc'][2])

  return geomstring,primstruct

######################################################################

class Cif2Crystal:
  _name_="Cif2Crystal"
  def run(self,job_record):
    #TODO: support  kmesh,charge
    if job_record['pseudopotential']!='BFD':
      print("ERROR: only support BFD pseudoptentials for now")
      quit()

    geomstring,primstruct=cif2geom(cStringIO.StringIO(job_record['cif']))
    basisstring=basis_section(primstruct,job_record['dft']['basis'],
                              job_record['dft']['initial_charges'])
    supercell="SUPERCEL\n"
    for row in job_record['supercell']:
      supercell+="%i %i %i\n"%(row[0],row[1],row[2])

    f=open("autogen.d12",'w')
    f.write("Generated by cif2crystal\n"+geomstring+supercell+"END\n"+\
          basisstring+"99 0\nCHARGED\nEND\n")
    f.write("""SHRINK
8 16
DFT\n""")
    if job_record['dft']['spin_polarized']:
      f.write("SPIN\n")
    f.write("""EXCHANGE\n"""+
job_record['dft']['functional']['exchange'] +
"\nCORRELAT\n"+
job_record['dft']['functional']['correlation'] +
"\nHYBRID\n"+
str(job_record['dft']['functional']['hybrid']) +
"""\nXLGRID
END
SCFDIR
BIPOSIZE
100000000
FMIXING
"""
+str(job_record['dft']['fmixing'])+"""
BROYDEN\n"""+
"%g %i %i\n"%(job_record['dft']['broyden'][0],job_record['dft']['broyden'][1],job_record['dft']['broyden'][2]) +
"""TOLINTEG 
12 12 12 12 20
""")
    if job_record['dft']['spin_polarized']:
      f.write("SPINLOCK\n") 
      f.write(str(job_record['total_spin'])+' 200 \n')
      if len(job_record['dft']['initial_spin']) > 0:
        ninit=len(job_record['dft']['initial_spin'])
        f.write("ATOMSPIN\n"+str(ninit)+'\n')
        for i,s in enumerate(job_record['dft']['initial_spin']):
          f.write(str(i+1)+" "+str(s)+" ")
        f.write('\n')
    f.write("END\n")
    f.close()
    return 'ok'
  def check_status(self,job_record):
    if os.path.isfile("autogen.d12"):
      return 'ok'
    else:
      return 'not_started'
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
