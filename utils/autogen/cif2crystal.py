from __future__ import print_function
from __future__ import absolute_import
from __future__ import division
from __future__ import unicode_literals

from xml.etree.ElementTree import ElementTree
from pymatgen.io.cif import CifParser
from pymatgen.core.periodic_table import Element
import os
from io import StringIO 
import sys
import shutil
import string
import numpy as np
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer


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
        list of lines of pseudopotential section (edit by Brian Busemeyer).
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
    strlist = []
    strlist.append('INPUT')
    strlist.append(' '.join(map(str,[eff_core_charge,nlocal,
                                     m[0],m[1],m[2],m[3],m[4]])))
    for local_component in local_list:
        exp_gaus = local_component.find('./Exp').text
        coeff_gaus = local_component.find('./Coeff').text
        r_to_n = local_component.find('./r_to_n').text
        strlist.append(' '.join([exp_gaus, coeff_gaus,r_to_n]))
    for non_local_component in non_local_list:
        exp_gaus = non_local_component.find('./Exp').text
        coeff_gaus = non_local_component.find('./Coeff').text
        r_to_n = non_local_component.find('./r_to_n').text
        strlist.append(' '.join([exp_gaus, coeff_gaus,r_to_n]))
    return strlist


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
    if symbol in initial_charges.keys():
      atom_charge-=initial_charges[symbol]
    basis_path = './Basis-set[@name="{}"]/Contraction'.format(basis_name)
    found_orbitals = []
    totcharge=0
    ret=[]
    ncontract=0
    for contraction in element.findall(basis_path):
        angular = contraction.get('Angular_momentum')
        if found_orbitals.count(angular) >= nangular[angular]:
            continue

        #Figure out which coefficients to print out based on the minimal exponent
        nterms = 0
        basis_part=[]
        for basis_term in contraction.findall('./Basis-term'):
            exp = basis_term.get('Exp')
            coeff = basis_term.get('Coeff')
            if float(exp) > cutoff:
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

    #Add in the uncontracted basis elements
    angular_uncontracted=['s','p']
    if symbol in transition_metals:
        angular_uncontracted.append('d')

    for angular in angular_uncontracted:
        for i in range(0,params[1]):
            exp=params[0]*params[2]**i
            line='{} {}'.format(exp,1.0)
            ret+=["0 %i %i %g 1"%(basis_index[angular],1,0.0),line]
            ncontract+=1

    return ["%i %i"%(Element(symbol).number+200,ncontract)] +\
            pseudopotential_section(symbol,xml_name) +\
            ret

######################################################################

def basis_section(struct,params=[0.2,2,3],initial_charges={},library_directory="../"):
  sites=struct.as_dict()['sites']
  elements=set()
  for s in sites:
    nm=s['species'][0]['element']
    elements.add(nm)
  basislines=[]
  elements = sorted(list(elements)) # Standardize ordering.
  for e in elements:
    basislines+=modify_basis(e,library_directory+"BFD_Library.xml",
            params=params,initial_charges=initial_charges)
  return basislines


######################################################################

def cif2geom(cif):
  parser=CifParser.from_string(cif)
  struct=parser.get_structures()[0]
  primstruct=struct
  lat=primstruct.as_dict()['lattice']
  sites=primstruct.as_dict()['sites']
  geomlines=["CRYSTAL","0 0 0","1","%g %g %g %g %g %g"%\
          (lat['a'],lat['b'],lat['c'],lat['alpha'],lat['beta'],lat['gamma'])]

  geomlines+=["%i"%len(sites)]
  for v in sites:
      nm=v['species'][0]['element']
      nm=str(Element(nm).Z+200)
      geomlines+=[nm+" %g %g %g"%(v['abc'][0],v['abc'][1],v['abc'][2])]

  return geomlines,primstruct

######################################################################
def cif2geom_sym(cif):
  def parse(s):
    return str(s).split('(')[0]
  parser=CifParser.from_string(cif)
  struct=parser.get_structures()[0]
  primstruct=struct.get_primitive_structure()  
  cifdict=parser.as_dict().popitem()[1]
  sym_num=int(cifdict['_symmetry_Int_Tables_number'])
  sg_name=cifdict['_symmetry_space_group_name_H-M']
  cell_parms=[]
  if sym_num in range(0,3):# Triclinic
    cell_parms=['_cell_length_a','_cell_length_b','_cell_length_c',
                '_cell_angle_alpha', '_cell_angle_beta', '_cell_angle_gamma']
  elif sym_num in range(3,16): #Monoclinic
    cell_parms=['_cell_length_a','_cell_length_b','_cell_length_c',
                '_cell_angle_alpha']
  elif sym_num in range(16,75): #Orthorhombic
    cell_parms=['_cell_length_a','_cell_length_b','_cell_length_c']
  elif sym_num in range(75,143): #Tetragonal
    cell_parms=['_cell_length_a','_cell_length_c']
  elif sym_num in range(143,168): #Trigonal
    cell_parms=['_cell_length_a','_cell_length_c']
    print("WARNING: Trigonal is untested and may have some bugs.")
  elif sym_num in range(168,195): #Hexagonal
    cell_parms=['_cell_length_a','_cell_length_c']
  elif sym_num in range(195,231): #Cubic
    cell_parms=['_cell_length_a']

  latconst=""
  for p in cell_parms:
    latconst+= parse(cifdict[p])+" "
  geomlines=["CRYSTAL","0 0 1",str(sym_num),latconst,str(len(cifdict['_atom_site_label']))]

  for nm,x,y,z in zip(cifdict['_atom_site_label'],
                      cifdict['_atom_site_fract_x'],
                      cifdict['_atom_site_fract_y'],
                      cifdict['_atom_site_fract_z'] ):
    nm=''.join(c for c in nm if not c.isdigit())
    #print(x,y,z)
    geomlines.append(str(Element(nm).Z+200)+" "+parse(x) + " " + parse(y) + " " + parse(z))
  return geomlines,primstruct


######################################################################

def cif2geom_sym2(cif):
  parser=CifParser.from_string(cif)
  struct=parser.get_structures()[0]
  sg = SpacegroupAnalyzer(struct)
  struct = sg.get_conventional_standard_structure()
  sg = SpacegroupAnalyzer(struct)

  geomlines=["CRYSTAL"]
  geomlines += ["0 0 1"]
  geomlines += [str(sg.get_spacegroup_number())]
  cry_sys = sg.get_crystal_system()
  lattice = struct.lattice

  if cry_sys == 'trigonal' or cry_sys == 'hexagonal' or cry_sys == 'tetragonal':
    geomlines += ["%s %s" %(lattice.a,lattice.c)]
  elif cry_sys == 'cubic':
    geomlines += ["%s" %(lattice.a)]
  elif cry_sys == 'triclinic':
    geomlines += ["%s %s %s %s %s %s" %(lattice.a,lattice.b,lattice.c,lattice.alpha,lattice.beta,lattice.gamma)]
  elif cry_sys == 'monoclinic':
    geomlines += ["%s %s %s %s" %(lattice.a,lattice.b,lattice.c,lattice.beta)]
  elif cry_sys == 'orthorhombic':
    geomlines += ["%s %s %s" %(lattice.a,lattice.b,lattice.c)]
  else:
    print('Error printing symmetrized structure.')
    quit()
  
  ds = sg.get_symmetry_dataset()
  eq_sites = np.unique(ds['equivalent_atoms'])
  geomlines += [str(len(eq_sites))]
  for eq_site in eq_sites:
    site = struct.sites[eq_site]
    geomlines += ["%s %s %s %s" %(site.specie.Z+200,site.a,site.b,site.c)]

  return geomlines,struct

######################################################################

class Cif2Crystal:
  _name_="Cif2Crystal"

  def __init__(self,library_directory="../"):
    self.library_directory=library_directory
  # Currently, check_status() and runcrystal requires user not to modify outfn. 
  # In the future, we should probably store name of d12 in record, 
  # so this is not needed.
  def run(self,job_record,outfn="autogen.d12"):
    if job_record['pseudopotential']!='BFD':
      print("ERROR: only support BFD pseudoptentials for now")
      quit()

    if job_record['dft']['symmetrized'] == False:
      geomlines,primstruct=cif2geom(job_record['cif'])
    elif job_record['dft']['symmetrized'] == True:
      geomlines,primstruct=cif2geom_sym2(job_record['cif'])
    else:
      print('Improper input for "symmetrized" keyword. Value must be True or False.')
      quit()
    basislines=basis_section(primstruct,job_record['dft']['basis'],
                              job_record['dft']['initial_charges'],
                              self.library_directory)
    supercell=["SUPERCEL"]
    for row in job_record['supercell']:
      supercell+=[' '.join(map(str,row))]

    modisym=[]
    if len(job_record['dft']['initial_spin']) > 0:
      ninit=len(job_record['dft']['initial_spin'])
      #modisym += ["MODISYMM",str(ninit)]
      count=1
      for i,s in enumerate(job_record['dft']['initial_spin']):
        if s!=0:
          modisym += [str(i+1)+" "+str(count)+" "]
          count+=1
      if count > 0:
        modisym.insert(0,str(count-1))
        modisym.insert(0,"MODISYMM")

    # List of lines will be joined by '\n'.
    outlines = ["Generated by cif2crystal"] +\
                geomlines +\
                supercell + \
                modisym + \
                ["END"] +\
                basislines +\
                ["99 0"] +\
                ["CHARGED"] +\
                ["END"]
    if min(job_record['dft']['kmesh']) != max(job_record['dft']['kmesh']):
      print("Non-cubic meshes not implemented")
      quit()
    else:
      val = job_record['dft']['kmesh'][0]
      outlines += [
        "SHRINK",
        str(val)+' '+str(2*val)
      ]
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
      "SAVEWF",
      "BIPOSIZE",
      "100000000",
      "EXCHSIZE",
      "10000000",
      "TOLDEE",
      str(job_record['dft']['edifftol']),
      "FMIXING",
      str(job_record['dft']['fmixing']),
      "TOLINTEG",
      ' '.join(map(str,job_record['dft']['tolinteg'])),
      "MAXCYCLE",
      str(job_record['dft']['maxcycle'])
    ]
    
    if job_record['dft']['levshift']!=None:
      outlines+=["LEVSHIFT",' '.join(map(str,job_record['dft']['levshift']))]
    else:
      outlines+=["BROYDEN",
                 ' '.join(map(str,job_record['dft']['broyden']))]
    

    if job_record['dft']['smear'] != None:
      outlines += ["SMEAR",str(job_record['dft']['smear'])]
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
    new = open("new.autogen.d12",'r').read().split()
    old = open("autogen.d12",'r').read().split()
    for doesntmatter in ["GUESSP","SAVEWF"]:
      if doesntmatter in new: new.remove(doesntmatter)
      if doesntmatter in old: old.remove(doesntmatter)
    if new != old:
      if not job_record['assert_nochanges']:
        print("Warning: job record inconsistent with past input")
        return 'ok'
      else:
        print("Error: job record inconsistent with past input")
        return 'failed'
    else:
      return 'ok'

  def retry(self,job_record):
    return self.run(job_record)
  def output(self,job_record):
    return job_record
        


######################################################################

if __name__ == "__main__":
  geomstring,primstruct=cif2geom(sys.argv[1])
  basisstring=basis_section(primstruct,'.')
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
