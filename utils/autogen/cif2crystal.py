from __future__ import print_function
from xml.etree.ElementTree import ElementTree
from collections import OrderedDict
import pymatgen as mg
from pymatgen.io.cifio import CifParser
from pymatgen.core.periodic_table import Element
import sys

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
def angular_momentum(string):
    """
    Returns a corresponding CRYSTAL input value and the maximum number of
    electrons in the shell of the orbital angular momentum represented by
    the string parameter.

    Args:
        string (str): The orbital angular momentum.

    Returns:
        str: The number represents the orbital as specified in the
            CRYSTAL manual.
        int: The maximum number of electrons in the given shells.
    """
    if string == 's':
        return '0', 2
    elif string == 'sp':
        return '1', 8
    elif string == 'p':
        return '2', 6
    elif string == 'd':
        return '3', 10
    elif string == 'f':
        return '4', 14
    elif string == 'g':
        return '5', -99
    elif string == 'h':
        return '6', -99
    elif string == 'i':
        return '7', -99
    else:
        return 'error', -99
######################################################################
def fill(charge_left, cap):
    """
    Fills shell with the effective core charges and returns the number of
    charges left and charges in the shell.
    """
    if charge_left <= cap:
        charge = charge_left
        charge_left = 0
    else:
        charge_left -= cap
        charge = cap
    return charge_left, charge
######################################################################

def basis_modified(symbol, xml_name, basis_name, params):
    """
    Author: "Kittithat (Mick) Krongchon" <pikkienvd3@gmail.com>
    
    Returns a string of the modified basis set section for D12 input file.

    Args:
        symbol (str): The symbol of the element to be specified
            in the D12 file
        xml_name (str): The name of the XML pseudopotential
            and basis set database
        basis_name (str): The name of the basis set (e.g. 'vtz')
        params (list): [coeff, n, base]
            coeff (double): aka alpha. The coefficient of the exponential
                operation.
            n (int): The exponent.
            base (double): aka c. The base of the exponential operation.

    Returns:
        str: The basis set section.
    """
    tree = ElementTree()
    tree.parse(xml_name)
    element = tree.find('./Pseudopotential[@symbol="{}"]'.format(symbol))
    eff_core_charge = int(element.find('./Effective_core_charge').text)
    basis_path = './Basis-set[@name="{}"]/Contraction'.format(basis_name)
    contractions = OrderedDict()
    found_orbitals = []
    for contraction in element.findall(basis_path):
        angular = contraction.get('Angular_momentum')
        if angular in found_orbitals or not angular in 'spdf':
            continue
        found_orbitals.append(angular)
        contractions[angular] = []
        indent = 0
        for basis_term in contraction.findall('./Basis-term'):
            exp = basis_term.get('Exp')
            coeff = basis_term.get('Coeff')
            if float(exp) > params[0]:
                line = '  {} {}\n'.format(exp, coeff)
                if indent < 2:
                    contractions[angular].append('  ' + line)
                else:
                    contractions[angular].append(line)
                indent = (indent + 1) % 4
    result = []
    charge_left = eff_core_charge
    ncontractions = 0

    maxorb=3
    if symbol in ["Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn"]:
      maxorb=4
    for angular, basis_terms in contractions.items():
        orbital, cap = angular_momentum(angular)
        charge_left, charge = fill(charge_left, cap)
        l = orbital, len(contractions[angular]), charge
        if basis_terms:
            result.append('0 {} {} {} 1\n'.format(l[0], l[1], l[2]))
            ncontractions += 1
            for basis_term in basis_terms:
                result.append(basis_term)
        else:
            break
        if int(orbital) < maxorb:
            for i in range(params[1]):
                result.append('0 {} 1 0 1\n'.format(orbital))
                ncontractions += 1
                exp = params[0] * params[2] ** i
                result.append('    {} 1.000000\n'.format(exp))
    atomic = '2{:02}'.format(mg.Element(symbol).number)
    result = (['{} {}\n'.format(atomic, ncontractions),
               pseudopotential_section(symbol, xml_name)] +
               result)
    return ''.join(result)

######################################################################

def basis_section(struct):
  sites=struct.as_dict()['sites']
  elements=set()
  for s in sites:
    nm=s['species'][0]['element']
    elements.add(nm)
  basisstring=""
  for e in elements:
    basisstring+=basis_modified(e,"BFD_Library.xml","vtz",[0.2,2,3])
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
BIPOSIZE
100000000
FMIXING
99
BROYDEN
.01 60 5
END
""")
