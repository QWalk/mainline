from pymatgen.io.cifio import CifParser,CifWriter
import numpy as np

def move_atom(input_cif, output_cif, atom_idx, displacement_vec):
  """ 
  Reads input cif (as file name), moves atom by displacement vector (in
  angstroms), and writes to output_cif (as file name).

  Currently, doesn't take into account the symmetry operations. This is
  theoretically possible with PyMatGen, although I found it to be buggy so
  decided it might not be worth it. 
  Useful things to look up for implementing this:

  pymatgen.symmetry.analyzer.SpacegroupAnalyzer.get_symmetrized_structure

  and its resulting method

  find_equivalent_sites()
  """
  parser    = CifParser(input_cif)
  struct    = parser.get_structures()[0]
  lat_const = np.array(struct.lattice.abc)
  frac_pos  = struct[atom_idx].frac_coords
  new_pos   = frac_pos + displacement_vec/lat_const
  struct[atom_idx] = struct[atom_idx].species_and_occu.elements[0], new_pos
  writer = CifWriter(struct)
  writer.write_file(output_cif)
