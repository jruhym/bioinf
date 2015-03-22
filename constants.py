from numpy import array

# An array with the aliases for the 20 common amino acids in proteins.
# Last element for each amino acid is it's group as designated in 
# Sec. 2.4.2 of Ray et al. Strucl. Bioinf. vol. 26 p. 3067 (2010) 
aminoacids = array([
     ['A', 'ALA', 'Alanine', 5], 
     ['V', 'VAL', 'Valine', 5], 
     ['L', 'LEU', 'Leucine', 5],
     ['I', 'ILE', 'Isoleucine', 5], 
     ['P', 'PRO', 'Proline', 6], 
     ['F', 'PHE', 'Phenylalanine', 3], 
     ['W', 'TRP', 'Tryptophan', 3], 
     ['M', 'MET', 'Methionine', 5], 
     ['G', 'GLY', 'Glycine', 6], 
     ['S', 'SER', 'Serine', 4], 
     ['T', 'THR', 'Threonine', 4], 
     ['C', 'CYS', 'Cysteine', 5], 
     ['Y', 'TYR', 'Tyrosine', 3], 
     ['N', 'ASN', 'Asparagine', 4], 
     ['Q', 'GLN', 'Glutamine', 4], 
     ['D', 'ASP', 'Aspartic Acid', 2], 
     ['E', 'GLU', 'Glutamic Acid', 2], 
     ['K', 'LYS', 'Lysine', 1], 
     ['R', 'ARG', 'Arginine', 1], 
     ['H', 'HIS', 'Histidine', 3]
])


class AminoAcidPrimitive(object):

     def __init__(self, letter, abbr, name, group):
          self._letter = letter   
          self._abbr = abbr
          self._name = name
          self._group = group


     letter = property(lambda self: self._letter)
     abbr = property(lambda self: self._abbr)
     name = property(lambda self: self._name)
     group = property(lambda self: self._group)

class AminoAcid(AminoAcidPrimitive):

     _all_values = (
          AminoAcidPrimitive('A', 'ALA', 'Alanine', 5),
          AminoAcidPrimitive('V', 'VAL', 'Valine', 5),
          AminoAcidPrimitive('L', 'LEU', 'Leucine', 5),
          AminoAcidPrimitive('I', 'ILE', 'Isoleucine', 5),
          AminoAcidPrimitive('P', 'PRO', 'Proline', 6),
          AminoAcidPrimitive('F', 'PHE', 'Phenylalanine', 3),
          AminoAcidPrimitive('W', 'TRP', 'Tryptophan', 3),
          AminoAcidPrimitive('M', 'MET', 'Methionine', 5),
          AminoAcidPrimitive('G', 'GLY', 'Glycine', 6),
          AminoAcidPrimitive('S', 'SER', 'Serine', 4),
          AminoAcidPrimitive('T', 'THR', 'Threonine', 4),
          AminoAcidPrimitive('C', 'CYS', 'Cysteine', 5),
          AminoAcidPrimitive('Y', 'TYR', 'Tyrosine', 3),
          AminoAcidPrimitive('N', 'ASN', 'Asparagine', 4),
          AminoAcidPrimitive('Q', 'GLN', 'Glutamine', 4),
          AminoAcidPrimitive('D', 'ASP', 'Aspartic Acid', 2),
          AminoAcidPrimitive('E', 'GLU', 'Glutamic Acid', 2),
          AminoAcidPrimitive('K', 'LYS', 'Lysine', 1),
          AminoAcidPrimitive('R', 'ARG', 'Arginine', 1),
          AminoAcidPrimitive('H', 'HIS', 'Histidine', 3),
     )

     _by_abbr = {}

     for acid in _all_values:
          _by_abbr[acid.abbr] = acid

     @staticmethod
     def by_abbr(abbr):
          assert isinstance(abbr, basestring)
          return AminoAcid._by_abbr[abbr.upper()]
