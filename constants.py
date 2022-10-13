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

     _by_abbr = { acid.abbr : acid for acid in _all_values }
     _by_letter = { acid.letter : acid for acid in _all_values }

     for acid in _all_values:
          _by_abbr[acid.abbr] = acid

     @staticmethod
     def by_abbr(abbr):
          assert isinstance(abbr, basestring)
          return AminoAcid._by_abbr[abbr.upper()]

     @staticmethod   
     def by_letter(letter):
        assert isinstance(letter, str)
        return AminoAcid._by_letter[letter.upper()]

class Residue:
    def __init__(self, sequence_number, amino_acid, protein_segment, bw=None):
        self._sequence_number = sequence_number
        self._amino_acid = amino_acid   
        self._bw = bw
        self._protein_segment = protein_segment

    def __init__(self, gpcrdb_dict):
        self._sequence_number = gpcrdb_dict.get('sequence_number', None)
        self._amino_acid = AminoAcid.by_letter(gpcrdb_dict.get('amino_acid', None))
#       https://stackoverflow.com/questions/2492087/how-to-get-the-nth-element-of-a-python-list-or-a-default-if-not-available
        try:
            self._bw = [a.get('label', None) for a in gpcrdb_dict.get('alternative_generic_numbers', []) if a.get('scheme', None) == 'BW'][0]
        except IndexError:
            self._bw = None
        self._protein_segment = gpcrdb_dict.get('protein_segment', None)

    def __str__(self):
        return "{abbr_seq_num:s} : {bw:s}".format(abbr_seq_num = self.abbr_seq_num, bw = (self.bw or self.protein_segment))

    sequence_number = property(lambda self: self._sequence_number)
    amino_acid = property(lambda self: self._amino_acid)
    bw = property(lambda self: self._bw)
    protein_segment = property(lambda self: self._protein_segment)
    lett_seq_num = property(lambda self: f'{self._amino_acid.letter}{self._sequence_number}')
    abbr_seq_num = property(lambda self: f'{self._amino_acid.abbr}{self._sequence_number}')
    helix = property(lambda self: self._bw.split('.')[0] if self._bw else None)

from pathlib import Path
import requests
import json

class Receptor:
    def __residues__(self, receptor):
        assert isinstance(receptor, str)
        directory = Path.home().joinpath(Path('GPCRDB/'))
        Path.mkdir(directory, exist_ok=True)
        file = directory.joinpath(f'{receptor}{".json"}')
        try:
            fp = file.open()
            responseJSON = json.load(fp)
        except:
            response = requests.get(f'{"https://gpcrdb.org/services/residues/extended/"}{receptor}')
            if response.ok:    
                fp = file.open('w')
                responseJSON = response.json()
                json.dump(responseJSON, fp)
        finally:
            fp.close()
            return { Residue(residueDict).abbr_seq_num : Residue(residueDict) for residueDict in responseJSON }

    def __init__(self, uniprot_name):
        self._uniprot_name_ = uniprot_name
        self._residues_ = self.__residues__(self._uniprot_name_)
    
    def residue_for_index(self, index):
        assert isinstance(index, str)
        a = self._residues_.get(index, None)
        if a: 
            return a
        ai = None
        a = index
        while not ai and a:
            try:
                ai = int(a)
            except:
                a = a[1:]
        for key,value in self._residues_.items():
            if key[3:] == a:
                return value
        return None
    
    uniprot_name = property(lambda self: self._uniprot_name_)
    residues = property(lambda self: self._residues_)
