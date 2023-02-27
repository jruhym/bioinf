import urllib.request as urlrequest
from urllib.error import URLError as urlerror
import os
import os.path
import re
from numpy import array, append # TODO: use np
import numpy as np
from bs4 import BeautifulSoup
import functools

class PDBHelixLine(object):
    @classmethod
    def parse_string(cls, string):
        assert type(string) == str 

        cleanstring = string.strip()
        # make sure string is long enough. If not, pad it on the right with ' ' 
        cleanstring += (76 - len(cleanstring)) * ' '

        assert cleanstring[:6].strip() == 'HELIX'

        return PDBHelixLine(
            cleanstring[7:10],
            cleanstring[11:14],
            cleanstring[15:18],
            cleanstring[19],
            cleanstring[21:25],
            cleanstring[25],
            cleanstring[27:30],
            cleanstring[31],
            cleanstring[33:37],
            cleanstring[37],
            cleanstring[38:40],
            cleanstring[40:70],
            cleanstring[71:]
        )

    def __init__(self, serNum, helixID, intResName, initChainID, initSeqNum, 
        initICode, endResName, endChainID, endSeqNum, endICode, helixClass, 
        comment, length
        ):
        self._serNum = serNum
        self._helixID = helixID
        self._intResName = intResName
        self._initChainID = initChainID
        self._initSeqNum = initSeqNum
        self._initICode = initICode
        self._endResName = endResName
        self._endChainID = endChainID
        self._endSeqNum = endSeqNum
        self._endICode = endICode
        self._helixClass = helixClass 
        self._comment = comment
        self._length = length

    serNum = property(lambda self: self._serNum)
    helixID = property(lambda self: self._helixID)
    intResName = property(lambda self: self._intResName)
    initChainID = property(lambda self: self._initChainID)
    initSeqNum = property(lambda self: self._initSeqNum)
    initICode = property(lambda self: self._initICode)
    endResName = property(lambda self: self._endResName)
    endChainID = property(lambda self: self._endChainID)
    endSeqNum = property(lambda self: self._endSeqNum)
    endICode = property(lambda self: self._endICode)
    helixClass = property(lambda self: self._helixClass)
    comment = property(lambda self: self._comment)
    length = property(lambda self: self._length)



def parse_pdb_HELIX_line(hlx_line):

    """Returns a dictionary of strings parsed from a HELIX line from a pdb 
    file

    Arguments:

    hlx_line: a string containing a HELIX line from a pdb file in string

    The keys and elements of the dictionary are
    'HELIX':       "HELIX  "
    'blank1':      A space
    'serNum':      Helix serial number
    'blank2':      A space
    'helixID':     Helix identifier
    'blank3':      A space
    'intResName':  Initial Residue name
    'blank4':      A space
    'initchainID': Chain identifier for chain containing this helix
    'blank5'       A space
    'initSeqNum':  Sequence number of the initial residue.
    'initICode':   Insertion code of the initial residue.
    'blank6':      A space
    'endResName':  Name of the terminal residue of the helix.
    'blank7':      A space
    'endChainID'   Chain identifier for the chain containing this helix.
    'blank8':      A space
    'endSeqNum':   Sequence number of the terminal residue.
    'endICode':    Insertion code of the terminal residue.
    'helixClass':  Helix class.
    'comment':     Comment about this helix.
    'blank9':      A space
    'length':      Length of this helix.

    Source:
    http://www.wwpdb.org/documentation/format33/sect5.html#HELIX

    """

    assert type(hlx_line) == str 

    parts = {}
    cleanstring = hlx_line.strip()
    # make sure string is long enough. If not, pad it on the right with ' ' 
    cleanstring += (76 - len(cleanstring)) * ' '

    parts['HELIX'] = cleanstring[:6]

    assert parts['HELIX'].strip() == 'HELIX'

    line = PDBHelixLine.parse_string(hlx_line)

    parts['blank1'] = ' '
    parts['serNum'] = line.serNum
    parts['blank2'] = ' '
    parts['helixID'] = line.helixID
    parts['blank3'] = ' '
    parts['intResName'] = line.intResName
    parts['blank4'] = ' '
    parts['initchainID'] = line.initChainID
    parts['blank5'] = ' '
    parts['initSeqNum'] = line.initSeqNum
    parts['initICode'] = line.initICode
    parts['blank6'] = ' '
    parts['endResName'] = line.endResName
    parts['blank7'] = ' '
    parts['endChainID'] = line.endChainID
    parts['blank8'] = ' '
    parts['endSeqNum'] = line.endSeqNum
    parts['endICode'] = line.endICode
    parts['helixClass'] = line.helixClass
    parts['comment'] = line.comment
    parts['blank9'] = ' '
    parts['length'] = line.length

    return parts


class PDBAtomLine(object):
    @classmethod
    def parse_string(cls, string):
        assert type(string) == str 

        cleanstring = string.strip()
        # make sure string is long enough. If not, pad it on the right with ' ' 
        cleanstring += (80 - len(cleanstring)) * ' '

        if not cleanstring[:6].strip() in ['ATOM', 'HETATM']:
            return None

        return PDBAtomLine(
            cleanstring[6:11],
            cleanstring[12:16],
            cleanstring[16],
            cleanstring[17:21].strip(), #4-letter residues
            cleanstring[21],
            cleanstring[22:26],
            cleanstring[26],
            cleanstring[30:38],
            cleanstring[38:46],
            cleanstring[46:54],
            cleanstring[54:60],
            cleanstring[60:66],
            cleanstring[76:78],
            cleanstring[78:],
            cleanstring[:6].strip()
        )
    
    def __init__(self, serial, name, altLoc, resName, chainID, resSeq, 
        iCode, x, y, z, occupancy, tempFactor, element, charge, kind='ATOM'):
        self._serial = serial.strip()
        self._name = name.strip()
        self._altLoc = altLoc.strip()
        self._resName = resName.strip()
        self._chainID = chainID.strip()
        self._resSeq = resSeq.strip()
        self._iCode = iCode.strip()
        self._x = x.strip()
        self._y = y.strip()
        self._z = z.strip()
        self._occupancy = occupancy.strip()
        self._tempFactor = tempFactor.strip()
        self._element = element.strip()
        self._charge = charge.strip()
        self._kind = kind

    @classmethod
    def for_dict(cls, dict):
        return PDBAtomLine(dict['serial'], dict['name'], dict['altLoc'], dict['resName'], 
            dict['chainID'], dict['resSeq'], dict['iCode'], dict['x'], dict['y'], dict['z'], 
            dict['occupancy'], dict['tempFactor'], dict['element'], dict['charge'], dict.get('ATOM', 'ATOM'))

    serial = property(lambda self: self._serial)
    name = property(lambda self: self._name)
    altLoc = property(lambda self: self._altLoc)
    resName = property(lambda self: self._resName)
    chainID = property(lambda self: self._chainID)
    resSeq = property(lambda self: self._resSeq)
    iCode = property(lambda self: self._iCode)
    x = property(lambda self: self._x)
    y = property(lambda self: self._y)
    z = property(lambda self: self._z)
    occupancy = property(lambda self: self._occupancy)
    tempFactor = property(lambda self: self._tempFactor)
    element = property(lambda self: self._element)
    charge = property(lambda self: self._charge)
    kind = property(lambda self: self._kind)

    def copy_with_serial(self, serial):
        return PDBAtomLine(serial, self._name, self._altLoc, self._resName, self._chainID, 
            self._resSeq, self._iCode, self._x, self._y, self._z, self._occupancy, 
            self._tempFactor, self._element, self._charge, self._kind)

    def copy_with_chainID(self, chainID):
        return PDBAtomLine(self._serial, self._name, self._altLoc, self._resName, chainID, 
            self._resSeq, self._iCode, self._x, self._y, self._z, self._occupancy, 
            self._tempFactor, self._element, self._charge, self._kind)

    def copy_with_name(self, name):
        return PDBAtomLine(self._serial, name, self._altLoc, self._resName, self._chainID, 
            self._resSeq, self._iCode, self._x, self._y, self._z, self._occupancy, 
            self._tempFactor, self._element, self._charge, self._kind)

    def copy_with(self, serial="<replace>", name="<replace>", altLoc="<replace>", resName="<replace>",
        chainID="<replace>", resSeq="<replace>", iCode="<replace>", x="<replace>", y="<replace>", z="<replace>",
        occupancy="<replace>", tempFactor="<replace>", element="<replace>", charge="<replace>",
        kind="<replace>", r=np.empty(0)):
        if len(r):
            assert(x == y and y == z and z == "<replace>")
            assert(len(r) == 3)
            xx, yy, zz = [f'{d:7.3f}' for d in r]
        else:
            xx, yy, zz = x, y, z
        return PDBAtomLine(self._serial if serial == "<replace>" else serial,
            self._name if name == "<replace>" else name,
            self._altLoc if altLoc == "<replace>" else altLoc,
            self._resName if resName == "<replace>" else resName,
            self._chainID if chainID == "<replace>" else chainID,
            self._resSeq if resSeq == "<replace>" else resSeq,
            self._iCode if iCode == "<replace>" else iCode,
            self._x if xx == "<replace>" else xx,
            self._y if yy == "<replace>" else yy,
            self._z if zz == "<replace>" else zz,
            self._occupancy if occupancy == "<replace>" else occupancy,
            self._tempFactor if tempFactor == "<replace>" else tempFactor,
            self._element if element == "<replace>" else element,
            self._charge if charge == "<replace>" else charge,
            self._kind if kind == "<replace>" else kind
        )

    def as_dict(self, atom='ATOM'):
        parts = {}
        parts['ATOM'] = self._kind
        parts['serial'] = self._serial
        parts['blank1'] = ' '
        parts['name'] = self._name
        parts['altLoc'] = self._altLoc
        parts['resName'] = self._resName
        parts['blank2'] = ' '
        parts['chainID'] = self._chainID
        parts['resSeq'] = self._resSeq
        parts['iCode'] = self._iCode
        parts['blank3'] = 3 * ' '
        parts['x'] = self._x
        parts['y'] = self._y
        parts['z'] = self._z
        parts['occupancy'] = self._occupancy
        parts['tempFactor'] = self._tempFactor
        parts['blank4'] = 10 * ' '
        parts['element'] = self._element
        parts['charge'] = self._charge
        return parts

    def print(self, atom='ATOM'): #TODO: give class an ATOM propery with default value and convert this to __string__ method
        return print_pdb_ATOM_line(self.as_dict(atom))

    def __str__(self):
        return print_pdb_ATOM_line(self.as_dict())

class PDBProtein(object):
    def __init__(self, atoms_lines, forces=np.array([])):
        does_have_forces = forces.size > 0
        if does_have_forces:
            assert(type(forces) == np.ndarray)
            assert(forces.shape == (len(atoms), 3))
        
        residues = {}
        for i, atom_line in enumerate(atoms_lines):
            atom = PDBAtom(atom_line, forces[i]) if does_have_forces else PDBAtom(atom_line)
            res_atoms = residues.get(atom.resSeq, [])
            res_atoms.append(atom)
            residues[atom.resSeq] = res_atoms
        self._residues = { k:PDBResidue(v) for (k, v) in residues.items()}
    
    residues = property(lambda self: list(self._residues.values()))
    residues_dict = property(lambda self: self._residues)

    def strip_hydrogens(self):
        for residue in self.residues:
            residue.strip_hydrogens()

    def rename_atoms(self, names_map):
        for residue in self.residues:
            residue.rename_atoms(names_map)

    def is_closer_than(self, d, residue):
        for this_residue in self.residues:
            if this_residue.is_closer_than(d, residue):
                return True
        return False

    def __str__(self):
        return '\n'.join(f"{residue_out}" for residue_out in [str(residue) for residue in self.residues])


@functools.total_ordering
class PDBResidue(object):
    def __init__(self, atoms, forces=np.array([])):
        assert(type(atoms[0]) == PDBAtom)
        if forces.size > 0:
            assert(type(forces) == np.ndarray)
            assert(forces.shape == (len(atoms), 3))
            for i,force in enumerate(forces):
                atoms[i].set_force(force)
        self._atoms = atoms

    atoms = property(lambda self: self._atoms)
    force = property(lambda self: np.sum(np.array([atom.force for atom in self._atoms]), axis=0))
    resSeq = property(lambda self: self._atoms[0].resSeq)
    resName = property(lambda self: self._atoms[0].resName)
    chainID = property(lambda self: self._atoms[0].chainID)

    def strip_hydrogens(self):
        self._atoms = [atom for atom in self._atoms if not atom.is_hydrogen()]

    def reorder_by_names(self, names):
        assert(len(names) == len(self._atoms))
        atom_by_name = {atom.name : atom for atom in self._atoms}
        try:
            self._atoms = [atom_by_name[name] for name in names]
        finally:
            return

    def rename_atoms(self, names_map):
        for atom in self._atoms:
            atom.set_name(names_map[atom.name])

    def is_closer_than(self, d, residue):
        for atom in self._atoms:
            if atom.is_hydrogen(): 
                continue
            for residue_atom in residue.atoms:
                if residue_atom.is_hydrogen():
                    continue
                if atom.distance_to(residue_atom) < d:
                    return True
        return False

    def __str__(self):
        return '\n'.join(f"{atom_out}" for atom_out in [str(atom) for atom in self._atoms])

    def __lt__(self, other):
        return int(self.resSeq) < int(other.resSeq)

    def __eq__(self, other):
        return self.resSeq == other.resSeq


class PDBAtom(object):
    def __init__(self, atom_line, force=np.zeros(3)):
        self._atom_line = atom_line
        self._force = force
        self._r = np.array([float(self._atom_line.x), float(self._atom_line.y), float(self._atom_line.z)])
    
    def set_force(self, force):
        assert(type(force) == np.ndarray)
        assert(force.dtype == np.float64)
        self._force = force
    
    name = property(lambda self: self._atom_line.name)
    force = property(lambda self: self._force, set_force)
    resSeq = property(lambda self: self._atom_line.resSeq)
    resName = property(lambda self: self._atom_line.resName)
    chainID = property(lambda self: self._atom_line.chainID)
    r = property(lambda self: self._r)
    serial = property(lambda self: self._atom_line.serial)

    def set_name(self, name):
        assert type(name) == str
        self._name = name
        self._atom_line = self._atom_line.copy_with_name(name)

    def is_hydrogen(self):
        return self._atom_line.element == 'H'

    def __str__(self):
        return str(self._atom_line)

    def distance_to(self, atom):
        return np.linalg.norm(self._r - atom.r)


def parse_pdb_ATOM_line(atm_line):

    """Returns a dictionary of strings parsed from an ATOM line from a pdb 
    file

    Arguments:

    atm_line: a string containing an ATOM line from a pdb file in string

    The keys and elements of the dictionary are
    'ATOM':       "ATOM  "
    'serial':     Atom serial number
    'blank1':     A space
    'name':       Atom name
    'altLoc':     Alternate location indicator
    'resName':    Residue name
    'blank2':     A space
    'chainID':    Chain identifier
    'resSeq':     Residue sequence number
    'iCode':      Code for insertion of residues
    'blank3':     3 spaces
    'x':          Orthogonal coordinates for X in Angstroms
    'y':          Orthogonal coordinates for Y in Angstroms
    'z':          Orthogonal coordinates for Z in Angstroms
    'occupancy':  Occupancy
    'tempFactor': Temperature factor
    'blank4':     10 spaces
    'element':    Element symbol, right-justified
    'charge':     Charge on the atom

    Source:
    http://www.wwpdb.org/documentation/file-format-content/format33/sect9.html#ATOM
    """

    assert type(atm_line) == str

    parts = {}
    cleanstring = atm_line.strip()
    cleanstring += (80 - len(cleanstring)) * ' '
    parts['ATOM'] = cleanstring[:6]

    atm = parts['ATOM'].strip()
    if atm not in ['ATOM', 'HETATM']:
        return None

    assert atm == 'ATOM' or atm == 'HETATM'

    line = PDBAtomLine.parse_string(atm_line)
    
    return line.as_dict(atm)


def print_pdb_ATOM_line(atm_dic):

    assert type(atm_dic) == dict

    atmLine = ''
    keys = ['ATOM','serial', 'blank1', 'name', 'altLoc', 'resName', 
                'chainID', 'resSeq', 'iCode', 'blank3', 'x', 'y', 'z', 
                'occupancy', 'tempFactor', 'blank4', 'element', 'charge'
            ]
    lengths = [6,5,1,4,1,4,1,4,1,3,8,8,8,6,6,10,2,2]
    justification = [0,1,0,0,0,0,0,1,0,0,1,1,1,1,1,0,1,1]
    for i, key in enumerate(keys):
        cur_val = atm_dic[key]
        if key == 'name' and cur_val[0] != ' ' and len(cur_val) < 4:
            cur_val = ' ' + cur_val
        if justification[i]:
            cur_val = (lengths[i] - len(cur_val)) * ' '  + cur_val
        else:
            cur_val += (lengths[i] - len(cur_val)) * ' '
        atmLine += cur_val
    return atmLine


def get_pdb_resolution_from_web(pdbid):
    url = f'https://www.rcsb.org/structure/{pdbid.upper()}'
    try:
        pdbWebPage_f = urlrequest.urlopen(url)
        pagecontents = pdbWebPage_f.read().decode('utf-8')
    except urlerror:
        print('error')
        resolution = 'N/F'
        return resolution
    soup = BeautifulSoup(pagecontents, features="lxml")
    tags = soup.find_all(id=re.compile('exp_header_.*_resolution'))
    resolution = (re.findall('\d+\.\d+', str(tags[0])) if len(tags) else ['N/F'])[0]
    if not len(resolution):
        resolution = 'N/A'
    return resolution


def pdb_rsln(pdbid):
    
    """Return resolution of pdbid

    Arguments: pdbid, a string with the four-character pdbid of the protein  
    for which you want the resoulution. This function seeks that resolution 
    on the protein databank website, rcsb.org
    
    """

    assert len(pdbid) > 3
    pdbid = pdbid.upper()[:4]
    db_path = os.path.expanduser("~/Proteins/misc/pdb_rsln.dat")
    (db_dir, db_file) = os.path.split(db_path)
    os.makedirs(db_dir, exist_ok=True)

    try:
        with open(db_path, 'r') as f:
            pdbids = array([])
            rslns = array([])
            for line in f:
                line_parts = line.split()
                pdbids = append(pdbids, line_parts[0])
                rslns = append(rslns, line_parts[1])
            if pdbid in pdbids:
                rsln = rslns[pdbids == pdbid][0]
                if rsln not in ['N/A', 'N/F']:
                    return rsln
    except IOError:
        pass
    rsln = '%.2f' % float(get_pdb_resolution_from_web(pdbid))
    with open(db_path, 'a') as f:
        f.write('%s %s\n' % (pdbid, rsln))
    return rsln
