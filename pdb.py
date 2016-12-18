import urllib2
from HTMLParser import HTMLParser
# import formatter
import os
import os.path
from numpy import array, append


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

        assert cleanstring[:6].strip() == 'ATOM'

        return PDBAtomLine(
            cleanstring[6:11],
            cleanstring[12:16],
            cleanstring[16],
            cleanstring[17:20],
            cleanstring[21],
            cleanstring[22:26],
            cleanstring[26],
            cleanstring[30:38],
            cleanstring[38:46],
            cleanstring[46:54],
            cleanstring[54:60],
            cleanstring[60:66],
            cleanstring[76:78],
            cleanstring[78:]
        )
    def __init__(self, serial, name, altLoc, resName, chainID, resSeq, 
        iCode, x, y, z, occupancy, tempFactor, element, charge):
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
    http://www.wwpdb.org/documentation/format32/sect9.html#ATOM

    """

    assert type(atm_line) == str

    parts = {}
    cleanstring = atm_line.strip()
    #short = 80 - len(cleanstring)
    #make sure the string is long enough. If not, pad it on the right with ' '
    #if short > 0:
    #    for i in range(short):
    #        cleanstring += ' '
    cleanstring += (80 - len(cleanstring)) * ' '
    parts['ATOM'] = cleanstring[:6]
    assert parts['ATOM'].strip() == 'ATOM'

    line = PDBAtomLine.parse_string(atm_line)

    parts['serial'] = line.serial
    parts['blank1'] = ' '
    parts['name'] = line.name
    parts['altLoc'] = line.altLoc
    parts['resName'] = line.resName
    parts['blank2'] = ' '
    parts['chainID'] = line.chainID
    parts['resSeq'] = line.resSeq
    parts['iCode'] = line.iCode
    parts['blank3'] = 3 * ' '
    parts['x'] = line.x
    parts['y'] = line.y
    parts['z'] = line.z
    parts['occupancy'] = line.occupancy
    parts['tempFactor'] = line.tempFactor
    parts['blank4'] = 10 * ' '
    parts['element'] = line.element
    parts['charge'] = line.charge
    return parts



def print_pdb_ATOM_line(atm_dic):

    assert type(atm_dic) == dict

    atmLine = ''
    keys = ['ATOM','serial', 'blank1', 'name', 'altLoc', 'resName', 'blank2', 
               'chainID', 'resSeq', 'iCode', 'blank3', 'x', 'y', 'z', 
                'occupancy', 'tempFactor', 'blank4', 'element', 'charge']
    lengths = [6,5,1,4,1,3,1,1,4,1,3,8,8,8,6,6,10,2,2]
    justification = [0,1,0,0,0,0,0,0,1,0,0,1,1,1,1,1,0,1,1]
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



# This class will be used by the pdbid_rsln
class PDBParser(HTMLParser):
    def __init__(self):

        """Inherited from HTMLParser.HTMLParser + init of outwithit and 
        resolution.

        """

        HTMLParser.__init__(self)
        self.outwithit = False
        self.resolution = ''
 
    def handle_data(self, data):
        if data == "Resolution":
            self.outwithit = True
        elif self.outwithit:
            self.resolution = data
            self.outwithit = False

def get_pdb_resolution_from_web(pdbid):

    url = 'http://www.pdb.org/pdb/explore/explore.do?structureId=%s' % \
        (pdbid.upper())
    pdbWebPage_f = urllib2.urlopen(url)
    pagecontents = pdbWebPage_f.read()
    pdbWebPage_f.close()
    if 'No results were found matching your query' in pagecontents:
        resolution = 'N/F'
    else:
        parse = PDBParser()
        parse.feed(pagecontents)
        resolution = parse.resolution.split(":")[-1].strip()
    if not len(resolution):
        resolution = 'N/A'
    return resolution


def pdb_rsln(pdbid):
    
    """Return resolution of pdbid

    Arguments: pdbid, a string with the four-character pdbid of the protein  
    for which you want the resoulution. This function seeks that resolution 
    on the protein databank website, pdb.org
    
    """

    assert len(pdbid) > 3
    pdbid = pdbid.upper()[:4]

    # access local db of resolution if available
    db_path = os.path.expanduser('~/Proteins/misc/pdb_rsln.dat')
    db_exists = os.path.exists(db_path)
    db_dir, db_file = os.path.split(db_path)

    # is it available
    if db_exists:
        # read in its contents
        f = open(db_path, 'r')
        contents = f.read()
        f.close()
        contentlines = contents.splitlines()
        pdbids = array([])
        rslns = array([])
        for line in contentlines:
            line_parts = line.split()
            pdbids = append(pdbids, line_parts[0])
            rslns = append(rslns, line_parts[1])
        # is pdbid of interest already in local db
        if pdbid in pdbids:
            rsln = rslns[pdbids == pdbid][0]
            if rsln not in ['N/A', 'N/F']:
                return rsln
    # if not
    else:
        # does the directory not exist
        if not os.path.exists(db_dir):
            os.makedirs(db_dir)

    get_pdb_resolution_from_web(pdbid)

    f = open(db_path, 'a')
    f.write('%s %s\n' % (pdbid, rsln))
    f.close()
    return rsln

    