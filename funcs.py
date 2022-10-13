# Library for bioinformatics stuff
# University of the Sciences 2012
# All rights reserved. To be used for academic purposes only.
# Andrew J. Heim <lastfirstMI at gmail>

'''This is a library for bioinformatics stuff.

imports sys, glob, os, os.path, urllib2, formatter, * from numpy, and htmllib
from HTMLParser

List of functions:

parse_pdb_HELIX_line
parse_pdb_ATOM_line
print_pdb_ATOM_line
parse_prot_cntct_line
float_list
get_ranges
intrcn_distro
unionize_distro_given_x
unionize_2_distros
pdb_rsln

Classes:
PDBParser

'''

import sys
import glob 
import os
import os.path
from numpy import array, where, append, unique, zeros, ndarray, empty, abs, \
    ndarray, pi, e
import formatter

from .constants import AminoAcid
   


def parse_prot_cntct_line(prot_cntct_line):

    """Returns a dictionary containing the info from a line of a protein
    contact report from MOE

    Usage:
    
    parse_prot_cntct_line(prot_cntct_line)
    where pc_line is a string from a protein contact line from a MOE 
    protein contact report.

    Returns a dictionary of the parts of the line, the keys of which are

    'line_num':        The contact line number
    'intrcn_typ':      The type of residue--residue interaction
    'chain1'/'2':      The chain of the first/second residue in the contact
    'pdb_res_dex1'/'2' The pdb index of the first/second residue
    'moe_dex1'/'2':    The MOE-assigned index of the first/second residue
    'atom1'/'2':       The atom designation w/n residue of 1st/2nd residue
    'res_typ1'/'2':    The type of residue 
    'net':             To which network the contact belongs

    """

    assert type(prot_cntct_line) == str
    parts = {}
    line_parts = prot_cntct_line.split()
    parts['line_num'] = line_parts[0]
    parts['intrcn_typ'] = line_parts[1]
    parts['chain1'] = line_parts[2]
    parts['moe_dex1'] = line_parts[3]
    contact_point = line_parts[4]
    res, parts['atom1'] = contact_point.split('.')
    parts['res_typ1'] = res[:3]
    parts['pdb_res_dex1'] = res[3:]
    parts['chain2'] = line_parts[5]
    parts['moe_dex2'] = line_parts[6]
    contact_point = line_parts[7]
    res, parts['atom2'] = contact_point.split('.')
    parts['res_typ2'] = res[:3]
    parts['pdb_res_dex2'] = res[3:]
    parts['net'] = line_parts[8]
    assert type(parts) == dict
    return parts



def float_list(listFloats):

    """Return a numpy array of floats from a list of strings of floats

    """
    if isinstance(listFloats, (str)):
        templistFloats = [listFloats]
        listFloats = templistFloats
    try:
        list_len = len(listFloats)
    except TypeError:
        listFloats = [listFloats]
        list_len = len(listFloats)
    list_blank = empty(list_len, dtype=float)
    for i in range(list_len):
        try:
            if isinstance(listFloats[i], complex):
                listFloats[i] = abs(listFloats[i])
            if isinstance(listFloats[i], str):
                if listFloats[i].lower().strip() == 'pi':
                    listFloats[i] = pi
                elif listFloats[i].lower().strip() == 'e':
                    listFloats[i] = e
            list_blank[i] = float(listFloats[i])
        except ValueError:
            raise ValueError(listFloats[i] + ' is not float-able.')
    return list_blank



def get_ranges(args, key='-range', flatten=True):

    """Return array of ints specified by form 'keyx1:y1[,x2:y2...]'.

    The resulting array looks like (x1,...y1[,x2...y2...]).

    Arguments:
    
    args: list of string(s) at least one of which in the form 
    keyx1:y1[,x2:y2...]
    
    key: string of the first part of args (default '-range')
    
    flatten: if True, returns a flattened array of the indecies contained
    in the ranges. Else returns a list of lists. The latter can have sublists
    of varying length. The former can use 'in' to test for the presence of a
    single element.

    The point is to determine ranges of residue indecies to be included 
    or excluded by other functions which take arguments of this form.

    """
    assert isinstance(args, (list, tuple, ndarray))
    assert isinstance(flatten, bool) 

    if flatten:
        vld_dex = array([], dtype=int) 
    else:
        vld_dex = []
    for i, arg in enumerate(args):
        if arg.startswith(key):
            assert isinstance(arg, str)
            option = args.pop(i)
            ranges = option[len(key):].split(',')
            for cur_rng in ranges:
                beg,end = cur_rng.split(':')  
                try:
                    a = [int(beg), int(end)]
                    a.sort()
                    ibeg, iend = a
                except ValueError:
                    sys.exit('Part of ' + cur_rng + ' != int')
                if flatten:
                    vld_dex = append(vld_dex, range(ibeg, iend + 1))
                else:
                    vld_dex.append(range(ibeg, iend + 1))
    if flatten:
        return vld_dex.flatten()
    else:
        return vld_dex



class ContactEnd(object):
    
    @classmethod
    def parse_string(cls, string):

        # A MOE contact report line looks like this
        #     Type   Chain     Pos Residue       Chain     Pos Residue      Net
        #1     HB    1:1C3WA_c   3 LEU13.O       1:1C3WA_c   7 THR17.OG1     16

        # First split off the bonding atoms in each residue, 
        # deliminated by a '.'
        parts = string.split('.')

        # this next line replaces the first element (a string) with a 
        # 2 elements two parts of a string separated at the third position 
        # (residues have 3-letter abbrs.)
        parts[:1] = [parts[0][:3], parts[0][3:]]
        return ContactEnd(parts[0], int(parts[1]), parts[2], )

    def __init__(self, residue1_name, UID1, atom1_name):
        assert isinstance(residue_name, str)
        assert isinstance(UID, int)
        assert isinstance(atom_name, str)

        self._residue_name = residue_name
        self._UID = UID
        self._atom_name = atom_name

    UID = property(lambda self: self._UID)
    atom_name = property(lambda self: self._atom_name)
    residue_name = property(lambda self: self._residue_name)
        
    def find_name_in_amino_acids(self, aminoacids):
        aadex = where(aminoacids[:, 1] == self.residue_name)
        if len(aadex[0]) == 0:
            raise LookupError('%s %s of %s is not in list' % (
                self.residue_name, 
                residue1UID, 
                report
            ))
        return aminoacids[aadex, 3][0, 0]



class InterresidueContact(object):

    def __init__(self, contact_line):
        assert isinstance(contact_line, str)
        parts = contact_line.split()
        self._contact_dex = int(parts[0])
        self._type = parts[1]
        res_and_UID_and_atom = parts[4]
        
        res_and_UID, self._atom1 = res_and_UID_and_atom.split('.')
        self._residue1_abbr = res_and_UID[:3]
        self._residue1_group = AminoAcid.by_abbr(self._residue1_abbr).group
        self._residue1_UID = int(res_and_UID[3:])
        
        res_and_UID_and_atom = parts[7]

        res_and_UID, self._atom2 = res_and_UID_and_atom.split('.')
        self._residue2_abbr = res_and_UID[:3]
        self._residue2_group = AminoAcid.by_abbr(self._residue2_abbr).group
        self._residue2_UID = int(res_and_UID[3:])
        interaction_pair = [self._residue1_group, self._residue2_group]
        interaction_pair.sort()
        self._interaction_pair_string = '%s %s %s' % (
            interaction_pair[0], 
            interaction_pair[1], self._type
        )

    atom1 = property(lambda self: self._atom1)
    residue1_UID = property(lambda self: self._residue1_UID)
    residue1_abbr = property(lambda self: self._residue1_abbr)
    residue1_group = property(lambda self: self._residue1_group)
    atom2 = property(lambda self: self._atom2)
    residue2_UID = property(lambda self: self._residue2_UID)
    residue2_name = property(lambda self: self._residue2_name)
    residue2_group = property(lambda self: self._residue2_group)
    contact_type = property(lambda self: self._type)
    interaction_pair_string = property(
        lambda self: self._interaction_pair_string
    )



def unionize_distro_given_x(cmbn_x, x, y):

    """Return union of y's given union of x's

    Arguments:
    cmbn_x: a numpy array of union of x's from first and second sets
    x: array/list of x values of first set
    y: array/list of y values of first set

    """

    assert(type(cmbn_x) == ndarray)
    assert(len(x) == len(y))

    # make a new y lists as long as the combination 
    lngth = len(cmbn_x)
    cmbn_y = zeros(lngth)

    # copy the y's to combination  array 
    # at the loctions at which the array of x corresponds to the 
    # array of the combination of x.
    # In other words copy the y's from the original set to y's of the 
    # combination set at the proper position as determined by comparing the 
    # original x positions to the combined x positions
    for i, j in enumerate(x):
        bool_ar = (j == cmbn_x)
        cmbn_y[bool_ar] = y[i]

    return cmbn_y
 
# Unionize the x axis of two distros. Inputs are 4 arrays. Output is three 
# arrays the new x, y1, and y2

def unionize_2_distros(x1, y1, x2, y2):

    """Return union of abscissae from two sets and corresponding ordinates
    
    Arguments:
    
    x1: abscissa of set 1
    y1: ordinate of set 1
    x2: abscissa of set 2
    y2: ordinate of set 2 

    """

    assert(len(x1) == len(y1))
    assert(len(x2) == len(y2))

    # Make a union of unique elements between the x1 and x2 sets
    cmbn_x = unique(append(x1, x2))

    cmbn_y1 = unionize_distro_given_x(cmbn_x, x1, y1)
    cmbn_y2 = unionize_distro_given_x(cmbn_x, x2, y2)
    
    return cmbn_x, cmbn_y1, cmbn_y2


