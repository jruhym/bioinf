import unittest
import sure
from .funcs import *
from .pdb import *
from numpy.testing import assert_array_equal
from io import StringIO
 


class ProteinDatabaseHelixParser(unittest.TestCase):

    def test_should_return_the_proper_keys_and_values(self):
        parse_pdb_HELIX_line.when.called_with(
            'HELIX   1   1A ARG A 1125T ALA B 1134W ' + \
            '2THIS_IS_A_COMMENT_LIKE_IT_HUH?   10    '
        ).should.return_value(
            { 
                'comment': 'THIS_IS_A_COMMENT_LIKE_IT_HUH?', 
                'blank1': ' ', 
                'blank2': ' ', 
                'blank3': ' ', 
                'blank4': ' ', 
                'blank5': ' ', 
                'blank6': ' ', 
                'blank7': ' ', 
                'serNum': ' 1 ', 
                'blank9': ' ', 
                'helixClass': ' 2', 
                'endICode': 'W', 
                'HELIX': 'HELIX ', 
                'blank8': ' ', 
                'endSeqNum': '1134', 
                'initICode': 'T', 
                'endResName': 'ALA', 
                'helixID': ' 1A', 
                'initchainID': 'A', 
                'initSeqNum': '1125', 
                'intResName': 'ARG', 
                'length': '  10 ', 
                'endChainID': 'B'
            }        
        )

    def test_should_return_a_dict(self):
        parse_pdb_HELIX_line(
            'HELIX   1   1A ARG A 1125T ALA B 1134W ' + \
            '2THIS_IS_A_COMMENT_LIKE_IT_HUH?   10    '
            ).should.be.a(dict)



class ProteinDataBankAtomParser(unittest.TestCase):

    def test_should_return_the_correct_dictionary(self):
        parse_pdb_ATOM_line.when.called_with(
            'ATOM     32  N2 AARG A  -3B     11.281  86.699  ' + \
            '94.383  0.50 35.88           N-12'
            ).should.return_value(
            {'ATOM': 'ATOM  ',
             'altLoc': 'A',
             'blank1': ' ',
             'blank2': ' ',
             'blank3': '   ',
             'blank4': '          ',
             'chainID': 'A',
             'charge': '-12',
             'element': 'N',
             'iCode': 'B',
             'name': 'N2',
             'occupancy': '0.50',
             'resName': 'ARG',
             'resSeq': '-3',
             'serial': '32',
             'tempFactor': '35.88',
             'x': '11.281',
             'y': '86.699',
             'z': '94.383'}
             )



class TestPrintAtomLine(unittest.TestCase):

    def test_should_return_proper_line(self):
        print_pdb_ATOM_line.when.called_with(
            {'ATOM': 'ATOM  ', 
            'altLoc': 'A', 
            'blank1': ' ', 
            'blank2': ' ', 
            'blank3': '   ', 
            'blank4': '          ', 
            'chainID': 'A', 
            'charge': '-12', 
            'element': ' N', 
            'iCode': 'B', 
            'name': ' N2 ', 
            'occupancy': '  0.50', 
            'resName': 'ARG', 
            'resSeq': '  -3', 
            'serial': '   32', 
            'tempFactor': ' 35.88', 
            'x': '  11.281', 
            'y': '  86.699', 
            'z': '  94.383'}
            ).should.return_value(
            'ATOM     32  N2 AARG A  -3B     ' +\
            '11.281  86.699  94.383  0.50 35.88           N-12'
            )

    def test_name_should_begin_with_space(self):
        print_pdb_ATOM_line.when.called_with(
            {'ATOM': 'ATOM  ', 
            'altLoc': 'A', 
            'blank1': ' ', 
            'blank2': ' ', 
            'blank3': '   ', 
            'blank4': '          ', 
            'chainID': 'A', 
            'charge': '-12', 
            'element': ' N', 
            'iCode': 'B',
            # name doesn't start with space
            'name': 'N2 ', 
            'occupancy': '  0.50', 
            'resName': 'ARG', 
            'resSeq': '  -3', 
            'serial': '   32', 
            'tempFactor': ' 35.88', 
            'x': '  11.281', 
            'y': '  86.699', 
            'z': '  94.383'}
            ).should.return_value(
            'ATOM     32  N2 AARG A  -3B     ' +\
            '11.281  86.699  94.383  0.50 35.88           N-12'
            )


        
class TestProtCntctParser(unittest.TestCase):

    def test_should_return_correct_dictionary(self):
        parse_prot_cntct_line.when.called_with('1     HB    1:1ORSC_c   ' + \
            '8 ALA34.O       1:1ORSC_c  12 SER38.OG       3'
            ).should.return_value(
            {'atom1': 'O',
             'atom2': 'OG',
             'chain1': '1:1ORSC_c',
             'chain2': '1:1ORSC_c',
             'intrcn_typ': 'HB',
             'line_num': '1',
             'moe_dex1': '8',
             'moe_dex2': '12',
             'net': '3',
             'pdb_res_dex1': '34',
             'pdb_res_dex2': '38',
             'res_typ1': 'ALA',
             'res_typ2': 'SER'}
        )       



class TestFloatList(unittest.TestCase):

    def test_correct_floats_from_list_o_str_o_floats(self):
        assert_array_equal([2.0, 3.0, 7.0], float_list(['2', '3E0', '7.0']))

    def test_correct_float_if_passed_a_single_str(self):
        assert_array_equal(float_list('5'), array([5.]))
    
    def test_correct_float_if_passed_a_single_nonstr_item(self):
        assert_array_equal(float_list(True), array([1.]))

    def test_correct_abs_if_item_is_complex(self):
        assert_array_equal(float_list(complex(1,1)), array(abs(complex(1,1))))



class TestGetRanges(unittest.TestCase):

    def test_should_rtrn_corct_flatnd_ary_o_ints_frm_lst_o_str_o_rngs(self):
        assert_array_equal(((1,2,3),(4,5,6),(7,8,9)), 
            get_ranges(['-range1:3,4:6,7:9'], flatten=False))
        
    def test_should_return_correct_matrix_from_list_of_strs_of_ranges(self):
        assert_array_equal((1,2,3,4,5,6,7,8,9), 
            get_ranges(['-range1:3,4:6,7:9']))
        
    def test_different_option_name_should_work_for_ranges(self):
        assert_array_equal((1,2,3,4,5,6,7,8,9), 
            get_ranges(['-carlos1:3,4:6,7:9'], key='-carlos'))
        


class TestUnionOfDistrosGivenX(unittest.TestCase):
 
    def setUp(self):
        self.unionOfXsToTest = array([1,2,3,4])
        self.distroXToTest = [2,3]
        self.distroYToTest = [8,9]
        self.unionizedYDistro = unionize_distro_given_x(self.unionOfXsToTest, 
            self.distroXToTest, self.distroYToTest)
        self.correctUnionOfYs = array([0,8,9,0])

    def test_gives_correct_y_union(self):
        self.assertTrue((self.unionizedYDistro == self.correctUnionOfYs).
            all())



class UnionOfTwoDistros(unittest.TestCase):

    def setUp(self):
        self.distro_1_X_ToTest = [0, 1, 2, 3]
        self.distro_1_Y_ToTest = [4, 5, 6, 7]
        self.distro_2_X_ToTest = [3, 4, 5, 6]
        self.distro_2_Y_ToTest = [6, 7, 8, 9]
        self.unionizedX, self.unionizedY1, self.unionizedY2 = \
            unionize_2_distros(self.distro_1_X_ToTest, self.distro_1_Y_ToTest,
                self.distro_2_X_ToTest, self.distro_2_Y_ToTest)
        self.correctUnionX = array([0., 1., 2., 3., 4., 5., 6.])
        self.correctUnionY1 = array([4., 5., 6., 7., 0., 0., 0.])
        self.correctUnionY2 = array([0., 0., 0., 6., 7., 8., 9.])

    def test_gives_correct_x_union(self):
        self.assertTrue((self.unionizedX == self.correctUnionX).all())

    def test_gives_correct_y1_union(self):
        self.assertTrue((self.unionizedY1 == self.correctUnionY1).all())

    def test_gives_correct_y2_union(self):
        self.assertTrue((self.unionizedY2 == self.correctUnionY2).all())

 

class ProteinDataBankResolutionSpec(unittest.TestCase):

    def test_should_resolve_from_the_web(self):
        get_pdb_resolution_from_web.when.called_with(
            '2rh1').should.return_value('2.4')

    def test_should_give_not_found_if_pdbID_does_not_exist(self):
        get_pdb_resolution_from_web.when.called_with(
            '1234').should.return_value('N/F')

    def test_should_resolve_from_the_web_or_from_the_cache(self):
        pdb_rsln.when.called_with('2rh1').should.return_value('2.40')


