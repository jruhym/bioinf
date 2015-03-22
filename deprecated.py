def what_res_pos(stringar, file_path, index=False):

    """Deprecated 3/1/13. Use parse_pdb_ATOM_line
    
    Return string of pdb index.

    The pdb files don't have a uniform number of whitespaces before the pdb 
    index. This function searches for and returns the index if the line has 
    been split into a list of strings.

    Arguments: 

    stringar: a list/array of strings parsed from an ATOM line 
    in a pdb file

    file_path: a string of the file path of that pdb file

    index: a bool. If true, return index in array at which the string refers
    to the pdb index. Else, return the string of the pdb index 
    (default = False).

    """

    if len(stringar) > 5:
        for l in range(3,6,1):
            if stringar[l].isdigit():
                if index == True:
                    return l
                return stringar[l]
            elif len(stringar[l]) > 1:
                if stringar[l][1:].isdigit():
                    if index == True:
                        print 'Warning: requesting index of the res position, '
                        print 'but value at that index starts with ' + \
                            'a nondigit'
                        return l
                    else: 
                        return stringar[l][1:] 
        dex = 'fidget'
        weregood = True
        while (weregood):
            print stringar
            if index:
                print 'Warning: requesting index of the res position, ' + \
                    'but value at that index might not be a digit.'
                dex = raw_input('What is the idex of the position of the ' + \
                                    'residue in ' + file_path + '\n')
            else:
                dex = raw_input('What is the position of the ' + \
                                    'residue in ' + file_path + '\n')
    
            if dex.isdigit():
                if index: 
                    indx = int(dex)
                    if indx < len(stringar):
                        print 'Returning %i' % indx + ', which points to ' + \
                            stringar[indx]
                        return indx
                    else: 
                        print 'Use an index < ', len(stringar), '\n'
                else:
                    print 'Returning ' + dex + ' as the position of res.'
                    return dex
            else:
                print 'Use an int.\n'
