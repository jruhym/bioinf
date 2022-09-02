import subprocess as subp

def add_group_to_index(groups, structure_file="md_0_1.tpr", input_index_file="index.ndx", output_index_file="index.ndx", should_backup=False):
    assert(type(groups) == dict)
    assert(len(groups) > 0)
    initial_last_dex = 0
    if input_index_file:
        with open(input_index_file) as f:
            for line in f:
                if line.startswith('[ '):
                    initial_last_dex += 1
        assert(initial_last_dex != 0)

    make_ndex_options = ['gmx', 'make_ndx', '-f', structure_file, '-o', output_index_file, '-backup', 'no'] + ([] if input_index_file == None else ['-n', input_index_file])
    for (name, group) in groups.items():
        make_ndx = subp.Popen(make_ndex_options, stdin=subp.PIPE, stdout=subp.PIPE)
        rename = f'\nname {initial_last_dex} {name}'
        make_ndx.communicate(f'{group}{rename}\nq\n'.encode()) #use q to close https://mailman-1.sys.kth.se/pipermail/gromacs.org_gmx-users/2011-July/063179.html
        initial_last_dex += 1