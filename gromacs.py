import subprocess as subp

def add_group_to_index(groups, structure_file="md_0_1.tpr", input_index_file="index.ndx", output_index_file="index.ndx", should_backup=False):
    assert(type(groups) == dict)
    assert(len(groups) > 0)
    assert(output_index_file != None)
    assert(output_index_file)

    initial_last_dex = 0
    do_backup = 'yes' if should_backup else 'no'

    make_ndex_options = ['gmx', 'make_ndx', '-f', structure_file, '-o', output_index_file, '-backup', do_backup]
    if input_index_file == None or not input_index_file:
        make_ndx = subp.Popen(make_ndex_options, stdin=subp.PIPE, stdout=subp.PIPE)
        make_ndx.communicate('\nq\n'.encode())
        input_index_file = output_index_file

    with open(input_index_file) as f:
        for line in f:
            if line.startswith('[ '):
                initial_last_dex += 1
    assert(initial_last_dex != 0)

    make_ndex_options += ['-n', input_index_file]
    for (name, group) in groups.items():
        make_ndx = subp.Popen(make_ndex_options, stdin=subp.PIPE, stdout=subp.PIPE)
        rename = f'\nname {initial_last_dex} {name}'
        make_ndx.communicate(f'{group}{rename}\nq\n'.encode()) #use q to close https://mailman-1.sys.kth.se/pipermail/gromacs.org_gmx-users/2011-July/063179.html
        initial_last_dex += 1
    return initial_last_dex
