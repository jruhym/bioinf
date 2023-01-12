import subprocess as subp
import asyncio
import signal

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

async def run_async_shell_cmd(cmd, inputs=[], do_stdout=False, do_stderr=False, requires_kill=False):
    assert(type(inputs)==list)
    proc = await asyncio.create_subprocess_shell(cmd.encode(), stdin=asyncio.subprocess.PIPE if inputs else None, stdout=asyncio.subprocess.PIPE, stderr=asyncio.subprocess.PIPE)
    stdout, stderr = await proc.communicate() if not inputs else await proc.communicate("\n".join(inputs).encode())
    print(f'[{cmd!r} exited with {proc.returncode}]')
    if do_stdout and stdout:
        print(f'[stdout]\n{stdout.decode()}')
    if do_stderr and stderr:
        print(f'[stderr]\n{stderr.decode()}')
    if requires_kill:
        print("Killing")
        stdin = proc.stdin
        if stdin.can_write_eof():
            print("Writing EOF")
            stdin.write_eof()
        stdin.close()
        await stdin.wait_closed()
        #proc.terminate()
    return
