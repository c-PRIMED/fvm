"""
build configuration utilities
"""
import os

# defaults

_config_pkgs = {}

_config_srcs = {
    'Testing': {
        'queue': 'standby',
        'cpus' : '8',
        'nodes': '1',
        'walltime': '1:00:00',
        },
    'lammps': {
        'build': 'openmpi',
        },
    'fvm': {
        'parallel': False,
        'version': 'debug',
        'compiler': 'gcc',
        },
}

def config(x,y):
    try:
        # print "GET %s %s -> %s" % (x, y, _config[x][y])
        return _config[x][y]
    except:
        return ''


section = ''

def set_section(sec):
    global section
    section = sec
    return True

def set_value(val):
    global section
    # print "SET %s %s" %(section, val)
    if section == '':
        print "Error: No section set."
        return False

    eq = val.find('=')
    if eq < 0:
        return False
    lval = val[:eq]
    rval = val[eq+1:]           
    try: 
        _config[section][lval] = rval
    except KeyError: 
        _config[section] = {lval:rval}
    return True

def read(srcpath, filename, sources, packages):
    global _config

    if filename == '': return False
    if sources:
        _config = _config_srcs.copy()
    else:
        _config = _config_pkgs.copy()

    filename = os.path.join(srcpath, "config", filename)
    filenames = []
    if sources:
        filenames.append(filename)
    if packages:
        filenames.append(filename+'-pkgs')

    for fn in filenames:
        lnum = 0
        f = open(fn, 'r')
        while True:
            line = f.readline()
            lnum = lnum + 1
            if not line:
                break
            line = line.rstrip()
            if line == '' or line[0] in '#;':
                continue
            if line[-1] == ':' and set_section(line[:-1]):
                continue
            if line[0].isspace():
                line = line.lstrip()
                if line[0] in '#;':
                    continue
                if set_value(line):
                    continue
            print "Cannot parse line %s in %s: %s" % (lnum, fn, line)
            return False
    return  True
        
