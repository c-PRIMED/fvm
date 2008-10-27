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
    'fvm': {
        'Build': '1',
        },
    'lammps': {
        'Build': '1',
        'build': 'openmpi',
        },
}

def config(x,y):
    try:
        #print "GET %s %s -> %s" % (x, y, _config[x][y])
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
    #print "SET %s %s" %(section, val)
    if section == '':
        print "Error: No section set."
        return False
    if section == 'before' or section == 'after':
        try:
            _config[section][0].append(val)
        except:
            _config[section] = {0: [val]}
    else:
        val = val.split('=')
        if val[0] == 'skip' and len(val) == 1:
            val.append(1)
        if len(val) != 2:
            return False
        if val[0] == 'before' or val[0] == 'after':
            try:
                _config[section][val[0]].append(val[1])
            except:
                try: 
                    _config[section][val[0]] = [val[1]]
                except KeyError: 
                    _config[section] = {val[0]:[val[1]]}
        else:
            try: 
                _config[section][val[0]] = val[1]
            except KeyError: 
                _config[section] = {val[0]:val[1]}
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
            if line[0].isspace() and set_value(line.lstrip()):
                continue
            print "Cannot parse line %s in %s: %s" % (lnum, fn, line)
            return False
    return  True
        
