"""
build configuration utilities
"""
import os, re

# defaults

_config = {
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
        'parallel': 'False',
        'version': 'release',
        'compiler': 'gcc',
        'coupling': 'False',
        },
    'MEMOSA': {
        'Build' : 'True',
        },
}

def config(x,y):
    # special hack for 'env'
    if y == 'env':
        new_env = []
        if  _config.has_key(x):
            for e in _config[x]:
                try:
                    new_env.append('%s=%s' % (re.compile(r'env\[(.*)\]').findall(e)[0], _config[x][e]))
                except IndexError:
                    pass
            # Keep compatibility with old format, for now.
            if  _config[x].has_key(y):            
                new_env.append(_config[x][y])
        return new_env
    try:
        # print "GET %s %s -> %s" % (x, y, _config[x][y])
        return _config[x][y]
    except KeyError:
        return ''


section = ''

def set_section(sec):
    global section
    section = sec
    if sec == 'ALL':
        print 50 * '-'
        print "WARNING: 'ALL' section is deprecated. Please see"
        print "https://memshub.org/infrastructure/memosa/wiki/Configuration"
        print 50 * '-'        
    return True

def set_value(val):
    global section
    #print "SET %s %s" %(section, val)
    if section == '':
        print "Error: No section set."
        return False

    if section == 'before':
        try: 
            _config['ALL']['before'].append(val)
        except KeyError: 
            _config['ALL'] = {'before':[val]}
        return True
    
    eq = val.find('=')
    if eq < 0:
        return False
    lval = val[:eq]
    rval = val[eq+1:]
    if lval == 'modules':
        print 50 * '-'        
        print "ERROR: 'modules' config variables has been removed. Modules"
        print "should all be loaded in the 'before' section. Please see"
        print "https://memshub.org/infrastructure/memosa/wiki/Configuration"
        print 50 * '-'          
        return False
    try: 
        _config[section][lval] = rval
    except KeyError: 
        _config[section] = {lval:rval}
    return True

def read(srcpath, cname):
    global _config

    if cname == '': return False

    filename = os.path.join(srcpath, cname)

    lnum = 0
    f = open(filename, 'r')
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
        print "Cannot parse line %s in %s: %s" % (lnum, filename, line)
        return False
    return  True
        
