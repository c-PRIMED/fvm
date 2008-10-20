"""
build utility functions.
"""
import sys, os, shutil
from config import config

colors = {
    'BOLD'  :'\033[1m',
    'RED'   :'\033[1;31m',
    'GREEN' :'\033[1;32m',
    'YELLOW':'\033[1;33m',
    'PINK'  :'\033[1;35m',
    'BLUE'  :'\033[1;34m',
    'CYAN'  :'\033[1;36m',
    'DRED'   :'\033[0;31m',
    'DGREEN' :'\033[0;32m',
    'DYELLOW':'\033[0;33m',
    'DPINK'  :'\033[0;35m',
    'DBLUE'  :'\033[0;34m',
    'DCYAN'  :'\033[0;36m',
    'NORMAL':'\033[0m'
    }
"colors used for printing messages"


maxlen = 20
verbose = False

def _reset_types():
    global g_type
    g_type = {
        'CONF':"%s%s%s"%(colors['BOLD'],"Configuring",colors['NORMAL']),
        'BUILD':"%s%s%s"%(colors['BOLD'],"Building",colors['NORMAL']),
        'INSTALL':"%s%s%s"%(colors['BOLD'],"Installing",colors['NORMAL']),
        'TEST':"%s%s%s"%(colors['BOLD'],"Testing",colors['NORMAL'])
        }

def clear_colors():
    for k in colors.keys():
        colors[k]=''
    _reset_types()

        
_reset_types()
if ('NOCOLOR' in os.environ) \
        or (os.environ.get('TERM', 'dumb') in ['dumb', 'emacs']) \
        or (not sys.stdout.isatty()):
    clear_colors()
    
def cprint(col, str):
    try: mycol=colors[col]
    except KeyError: mycol=''
    print "%s%s%s" % (mycol, str, colors['NORMAL'])

def _niceprint(msg, type=''):
    def print_pat(color):
        print '\n%s: %s%s%s' % (type, colors[color], msg, colors['NORMAL'])
    if type == 'ERROR' or type == 'WARNING':
        print_pat('RED')
    elif type=='DEBUG':
        print_pat('DCYAN')
    else:
        print_pat('NORMAL')
                      
def debug(msg):
    global verbose
    if not verbose:
        return
    _niceprint(msg, 'DEBUG')

def warning(msg):
    _niceprint(msg, 'WARNING')

def error(msg):
    _niceprint(msg, 'ERROR')

def fatal(msg, ret=1, trace=1):
    import traceback
    cprint('RED', '%s' % msg)
    if trace: traceback.print_stack()
    sys.exit(ret)

def pmess(type, pkg, dir):
    "print a build message"
    sr = '%s %s%s%s in %s' % (g_type[type], colors['DCYAN'], pkg, colors['NORMAL'], dir)
    global maxlen
    maxlen = max(maxlen, len(sr))
    print "%s :" % sr.ljust(maxlen),
    sys.stdout.flush()

def remove_file(name):
    try:
        os.remove(name)
    except OSError:
        pass

def do_env(c):
    a = c.split('=')
    if len(a) == 1 or a[1] == '':
        debug("unset "+a[0])
        os.unsetenv(a[0])
    elif len(a) > 2:
        fatal("Cannot parse command: env",c)
    else:
        debug("%s=%s" % (a[0], a[1]))
        os.putenv(a[0],a[1])

def fix_path(k, v, prepend, unload):
    if unload:
        s1 = "remove"
        s2 = "from"
    else:
        s2 = "to"
        if prepend: s1 = "prepend"
        else: s1 = "append"

    debug("%s '%s' %s %s" % (s1, v, s2, k))

    try:
        e = os.environ[k]
    except KeyError: 
        e = ''
        
    if unload:
        os.environ[k] = e.replace(v,'').strip(':')
    else:
        if prepend:
            os.environ[k] = (v + ':' + e).strip(':')
        else:
            os.environ[k] = (e + ':' + v).strip(':')
    debug("%s=%s" % (k, os.environ[k]))


# python handler for Environment Modules
# http://modules.sourceforge.net/
# We can't just load modules with system()
# or popen() because those run in subshells and do not
# modify the current environment.
def module_load(m, unload=False):
    if unload:
        debug("module_unload: "+m)
    else:
        debug("module_load: "+m)
        
    for line in os.popen("/bin/bash -l -c 'module show %s 2>&1'" % m):
        x = line.split()
        if not x: continue
        if x[0] == 'setenv':
            val = ' '.join(x[2:])
            if unload:
                debug("unsetenv %s" % x[1])
                os.unsetenv(x[1])
            else:
                debug("%s=%s" % (x[1], val))
                os.putenv(x[1],val)
        elif x[0] == 'prepend-path':
            val = ' '.join(x[2:])
            fix_path(x[1], val, 1, unload)
        elif x[0] == 'append-path':
            val = ' '.join(x[2:])
            fix_path(x[1], val, 0, unload)
        elif x[0] == 'module-whatis':
            val = ' '.join(x[1:])
            debug ("[un]loading %s: %s" % (m, val))
        elif x[0] == 'module':
            if x[1] == 'load':
                val = ' '.join(x[2:])
                debug("need to load module " + val)
                module_load(val, unload)
            else:
                warning ("Unknown module command " + x[1])
        else:
            if x[0].find('ERROR') > 0:
                fatal("ERROR loading module '%s'" %  m, -1, 0)

def run_commands(section, pkg):
    if pkg == 0:
        c = config(section, 0)
    else:
        c = config(pkg, section)
    if c:
        for cmd in c:
            cargs = cmd.split(' ')
            if cargs[0] == "module":
                if cargs[1] == "load":
                    for m in cargs[2:]:
                        module_load(m)
                elif cargs[1] == "unload":
                    for m in cargs[2:]:
                        module_load(m, 1)
            elif cargs[0] == 'env':
                do_env(''.join(cargs[1:]))
            else:
                debug("executing: "+cmd)
                s = os.system(cmd)
                if s:
                    error("While running 'before' commands. Execution of")
                    print cmd
                    print "failed."
                    sys.exit(-1)

def copytree(src, dst, ctype):
    names = os.listdir(src)
    os.makedirs(dst)
    if ctype == 0:
        return
    errors = []
    for name in names:
        srcname = os.path.join(src, name)
        dstname = os.path.join(dst, name)
        try:
            if os.path.islink(srcname):
                linkto = os.readlink(srcname)
                os.symlink(linkto, dstname)
            elif os.path.isdir(srcname):
                copytree(srcname, dstname, ctype)
            elif ctype == 2:
                shutil.copy2(srcname, dstname)
            else:
                os.symlink(srcname, dstname)
        except (IOError, os.error), why:
            errors.append((srcname, dstname, str(why)))
    if errors:
        raise Error, errors


if __name__ == '__main__':
    for c in colors:
        cprint(c,c)
    print "\n"

    def pstatus(state, option=''):
        if state == 0 or state == None:
            cprint('GREEN', 'ok ' + option)
        elif isinstance(state, int):
            cprint('RED', 'failed ' + option)
        else:
            cprint('YELLOW', state)


    pmess("CONF","xyzzy", "/tmp")
    os.system("sleep 1")
    pstatus(0)
    pmess("BUILD","gmsh", "/tmp/foo/bar")
    os.system("sleep 1")
    pstatus(1)

    pmess("BUILD","foo", "/tmp/foo")
    pstatus("skipped")
    pmess("INSTALL","foobar", "/tmp/foobar")
    pstatus(0,"(excellent in fact)")

    debug("unseen debug message")
    verbose = True
    debug("debug message")
    warning("reactor overheating")
    error("reactor exploding")
    fatal("BOOM")
