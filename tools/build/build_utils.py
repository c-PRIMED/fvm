"""
build utility functions.
"""
import sys, os
from config import config

colors = {
    'BOLD'  :'\033[1m',
    'RED'   :'\033[1;31m',
    'GREEN' :'\033[1;32m',
    'YELLOW':'\033[1;33m',
    'PINK'  :'\033[1;35m',
    'BLUE'  :'\033[1;34m',
    'CYAN'  :'\033[1;36m',
    'NORMAL':'\033[0m'
    }
"colors used for printing messages"


maxlen = 20
verbose = False

def _clear_colors():
    for k in colors.keys():
        colors[k]=''
        
if (sys.platform=='win32') or ('NOCOLOR' in os.environ) \
        or (os.environ.get('TERM', 'dumb') in ['dumb', 'emacs']) \
        or (not sys.stdout.isatty()):
    _clear_colors()

g_type = {
    'CONF':"%s%s%s"%(colors['BLUE'],"Configuring",colors['NORMAL']),
    'BUILD':"%s%s%s"%(colors['BLUE'],"Building",colors['NORMAL']),
    'INSTALL':"%s%s%s"%(colors['BLUE'],"Installing",colors['NORMAL'])
    }
    
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
        print_pat('CYAN')
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

def fatal(msg, ret=1):
    import traceback
    cprint('RED', '%s' % msg)
    traceback.print_stack()
    sys.exit(ret)

def pmess(type, pkg, dir):
    "print a build message"
    sr = '%s %s%s%s in %s' % (g_type[type], colors['PINK'], pkg, colors['NORMAL'], dir)
    global maxlen
    maxlen = max(maxlen, len(sr))
    print "%s :" % sr.ljust(maxlen),
    sys.stdout.flush()

def run_commands(section, pkg):
    if pkg == 0:
        c = config(section, 0)
    else:
        c = config(pkg, section)
    if c:
        for cmd in c:
            s = os.system(cmd)
            if s:
                error("While running 'before' commands. Execution of")
                print cmd
                print "failed."
                sys.exit(-1)


if __name__ == '__main__':
    for c in colors:
        cprint(c,c)
    print "\n"


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
