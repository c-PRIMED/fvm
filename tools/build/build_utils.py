"""
build utility functions.
"""
import sys, os

g_colors = {
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


g_maxlen = 20
g_verbose = 0

def _clear_colors():
    for k in g_colors.keys():
        g_colors[k]=''
        
if (sys.platform=='win32') or ('NOCOLOR' in os.environ) \
        or (os.environ.get('TERM', 'dumb') in ['dumb', 'emacs']) \
        or (not sys.stdout.isatty()):
    _clear_colors()

g_type = {
    'CONF':"%s%s%s"%(g_colors['BLUE'],"Configuring",g_colors['NORMAL']),
    'BUILD':"%s%s%s"%(g_colors['BLUE'],"Building",g_colors['NORMAL']),
    'INSTALL':"%s%s%s"%(g_colors['BLUE'],"Installing",g_colors['NORMAL'])
    }
    
def cprint(col, str):
    try: mycol=g_colors[col]
    except KeyError: mycol=''
    print "%s%s%s" % (mycol, str, g_colors['NORMAL'])

def _niceprint(msg, type=''):
    def print_pat(color):
        print '%s: %s%s%s' % (type, g_colors[color], msg, g_colors['NORMAL'])
    if type == 'ERROR' or type == 'WARNING':
        print_pat('RED')
    elif type=='DEBUG':
        print_pat('CYAN')
    else:
        print_pat('NORMAL')
                      
def debug(msg):
    global g_verbose
    if g_verbose < 2:
        return
    _niceprint(msg, 'DEBUG')

def warning(msg, zone=0):
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
    sr = '%s %s%s%s in %s' % (g_type[type], g_colors['PINK'], pkg, g_colors['NORMAL'], dir)
    global g_maxlen
    g_maxlen = max(g_maxlen, len(sr))
    print "%s :" % sr.ljust(g_maxlen),
    sys.stdout.flush()

if __name__ == '__main__':
    for c in g_colors:
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
    g_verbose = 2
    debug("debug message")
    warning("reactor overheating")
    error("reactor exploding")
    fatal("BOOM")
