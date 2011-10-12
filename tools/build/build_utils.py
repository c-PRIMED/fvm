"""
build utility functions.
"""
import sys, os, shutil, re, subprocess
from config import config

# define some exceptions
class FatalException: pass # Does a trace
class CompileException: pass

#colors used for printing messages
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


maxlen = 20
opt_verbose = False
opt_debug = False
myenv = {}

def set_options(options):
    global opt_verbose, opt_debug, opt_jobs
    if options.nocolor:
        clear_colors()
    opt_verbose = options.verbose or options.debug
    opt_debug = options.debug
    if options.jobs:
        opt_jobs = options.jobs
    else:
        try:
            import multiprocessing
            opt_jobs = multiprocessing.cpu_count()
        except ImportError:
            opt_jobs = 1

# number of jobs for make -j argument
def jobs(pkg):
    if config(pkg, 'jobs'):
        return config(pkg, 'jobs')
    return opt_jobs

def _reset_types():
    global g_type
    g_type = {
        'CONF':"%s%s%s" % (colors['BOLD'], "Configuring", colors['NORMAL']),
        'BUILD':"%s%s%s" % (colors['BOLD'], "Building", colors['NORMAL']),
        'INSTALL':"%s%s%s" % (colors['BOLD'], "Installing", colors['NORMAL']),
        'TEST':"%s%s%s" % (colors['BOLD'], "Testing", colors['NORMAL']),
        'CLEAN':"%s%s%s" % (colors['BOLD'], "Cleaning", colors['NORMAL'])
        }

def clear_colors():
    for k in colors.keys():
        colors[k] = ''
    _reset_types()


_reset_types()
if ('NOCOLOR' in os.environ) \
        or (os.environ.get('TERM', 'dumb') in ['dumb', 'emacs']) \
        or (not sys.stdout.isatty()):
    clear_colors()

def cprint(col, str, newline=True):
    try: mycol = colors[col]
    except KeyError: mycol = ''
    if newline:
        print "%s%s%s" % (mycol, str, colors['NORMAL'])
    else:
        print "%s%s%s" % (mycol, str, colors['NORMAL']),

def _niceprint(msg, type=''):
    print_type = True
    def print_pat(color):
        global opt_verbose
        if not opt_verbose:
            print '\n'
        if print_type:
            print '%s: %s%s%s' % (type, colors[color], msg, colors['NORMAL'])
        else:
            print '%s%s%s' % (colors[color], msg, colors['NORMAL'])
    if type == 'ERROR' or type == 'WARNING':
        print_pat('RED')
    elif type == 'DEBUG':
        print_pat('DCYAN')
    elif type == 'VERBOSE':
        print_type = False
        print_pat('DYELLOW')
    else:
        print_pat('NORMAL')

def debug(msg):
    global opt_debug
    if not opt_debug:
        return
    _niceprint(msg, 'DEBUG')

def verbose(level, msg):
    global opt_verbose
    if level <= opt_verbose:
        _niceprint(msg, 'VERBOSE')

def warning(msg):
    _niceprint(msg, 'WARNING')

def verbose_warn(msg):
    global opt_verbose
    if not opt_verbose:
        return
    _niceprint(msg, 'WARNING')

def error(msg):
    _niceprint(msg, 'ERROR')

def fatal(msg, ret=1, trace=1):
    cprint('RED', '%s' % msg)
    raise FatalException

def pmess(type, pkg, dir):
    "print a build message"
    sr = '%s %s%s%s in %s' % (g_type[type], colors['DCYAN'], pkg, colors['NORMAL'], dir)
    global maxlen
    maxlen = max(maxlen, len(sr))
    if opt_verbose:
        print "%s :" % sr.ljust(maxlen),
    else:
        print "%s :" % sr.ljust(maxlen),
        sys.stdout.flush()

def remove_file(name):
    try:
        os.remove(name)
    except OSError:
        pass

def do_env(c, unload=False):
    try:
        a, b = c.split('=')
    except:
        fatal("Cannot parse command: env", c)
    if unload:
        try:
            old = myenv[a].pop()
            os.environ[a] = old
            if not myenv[a]:
                del myenv[a]
            debug("Set %s back to %s" % (a, old))
        except:
            debug("Unset " + a)
            if os.environ.has_key(a):
                os.environ.pop(a)
    else:
        debug("Set %s=%s" % (a, b))
        if not myenv.has_key(a):
            myenv[a] = []
        if os.environ.has_key(a):
            myenv[a].append(os.environ[a])
            if b == '':
                del os.environ[a]
            else:
                os.environ[a] = b
        elif b:
            os.environ[a] = b

def fix_path(varname, val, prepend, unload):
    if unload:
        s1 = "remove"
        s2 = "from"
    else:
        s2 = "to"
        if prepend: s1 = "prepend"
        else: s1 = "append"

    debug("%s '%s' %s %s" % (s1, val, s2, varname))

    try:
        e = os.environ[varname]
    except KeyError:
        e = ''

    if unload:
        p = e.split(':')

        for v in val.split(':'):
            try:
                p.remove(v)
            except ValueError:
                fatal("While trying to remove '%s' from '%s'\n[%s]"
                      % (v, p, e))
        os.environ[varname] = ':'.join(p)
    else:
        if prepend:
            os.environ[varname] = (val + ':' + e).strip(':')
        else:
            os.environ[varname] = (e + ':' + val).strip(':')
    debug("%s=%s" % (varname, os.environ[varname]))


# python handler for Environment Modules
# http://modules.sourceforge.net/
# We can't just load modules with system()
# or popen() because those run in subshells and do not
# modify the current environment.
def module_load(m, unload=False):
    if unload:
        debug("Unoading module " + m)
    else:
        debug("Loading module " + m)

    for line in os.popen("/bin/bash -l -c 'module show %s 2>&1'" % m):
        x = line.split()
        if not x: continue
        if x[0] == 'setenv':
            val = ' '.join(x[2:])
            do_env("%s=%s" % (x[1], val), unload)
        elif x[0] == 'prepend-path':
            val = ' '.join(x[2:])
            fix_path(x[1], val, 1, unload)
        elif x[0] == 'append-path':
            val = ' '.join(x[2:])
            fix_path(x[1], val, 0, unload)
        elif x[0] == 'module-whatis':
            val = ' '.join(x[1:])
        elif x[0] == 'module':
            if x[1] == 'load':
                val = ' '.join(x[2:])
                debug("need to load module " + val)
                module_load(val, unload)
            else:
                warning ("Unknown module command " + x[1])
        else:
            if x[0].find('ERROR') > 0:
                fatal("ERROR loading module '%s' while processing\n'%s'"
                      % (m, line), -1, 0)

# pkg = 'ALL' or package name
# section = 'before' or 'after'
# Runs commands and [un]loads modules and [un]sets env variables
def run_commands(pkg, section):

    # Optionally run a command
    cmd = config(pkg, section)
    if cmd:
        debug("executing: " + cmd)
        s = os.system(cmd)
        if s:
            error("While running 'before' commands. Execution of")
            print cmd
            print "failed."
            sys.exit(-1)

    # Optionally set an environment variable
    for env in config(pkg, 'env'):
        do_env(env, section == 'after')

    # First load modules
    mods = config(pkg, 'modules')
    if mods:
        for m in mods.split():
            module_load(m, section == 'after')


def copytree(src, dst, ctype):
    os.makedirs(dst)
    if ctype <= 0:
        return
    for name in os.listdir(src):
        # filter out things we don't want to copy
        if name in ['.svn', 'CVS'] or name.endswith('.o'):
            continue
        srcname = os.path.join(src, name)
        dstname = os.path.join(dst, name)
        if os.path.islink(srcname):
            linkto = os.readlink(srcname)
            os.symlink(linkto, dstname)
        elif os.path.isdir(srcname):
            copytree(srcname, dstname, ctype)
        elif ctype == 2:
            shutil.copy2(srcname, dstname)
        else:
            os.symlink(srcname, dstname)

def set_python_path(dir, reset=False):
    #py = find_executable('python')
    ver = os.popen("/bin/bash -c 'python -V 2>&1'").readline()
    try:
        a, b = re.compile(r'Python ([^.]*).([^.\n]*)').findall(ver)[0]
    except:
        print "Unable to determine python version."
        print "ver=%s" % ver
        sys.exit(2)

    libpath = os.path.join(dir, 'lib')
    pypath1 = os.path.join(dir, 'lib64', 'python%s.%s' % (a, b), 'site-packages')
    pypath2 = os.path.join(dir, 'lib', 'python%s.%s' % (a, b), 'site-packages')
    pypath = pypath1 + ':' + pypath2

    for p in [pypath1, pypath2]:
        try:
            os.makedirs(p)
        except OSError:
            pass
    if 'PYTHONPATH' in os.environ and not reset:
        orig = os.environ['PYTHONPATH']
    else:
        orig = ''
    os.environ['PYTHONPATH'] = libpath + ':' + pypath
    os.environ['PYTHONPATH'] += ':' + os.path.join(dir, 'bin')
    if orig:
        os.environ['PYTHONPATH'] += ':' + orig

def write_env(bld, cwd, cname):
    idir = os.path.join(bld.blddir, "include")
    # write out env.csh for people who haven't yet learned bash
    env_name = os.path.join(cwd, 'env.csh')
    c_name = os.path.join(cwd, 'env-%s.csh' % cname)
    f = open(env_name, 'w')
    for cmd in config('ALL', 'before'):
        exp = re.findall(r'export (\S+)=(\S+)', cmd)
        if exp:
            f.write('setenv %s %s\n' % (exp[0][0], exp[0][1]))
        else:
            f.write('%s\n' % cmd)

    print >> f, "if ${?LD_LIBRARY_PATH} then"
    print >> f, "  setenv LD_LIBRARY_PATH " + bld.libdir + ":${LD_LIBRARY_PATH};"
    print >> f, "else"
    print >> f, "  setenv LD_LIBRARY_PATH " + bld.libdir + ";"
    print >> f, "endif"

    try:
        if os.environ['PYTHONPATH']:
            print >> f, "setenv PYTHONPATH " + os.environ['PYTHONPATH']
    except:
        pass
    print >> f, "# Add include and lib directories to compiler paths"
    print >> f, "if ${?C_INCLUDE_PATH} then"
    print >> f, "  setenv C_INCLUDE_PATH %s:${C_INCLUDE_PATH};" % idir
    print >> f, "else"
    print >> f, "  setenv C_INCLUDE_PATH %s;" % idir
    print >> f, "endif"

    print >> f, "if ${?CPLUS_INCLUDE_PATH} then"
    print >> f, "  setenv CPLUS_INCLUDE_PATH %s:${CPLUS_INCLUDE_PATH};" % idir
    print >> f, "else"
    print >> f, "  setenv CPLUS_INCLUDE_PATH %s;" % idir
    print >> f, "endif"

    print >> f, "setenv PATH %s:$PATH" % bld.bindir
    print >> f, "\n# Need this to recompile MPM in its directory."
    print >> f, "setenv MEMOSA_CONFNAME %s" % cname
    f.close()
    shutil.copy2(env_name, c_name)

    # write out env.sh
    env_name = os.path.join(cwd, 'env.sh')
    c_name = os.path.join(cwd, 'env-%s.sh' % cname)
    f = open(env_name, 'w')
    for cmd in config('ALL', 'before'):
        f.write('%s\n' % cmd)

    print >> f, "export LD_LIBRARY_PATH=" + bld.libdir + ":$LD_LIBRARY_PATH"
    try:
        if os.environ['PYTHONPATH']:
            print >> f, "export PYTHONPATH=" + os.environ['PYTHONPATH']
    except:
        pass

    print >> f, "# Add include and lib directories to compiler paths"
    print >> f, "export C_INCLUDE_PATH=%s:$C_INCLUDE_PATH" % idir
    print >> f, "export CPLUS_INCLUDE_PATH=%s:$CPLUS_INCLUDE_PATH" % idir
    print >> f, "export PATH=%s:$PATH" % bld.bindir
    print >> f, "\n# Need this to recompile MPM in its directory."
    print >> f, "export MEMOSA_CONFNAME=%s" % cname
    f.close()
    shutil.copy2(env_name, c_name)
    return c_name


def find_executable(executable, path=None):
    """Try to find 'executable' in the directories listed in 'path' (a
    string listing directories separated by 'os.pathsep'; defaults to
    os.environ['PATH']).  Returns the complete filename or None if not
    found
    """
    if path is None:
        path = os.environ['PATH']
    paths = path.split(os.pathsep)
    extlist = ['']
    if sys.platform == 'win32':
        pathext = os.environ['PATHEXT'].lower().split(os.pathsep)
        (base, ext) = os.path.splitext(executable)
        if ext.lower() not in pathext:
            extlist = pathext
    for ext in extlist:
        execname = executable + ext
        if os.path.isfile(execname):
            return execname
        else:
            for p in paths:
                f = os.path.join(p, execname)
                if os.path.isfile(f):
                    return f
    else:
        return None

def python_path(name):
    cmd = "python -c 'import %s; print %s.__path__[0]'" % (name, name)
    try:
        path = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE).stdout.read()
    except:
        path = ''
    return path

def python_package(name, version):
    cmd = "python -c 'import %s; print %s.__version__'" % (name, name)
    try:
        ver = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE).stdout.read()
        # parse version string and convert to numbers
        ver = re.findall(r'([\d]+)', ver)
        if not ver:
            return False

        for i,v in enumerate(ver):
            try:
                ver[i] = int(v)
            except:
                ver[i] = 0

        for v1,v2 in zip(version, ver):
            if v1 > v2:
                return False
            if v2 > v1:
                return True
        return True
    except:
        pass
    return False

if __name__ == '__main__':
    for c in colors:
        cprint(c, c)
    print "\n"

    def pstatus(state, option=''):
        if state == 0 or state == None:
            cprint('GREEN', 'ok ' + option)
        elif isinstance(state, int):
            cprint('RED', 'failed ' + option)
        else:
            cprint('YELLOW', state)


    pmess("CONF", "xyzzy", "/tmp")
    os.system("sleep 1")
    pstatus(0)
    pmess("BUILD", "gmsh", "/tmp/foo/bar")
    os.system("sleep 1")
    pstatus(1)

    pmess("BUILD", "foo", "/tmp/foo")
    pmess("INSTALL", "foobar", "/tmp/foobar")
    pstatus(0, "(excellent in fact)")

    debug("unseen debug message")
    verbose("unseen verbose message")
    opt_verbose = True
    verbose("verbose message")
    debug("unseen debug message")
    opt_debug = True
    debug("debug message")
    warning("reactor overheating")
    error("reactor exploding")
    fatal("BOOM")
