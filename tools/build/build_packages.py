# Build packages
# Martin Hunt <mmh@purdue.edu>

"""
Build package definitions.
"""

import testing, string, time
from build_utils import *
from config import *
from build import Build

# superclass for all packages
class BuildPkg(Build):
    # Some packages require building in the source directory.
    # Others require separate build directories.
    # We don't want to build in our subversion sources. So we define
    # 'copy_sources' as follows:
    #     0 - no copy. Build directory is separate from sources.
    #     1 - copy with symlinks
    #     2 - full copy
    # For tarballs, copy_sources of nonzero means source and build directories are not the same.    
    copy_sources = 0

    def __init__(self, bld, sdir):
        debug('BuildPkg Init (%s,%s,%s)' % (self, bld, sdir))
        if not hasattr(self, 'name'):
            self.name = string.lower(self.__class__.__name__)
        self.__class__.name = self.name
        self.is_pkg = (sdir.split('/')[0] == 'pkgs')
        self.psdir = self.sdir = os.path.join(bld.topdir, sdir)
        self.bdir = os.path.join(bld.blddir, "build", self.name)
        self.blddir = bld.blddir
        self.logdir = bld.logdir
        self.libdir = bld.libdir
        self.bindir = bld.bindir
        self.bld = bld

    def pstatus(self, state, option=''):
        "update build status"
        if state == 0 or state == None:
            cprint('GREEN', 'ok ' + option)
        elif isinstance(state, int):
            cprint('RED', 'failed ' + option)
            print "Check contents of %s\n" % self.logfile
            os.system("tail -40 " + self.logfile)
            sys.exit(state)
        else:
            cprint('YELLOW', state)

    def copy_srcs(self):
        copytree(self.sdir, self.bdir, self.copy_sources)

    # unpack tarball
    def unpack_srcs(self):
        dst = self.bdir
        src = self.sdir
        suffix = src.split('.')[-1]

        compress = None
        if suffix == 'tgz' or suffix == 'gz':
            compress = 'z'
        elif suffix == 'bz2':
            compress = 'j'
        elif suffix == 'tar':
            compress = ''

        if compress == None:
            # it wasn't a tarball
            self.copy_srcs()
            return

        if self.copy_sources:
            self.sdir = os.path.join(dst, 'SOURCES')
            os.makedirs(self.sdir)
            self.bdir = os.path.join(dst, 'BUILD')
            os.makedirs(self.bdir)
        else:
            self.sdir = self.bdir
            os.makedirs(self.bdir)

        # If tarballs have everything in a directory with the same name as the tarball,
        # strip out that directory name.
        name = os.path.basename(src.rstrip('.' + suffix).rstrip('.tar'))
        dir = os.popen("/bin/bash -c 'tar -%stf %s 2> /dev/null'" % (compress, src)).readline().split()[0]
        if dir.startswith('%s/' % name):
            strip = '--strip-components 1'
        else:
            strip = ''
        os.system('tar -C %s %s -%sxf %s' % (self.sdir, strip, compress, src))

    def configure(self):
        self.state = 'configure'
        self.bld.database[self.name] = 0 # mark this package as not built        
        self.logfile = os.path.join(self.logdir, self.name + "-conf.log")
        remove_file(self.logfile)
        pmess("CONF", self.name, self.bdir)
        # remove any old sources
        os.system("/bin/rm -rf %s" % self.bdir)
        # get new sources
        self.unpack_srcs()

        os.chdir(self.bdir)
        run_commands(self.name, 'before')
        self.build_start_time = time.time()
        self.pstatus(self._configure())

    def clean(self):
        self.state = 'clean'
        self.logfile = ''
        self._clean()

    def build(self):
        self.state = 'build'
        self.logfile = os.path.join(self.logdir, self.name + "-build.log")
        remove_file(self.logfile)
        pmess("BUILD", self.name, self.bdir)
        os.chdir(self.bdir)
        self.pstatus(self._build())

    def install(self):
        self.state = 'install'
        self.logfile = os.path.join(self.logdir, self.name + "-install.log")
        remove_file(self.logfile)
        pmess("INSTALL", self.name, self.blddir)
        os.chdir(self.bdir)
        self.pstatus(self._install())
        debug('database[%s] = %s' % (self.name, self.build_start_time))
        self.bld.database[self.name] = self.build_start_time
        run_commands(self.name, 'after')

    def test(self):
        self.state = 'testing'
        self.logfile = os.path.join(self.logdir, self.name + "-test.xml")
        remove_file(self.logfile)
        pmess("TEST", self.name, self.blddir)
        ok, errs = self._test()
        if errs:
            cprint('YELLOW', "%s OK, %s FAIL" % (ok, errs))
        else:
            cprint('GREEN', "%s OK" % ok)
        return ok, errs

    def sys_log(self, cmd):
        "Execute a system call and log the result."
        # get configuration variable
        e = config(self.name, self.state)
        e = e.replace('BUILDDIR', self.blddir)
        e = e.replace('SRCDIR', self.sdir)
        e = e.replace('TMPBDIR', self.bdir)
        e = e.replace('LOGDIR', self.logdir)
        cmd = cmd + " " + e
        debug(cmd)
        if self.logfile != '':
            f = open(self.logfile, 'a')
            print >> f, "EXECUTING:", cmd
            f.close()
            cmd = cmd + " >>" + self.logfile + " 2>&1"
        return os.system("/bin/bash -c '%s'" % cmd)

    # subclasses must redefine these
    def _configure(self):
        pass
    def _clean(self):
        pass
    def _build(self):
        pass
    def _install(self):
        pass
    def _test(self):
        return testing.do_tests(self.name, self.sdir, self.logfile)

