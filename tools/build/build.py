# Build Class
# Martin Hunt <mmh@purdue.edu>

"""
Build definitions.
"""
import os, glob, shelve, dirwalk, time, config
from build_utils import fatal, debug, verbose
from tsort import topological_sort

class Build:

    def __init__(self, cname, srcpath):

        # create build directories
        self.topdir = srcpath
        self.blddir = os.path.join(os.getcwd(), "build-%s" % cname)
        self.logdir = os.path.join(self.blddir, "log")
        self.bindir = os.path.join(self.blddir, "bin")
        self.libdir = os.path.join(self.blddir, "lib")
        for p in [self.blddir, self.logdir,
                  self.libdir, self.bindir ]:
            if not os.path.isdir(p):
                try:
                    os.makedirs(p)
                except:
                    fatal("error creating directory " + p)

        # load build database
        self.database = shelve.open(os.path.join(self.logdir, 'PACKAGES'))
        cwd = os.getcwd()

        # import modules from 'packages' directory
        os.chdir(os.path.join(srcpath, 'tools', 'build', 'packages'))
        for m in glob.glob('*.py'):
            m = m.split('.')[0]
            if m == '__init__':
                continue
            debug("importing %s" % m)
            temp = __import__('packages.' + m, globals(), locals(), [m])
            exec('global %s; %s=temp.%s; self.%s=%s' % (m, m, m, m, m))

        # create package classes
        self.all_packages = []
        for line in open('SOURCES'):
            cls, sdir = line.split()
            acls = []
            exec("acls = %s(self, '%s')" % (cls, sdir))
            self.all_packages.append(acls)
            if cls == 'Scons':
                self.Scons = acls

        # load PACKAGES database
        for pname in [p.name for p in self.all_packages]:
            try:
                tmp = self.database[pname]
                debug('database[%s] = %s' % (pname, tmp))
            except KeyError:
                self.database[pname] = 0

        # cache package instances by name key
        self.pkglist = {}
        for p in self.all_packages:
            self.pkglist[p.name] = p

        # now built dependency list for all packages
        for p in self.all_packages:
            deps = []
            p.deps = []
            if hasattr(p, 'requires'): deps += p.requires
            for dep in deps:
                try:
                    p.deps.append(self.pkglist[dep])
                except:
                    fatal('Dependency %s for %s not recognized' % (dep, p.name))

        self.packages = []
        self.build_pkg_list(True)
        if self.packages:
            self.reorder_pkgs()
        os.chdir(cwd)

    def build_pkg_list(self, check_timestamp):
        ''' Build list of packages to build from config file. '''
        for p in self.all_packages:
            x = config.config(p.name, 'Build')
            if x != '' and eval(x):
                self.add_build(p, check_timestamp)

    def add_build(self, pkg, check_timestamp):
        ''' Add a package to the build list if it or its dependencies have changed. '''
        debug ('add_build %s [%s]' % (pkg, pkg.deps))
        deps_needed = 0
        for p in pkg.deps:
            deps_needed += self.add_build(p, check_timestamp)
        if check_timestamp and deps_needed == 0:
            f, t = dirwalk.DirWalk(pkg.psdir).find_newest()
            debug('%s: %s\t%s' % (pkg.psdir, f, time.strftime("%b %d %Y %H:%M %Z", time.localtime(t))))
            if t < self.database[pkg.name]:
                return 0
        if not pkg in self.packages:
            if check_timestamp:
                verbose(1, '%s needs rebuilt' % pkg.name)
            self.packages.append(pkg)
        return 1

    def reorder_pkgs(self):
        ''' Sort packages so they are built in the correct order. '''
        deplist = []
        for pkg in self.packages:
            for dep in pkg.deps:
                if dep in self.packages:
                    deplist.append((dep, pkg))

        #print "PACKAGES=%s" % [p.name for p in self.packages]
        #print "DEPLIST=%s" % deplist
        self.packages = topological_sort(self.packages, deplist)


    def done(self):
        self.database.close()

