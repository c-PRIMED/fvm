# Build Class
# Martin Hunt <mmh@purdue.edu>

"""
Build definitions.
"""
import os, glob, shelve, dirwalk, time, config, sys
from build_utils import fatal, debug, verbose, python_path
from tsort import topological_sort

class Build:

    def __init__(self, cname, topdir, make_path):
        # create build directories
        cwd = os.getcwd()
        topdir = os.path.abspath(topdir)
        self.topdir = topdir
        self.blddir = os.path.join(os.getcwd(), "build-%s" % cname)
        self.logdir = os.path.join(self.blddir, "log")
        self.bindir = os.path.join(self.blddir, "bin")
        self.libdir = os.path.join(self.blddir, "lib")
        self.incdir = os.path.join(self.blddir, "include")
        for p in [self.blddir, self.logdir,
                  self.libdir, self.bindir, self.incdir]:
            if not os.path.isdir(p):
                try:
                    os.makedirs(p)
                except:
                    fatal("error creating directory " + p)

        # load build database
        self.database = shelve.open(os.path.join(self.logdir, 'PACKAGES'))

        # import modules from 'packages' directory
        for path in (os.path.join(topdir, 'config', 'packages'),
                     os.path.join(make_path, 'packages')):
            if os.path.isdir(path):
                debug("importing package modules from %s" % path)
                try:
                    os.chdir(path)
                    sys.path = [path] + sys.path
                except:
                    fatal('Cannot read modules from %s' % path)
                for m in glob.glob('*.py'):
                    m = m.split('.')[0]
                    if m == '__init__':
                        continue
                    debug("importing '%s'" % m)
                    try:
                        loaded = eval(m)
                    except NameError:
                        loaded = False
                    if not loaded:
                        _temp = __import__(m, globals(), locals(), [m])
                        exec('global %s; %s=_temp.%s; self.%s=%s' % (m, m, m, m, m))

        # create package class instances
        self.all_packages = []
        path = os.path.join(topdir, 'config', 'packages')
        os.chdir(path)

        for line in open('SOURCES'):
            cls, sdir = line.split()
            acls = []
            try:
                exec("acls = %s(self, '%s')" % (cls, sdir))
                self.all_packages.append(acls)
                if cls == 'Scons':
                    self.Scons = acls
            except NameError:
                Warning('Package %s in SOURCES, but no build script found.' % cls)

        # load PACKAGES database
        for pname in [p.name for p in self.all_packages]:
            try:
                tmp = self.database[pname]
                debug('database[%s] = %s' % (pname, tmp))
            except KeyError:
                self.database[pname] = 0
            try:
                tmp = self.database[pname + '-status']
                debug('database[%s] = %s' % (pname + '-status', tmp))
            except KeyError:
                self.database[pname + '-status'] = ''

        # cache package instances by name key
        self.pkglist = {}
        for p in self.all_packages:
            self.pkglist[p.name] = p

        # now build dependency list for all packages
        for p in self.all_packages:
            deps = []
            p.deps = []
            if hasattr(p, 'requires'): deps += p.requires
            deps = p.add_required(deps)
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
                if check_timestamp:
                    verbose(1, 'Checking required packages for %s: ' % p.name)
                if p.name == 'python':
                    # FIXME: major hack
                    self.add_build(p, check_timestamp, True)
                    self.add_build([p for p in self.all_packages if p.name == 'matplotlib'][0], check_timestamp, True)
                    self.add_build([p for p in self.all_packages if p.name == 'scipy'][0], check_timestamp, True)
                    self.add_build([p for p in self.all_packages if p.name == 'nose'][0], check_timestamp, True)
                    self.add_build([p for p in self.all_packages if p.name == 'sympy'][0], check_timestamp, True)
                else:
                    self.add_build(p, check_timestamp)
                    
    def add_build(self, pkg, check_timestamp, force=False):
        ''' Add a package to the build list if it or its dependencies have changed. '''
        debug ('add_build %s [%s] %s' % (pkg, pkg.deps, force))
        if self.database[pkg.name + '-status'] == 'installed':
            return 0
        if not force and pkg.installed():
            path = python_path(pkg.name)
            if os.path.commonprefix([path,self.libdir]) != self.libdir:
                return 0
        deps_needed = 0
        for p in pkg.deps:
            deps_needed += self.add_build(p, check_timestamp, force)
        if check_timestamp and deps_needed == 0:
            f, t = dirwalk.DirWalk(pkg.psdir).find_newest()
            debug('%s: %s\t%s' % (pkg.psdir, f, time.strftime("%b %d %Y %H:%M %Z", time.localtime(t))))
            if t < self.database[pkg.name] and self.database[pkg.name + '-status'] == pkg.status():
                return 0
        if not pkg in self.packages:
            if check_timestamp:
                if self.database[pkg.name + '-status'] == pkg.status():
                    verbose(1, '\t%s needs rebuilt' % pkg.name)
                elif self.database[pkg.name + '-status'] == '':
                    verbose(1, '\t%s needs built' % pkg.name)
                else:
                    verbose(1, "\t%s needs rebuilt. Status changed from '%s' to '%s.'"
                            % (pkg.name, self.database[pkg.name + '-status'], pkg.status()))
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

