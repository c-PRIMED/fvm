from build_packages import *
import config

class Cgal(BuildPkg):

    def find_cgal_inc(self, installed=False, in_build=False):
        f = ''
        if installed or not in_build:
            pathlist=['/usr', '/usr/local', '/opt/local']
            if os.environ.has_key('CGAL_DIR'):
                p = os.path.join(os.environ['CGAL_DIR'], '..', '..')
                pathlist.insert(0, os.path.abspath(p))
        else:
            pathlist=[self.blddir]

        for path in pathlist:
            verbose(2,'Checking for CGAL headers in %s' % path)
            try:
                f = open(os.path.join(path, 'include', 'CGAL', 'version.h'), 'r')
            except:
                pass
            if f:
                verbose(2,'Found CGAL headers.')
                #if path != '/usr':
                #    do_env('HDF5_DIR=%s' % path)
                return f
        return f

    def find_cgal_vers(self, installed=False, in_build=False):
        v = ''
        f = self.find_cgal_inc(installed, in_build)
        if f:
            while True:
                line = f.readline()
                if not line:
                    break
                v = re.findall(r'#define CGAL_VERSION ([^\r\n]*)', line)
                if v:
                    v = v[0]
                    break
            f.close()
        verbose(2, 'CGAL version=%s' % v)
        return v

    def _installed(self):
        # Do we really need to build CGAL or can we use the installed one?
        c = self.find_cgal_inc(in_build=True)
        if c:
            # includes were found in build directory
            c.close()
        else:
            # try system includes
            v = self.find_cgal_vers(installed=True)
            if v:
                v = v.split('.')
                if int(v[0]) > 3 or (int(v[0]) == 3 and int(v[1]) >= 5):
                    return True
        return False

    def _configure(self):
        return self.sys_log("cmake %s -DCMAKE_INSTALL_PREFIX=%s -DWITH_CGAL_Qt3=0 -DWITH_CGAL_Qt4=0 -DCMAKE_BUILD_TYPE=Release" % (self.sdir, self.blddir))
    def _build(self):
        return self.sys_log("make -j%s" % jobs(self.name))
    def _install(self):
        return self.sys_log("make install")


