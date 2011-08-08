from build_packages import *

# The FreeType Project
# http://freetype.sourceforge.net

class Freetype(BuildPkg):
    def _installed(self):
        for path in [os.path.join(self.blddir, 'include'), '/usr/include', '/usr/local/include']:
            verbose(2,'Checking for freetype headers in %s' % path)
            f = ''
            try:
                f = open(os.path.join(path, 'ft2build.h'), 'r')
                f.close()
                verbose(2,'Found freetype headers.')
                return True
            except:
                pass
        return False

    def _configure(self):
        return self.sys_log("%s/configure --program-prefix=g -prefix=%s" % (self.sdir, self.blddir))
    def _build(self):
        return self.sys_log("make -j%s" % jobs(self.name))
    def _install(self):
        return self.sys_log("make install")
