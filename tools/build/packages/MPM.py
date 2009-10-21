from build_packages import *

class MPM(BuildPkg):
    name = "MPM"
    requires = ['netcdf', 'netCDF4', 'h5py']

    # just copy the sources we need
    def copy_srcs(self):
        os.makedirs(self.bdir)
        shutil.copy2(os.path.join(self.sdir,'Makefile'), self.bdir)
        for dirname in ['config', 'F95', 'py']:
            src = os.path.join(self.sdir, dirname)
            dst = os.path.join(self.bdir, dirname)
            copytree(src, dst, 1)

    def _configure(self):
        pass
        e = config(self.name,'configname')
        bfile = os.path.join(self.sdir, "config", e)
        if os.path.isfile(bfile):
            os.chdir(os.path.join(self.sdir, "config"))
            self.sys_log("/bin/ln -fs %s CURRENT" % bfile)
            os.chdir(os.path.join(self.bdir, "config"))
            bfile = os.path.join(self.bdir, "config", e)            
            return self.sys_log("/bin/ln -fs %s CURRENT" % bfile)
        else:
            f = open(self.logfile, 'a')
            print >> f, "Cannot open config file %s." % bfile
            f.close()        
        return False
    def _build(self):
        # Don't use parallel make here. Breaks on some systems.
        return self.sys_log("make")
    def _install(self):
        return self.sys_log("make install")
