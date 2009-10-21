from build_packages import *

class NetCDF4(BuildPkg):
    requires = ['netcdf','numpy']
    name = "netCDF4"
    def _install(self):
        return self.sys_log("NETCDF3_DIR=%s python setup-nc3.py install --prefix=%s" % (self.blddir, self.blddir))
    

