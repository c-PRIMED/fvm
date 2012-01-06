from build_packages import *
import config

class Pcl(BuildPkg):
    requires = ['flann', 'eigen', 'boost', 'qhull']
    
    def _configure(self):
        cmdline = "cmake %s -DBUILD_visualization=OFF -DCMAKE_INSTALL_PREFIX=%s -DCMAKE_BUILD_TYPE=Release" % (self.sdir, self.blddir)
        if os.path.isfile(os.path.join(self.libdir, 'libboost_date_time.so')):
            cmdline += " -DBoost_DATE_TIME_LIBRARY=%s" % os.path.join(self.libdir, 'libboost_date_time.so')
            cmdline += " -DBoost_DATE_TIME_LIBRARY_RELEASE=%s" % os.path.join(self.libdir, 'libboost_date_time.so')
            cmdline += " -DBoost_FILESYSTEM_LIBRARY=%s" % os.path.join(self.libdir, 'libboost_filesystem.so')
            cmdline += " -DBoost_FILESYSTEM_LIBRARY_RELEASE=%s" % os.path.join(self.libdir, 'libboost_filesystem.so')
            cmdline += " -DBoost_IOSTREAMS_LIBRARY=%s" % os.path.join(self.libdir, 'libboost_iostreams.so')
            cmdline += " -DBoost_IOSTREAMS_LIBRARY_RELEASE=%s" % os.path.join(self.libdir, 'libboost_iostreams.so')
            cmdline += " -DBoost_SYSTEM_LIBRARY=%s" % os.path.join(self.libdir, 'libboost_system.so')
            cmdline += " -DBoost_SYSTEM_LIBRARY_RELEASE=%s" % os.path.join(self.libdir, 'libboost_system.so')
            cmdline += " -DBoost_THREAD_LIBRARY=%s" % os.path.join(self.libdir, 'libboost_thread.so')            
            cmdline += " -DBoost_THREAD_LIBRARY_RELEASE=%s" % os.path.join(self.libdir, 'libboost_thread.so')            
        return self.sys_log(cmdline)
    def _build(self):
        return self.sys_log("make -j%s" % jobs(self.name))
    def _install(self):
        return self.sys_log("make install")


