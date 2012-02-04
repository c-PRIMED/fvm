from build_packages import *
import re, os

class Pcl(BuildPkg):
    requires = ['flann', 'eigen', 'boost', 'qhull']
    
    def _installed(self):
        return find_lib('pcl_common', '1.4')

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
            cmdline += " -DBoost_INCLUDE_DIR=%s" % self.incdir
        if os.path.isfile(os.path.join(self.libdir, 'libqhull6.so')):
            cmdline += " -DQHULL_LIBRARY=%s" % os.path.join(self.libdir, 'libqhull6.so')

        # gcc 4.6.2 does not compile PCL with -O3
        try:
            line = os.popen("/bin/bash -c 'gcc --version 2>&1'").readline()
            version = re.compile(r'[^(]*[^)]*\) ([^\n ]*)').findall(line)[0]
            if version == '4.6.2':
                cmdline += " -DCMAKE_CXX_FLAGS_RELEASE='-O -DNDEBUG'"
        except:
            pass
        return self.sys_log(cmdline)
    def _build(self):
        return self.sys_log("make -j%s" % jobs(self.name))
    def _install(self):
        return self.sys_log("make install")


