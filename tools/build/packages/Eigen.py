from build_packages import *
import shutil

class Eigen(BuildPkg):
    def _install(self):
        for idir in ['Eigen', 'unsupported']:
            dname = os.path.join(self.blddir, "include", idir)
            shutil.rmtree(dname, ignore_errors=1)
            shutil.copytree(idir, dname)
        return 0
