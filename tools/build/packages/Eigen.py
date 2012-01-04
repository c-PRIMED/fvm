from build_packages import *
import shutil

class Eigen(BuildPkg):
    def _install(self):
        for idir in ['Eigen', 'unsupported']:
            shutil.copytree(idir, os.path.join(self.blddir, "include", idir))
        return 0
