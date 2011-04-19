from build_packages import *
import shutil

class Ufconfig(BuildPkg):
    def _install(self):
        # Need to rename the unpack directory because UMFPACK
        # is hardcoded to look for it with a certain name
        newdir = os.path.join(os.path.split(self.sdir)[0], 'UFconfig')
        shutil.move(self.sdir, newdir)
        self.sdir = newdir
        self.bdir = newdir