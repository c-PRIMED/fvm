from build_packages import *
import glob

class Scons(BuildPkg):
    name='scons'
    requires = ['python']
    def set_path(self, unload):
        path = os.path.join(self.blddir, 'build', 'scons')
        sc = glob.glob(path+'/scons-local*/SCons')
        path = '/'.join(sc[0].split('/')[:-1])
        if os.path.isdir(path):
            fix_path('PYTHONPATH', path, 1,unload)
