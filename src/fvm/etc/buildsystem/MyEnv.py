import sys, os, os.path
from SCons.Script.SConscript import SConsEnvironment
from SCons.Util import CLVar
import SCons.Tool
import SCons.Errors
import Arch
import MyTools
import Packages
import distutils.sysconfig
import types
import re

optionsList = [
  ('VERSION', 'debug', 'Build version (default debug, profile)'),
  ('OPENMP', False, 'Compile with OpenMP'),
  ('COMPILER', None, 'C++ compiler (default is the one found in path)'),
  ('PARALLEL', False, 'Build using MPI'),
  ('ATYPES', ['double'],'The atypes to build (e.g., float,double,tangent)'),
  ('BUILDDIR', None, 'build directory (default build)'),
  ('PACKAGESDIR', None, 'packages tree (default packages)'),
  ('COMPACTOUTPUT', True, 'Use compacted form of build command lines in ouotput'),
 ]

swigre = re.compile('(.*).i')

def setBaseVars(env,cla):

    for o in optionsList:
        env[o[0]] = o[1]

    
    for k,v in cla.iteritems():
        if v in ['True','true']:
            v = True
        elif v in ['False', 'false']:
            v = False
            
        env[k] = v

    env['DEBUG'] = env['VERSION'] == 'debug'
        
    env['ARCH'] = Arch.getArch()
    sconstructDir = os.path.abspath(str(env.fs.SConstruct_dir))
    #print sconstructDir
    env['TOPDIR'] = sconstructDir
    env['SRCDIR'] = '$TOPDIR/src'

    if env['BUILDDIR'] is None:
        env['BUILDDIR'] = '$TOPDIR/build'

    if isinstance(env['ATYPES'],types.StringType):
        env['ATYPES'] = [atype.lower() for atype in env['ATYPES'].split(',')]

    # Third-party packages directory variables
    if env['PACKAGESDIR'] is None:
        env['PACKAGESDIR'] = '$TOPDIR/packages/$ARCH'
    env['PACKAGESLIBDIR'] = '$PACKAGESDIR/lib'

    # Propagate the following env variables
    for var in ['TMP', 'TEMP', 'LD_LIBRARY_PATH', 'PATH']:
        if var in os.environ:
            env['ENV'][var] = os.environ[var]

    # Redefine command string variables for compact output.
    if env['COMPACTOUTPUT']:
        for var in ['SHCXXCOMSTR', 'CXXCOMSTR', 'CCCOMSTR', 'SHCCCOMSTR']:
            env[var] = 'Compiling $SOURCE'

        for var in ['SHLINKCOMSTR', 'LINKCOMSTR', 'LDMODULECOMSTR']:
            env[var] = 'Linking $TARGET'

        env['SWIGCOMSTR'] = 'Building swig module for $SOURCE'

def setVersionPaths(env):
    env['BUILDVERSION'] = '$COMPILER/$VERSION'
    env['BUILDVERSIONDIR'] = '$BUILDDIR/$ARCH/$BUILDVERSION'
    env['AUTOGENDIR'] = '$BUILDVERSIONDIR/autogen'
    buildvdir =  os.path.join(env['TOPDIR'],
                              env.GetBuildPath('$BUILDVERSIONDIR'))
    if not os.access(buildvdir,os.F_OK):
        os.makedirs(buildvdir)
        
class MyEnvironment(SConsEnvironment):

    ### mapping between components and their directories. a component
    ### is a directory anywhere in the src tree that has a
    ### directory.scons file. Each such .scons can specify multiple targets.
    _compDirMap = {}

    def __init__(self, platform=None, **kwargs):
        SConsEnvironment.__init__(self, platform=platform, tools=[])

        self._targets = []

        ### our build process requires that we run everything from the top level dir
        self.SConscriptChdir(False)

        ### extract our specific command line variables 
        commandLineOptions = self.getCommandLineOptions()

        ### set the arch and compiler independent options, overriding
        ### defauls with values specified on the command line
        
        setBaseVars(self,commandLineOptions)

        ### now load the architecture specific defaults
        Arch.updateEnv(self)

        ### should allow for a mechanism for user to override
        ### architecture specific defaults. currently we can only
        ### handle specification of the compiler
        

        ### since our build paths depend on the arch/compiler settings
        ### we can only set these paths now
        
        setVersionPaths(self)

        ### keep the sconsign database local to each build rather
        ### than polluting the source tree
        self.SConsignFile('$BUILDVERSIONDIR/sconsign')

        ### for emacs compile mode
        SCons.Util.display("scons: Entering directory `%s'" %  os.getcwd())


    def constructComponents(self):
        for component, dir in self._compDirMap.iteritems():
            env = self.myClone(newTargets=True)

            env['COMPONENT'] = component
            env['COMPONENTSRCDIR'] = dir
            env['COMPONENTBUILDDIR'] = '$BUILDVERSIONDIR/obj/%s' % component
            env['COMPONENTTARGETDIR'] = '$BUILDVERSIONDIR/bin'

            env.AppendUnique(CPPPATH=['$COMPONENTSRCDIR'])
            env.AppendUnique(CPPPATH=['$COMPONENTBUILDDIR'])
            env.AppendUnique(LIBPATH=['$COMPONENTTARGETDIR'])

            script = '$COMPONENTSRCDIR/${COMPONENT}.scons'
            buildDir = '$COMPONENTBUILDDIR'
            srcDir = '$COMPONENTSRCDIR'
            exports = {'env': env}
            
            env.SConscript(script, build_dir=buildDir, src_dir=srcDir,
                           duplicate=False, exports=exports)

            self.Alias(component, env._targets)
            self.Default(component)


    def Copy(self, tools=None, toolpath=[], **kwargs):
        return self


    def createExe(self, target, source=[], deplibs=[], dynlibs=[], **kwargs):
        env = self.myClone()
        env._updateEnv(target, deplibs)

        nodes = env.Program('$COMPONENTTARGETDIR/'+target, source, **kwargs)

        self._addTarget(nodes)
        #self.Alias(target, nodes)
        for lib in dynlibs:
            atypedLib = '_atyped' in lib
            if atypedLib:
                for atype in self['ATYPES']:
                    self.Alias(target,'%s_%s' % (lib, atype))
            else:
                self.Alias(target, lib)
        return nodes


    def createSharedLibrary(self, target, sources=[],
                            deplibs=[],
                            createPublicHeader=True,
                            targetAType=None,
                            **kwargs):

        env = self.myClone()
        env._updateEnv(target, deplibs,targetAType)
        env.AppendUnique(CPPDEFINES='RLOG_COMPONENT=%s' % target)

        env.Append(CPPPATH=['$AUTOGENDIR'])
        env['CURRENTLIB'] = target.upper()
        if Arch.isWindows():
            env.AppendUnique(CPPDEFINES= CLVar('BUILDING_' + env['CURRENTLIB']))

        nodes = env.SharedLibrary('$COMPONENTTARGETDIR/'+target, sources, **kwargs)

        self._addTarget(nodes)
#        if target not in self._compDirMap.keys():
#            self.Alias(target, nodes)
        return nodes

    def createATypedSharedLibrary(self, target, sources=[], deplibs=[],
                                  skip=[],**kwargs):
        env = self.myClone()

        nodes = []
        for atype in self['ATYPES']:
            if atype not in skip:
                libEnv = env.myClone()
                if Arch.isWindows():
                    libEnv.AppendUnique(CPPDEFINES= CLVar('BUILDING_' + target.upper()))
                libTarget = '%s_%s' % (target, atype)
                libEnv['SHOBJSUFFIX'] = '_%s%s' % (atype, self['SHOBJSUFFIX'])
                nodes += libEnv.createSharedLibrary(libTarget, sources, deplibs,
                                                    createPublicHeader=False,
                                                    targetAType=atype,
                                                    **kwargs)

        return nodes

 
    def createSwigModule(self, target, sources=[], deplibs=[], atype=None,**kwargs):

        env = self.myClone()

        env.loadTools(['swig'])
        
        env.Append(CPPPATH=['$AUTOGENDIR'])
        env['SHLIBPREFIX'] = ""
        if env['PARALLEL']:
            env.AppendUnique(SWIGFLAGS = CLVar('-DFVM_PARALLEL -c++ -python -module %s' % target))
        else:
            env.AppendUnique(SWIGFLAGS = CLVar('-c++ -python -module %s' % target))
        env['SWIGOUTDIR'] = '$COMPONENTTARGETDIR'


        env.AppendUnique(CPPDEFINES='RLOG_COMPONENT=%s' % target[1:-3])
        deplibs += ['python' ]
        env._updateEnv(target, deplibs,atype,True)

#        wrapSources = [ swigre.sub(r'\1_wrap.cc', source) for source in sources]

        nodes = env.SharedLibrary('$COMPONENTTARGETDIR/_'+target, sources, targetAType=atype,**kwargs)

        #self._addTarget(wrapSources)
        self._addTarget(nodes)
        return nodes

    def createATypedSwigModule(self, target, sources=[], deplibs=[],
                                  skip=[],**kwargs):
        env = self.myClone()
        env.loadTools(['swig'])

        nodes = []
        for atype in self['ATYPES']:
            if atype not in skip:
                libEnv = env.myClone()
                libTarget = '%s_%s' % (target, atype)
                libEnv.Append(CPPPATH=['$SRCDIR/modules/atypes/%s' % atype])
                libEnv.AppendUnique(SWIGPATH=['$SRCDIR/modules/atypes/%s' % atype])
                libEnv.AppendUnique(SWIGFLAGS=CLVar('-module %s_atyped_%s' % (target,atype)))
                
                libEnv['SHOBJSUFFIX'] = '_%s%s' % (atype, self['SHOBJSUFFIX'])
                nodes += libEnv.createSwigModule(libTarget, sources, deplibs,
                                                 createPublicHeader=False,
                                                 atype=atype,
                                                 **kwargs)

        return nodes

 
    def getAbsPath(self, pathname):
        pathname = self.GetBuildPath(pathname)
        if not os.path.isabs(pathname):
            pathname = os.path.join(self['TOPDIR'], pathname)
        return pathname

    def locateComponents(self):
        env = self.myClone()
        def findComponents(dir):
            for entry in os.listdir(dir):
                if entry in ['.svn']:
                    continue
                subdir = os.path.join(dir, entry)
                if os.path.isdir(subdir):
                    # skip the parallel module for non-parallel builds
                    if entry != 'parallel' or env['PARALLEL']:
                        sconsFile = os.path.join(subdir,'%s.scons' % entry)
                        if os.path.isfile(sconsFile):
                            yield entry, subdir
                        else:
                            for comp in findComponents(subdir):
                                yield comp
        srcDir = self.getAbsPath('$SRCDIR')
        for comp in findComponents(srcDir):
            self._compDirMap[comp[0]] = comp[1]

        #print self._compDirMap


    def loadTools(self, tools):
        for tool in tools:
            SConsEnvironment.Tool(self, tool, MyTools.__path__)

    def _addTarget(self, nodes):
        self._targets.extend(nodes)

        
    def myClone(self, newTargets=False):
        newEnv = SConsEnvironment.Clone(self)
        if newTargets:
            newEnv._targets = []
        else:
            newEnv._targets = self._targets
        return newEnv

    def getCommandLineOptions(self):
        from SCons.Script.SConscript import GlobalDict
        if GlobalDict is None:
            SCons.Script.SConscript.BuildDefaultGlobals()

        validOptions = [opt[0] for opt in optionsList]
        args = {}
        for arg, value in GlobalDict['ARGUMENTS'].items():
            arg = arg.upper()
            if arg not in validOptions:
                raise SCons.Errors.UserError, 'Invalid option: %s' % arg
            args[arg] = value

        ###options = SCons.Options.Options(args=args)
        ###options.AddOptions(*CLOptions.optionsList)
        return args

    def _updateEnv(self, target, deplibs, targetAType=None, updateSwigPath=False):

        if targetAType is not None:
            atypeDir = '$SRCDIR/modules/atypes/%s' % targetAType
            self.AppendUnique(CPPPATH=[atypeDir])
            if updateSwigPath:
                self.AppendUnique(SWIGPATH=[atypeDir])

        for dep in deplibs:
            if targetAType is not None and '_atyped' in dep:
                dep += '_' + targetAType

            if hasattr(Packages, dep):
                getattr(Packages, dep).generate(self)
            else:
                if dep in self._compDirMap.keys():
                    depDir = self._compDirMap[dep]
                    self.AppendUnique(CPPPATH=[depDir])
                    if updateSwigPath:
                        self.AppendUnique(SWIGPATH=[depDir])
                        if targetAType is not None:
                            self['SWIGCFILESUFFIX']   = '_wrap_%s$CFILESUFFIX' % targetAType
                            self['SWIGCXXFILESUFFIX'] = '_wrap_%s$CXXFILESUFFIX' % targetAType
                self.AppendUnique(LIBS=[dep])

        # Update env with the system libs.
        self.Append(LIBS=self['SYSLIBS'])


SCons.Environment.Environment = MyEnvironment
