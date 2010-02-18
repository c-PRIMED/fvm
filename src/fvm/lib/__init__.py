# import this for all fvm scripts
import sys
sys.setdlopenflags(0x100|0x2)

class AtypeException: pass

global atype, models, exporters

def set_atype(t):
    global atype, models, exporters
    atype = t
    if atype == 'double':
        import fvm.models_atyped_double as models
        import fvm.exporters_atyped_double as exporters
    elif atype == 'tangent':
        import fvm.models_atyped_tangent_double as models
        import fvm.exporters_atyped_tangent_double as exporters
    else:
        raise AtypeException

def dump_hdf5(name, v, desc=''):
    print 'HDF5:%s:5FDH' % repr({'name': name, 'desc': desc, 'value':v})
