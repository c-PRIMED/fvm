
from enthought.traits.api import      List, Dict, Tuple, Trait, CInt, Float


from HasCImpl import CMethod
from MultiField import MultiField, ArrayIndexTrait
from StorageSite import StorageSite
from Matrix import Matrix
from Array import Array
from Connectivity import Connectivity

from CDict import CDictTrait
ArrayIndexTType = 'pair<const Field*,const StorageSite*>'
MultiFieldMatrixIndexTType = 'pair<' + ArrayIndexTType + ',' + ArrayIndexTType + '>'

class MultiFieldMatrix(Matrix):

    templateType = ""
    matrices = CDictTrait(MultiFieldMatrixIndexTType+',Matrix*')

    coarseSizes = CDictTrait(ArrayIndexTType + ',int')

    coarseGhostSizes = CDictTrait(ArrayIndexTType + ',int')
    
    coarseMappers = CDictTrait(MultiFieldMatrixIndexTType + ',OneToOneIndexMap*')

    coarseSites = CDictTrait(ArrayIndexTType + ',const StorageSite*')

    coarseToFineMappings = CDictTrait(ArrayIndexTType +
                                      ',const CRConnectivity*')

    coarseConnectivities = CDictTrait(MultiFieldMatrixIndexTType+
                                      ',const Connectivity*')

    coarseMatrices = CDictTrait(MultiFieldMatrixIndexTType+',Matrix*')

    CMethod('_createCoarsening',None,MultiField,CInt,Float)
    CMethod('_syncGhostCoarsening',None,MultiField)
    CMethod('_createCoarseToFineMapping',None,MultiField)
    CMethod('_createCoarseConnectivity',None,MultiField)
    CMethod('_createCoarseMatrices',None,MultiField)
    CMethod('injectResidual',None,MultiField,MultiField,MultiField)
    CMethod('correctSolution',None,MultiField,MultiField,MultiField)
    
            
    def addMatrix(self,rowField,rowSite,colField,colSite,m):
        self.matrices[((rowField,rowSite),(colField,colSite))] = m

    def __setitem__(self, p, m):
        self.matrices[p] = m

    def __getitem__(self, p):
        return self.matrices[p]

    def __contains__(self, p):
        return p in self.matrices
    
    def initAssembly(self):
        for m in self.matrices.values():
            m.initAssembly()

    def getSize(self):
        size=0
        for k,m in self.matrices.iteritems():
            if (k[0][0],k[0][1]) == (k[1][0],k[1][1]):
                size += k[0][1].selfCount
        return size
        
MultiFieldMatrixTrait  = Trait(MultiFieldMatrix)
