import numpy

import fvm.fvmbaseExt as fvmbaseExt
import fvm.models_atyped_double as models


def generate2DMesh(la, lb, ta, tb, tc, nxa, nxb, nxc, nya, nyb):


    # number of cells
    imax = 2*(nxa + nxc) + nxb
    jmax = nya + nyb

    # cell lengths in x and y 
    dx = numpy.zeros(shape=(imax,), dtype='double')
    dy = numpy.zeros(shape=(jmax,), dtype='double')

    dxa = la/nxa
    dxb = lb/nxb
    dxc = tc/nxc

    dx[0:nxa] = dxa
    dx[nxa:nxa+nxc] = dxc
    dx[nxa+nxc:nxa+nxc+nxb] = dxb
    dx[nxa+nxc+nxb:nxa+2*nxc+nxb] = dxc
    dx[nxa+2*nxc+nxb:] = dxa
    
    dya = ta/nya
    dyb = tb/nyb

    dy[0:nya] = dya
    dy[nya:] = dyb
    
    xcoord = numpy.zeros(shape=(imax+1,), dtype='double')
    ycoord = numpy.zeros(shape=(jmax+1,), dtype='double')

    xcoord[0] = ycoord[0] = 0
    
    for i in range(0,imax):
        xcoord[i+1] = xcoord[i] + dx[i]
    for j in range(0,jmax):
        ycoord[j+1] = ycoord[j] + dy[j]

    ## move the x centroid to center of beam
    xmax = xcoord[imax]
    xcoord[:] -= xmax/2.0
    
    # this will store the node index in the mesh for the cartesian node location
    nodeIndex = numpy.zeros(shape=(imax+1,jmax+1), dtype='int')
    nodeIndex[:,:] = -1

    nNodes = 0
    ## mark nodes in left anchor
    for j in range(0, nya+1):
        for i in range(0, nxa+nxc+1):
            if nodeIndex[i,j] == -1:
                nodeIndex[i,j] = nNodes
                nNodes += 1

    ## mark nodes in beam
    for j in range(nya, jmax+1):
        for i in range(nxa, nxa+2*nxc+nxb+1):
            if nodeIndex[i,j] == -1:
                nodeIndex[i,j] = nNodes
                nNodes += 1
                
    ## mark nodes in right anchor
    for j in range(0, nya+1):
        for i in range(nxa+nxc+nxb, imax+1):
            if nodeIndex[i,j] == -1:
                nodeIndex[i,j] = nNodes
                nNodes += 1

    nodeCoordsN = models.newVec3Array(nNodes)
    nodeCoordsA = nodeCoordsN.asNumPyArray()

    nodeCoordsA[:,2] = 0.
    for i in range(0, imax+1):
        for j in range(0, jmax+1):
            ni = nodeIndex[i,j]
            if (ni >= 0):
                nodeCoordsA[ni,0] = xcoord[i]
                nodeCoordsA[ni,1] = ycoord[j]
                

    # this will store the cell index in the mesh for the cartesian cell location
    cellIndex = numpy.zeros(shape=(imax,jmax), dtype='int')
    cellIndex[:,:] = -1

    nCells = 0
    ## mark cells in left anchor
    for j in range(0, nya):
        for i in range(0, nxa+nxc):
            if cellIndex[i,j] == -1:
                cellIndex[i,j] = nCells
                nCells += 1

    ## mark cells in beam
    for j in range(nya, jmax):
        for i in range(nxa, nxa+2*nxc+nxb):
            if cellIndex[i,j] == -1:
                cellIndex[i,j] = nCells
                nCells += 1
                
    ## mark cells in right anchor
    for j in range(0, nya):
        for i in range(nxa+nxc+nxb, imax):
            if cellIndex[i,j] == -1:
                cellIndex[i,j] = nCells
                nCells += 1
                

    nFacesInterior = 0

    ## interior faces of left and right anchors

    nFacesInterior += (nxa+nxc-1)*nya + (nya-1)*(nxa+nxc)
    nFacesInterior += (nxa+nxc-1)*nya + (nya-1)*(nxa+nxc)

    ## interior faces of beam

    nFacesInterior += (nxb+2*nxc-1)*nyb + (nyb-1)*(nxb+2*nxc)

    ## interior faces between anchor and beam

    nFacesInterior += 2*nxc

    nFaceZones = 8

    faceGroupCountN = fvmbaseExt.newIntArray(nFaceZones)
    faceGroupCount = faceGroupCountN.asNumPyArray()

    faceGroupCount[0] = nFacesInterior 

    # bottom of anchors
    faceGroupCount[1] = 2*(nxa+nxc)


    #outer sides of anchors
    faceGroupCount[2] = 2*nya

    #inner sides of anchors
    faceGroupCount[3] = 2*nya

    #beam bottom
    faceGroupCount[4] = nxb

    #beam top
    faceGroupCount[5] = nxb + 2*nxc

    #beam sides
    faceGroupCount[6] = 2*nyb

    # top of anchors
    faceGroupCount[7] = 2*nxa

    nFaces = int(faceGroupCount.sum())

    ## allocate arrays for face nodes 
    faceNodeCountN = fvmbaseExt.newIntArray(nFaces)
    faceNodesN = fvmbaseExt.newIntArray(nFaces*2)
    faceCellsN = fvmbaseExt.newIntArray(nFaces*2)

    faceNodesA = faceNodesN.asNumPyArray()
    faceCellsA = faceCellsN.asNumPyArray()

    faceNodeCountA = faceNodeCountN.asNumPyArray()
    faceNodeCountA[:] = 2
    
    ## reshape for convenience

    faceNodes = faceNodesA.reshape((nFaces,2))
    faceCells = faceCellsA.reshape((nFaces,2))
    
    # interior x faces of left anchor

    nf = 0
    for j in range(0, nya):
        for i in range(0, nxa+nxc-1):
            faceNodes[nf,0] = nodeIndex[i+1,j]
            faceNodes[nf,1] = nodeIndex[i+1,j+1]

            faceCells[nf,0] = cellIndex[i,j]
            faceCells[nf,1] = cellIndex[i+1,j]

            nf += 1

    # interior y faces of left anchor

    for j in range(0, nya-1):
        for i in range(0, nxa+nxc):
            faceNodes[nf,0] = nodeIndex[i+1,j+1]
            faceNodes[nf,1] = nodeIndex[i,j+1]

            faceCells[nf,0] = cellIndex[i,j]
            faceCells[nf,1] = cellIndex[i,j+1]

            nf += 1


    # interior x faces of beam

    for j in range(nya, jmax):
        for i in range(nxa, nxa+2*nxc+nxb-1):
            faceNodes[nf,0] = nodeIndex[i+1,j]
            faceNodes[nf,1] = nodeIndex[i+1,j+1]

            faceCells[nf,0] = cellIndex[i,j]
            faceCells[nf,1] = cellIndex[i+1,j]

            nf += 1

    # interior y faces of beam

    for j in range(nya, jmax-1):
        for i in range(nxa, nxa+2*nxc+nxb):
            faceNodes[nf,0] = nodeIndex[i+1,j+1]
            faceNodes[nf,1] = nodeIndex[i,j+1]

            faceCells[nf,0] = cellIndex[i,j]
            faceCells[nf,1] = cellIndex[i,j+1]

            nf += 1
            
    ## interior x faces in right anchor
    for j in range(0, nya):
        for i in range(nxa+nxc+nxb, imax-1):
            faceNodes[nf,0] = nodeIndex[i+1,j]
            faceNodes[nf,1] = nodeIndex[i+1,j+1]

            faceCells[nf,0] = cellIndex[i,j]
            faceCells[nf,1] = cellIndex[i+1,j]

            nf += 1

    ## interior y faces in right anchor
    for j in range(0, nya-1):
        for i in range(nxa+nxc+nxb, imax):
            faceNodes[nf,0] = nodeIndex[i+1,j+1]
            faceNodes[nf,1] = nodeIndex[i,j+1]

            faceCells[nf,0] = cellIndex[i,j]
            faceCells[nf,1] = cellIndex[i,j+1]

            nf += 1

    ## left interior faces between anchor and beam

    j = nya

    for i in range(nxa, nxa+nxc):
        faceNodes[nf,0] = nodeIndex[i+1,j]
        faceNodes[nf,1] = nodeIndex[i,j]
        
        faceCells[nf,0] = cellIndex[i,j-1]
        faceCells[nf,1] = cellIndex[i,j]
        
        nf += 1
        
    ## right interior faces between anchor and beam

    for i in range(nxa+nxc+nxb, nxa+2*nxc+nxb):
        faceNodes[nf,0] = nodeIndex[i+1,j]
        faceNodes[nf,1] = nodeIndex[i,j]
        
        faceCells[nf,0] = cellIndex[i,j-1]
        faceCells[nf,1] = cellIndex[i,j]
        
        nf += 1
        
    
    
    nb = 0

    # left anchor bottom
    j = 0
    for i in range(0, nxa+nxc):
        faceNodes[nf,0] = nodeIndex[i,j]
        faceNodes[nf,1] = nodeIndex[i+1,j]
        
        faceCells[nf,0] = cellIndex[i,j]
        faceCells[nf,1] = nCells + nb
        
        nf += 1
        nb += 1

    # right anchor bottom

    for i in range(nxa+nxc+nxb,imax):
        faceNodes[nf,0] = nodeIndex[i,j]
        faceNodes[nf,1] = nodeIndex[i+1,j]
        
        faceCells[nf,0] = cellIndex[i,j]
        faceCells[nf,1] = nCells + nb
        
        nf += 1
        nb += 1
        
    # left anchor outer side
    i = 0
    for j in range(0, nya):
        faceNodes[nf,0] = nodeIndex[i,j+1]
        faceNodes[nf,1] = nodeIndex[i,j]
        
        faceCells[nf,0] = cellIndex[i,j]
        faceCells[nf,1] = nCells + nb

        nf += 1
        nb += 1
        
    # right anchor outer side
    i = imax
    for j in range(0, nya):
        faceNodes[nf,0] = nodeIndex[i,j]
        faceNodes[nf,1] = nodeIndex[i,j+1]
        
        faceCells[nf,0] = cellIndex[i-1,j]
        faceCells[nf,1] = nCells + nb

        nf += 1
        nb += 1
        
    # left anchor inner side
    i = nxa+nxc
    for j in range(0, nya):
        faceNodes[nf,0] = nodeIndex[i,j]
        faceNodes[nf,1] = nodeIndex[i,j+1]
        
        faceCells[nf,0] = cellIndex[i-1,j]
        faceCells[nf,1] = nCells + nb

        nf += 1
        nb += 1
        

    # right anchor inner side
    i = nxa+nxc+nxb
    for j in range(0, nya):
        faceNodes[nf,0] = nodeIndex[i,j+1]
        faceNodes[nf,1] = nodeIndex[i,j]
        
        faceCells[nf,0] = cellIndex[i,j]
        faceCells[nf,1] = nCells + nb

        nf += 1
        nb += 1

    # beam bottom
    j = nya
    for i in range(nxa+nxc, nxa+nxc+nxb):
        faceNodes[nf,0] = nodeIndex[i,j]
        faceNodes[nf,1] = nodeIndex[i+1,j]
        
        faceCells[nf,0] = cellIndex[i,j]
        faceCells[nf,1] = nCells + nb
        
        nf += 1
        nb += 1

    # beam top
    j = nya+nyb
    for i in range(nxa, nxa+2*nxc+nxb):
        faceNodes[nf,0] = nodeIndex[i+1,j]
        faceNodes[nf,1] = nodeIndex[i,j]
        
        faceCells[nf,0] = cellIndex[i,j-1]
        faceCells[nf,1] = nCells + nb
        
        nf += 1
        nb += 1
        
    # left beam side
    i = nxa
    for j in range(nya, nya+nyb):
        faceNodes[nf,0] = nodeIndex[i,j+1]
        faceNodes[nf,1] = nodeIndex[i,j]
        
        faceCells[nf,0] = cellIndex[i,j]
        faceCells[nf,1] = nCells + nb

        nf += 1
        nb += 1
        
    # right beam side
    i = nxa+2*nxc+nxb
    for j in range(nya, nya+nyb):
        faceNodes[nf,0] = nodeIndex[i,j]
        faceNodes[nf,1] = nodeIndex[i,j+1]
        
        faceCells[nf,0] = cellIndex[i-1,j]
        faceCells[nf,1] = nCells + nb

        nf += 1
        nb += 1
        
    # left anchor top
    j = nya
    for i in range(0, nxa):
        faceNodes[nf,0] = nodeIndex[i+1,j]
        faceNodes[nf,1] = nodeIndex[i,j]
        
        faceCells[nf,0] = cellIndex[i,j-1]
        faceCells[nf,1] = nCells + nb
        
        nf += 1
        nb += 1
        
    # right anchor top
    for i in range(nxa+2*nxc+nxb, imax):
        faceNodes[nf,0] = nodeIndex[i+1,j]
        faceNodes[nf,1] = nodeIndex[i,j]
        
        faceCells[nf,0] = cellIndex[i,j-1]
        faceCells[nf,1] = nCells + nb
        
        nf += 1
        nb += 1

    #numpy.set_printoptions(threshold=1e8)

    #print faceCellsA[:]
    #print faceNodesA
    mesh = fvmbaseExt.Mesh(2, nCells, nodeCoordsN, faceCellsN,
                           faceNodesN, faceNodeCountN, faceGroupCountN)

    print '2d mesh nnodes = %s' % nNodes
    print 'midpoint node index is %s' % nodeIndex[imax/2,nya]

    return mesh
