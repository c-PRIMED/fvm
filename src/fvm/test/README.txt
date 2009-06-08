DESCRIPTION OF TESTS FOR FVM

Tests in TESTS

Fvm001 – testLinearSolver  ???
Fvm002 – FvmTestFlowModel.py  ???


Tests in Octree/TESTS

General purpose: given a list of points with Cartesian coordinates (x, y, z),  build up an Octree structure of the points.
For any given point (x0, y0, z0), search the Octree to find out its nearest neighbor or neighbors within a range.
In the test, an Octree structure is built on the cell centroids of a 32x32 Cartesian mesh.  Then a number of points are given representing different locations in the domain, like external, internal, edge, center and so on. The octree search function should return the correct nearest neighbor for all these points in a fast way. Similarly, all the neighbors within a prescribed range are found out using the search function.



Tests in CellMark/TESTS

Test the full immersed boundary cell marking scheme
step1: import the mesh (fluid) and material points (solid)
step2: build octree structure for fluid mesh cells
step3: for each material points, find out which cell it falls into using fast octree search function
step4: create connectivity for cell to particle and particle to cell
step5: mark the cells which contain particles as Solid type, mark the cells which does not contain any particles as Fluid type
step5:  for Solid cells,  if it has at least Fluid cell as neighbor, mark it as immersed boundary (IB) type
step6: mark the faces shared by a Fluid cell and a IB cell  as IB face
step7: create the connectivity of IB face to particles
step8: create the connectivity of IB face to neighboring fluid cells
step9: store all the connectivity and marking information for further use
In the tests,  three different geometries are tested;  beam, cylinder and sphere. The first two are 2D cases, the last one is 3D.


Tests in Grid/TESTS

CantileverVelocityInterpolation
From experiments, the velocity profile of the cantilever is obtained by measuring on 21 discrete grids on the cantilever surface. This velocity profile will be used as input for the simulation. An interpolation and extrapolation needs to be applied to get the full velocity field on the cantilever mesh using the available 21 grids data. This test is to make sure the interpolation and extrapolation work properly.




Tests in FVMParticleMark/TESTS

 CAVITY_TRI44_PROCS1
 CAVITY_TRI44_PROCS2
 CAVITY_TRI44_PROCS4
 CAVITY_TRI44_PROCS8
 CAVITY_TRI44_PROCS12
 
 CAVITY_TRI26_PROCS1
 CAVITY_TRI26_PROCS2
 CAVITY_TRI26_PROCS4
 CAVITY_TRI26_PROCS6
 CAVITY_TRI26_PROCS8
 
 
 CAVITY_TRI894_PROCS1
 CAVITY_TRI894_PROCS2
 CAVITY_TRI894_PROCS4
 CAVITY_TRI894_PROCS8
 CAVITY_TRI894_PROCS12
 CAVITY_TRI894_PROCS24

The above cases design to test PartMesh class for Triangular Meshes (44, 26 and 894 triangular elements). These meshes is decomposed into 1, 2, 4, 8, 12, 24 submeshes. Specifically, I am testing many small test cases 

Checking following items:
distribution of number of elements
GlobalIndx to hold where first (manual) partitioning starts
Partmetis "eptr" array
Parmetis  "elm" array
Parmetis  "eind" array
Parmetis  "elmwgt" array
Parmetis  "wgtflag"  parameter
Parmetis  "numflag" parameter
Parmetis  "ncon" parameter
Parmetis  "ncommonodes" parameter
Parmetis   "tpwgts" array
Parmetis   "ubvec" array
Parmetis   "option" parameter
Parmetis   "edgecut" parameter"
checking    map between Partition and elements
checking    total elements of partitioning
checking    connectivity returned by Parmetis
checking   total elemnts including ghost cells
checking   cellPartitionin CRConnectivity
checking   Boundary Map (map between boundary ID and boundary type)
checking   faceCells CRConnectivity
checking   globalToLocalMap 
checking   localToGlobalMap 
checking   faceNodes CRConnectivity
checking   cellNodes CRConnectivity
checking   cellCells CRConnectivity
checking   Node coordinates
checking   InterfaceMap( interface ID and its neighbouring mesh_id)
checking   interior face counts
checking   non-interior cells
checking   boundary grouID and offsets
checking   interface ID and offsets
checking   faceCellsOrdered  CRConnectivity (now face ordered such that first
interior then boundary   and interfaces
checking   faceNodesOrdered  CRConnectivity (now face ordered such that first
interior then boundary   and interfaces
checking total number of meshes surrounding local mesh
checking  offset for ghost Cells from adjacent meshes
checking  cells looking interior domain from interfaces (global numbering)
checking  cells looking interior domain from interfaces (local numbering)
checking  mesh formats (ASCII) dumped for Tecplot (visualization program) 
checking  mapping information between meshes if any


The following cases (five of them ) test parallel solver (AMG).
 CAVITY_TRI894_PROCS1_THERMALSOLVER
 CAVITY_TRI894_PROCS2_THERMALSOLVER
 CAVITY_TRI894_PROCS4_THERMALSOLVER
 CAVITY_TRI894_PROCS8_THERMALSOLVER
 CAVITY_TRI894_PROCS12_THERMALSOLVER

Global Mesh size of the above cases is 894 triangular elements for Cavity
problem. This test Algebric Multi Grid Solver (AMG) (communication happens for
each iteration in AMG). The final result is written in tecplot format. The
global mesh has been decomposed in 1, 2, 4, 8, 12 submeshes


The following cases is part of MPM and FVM coupling. We  need to send 
fluid particles position, volume, and stress information. First of all, we 
have to identify fluid particles. The fluid particles are the fluid cells in FVM
mesh just around the immersed boundary cell (IBC). This identification is 
done in FVM_Particles class in FVM code.  We have used two meshes containing
Quad and triangular elements. The class is designed to control number of fluid
particles with "sweep" argument. If you choose sweep "1", the algorithm will
pick only first layer fluid cells around solid body. If you choose "3", you will
get three layer fluid cells around the body. 

 FVMPARTICLE_CAV32_QUAD_SWEEP1
 FVMPARTICLE_CAV32_QUAD_SWEEP2
 FVMPARTICLE_CAV32_QUAD_SWEEP3
 FVMPARTICLE_CAV32_QUAD_SWEEP4
 FVMPARTICLE_CAV32_QUAD_SWEEP5
 FVMPARTICLE_TRI_100_100_SWEEP1
 FVMPARTICLE_TRI_100_100_SWEEP2
 FVMPARTICLE_TRI_100_100_SWEEP3
 FVMPARTICLE_TRI_100_100_SWEEP4
 FVMPARTICLE_TRI_100_100_SWEEP5


