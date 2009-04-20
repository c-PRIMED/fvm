//
// C++ Implementation: PartMesh
//
// Description: 
//
//
// Author: yildirim,,, <yildirim@cfm>, (C) 2008
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include<iostream>
#include <cassert>
#include <fstream>
#include <string>
#include <sstream>
#include <iomanip>
#include <algorithm>
#include "PartMesh.h"
#include "Mesh.h"
#include "Array.h"
#include "CRConnectivity.h"
#include "StorageSite.h"
#include "Vector.h"
#include "OneToOneIndexMap.h"
#include <parmetis.h>



PartMesh::PartMesh( const MeshList &mesh_list, vector<int> nPart,  vector<int> eType ):
_meshList(mesh_list), _nPart(nPart), _eType(eType), _options(0)
{
   if ( !MPI::Is_initialized() )  MPI::Init();
   init();
}


PartMesh::~PartMesh()
{

    //releae memory of vector elements dynamically allocated memory
    vector< Array<int>* > ::iterator it_array;
    for ( it_array = _elemDist.begin(); it_array != _elemDist.end(); it_array++)
       delete *it_array;

    for ( it_array = _globalIndx.begin(); it_array != _globalIndx.end(); it_array++)
       delete *it_array;

    vector< int* > ::iterator it_int;
    for ( it_int = _ePtr.begin(); it_int != _ePtr.end(); it_int++)
       delete [] *it_int;

    for ( it_int = _eInd.begin(); it_int != _eInd.end(); it_int++)
       delete [] *it_int;

    for ( it_int = _eElm.begin(); it_int != _eElm.end(); it_int++)
       delete [] *it_int;


    for ( it_int = _elmWght.begin(); it_int != _elmWght.end(); it_int++)
       delete [] *it_int;

    vector< float* > ::iterator it_float;
    for ( it_float = _tpwgts.begin(); it_float != _tpwgts.end(); it_float++)
       delete [] *it_float;

    for ( it_float = _ubvec.begin(); it_float != _ubvec.end(); it_float++)
       delete [] *it_float;

    for ( it_int = _part.begin(); it_int != _part.end(); it_int++)
       delete [] *it_int;

    for ( it_int = _row.begin(); it_int != _row.end(); it_int++)
       delete [] *it_int;

    for ( it_int = _col.begin(); it_int != _col.end(); it_int++)
       delete [] *it_int;
   
   for ( it_int  = _elem.begin(); it_int != _elem.end(); it_int++)
       delete [] *it_int;

   for ( it_int  = _elemWithGhosts.begin(); it_int != _elemWithGhosts.end(); it_int++)
       delete [] *it_int;

   vector< Mesh* >::iterator it_mesh;
   for ( it_mesh = _meshListLocal.begin(); it_mesh != _meshListLocal.end(); it_mesh++)
        delete *it_mesh;
   
   //if ( !MPI::Is_finalized() )  MPI::Finalize();
}



void 
PartMesh::mesh()
{ 

    CRConnectivity_cellParts();
    CRConnectivity_faceParts();
    interfaces();
    faceCells_faceNodes();
    local_number_elems();
    non_interior_cells();
    order_faceCells_faceNodes();
    coordinates();
    exchange_interface_meshes();
    mesh_setup();
    mappers();

    debug_print();

}


void
PartMesh::partition()
{
 
   compute_elem_dist();
    elem_connectivity();
    parmetis_mesh();
    map_part_elms();
    count_elems_part();
    exchange_part_elems();
  

}


void 
PartMesh::dumpTecplot()
{

//just for mesh 0;
    MPI::COMM_WORLD.Barrier();

    const Mesh& mesh = *(_meshList.at(0));
    const Array<Mesh::VecD3>&  coord = mesh.getNodeCoordinates();

   int tot_nodes = coord.getLength();

    if ( _procID == 0 ) {
       MPI::Status status_row;
       MPI::Status status_col;
       int tag_row = 9;
       int tag_col = 99;
       int tot_elems = _totElems.at(0);
       ofstream part_file( "partiton.dat" );
       
       int *row = new int[tot_elems];  //open big array
       int *col = new int[8*tot_elems]; //maximum hexa 

      
        part_file << "title = \" tecplot file for partionoing \" "  << endl;
        part_file << "variables = \"x\",  \"y\", \"z\", \"partID\" " << endl;
        part_file << "zone N = " << tot_nodes << " E = " << tot_elems <<
           " DATAPACKING = BLOCK,  VARLOCATION = ([4]=CELLCENTERED), ZONETYPE=FETRIANGLE " << endl;

        //x 
        for ( int n = 0; n < tot_nodes; n++){
           part_file << scientific  << coord[n][0] << "     " ;
           if ( n % 5 == 0 ) part_file << endl;
        }
        part_file << endl;
        //y
        for ( int n= 0; n < tot_nodes; n++){
           part_file << scientific << coord[n][1] << "     ";
           if ( n % 5 == 0 ) part_file << endl;
        }
        part_file << endl;
        //z
        for ( int n = 0; n < tot_nodes; n++){
           part_file << scientific << coord[n][2] << "     ";
           if ( n % 5 == 0) part_file << endl;
        }
        part_file <<  endl;

      for ( int elem = 0; elem < _nelems.at(0);  elem++){
               part_file << 0  << "     ";
               if (  elem % 10 == 0 ) part_file << endl;
       }

       vector< shared_ptr< Array<int> > >  row_root;
       vector< shared_ptr< Array<int> > >  col_root;

       //copy row and col values from root
       row_root.push_back( shared_ptr< Array<int> > (new Array<int> (_nelems.at(0)+1) ) );
       col_root.push_back( shared_ptr< Array<int> > (new Array<int> (_colDim.at(0)  ) ) );

       for ( int r = 0; r < _nelems.at(0)+1; r++)
          (*row_root.at(0))[r] = _row.at(0)[r];

       for ( int c = 0; c < _colDim.at(0); c++)
          (*col_root.at(0))[c] = _col.at(0)[c];


       for ( int p = 1; p < _nPart.at(0); p++){
           MPI::COMM_WORLD.Recv( row, tot_elems, MPI::INT, p, tag_row, status_row);
           MPI::COMM_WORLD.Recv( col, 8*tot_elems, MPI::INT, p, tag_col, status_col);

           int nelems = status_row.Get_count(MPI::INT) - 1;
           int col_dim = status_col.Get_count(MPI::INT);
           row_root.push_back( shared_ptr< Array<int> > (new Array<int> (nelems+1) ) );
           col_root.push_back( shared_ptr< Array<int> > (new Array<int> (col_dim)  ) );
           //fill row_root
           for ( int r = 0; r < nelems+1; r++)
              (*row_root.at(p))[r] = row[r];
           //fill col_root
           for ( int c = 0; c < col_dim; c++)
              (*col_root.at(p))[c] = col[c];

           for ( int elem = 0; elem < nelems;  elem++){
               part_file << p  << "     ";
               if (  elem % 10 == 0 ) part_file << endl;
           }
       }
      part_file << endl;


      for ( int p = 0; p < _nPart.at(0); p++){
         //connectivity
        int nelems = row_root.at(p)->getLength() - 1;
         for ( int elem = 0; elem < nelems;  elem++){
             int node_start = (*row_root.at(p))[elem];
             int node_end   = (*row_root.at(p))[elem+1];
             for ( int node = node_start; node < node_end; node++)
                 part_file << setw(6) <<  (*col_root.at(p))[node]+1 << "     ";
              part_file <<endl;
         }
       }

       part_file.close();
       delete [] row;
       delete [] col;

    } else { 
       int tag_row = 9;
       int tag_col = 99;

        MPI::COMM_WORLD.Send(_row.at(0), _nelems.at(0)+1, MPI::INT, 0, tag_row);
        MPI::COMM_WORLD.Send(_col.at(0), _colDim.at(0)  , MPI::INT, 0, tag_col);

    }   
      

 


}




void 
PartMesh::debug_print()
{

  stringstream ss;
  ss << "proc" << _procID << "_debug_print.dat";
  ofstream  debug_file( (ss.str()).c_str() );

//    if ( _procID == 0 ){
       for ( int id = 0; id < _nmesh; id++){
            debug_file << " procID = " << _procID << endl;
            debug_file << " npart  = " << _nPart.at(id) << endl;

            debug_file << endl;

            //elemDist
            for ( int n = 0; n < _nPart.at(id); n++)
               debug_file << " elemDist[" << n << "] = " << (*_elemDist.at(id))[n] << endl;

            debug_file << endl;

           //globalIndx
            for ( int n = 0; n <= _nPart.at(id); n++)
               debug_file << " n = " << n << " globalIndx[" << n << "] = " << (*_globalIndx.at(id))[n]  << endl;

            debug_file << endl;

          //eptr 
          int mesh_nlocal = (*_elemDist.at(id))[_procID];
          for ( int i = 0; i <= mesh_nlocal; i++)
              debug_file << " eptr[" << i << "] = " << _ePtr.at(id)[i] << endl;

            debug_file << endl;
         
         //eelm 
          for ( int i = 0; i < mesh_nlocal; i++)
              debug_file << " eelm[" << i << "] = " << _eElm.at(id)[i] << endl;

            debug_file << endl;
 
          //eind
          int nlocal_elem = (*_elemDist.at(id))[_procID];
          int indx = 0;
          for ( int i = 0; i < nlocal_elem; i++){
              int nnodes = _ePtr.at(id)[i+1] - _ePtr.at(id)[i];
             debug_file << " elemID  = " << i << ",  ";
              for ( int j = 0; j < nnodes; j++){
                 debug_file << " eind[" << indx << "]=" << _eInd.at(id)[indx] <<  "   ";
                 indx++;
              }
              debug_file << endl;
         }

        debug_file << endl;

         //elmwgt
          for ( int i = 0; i < nlocal_elem; i++){
                 debug_file << " elmwgt[" << i << "]=" << _elmWght.at(id)[i] <<  endl;
         }

         debug_file << endl;

         //wghtflag
          debug_file << " wgtflag = " << _wghtFlag.at(id) << endl;

          debug_file << endl;

         //numflag
        debug_file << " numflag  = " << _numFlag.at(id) << endl;
          
         debug_file << endl;

         //ncon
        debug_file << " ncon = " << _ncon.at(id) << endl;
         
        debug_file << endl;
        //ncommonnodes
        debug_file << " ncommonnodes = " << _ncommonNodes.at(id) << endl;
        debug_file << endl;         //nparts
        debug_file << " nparts = " << _nPart.at(id) << endl;
        debug_file << endl;
         //tpwgts
         for ( int i = 0; i < _nPart.at(id); i++)
            debug_file << "tpwgts[" << i << "] = " << _tpwgts.at(id)[i] << endl;
        debug_file << endl;        //ubvec
        for ( int i = 0; i < _ncon.at(id); i++)
           debug_file << " ubvec = " << _ubvec.at(id)[i] << endl;
        debug_file << endl;
       //options
        debug_file << " options = " << _options << endl;
        debug_file << endl;
       //edgecut (output)
        debug_file << " edgecut = " << _edgecut.at(id) << endl;
        debug_file << endl;
       //part(output)
       for ( int i = 0; i < nlocal_elem; i++){
          debug_file << " elem = " << (*_globalIndx.at(id))[_procID] + i << " partion = " <<  _part.at(id)[i] << endl;
       }
       debug_file << endl;
       for ( int p = 0; p < _nPart.at(id); p++){
          multimap<int,int>::iterator it, itlow, itup;
          itlow = _mapPartAndElms.at(id).lower_bound(p);
          itup  = _mapPartAndElms.at(id).upper_bound(p);
          for( it = itlow; it != itup; it++)
             debug_file << " partID = " << it->first << " elemID = " << it->second << endl;
       }
       debug_file << endl;
        //total elems
       debug_file << " total elements  = " << _nelems.at(id)  <<  endl;

       //total dimension for col array
       debug_file << " total dim of col = " << _colDim.at(id) << endl;
       debug_file << endl;
       //row
       for ( int n = 0; n < _nelems.at(id)+1; n++){
           debug_file << " _row[" << n << "] = " << _row.at(id)[n] << endl;
       }
       debug_file << endl;
      
       //elem
      for ( int n = 0; n < _nelems.at(id); n++){
          debug_file << " _elem[" << n << "] = " << _elem.at(id)[n] << endl;
      }
      debug_file << endl;

       //elemWithGhosts
      for ( int n = 0; n < _nelemsWithGhosts.at(id); n++){
          debug_file << " _elemWithGhosts[" << n << "] = " << _elemWithGhosts.at(id)[n] << endl;
      }
      debug_file << endl;


       //col
       for ( int n = 0; n < _colDim.at(id); n++){
           debug_file << " _col[" << n << "] = " << _col.at(id)[n] << endl;
       }

       //cellParts
       debug_file << " _cellParts : " << endl;
       debug_file << " _cellParts->getRowDim() = " << _cellParts.at(id)->getRowDim() << endl;
       debug_file << " _cellParts->getColDim() = " << _cellParts.at(id)->getColDim() << endl;
       const Array<int>&   rowCellParts = _cellParts.at(id)->getRow();
       const Array<int>&   colCellParts = _cellParts.at(id)->getCol();

       for ( int cell = 0; cell < _cellParts.at(id)->getRowDim(); cell++){
          debug_file << " row[" << cell<<"] = " << rowCellParts[cell] << "    ";
          int nnodes = rowCellParts[cell+1] - rowCellParts[cell];
          for ( int node = 0; node < nnodes; node++){
             debug_file << colCellParts[ rowCellParts[cell] + node ] << "    ";
          }
          debug_file << endl;
       }


       debug_file << endl;
       multimap<int,int>::iterator it_multimap;
       for ( it_multimap  = _mapBounIDAndCell.at(id).begin(); 
             it_multimap != _mapBounIDAndCell.at(id).end(); it_multimap++)
           debug_file << "Boundary multimap = " << it_multimap->first << "    " << it_multimap->second << endl;
        multimap<int,string>::iterator it_multimapS;
        for ( it_multimapS  = _mapBounIDAndBounType.at(id).begin(); 
             it_multimapS != _mapBounIDAndBounType.at(id).end(); it_multimapS++)
           debug_file << "Boundary multimap = " << it_multimapS->first << "    " << it_multimapS->second << endl;
  
        debug_file << endl;
       //faceParts
       debug_file << " _faceParts : " << endl;
       debug_file << " _faceParts->getRowDim() = " << _faceParts.at(id)->getRowDim() << endl;
       debug_file << " _faceParts->getColDim() = " << _faceParts.at(id)->getColDim() << endl;
       const Array<int>&   rowFaceParts = _faceParts.at(id)->getRow();
       const Array<int>&   colFaceParts = _faceParts.at(id)->getCol();

       for ( int cell = 0; cell < _faceParts.at(id)->getRowDim(); cell++){
          debug_file << " row[" << cell<<"] = " << rowFaceParts[cell] << "    ";
          int nnodes = rowFaceParts[cell+1] - rowFaceParts[cell];
          for ( int node = 0; node < nnodes; node++){
             debug_file << colFaceParts[ rowFaceParts[cell] + node ] << "    ";
          }
          debug_file << endl;
       }
       
       debug_file << endl;

      //faceCells 
      debug_file << " _faceCells :  " << endl;
      debug_file << " _faceCells->getRowDim() = " << _faceCells.at(id)->getRowDim() << endl;
      debug_file << " _faceCells->getColDim() = " << _faceCells.at(id)->getColDim() << endl;
      const Array<int>&   rowFaceCells = _faceCells.at(id)->getRow();
      const Array<int>&   colFaceCells = _faceCells.at(id)->getCol();
      const Array<int>&   globalToLocalMap = _faceCells.at(id)->getGlobalToLocalMap();
      const Array<int>&   localToGlobalMap = _faceCells.at(id)->getLocalToGlobalMap();
      debug_file << " globalToLocalMap.length() = " << globalToLocalMap.getLength() << endl;
      
      for ( int n = 0; n < globalToLocalMap.getLength(); n++)
          debug_file << " globalToLocalMap[" << n << "] = " << globalToLocalMap[n] << endl;
      debug_file << endl;
      debug_file << " localToGlobalMap.length() = " << localToGlobalMap.getLength() << endl;
      for ( int n = 0; n < localToGlobalMap.getLength(); n++)
          debug_file << " localToGlobalMap[" << n << "] = " << localToGlobalMap[n] << endl;
        



       for ( int face = 0; face < _faceCells.at(id)->getRowDim(); face++){
          debug_file << " row[" << face+1 <<"] = " << (*_partFaces.at(id))(_procID,face)+1 << "    ";
          int ncells = _faceCells.at(id)->getCount(face);
          for ( int cell = 0; cell < ncells; cell++){
             debug_file << colFaceCells[ rowFaceCells[face] + cell ] + 1 << "    ";
          }
          debug_file << endl;
       }

       debug_file << endl;

      //faceNodes
      debug_file << " _faceNodes :  " << endl;
      debug_file << " _faceNodes->getRowDim() = " << _faceNodes.at(id)->getRowDim() << endl;
      debug_file << " _faceNodes->getColDim() = " << _faceNodes.at(id)->getColDim() << endl;
      const Array<int>&   rowFaceNodes = _faceNodes.at(id)->getRow();
      const Array<int>&   colFaceNodes = _faceNodes.at(id)->getCol();

       for ( int face = 0; face < _faceNodes.at(id)->getRowDim(); face++){
          debug_file << " row[" << face+1<<"] = " << (*_partFaces.at(id))(_procID,face)+1 << "    ";
          int nnodes = _faceNodes.at(id)->getCount(face);
          for ( int node = 0; node < nnodes; node++){
             debug_file << colFaceNodes[ rowFaceNodes[face] + node ]+1 << "    ";
          }
          debug_file << endl;
       }
       
       debug_file << endl;


      //faceNodes
      debug_file << " _cellNodes(Local Numbering) :  " << endl;
      debug_file << " _cellNodes->getRowDim() = " << _cellNodes.at(id)->getRowDim() << endl;
      debug_file << " _cellNodes->getColDim() = " << _cellNodes.at(id)->getColDim() << endl;
      const Array<int>&   rowCellNodes = _cellNodes.at(id)->getRow();
      const Array<int>&   colCellNodes = _cellNodes.at(id)->getCol();

       for ( int cell = 0; cell < _cellNodes.at(id)->getRowDim(); cell++){
          debug_file << " row[" << cell+1 << "]  = " ;
          int nnodes = _cellNodes.at(id)->getCount(cell);
          for ( int node = 0; node < nnodes; node++){
             debug_file << colCellNodes[ rowCellNodes[cell] + node ]+1 << "    ";
          }
          debug_file << endl;
       }
       
       debug_file << endl;


      //cellCells
      debug_file << " _cellCells :  " << endl;
      debug_file << " _cellCells->getRowDim() = " << _cellCells.at(id)->getRowDim() << endl;
      debug_file << " _cellCells->getColDim() = " << _cellCells.at(id)->getColDim() << endl;
      const Array<int>&   rowCellCells = _cellCells.at(id)->getRow();
      const Array<int>&   colCellCells = _cellCells.at(id)->getCol();

       for ( int cell = 0; cell < _cellCells.at(id)->getRowDim(); cell++){
          debug_file << " row[" << cell+1<<"] = "  << "    ";
          int nnodes = _cellCells.at(id)->getCount(cell);
          for ( int node = 0; node < nnodes; node++){
             debug_file << colCellCells[ rowCellCells[cell] + node ]+1 << "    ";
          }
          debug_file << endl;
       }
       
       debug_file << endl;

       //node coordinates
        int node_count = _partNodes.at(id)->getCount( _procID );
       // const Array<int>&    rowPartNodes = _partNodes.at(id)->getRow();
        //const Array<int>&    colPartNodes = _partNodes.at(id)->getCol();
        for ( int node = 0; node < node_count; node++){
 //           int nodeID = colPartNodes[rowPartNodes[_procID]+node];
            debug_file << fixed;
            debug_file << " node ID = " <<setw(10)<< node+1 << setprecision(7) << ",  x = " << (*_coord.at(id))[node][0] << 
                          setprecision(7) << ",  y = " << (*_coord.at(id))[node][1] << 
                          setprecision(7) << ",  z = " << (*_coord.at(id))[node][2] << endl; 
       }
       
       debug_file << endl;



       //interfaceMap
       pair<multimap<int,int>::iterator,multimap<int,int>::iterator> ret;
       debug_file << " _interfaceMap.size() = " <<  _interfaceMap.at(id).size() << endl;
       for ( int part = 0; part < _nPart.at(id); part++){
//           int tot_interfaces = _interfaceMap.at(id).count(part);
             
              ret = _interfaceMap.at(id).equal_range(part);
              debug_file  << " interface ID =  "  << part << "  =>  ";
              for (it_multimap=ret.first; it_multimap!=ret.second; ++it_multimap)
                 debug_file  <<  (*_partFaces.at(id))(_procID,  (*it_multimap).second) + 1  << "  ";
              debug_file << endl;
       }
       debug_file << endl;
       
       //interior face counts
       debug_file << " interior face counts = " << count_interior_faces( id ) << endl;
       debug_file << endl;

     
      //interior cells
      set<int>::const_iterator it_set;
      debug_file << " total interior cells = " << _elemLocal.at(id).size() << endl;
      for ( it_set = _elemLocal.at(id).begin(); it_set != _elemLocal.at(id).end(); it_set++ )
          debug_file <<  "      " <<  *it_set  << endl;

      debug_file << endl;

       //non-interior cells (only for last mesh)
       debug_file << " total non-interior cells = " << _nonInteriorCells.at(id).size() << endl;
       for ( it_set = _nonInteriorCells.at(id).begin(); it_set != _nonInteriorCells.at(id).end(); it_set++ )
          debug_file <<  "      " <<  *it_set  << endl;

      debug_file << endl;
       
        for ( it_multimap = _bndryOffsets.at(id).begin(); it_multimap != _bndryOffsets.at(id).end(); it_multimap++ )
           debug_file  << "   bndry group ID = " << it_multimap->first << " offsets = " << it_multimap->second << endl;

      debug_file << endl;
       
        for ( it_multimap = _interfaceOffsets.at(id).begin(); it_multimap != _interfaceOffsets.at(id).end(); it_multimap++ )
           debug_file  << "   interface ID = " << it_multimap->first << " offsets = " << it_multimap->second << endl;

        debug_file << endl;
    
      //faceCellsOrdered
      debug_file << " _faceCellsOrdered :  " << endl;
      debug_file << " _faceCellsOrdered->getRowDim() = " << _faceCellsOrdered.at(id)->getRowDim() << endl;
      debug_file << " _faceCellsOrdered->getColDim() = " << _faceCellsOrdered.at(id)->getColDim() << endl;
      const Array<int>& rowFaceCellsOrdered = _faceCellsOrdered.at(id)->getRow();
      const Array<int>& colFaceCellsOrdered = _faceCellsOrdered.at(id)->getCol();
      
       for ( int face = 0; face < _faceCellsOrdered.at(id)->getRowDim(); face++){
          debug_file << " row[" << face <<"] = " ;
          int ncells = _faceCellsOrdered.at(id)->getCount(face);
          for ( int cell = 0; cell < ncells; cell++){
             debug_file << colFaceCellsOrdered[ rowFaceCellsOrdered[face] + cell ] + 1 << "    ";
          }
          debug_file << endl;
       }

       debug_file << endl;

      //faceNodes
      debug_file << " _faceNodesOrdered :  " << endl;
      debug_file << " _faceNodesOrdered->getRowDim() = " << _faceNodesOrdered.at(id)->getRowDim() << endl;
      debug_file << " _faceNodesOrdered->getColDim() = " << _faceNodesOrdered.at(id)->getColDim() << endl;
      const Array<int>& rowFaceNodesOrdered = _faceNodesOrdered.at(id)->getRow();
      const Array<int>& colFaceNodesOrdered = _faceNodesOrdered.at(id)->getCol();

      for ( int face = 0; face < _faceNodesOrdered.at(id)->getRowDim(); face++){
          debug_file << " row[" << face<<"] = " ;
          int nnodes = _faceNodesOrdered.at(id)->getCount(face);
          debug_file << nnodes << "        ";
          for ( int node = 0; node < nnodes; node++){
             debug_file << colFaceNodesOrdered[ rowFaceNodesOrdered[face] + node ]+1 << "    ";
          }
          debug_file << endl;
      }
       
       debug_file << endl;


       //interfaceMesheCounts
      for ( int proc = 0; proc < _nPart.at(id); proc++)
          debug_file << " total mesh surrounding = " << (*_interfaceMeshCounts.at(id))[proc] << endl;

       debug_file << endl;

      //ofsest
       debug_file << " offset for ghost Cells from adjacent meshes to read data from _ghostCellsGlobal : "  << endl;
       for ( int n = 0; n < _offsetInterfaceCells.at(id)->getLength(); n++ )
          debug_file << "    n  =  " << n << " offsetInterfaceCells = " << (*_offsetInterfaceCells.at(id))[n] << endl;
    
       debug_file << endl;
     //interfaceMeshIDs
       debug_file << " neightboorhood cell IDs : "  << endl;
       for ( int n = 0; n < _interfaceMeshIDs.at(id)->getLength(); n++ )
          debug_file << "    n  =  " << n << "  interfaced Mesh ID = " <<  (*_interfaceMeshIDs.at(id))[n] << endl;

       debug_file << endl;

     //global Interface cells (interior ones, global numbering)
       debug_file << "interface cells looking interior domain (global numbering)  : " << endl;
       for ( int n = 0; n < _ghostCellsGlobal.at(id)->getLength(); n++ )
          debug_file << "    n  =  " << n << "  cell ID = " <<  (*_ghostCellsGlobal.at(id))[n] << endl;

     //global Interface cells (interior ones, global numbering)
       debug_file << "interface cells looking interior domain (local numbering)  : " << endl;
       for ( int n = 0; n < _ghostCellsLocal.at(id)->getLength(); n++ )
          debug_file << "    n  =  " << n << "  interfaced Mesh ID = " <<  (*_ghostCellsLocal.at(id))[n] << endl;



    }

       

   


 

//  }

  debug_file.close();
}


             //SET PROPERTIES METHODS
void 
PartMesh::setWeightType( int weight_type )
{
   for ( int id = 0; id < _nmesh; id++)
      _wghtFlag.at(id) =  weight_type;

}

void 
PartMesh::setNumFlag( int num_flag )
{
  for ( int id = 0; id < _nmesh; id++)
      _numFlag.at(id) = num_flag;

}

           // PRIVATE METHODS

void
PartMesh::init()
{

    _procID = MPI::COMM_WORLD.Get_rank();
    _nmesh  = _meshList.size();

   _totElems.resize  ( _nmesh );
   _totElemsAndGhosts.resize( _nmesh );
   _wghtFlag.resize( _nmesh );
   _numFlag.resize ( _nmesh );
   _ncon.resize    ( _nmesh ); 
   _ncommonNodes.resize( _nmesh );
   _mapPartAndElms.resize( _nmesh );
   _boundarySet.resize( _nmesh );
   _interfaceSet.resize( _nmesh );
   _mapBounIDAndCell.resize( _nmesh );
   _mapBounIDAndBounType.resize( _nmesh );
   _nelems.resize( _nmesh );
   _nelemsWithGhosts.resize( _nmesh );
   _colDim.resize( _nmesh );
   _edgecut.resize( _nmesh );
   _elemLocal.resize(_nmesh);  //local numbering
   _nonInteriorCells.resize(_nmesh ); //local numbering
   _bndryOffsets.resize( _nmesh );
   _interfaceOffsets.resize( _nmesh );
   _cellToOrderedCell.resize( _nmesh );
   _isOneToOneMap.resize( _nmesh );
   _globalToLocalMappers.resize( _nmesh );
   _localToGlobalMappers.resize( _nmesh );
   _windowSize.resize( _nmesh );
   _fromIndices.resize( _nmesh );
   _toIndices.resize( _nmesh );

    for ( int id = 0; id < _nmesh; id++){
        StorageSite& site = _meshList[id]->getCells();
       _totElems.at(id)   = site.getSelfCount();
       _totElemsAndGhosts.at(id) = site.getCount();
       _wghtFlag.at(id) = int( NOWEIGHTS ); //No Weights : default value
       _numFlag.at(id)  = int( C_STYLE );   //C Style numbering :: default_value
       _ncon.at(id)     =  1;               //number of specified weights : default value

       //assign ubvec
       _ubvec.push_back( new float[_ncon.at(id) ] );
       for ( int n = 0; n < _ncon.at(id); n++)
         _ubvec.at(id)[n] = 1.05f; //1.05 suggested value from parMetis manual


        //assign elementy type
        switch (_eType.at(id) ){
            case   TRI :
                   _ncommonNodes.at(id) = 2;
                   break;
            case   TETRA :
                   _ncommonNodes.at(id) = 3;
                   break;
            case   HEXA :
                   _ncommonNodes.at(id) = 4;
                   cout << " NO SUPPORT FOR NOW " << endl;
                    abort();
                   break;                
            case   QUAD :
                   _ncommonNodes.at(id) = 2;
                   cout << " NO SUPPORT FOR NOW " << endl;
                   abort();
                    break;
            default :
                    cout << " ONLY TRIANGLE, TETRAHEDRAL, HEXAHEDRAL and QUADRILATERAL elements must be chose " <<
                    endl;     abort(); 
        } 

        //get tpwgts
        int ncon_by_nparts = _ncon.at(id) * _nPart.at(id);
        _tpwgts.push_back( new float[ ncon_by_nparts ] );
        for (int n = 0; n < ncon_by_nparts; n++){
            _tpwgts.at(id)[n] = 1.0f / float( ncon_by_nparts );
        }

       //edgecut
        _edgecut.at(id) = -1;

      _interfaceMeshCounts.push_back (  ArrayIntPtr( new Array<int>(_nPart.at(id)) )  );
      _procTotalInterfaces.push_back (  ArrayIntPtr( new Array<int>(_nPart.at(id)) )  );

    }


}



void 
PartMesh::compute_elem_dist()
{

   for (int id = 0; id < _nmesh; id++){
      
      int nelems = _totElems[id];
      int npart  = _nPart[id];
      int nremainder = nelems % npart;
      int init_dist = (nelems - nremainder) / npart;
      _elemDist.push_back( new Array<int>( npart) );
      _globalIndx.push_back( new Array<int>(npart+1) );
      *_elemDist.at(id) = init_dist;

       int p = 0;
       while ( nremainder != 0 ){
           (*_elemDist.at(id))[p % npart]++;
            nremainder--;
            p++;
       }

       (*_globalIndx[id]) = 0;
       int sum_elem = 0;
       for ( int n = 1; n <= npart; n++){
           sum_elem += (*_elemDist.at(id))[n-1];
           ((*_globalIndx.at(id))[n]) =  sum_elem;
       }
   }

}

void
PartMesh::elem_connectivity()
{

   for (int id = 0; id < _nmesh; id++){
       //allocate local ePtr for parMetis
       int mesh_nlocal = (*_elemDist.at(id))[_procID];
      _ePtr.push_back( new int[mesh_nlocal+1] );
      _eElm.push_back( new int[mesh_nlocal+1] );
       //element weights  
      _elmWght.push_back( new int[mesh_nlocal] );
       for ( int n = 0; n < mesh_nlocal; n++)
          _elmWght.at(id)[n] = 1;

       //allocate local eInd for ParMETIS
      _eInd.push_back( new int[local_nodes(id)] );

 
//       //setting ePtr and eInd for ParMETIS
       set_eptr_eind(id);

     _part.push_back( new  int[mesh_nlocal] );
      for ( int n = 0; n < mesh_nlocal; n++)
         _part.at(id)[n] = -1;
   }

}

int
PartMesh::local_nodes( int id )
{

      const Mesh* mesh = _meshList[id];
      const CRConnectivity& cellNodes = mesh->getCellNodes();

      //get local nodes
      int local_nodes = 0;
      int nstart = (*_globalIndx.at(id))[_procID];
      int npart  = nstart + (*_elemDist[id])[_procID];
      for ( int n = nstart; n < npart; n++)
           local_nodes += cellNodes.getCount(n);

       return local_nodes;
}

void
PartMesh::set_eptr_eind( int id )
{
      const Mesh* mesh = _meshList[id];
      const CRConnectivity& cellNodes = mesh->getCellNodes();

      int elem_start   = (*_globalIndx.at(id))[_procID];
      int elem_finish  = elem_start + (*_elemDist[id])[_procID];
      int indxInd   = 0;
      int indxPtr   = 0;
      _ePtr.at(id)[indxPtr] = 0;
      for ( int elem = elem_start; elem < elem_finish; elem++ ){
          _eElm.at(id)[indxPtr] = elem;
           indxPtr++;
          _ePtr.at(id)[indxPtr] = _ePtr.at(id)[indxPtr-1] +  cellNodes.getCount(elem);
          for (  int node = 0; node < cellNodes.getCount(elem); node++)
             _eInd.at(id)[indxInd++] = cellNodes(elem,node);
      }

}


void
PartMesh::parmetis_mesh()
{
   MPI_Comm comm_world = MPI::COMM_WORLD;
   for ( int id = 0; id < _nmesh; id++){
    
       ParMETIS_V3_PartMeshKway( &(*_globalIndx.at(id))[0], _ePtr.at(id), _eInd.at(id),
        _elmWght.at(id), &_wghtFlag.at(id), &_numFlag.at(id), &_ncon.at(id),  &_ncommonNodes.at(id),
        &_nPart.at(id), _tpwgts.at(id), _ubvec.at(id), &_options, &_edgecut.at(id), _part.at(id), &comm_world );

   } 
}   



void 
PartMesh::map_part_elms()
{

   for (int  id = 0; id < _nmesh; id++){
      int nlocal_elem = (*_elemDist.at(id))[_procID];

      for ( int elm = 0; elm < nlocal_elem; elm++){
         int partID = _part.at(id)[elm];
         _mapPartAndElms.at(id).insert(pair<int,int>(partID,elm));
      }
  
   }
}

void 
PartMesh::count_elems_part()
{
   for ( int id = 0; id < _nmesh; id++){
       for ( int partID = 0; partID < _nPart.at(id); partID++){
          int nelems_local = _mapPartAndElms.at(id).count(partID);

          //sum to find size of ncol
          multimap<int,int>::const_iterator it = _mapPartAndElms.at(id).find(partID);
          multimap<int,int>::const_iterator itlow = _mapPartAndElms.at(id).lower_bound(partID);
          multimap<int,int>::const_iterator itup = _mapPartAndElms.at(id).upper_bound(partID);
          int ncol_local = 0;
          for ( it = itlow; it != itup; it++){
              int pos = it->second;  // element number
              ncol_local += _ePtr.at(id)[pos+1] - _ePtr.at(id)[pos];
          }

          MPI::COMM_WORLD.Reduce(&nelems_local, &_nelems.at(id), 1, MPI::INT, MPI::SUM, partID);
          MPI::COMM_WORLD.Reduce(&ncol_local, &_colDim.at(id), 1, MPI::INT, MPI::SUM, partID);

       }
       //now each processor now how many elements and nodes
      _row.push_back ( new int[_nelems.at(id)+1] );
      _elem.push_back( new int[_nelems.at(id) ] );
      _col.push_back ( new int[_colDim.at(id)]   );
       for ( int n = 0; n < _nelems.at(id)+1; n++ )
            _row.at(id)[n] = -1;
       for ( int n = 0; n < _colDim.at(id); n++)
            _col.at(id)[n] = -1;
       for ( int n = 0; n < _nelems.at(id); n++)
            _elem.at(id)[n] = -1;

   }


}

void
PartMesh::exchange_part_elems()
{
  
   for ( int id = 0; id < _nmesh; id++){

       int *countsRow  =  new int[_nPart.at(id)];
       int *countsCol  =  new int[_nPart.at(id)];
       int *offsetsRow =  new int[_nPart.at(id)];
       int *offsetsCol =  new int[_nPart.at(id)];


       for ( int partID = 0; partID < _nPart.at(id); partID++){
          int nelems_local = _mapPartAndElms.at(id).count(partID);
          int *row_local  = new int[nelems_local];
          int *elem_local = new int[nelems_local];

           multimap<int,int>::const_iterator it = _mapPartAndElms.at(id).find(partID);
           multimap<int,int>::const_iterator itlow = _mapPartAndElms.at(id).lower_bound(partID);
           multimap<int,int>::const_iterator itup = _mapPartAndElms.at(id).upper_bound(partID);

           //fill row array
           int indx = 0;
           int ncol_local = 0; 
           for ( it = itlow; it != itup; it++){
              int pos         = it->second;  // element number
              ncol_local     += _ePtr.at(id)[pos+1] - _ePtr.at(id)[pos];
              row_local[indx]  = _ePtr.at(id)[pos+1]-_ePtr.at(id)[pos]; //aggregation before shipping
              elem_local[indx] = _eElm.at(id)[pos];
              indx++;
           }

           //fill col array
           int *col_local = new int[ncol_local];
           indx = 0;
           for ( it = itlow; it != itup; it++ ){
              int elID = it->second;
              int node_start = _ePtr.at(id)[elID];
              int node_end   = _ePtr.at(id)[elID+1];
              for ( int node = node_start; node < node_end; node++){
                 col_local[indx++] = _eInd.at(id)[node];
              }
           }

         //forming counts
         int nrow_local = nelems_local;
         MPI::COMM_WORLD.Allgather(&nrow_local, 1, MPI::INT, countsRow, 1, MPI::INT);
         MPI::COMM_WORLD.Allgather(&ncol_local, 1, MPI::INT, countsCol, 1, MPI::INT);


        //form offsets
        offsetsRow[0]  = 0;
        offsetsCol[0]  = 0;
        for ( int p = 1; p < int(_nPart.at(id)); p++ ){
           offsetsRow[p]  = countsRow[p-1]  + offsetsRow[p-1];
           offsetsCol[p]  = countsCol[p-1]  + offsetsCol[p-1];
        } 

        //gathering partial partions for _row and _col
        MPI::COMM_WORLD.Gatherv(row_local, countsRow[_procID], MPI::INT, _row.at(id), 
                                countsRow, offsetsRow, MPI::INT, partID);

        MPI::COMM_WORLD.Gatherv(col_local, countsCol[_procID], MPI::INT, _col.at(id), 
                                countsCol, offsetsCol, MPI::INT, partID);
 
       MPI::COMM_WORLD.Gatherv(elem_local, countsRow[_procID], MPI::INT, _elem.at(id),
                                countsRow, offsetsRow, MPI::INT, partID);


        delete [] row_local;
        delete [] col_local;
        delete [] elem_local;

      
       }  // for::partID

       delete [] countsRow ;
       delete [] countsCol ;
       delete [] offsetsRow;
       delete [] offsetsCol;

    }  // for::meshID


   shift_sum_row();

}


void
PartMesh::shift_sum_row()
{
   for ( int id = 0; id < _nmesh; id++){
       //shift [0,n] to [1,n+1]
       for ( int n = _nelems.at(id); n > 0; n--)
           _row.at(id)[n] = _row.at(id)[n-1];
      _row.at(id)[0] = 0;
      //summing row ex: row = {0,3,6,9,...} for triangle (three nodes)
       for ( int n = 1; n < _nelems.at(id)+1; n++ )
          _row.at(id)[n] += _row.at(id)[n-1];

   }


}



void 
PartMesh::mesh_setup()
{
    for ( int id = 0; id < _nmesh; id++){
        //interior faces
        _meshListLocal.at(id)->createInteriorFaceGroup( count_interior_faces(id) );

         //boundary faces
        set<int>::const_iterator it_set;
        for ( it_set = _boundarySet.at(id).begin(); it_set != _boundarySet.at(id).end(); it_set++){
           int bndryID = *it_set;
           int size   =  _mapBounIDAndCell.at(id).count(bndryID);
           if ( size > 0 ){
              int offset =  _bndryOffsets.at(id)[ bndryID ] ;
              string boundaryType = _mapBounIDAndBounType.at(id)[bndryID];
              _meshListLocal.at(id)->createBoundaryFaceGroup( size, offset, bndryID, boundaryType);
           }
         }


        //then interface faces
        for ( it_set = _interfaceSet.at(id).begin(); it_set != _interfaceSet.at(id).end(); it_set++){
           int interfaceID = *it_set;
           int size = int(_interfaceMap.at(id).count( interfaceID ) );
           int offset = _interfaceOffsets.at(id)[ interfaceID ];
           _meshListLocal.at(id)->createInterfaceGroup( size, offset, interfaceID );
            shared_ptr<StorageSite> site( new StorageSite(size) );
            site->setScatterProcID( _procID );
            site->setGatherProcID ( interfaceID );
           _meshListLocal.at(id)->createGhostCellSite( interfaceID,  site );
         }


       _meshListLocal.at(id)->setCoordinates( _coord.at(id) );
       _meshListLocal.at(id)->setFaceNodes  ( _faceNodesOrdered.at(id) );
       _meshListLocal.at(id)->setFaceCells  ( _faceCellsOrdered.at(id) );


    }

}


//construct CRConnectivity faceParts
void
PartMesh::CRConnectivity_faceParts()
{

     for ( int id = 0; id < _nmesh; id++){
          _faceCellsGlobal.push_back( &_meshList.at(id)->getAllFaceCells() );
          _faceNodesGlobal.push_back( &_meshList.at(id)->getAllFaceNodes() );
          _faceParts.push_back( _faceCellsGlobal.at(id)->multiply( *_cellParts.at(id), false) );
          _partFaces.push_back( _faceParts.at(id)->getTranspose() );
          _partNodes.push_back( _partFaces.at(id)->multiply( *_faceNodesGlobal.at(id), false) );
    }

}


//construct CRConnectivity cellParts
void
PartMesh::CRConnectivity_cellParts()
{
    vector< int* >  elemGlobal;
    vector< int* >  distGlobal;  //total partition + 1 suche that 0, 5, 10, 15 

    for ( int id = 0; id < _nmesh; id++){
        mapBounIDAndCell(id);
        resize_elem(id);

        elemGlobal.push_back( new int [_totElemsAndGhosts.at(id) ]   ); //global array to aggregation
        distGlobal.push_back( new int [ _nPart.at(id) + 1 ] );
        int *offsets = new int [ _nPart.at(id) ];

        MPI::COMM_WORLD.Allgather(&_nelemsWithGhosts.at(id), 1, MPI::INT, distGlobal.at(id), 1, MPI::INT);
        //form offsets
        offsets[0]  = 0;
        for ( int p = 1; p < int(_nPart.at(id)); p++ ){
           offsets[p]  = distGlobal.at(id)[p-1]  + offsets[p-1];
        } 

        //gathering partial partions for _row and _col
        MPI::COMM_WORLD.Allgatherv(_elemWithGhosts.at(id), _nelemsWithGhosts.at(id), MPI::INT, elemGlobal.at(id), 
                                distGlobal.at(id), offsets, MPI::INT);

       //shift  distGlobal to one right
        for ( int i = int(_nPart.at(id)); i > 0; i--)
           distGlobal.at(id)[i] = distGlobal.at(id)[i-1];

        distGlobal.at(id)[0] = 0;

       //summing distGlobal
        for ( int i = 1; i < _nPart.at(id)+1; i++ ){
             distGlobal.at(id)[i] += distGlobal.at(id)[i-1];
        }
        delete [] offsets;

      //forming CRConnectivity for cellPart
      _cellSiteGlobal.push_back( StorageSitePtr(new StorageSite(_totElemsAndGhosts.at(id) )) );
      _partSite.push_back( StorageSitePtr(new StorageSite(_nPart.at(id)) )  );

       _cellParts.push_back( CRConnectivityPtr(new CRConnectivity( *_cellSiteGlobal.at(id), *_partSite.at(id)) ) );

      _cellParts.at(id)->initCount();

      for ( int indx = 0; indx < _totElemsAndGhosts.at(id); indx++)
          _cellParts.at(id)->addCount(indx,1);

      _cellParts.at(id)->finishCount();

      int index = 0;
      while ( index < _nPart.at(id) ){
         for ( int n = distGlobal.at(id)[index]; n < distGlobal.at(id)[index+1]; n++){

             _cellParts.at(id)->add(elemGlobal.at(id)[n],index);
         }
         index++;
      }

      _cellParts.at(id)->finishAdd();

   }

    //deleting allocated arrays in this method
    vector< int* > ::iterator it_int;
    for ( it_int = elemGlobal.begin(); it_int != elemGlobal.end(); it_int++)
       delete [] *it_int;

    for ( it_int = distGlobal.begin(); it_int != distGlobal.end(); it_int++)
       delete [] *it_int;



}


//get boundary information for process
void
PartMesh::mapBounIDAndCell(int id)
{

   //mapBounIDAndCell  store global information
   //_mapBounIDAndCell store  local (process) information
    multimap<int,int>  mapBounIDAndCell;

    //boundary information has been stored
    const FaceGroupList&  boundaryFaceGroups = _meshList.at(id)->getBoundaryFaceGroups();
    int indx = _totElems.at(id);

    for ( int bounID = 0; bounID < int(boundaryFaceGroups.size()); bounID++){
        int group_id = boundaryFaceGroups.at(bounID)->id;
        _boundarySet.at(id).insert( group_id );
        string boun_type( boundaryFaceGroups.at(bounID)->groupType );
        _mapBounIDAndBounType.at(id).insert( pair<int,string>(group_id, boun_type) );

        int nBounElm = boundaryFaceGroups.at(bounID)->site.getCount();
        for ( int n = 0; n < nBounElm; n++){
           mapBounIDAndCell.insert( pair<int,int>(group_id, indx) );
           indx++;
        }
     }

    //loop over local elements of process
    multimap<int,int>::const_iterator  it;
    const CRConnectivity& cellCells = _meshList.at(id)->getCellCells();
    for ( int n = 0; n < _nelems.at(id); n++){
       int ncells = cellCells.getCount( _elem.at(id)[n] );
       for ( int m = 0; m < ncells; m++){
           int elem_id = cellCells( _elem.at(id)[n], m );
           //loop over boundary elements
           for ( it = mapBounIDAndCell.begin(); it != mapBounIDAndCell.end(); it++){
               if ( elem_id == it->second )
                   _mapBounIDAndCell.at(id).insert( pair<int,int>(it->first, it->second) );
           }
       } 
    }


    

}


//adding ghost cell elements to local elements
void
PartMesh::resize_elem(int id)
{
   int  tot_cells = _mapBounIDAndCell.at(id).size() + _nelems.at(id);
   _nelemsWithGhosts.at(id) = tot_cells;
   _elemWithGhosts.push_back(  new int[ tot_cells] );
   //assign old values
   for ( int n = 0;  n < _nelems.at(id); n++)
       _elemWithGhosts.at(id)[n] = _elem.at(id)[n];
   //ghost part assigned
    multimap<int,int>::const_iterator it;
    int indx = _nelems.at(id);
   for ( it = _mapBounIDAndCell.at(id).begin(); it != _mapBounIDAndCell.at(id).end(); it++){
       _elemWithGhosts.at(id)[indx] = it->second;
       indx++;
   }


}



//forming faceCells and faceNodes on each local mesh
void 
PartMesh::faceCells_faceNodes()
{
    vector< ArrayIntPtr  >   indices;
    for ( int id = 0; id < _nmesh; id++){
        //form site 
        int face_count = _partFaces.at(id)->getCount( _procID );
        int node_count = _partNodes.at(id)->getCount( _procID );

         _faceSite.push_back( StorageSitePtr(new  StorageSite(face_count)) );
        _nodeSite.push_back( StorageSitePtr(new  StorageSite(node_count)) );

        const Array<int>&  row = _partFaces.at(id)->getRow();
        const Array<int>&  col = _partFaces.at(id)->getCol();
        //forming indices
        indices.push_back( ArrayIntPtr(new Array<int>( face_count )) );
        int n_start = row[_procID];
        int indx = 0;
        for ( int n = n_start; n < n_start + face_count; n++){
            (*indices.at(id))[indx] = col[n];
            indx++;
        }
        //getting subset from global _faceCellsGlobal and _faceNodesGlobal
        int cell_count = _nelemsWithGhosts.at(id) + _interfaceMap.at(id).size();
       _cellSite.push_back( StorageSitePtr(new StorageSite(cell_count)) );
      _faceCells.push_back( _faceCellsGlobal.at(id)->getLocalizedSubset( *_faceSite.at(id), *_cellSite.at(id), *indices.at(id) )  );
      _faceNodes.push_back( _faceNodesGlobal.at(id)->getLocalizedSubset( *_faceSite.at(id), *_nodeSite.at(id), *indices.at(id) )  );
       _cellCells.push_back( (_faceCells.at(id)->getTranspose())->multiply(*_faceCells.at(id), true) );
       _cellNodes.push_back( (_faceCells.at(id)->getTranspose())->multiply(*_faceNodes.at(id), false) );

    }



}


//form interfaces
void
PartMesh::interfaces()
{
    _interfaceMap.resize( _nmesh );
    for ( int id = 0; id < _nmesh; id++){
       int nface =   _partFaces.at(id)->getCount( _procID );
        for ( int face = 0; face < nface; face++ ){
             int face_globalID = (*_partFaces.at(id))(_procID,face);

            if (_faceParts.at(id)->getCount(face_globalID) == 2 ){  // ==2 means sharing interface
                int neighPart = (*_faceParts.at(id))(face_globalID,0) +           
                                (*_faceParts.at(id))(face_globalID,1) - _procID; 
                _interfaceSet.at(id).insert( neighPart );
                _interfaceMap.at(id).insert(  pair<int,int>(neighPart,face) );
            }
        }
    }

}


//form coordinates 
void
PartMesh::coordinates()
{

  for ( int id = 0; id < _nmesh; id++){
      const Mesh& mesh = *(_meshList.at(id));
        const Array<Mesh::VecD3>&   global_coord = mesh.getNodeCoordinates();
        int node_count = _partNodes.at(id)->getCount( _procID );
        _coord.push_back( ArrayVecD3Ptr(new Array<Mesh::VecD3>(node_count)) );
        const Array<int>&    rowPartNodes = _partNodes.at(id)->getRow();
        const Array<int>&    colPartNodes = _partNodes.at(id)->getCol();

        for ( int node = 0; node < node_count; node++)
            (*_coord.at(id))[node] = global_coord[ colPartNodes[rowPartNodes[_procID]+node] ];
    }


   
}


int
PartMesh::count_interior_faces( int id )
{
   return _partFaces.at(id)->getCount(_procID) - (_nelemsWithGhosts.at(id) - _nelems.at(id)) 
          - _interfaceMap.at(id).size();

}

void
PartMesh::local_number_elems()
{
    for ( int id = 0; id < _nmesh; id++){
        const Array<int>& globalToLocal = _faceCells.at(id)->getGlobalToLocalMap();
         for ( int n = 0; n < _nelems.at(id); n++){
               int local_elem_id = globalToLocal[ _elem.at(id)[n] ];
               _elemLocal.at(id).insert(local_elem_id);
         }
    }

}
//store all boundary and interface cells into vector for using "find" in algorithm
void
PartMesh::non_interior_cells()
{
    multimap<int,int>::const_iterator it;
    for ( int id = 0; id < _nmesh; id++){
        //all boundary and interface  cells
        for ( it = _mapBounIDAndCell.at(id).begin(); it != _mapBounIDAndCell.at(id).end(); it++){
            int local_cell_id = _faceCells.at(id)->getGlobalToLocalMap()[it->second];
           _nonInteriorCells.at(id).insert( local_cell_id );
        }

        for ( it = _interfaceMap.at(id).begin(); it != _interfaceMap.at(id).end(); it++ ){
            int cell_0 = (*_faceCells.at(id))(it->second,0);
            int cell_1 = (*_faceCells.at(id))(it->second,1);
            if ( _elemLocal.at(id).count(cell_0) == 0 )
               _nonInteriorCells.at(id).insert(cell_0);

            if ( _elemLocal.at(id).count(cell_1) == 0 )
               _nonInteriorCells.at(id).insert(cell_1);

        }

    }

}

//faceCells and faceNodes are order such that interior face come first and 
//then  boundary faces or interfaces follows after that
//interface and boundary cells which are always stored as second element 
//faceCellsOrdered(face,0) => interior cells, faceCellsOrdered(face,1)=>boundary or interface cells
void 
PartMesh::order_faceCells_faceNodes()
{

     for ( int id = 0; id < _nmesh; id++ ){
        int tot_cells = _nelemsWithGhosts.at(id) + _interfaceMap.at(id).size();
        construct_mesh( id );

        _faceCellsOrdered.push_back( CRConnectivityPtr( new  CRConnectivity(_meshListLocal.at(id)->getFaces(), _meshListLocal.at(id)->getCells() ) ) );
        _faceNodesOrdered.push_back( CRConnectivityPtr( new  CRConnectivity(_meshListLocal.at(id)->getFaces(), _meshListLocal.at(id)->getNodes() ) ) );
        _cellToOrderedCell[id].assign(tot_cells, -1);
        _isOneToOneMap[id].assign(tot_cells, true); 

          //faceCells 
         _faceCellsOrdered.at(id)->initCount();
         _faceNodesOrdered.at(id)->initCount();

         int nface = _partFaces.at(id)->getCount(_procID);
         int count_node = _faceNodes.at(id)->getRow()[1] - _faceNodes.at(id)->getRow()[0];
         int count_cell = _faceCells.at(id)->getRow()[1] - _faceCells.at(id)->getRow()[0];
         for ( int face = 0; face < nface; face++){
            _faceCellsOrdered.at(id)->addCount(face,count_cell);  //two cells (always)
            _faceNodesOrdered.at(id)->addCount(face,count_node);  //two, three or four nodes
         }

         _faceCellsOrdered.at(id)->finishCount();
         _faceNodesOrdered.at(id)->finishCount();

         //start with interior faces
          int array_length = _faceCells.at(id)->getGlobalToLocalMap().getLength();
         _globalToLocalMap.push_back(  ArrayIntPtr( new Array<int>(array_length) ) );
         array_length = _faceCells.at(id)->getLocalToGlobalMap().getLength();
         _localToGlobalMap.push_back(  ArrayIntPtr( new Array<int>(array_length) ) );

         *_globalToLocalMap.at(id) = -1; 
         *_localToGlobalMap.at(id) = -1; 

         int cellID = 0;
         int face_track = 0;
         int nface_local = _partFaces.at(id)->getCount( _procID );
         for ( int face = 0; face < nface_local; face++){
             int cell_0 = (*_faceCells.at(id))(face,0);
             int cell_1 = (*_faceCells.at(id))(face,1);
             bool is_interior = _nonInteriorCells.at(id).count(cell_0) == 0 &&
                                _nonInteriorCells.at(id).count(cell_1) == 0; 
              if ( is_interior ) {
                  bool is_Counted =  ( _cellToOrderedCell[id][cell_0] != -1 );
                  if (  !is_Counted ){
                      _cellToOrderedCell[id][cell_0] = cellID;
                      _faceCellsOrdered.at(id)->add(face_track,cellID);
                       int global_id = _faceCells.at(id)->getLocalToGlobalMap()[cell_0];
                       (*_globalToLocalMap.at(id))[global_id] = cellID;
                       (*_localToGlobalMap.at(id))[cellID] = global_id;
                       _globalToLocalMappers.at(id).insert( pair<int,int>(global_id,cellID)  );
                       _localToGlobalMappers.at(id).insert( pair<int,int>(cellID, global_id) );
                       cellID++;
                  } else {
                      _faceCellsOrdered.at(id)->add(face_track,_cellToOrderedCell[id][cell_0]);
                  } 

                  is_Counted = (_cellToOrderedCell[id][cell_1] != -1);
                  if (  !is_Counted ){
                      _cellToOrderedCell[id][cell_1] = cellID;
                      _faceCellsOrdered.at(id)->add(face_track,cellID);
                       int global_id = _faceCells.at(id)->getLocalToGlobalMap()[cell_1];
                       (*_globalToLocalMap.at(id))[global_id] = cellID;
                       (*_localToGlobalMap.at(id))[cellID] = global_id;
                       _globalToLocalMappers.at(id).insert( pair<int,int>(global_id,cellID)  );
                       _localToGlobalMappers.at(id).insert( pair<int,int>(cellID, global_id) );
                      cellID++;
                  } else {
                      _faceCellsOrdered.at(id)->add(face_track,_cellToOrderedCell[id][cell_1]);
                  } 


                 for  ( int node = 0; node < count_node; node++ )
                    _faceNodesOrdered.at(id)->add( face_track, (*_faceNodes.at(id))( face, node ) );

                 face_track++;
              }
         }

        //then boundary faces
       multimap<int,int>::const_iterator it_cell;
       pair<multimap<int,int>::const_iterator,multimap<int,int>::const_iterator> it;
       set<int> ::const_iterator it_set;
       int offset = face_track;
       //loop over boundaries
       for ( it_set = _boundarySet.at(id).begin(); it_set != _boundarySet.at(id).end(); it_set++){
          int bndryID = *it_set;

          it = _mapBounIDAndCell.at(id).equal_range(bndryID);
          //if it is not empty
          if ( _mapBounIDAndCell.at(id).count( bndryID ) > 0 )
             _bndryOffsets.at(id).insert( pair<int,int>(bndryID, offset) );
 
          for ( it_cell = it.first; it_cell != it.second; it_cell++ ){

               int elem_0 =  _faceCells.at(id)->getGlobalToLocalMap()[it_cell->second];
               int elem_1 =  (*_cellCells.at(id))(elem_0, 0);
               int inner_elem = _cellToOrderedCell[id][elem_1];
               int outer_elem = cellID;

               //update globalToLocal and localToGlobalMaps
               (*_globalToLocalMap.at(id))[it_cell->second] = cellID;
               (*_localToGlobalMap.at(id))[cellID] = it_cell->second;
               _globalToLocalMappers.at(id).insert( pair<int,int>(it_cell->second,cellID ) );
               _localToGlobalMappers.at(id).insert( pair<int,int>(cellID, it_cell->second) );

               _cellToOrderedCell[id][elem_0] = cellID;

              _faceCellsOrdered.at(id)->add(face_track, inner_elem);
              _faceCellsOrdered.at(id)->add(face_track, outer_elem);

               int count_node = _faceNodes.at(id)->getRow()[1] - _faceNodes.at(id)->getRow()[0];
               for  ( int node = 0; node < count_node; node++)
                    _faceNodesOrdered.at(id)->add( face_track, (*_cellNodes.at(id))(elem_0,node) );

              face_track++;
              offset++;
              cellID++;
          }
       }

       //then interface faces
       multimap<int,int>::const_iterator it_face;
       for ( it_set = _interfaceSet.at(id).begin(); it_set != _interfaceSet.at(id).end(); it_set++){
          int interfaceID = *it_set;
          it = _interfaceMap.at(id).equal_range( interfaceID );
          _interfaceOffsets.at(id).insert( pair<int,int>(interfaceID,offset) ) ;
          for ( it_face = it.first; it_face != it.second; it_face++ ){
              int face_id = it_face->second;
              int elem_0 =  (*_faceCells.at(id))(face_id,0);
              int elem_1 =  (*_faceCells.at(id))(face_id,1);
              int outer_elem_id = -1;

             if ( _nonInteriorCells.at(id).count( elem_1 ) > 0 ){ //if elem_1 is non-interior cell
                _faceCellsOrdered.at(id)->add(face_track,_cellToOrderedCell[id][elem_0]);
                 outer_elem_id = elem_1;
             } else { 
                _faceCellsOrdered.at(id)->add(face_track,_cellToOrderedCell[id][elem_1]);
                outer_elem_id = elem_0;
             }
            _faceCellsOrdered.at(id)->add(face_track,cellID);

             //update maps
             int global_id = _faceCells.at(id)->getLocalToGlobalMap()[outer_elem_id];
             (*_globalToLocalMap.at(id))[global_id] = cellID;
             (*_localToGlobalMap.at(id))[cellID] = global_id;
             _globalToLocalMappers.at(id).insert( pair<int,int>(global_id, cellID) );
             _localToGlobalMappers.at(id).insert( pair<int,int>(cellID, global_id) );
             _cellToOrderedCell[id][outer_elem_id] = cellID;

             int count_node = _faceNodes.at(id)->getRow()[1] - _faceNodes.at(id)->getRow()[0];

             if ( outer_elem_id == elem_1 ) {
                 for  ( int  node = 0; node < count_node; node++)
                    _faceNodesOrdered.at(id)->add( face_track, (*_faceNodes.at(id))( face_id, node ) );
             } else {
                 for  ( int  node = count_node-1; node >= 0; node--)
                    _faceNodesOrdered.at(id)->add( face_track, (*_faceNodes.at(id))( face_id, node ) );
             }

              face_track++;
              offset++;
              cellID++;
          }
       }

        _faceCellsOrdered.at(id)->finishAdd(); 
        _faceNodesOrdered.at(id)->finishAdd();

//        //checking _cellToOrderedCell
//        double sum = 0;
//        for ( int n = 0; n < _cellToOrderedCell.at(id).size(); n++ )
//           sum += _cellToOrderedCell[id][n] - n;
// 
//        if ( _procID ==  0 )
//           for ( int n = 0; n < _cellToOrderedCell.at(id).size(); n++ )
//               cout << " celToOrededcell[" << n+1 << "] = " << _cellToOrderedCell[id][n]+1 << endl;
// 
//        cout << " proc id = " << _procID << "  sum ======= " << sum << endl;
//        assert ( !sum );  //sum should be zero 
//          if ( _procID == 0 ){ 
//            multimap<int,int>::iterator it_test;
//            pair<multimap<int,int>::iterator, multimap<int,int>::iterator> ret;
//            for ( int key = 0; key < 38; key++){
//                 cout << key+1 << " ==>";
//                 ret = _globalToLocalMappers.at(id).equal_range(key);
//                 for ( it_test = ret.first; it_test != ret.second; it_test++)
//                     cout << " " << it_test->second+1;
//                cout << endl;
//            }
//         }
// 

     }

}


void 
PartMesh::construct_mesh( int id )
{
        int dim = _meshList.at(id)->getDimension();
        _meshListLocal.push_back(  new Mesh( dim, _procID)   );

        StorageSite& faceSite = _meshListLocal.at(id)->getFaces();
        StorageSite& cellSite = _meshListLocal.at(id)->getCells();
        StorageSite& nodeSite = _meshListLocal.at(id)->getNodes();
        int nface_local = _partFaces.at(id)->getCount( _procID );
        int tot_cells = _nelemsWithGhosts.at(id) + _interfaceMap.at(id).size();
        int nGhostCell_local =  tot_cells - _nelems.at(id);
        int nnode_local =_partNodes.at(id)->getCount( _procID );
        //Storage sites
        faceSite.setCount( nface_local );
        cellSite.setCount( _nelems.at(id), nGhostCell_local );
        nodeSite.setCount( nnode_local );
}

//collect each mesh neightbourhood interface cells, ...
void 
PartMesh::exchange_interface_meshes()
{
    vector<int>  offset;
    vector<int>  interfaceMeshIDs;
    int *recv_counts = NULL;
    int *displ       = NULL;
    for ( int id = 0; id < _nmesh; id++){
        recv_counts = new int[ _nPart.at(id) ];
        displ       = new int[ _nPart.at(id) ];

        int total_interface_mesh = int( _interfaceSet.at(id).size() );
        int total_faces = int( _interfaceMap.at(id).size() );

        MPI::COMM_WORLD.Allgather(&total_interface_mesh, 1, MPI::INT, _interfaceMeshCounts.at(id)->getData(), 1, MPI::INT);
        MPI::COMM_WORLD.Allgather(&total_faces, 1, MPI::INT, _procTotalInterfaces.at(id)->getData(), 1, MPI::INT);

        //now find offsets for ghostCells 
        int total_interface_local = _interfaceSet.at(id).size();
        int total_interface_global = -1;

        MPI::COMM_WORLD.Allreduce( &total_interface_local, &total_interface_global, 1, MPI::INT, MPI::SUM );
        MPI::COMM_WORLD.Allgather( &total_interface_local, 1, MPI::INT, recv_counts, 1, MPI::INT );
        MPI::COMM_WORLD.Allreduce( &total_faces, &_windowSize.at(id), 1, MPI::INT, MPI::MAX);

      //enough space for gathering
       _offsetInterfaceCells.push_back(  ArrayIntPtr( new Array<int>(total_interface_global) )  );
       _interfaceMeshIDs.push_back    (  ArrayIntPtr( new Array<int>(total_interface_global) )  ); 

       _ghostCellsGlobal.push_back    (  ArrayIntPtr( new Array<int>(total_faces ) )  );
       _ghostCellsLocal.push_back     (  ArrayIntPtr( new Array<int>(total_faces ) )  );

        //local offset and interfaceMeshID are stored in contigous memory
        int index = 0;
        set<int>::const_iterator it_set;
        for ( it_set = _interfaceSet.at(id).begin(); it_set != _interfaceSet.at(id).end(); it_set++){
            int neighMeshID = *it_set;
            interfaceMeshIDs.push_back( neighMeshID );
            //loop over interface
             int  nstart = _interfaceOffsets.at(id)[neighMeshID];
             offset.push_back( nstart );

             int  nend  = nstart + _interfaceMap.at(id).count( neighMeshID );
            for ( int n = nstart;  n < nend; n++){
                  int elem_local_id  = (*_faceCellsOrdered.at(id))(n,0);
                  int elem_global_id = (*_localToGlobalMap.at(id))[ elem_local_id ];

                 (*_ghostCellsLocal.at(id))[index]  = elem_local_id;
                 (*_ghostCellsGlobal.at(id))[index] = elem_global_id;
                 index++;
            }

        }

        displ[0] = 0;
        for ( int i = 1; i < _nPart.at(id); i++)
            displ[i] = recv_counts[i-1] + displ[i-1]; 


         //now gather for _interface...
         MPI::COMM_WORLD.Allgatherv( &offset[0], total_interface_local, MPI::INT, 
                  _offsetInterfaceCells.at(id)->getData(), recv_counts, displ, MPI::INT); 

         MPI::COMM_WORLD.Allgatherv( &interfaceMeshIDs[0], total_interface_local, MPI::INT, 
                  _interfaceMeshIDs.at(id)->getData(), recv_counts, displ, MPI::INT); 

       offset.clear();
       interfaceMeshIDs.clear();
 
       delete [] recv_counts;
       delete [] displ;

    }
}


void
PartMesh::mappers()
{


    for ( int id = 0; id < _nmesh; id++){
 
        create_window( id );
        fence_window();      

        StorageSite::ScatterMap & cellScatterMap = _meshListLocal.at(id)->getCells().getScatterMap();
        StorageSite::GatherMap  & cellGatherMap  = _meshListLocal.at(id)->getCells().getGatherMap();

       //getting data
        set<int>::const_iterator it_set;
        int interfaceIndx = 0;
        for ( it_set = _interfaceSet.at(id).begin(); it_set != _interfaceSet.at(id).end(); it_set++){

            int neighMeshID = *it_set;
            int size = int(_interfaceMap.at(id).count(neighMeshID) );
            _fromIndices.at(id).push_back( ArrayIntPtr( new Array<int>(size) ) );
              _toIndices.at(id).push_back( ArrayIntPtr( new Array<int>(size) ) );
           *_fromIndices.at(id).at(interfaceIndx)  = -1;
              *_toIndices.at(id).at(interfaceIndx) = -1;

             int window_displ = -1;
             window_displ = get_window_displ( id, neighMeshID );
            _winLocal.Get ( _fromIndices.at(id).at(interfaceIndx)->getData(), size, MPI::INT, neighMeshID, window_displ, size, MPI::INT );
            _winGlobal.Get( _toIndices.at(id).at(interfaceIndx)->getData()  , size, MPI::INT, neighMeshID, window_displ, size, MPI::INT );
             interfaceIndx++;

        }

       fence_window();
       free_window();

       interfaceIndx = 0;
       for ( it_set = _interfaceSet.at(id).begin(); it_set != _interfaceSet.at(id).end(); it_set++ ){

           int neighMeshID = *it_set;
           int size = int(_interfaceMap.at(id).count(neighMeshID) );
           map<int, int>  mapKeyCount;        //map between key and count of that key
           for ( int n = 0; n < size; n++){

               int  key     = (*_toIndices.at(id).at(interfaceIndx))[n];
               int  count   = _globalToLocalMappers.at(id).count( key ); 

               if ( mapKeyCount.count( key ) > 0 ) { //it has elements
                  mapKeyCount[key] = mapKeyCount[key] + 1; //increase one
               } else {  //if it is empty  
                  mapKeyCount.insert(pair<int,int>(key,0));
               }

               multimap<int,int>::const_iterator it;
               it = _globalToLocalMappers.at(id).lower_bound( key );
               for ( int n_iter = 0; n_iter < mapKeyCount[key]; n_iter++)
                   it++;

               int elem_id =  it->second;
               (*_toIndices.at(id).at(interfaceIndx))[n] =  elem_id;

//                if ( _procID == 0 ) {
//                   cout << " neighMeshID = " << neighMeshID << endl;
//                   cout << " size        = " << size        << endl;
//                   cout << " key         = " << key         << endl;
//                   cout << " count       = " << count       << endl;
//                   cout << " elem_id     = " << elem_id     << endl;
//                   cout << " mapKeyCount[" << key << "] = " << mapKeyCount[key] << endl;
//                   cout << " _fromIndices= " << (*_fromIndices.at(id).at(interfaceIndx))[n] << endl;
//                   cout << " _toIncides  = " << (*_toIndices.at(id).at(interfaceIndx))[n]   << endl;
//                   cout << endl;
//                }


           }
          //from indices seems useless for now but we need to find scatterCells = cellCells(toIndices) and
          //use fromindices as storage Array
          for ( int i = 0; i < _fromIndices.at(id).at(interfaceIndx)->getLength(); i++){
               int elem_id = (*_toIndices.at(id).at(interfaceIndx))[i];
                (*_fromIndices.at(id).at(interfaceIndx))[i] = _meshListLocal.at(id)->getCellCells()(elem_id,0);
          }

          cellScatterMap[ _meshListLocal.at(id)->getGhostCellSite( neighMeshID ) ] = _fromIndices.at(id).at(interfaceIndx);
          cellGatherMap [ _meshListLocal.at(id)->getGhostCellSite( neighMeshID ) ] = _toIndices.at(id).at(interfaceIndx);
 
          interfaceIndx++;

        }

    }

}


int
PartMesh::get_window_displ( int id, int neigh_mesh_id )
{
      int loc = 0;
      int window_displ = 0;
      for ( int i = 0; i < neigh_mesh_id; i++)
          loc += (*_interfaceMeshCounts.at(id))[i];
         
      while ( (*_interfaceMeshIDs.at(id))[loc] != _procID){
          window_displ += (*_offsetInterfaceCells.at(id))[loc+1] - (*_offsetInterfaceCells.at(id))[loc];
         loc++;
      }


    return window_displ;
}

void
PartMesh::create_window( int id )
{
           int int_size = MPI::INT.Get_size();
           MPI::Aint lb, sizeofAint;
           MPI::INT.Get_extent(lb,sizeofAint);

           int window_size = _windowSize.at(id);  //already MPI::MAX  has taken  maximum size
           _winLocal  = MPI::Win::Create(_ghostCellsLocal.at(id)->getData(), window_size*sizeofAint, int_size,
                                         MPI_INFO_NULL, MPI::COMM_WORLD);
           _winGlobal = MPI::Win::Create(_ghostCellsGlobal.at(id)->getData(), window_size*sizeofAint, int_size,
                                         MPI_INFO_NULL, MPI::COMM_WORLD);

}

void
PartMesh::free_window( )
{
     _winLocal.Free();
     _winGlobal.Free();
}

void 
PartMesh::fence_window()
{

//     _winLocal.Fence ( 0 );
//     _winGlobal.Fence( 0 );

      _winLocal.Fence(MPI::MODE_NOPUT);
      _winGlobal.Fence(MPI::MODE_NOPUT);

}

//dump all mesh information 
void
PartMesh::mesh_debug()
{
    mesh_file();
    mesh_tecplot();


}

void 
PartMesh::mesh_file()
{
     stringstream ss;
     ss << "mesh_proc" << _procID << "_info.dat";
     ofstream  mesh_file( (ss.str()).c_str() );
     for ( int id = 0; id < _nmesh; id++ ){

          const StorageSite::ScatterMap& cellScatterMap  = _meshListLocal.at(id)->getCells().getScatterMap();
          const StorageSite::GatherMap& cellGatherMap   = _meshListLocal.at(id)->getCells().getGatherMap();
          const Mesh::GhostCellSiteMap& ghostCellSiteMap = _meshListLocal.at(id)->getGhostCellSiteMap();

           Mesh       ::GhostCellSiteMap::const_iterator it_ghost;
           StorageSite::ScatterMap::const_iterator it_mapper;
           it_ghost = ghostCellSiteMap.begin();
            //loop over interfaces
           for ( it_mapper = cellScatterMap.begin(); it_mapper != cellScatterMap.end(); it_mapper++){
               const StorageSite *site = it_mapper->first;
               const Array<int>&  scatterArray = *(it_mapper->second);
               const Array<int>&  gatherArray  = *(cellGatherMap.find(site)->second);
               for ( int i = 0; i < site->getCount(); i++){
                     mesh_file <<   "  neightMeshID = " <<  it_ghost->first  << "        "
                               << gatherArray[i]  + 1  << "    ===>    " 
                               << scatterArray[i] + 1  << endl;
               }
               it_ghost++;
          }

      }

      mesh_file.close();

}




//need modification for quad, hexa, tetra
void
PartMesh::mesh_tecplot()
{
     stringstream ss;
     ss << "mesh_proc" << _procID << ".dat";
     ofstream  mesh_file( (ss.str()).c_str() );

      const Mesh& mesh = *(_meshListLocal.at(0));
      const CRConnectivity&  cellNodes = mesh.getCellNodes();

     const Array<Mesh::VecD3>&  coord = mesh.getNodeCoordinates();
     int tot_elems = cellNodes.getRowDim();
     int tot_nodes = cellNodes.getColDim();

     mesh_file << "title = \" tecplot file for process Mesh \" "  << endl;
     mesh_file << "variables = \"x\",  \"y\", \"z\", \"cell_type\" " << endl;
     mesh_file << "zone N = " << tot_nodes << " E = " << tot_elems <<
    " DATAPACKING = BLOCK,  VARLOCATION = ([4]=CELLCENTERED), ZONETYPE=FETRIANGLE " << endl;

     //x 
     for ( int n = 0; n < tot_nodes; n++){
         mesh_file << scientific  << coord[n][0] << "     " ;
         if ( n % 5 == 0 ) mesh_file << endl;
     }

     mesh_file << endl;
     //y
     for ( int n= 0; n < tot_nodes; n++){
         mesh_file << scientific << coord[n][1] << "     ";
         if ( n % 5 == 0 ) mesh_file << endl;
     }


     mesh_file << endl;
     //z
     for ( int n = 0; n < tot_nodes; n++){
          mesh_file << scientific << coord[n][2] << "     ";
          if ( n % 5 == 0) mesh_file << endl;
     }

     mesh_file <<  endl;
     mesh_file << endl;  
     //cell type
      int cell_type = -1;
      for ( int n = 0; n < tot_elems;  n++){
          int elem_id = _cellToOrderedCell[0][n];
          if ( _nonInteriorCells.at(0).count(elem_id) == 0 ){
             cell_type = 0;
          } else {
             cell_type = 1;
          }

          mesh_file << cell_type  << "      ";
          if (  n % 10 == 0 ) mesh_file << endl;

       }

      mesh_file << endl;
      mesh_file << endl;
     //connectivity
    for (int n = 0; n < tot_elems; n++){
         int nnodes =  cellNodes.getCount(n);
         if (  n < _nelems.at(0) ){
            for ( int node = 0; node < nnodes; node++)
                mesh_file << cellNodes(n,node)+1 << "      ";
         } else {
                mesh_file << cellNodes(n,0)+1 << "      " << cellNodes(n,1)+1 <<
                "       " << cellNodes(n,0)+1 << "      ";
         }
        mesh_file << endl;
     }


    mesh_file.close();

}
