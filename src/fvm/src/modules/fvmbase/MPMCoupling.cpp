// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#ifdef FVM_PARALLEL
#include <mpi.h>
#include <iostream>
#include <cassert>
#include "MPMCoupling.h"
#include "CRConnectivity.h"


using namespace std;

MPMCoupling::MPMCoupling(MPM& mpm, Field& coordinate, Field& velocity, StorageSite& cell_site  )
:_mpm(mpm), _coordinate(coordinate), _velocity(velocity), _cellSite(cell_site), MPM_PORTNAME_TAG(7777),
 _dtTAG(1122), _timeTAG(2122), _totParticlesTAG(3122), _isContinueTAG(4122), _particlesPosTAG(5122), _particlesVelTAG(6122)
{
    _procID = MPI::COMM_WORLD.Get_rank();
    PARENT = MPI::Comm::Get_parent();
    char port_name[ MPI::MAX_PORT_NAME];
    assert( PARENT != MPI::COMM_NULL );

    // get portName from Parent ( which is sent to parent by MPM) 
    if ( _procID == 0 )
         PARENT.Recv( port_name, MPI::MAX_PORT_NAME, MPI::CHAR, 0, MPM_PORTNAME_TAG, _portNameStatus);

    FVM_COMM_MPM = MPI::COMM_WORLD.Connect( port_name, MPI::INFO_NULL, 0);

    cout << " max-pot -name = " << MPI::MAX_PORT_NAME << " port_name = " << port_name << endl;
   _pType =  shared_ptr< Array<int> > ( new Array<int>(1) );

}


MPMCoupling::~MPMCoupling()
{    
    FVM_COMM_MPM.Barrier();
    FVM_COMM_MPM.Disconnect();
    PARENT.Barrier();
    PARENT.Disconnect();
   
}




void 
MPMCoupling::updateMPM()
{
    
}


void 
MPMCoupling::acceptMPM( )
{
    cout << " isContine  = " << _isContinueMPM << endl;
       FVM_COMM_MPM.Recv( &_isContinueMPM, 1, MPI::BOOL, 0, _isContinueTAG );
        cout << " isContine  = " << _isContinueMPM << endl;
       if ( _isContinueMPM ) {
           FVM_COMM_MPM.Recv( &_dtMPM, 1, MPI::DOUBLE, 0, _dtTAG   );
           FVM_COMM_MPM.Recv( &_timeMPM , 1, MPI::DOUBLE, 0, _timeTAG );
           FVM_COMM_MPM.Recv( &_totParticlesMPM, 1, MPI::INT, 0, _totParticlesTAG );
           _px = _coordinate[_cellSite].newSizedClone( _totParticlesMPM );
           _pv = _coordinate[_cellSite].newSizedClone( _totParticlesMPM );
           _pType->resize( _totParticlesMPM );
           *_pType = 1;
           FVM_COMM_MPM.Recv( _px->getData(), _px->getDataSize(), MPI::BYTE, 0, _particlesPosTAG );
           FVM_COMM_MPM.Recv( _pv->getData(), _pv->getDataSize(), MPI::BYTE, 0, _particlesVelTAG );
           const StorageSite& particleSite = _mpm.getParticles( _totParticlesMPM );
           _coordinate[particleSite] = *_px;
           _velocity  [particleSite] = *_pv;
           _mpm.setCoordinates( _px );
           _mpm.setVelocities ( _pv );
           _mpm.setTypes      ( _pType );

       }



}
#endif
