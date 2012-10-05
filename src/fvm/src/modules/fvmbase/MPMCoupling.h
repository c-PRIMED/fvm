// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#ifndef MPMCoupling_H
#define MPMCoupling_H
#ifdef FVM_PARALLEL

#include <mpi.h>
#include "Mesh.h"
#include "Array.h"
#include <vector>
#include <set>

#include "MPM_Particles.h"
#include "Field.h"

class MPMCoupling{

    public:

       MPMCoupling( MPM& mpm, Field& coordinate, Field& velocity, StorageSite& cell_site );
      ~MPMCoupling();
       
       void updateMPM();
       void acceptMPM();



    private:

       MPM& _mpm;
       Field& _coordinate;
       Field& _velocity;
       StorageSite& _cellSite;

       int _procID;  
       MPI::Intercomm  PARENT;
       MPI::Intercomm  FVM_COMM_MPM;
       MPI::Status  _portNameStatus;
       MPI::Info    _connectInfo;
       int MPM_PORTNAME_TAG;

       bool     _isContinueMPM;
       double  _dtMPM;
       double  _timeMPM;
       int  _totParticlesMPM;

       int  _dtTAG;
       int  _timeTAG;
       int  _totParticlesTAG;
       int  _isContinueTAG;
       int  _particlesPosTAG;
       int  _particlesVelTAG;
       shared_ptr< ArrayBase > _px;
       shared_ptr< ArrayBase > _pv;
       shared_ptr< Array<int> > _pType;
        

};

#endif
#endif
