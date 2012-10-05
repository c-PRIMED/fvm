// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#ifndef _SHOCKTUBE_H_
#define _SHOCKTUBE_H_

#include "Array.h"
#include "Vector.h"

template<class T>
class ShockTube
{
public:
   
  typedef Vector<T,3> Vec3;

  ShockTube(const int nCells) :
    pL(10.0),
    pR(1.0),
    rhoL(8.0),
    rhoR(1.0),
    uL(0),
    uR(0),
    _nCells(nCells),
    _gmm(0.4),
    _nStages(5),
    _stageCoeffs(5),
    q(_nCells),
    pressure(_nCells),
    rho(nCells),
    u(nCells),
    M(nCells)
  {
    _stageCoeffs[0] = 0.25;
    _stageCoeffs[1] = 1.0/6.0;
    _stageCoeffs[2] = 0.375;
    _stageCoeffs[3] = 0.5;
    _stageCoeffs[4] = 1;
  }

  T getPressure(const Vec3& q)
  {
    const T u = q[1]/q[0];
    return (q[2] - 0.5*(q[0]*u*u))*_gmm;
  }

  Vec3 getConservedVars(const T rho, const T p, const T u)
  {
    Vec3 q;
    q[0] = rho;
    q[1] = rho*u;
    q[2] = p/_gmm + 0.5*rho*u*u;
    return q;
  }
   
  Vec3 getFlux(const Vec3& q)
  {
    const T u = q[1]/q[0];
    const T p = getPressure(q);
    Vec3 flux;
    flux[0] = q[1];
    flux[1] = q[0]*u*u + p;
    flux[2] = (q[2]+p)*u;
    return flux;
  }

  Vec3 getRoeFlux_Frink(const Vec3& q0, const Vec3& q1)
  {
    const T half(0.5);
    
    const T sr0 = sqrt(q0[0]);
    const T sr1 = sqrt(q1[0]);

    const T sr0psr1 = sr0 + sr1;
    const T u0 = q0[1]/q0[0];
    const T u1 = q1[1]/q1[0];

    const T p0 = getPressure(q0);
    const T p1 = getPressure(q1);

    const T h0 = (q0[2]+p0) / q0[0];
    const T h1 = (q1[2]+p1) / q1[0];
    
    const T rho_f = sr0 * sr1;
    const T u_f = (sr0*u0 + sr1*u1) / sr0psr1;
    const T h_f = (sr0*h0 + sr1*h1) / sr0psr1;

    const T c_f = sqrt( _gmm * (h_f - half*u_f*u_f));

    const T abs_uf = fabs(u_f);
    const T abs_upc = fabs(u_f + c_f);
    const T abs_umc = fabs(u_f - c_f);

    
    const Vec3 flux0 = getFlux(q0);
    const Vec3 flux1 = getFlux(q1);

    const T dr = q1[0]-q0[0];
    const T rdu = rho_f*(u1-u0);

    const T dp = p1-p0;

    const T dpc2 = dp*c_f*c_f;

    const T a0 = half*(dpc2 + rdu*c_f)*abs_upc;
    const T a4 = half*(dpc2 - rdu*c_f)*abs_umc;

    const T drmdpc2 = dr-dpc2;
    Vec3 dFlux;

    dFlux[0] = abs_uf*drmdpc2 + a0 + a4;
    dFlux[1] = abs_uf*drmdpc2*u_f + a0*(u_f+c_f) + a4*(u_f-c_f);
    dFlux[2] = abs_uf*(half*drmdpc2*u_f*u_f) + a0*(h_f+u_f*c_f) + a4*(h_f-u_f*c_f);
    
    Vec3 flux(flux0+flux1-dFlux);
    return half*flux;
  }
  
  Vec3 getRoeFlux(const Vec3& q0, const Vec3& q1)
  {
    const T half(0.5);
    
    const T sr0 = sqrt(q0[0]);
    const T sr1 = sqrt(q1[0]);

    const T sr0psr1 = sr0 + sr1;
    const T u0 = q0[1]/q0[0];
    const T u1 = q1[1]/q1[0];

    const T p0 = getPressure(q0);
    const T p1 = getPressure(q1);

    const T h0 = (q0[2]+p0) / q0[0];
    const T h1 = (q1[2]+p1) / q1[0];
    
    const T rho_f = sr0 * sr1;
    const T u_f = (sr0*u0 + sr1*u1) / sr0psr1;
    const T h_f = (sr0*h0 + sr1*h1) / sr0psr1;

    const T c_f = sqrt( _gmm * (h_f - half*u_f*u_f));

    const T abs_uf = fabs(u_f);
    const T abs_upc = fabs(u_f + c_f);
    const T abs_umc = fabs(u_f - c_f);

    
    const Vec3 flux0 = getFlux(q0);
    const Vec3 flux1 = getFlux(q1);

    const T dr = q1[0]-q0[0];
    const T rdu = rho_f*(u1-u0);

    const T dp = p1-p0;

    const T dpbyc2 = dp/(c_f*c_f);

    const T a0 = half*(dpbyc2 + rdu/c_f)*abs_upc;
    const T a4 = half*(dpbyc2 - rdu/c_f)*abs_umc;

    const T drmdpbyc2 = dr-dpbyc2;
    Vec3 dFlux;

    dFlux[0] = abs_uf*drmdpbyc2 + a0 + a4;
    dFlux[1] = abs_uf*drmdpbyc2*u_f + a0*(u_f+c_f) + a4*(u_f-c_f);
    dFlux[2] = abs_uf*(half*drmdpbyc2*u_f*u_f) + a0*(h_f+u_f*c_f)
      + a4*(h_f-u_f*c_f);
    
    Vec3 flux(flux0+flux1-dFlux);
    return half*flux;
  }
  
  
  void init()
  {
    Vec3 qL(getConservedVars(rhoL, pL, uL));
    Vec3 qR(getConservedVars(rhoR, pR, uR));

    for(int i=0; i<_nCells/2; i++)
      q[i] = qL;
    
    for(int i=_nCells/2; i<_nCells; i++)
      q[i] = qR;
  }
  
  void solve(const double dt, const int nsteps)
  {
    const T dx = T(1.0)/_nCells;

    Array<Vec3> R(_nCells);

    Array<Vec3> qN(_nCells);

    for(int n=0; n<nsteps; n++)
    {
        cout << "step " << n << endl;
        qN = q;
        for(int s=0; s<_nStages; s++)
        {
            R.zero();

            for(int f=1; f<_nCells; f++)
            {
                int c0 = f-1;
                int c1 = f;
                const Vec3 flux = getRoeFlux(q[c0], q[c1]);
                R[c0] += flux;
                R[c1] -= flux;
            }

            R[0] -= getRoeFlux(q[0],q[0]);
            R[_nCells-1] += getRoeFlux(q[_nCells-1],q[_nCells-1]);
            
            const T sCoeff = _stageCoeffs[s]*dt/dx;
            for(int c=0; c<_nCells; c++)
            {
                q[c] = qN[c] - sCoeff*R[c];
            }
        }
    }

    for(int c=0; c<_nCells; c++)
    {
        rho[c] = q[c][0];
        pressure[c] = getPressure(q[c]);
        u[c] = q[c][1]/q[c][0];
        T soundSpeed = sqrt((_gmm+1)*pressure[c]/rho[c]);
        M[c] = u[c]/soundSpeed;
    }
  }
  
  ArrayBase& getSolution() {return q;}
  ArrayBase& getPressure() {return pressure;}
  ArrayBase& getDensity() {return rho;}
  ArrayBase& getVelocity() {return u;}
  ArrayBase& getMachNumber() {return M;}
  

  T pL;
  T pR;
  T rhoL;
  T rhoR;
  T uL;
  T uR;
  
  const int _nCells;
  T _gmm;
  const int _nStages;
  Array<double> _stageCoeffs;
  Array<Vec3> q;
  Array<T> pressure;
  Array<T> rho;
  Array<T> u;
  Array<T> M;
};
#endif
