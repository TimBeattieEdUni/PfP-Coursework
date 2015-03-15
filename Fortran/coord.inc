C $Id: coord.inc,v 1.1 2002/01/08 12:30:07 spb Exp spb $
C This file defines a common block that contains the primary coordinates
C of the particles,
C
C  Nbody	Number of particles
C  Npair	Number of particle pairs
C  pos		Position of the particles
C  r            Distance of particle from central mass
C  vel		velocity of the particles
C  f		Forces acting on each particle
C  visc		viscosity coefficient for each particle
C  mass		mass of each particle
C  delta_pos	seperation vector for each particle pair
C  delta_r		seperation for each particle pair
C
      INTEGER Nbody, Npair
      PARAMETER( Nbody=4*1024, Npair=(Nbody*(Nbody-1))/2 )
      INTEGER Xcoord, Ycoord, Zcoord, Ndim
      PARAMETER( Xcoord=1, Ycoord=2, Zcoord=3, Ndim=3 )
      DOUBLE PRECISION r(Nbody)
      DOUBLE PRECISION pos(Nbody,Ndim), vel(Nbody,Ndim)
      DOUBLE PRECISION f(Nbody,Ndim), visc(Nbody), mass(Nbody)
      DOUBLE PRECISION delta_pos(Nbody*Nbody,Ndim)
      DOUBLE PRECISION delta_r(Nbody*Nbody)
      DOUBLE PRECISION wind(Ndim)
      INTEGER collisions

      COMMON /coord/ f, pos, visc, mass, vel, delta_pos, delta_r,
     $         r,collisions,wind

C
C Nominal diameter of the particles.
C
      DOUBLE PRECISION Size
      PARAMETER( Size = 1d0)

C
C Size of central mass.
C
      DOUBLE PRECISION M_central
      PARAMETER( M_central = 1000d0)

C
C Interaction strength for pairwise interactions.
C
      DOUBLE PRECISION G
      PARAMETER( G = 2d0 )


