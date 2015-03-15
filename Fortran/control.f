C
C $Id: control.f,v 1.1 2002/01/08 12:30:07 spb Exp spb $
C
C Control program for the MD update
C

      PROGRAM MD
      IMPLICIT NONE
      INCLUDE 'coord.inc'
      INTEGER i,j
      REAL start, stop
      REAL tarray(2)
      REAL etime

C  timestep value
      DOUBLE PRECISION dt
      PARAMETER( dt = 0.02 )

C  number of timesteps to use.
      INTEGER Nstep
      INTEGER Nsave
      PARAMETER(Nstep=100,Nsave=5)
      CHARACTER*80 name

C read the initial data from a file

      wind(XCoord) = 0.9
      wind(YCoord) = 0.4
      wind(ZCoord) = 0.0
      collisions=0
      OPEN(unit=1,file='input.dat',status='OLD')
      DO i=1,Nbody
        READ(1,10) mass(i),visc(i),
     $    pos(i,Xcoord),pos(i,Ycoord),pos(i,Zcoord),
     $    vel(i,Xcoord),vel(i,Ycoord),vel(i,Zcoord)
      END DO
      CLOSE(unit=1)

C
C Run 20 timesteps and time how long it took
C


      start = etime(tarray)
      DO j=1,Nsave
        CALL evolve(Nstep,dt)
        stop = etime(tarray)

        WRITE(*,*) (j*Nstep), ' timesteps took ',(stop-start)
        write(*,*) collisions, ' collisions'


C write result to a file
        WRITE(name,12) 'output.dat', j*Nstep
        OPEN(unit=1,file=name)
        DO i=1,Nbody
          WRITE(1,10) mass(i),visc(i),
     $      pos(i,Xcoord),pos(i,Ycoord),pos(i,Zcoord),
     $      vel(i,Xcoord),vel(i,Ycoord),vel(i,Zcoord)
        END DO
        CLOSE(unit=1)
      END DO

10    FORMAT(8E16.8)
12    FORMAT(A10,I3.3)

      END

