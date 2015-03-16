
      SUBROUTINE visc_force(N,f,visc,vel)
      IMPLICIT NONE
      INTEGER N, i
      DOUBLE PRECISION f(N),visc(N),vel(N)   
          DO i=1,N
            f(i) = -visc(i) * vel(i)
          END DO
      END

      SUBROUTINE wind_force(N,f,visc,vel)
      IMPLICIT NONE
      INTEGER N, i
      DOUBLE PRECISION f(N),visc(N),vel
          DO i=1,N
            f(i) = f(i) -visc(i) * vel
          END DO
      END

      SUBROUTINE add_norm(N,r,delta)
      IMPLICIT NONE
      INTEGER N, k
      DOUBLE PRECISION r(N),delta(N)
        DO k=1,N
          r(k) = r(k) + (delta(k) * delta(k))
        END DO
      END

      FUNCTION force(W,delta,r)
      IMPLICIT NONE
      DOUBLE PRECISION force,W,delta,r
        force=W*delta/(r**3)
        RETURN
      END

