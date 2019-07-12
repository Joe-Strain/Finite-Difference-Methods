MODULE fin_diff
! *** Module containing subroutines and functions required to run finite
! *** difference methods to approximate PDEs

  IMPLICIT NONE 
  CONTAINS

! *****************************************************************************

  SUBROUTINE timestep(uu,dt,alpha,dy,ny)
! *** Subroutine to calc the solution to PDE as we step along in time
! *** as we step along in time
    REAL, DIMENSION(:), INTENT(INOUT) :: uu
    REAL, INTENT(IN) :: dt, dy, alpha
    INTEGER, INTENT(IN) :: ny

    REAL, DIMENSION(SIZE(uu)) :: du

    CALL pde(alpha,uu,du,dy,ny)

! *** Timestep the internal points of the spatial grid
    uu(2:ny−1) = uu(2:ny−1) + dt*du(2:ny−1)
	
  END SUBROUTINE timestep

! *****************************************************************************

  SUBROUTINE pde(alpha,uu,du,dy,ny)
! *** Subroutine to calculate the right−hand side of out original PDE

    REAL, DIMENSION(:), INTENT(IN) :: uu 
    REAL, DIMENSION(:), INTENT(OUT) :: du
    REAL, INTENT(IN) :: alpha, dy
    INTEGER, INTENT(IN) :: ny

    du = alpha*der2(uu,dy,ny)

  END SUBROUTINE pde

! *****************************************************************************

  FUNCTION der2(uu,dy,ny)
! ***  This calculates the  2nd order central difference of the second
! *** derivative of our function
    REAL, DIMENSION(:), INTENT(IN) :: uu
    REAL, INTENT(IN) :: dy
    INTEGER, INTENT(IN) :: ny

    REAL, DIMENSION(SIZE(uu)) :: der2

    der2(2:ny−1) = (uu(3:ny) − 2*uu(2:ny−1) + uu(1:ny−2))/dy**2
	
  END FUNCTION der2

! *****************************************************************************

  SUBROUTINE dftimestep(uu,dt,alpha,dy,ny,tempuu)
! *** Updates the uu values by one time step using the DuFort−Frankel method

    REAL, DIMENSION(:) , INTENT(INOUT) :: uu
    REAL, DIMENSION(:), INTENT(IN) :: tempuu
    REAL, INTENT(IN) :: dt, alpha, dy
    INTEGER, INTENT(IN) :: ny

    REAL :: cc

    cc = (2*alpha*dt)/dy**2
    uu(2:ny−1) = ((1−cc)*tempuu(2:ny−1) + cc*(uu(3:ny) + uu(1:ny−2))) / (1+cc)

  END SUBROUTINE dftimestep

!******************************************************************************

  SUBROUTINE laastimestep(uu,ny,inuu,k)
! *** Updates the uu values by one timestep using the Laasonen implicit method

    REAL, DIMENSION(:), INTENT(INOUT) :: uu
    REAL, INTENT(IN) :: inuu, k
    INTEGER, INTENT(IN) :: ny

    REAL, DIMENSION(SIZE(uu)) :: dd
    REAL :: a, b, c

! *** Set up our values for the Thompson Algorithm
    dd(1:ny) = −uu(1:ny)
    b = −(2*k+1)
    a = k 
    c = k

    CALL thompson(uu,a,b,c,dd,inuu,ny)

  END SUBROUTINE laastimestep

!******************************************************************************

  SUBROUTINE cranktimestep(uu,ny,inuu,k)
! *** Updates the uu values by one timestep using the Crank−Nicolsen implicit
! *** method

    REAL, DIMENSION(:), INTENT(INOUT) :: uu
    REAL, INTENT(IN) :: inuu, k
    INTEGER, INTENT(IN) :: ny

    REAL, DIMENSION(SIZE(uu)) :: dd
    REAL :: a, b, c

! *** Set up our values for the Thompson Algorithm
    dd(2:ny−1) = (k/2)*(uu(3:ny) − 2*uu(2:ny−1) + uu(1:ny−2)) + uu(2:ny−1)
    dd(ny) = 0

    a = −k/2 
    b = 1+k
    c = −k/2
    CALL thompson (uu,a,b,c,dd,inuu,ny)

  END SUBROUTINE cranktimestep

! *****************************************************************************

  SUBROUTINE thompson(uu,a,b,c,dd,inuu,ny)
! *** A subroutine to perform Thompson Algorithm for tridiagonal matrices

    REAL, DIMENSION(:), INTENT(INOUT) :: uu 
    REAL, DIMENSION(:), INTENT(IN) :: dd
    REAL, INTENT(IN) ::  a, b, c, inuu 
    INTEGER, INTENT(IN) :: ny

    INTEGER :: i
    REAL, DIMENSION(ny) :: hh, gg

    hh(1) = 0
    gg(1) = inuu

    hh(2:ny) = c/(b − a*hh(1:ny−1))
    gg(2:ny) = (dd(2:ny) − a*gg(1:ny−1)) / (b − a*hh(1:ny−1))
    DO i=1, ny−2
      uu(ny−i) = −hh(ny−i) * uu(ny−i+1) + gg(ny−i)
    END DO

  END SUBROUTINE thompson

!***********************************************************

END MODULE fin_diff
