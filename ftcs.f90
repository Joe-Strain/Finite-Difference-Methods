PROGRAM ftcs

! *** Program: First time central space finite difference method to
! *** approximate  a PDE,  in  this   case  the  heat   diffusion equation

  USE fin_diff

  IMPLICIT NONE

  REAL, PARAMETER :: inuu=40, dt=0.002, alpha=0.000217, ymin=0, ymax=0.04
  REAL :: dy, tt, yr
  INTEGER, PARAMETER :: ny=41, nt=540, save_iter=1
  INTEGER :: n
  REAL, DIMENSION(ny) :: yy, uu
  CHARACTER(LEN=*), PARAMETER :: fname=”data/ftcs.dat”

! *** Set up values needed
  tt = 0 
  yr = ymax − ymin

! *** Create yy and give it values
  DO n=1, ny
    yy(n) = yr * (n−1)
  END DO
  
  yy = yy / (ny−1) + ymin
  dy = yy(2) − yy(1)

! *** Intial Conditions 
  uu(1) = inuu
  uu(2:ny) = 0

! *** Perform FTCS and save to ’fname’
  OPEN(UNIT=1,FILE=fname,FORM=”UNFORMATTED”,STATUS=”REPLACE”,ACTION=”WRITE”)

  WRITE(UNIT=1) uu, tt

  DO n=1, nt
    CALL timestep(uu,dt,alpha,dy,ny) 
    tt = tt + dt
    IF (MOD(n,save_iter)==0) THEN
      WRITE(UNIT=1) uu, tt
    END IF
  END DO 
  
  CLOSE(UNIT=1)

END PROGRAM ftcs
