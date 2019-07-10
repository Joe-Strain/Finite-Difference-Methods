PROGRAM dufort

! *** Program: DuFort−Frankel method to
! *** approximate a PDE, in this case, the heat diffusion equation

  USE fin_diff

  IMPLICIT NONE

  REAL, PARAMETER :: inuu=40, dt=0.002, alpha=0.000217, ymin=0, ymax=0.04
  REAL :: dy, tt, yr
  INTEGER, PARAMETER :: ny=41, nt=540, save_iter=1
  INTEGER :: n
  REAL, DIMENSION(ny) :: yy, uu, tempuu, temp2uu
  CHARACTER(LEN=*), PARAMETER :: fname=”data/dufort.dat”

! *** Set up values needed
  tt = 0 
  yr = ymax − ymin

! *** Create array yy and give it appropriate values
  DO n=1, ny
    yy(n) = yr * (n − 1)
  END DO
  yy = yy / (ny − 1) + ymin 
  dy = yy(2) − yy(1)

! *** Initial conditions 
  uu(1) = inuu 
  uu(2:ny) = 0

! *** Perform one step of the FTCS method and save results to ’fname’
  OPEN(UNIT=1,FILE=fname,FORM=”UNFORMATTED”,STATUS=”REPLACE”,ACTION=”WRITE”)

  WRITE(UNIT=1) uu, tt

! *** Temporary array to hold values of uu one time−step behin
  tempuu = uu

  CALL timestep(uu,dt,alpha,dy,ny)
  tt = tt + dt
  IF (MOD(1,saveiter)==0) THEN
    WRITE(UNIT=1) uu, tt
  END IF

! *** Perform the DuFort−Frankel method and save results to ’fname’
! *** Only (nt−1) steps, as one step of FTCS already completed
  DO  n=1, nt−1	
    temp2uu = uu
    CALL dftimestep(uu,dt,alpha,dy,ny,tempuu)

    tempuu = temp2uu
    tt = tt + dt

    IF (MOD(n+1,saveiter)==0) THEN
      WRITE(UNIT=1) uu, tt
    END IF 
  END DO

  CLOSE(UNIT=1)

END PROGRAM dufort
