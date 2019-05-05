C   compile and link with SLATEC and LAPACK libraries
C   -L /opt/local/lib -lslatec -llapack

      SUBROUTINE BLADESETUP(PARAM)
C   setup the grid and initial blade shape
C   infinite number of thin blades, e.g. tangentially averaged flow
C   GRID: NGRID, XGRID(100); additional IGRID, BLADESEGMENT
C   BLADE:  BLADEOLD/NEW(2*(I-1)+1) = blade shape
C           BLADEOLD/NEW(2*(I-1)+2) = blade slope
C           I=1,...,NGRID
C           BLADEOLD/NEW(200)
C   PARAMS: PAR(3) loading function parameters, optimisation space
C   FLOWDIR: upstream (1) and downstream (2) flow direction VY/VX
C           VX=1 (reference velocity), thus we have VY1 and VY2
      IMPLICIT NONE
      DOUBLEPRECISION PARAM
      DIMENSION PARAM(*)
      INTEGER NMAX
      PARAMETER (NMAX=100)
      DOUBLEPRECISION XGRID(NMAX),
     &  BLADESHAPEOLD(NMAX),BLADESLOPEOLD(NMAX),
     &  BLADESHAPENEW(NMAX),BLADESLOPENEW(NMAX),
     &  BLADEINILENGTH,BLADESEGMENT
      DOUBLEPRECISION VY1,VY2
      DOUBLEPRECISION FROM,UPTO,GUESS
      INTEGER NGRID,IGRID,IERR
      DOUBLEPRECISION PAR(3)
      DOUBLEPRECISION DX,DL
      DOUBLEPRECISION LOADING1,LOADING2
      COMMON /GRID/ XGRID,BLADESEGMENT,NGRID,IGRID
      COMMON /BLADE/BLADESHAPEOLD,BLADESLOPEOLD,
     &              BLADESHAPENEW,BLADESLOPENEW
      COMMON /PARAMS/ PAR
      COMMON /FLOWDIR/ VY1,VY2
      EXTERNAL DX,DL

C   copy PARAMs in COMMON block
      PAR(1)=PARAM(1)
      PAR(2)=PARAM(2)
      PAR(3)=PARAM(3)
C   compute the initial blade length (for grid generation)
C   check the DL function and compute the BLADESEGMENT length
      CALL DGAUS8(DL,0.D0,1.D0,1.D-12,BLADEINILENGTH,IERR)
      BLADESEGMENT=BLADEINILENGTH/DFLOAT(NGRID-1)
      WRITE(*,'('' initial blade length '',G15.7)') BLADEINILENGTH
C   generate a fancy grid from 0 to 1
C   equal segments along the initial blade (average streamline)
      XGRID(1)=0.D0
      XGRID(NGRID)=1.D0
      DO IGRID=2,NGRID-1
        FROM=XGRID(IGRID-1)
C       careful with UPTO !
        UPTO=XGRID(IGRID-1)+4.D0/DFLOAT(NGRID)
        GUESS=FROM
        CALL DFZERO(DX,FROM,UPTO,GUESS,1.D-10,1.D-10,IERR)
        IF (IERR.NE.1) WRITE(*,'(''node'',I3,'' status FZERO'',I3)')
     &  IGRID,IERR
        XGRID(IGRID)=FROM
      END DO
C   compute the 'streamline' geometry and slope
C      BLADEOLD(2*(1-1)+1)=0.D0  ! blade leading edge at origin
C      BLADEOLD(2*(1-1)+2)=VY1   ! blade slope at leading edge
      DO IGRID=1,NGRID
        BLADESLOPEOLD(IGRID)=VY1+(VY2-VY1)*LOADING1(XGRID(IGRID))
        BLADESHAPEOLD(IGRID)=XGRID(IGRID)*VY1
     &  +(VY2-VY1)*LOADING2(XGRID(IGRID))
      END DO
      RETURN
      END SUBROUTINE BLADESETUP

      DOUBLEPRECISION FUNCTION DL(X)
      IMPLICIT NONE
      DOUBLEPRECISION X,VY1,VY2,LOADING1
      COMMON /FLOWDIR/ VY1,VY2
C      EXTERNAL LOADING1
      DL=VY1+(VY2-VY1)*LOADING1(X)
      DL=DSQRT(1.D0+DL*DL)
      RETURN
      END FUNCTION DL

      DOUBLEPRECISION FUNCTION DX(X)
C   for finding X such that DX=BLADESEGMENT
      IMPLICIT NONE
      INTEGER NGRID,IGRID,IERR
      DOUBLEPRECISION XGRID(100),BLADESEGMENT
      DOUBLEPRECISION X,DL
      COMMON /GRID/ XGRID,BLADESEGMENT,NGRID,IGRID
      EXTERNAL DL
      CALL DGAUS8(DL,XGRID(IGRID-1),X,1.D-12,DX,IERR)
      DX=DX-BLADESEGMENT ! this is the equation
      RETURN
      END FUNCTION DX

      DOUBLEPRECISION FUNCTION LOADING2(X)
      IMPLICIT NONE
      DOUBLEPRECISION X,LOADING1
      INTEGER IERR
      EXTERNAL LOADING1
      CALL DGAUS8(LOADING1,0.D0,X,1.D-12,LOADING2,IERR)
      RETURN
      END FUNCTION LOADING2

      DOUBLEPRECISION FUNCTION LOADING1(X)
      IMPLICIT NONE
      DOUBLEPRECISION X,LOADING
      INTEGER IERR
      EXTERNAL LOADING
      CALL DGAUS8(LOADING,0.D0,X,1.D-12,LOADING1,IERR)
      RETURN
      END FUNCTION LOADING1

      DOUBLEPRECISION FUNCTION LOADING(X)
C   three-parameter loading shape
C   normalized to unit integral
      IMPLICIT NONE
      DOUBLEPRECISION X,LOADING0,LOADNORM
      INTEGER IERR
      EXTERNAL LOADING0
      CALL DGAUS8(LOADING0,0.D0,1.D0,1.D-12,LOADNORM,IERR)
      LOADING=LOADING0(X)/LOADNORM
      RETURN
      END FUNCTION LOADING

      DOUBLEPRECISION FUNCTION LOADING0(X)
C   three-parameter loading shape
C   this is not normalized to unit integral
      IMPLICIT NONE
      DOUBLEPRECISION X,PAR(3)
      COMMON /PARAMS/ PAR
C   three-slopes loading function
      LOADING0=DERF(PAR(1)*X)*DERF(PAR(2)*(1.D0-X))
     &  *(1.D0+PAR(3)*X)
      RETURN
      END FUNCTION LOADING0
