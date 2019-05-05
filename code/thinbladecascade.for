C   compile and link with SLATEC and LAPACK libraries
C   -L /opt/local/lib -lslatec -llapack

C   thin blade cascade design for minimum maximum velocity

      PROGRAM THINBLADECASCADE
C   use the BOBYQA algorithm for bound constrained
C       optimization without derivatives (DAMTP 2009.NA06)
C   SUBROUTINE BOBYQA (N,NPT,X,XL,XU,RHOBEG,RHOEND,IPRINT,MAXFUN,W)
C   IMPLICIT REAL*8 (A-H,O-Z)
C   DIMENSION X(*),XL(*),XU(*),W(*)
C   This subroutine seeks the least value of a function of many variables,
C   by applying a trust region method that forms quadratic models by
C   interpolation. There is usually some freedom in the interpolation
C   conditions, which is taken up by minimizing the Frobenius norm of
C   the change to the second derivative of the model, beginning with the
C   zero matrix. The values of the variables are constrained by upper and
C   lower bounds.
C   SUBROUTINE CALFUN (N,X,F) has to be provided by the user. It must set
C   F to the value of the objective function for the current values of the
C   variables X(1),X(2),...,X(N), which are generated automatically in a
C   way that satisfies the bounds given in XL and XU.
      IMPLICIT NONE
      INTEGER NMAX
      PARAMETER (NMAX=100)
      INTEGER NGRID,IGRID,ITER
      DOUBLEPRECISION XGRID(NMAX),BLADESEGMENT,
     &  BLADESHAPEOLD(NMAX),BLADESLOPEOLD(NMAX),
     &  BLADESHAPENEW(NMAX),BLADESLOPENEW(NMAX),
     &  PAR(3),VY1,VY2,S,
     &  TWOPI,LOADING,LOADING1,LOADING2
      PARAMETER (TWOPI=6.283185307179586D0)
      DOUBLEPRECISION RELERR
      DOUBLEPRECISION VBL(NMAX),VBLUP(NMAX),VBLLO(NMAX)

      INTEGER NPARAM,IPRINT,NPT,MAXFUN
      PARAMETER (NPARAM=3)
      DOUBLEPRECISION PARAM(NPARAM),VELMAX,
     &  PARAMLO(NPARAM),PARAMUP(NPARAM)
      DOUBLEPRECISION RHOBEG,RHOEND,W(500)
C   length of W must be at least (NPT+5)*(NPT+N)+3*N*(N+5)/2


      COMMON /GRID/ XGRID,BLADESEGMENT,NGRID,IGRID
      COMMON /BLADE/BLADESHAPEOLD,BLADESLOPEOLD,
     &              BLADESHAPENEW,BLADESLOPENEW
      COMMON /PARAMS/ PAR
      COMMON /FLOWDIR/ VY1,VY2
      COMMON /SPACING/ S ! cascade spacing
      COMMON /VELOCITY/ VBL,VBLLO,VBLUP

C   AXENT stator blades problem setup
c      VY1=0.0D0 ! AXENT stator inlet
c      VY2=14.7D0/6.8D0 ! AXENT stator outlet
c      S=1.05D0 ! from ZWEIFEL criterion

C   AXENT runner blades problem setup
      VY1=(14.7D0-16.0D0)/6.8D0 ! AXENT runner inlet
      VY2=-16.0D0/6.8D0 ! AXENT runner outlet
      S=1.03D0 ! from ZWEIFEL criterion

C   grid refinement
      NGRID=76 ! grid refinement

C   solution for stator blade
C   PARAM(1)= 1.821132D1 ! slope LE
C   PARAM(2)= 3.202183D0 ! slope TE
C   PARAM(3)= 1.459788D0 ! slope middle
C   VELMAX  = 2.509556D0 ! maximum velocity on the suction side

C   initialize parameters
c      PARAM(1)= 1.820242D1 ! slope LE
c      PARAM(2)= 3.201883D0 ! slope TE
c      PARAM(3)= 1.461346D0 ! slope middle
      PARAM(1)= 1.530426D1 ! slope LE
      PARAM(2)= 2.219675D0 ! slope TE
      PARAM(3)= 1.227281D0 ! slope middle

C   upper and lower limits
      PARAMLO(1) = 1.0D1
      PARAMUP(1) = 2.0D1
      PARAMLO(2) = 1.D0
      PARAMUP(2) = 5.D1
      PARAMLO(3) = 1.0D0
      PARAMUP(3) = 2.0D0

C   setup gridpoints and blade (shape & slope)
      CALL BLADESETUP(PARAM)
      WRITE(*,'('' done setup problem '')')

C   test CALFUN
      CALL CALFUN(NPARAM,PARAM,VELMAX)
      WRITE(*,'('' maximum velocity '',G15.7)') VELMAX

      IPRINT=3
C   NPT must be in the interval [N+2,(N+1)(N+2)/2], i.e. [5,10]
C   choices that exceed 2*N+1 are not recommended
      NPT=2*NPARAM+1 ! NPT=7 in our case
C   RHOBEG and RHOEND must be set to the initial and final values
C   of a trust region radius, so both must be positive with
C   RHOEND no greater than RHOBEG.
C   Typically, RHOBEG should be about one tenth of the greatest
C   expected change to a variable, while RHOEND should indicate the
C   accuracy that is required in the final values of the variables.
C   An error return occurs if any of the differences XU(I)-XL(I),
C   I=1,...,N, is less than 2*RHOBEG.
      RHOBEG=1.D-2
      RHOEND=1.D-4
      MAXFUN=200
      CALL BOBYQA (NPARAM,NPT,PARAM,PARAMLO,PARAMUP,
     &  RHOBEG,RHOEND,IPRINT,MAXFUN,W)
      CALL CALFUN(NPARAM,PARAM,VELMAX)


      WRITE(*,'('' maximum velocity '',G15.7)') VELMAX
      WRITE(*,'('' slope LE  '',G15.7)') PAR(1)
      WRITE(*,'('' slope TE  '',G15.7)') PAR(2)
      WRITE(*,'('' slope mid '',G15.7)') PAR(3)

      OPEN(1,FILE='loading.dat')
      OPEN(2,FILE='blade.dat')
      OPEN(3,FILE='velocity.dat')
      OPEN(4,FILE='cp_lo.dat')
      OPEN(5,FILE='cp_up.dat')
      DO IGRID=1,NGRID
        WRITE(1,'(2G15.7)') XGRID(IGRID),LOADING(XGRID(IGRID))
C   plot the initial and the final blade
        BLADESHAPEOLD(IGRID)=XGRID(IGRID)*VY1
     &  +(VY2-VY1)*LOADING2(XGRID(IGRID))
        BLADESLOPEOLD(IGRID)=VY1+(VY2-VY1)*LOADING1(XGRID(IGRID))
        WRITE(2,'(5G15.7)') XGRID(IGRID),
     &  BLADESHAPEOLD(IGRID), ! vanishing spacing f(x)
     &  BLADESHAPENEW(IGRID), ! final blade shape f(x)
     &  BLADESLOPEOLD(IGRID), ! vanishing spacing f'(x)
     &  BLADESLOPENEW(IGRID)  ! final blade slope f'(x)
        WRITE(3,'(3G15.7)') XGRID(IGRID),
     &  VBLLO(IGRID), ! velocity on the lower side of the blade
     &  VBLUP(IGRID)  ! velocity on the upper side of the blade
C   pressure coefficient with respect to downstream conditions
C   downstream velocity V2**2=1+VY2**2
C   cp = 1 - V**2/V2**2
        WRITE(4,'(2G15.7)') XGRID(IGRID),
     &  (1.D0-VBLLO(IGRID)**2/(1.D0+VY2**2))
        WRITE(5,'(2G15.7)') XGRID(IGRID),
     &  (1.D0-VBLUP(IGRID)**2/(1.D0+VY2**2))
      END DO
      CLOSE(1)
      CLOSE(2)
      CLOSE(3)
      CLOSE(4)
      CLOSE(5)

      CALL GEOMETRY()

C   miscelaneous
C   write and compare
      DO IGRID=1,NGRID
        WRITE(*,'(5G15.7)') XGRID(IGRID),
     &  BLADESHAPENEW(IGRID),BLADESLOPENEW(IGRID),
     &  VBLLO(IGRID),VBLUP(IGRID)
      END DO
C   check geometrical approximation for leading edge slope
      WRITE(*,'('' LE slope exact = '',G15.7,'' LE aprox = '',G15.7)')
     &  BLADESLOPENEW(1),2.D0*BLADESHAPENEW(2)/XGRID(2)-BLADESLOPENEW(2)

      STOP
      END PROGRAM THINBLADECASCADE

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE GEOMETRY()
C   generate lines for FLUEND/GAMBIT geometry
      IMPLICIT NONE
      INTEGER NMAX,NGRID,IGRID
      PARAMETER (NMAX=100)
      DOUBLEPRECISION XGRID(NMAX),BLADESEGMENT,
     &  BLADESHAPEOLD(NMAX),BLADESLOPEOLD(NMAX),
     &  BLADESHAPENEW(NMAX),BLADESLOPENEW(NMAX),S,
     &  THICKNESS,LOADING1,LOADING2,VY1,VY2
      PARAMETER (THICKNESS=5.D-4)
      DOUBLEPRECISION X,FX,DFX,XPS,YPS,XSS,YSS
      COMMON /GRID/ XGRID,BLADESEGMENT,NGRID,IGRID
      COMMON /BLADE/BLADESHAPEOLD,BLADESLOPEOLD,
     &              BLADESHAPENEW,BLADESLOPENEW
      COMMON /SPACING/ S ! cascade spacing
      COMMON /FLOWDIR/ VY1,VY2

      OPEN (1,FILE='blade_ps.dat')
      OPEN (2,FILE='blade_ss.dat')
      OPEN (3,FILE='periodic_up.dat')
      OPEN (4,FILE='periodic_lo.dat')
C   write leading edge
      WRITE(1,'(2G15.7)') XGRID(1),BLADESHAPENEW(1)
      WRITE(2,'(2G15.7)') XGRID(1),BLADESHAPENEW(1)
C   write interior points
      DO IGRID=2,NGRID-1
        X=XGRID(IGRID)
        FX=BLADESHAPENEW(IGRID)
        DFX=BLADESLOPENEW(IGRID)
        XPS=X -THICKNESS*DFX/DSQRT(1.D0+DFX*DFX)
        YPS=FX+THICKNESS/DSQRT(1.D0+DFX*DFX)
        XSS=X+THICKNESS*DFX/DSQRT(1.D0+DFX*DFX)
        YSS=FX-THICKNESS/DSQRT(1.D0+DFX*DFX)
        WRITE(1,'(2G15.7)') XPS,YPS
        WRITE(2,'(2G15.7)') XSS,YSS
      END DO
C   write trailing edge
      WRITE(1,'(2G15.7)') XGRID(NGRID),BLADESHAPENEW(NGRID)
      WRITE(2,'(2G15.7)') XGRID(NGRID),BLADESHAPENEW(NGRID)
      CLOSE(1)
      CLOSE(2)
      DO IGRID=1,NGRID
        X=XGRID(IGRID)
        FX=X*VY1+(VY2-VY1)*LOADING2(X)
        WRITE(3,'(2G15.7)') X,FX+S/2.D0
        WRITE(4,'(2G15.7)') X,FX-S/2.D0
      END DO
      CLOSE(3)
      CLOSE(4)

C   write *.jou file for GAMBIT
      OPEN (1,FILE='thinbladecascade.jou')
C   write upper periodic boundary
      DO IGRID=1,NGRID
        X=XGRID(IGRID)
        FX=X*VY1+(VY2-VY1)*LOADING2(X)
        WRITE(1,
     &  '(''vertex create "p'',I3.3,''" coordinates'',2G15.7,F4.1)')
     &  IGRID, X,FX+0.5D0*S,0.D0 ! upper periodic
      END DO
      WRITE(1,'(''edge create "periodic-up" nurbs \'')')
      DO IGRID=1,NGRID
        WRITE(1,'(''"p'',I3.3,''" \'')') IGRID
      END DO
      WRITE(1,'(''interpolate'')')
C   write lower periodic boundary
      DO IGRID=1,NGRID
        X=XGRID(IGRID)
        FX=X*VY1+(VY2-VY1)*LOADING2(X)
        WRITE(1,
     &  '(''vertex create "p'',I3.3,''" coordinates'',2G15.7,F4.1)')
     &  IGRID+NGRID, X,FX-0.5D0*S,0. ! upper periodic
      END DO
      WRITE(1,'(''edge create "periodic-lo" nurbs \'')')
      DO IGRID=1,NGRID
        WRITE(1,'(''"p'',I3.3,''" \'')') IGRID+NGRID
      END DO
      WRITE(1,'(''interpolate'')')
C   write blade LE and TE
      WRITE(1,
     &  '(''vertex create "pLE" coordinates'',2G15.7,F4.1)')
     &  XGRID(1),BLADESHAPENEW(1),0. ! leading edge
      WRITE(1,
     &  '(''vertex create "pTE" coordinates'',2G15.7,F4.1)')
     &  XGRID(NGRID),BLADESHAPENEW(NGRID),0. ! leading edge
C   write blade upper side
      DO IGRID=2,NGRID-1
        X=XGRID(IGRID)
        FX=BLADESHAPENEW(IGRID)
        DFX=BLADESLOPENEW(IGRID)
        XPS=X -THICKNESS*DFX/DSQRT(1.D0+DFX*DFX)
        YPS=FX+THICKNESS/DSQRT(1.D0+DFX*DFX)
        WRITE(1,
     &  '(''vertex create "p'',I3.3,''" coordinates'',2G15.7,F4.1)')
     &  IGRID+300,XPS,YPS,0. ! blade upper side
      END DO
      WRITE(1,'(''edge create "blade-up" nurbs \'')')
      WRITE(1,'(''"pLE" \'')')
      DO IGRID=2,NGRID-1
        WRITE(1,'(''"p'',I3.3,''" \'')') IGRID+300
      END DO
      WRITE(1,'(''"pTE" \'')')
      WRITE(1,'(''interpolate'')')
C   write blade lower side
      DO IGRID=2,NGRID-1
        X=XGRID(IGRID)
        FX=BLADESHAPENEW(IGRID)
        DFX=BLADESLOPENEW(IGRID)
        XSS=X +THICKNESS*DFX/DSQRT(1.D0+DFX*DFX)
        YSS=FX-THICKNESS/DSQRT(1.D0+DFX*DFX)
        WRITE(1,
     &  '(''vertex create "p'',I3.3,''" coordinates'',2G15.7,F4.1)')
     &  IGRID+400,XSS,YSS,0. ! blade lower side
      END DO
      WRITE(1,'(''edge create "blade-lo" nurbs \'')')
      WRITE(1,'(''"pLE" \'')')
      DO IGRID=2,NGRID-1
        WRITE(1,'(''"p'',I3.3,''" \'')') IGRID+400
      END DO
      WRITE(1,'(''"pTE" \'')')
      WRITE(1,'(''interpolate'')')
      CLOSE(1)

      RETURN
      END SUBROUTINE GEOMETRY


      SUBROUTINE CALFUN(N,X,F)
C   set F to the value of the objective function for the
C   current values of the variables X(1),X(2),X(3) which
C   are generated automatically in a way that satisfies the
C   bounds given in XL and XU
      IMPLICIT NONE
      INTEGER N
      DOUBLEPRECISION X,F
      DIMENSION X(*)
      INTEGER NMAX
      PARAMETER (NMAX=100)
      DOUBLEPRECISION PAR,VBL,VBLLO,VBLUP
      DIMENSION PAR(3),VBL(NMAX),VBLLO(NMAX),VBLUP(NMAX)
      INTEGER NGRID,IGRID,ITER
      DOUBLEPRECISION XGRID,BLADESEGMENT,RELERR,ERRMAX
      PARAMETER (ERRMAX=1.D-5)
      DIMENSION XGRID(NMAX)
      DOUBLEPRECISION BLADESHAPEOLD,BLADESLOPEOLD,
     &  BLADESHAPENEW,BLADESLOPENEW
      DIMENSION BLADESHAPEOLD(NMAX),BLADESLOPEOLD(NMAX),
     &  BLADESHAPENEW(NMAX),BLADESLOPENEW(NMAX)
      COMMON /GRID/ XGRID,BLADESEGMENT,NGRID,IGRID
      COMMON /PARAMS/ PAR
      COMMON /VELOCITY/ VBL,VBLLO,VBLUP
      COMMON /BLADE/ BLADESHAPEOLD,BLADESLOPEOLD,
     &  BLADESHAPENEW,BLADESLOPENEW
C   copy the parameters in COMMON block
      PAR(1)=X(1)
      PAR(2)=X(2)
      PAR(3)=X(3)
C   compute the grid and initialize blade shape & slope
C      CALL BLADESETUP() ! set XGRID & initial BLADESHAPE BLADESLOPE
C   iteratively correct the BLADESHAPE and BLADESLOPE
      RELERR=1.D0
      ITER=0
      DO WHILE(RELERR.GT.ERRMAX)
        ITER=ITER+1
        CALL BLADEUPDATE(RELERR)
C        WRITE(*,'('' RELERR = '',I4,G15.7)') ITER,RELERR
        DO IGRID=1,NGRID
            BLADESHAPEOLD(IGRID)=BLADESHAPENEW(IGRID)
            BLADESLOPEOLD(IGRID)=BLADESLOPENEW(IGRID)
        END DO
      IF (ITER.GT.30) GO TO 1
      END DO
    1 CALL BLADEUPDATE(RELERR) ! one more final update
      WRITE(*,'('' RELERR = '',I4,G15.7)') ITER,RELERR
C   compute the velocity on the blade
      CALL BLADEVELOCITY()
C   find the maximum velocity (objective function)
      F=0.D0 ! initialize maximum velocity value
      DO IGRID=1,NGRID
        IF (VBLLO(IGRID).GT.F) F=VBLLO(IGRID)
        IF (VBLUP(IGRID).GT.F) F=VBLUP(IGRID)
      END DO
      RETURN
      END SUBROUTINE CALFUN
