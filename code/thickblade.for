C   add thickness to a camberline (thin blade)
C   use airfoil thickness from airfoiltools.com
C   e.g. NACA 16-018,
C   http://airfoiltools.com/airfoil/details?airfoil=naca16018-il
      PROGRAM THICKBLADE
C   first read the thinblade data file,
C   XGRID(IGRID), ! axial coordinate from 0 to 1
C   BLADESHAPEOLD(IGRID), ! vanishing spacing f(x)
C   BLADESHAPENEW(IGRID), ! final blade shape f(x)
C   BLADESLOPEOLD(IGRID), ! vanishing spacing f'(x)
C   BLADESLOPENEW(IGRID)  ! final blade slope f'(x)

      IMPLICIT NONE
      INTEGER IGRID,NGRID,IO,NFOIL,IFOIL,IBLADE
      INTEGER NMAX,NBLADE,INCFD,IERR
      PARAMETER (NMAX=100,NBLADE=150,INCFD=1)
      LOGICAL SKIP
      DOUBLEPRECISION XGRID(NMAX),
     &  BLADESHAPEOLD(NMAX),BLADESLOPEOLD(NMAX),
     &  BLADESHAPENEW(NMAX),BLADESLOPENEW(NMAX)
      DOUBLEPRECISION XFOIL(NMAX),YFOIL(NMAX),DYFOIL(NMAX)
      DOUBLEPRECISION THINBLADESHAPE(INCFD,NMAX),
     &  THINBLADESLOPE(INCFD,NMAX)
      DOUBLEPRECISION XBLADE(4*NMAX),
     &  CAMBERSHAPE(4*NMAX),CAMBERSLOPE(4*NMAX),
     &  THICKSHAPE(4*NMAX),THICKSLOPE(4*NMAX),
     &  XBLADE1(4*NMAX),YBLADE1(4*NMAX),DYBLADE1(4*NMAX),
     &  XBLADE2(4*NMAX),YBLADE2(4*NMAX),DYBLADE2(4*NMAX),
     &  THICKU1(4*NMAX),THICKU2(4*NMAX)

      DOUBLEPRECISION PI
      PARAMETER (PI=4.D0*DATAN(1.D0))
      DOUBLEPRECISION SCALING
      PARAMETER (SCALING=1.D-1/9.D-2)
      WRITE (*,*) 'PI=',PI

      OPEN(1,FILE='AXENT-stator-thinblade.data',STATUS='OLD')
      NGRID=1 ! initialize the number of points
      DO
        READ(1,*,IOSTAT=IO) XGRID(NGRID),
     &  BLADESHAPEOLD(NGRID),BLADESHAPENEW(NGRID),
     &  BLADESLOPEOLD(NGRID),BLADESLOPENEW(NGRID)
        IF (IO.GT.0) THEN
            WRITE(*,*) 'Check input. Something was wrong'
            EXIT
        ELSE IF (IO.LT.0) THEN
C           there was nothing to read after last increment
            NGRID=NGRID-1
            WRITE(*,*) 'number of camber points', NGRID
            GO TO 100
        ELSE
            NGRID=NGRID+1 ! increment the number of points
        END IF
      END DO
      CLOSE(1)
100   CONTINUE
C   check file reading: camberline shape and slope
      DO IGRID=1,NGRID
        THINBLADESHAPE(INCFD,IGRID)=BLADESHAPENEW(IGRID)
        THINBLADESLOPE(INCFD,IGRID)=BLADESLOPENEW(IGRID)
        WRITE(*,'(I3.3,4G15.7)') IGRID,XGRID(IGRID),
     &  BLADESHAPENEW(IGRID),BLADESLOPENEW(IGRID),
     &  1.D0/DSQRT(1.0D0+BLADESLOPENEW(IGRID)**2)
      END DO

C   read thickness distribution
      OPEN(1,FILE='YS-900.data',STATUS='OLD')
C   http://airfoiltools.com/airfoil/lednicerdatfile?airfoil=naca16015-il
C   http://airfoiltools.com/airfoil/lednicerdatfile?airfoil=naca16018-il
C   http://airfoiltools.com/airfoil/details?airfoil=ys900-il
C   R. Eppler and Y. T. Shen
C   WING SECTIONS FOR HYDROFOILS--PART 1: SYMMETRICAL PROFILES
C   Journal of Ship Research, vol.23, no.3, pp.209-2017, 1979
      NFOIL=1 ! initialize number of points on the foil
C   Selig data file TE-upperside-LE-lowerside-TE
      DO
        READ(1,*,IOSTAT=IO) XFOIL(NFOIL),YFOIL(NFOIL)
        IF (IO.GT.0) THEN
            WRITE(*,*) 'Check input. Something was wrong'
            EXIT
        ELSE IF (IO.LT.0) THEN
C           there was nothing to read after last increment
            NFOIL=NFOIL-1
            WRITE(*,*) 'number of thickness points', NFOIL
            GO TO 200
        ELSE
            NFOIL=NFOIL+1 ! increment the number of points
        END IF
      END DO
200   CONTINUE
C   check file reading: XFOIL and YFOIL
      DO IFOIL=1,NFOIL
        WRITE(*,'(I3.3,2G15.7)') IFOIL,XFOIL(IFOIL),YFOIL(IFOIL)
      END DO
      CLOSE(1)

C   SUBROUTINE PCHFD (N, X, F, D, INCFD, SKIP, NE, XE, FE, DE, IERR)
C   Evaluate a piecewise cubic Hermite function and its first
C   derivative at an array of points. Used for CAMBERLINE interpolation

C   build a grid for blade geometry
      DO IBLADE=1,NBLADE
        XBLADE(IBLADE)=DCOS(PI*DFLOAT(IBLADE-1)/DFLOAT(NBLADE-1))
        XBLADE(IBLADE)=0.5D0*(1.D0-XBLADE(IBLADE))
C        WRITE(*,'(I4.3,G15.7)') IBLADE,XBLADE(IBLADE)
      END DO

      SKIP=.TRUE.
C   SKIP -- (input/output) logical variable which should be set to
C           .TRUE. if the user wishes to skip checks for validity of
C           preceding parameters, or to .FALSE. otherwise.
C           This will save time in case these checks have already
C           been performed (say, in PCHIM or PCHIC).
C           SKIP will be set to .TRUE. on normal return.
      CALL DPCHFD(NGRID,XGRID,THINBLADESHAPE,THINBLADESLOPE,
     &  INCFD,SKIP,NBLADE,XBLADE,CAMBERSHAPE,CAMBERSLOPE,IERR)
      WRITE(*,*) 'IERR from PCHFD camberline ',IERR
      OPEN(1,FILE='camberline.dat')
      DO IBLADE=1,NBLADE
        WRITE(1,'(3G15.7)') XBLADE(IBLADE),
     &  CAMBERSHAPE(IBLADE),CAMBERSLOPE(IBLADE)
      END DO
      CLOSE(1)

C   interpolate the thickness distribution
      CALL DPCHIM(NFOIL,XFOIL,YFOIL,DYFOIL,INCFD,IERR)
      WRITE(*,*) 'IERR from PCHIM thickness ',IERR
C   IERR.GT.0  means that IERR switches in the direction
C   of monotonicity were detected.
      SKIP=.TRUE.
      CALL DPCHFD(NFOIL,XFOIL,YFOIL,DYFOIL,
     &  INCFD,SKIP,NBLADE,XBLADE,THICKSHAPE,THICKSLOPE,IERR)
      WRITE(*,*) 'IERR from PCHFD thickness ',IERR
      OPEN(1,FILE='thickness.dat')
      DO IBLADE=1,NBLADE
        WRITE(1,'(2G15.7)') XBLADE(IBLADE),THICKSHAPE(IBLADE)
     &  /DSQRT(1.D0+CAMBERSLOPE(IBLADE)**2)
      END DO
      CLOSE(1)

C   add thickness to the camberline
      DO IBLADE=1,NBLADE
        XBLADE1(IBLADE)=XBLADE(IBLADE)
     &  -SCALING*THICKSHAPE(IBLADE)*CAMBERSLOPE(IBLADE)
     &  /(1.D0+CAMBERSLOPE(IBLADE)**2)
        YBLADE1(IBLADE)=CAMBERSHAPE(IBLADE)
     &  +SCALING*THICKSHAPE(IBLADE)
     &  /(1.D0+CAMBERSLOPE(IBLADE)**2)
        XBLADE2(IBLADE)=XBLADE(IBLADE)
     &  +SCALING*THICKSHAPE(IBLADE)*CAMBERSLOPE(IBLADE)
     &  /(1.D0+CAMBERSLOPE(IBLADE)**2)
        YBLADE2(IBLADE)=CAMBERSHAPE(IBLADE)
     &  -SCALING*THICKSHAPE(IBLADE)
     &  /(1.D0+CAMBERSLOPE(IBLADE)**2)
        WRITE(*,'(4G15.7)') XBLADE1(IBLADE),YBLADE1(IBLADE),
     &  XBLADE2(IBLADE),YBLADE2(IBLADE)
      END DO
      OPEN(1,FILE='blade1.dat')
      OPEN(2,FILE='blade2.dat')
      DO IBLADE=1,NBLADE
        WRITE(1,'(2G15.7)') XBLADE1(IBLADE),YBLADE1(IBLADE)
        WRITE(2,'(2G15.7)') XBLADE2(IBLADE),YBLADE2(IBLADE)
      END DO
      CLOSE(1)
      CLOSE(2)

C   check tangential thickness vs. symmetrical hydrofoil
      CALL DPCHIM(NBLADE,XBLADE1,YBLADE1,DYBLADE1,INCFD,IERR)
      WRITE(*,*) 'IERR from PCHIM BLADE1 ',IERR
      CALL DPCHIM(NBLADE,XBLADE2,YBLADE2,DYBLADE2,INCFD,IERR)
      WRITE(*,*) 'IERR from PCHIM BLADE2 ',IERR
      SKIP=.FALSE.
      CALL DPCHFE(NBLADE,XBLADE1,YBLADE1,DYBLADE1,
     &  INCFD,SKIP,NBLADE,XBLADE,THICKU1,IERR)
      WRITE(*,*) 'IERR from PCHFE BLADE1 ',IERR
      SKIP=.FALSE.
      CALL DPCHFE(NBLADE,XBLADE2,YBLADE2,DYBLADE2,
     &  INCFD,SKIP,NBLADE,XBLADE,THICKU2,IERR)
      WRITE(*,*) 'IERR from PCHFE BLADE2 ',IERR
      OPEN(1,FILE='thicktang.dat')
      DO IBLADE=1,NBLADE
        WRITE(1,'(3G15.7)') XBLADE(IBLADE),
     &  SCALING*2.D0*THICKSHAPE(IBLADE), ! prescribed thickness
     &  THICKU1(IBLADE)-THICKU2(IBLADE) ! actual thickness
      END DO
      CLOSE(1)

      CALL GAMBITJOURNAL('thickstator.jou',
     &  NGRID,XGRID,BLADESHAPEOLD,
     &  NBLADE,XBLADE1,YBLADE1,XBLADE2,YBLADE2)

      END PROGRAM THICKBLADE

      SUBROUTINE GAMBITJOURNAL(FILENAME,
     &  NGRID,XGRID,BLADESHAPEOLD,
     &  NBLADE,XBLADE1,YBLADE1,XBLADE2,YBLADE2)
C   write journal file for GAMBIT
      IMPLICIT NONE
      CHARACTER(*) FILENAME ! journal file name
      INTEGER NGRID,NBLADE,I
      DOUBLEPRECISION XGRID(*),BLADESHAPEOLD(*),
     &  XBLADE1(*),YBLADE1(*),XBLADE2(*),YBLADE2(*)

      OPEN(1,FILE=FILENAME)

C   periodic up
      DO I=1,NGRID
        WRITE(1,
     &  '(''vertex create "ps'',I4.4,''" coordinates'',2G15.7,F4.1)'
     &  ) I+1000,XGRID(I),BLADESHAPEOLD(I)+0.5D0,0.D0
      END DO
      WRITE(1,'(''edge create "periodic-up-s" nurbs \'')')
      DO I=1,NGRID
        WRITE(1,'(''"ps'',I4.4,''" \'')') I+1000
      END DO
      WRITE(1,'(''interpolate'')')
C   periodic lo
      DO I=1,NGRID
        WRITE(1,
     &  '(''vertex create "ps'',I4.4,''" coordinates'',2G15.7,F4.1)'
     &  ) I+2000,XGRID(I),BLADESHAPEOLD(I)-0.5D0,0.D0
      END DO
      WRITE(1,'(''edge create "periodic-lo-s" nurbs \'')')
      DO I=1,NGRID
        WRITE(1,'(''"ps'',I4.4,''" \'')') I+2000
      END DO
      WRITE(1,'(''interpolate'')')
C   write LE and TE
      WRITE(1,
     &  '(''vertex create "pLEs" coordinates'',2G15.7,F4.1)')
     &  XBLADE1(1),YBLADE1(1),0.D0 ! leading edge
      WRITE(1,
     &  '(''vertex create "pTEs" coordinates'',2G15.7,F4.1)')
     &  XBLADE1(NBLADE),YBLADE1(NBLADE),0.D0 ! trailing edge
C   generate points on the blade upper and lower sides
      DO I=2,NBLADE-1
        WRITE(1,
     &  '(''vertex create "ps'',I4.4,''" coordinates'',2G15.7,F4.1)'
     &  ) I+3000,XBLADE1(I),YBLADE1(I),0.D0
        WRITE(1,
     &  '(''vertex create "ps'',I4.4,''" coordinates'',2G15.7,F4.1)'
     &  ) I+4000,XBLADE2(I),YBLADE2(I),0.D0
      END DO
C   generate blade upper side
      WRITE(1,'(''edge create "blade-up-s" nurbs \'')')
      WRITE(1,'(''"pLEs" \'')')
      DO I=2,NBLADE-1
        WRITE(1,'(''"ps'',I4.4,''" \'')') I+3000
      END DO
      WRITE(1,'(''"pTEs" \'')')
      WRITE(1,'(''interpolate'')')
C   generate blade lower side
      WRITE(1,'(''edge create "blade-lo" nurbs \'')')
      WRITE(1,'(''"pLEs" \'')')
      DO I=2,NBLADE-1
        WRITE(1,'(''"ps'',I4.4,''" \'')') I+4000
      END DO
      WRITE(1,'(''"pTEs" \'')')
      WRITE(1,'(''interpolate'')')

      CLOSE(1)
      END SUBROUTINE GAMBITJOURNAL
