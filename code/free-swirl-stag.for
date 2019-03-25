        PROGRAM PIPEFREESWIRL
C solve cubic equation for stagnant region radius
C R.F. Susan-Resiga A variational model for swirling flow
C states with stagnant region Eur J Mech B/Fluids 55 104-115 (2016)
        IMPLICIT NONE
        DOUBLEPRECISION A0,A1,A2,A3,X1,X2,X3,S,Ssq
        INTEGER L,I
        OPEN (1,FILE='free-swirl-stag.dat')
        DO I=1,100
            S=DFLOAT(I)*1.D-2
C           S**2/X**2-2/(1-X)**3=0, Eq.(22)
            Ssq=S*S
            A3=Ssq
            A2=2.D0-3.D0*Ssq
            A1=3.D0*Ssq
            A0=-Ssq
            CALL TC(A3,A2,A1,A0,X1,X2,X3,L)
C           use only X1
            WRITE(1,'(2G15.7)') S,DSQRT(X1)
            PRINT '(4G15.7,I2)',S,X1,X2,X3,L
        END DO
        CLOSE(1)
        OPEN (1,FILE='free-swirl-stag-small.dat')
        DO I=0,100
            S=DFLOAT(I)*1.D-3
            Ssq=S*S
            A2=2.D0-3.D0*Ssq
            A1=3.D0*Ssq
            A0=-Ssq
C           quadratic approximation
            X1=(-A1+DSQRT(A1*A1-4.D0*A2*A0))/(2.D0*A2)
            X2=(-A1-DSQRT(A1*A1-4.D0*A2*A0))/(2.D0*A2)
            WRITE(1,'(2G15.7)') S,DSQRT(X1)
            PRINT '(3G15.7)',S,X1,X2
        END DO
        CLOSE(1)
        END PROGRAM PIPEFREESWIRL

        SUBROUTINE TC(A,B,C,D,X1,X2,X3,L)
C Lebedev V.I., On formulae for roots of cubic equation
C Sov.J.Numer.Anal.Math Modeling 6(4) pp.315-324 (1991)
        IMPLICIT NONE
        INTEGER L
        DOUBLEPRECISION A,B,C,D,X1,X2,X3,
     *  T,S,T1,T2,T3,T4,P
        T=DSQRT(3.D0)
        S=1.D0/3.D0
        T2=B*B
        T3=3.D0*A
        T4=T3*C
        P=T2-T4
        X3=DABS(P)
        X3=DSQRT(X3)
        X1=B*(T4-P-P)-3.D0*T3*T3*D
        X2=DABS(X1)
        X2=X2**S
        T2=1.D0/T3
        T3=B*T2
        IF(X3.GT.1.D-32*X2) GO TO 1
        IF(X1.LT.0.D0) X2=-X2
        X1=X2*T2
        X2=-.5D0*X1
        X3=-T*X2
        IF(DABS(X3).GT.1.D-32) GO TO 15
        X3=X2
        GO TO 2
    1   T1=.5D0*X1/(P*X3)
        X2=DABS(T1)
        T2=X3*T2
        T=T*T2
        T4=X2*X2
        IF(P.LT.0.D0) GO TO 7
        X3=DABS(1.D0-T4)
        X3=DSQRT(X3)
        IF(T4.GT.1.D0) GO TO 5
        T4=DATAN2(X3,T1)*S
        X3=DCOS(T4)
        T4=DSQRT(1.D0-X3*X3)*T
        X3=X3*T2
        X1=X3+X3
        X2=T4-X3
        X3=-(T4+X3)
        IF(X2.GT.X3) GO TO 2
        T2=X2
        X2=X3
        X3=T2
    2   L=0
        IF(X1.GT.X2) GO TO 3
        T2=X1
        X1=X2
        X2=T2
        IF(T2.GT.X3) GO TO 3
        X2=X3
        X3=T2
    3   X3=X3-T3
        GO TO 20
    5   P=(X2+X3)**S
        T4=1.D0/P
        IF(T1) 11,13,13
    7   P=X2+DSQRT(T4+1.D0)
        P=P**S
        T4=-1.D0/P
        IF(T1.LT.0.D0) GO TO 13
   11   T2=-T2
   13   X1=(P+T4)*T2
        X2=-.5D0*X1
        X3=.5D0*T*(P-T4)
   15   L=2
   20   X1=X1-T3
        X2=X2-T3
        RETURN
        END SUBROUTINE TC
