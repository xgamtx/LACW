      SUBROUTINE BESS_NUL(N,M,SOL,NDIM,iBesselFamily)
c       Here N - count of values, M - power of BESSEL J function
      INTEGER N,M,AM,I,Md,NDIM
      integer iBesselFamily	! 1==J, 2==Y
      DOUBLE PRECISION STEP,FIND0,Eps
      PARAMETER(Eps=1.0D-7)
      EXTERNAL FIND0
      DOUBLE PRECISION SOL(NDIM),XA,XB,BJA,BJB,DUMMY
      double precision FSTSOL(20,2)
      DATA FSTSOL/2.4048255577D0,5.5200781103D0,8.6537279129D0,
     .		 11.7915344391D0,14.9309177086D0,18.0710639679D0,
     .		 21.2116366299D0, 24.3524715308D0,27.4934791320D0,
     .		 30.6346064684D0,33.7758202136D0,36.9170983537D0,
     .		 40.0584257646D0, 43.1997917132D0,46.3411883717D0,
     .		 49.4826098974D0,52.6240518411D0,55.7655107550D0,
     .		 58.9069839261D0,62.0484691902D0,
     . 0.89357697D0,3.95767842D0,7.08605106D0,10.22234504D0,
     . 13.36109747D0,16.50092244D0,19.64130970D0,22.78202805D0,
     . 25.92295765D0,29.06403025D0,32.20520412D0,35.34645231D0,
     . 38.48775665D0,41.62910447D0,44.77048661D0,47.91189633D0,
     . 51.05332855D0,54.19477936D0,57.33624570D0,60.47772516D0/
      DATA STEP/.5D0/
*
      if (iBesselFamily.ne.1.and.iBesselFamily.ne.2)
     . stop 'Bess_nul.for: Unknown family of Bessel function'
*
d      print *,'	m	J0 roots	Y0 roots'
d      do i=1,20
d      print *,i,FSTSOL(i,1),FSTSOL(i,2)
d      enddo
d      stop '*TEST* Array traverse'
*
*      IF(M.EQ.0) STOP '**BESS_NUL**> Nothing to do!'
	AM=IABS(M)
      IF(N+AM.GT.NDIM) THEN
	 WRITE(*,1) NDIM,N+AM
 	 STOP
      ENDIF
      DO I=1,20
        IF(I.GT.AM+N) EXIT
        SOL(I)=FSTSOL(I,iBesselFamily)
      END DO
      XA=SOL(I-1)+Eps
      XB=XA+STEP
      DO WHILE(I.LE.AM+N)
        if (iBesselFamily.eq.1) then
	  CALL BESS(0,XA,BJA,DUMMY)
	  CALL BESS(0,XB,BJB,DUMMY)
	else
	  CALL BESSYPROC(0,XA,Eps,BJA,DUMMY)
	  CALL BESSYPROC(0,XB,Eps,BJB,DUMMY)
d	  BJA=BESSY(0,XA,Eps)
d	  BJB=BESSY(0,XB,Eps)
	endif
	  IF(BJA*BJB.LT.0D0) THEN
	    SOL(I)=FIND0(XA,XB,0,0,iBesselFamily)
	    I=I+1
	  ENDIF
	  XA=XB+Eps
	  XB=XB+STEP-Eps
      END DO
      Md=0
      DO I=1,AM
	 IF(I.EQ.AM) Md=1
	  DO J=1,AM+N-I
	    VAL=FIND0(SOL(J),SOL(J+1),I*ISIGN(1,M),Md,iBesselFamily)
	    SOL(J)=VAL
	  END DO
      END DO
      RETURN
   1  FORMAT(' **BESS_NUL**> max. dimension of array of solutions',I3,
     .       ', VALUE:',i3)
      END
*
      DOUBLE PRECISION FUNCTION FIND0(A,B,N,Mode,iBesselFamily)
      INTEGER N,Mode
      integer iBesselFamily	! 1==J, 2==Y
c       Here N - power of BESSEL J function
      DOUBLE PRECISION A,B,VAL,F,F1,F2,DUMMY,X,XA,XB
      DOUBLE PRECISION EpsHigh,EpsLow,Eps
      PARAMETER(EpsLow=1.D-4,EpsHigh=1.D-11)
      XA=A
      XB=B
      IF(Mode.EQ.0) THEN
	Eps=EpsLow
		    ELSE
	Eps=EpsHigh
      ENDIF
      if (iBesselFamily.eq.1) then
        CALL BESS(N,XA,F1,DUMMY)
        CALL BESS(N,XB,F2,DUMMY)
      else
	CALL BESSYPROC(N,XA,Eps,F1,DUMMY)
	CALL BESSYPROC(N,XB,Eps,F2,DUMMY)
d	F1=BESSY(N,XA,Eps)
d	F2=BESSY(N,XB,Eps)
	endif
      X=XA
      VAL=XB
      DO WHILE (DABS(X-VAL).GT.Eps)
	X=VAL
	VAL=XA-(XB-XA)/(F2-F1)*F1
        if (iBesselFamily.eq.1) then
	  CALL BESS(N,VAL,F,DUMMY)
        else
	  CALL BESSYPROC(N,VAL,Eps,F,DUMMY)
d	  F=BESSY(N,VAL,Eps)
        endif
	IF(F*F1.LT.0D0) THEN
	  XB=VAL
	  F2=F
        ELSE
	  XA=VAL
	  F1=F
	ENDIF
      END DO
      FIND0=(VAL+X)/2.D0
      RETURN
      END
