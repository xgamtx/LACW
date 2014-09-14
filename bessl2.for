      DOUBLE PRECISION FUNCTION BESSY(N,X,Eps)
c      implicit double precision(a-h,o-z)
c   -- Calculating of Y-Bessel function --
      INTEGER I,N,NMAX
      PARAMETER (NMAX=2048)
      DOUBLE PRECISION PI,X,Eps,DBESSY(3),GAMMA,BESSJ0,S,TMP,TMP1
      PARAMETER( GAMMA=.577215664901532860606512D0 )
*
      IF(DABS(X).LT.Eps) THEN
        STOP '*ERR* BESSY: argument must be not zero.'
      ENDIF
      PI=4.D0*DATAN(1.D0)
      S=0.D0
      DO I=1, NMAX
        CALL BESS(2*I,X,TMP,TMP1)
        TMP=TMP/DBLE(I)
        IF(MOD(I,2).EQ.0) THEN
          S=S+TMP
         ELSE
          S=S-TMP
        ENDIF
        IF(DABS(TMP).LT.Eps) EXIT
      END DO
      IF(I.GE.NMAX) THEN
        STOP '*ERR* Max. number of points exceeded.'
      ENDIF
      CALL BESS(0,X,BESSJ0,TMP1)
      DBESSY(1)=2.D0/PI*(DLOG(X/2.D0)+GAMMA)*BESSJ0-4.D0/PI*S   ! Y0(X)
      CALL BESS(1,X,TMP,TMP1)
      DBESSY(2)=(TMP*DBESSY(1)-2.D0/PI/X)/BESSJ0                ! Y1(X)
      DO I=3, N+1
        DBESSY(3)=DBLE(2*(I-2))*DBESSY(2)/X-DBESSY(1)
        DBESSY(1)=DBESSY(2)
        DBESSY(2)=DBESSY(3)
      END DO
      BESSY=DBESSY(3)
      RETURN
      END

*-------------------------------------------------------------------------

      SUBROUTINE BESSYPROC(N,X,Eps,Y,Y1)
      implicit double precision(a-h,o-z)
c   -- Calculating of Y-Bessel function --
      INTEGER I,N,NMAX
      double precision Y,Y1
* no external array to have NMAX as max index
*$ old fragment:      PARAMETER (NMAX=101)
*$ new fragment:
      PARAMETER (NMAX=2048)
*$ end fragment
      DOUBLE PRECISION PI,X,Eps,DBESSY(3),GAMMA,BESSJ0,S,TMP,TMP1
      PARAMETER( GAMMA=.577215664901532860606512D0 )
*
      IF(DABS(X).LT.Eps) THEN
d        print *,'BessYproc: n,x,eps:',n,x,eps
        STOP '*ERR* BESSYPROC: argument must be not zero.'
      ENDIF
      PI=4.D0*DATAN(1.D0)
      S=0.D0
      DO I=1, NMAX
        CALL BESS(2*I,X,TMP,TMP1)
        TMP=TMP/DBLE(I)
        IF(MOD(I,2).EQ.0) THEN
          S=S+TMP
         ELSE
          S=S-TMP
        ENDIF
        IF(DABS(TMP).LT.Eps) EXIT
      END DO
* no external array to have NMAX as max index
      IF(I.GE.NMAX) THEN
        STOP '*ERR* Max. number of points exceeded.'
      ENDIF
      CALL BESS(0,X,BESSJ0,TMP1)
      DBESSY(2)=2.D0/PI*(DLOG(X/2.D0)+GAMMA)*BESSJ0-4.D0/PI*S   ! Y0(X)
      CALL BESS(1,X,TMP,TMP1)
      DBESSY(3)=(TMP*DBESSY(2)-2.D0/PI/X)/BESSJ0                ! Y1(X)
      DBESSY(1)=-DBESSY(3)
      DO I=1,iabs(N)
        DBESSY(1)=DBESSY(2)
        DBESSY(2)=DBESSY(3)
        DBESSY(3)=DBLE(2*I)*DBESSY(2)/X-DBESSY(1)
      END DO
      IF((ISIGN(1,N).LT.0) .AND. MOD(N,2).NE.0) THEN
         do i=1,3
            DBESSY(i)=-DBESSY(i)        ! Y(-m,x)=Y(m,x)*(-1)^m
         enddo
      ENDIF
      Y=DBESSY(2)
      Y1=(DBESSY(1)-DBESSY(3))/2.0d0
      RETURN
      END
