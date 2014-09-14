      SUBROUTINE BESS(N,X,BESSLJ,BESSLJ1)
*   Calculating J-Bessel function and first derivative: BESSLJ & BESSLJ1
      INTEGER I,N,M,NM,NF,NMAX
      PARAMETER (NMAX=2048)
      DOUBLE PRECISION X,BESSJ(NMAX),PMIN,S,BESSLJ,BESSLJ1
      PARAMETER(PMIN=1.D-12)
*
d      print *,'Bess(): n,x:',n,x
      M=IABS(N)+1
      IF(DABS(X).LT.PMIN) THEN
	  NF=1
	  DO I=1,M-1
	    NF=NF*I
	  END DO
        DO I=1, M+2
	 IF(M-1.NE.0) THEN
	    BESSJ(i)=(X/2.D0)**(M-1)/DBLE(NF)
	   ELSE
	    BESSJ(I)=0.D0
	  ENDIF
        END DO
      ELSE
        NM=MAX0(M*2,2*IDNINT(X),20)
        IF(NM.GT.NMAX) THEN
	    WRITE(*,99) NM,NMAX
	    STOP
        ENDIF
        S=0.D0
        BESSJ(NM)=1	!  = MAX. value
        BESSJ(NM+1)=0	!  = MIN. value
        DO I=NM, 2, -1
	    BESSJ(I-1)=DBLE(2*(I-1))*BESSJ(I)/X-BESSJ(I+1)
        END DO
        DO I=3, NM+1, 2
	    S=S+2.D0*BESSJ(I)
        END DO
        S=S+BESSJ(1)
        do i=1, M+2
	    BESSJ(I)=BESSJ(i)/S
c            print *,'b1',i,bessj(i),s
        end do
      ENDIF
      if(dabs(x).lt.PMIN.and.mod(M-1,2).eq.1) besslj=-besslj  !by IMSL
      IF(M.EQ.1) THEN
	BESSLJ1=-BESSJ(2)
       ELSE
	BESSLJ1=(BESSJ(M-1)-BESSJ(M+1))/2.D0
      ENDIF
      BESSLJ=BESSJ(M)
      if(dabs(x).lt.PMIN .AND. N.EQ.0) BESSLJ=1.D0
      IF((ISIGN(1,N).LT.0) .AND. MOD(N,2).NE.0) THEN
         BESSLJ=-BESSLJ        ! (-1)**M for Bessel J-function ...
	 BESSLJ1=-BESSLJ1      ! and its derivative
      ENDIF
      RETURN
  99  FORMAT(3X,'**BESS**> WARNING! MAX. DIM. INDEX=',I5,' > ',
     .	     'MAX=',I3)
      END
