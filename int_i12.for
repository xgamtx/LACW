      DOUBLE PRECISION FUNCTION FI1(RAL,VKKP,KPNM,TH,M,L,I1_2,R_I)
*         Sub-integral values for I1 and I2
      LOGICAL R_I
      DOUBLE PRECISION RAL,VKKP,KPNM,TH,DUMMY,BJM,BJM1,TMP,Z,X,Y,PLM_S,
     .                 TMP1
      INTEGER M,L,I1_2
      IF(I1_2.EQ.1) THEN
*     Sub_integral I1 - function value (real and img. parts)
        TMP=RAL*VKKP*DCOS(TH)
        IF(R_I) THEN
           TMP1=DCOS(TMP)
         ELSE
           TMP1=DSIN(TMP)
        ENDIF
        CALL BESS(M,KPNM*RAL*DSIN(TH),BJM,DUMMY)
        FI1=TMP1*BJM*plgndr(L,IABS(M),DCOS(TH))*DSIN(TH)
      ENDIF
      IF(I1_2.EQ.2) THEN
c   -- Subintegral I2 - function value (real and img. parts) --
        TMP=KPNM*RAL*DSIN(TH)
        CALL BESS(M,TMP,BJM,BJM1)
        Z=KPNM*DSIN(TH)*BJM1
        Y=VKKP*DCOS(TH)
        X=Y*RAL
        PLM_S=plgndr(L,IABS(M),DCOS(TH))*DSIN(TH)
        IF(R_I) THEN
          FI1=PLM_S*(Z*DCOS(X)-BJM*Y*DSIN(X)) ! Real
         ELSE
          FI1=PLM_S*(Z*DSIN(X)+BJM*Y*DCOS(X)) ! Imag.
        ENDIF
        RETURN
      ENDIF
      END
*
      DOUBLE PRECISION FUNCTION SIMP_I1(Eps,RAL,VKP,KP,M,L,A,B,
     .                                  I1_2,R_I)
c -- Integr. I1 - calc. by Simpson's method. Uses FI1 as sub-inegral function --
      INTEGER JMAX,L,M,I1_2,I
      LOGICAL R_I
      DOUBLE PRECISION RAL,VKP,KP,Eps,A,B,s1,S_Old,H,Hprev,X,FI1,TMP
      EXTERNAL FI1
      PARAMETER (JMAX=100)
*
      H=(B-A)/2.D0
      TMP=2.d0*FI1(RAL,VKP,KP,a+h,M,L,I1_2,R_I)
      S1 = FI1(RAL,VKP,KP,A,M,L,I1_2,R_I)+TMP+
     .     FI1(RAL,VKP,KP,B,M,L,I1_2,R_I)
      SIMP_I1 = S1+TMP
      DO I=1, JMAX
        TMP=0.D0
        S_Old=SIMP_I1
        Hprev=H
        H=Hprev/2.d0
        DO X=A+H, B, Hprev
          TMP=TMP+FI1(RAL,VKP,KP,X,M,L,I1_2,R_I)
        END DO
        TMP=2.D0*TMP
        S1=S1+TMP       ! Saved part of integral for use on next iteration
        SIMP_I1=S1+TMP  ! Value of integral on current iteration
        IF(DABS(2.D0*S_Old-SIMP_I1)*H/3.D0 .LT. Eps) EXIT
      END DO
      IF(I.GE.JMAX) THEN
         WRITE(*,'(1X,A32,I5)') 'SIMP_I1---> limit of iterations',
     +         JMAX
         STOP
      ENDIF
      SIMP_I1=SIMP_I1*H/3.D0
      RETURN
      END
*
      DOUBLE COMPLEX FUNCTION SIMP_I12(Eps,RAL,VKKP,KPN,M,L,I1_2)
c      Integr. I1 - calc. by Simpson's method. Uses FI1,
c      calls SIMP_I1 and convert result to complex*16.
      LOGICAL RorI
      INTEGER L, M, I1_2
      DOUBLE PRECISION RAL,VKKP,KPN,Eps,FNC,PI
      DOUBLE PRECISION SIMP_I1
      EXTERNAL SIMP_I1
      PARAMETER(PI=3.14159265358979323846264338328D0)
      RorI=MOD(L+M,2).EQ.0      ! odd or even?
      IF(.NOT.RorI .AND. DABS(VKKP).LT.Eps) THEN
       FNC=0.D0
      ELSE
       FNC=SIMP_I1(Eps,RAL,VKKP,KPN,M,L,0.D0,PI/2.D0,I1_2,RorI)
      ENDIF
      IF(RorI) THEN
        SIMP_I12=DCMPLX(2.d0*FNC,0.D0)
      ELSE
        SIMP_I12=DCMPLX(0.d0,2.d0*FNC)
      ENDIF
      RETURN
      END
