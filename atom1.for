      SUBROUTINE POTAT(NCY,UNIT,NERRO,NFLOW)
* ---------------------------------------------------------
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      CHARACTER*10 XNAME
      INTEGER UNIT
* ----------------------------------------------------------
      COMMON  / B L A T /
     *  VR(300),RH(300),A(300),RAD(300),Y(300),B(300),
     *  A0(50),B0(50),XN(50),XL(50),
     *  XJ(50),XE(50),XZ(50),XDELL(50),DA(5),DB(5),VOC(5),
     *  RN,H,ZEFF,ZN,XION,PHI,EPS,DEL,DELRV,XALPH,XLATT,RNUC,
     *  CONVR,FN,FL,FJ,E,Q,EV,EZ,EY,V0,ET3,SUMMA,ANUC,FDLL,RBAR,
     *  VBAR,H3,FK,CS,G,Q11,Q22,HA1,HA2,HA3,HA4,HA5,HA6,HA7,HA8,
     *  HA9,TCS,ITSCF,INDSCF,NPRIN,N,J,NC1,IDIRC,NXT,NDEBU
      COMMON /TEXT1/ XNAME
      COMMON /TIP/ CORDNS1(50),EPSCOR1(50),TOP(300),EPEDE1,EPEDE2,
     * V_old(300),V_real(300),ICORE
* -------------------------------------------------------
      DATA CL/137.037D0/,D1/1.D0/,D12/0.5D0/,DO/0.0D0/
* -------------------------------------------------------
      NERRO=0
      ZEFF=XION+XLATT
      EPSLO=DO
*
      DO 1 K=1,NFLOW
      R=RAD(K)
      IF(R.LT.RNUC) THEN
         X=R/RNUC
         RVN=-ZN*X*(1.5D0-D12*X*X)
      ELSE
         RVN=-ZN
      ENDIF
*
      FLOW=A(K)*DSQRT(R)*D12
      TOP(K)=FLOW
      RV=FLOW+Y(K)+RVN
*
      IF(R.LE.RBAR) THEN
         RV=RV+R*VBAR
      ELSE
         RV=RV+RBAR*VBAR
         IF(RV+ZEFF.GT.DO) RV=-ZEFF
      ENDIF
      ERROR=DABS(RV-VR(K))
      IF(ERROR.GT.EPSLO) THEN
         EPSLO = ERROR
         NERRO = K
      ENDIF
      VR(K) = RV
    1 CONTINUE
*
      CONVR=DABS(EPSLO)
      RA=RAD(1)
      RB=RAD(2)/RA
      RC=RAD(3)/RA
      IF(RNUC.LE.DO) THEN
         VOC(1)=-ZN
         TA=(VR(1)+ZN)
         TB=(VR(2)+ZN)/RB
         TC=(VR(3)+ZN)/RC
      ELSE
         VOC(1)=DO
         TA=VR(1)
         TB=VR(2)/RB
         TC=VR(3)/RC
      ENDIF
      RA=D1
      DETA=RB*RC*(RC-RB)
      DETB=RA*RC*(RA-RC)
      DETC=RA*RB*(RB-RA)
      DET=DETA+DETB+DETC
      VOC(5)=(TA*DETA+TB*DETB+TC*DETC)/DET
      DETA=RA*RA
      DETB=RB*RB
      DETC=RC*RC
      VOC(3)=(TA*(DETB-DETC)+TB*(DETC-DETA)+TC*(DETA-DETB))/DET
      VOC(4)=(TA*(RC-RB)+TB*(RA-RC)+TC*(RB-RA))/DET
      IF(NDEBU.EQ.2)
     A   WRITE (UNIT,101) VOC(1),VOC(5),VOC(3),VOC(4)
  101 FORMAT(1X,'* * POTENTIAL EXPANSION COEFFICIENTS * *',/,2X,
     * ' VOC(1)=',D15.4,' VOC(2)=',D15.4,/,2X,
     * ' VOC(3)=',D15.4,' VOC(4)=',D15.4,/)
      IF(IDIRC.NE.0) THEN
        DO 2 I=1,5
    2   VOC(I)=-VOC(I)/CL
      ENDIF
      RETURN
      END
*
      SUBROUTINE HSDAT(UNIT,XNORM)
* ---------------------------------------------------------
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      CHARACTER*10 XNAME
      INTEGER UNIT
* ---------------------------------------------------------
      COMMON  / B L A T /
     *  VR(300),RH(300),A(300),RAD(300),Y(300),B(300),
     *  A0(50),B0(50),XN(50),XL(50),
     *  XJ(50),XE(50),XZ(50),XDELL(50),DA(5),DB(5),VOC(5),
     *  RN,H,ZEFF,ZN,XION,PHI,EPS,DEL,DELRV,XALPH,XLATT,RNUC,
     *  CONVR,FN,FL,FJ,E,Q,EV,EZ,EY,V0,ET3,SUMMA,ANUC,FDLL,RBAR,
     *  VBAR,H3,FK,CS,G,Q11,Q22,HA1,HA2,HA3,HA4,HA5,HA6,HA7,HA8,
     *  HA9,TCS,ITSCF,INDSCF,NPRIN,N,J,NC1,IDIRC,NXT,NDEBU
      COMMON /TEXT1/ XNAME
* ---------------------------------------------------------
      DATA DELT/0.1D0/,DELK/0.1D-01/
      DATA D1/1.D0/,D12/0.5D0/,ETOP/1.D-05/,
     *     D2/2.D0/,SMALL/.1D-18/,D14/0.25D0/,
     *     DO/0.0D0/,DD18/.1125D1/,D18/.125D0/
      DATA FMODC/.92561983D00/,FCORC/.07438017D00/
* ---------------------------------------------------------
      DE=DO
      EOLD=E
      NTEST=0
      FDLL=DELT*FDLL
*     DELS=DMAX1(FDLL,(DEL*(D1-E)))
      IF(FDLL.GE.(DEL*(D1-E))) THEN
          DELS=FDLL
      ELSE
          DELS=(DEL*(D1-E))
      ENDIF
*
      EMIN=-D12*(ZN/FN)**2
      EMAX=-1.D-5
      FLP=FL+D1
      TLPO=FLP+FL
      TFLP=FLP+FLP
      NOMAX=IDNINT(FN-FLP)
      IF(EMAX.LT.E.OR.E.LT.EMIN) E=D12*(EMAX+EMIN)
      A0(1)=D1
      DG=DEXP(H*FLP)
      RST=RN**FLP/DG**N
*
    1 RA=RAD(1)
      VOC(2)=VOC(5)-E*RA
      DO 2 I=1,4
    2 B0(I)=D2*VOC(I)
*
      FIS=D1
      FFL=TFLP
      DO 3 IS=1,10
      SUMA=DO
      DO 4 IT=1,(MIN0(IS,4))
    4 SUMA=SUMA+B0(IT)*A0(IS-IT+1)
      A0(IS+1)=RA*SUMA/(FIS*FFL)
      FIS=FIS+D1
    3 FFL=FFL+D1
*
      DO 5 I=1,N
    5 B(I)=D2*RAD(I)*(VR(I)-RAD(I)*E)
*
      DO 6 I=1,4
      SUMA=DO
      SUMB=DO
      FIS=DO
      TPOW=D1
      DO 7 IS=1,11
      SUMA=SUMA+A0(IS)*TPOW
      SUMB=SUMB+FIS*A0(IS)*TPOW
      FIS=FIS+D1
    7 TPOW=TPOW*(RAD(I)/RA)
      A(I)=SUMA
      DA(I)=SUMB
    6 DB(I)=-TLPO*DA(I)+B(I)*A(I)
*
      NODES=0
      DO 8 K=11,(N-10)
      KM=N-K
      IF(B(KM).LT.DO) GO TO 9
    8 CONTINUE
*
    9 KI=KM+1
      ADIF=DO
      BDIF=DO
      DO 10 K=5,KI
      R=RAD(K)
      BPRED=DA(1)+HA1*(DB(4)-D12*DB(3)+DB(2))
      APRED=A(K-4)+HA1*(DA(4)-D12*DA(3)+DA(2))
      BMODE=BPRED+FMODC*BDIF
      AMODE=APRED+FMODC*ADIF
      FMODE=-TLPO*BMODE+B(K)*AMODE
      BCORR=DD18*DA(4)-D18*DA(2)+HA2*(FMODE+D2*DB(4)-DB(3))
      ACORR=DD18*A(K-1)-D18*A(K-3)+HA2*(BMODE+D2*DA(4)-DA(3))
      ADIF=ACORR-APRED
      BDIF=BCORR-BPRED
      A(K)=ACORR-FCORC*ADIF
      DA(5)=BCORR-FCORC*BDIF
      DB(5)=-TLPO*DA(5)+B(K)*A(K)
      DA(1)=DA(2)
      DA(2)=DA(3)
      DA(3)=DA(4)
      DA(4)=DA(5)
      DB(2)=DB(3)
      DB(3)=DB(4)
      DB(4)=DB(5)
      IF(A(K)*A(K-1).LT.DO) NODES=NODES+1
   10 CONTINUE
*
      IF(NOMAX.NE.NODES) THEN
         IF(NOMAX.LT.NODES) THEN
            IF(E.LT.EMAX) EMAX=E
*           E=E+D12*DMAX1((E+E),(EMIN-E))
            IF((E+E).GE.(EMIN-E)) THEN
                 E=(E+E)
            ELSE
                 E=E+D12*(EMIN-E)
            ENDIF
*
            DL1=DABS(EMAX-EMIN)/(DABS(E)+D1)
            IF(DL1.LT.DELK) THEN
               WRITE (UNIT,101) NODES,NOMAX
               EMIN=1.2D0*EMIN
            ENDIF
         ELSE
            IF(E.GT.EMIN) EMIN=E
            E=D12*(E+EMAX)
            DL1=DABS(EMAX-EMIN)/(DABS(E)+D1)
            IF(DL1.GE.DELK) GO TO 1
            WRITE (UNIT,101) NODES,NOMAX,FN,FL
            IF(NTEST.NE.0) THEN
               WRITE(UNIT,102) NTEST,FN,FL
               STOP
            ENDIF
            NTEST=1
            EMAX=1.D-5
            E=D12*E+EMAX
         ENDIF
         GO TO 1
      ENDIF
  101 FORMAT(5X,'** WARNING - NUMBER OF NODES FOUND =' ,I3,/,5X,
     * 20X,' REQUIRED=',I3,/,5X,'SHELL N=',F3.1,' L=',F3.1,/)
  102 FORMAT(5X,'** AUTO STOP **',/,5X,'ERROR NTEST=',I2,/,5X,
     *   'SHELL N=',F3.1,' L=',F3.1,/)
*
      RA=A(KI)
      RB=DA(5)
   11 DO 12 K=1,N
      IF(EPS.LT.B(K)) GO TO 14
   12 KJ=K
*
   14 K=KJ
      DO 15 I=1,4
      R=RAD(K)
      A(K)=R**(IDNINT(FN))*DEXP(-R*DSQRT(D2*DABS(VR(K)/R-E)))
      IF(A(K).LT.SMALL) THEN
         EPS=EPS*D14
         GO TO 11
      ENDIF
      DA(I)=A(K)*(FN/R-D1)
      DB(I)=-TLPO*DA(I)+B(K)*A(K)
   15 K=K-1
*
      ADIF=DO
      BDIF=DO
   16 R=RAD(K)
      BPRED=DA(1)-HA1*(DB(4)-D12*DB(3)+DB(2))
      APRED=A(K+4)-HA1*(DA(4)-D12*DA(3)+DA(2))
      BMODE=BPRED+FMODC*BDIF
      AMODE=APRED+FMODC*ADIF
      FMODE=-TLPO*BMODE+B(K)*AMODE
      BCORR=DD18*DA(4)-D18*DA(2)-HA2*(FMODE+D2*DB(4)-DB(3))
      ACORR=DD18*A(K+1)-D18*A(K+3)-HA2*(BMODE+D2*DA(4)-DA(3))
      ADIF=ACORR-APRED
      BDIF=BCORR-BPRED
      A(K)=ACORR-FCORC*ADIF
      DA(5)=BCORR-FCORC*BDIF
      DB(5)=-TLPO*DA(5)+B(K)*A(K)
      IF(K.GT.KI) THEN
         K=K-1
         DA(1)=DA(2)
         DA(2)=DA(3)
         DA(3)=DA(4)
         DA(4)=DA(5)
         DB(2)=DB(3)
         DB(3)=DB(4)
         DB(4)=DB(5)
         GO TO 16
      ENDIF
      RA=RA/A(KI)
      DO 17 K=KI,KJ
   17 A(K)=A(K)*RA
      RC=DA(5)*RA
      RG=RST
      DO 18 K=1,KJ
      RG=RG*DG
   18 A(K)=A(K)*RG
      IF(KJ.LT.N) THEN
         RJ=RAD(KJ)
         RT2E=-DSQRT(-D2*E)
*        RT2E=-((-D2*E)**0.5D0)
         DO 19 K=KJ,N
         IF(((RT2E*(RAD(K)-RJ))+DLOG(DABS(A(KJ)))).LT.-75.D0) THEN
              DO 23 K1=K,N
   23         A(K1)=DO
              GO TO 24
         ENDIF
         A(K)=A(KJ)*DEXP(RT2E*(RAD(K)-RJ))
   19    CONTINUE
      ENDIF
*
   24 DO 20 K=1,N
   20 Y(K)=A(K)*A(K)
*
      W=ATINT(RAD,Y,N,H,3,D2*FLP+D1)
      DE=-D12*RAD(KI)**FL*A(KI)*(RC-RB)/W
      DL1=DABS(EMAX-EMIN)/(DABS(E)+D1)
      DL=DABS(DE)
*     DLL=DMIN1(DL,(-D12*E))
      IF(DL.GE.(-D12*E)) THEN
          DLL=(-D12*E)
      ELSE
          DLL=DL
      ENDIF
*
      IF(DL.GT.DELK.AND.DL1.LT.DELK) THEN
         WRITE (UNIT,103) FN,FL,EMIN,E,EMAX,DE
         GO TO 21
      ENDIF
  103 FORMAT(5X,'** WARNING - <DE> TOO LARGE <SHELL N=',F3.1,' L=',
     * F3.1,/,5X,' EMIN=',F10.5,' E=',F10.5,' EMAX=',F10.5,
     * ' DE=',F10.5,/)
*
      IF(DE.NE.DO) THEN
         IF(DE.GT.DO) THEN
            EMIN=E
            DE=DLL
         ELSE
            EMAX=E
            DE=-DLL
         ENDIF
      ENDIF
      E=E+DE
*     DELL=DMAX1((DELT*DABS(EOLD-E)),DELS)
      IF((DELT*DABS(EOLD-E)).GE.DELS) THEN
            DELL=(DELT*DABS(EOLD-E))
      ELSE
            DELL=DELS
      ENDIF
*
      IF(DL.GE.DELL) THEN
         IF(E.GE.EMAX) E=D12*(E-DE+EMAX)
         IF(E.LE.EMIN) E=D12*(E-DE+EMIN)
         GO TO 1
      ENDIF
*
   21 XNORM=D1/DSQRT(W)
      DO 22 K=1,N
   22 A(K)=A(K)*XNORM
*
      FDLL=DELU/DELT
      RETURN
      END
*
      SUBROUTINE DIRAT(UNIT,XNORM,FZ)
* ----------------------------------------------------------
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      CHARACTER*10 XNAME
      INTEGER UNIT
* ----------------------------------------------------------
      COMMON  / B L A T /
     *  VR(300),RH(300),A(300),RAD(300),Y(300),B(300),
     *  A0(50),B0(50),XN(50),XL(50),
     *  XJ(50),XE(50),XZ(50),XDELL(50),DA(5),DB(5),VOC(5),
     *  RN,H,ZEFF,ZN,XION,PHI,EPS,DEL,DELRV,XALPH,XLATT,RNUC,
     *  CONVR,FN,FL,FJ,E,Q,EV,EZ,EY,V0,ET3,SUMMA,ANUC,FDLL,RBAR,
     *  VBAR,H3,FK,CS,G,Q11,Q22,HA1,HA2,HA3,HA4,HA5,HA6,HA7,HA8,
     *  HA9,TCS,ITSCF,INDSCF,NPRIN,N,J,NC1,IDIRC,NXT,NDEBU
      COMMON /TEXT1/ XNAME
* ---------------------------------------------------------
      DATA GA1/251.D0/,GA2/646.D0/,GA3/264.D0/,GA4/106.D0/,GA5/19.D0/,
     *  GA6/116.D0/,GA7/496.D0/,GA8/96.D0/,GA9/16.D0/,GA10/4.D0/,
     *  GA11/81.D0/,GA12/306.D0/,GA13/216.D0/,GA14/126./,GA15/9.D0/,
     *  GA16/56.D0/,GA17/256.D0/,GA18/96.D0/,GA19/256.D0/,GA20/56.D0/
      DATA C/137.037D0/
      DATA DELT/0.1D0/,DELK/0.01D0/,EFAC/0.6D0/
      DATA DD/.12D1/,D12/.5D0/,SMALL1/1.D-5/,D2/2.D0/,
     * D1/1.D0/,DO/0.0D0/,DMUL/1.D20/,SMALL/1.D-20/
* ---------------------------------------------------------
      FLOP=D12*FL*(FL+D1)
      EOLD=E
      NTEST=0
*     DELS=DMAX1((DELT*FDLL),(DEL*(D1-E)) )
      IF((DELT*FDLL).GE.(DEL*(D1-E))) THEN
          DELS=(DELT*FDLL)
      ELSE
          DELS=(DEL*(D1-E))
      ENDIF
*
      DE=DO
      EMIN=-D12*(D1+ZN/C)*(ZN/FN)**2+VBAR
      EMAX=-SMALL1
      IF(E.GT.EMAX.OR.E.LT.EMIN) E=D12*(EMAX+EMIN)
      FK=D2*(FL-FJ)*(FJ+D12)
      CS=C
      TCS=C+C
      Q21=VOC(1)
*     G=DSQRT(FK*FK-Q21*Q21)
      G=(FK*FK-Q21*Q21)**0.5D0
      Q11=-G-FK
      Q22=-G+FK
      IF(FK.LT.DO) THEN
         A0(1)=G-FK
         B0(1)=Q21
      ELSE
         A0(1)=-Q21
         B0(1)=G+FK
      ENDIF
      DG=DEXP(H*G)
      RST=RN**G/DG**N
      IFN=IDNINT(FN)
      ISTEP=0
*
    1 CALL DIRAT1(KI,IFN1,FLOP,UNIT)
      ISTEP=ISTEP+1
      IF(ISTEP.GT.1000) THEN
         WRITE(UNIT,105) FN,FL,FJ,FZ,ICASE
         WRITE(3,105) FN,FL,FJ,FZ,ICASE
         IF(ICASE.EQ.1) GOTO 17
         IF(ICASE.EQ.2) GOTO 17
         IF(ICASE.EQ.3) GOTO 18
         IF(ISTEP.GT.1001) STOP '<DIRAT>: Check ISTEP'
      ENDIF
  105 FORMAT(5X,'+++ DIRAT emergency ESCAPE',/,5X,'+++ Shell',
     * ' N',F4.1,' L',F4.1,' J',F4.1,' occ',F5.2,/,5X,
     * '+++ Case No',I3,' Keep Running ...',/)
*
      IF(IFN.NE.IFN1) THEN
         IF(IFN.GT.IFN1) THEN
            IF(E.GT.EMIN) EMIN=E
            E=D12*(E+EMAX)
            DL1=DABS(EMAX-EMIN)/(DABS(E)+D1)
            ICASE=1
            IF(DL1.GE.DELK) GOTO 1
            WRITE (UNIT,101) IFN1,FN,FL,FJ,FZ
            IF(NTEST.NE.0) STOP '<DIRAT>: check NTEST'
            NTEST=1
            EMAX=-SMALL1
            E=D12*(E+EMAX)
         ELSE
*           ( FN < FN1 )
            IF(E.LT.EMAX) EMAX=E
*           E=E+D12*(DMAX1(E,(EMIN-E)))
            IF(E.GE.(EMIN-E)) THEN
                E=1.5D0*E
            ELSE
                E=E+D12*(EMIN-E)
            ENDIF
*
            DL1=DABS(EMAX-EMIN)/(DABS(E)+D1)
            ICASE=2
            IF(DL1.GE.DELK) GOTO 1
            WRITE (UNIT,102) IFN1,FN,FL,FJ,FZ
            EMIN=DD*EMIN
         ENDIF
         STOP '<DIRAT>: Internal'
      ENDIF
  101 FORMAT(5X,'** Warning <DIRAT> N1=',I3,' < N=',F3.1,
     *  ' L=',F3.1,' J=',F3.1,' Occ=',F4.2)
  102 FORMAT(5X,'** Warning <DIRAT> N1=',I3,' > N=',F3.1,
     *  ' L=',F3.1,' J=',F3.1,' Occ=',F4.2)
*
   17 CONTINUE
      RA=A(KI)
      RB=B(KI)
      DO 2 I=1,KI
      A(I)=A(I)/RA
    2 B(I)=B(I)/RA
      R=RAD(N)
      KJ=N
      DO 3 K=KI,N
      IF((EPS+RAD(K)*(E*RAD(K)-VR(K))).LT.DO) THEN
         KJ=K-1
         R=RAD(K)
         GOTO 4
      ENDIF
    3 CONTINUE
    4 RZ=-VR(KJ)/R
      RL=(FL+D12)/R
      RK=FK/R
      P=-D2*(E+RZ)+RL*RL
      IF(P.LT.DO) WRITE (UNIT,103) FN,FL,FJ,FZ
  103 FORMAT(5X,'Some <DIRAT> Trouble: N',F4.1,
     * ' L',F4.1,' J',F4.1,' occ',F5.2,/,5X,'Keep Running ...')
      A(KJ)=D1
      B(KJ)=CS*(RK-DSQRT(P)+D12*(RZ-RL*RL)/(R*P))/(-CS*TCS-E-RZ)
*     B(KJ)=CS*(RK-(P)**0.5D0+D12*(RZ-RL*RL)/(R*P))/(-CS*TCS-E-RZ)
      DO 5 L=1,4
      A(KJ-L)=A(KJ)
    5 B(KJ-L)=B(KJ)
      DO 6 I=1,4
      K=KJ+1
      DO 7 L=1,5
      K=K-1
      RP21=(E*RAD(K)-VR(K))/CS
      RP12=-RP21-TCS*RAD(K)
      DA(L)=Q11*A(K)+RP12*B(K)
    7 DB(L)=Q22*B(K)+RP21*A(K)
      A(KJ-1)=A(KJ)-(GA1*DA(1)+GA2*DA(2)-GA3*DA(3)+GA4*DA(4)-
     * GA5*DA(5))*HA6
      B(KJ-1)=B(KJ)-(GA1*DB(1)+GA2*DB(2)-GA3*DB(3)+GA4*DB(4)-
     *  GA5*DB(5))*HA6
      A(KJ-2)=A(KJ)-(GA6*DA(1)+GA7*DA(2)+GA8*DA(3)+GA9*DA(4)-
     * GA10*DA(5))*HA7
      B(KJ-2)=B(KJ)-(GA6*DB(1)+GA7*DB(2)+GA8*DB(3)+GA9*DB(4)-
     * GA10*DB(5))*HA7
      A(KJ-3)=A(KJ)-(GA11*DA(1)+GA12*DA(2)+GA13*DA(3)+GA14*DA(4)-
     * GA15*DA(5))*HA8
      B(KJ-3)=B(KJ)-(GA11*DB(1)+GA12*DB(2)+GA13*DB(3)+GA14*DB(4)-
     * GA15*DB(5))*HA8
      A(KJ-4)=A(KJ)-(GA16*DA(1)+GA17*DA(2)+GA18*DA(3)+GA19*DA(4)+
     * GA20*DA(5))*HA9
      B(KJ-4)=B(KJ)-(GA16*DB(1)+GA17*DB(2)+GA18*DB(3)+GA19*DB(4)+
     * GA20*DB(5))*HA9
    6 CONTINUE
*
      K=KJ-3
      DA(1)=DA(2)
      DB(1)=DB(2)
      DA(2)=DA(3)
      DB(2)=DB(3)
      DA(3)=DA(4)
      DB(3)=DB(4)
    8 K=K-1
      RP21=(E*RAD(K)-VR(K))/CS
      RP12=-RP21-TCS*RAD(K)
      AKK=A(K+4)-HA1*(DA(3)-D12*DA(2)+DA(1))
      BKK=B(K+4)-HA1*(DB(3)-D12*DB(2)+DB(1))
      DA(4)=Q11*AKK+RP12*BKK
      DB(4)=Q22*BKK+RP21*AKK
      A(K)=A(K+1)-HA2*DA(4)-HA3*DA(3)+HA4*DA(2)-HA5*DA(1)
      B(K)=B(K+1)-HA2*DB(4)-HA3*DB(3)+HA4*DB(2)-HA5*DB(1)
      DA(1)=DA(2)
      DB(1)=DB(2)
      DA(2)=DA(3)
      DB(2)=DB(3)
      DA(3)=Q11*A(K)+RP12*B(K)
      DB(3)=Q22*B(K)+RP21*A(K)
      IF(K.GT.KI) GO TO 8
      RC=A(KI)
      RB=RB/B(KI)
      DO 9 K=KI,KJ
      A(K)=A(K)/RC
    9 B(K)=B(K)/RC
      RG=RST
      DO 10 K=1,KJ
      RG=RG*DG
      A(K)=A(K)*RG
   10 B(K)=B(K)*RG
      IF(KJ.EQ.N) GOTO 15
      DO 11 K=KJ,N
      W=(-DSQRT(-D2*E))*(RAD(K)-RAD(KJ))
*     W=(-((-D2*E)**0.5D0))*(RAD(K)-RAD(KJ))
      IF((W+DLOG(DABS(A(KJ)))).LT.-75.0D0) THEN
         DO 16 K1=K,N
         A(K1)=DO
   16    B(K1)=DO
         GOTO 15
      ENDIF
      A(K)=A(KJ)*DEXP(W)
   11 B(K)=B(KJ)*DEXP(W)
*
   15 DO 12 K=1,N
   12 Y(K)=A(K)*A(K)+B(K)*B(K)
      W=ATINT(RAD,Y,N,H,3,D2*G+D1)
      DE=CS*A(KI)*B(KI)*(D1-RB*RC/RA)/W
      DL1=DABS(EMAX-EMIN)/(DABS(E)+D1)
      DL=DABS(DE)
*     DLL=DMIN1(DL,(-D12*E))
      IF(DL.GE.(-D12*E)) THEN
           DLL=-(D12*E)
      ELSE
           DLL=DL
      ENDIF
*
      IF(DL.LE.DELK.OR.DL1.GE.DELK) THEN
         IF(DE.GT.DO) THEN
            EMIN=E
            DE=DLL
         ELSE
            EMAX=E
            DE=-DLL
         ENDIF
         E=E+DE
         DELU=DELT*DABS(EOLD-E)
*        DELL=DMAX1(DELU,DELS)
         IF(DELU.GE.DELS) THEN
            DELL=DELU
         ELSE
            DELL=DELS
         ENDIF
*
         IF(DL.GE.DELL) THEN
            IF(E.GT.EMAX) E=D12*(E-DE+EMAX)
            IF(E.LT.EMIN) E=D12*(E-DE+EMIN)
            ICASE=3
            GOTO 1
         ENDIF
      ELSE
         WRITE (UNIT,104) FN,FL,FJ,FZ
      ENDIF
  104 FORMAT(5X,'Shell Trouble: N',F4.1,' L',F4.1,' J',F4.1,
     * ' occ',F5.2,' Keep Running ...')
*
   18 CONTINUE
      XNORM=D1/DSQRT(W)
*
      DO 14 K=1,N
      A(K)=A(K)*XNORM
   14 B(K)=B(K)*XNORM
      FDLL=DELU/DELT
      RETURN
      END
*
      SUBROUTINE DIRAT1(KI,IFN1,FLOP,UNIT)
* ---------------------------------------------------------
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      CHARACTER*10 XNAME
      INTEGER UNIT
* ---------------------------------------------------------
      COMMON  / B L A T /
     *  VR(300),RH(300),A(300),RAD(300),Y(300),B(300),
     *  A0(50),B0(50),XN(50),XL(50),
     *  XJ(50),XE(50),XZ(50),XDELL(50),DA(5),DB(5),VOC(5),
     *  RN,H,ZEFF,ZN,XION,PHI,EPS,DEL,DELRV,XALPH,XLATT,RNUC,
     *  CONVR,FN,FL,FJ,E,Q,EV,EZ,EY,V0,ET3,SUMMA,ANUC,FDLL,RBAR,
     *  VBAR,H3,FK,CS,G,Q11,Q22,HA1,HA2,HA3,HA4,HA5,HA6,HA7,HA8,
     *  HA9,TCS,ITSCF,INDSCF,NPRIN,N,J,NC1,IDIRC,NXT,NDEBU
      COMMON /TEXT1/ XNAME
* ---------------------------------------------------------
      DATA DO/.0D0/,D12/0.5D0/,D1/1.D0/
* ---------------------------------------------------------
      VOC(2)=VOC(5)+E*RAD(1)/CS
      FIS=DO
      DO 1 IS=2,11
      FIS=FIS+D1
      SUMA=DO
      SUMB=DO
      DO 2 IT=2,MIN0(IS,4)
      ITM=IS-IT+1
      SUMA=SUMA-VOC(IT)*(B0(ITM)*(FIS+G-FK)-VOC(1)*A0(ITM))
    2 SUMB=SUMB+VOC(IT)*(A0(ITM)*(FIS+G+FK)-VOC(1)*B0(ITM))
      A0(IS)=((-TCS*RAD(1))*(FIS+G-FK)*B0(IS-1)+SUMA)/(FIS*(FIS+G+G))
      B0(IS)=((-TCS*RAD(1))*VOC(1)*B0(IS-1)+SUMB)/(FIS*(FIS+G+G))
    1 CONTINUE
*
      DO 3 K=1,5
      A(K)=DO
      B(K)=DO
      TPOW=D1
      DO 4 IS=1,11
      A(K)=A(K)+A0(IS)*TPOW
      B(K)=B(K)+B0(IS)*TPOW
    4 TPOW=TPOW*(RAD(K)/RAD(1))
      RP21=(E*RAD(K)-VR(K))/CS
      RP12=-RP21-TCS*RAD(K)
      IF(K.GT.1) THEN
         DA(K-1)=Q11*A(K)+RP12*B(K)
         DB(K-1)=Q22*B(K)+RP21*A(K)
      ENDIF
    3 CONTINUE
*
      IFN1=1+IDNINT(FL)
      DO 5 KM=(N-11),10, -1
      IF((E*RAD(KM)-VR(KM)-FLOP/RAD(KM)).GT.DO) THEN
         KI=KM+1
         DO 6 K=5,KI
         RP21=(E*RAD(K)-VR(K))/CS
         RP12=-RP21-TCS*RAD(K)
         AKK=A(K-4)+HA1*(DA(3)-D12*DA(2)+DA(1))
         BKK=B(K-4)+HA1*(DB(3)-D12*DB(2)+DB(1))
         DA(4)=Q11*AKK+RP12*BKK
         DB(4)=Q22*BKK+RP21*AKK
         A(K)=A(K-1)+HA2*DA(4)+HA3*DA(3)-HA4*DA(2)+HA5*DA(1)
         B(K)=B(K-1)+HA2*DB(4)+HA3*DB(3)-HA4*DB(2)+HA5*DB(1)
         DA(1)=DA(2)
         DB(1)=DB(2)
         DA(2)=DA(3)
         DB(2)=DB(3)
         DA(3)=Q11*A(K)+RP12*B(K)
         DB(3)=Q22*B(K)+RP21*A(K)
         IF((A(K)*A(K-1)).LT.DO) IFN1=IFN1+1
    6    CONTINUE
         RETURN
      ENDIF
    5 CONTINUE
      WRITE(UNIT,101) FN,FL,FJ
  101 FORMAT(5X,'+ Can not glue solutions: N',F4.1,' L',F4.1,
     * ' J',F4.1,' Keep Running ...')
      RETURN
      END
*
      DOUBLE PRECISION FUNCTION ATINT(X,Y,N,H,NMX,DX)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION X(N),Y(N)
      DATA  D2/2.D0/,D4/4.D0/,DO/0.D0/,D3/3.D0/
*
      START=DO
      IF(NMX.EQ.2) DX=DLOG(X(2)*Y(2)/(X(1)*Y(1)))/H
      IF(NMX.GT.1) START=X(1)*Y(1)/DX
      SUM1=DO
      SUM2=DO
      DO 1 K=2,N,2
      SUM1=SUM1+X(K-1)*Y(K-1)
    1 SUM2=SUM2+X(K)*Y(K)
      ATINT=H*(D2*SUM1+D4*SUM2-X(1)*Y(1))/D3+START
      RETURN
      END
*
      DOUBLE PRECISION FUNCTION ATINT1(X,Y,N,H,NMX,DX)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION X(N),Y(N)
      DATA  D2/2.D0/,D4/4.0D0/,DO/0.0D0/,D3/3.0D0/
      NTEST=MOD(N,2)
      IF(NTEST.EQ.1) THEN
         ISTART=1
         BEG=DO
      ELSE
         ISTART=2
         BEG=(X(1)*Y(1)+X(2)*Y(2))*H/D2
      ENDIF
      SA=DO
      SB=DO
      DO 1 K=ISTART+2,N,2
      SA=SA+X(K)*Y(K)
    1 SB=SB+X(K-1)*Y(K-1)
      ATINT1=H*(D2*SA+D4*SB-X(N)*Y(N)+X(ISTART)*Y(ISTART))/D3+BEG
      RETURN
      END
*
      SUBROUTINE POISAT(PSQ,J,W,Z,XO,H)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
*
      DIMENSION PSQ(J),W(J),E(300),F(300)
      DATA DO/0.D0/,D1/1.D0/,D2/2.D0/,D3/3.D0/,DA/10.D0/,DB/48.D0/
*
      A=D1-H*H/DB
      B=-D2-H*H*DA/DB
      EDL=DEXP(H/D2)
      C=H*H/(D2*D3)
      C2=-B/A
      E(1)=DO
      F(1)=DEXP(H/D2)
      X=XO+H
      ITOP=J-1
      DO 1 I=2,ITOP
      G=C*DEXP(X/D2)*(EDL*PSQ(I+1)+DA*PSQ(I)+PSQ(I-1)/EDL)
      F(I)=C2-D1/F(I-1)
      E(I)=(G/A+E(I-1))/F(I)
   1  X=X+H
      W(J)=D2*Z*DEXP(-X/D2)
      DO 2 I=1,ITOP
      JV=J-I
   2  W(JV)=E(JV)+W(JV+1)/F(JV)
      RETURN
      END
*
      SUBROUTINE AVER(XNORM,I,HSTEP)
* ------------------------------------------------
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER UNIT
* ------------------------------------------------
      COMMON  / B L A T /
     *  VR(300),RH(300),A(300),RAD(300),Y(300),B(300),
     *  A0(50),B0(50),XN(50),XL(50),
     *  XJ(50),XE(50),XZ(50),XDELL(50),DA(5),DB(5),VOC(5),
     *  RN,H,ZEFF,ZN,XION,PHI,EPS,DEL,DELRV,XALPH,XLATT,RNUC,
     *  CONVR,FN,FL,FJ,E,Q,EV,EZ,EY,V0,ET3,SUMMA,ANUC,FDLL,RBAR,
     *  VBAR,H3,FK,CS,G,Q11,Q22,HA1,HA2,HA3,HA4,HA5,HA6,HA7,HA8,
     *  HA9,TCS,ITSCF,INDSCF,NPRIN,N,J,NC1,IDIRC,NXT,NDEBU
      COMMON /ADD/ AMTC(300),CORE(300),FINT(50,5),
     *             RMT,STEP,RSTART,JPR(50),JRI,UNIT
      DIMENSION IPOW(5)
* ------------------------------------------------
      DATA IPOW /-3,-1,1,2,4/,DO/0.D0/,D12/0.5D0/
* ------------------------------------------------
      DO 1 K=1,N
    1 Y(K)=Y(K)*XNORM*XNORM
      IF(XJ(I).EQ.D12) THEN
         ISTART=3
         FINT(I,1)=DO
         FINT(I,2)=DO
      ELSE
         ISTART=1
      ENDIF
      DO 2 IMAIN=ISTART,5
      IUP=IPOW(IMAIN)
         DO 3 K=1,N
    3    A(K)=Y(K)*RAD(K)**IUP
    2 FINT(I,IMAIN)=ATINT1(A,RAD,N,HSTEP,1,XNORM)
      RETURN
      END
*
      SUBROUTINE TOMAS
* -------------------------------------------------
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
* -------------------------------------------------
      COMMON  / B L A T /
     *  VR(300),RH(300),A(300),RAD(300),Y(300),B(300),
     *  A0(50),B0(50),XN(50),XL(50),
     *  XJ(50),XE(50),XZ(50),XDELL(50),DA(5),DB(5),VOC(5),
     *  RN,H,ZEFF,ZN,XION,PHI,EPS,DEL,DELRV,XALPH,XLATT,RNUC,
     *  CONVR,FN,FL,FJ,E,Q,EV,EZ,EY,V0,ET3,SUMMA,ANUC,FDLL,RBAR,
     *  VBAR,H3,FK,CS,G,Q11,Q22,HA1,HA2,HA3,HA4,HA5,HA6,HA7,HA8,
     *  HA9,TCS,ITSCF,INDSCF,NPRIN,N,J,NC1,IDIRC,NXT,NDEBU
* --------------------------------------------------
      DATA DO/0.0D0/,SMALL1/1.0D-07/,D1/1.0D0/,D8/8.0D0/,D2/2.D0/
* --------------------------------------------------
      NFLOW=0
      DO 1  K=1,J
*     B(K)=DSQRT(-D8*XE(K))
      B(K)=(-D8*XE(K))**0.5D0
      NFLOW=IDNINT(D2*(XN(K)-XL(K)))
      DFAC=D1
      DO 2 IK=1,NFLOW
    2 DFAC=DFAC*DBLE(IK)
      A(K)=B(K)**(NFLOW+1)/DFAC
    1 CONTINUE
      DO 3 I=1,N
      SUMDEN=DO
      DO 4 K=1,J
      VAL=RAD(I)*B(K)
      IF(VAL.LT.32.D0) THEN
         SUMDEN=SUMDEN+XZ(K)*A(K)*
     *       DEXP(-VAL)*RAD(I)**(IDNINT(D2*(XN(K)-XL(K))))
      ENDIF
    4 CONTINUE
      RH(I)=SUMDEN
    3 CONTINUE
      RETURN
      END
*
      SUBROUTINE READA1(IN,IOUT)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER UNIT,UNIT1
      CHARACTER*1 TY(5),TFLOW,TDUMP1
      CHARACTER*5 TXT1*30,TXT2,TEND
* ----------------------------------------------------------
      COMMON  / B L A T /
     *  VR(300),RH(300),A(300),RAD(300),Y(300),B(300),
     *  A0(50),B0(50),XN(50),XL(50),
     *  XJ(50),XE(50),XZ(50),XDELL(50),DA(5),DB(5),VOC(5),
     *  RN,H,ZEFF,ZN,XION,PHI,EPS,DEL,DELRV,XALPH,XLATT,RNUC,
     *  CONVR,FN,FL,FJ,E,Q,EV,EZ,EY,V0,ET3,SUMMA,ANUC,FDLL,RBAR,
     *  VBAR,H3,FK,CS,G,Q11,Q22,HA1,HA2,HA3,HA4,HA5,HA6,HA7,HA8,
     *  HA9,TCS,ITSCF,INDSCF,NPRIN,N,J,NC1,IDIRC,NXT,NDEBU
      COMMON /ADD/ AMTC(300),CORE(300),FINT(50,5),
     *             RMT,STEP,RSTART,JPR(50),JRI,UNIT
*-------------------------------------
      DATA TY /'S','P','D','F','G'/, TEND /'  END'/
      UNIT1=16
*-------------------------------------
      IOPEN=0
      DO 1 I=(IN+1),(IN+10)
      IF(IOPEN.EQ.0) THEN
         READ(UNIT1,101)
     *     TXT1,NFLOW,TFLOW,JF1,TDUMP1,JF2,ZFLOW,TXT2
         L1=6
         DO 2 L=1,5
         IF(TFLOW.EQ.TY(L)) L1=L-1
    2    CONTINUE
         IF(L1.EQ.6) THEN
              WRITE(UNIT,102)
              WRITE(UNIT,103) NFLOW,TFLOW,
     *              JF1,TDUMP1,JF2, ZFLOW, TXT2
              STOP
         ENDIF
         XN(I)=DBLE(NFLOW)
         XL(I)=DBLE(L1)
         XJ(I)=DBLE(JF1)/DBLE(JF2)
         XZ(I)=ZFLOW
         JPR(I)=1
         IF(TXT2.EQ.TEND) THEN
             IOPEN=1
             IOUT=I
             RETURN
         ENDIF
      ENDIF
    1 CONTINUE
  101 FORMAT(A30,I1,A1,I1,A1,I1,F5.2,A5)
  102 FORMAT(3X,'<AUTOSTOP IN <READA1> - READING DATA FOR <ATOM>>')
  103 FORMAT(10X,I1,A1,I1,A1,I1,F5.2,A5)
      RETURN
      END
*
      SUBROUTINE HABAT(NXT,N,RAD,IENUN)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      CHARACTER*1 YL(5)
      INTEGER UNIT
*
      COMMON /ADD/ AMTC(300),CORE(300),FINT(50,5),
     *             RMT,STEP,RSTART,JPR(50),JRI,UNIT
      COMMON /HUBB/ DDENS(300,4)
      DIMENSION RAD(N),PPOT(300,4),QFF(4),HABB(4,4),IOPEN(4)
*
      DATA DO/0.D0/,D2/2.D0/,NS/1/,D1/1.D0/,DMIN/1.D-3/
      DATA YL/'S','P','D','F','A'/,EVOLTS/13.6058/
*
      WRITE(UNIT,101)
      WRITE(3,101)
  101 FORMAT(2X,70('-'),/,5X,'ENTER IN <HUBBARD>-<ATOM> BLOCK',/,
     A  5X,'COLOMN INTRACENTER PARAMETERS CALCULATION')
*
      DO 1 IL=1,4
      DO 2 IR=1,N
    2 AMTC(IR)=DDENS(IR,IL)
*
      QFF(IL)=ATINT1(AMTC,RAD,N,STEP,1,STEP)
      IF(QFF(IL).GT.DMIN) THEN
         IOPEN(IL)=1
      ELSE
         IOPEN(IL)=0
         QFF(IL)=DO
         GO TO 1
      ENDIF
      CALL POISAT(AMTC,N,CORE,(QFF(IL)),RSTART,STEP)
*
      DO 3 IR=1,N
    3 PPOT(IR,IL)=CORE(IR)
*
    1 CONTINUE
*
      DO 4 IL=1,4
      DO 4 JL=1,4
      IF((IOPEN(IL)*IOPEN(JL)).NE.0) THEN
*
         DO 5 IR=1,N
*   5    AMTC(IR)=DDENS(IR,JL)*PPOT(IR,IL)/DSQRT(RAD(IR))
    5    AMTC(IR)=DDENS(IR,JL)*PPOT(IR,IL)/(RAD(IR)**0.5D0)
         IF(IENUN.EQ.0) THEN
           HABB(JL,IL)=
     *       ATINT1(AMTC,RAD,N,STEP,1,STEP)/(QFF(IL)*QFF(JL))
         ELSE
           HABB(JL,IL)= EVOLTS *
     *       ATINT1(AMTC,RAD,N,STEP,1,STEP)/(QFF(IL)*QFF(JL))
         ENDIF
      ELSE
         HABB(JL,IL)=DO
      ENDIF
*
    4 CONTINUE
*
      WRITE(UNIT,102) NXT
      WRITE(3,102) NXT
      WRITE(UNIT,103) (YL(JL),JL=1,4)
      WRITE(3,103) (YL(JL),JL=1,4)
      DO 6 IL=1,4
      IF(IENUN.EQ.0) THEN
         WRITE(UNIT,104)  YL(IL),(HABB(IL,JL),JL=1,4)
         WRITE(3,104)  YL(IL),(HABB(IL,JL),JL=1,4)
      ELSE
         WRITE(UNIT,107)  YL(IL),(HABB(IL,JL),JL=1,4)
         WRITE(3,107)  YL(IL),(HABB(IL,JL),JL=1,4)
      ENDIF
      WRITE(UNIT,105)   ((QFF(IL)*QFF(JL)),JL=1,4)
    6 WRITE(3,105)   ((QFF(IL)*QFF(JL)),JL=1,4)
*
  102 FORMAT(1X,/,5X,' L-L AND <OVERLAP> MATRIX FOR <NXT>=',I2)
  103 FORMAT(4X,' :',3X,'!',5(4X,A1,4X,' !'),/,1X,64('-'))
  104 FORMAT(4X,A1,4X,'!',5(F9.6,' !'))
  107 FORMAT(4X,A1,4X,'!',5(F9.2,' !'))
  105 FORMAT(' OVERLAP !',5(F9.2,' !'),/,1X,64('-'))
      WRITE(UNIT,106)
      WRITE(3,106)
  106 FORMAT(5X,' END OF <HUBBARD> BLOCK',/)
      RETURN
      END
*
      SUBROUTINE READA(PH,ICORE,INEW,ITERX,IENUN,IEXC)
*-------------------------------------------------------------------*
*             Subroutine to READ Input Information                  *
*-------------------------------------------------------------------*
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      CHARACTER*10 XNAME,TXT2,TXT4,TXT7,TXTS,TVAL,TCOR,TFLOW
      CHARACTER*5 TXT3,TXT5,TXTA6,TXT6,TXT8,TXT9,TXTY,TXTD,TEND
      CHARACTER*30 TXT1,INTXT*8,TY(5)*1,TXTN*5,TEXCH*17
      CHARACTER*11 TXTFAL,TXTFCO,TXTFL1,TXT1A*12,TXT1R*12
      INTEGER UNIT,UNIT1
* ----------------------------------------------------------
      COMMON  / B L A T /
     *  VR(300),RH(300),A(300),RAD(300),Y(300),B(300),
     *  A0(50),B0(50),XN(50),XL(50),
     *  XJ(50),XE(50),XZ(50),XDELL(50),DA(5),DB(5),VOC(5),
     *  RN,H,ZEFF,ZN,XION,PHI,EPS,DEL,DELRV,XALPH,XLATT,RNUC,
     *  CONVR,FN,FL,FJ,E,Q,EV,EZ,EY,V0,ET3,SUMMA,ANUC,FDLL,RBAR,
     *  VBAR,H3,FK,CS,G,Q11,Q22,HA1,HA2,HA3,HA4,HA5,HA6,HA7,HA8,
     *  HA9,TCS,ITSCF,INDSCF,NPRIN,N,J,NC1,IDIRC,NXT,NDEBU
      COMMON /ADD/ AMTC(300),CORE(300),FINT(50,5),
     *             RMT,STEP,RSTART,JPR(50),JRI,UNIT
      COMMON /TEXT1/ XNAME
      DIMENSION IXN1(24),IXL1(24),XJ1(24),IXZ1(24),INERT(6),INERTS(6)
      DIMENSION IXN2(15),IXL2(15),IXZ2(15),INTXT(6),INERTD(6),TEXCH(4)
      DIMENSION TXTFAL(5),TXTFCO(5)
* ----------------------------------------------------------
      DATA DO/0.0D0/,D2/2.D0/,D1/1.0D0/,D3/3.0D0/,D8/8.0D0/
      DATA SMALL1/1.0D-07/,SMALL2/1.0D-05/,TXT1R/'RYDBERG     '/
      DATA TXTY /'YES  '/,TXTD /'DIRAC'/,TEND /'END  '/
      DATA TXTS /'STANDART  '/,TXTN/'NO   '/
      DATA TXTFAL /'outatm1.all','outatm2.all','outatm3.all',
     *             'outatm4.all','outatm5.all'/
      DATA TXTFCO /'outatm1.dat','outatm2.dat','outatm3.dat',
     *             'outatm4.dat','outatm5.dat'/
      DATA TEXCH /'Barth-Hedin','Slater X-alpha:','Barth-Hedin-Janak',
     *            'Vosko-Wilk-Nusair'/
* ---------------------------------------------------------
      DATA IXN1/1, 2,2,2, 3,3,3, 4,3,3,4,4, 5,4,4,5,5, 6,4,4,5,5,6,6/,
     *     IXL1/0, 0,1,1, 0,1,1, 0,2,2,1,1, 0,2,2,1,1, 0,3,3,2,2,1,1/,
     *     IXZ1/2, 2,2,4, 2,2,4, 2,4,6,2,4, 2,4,6,2,4, 2,6,8,4,6,2,4/,
     * XJ1/0.5D0, 2*0.5D0,1.5D0, 2*0.5D0,1.5D0, 0.5D0,1.5D0,2.5D0,
     * 0.5D0,1.5D0, 0.5D0,1.5D0,2.5D0,0.5D0,1.5D0, 0.5D0,2.5D0,3.5D0,
     * 1.5D0,2.5D0,0.5D0,1.5D0/
* ---------------------------------------------------------
      DATA IXN2/1, 2,2, 3,3, 4,3,4, 5,4,5, 6,4,5,6/,
     *     IXL2/0, 0,1, 0,1, 0,2,1, 0,2,1, 0,3,2,1/,
     *     IXZ2/2, 2,6, 2,6, 2,10,6, 2,10,6, 2,14,10,6/
* ---------------------------------------------------------
      DATA INERT/2,10,18,36,54,86/, INTXT /' HE(2) +',' NE(10)+',
     * ' AR(18)+','KR(36)+',' XE(54)+',' RN(86)+'/
      DATA INERTD/1,4,7,12,17,24/,INERTS/1,3,5,8,11,15/
      DATA TCOR/'<SEMICORE>'/, TVAL/'  <VALENT>'/
      DATA TY/'S','P','D','F','G'/
* ---------------------------------------------------------
      UNIT1=16
* ----------------------------------------------------------
  102 FORMAT(A30,A10)
  103 FORMAT(A30,F10.5)
  109 FORMAT(A30,A5)
*
      N=300
      EPS=2.5D3
      DEL=SMALL1
      RBAR=DO
      VBAR=DO
      NDEBU=0
*
* the TITLE and ATOMIC number
*
      READ(UNIT1,101) TXT1,NXT,TXT3,XNAME
  101 FORMAT(A30,I5,A3,A10)
      IF(NXT.LE.0) RETURN
*
* the ENERGY UNITs determination (EV or RY)
*
      READ(UNIT1,126) TXT1,TXT1A
  126 FORMAT(A30,A12)
      IENUN=-1
      IF(TXT1A.EQ.TXT1R) THEN
         IENUN=0
      ELSE
         IENUN=1
      ENDIF
*
      READ(UNIT1,102) TXT1,TXT2
      IF(TXT2.EQ.TXTS) THEN
         NPRIN=0
         ITSCF=10
      ELSE
         NPRIN=1
         ITSCF=0
      ENDIF
*
      READ(UNIT1,103) TXT1,RN
      IF(RN.EQ.DO) RN=60.D0
*
      READ(UNIT1,103) TXT1, H
      IF(H.EQ.DO) H=32.0
      H=D1/H
      STEP=H
      H3=H/D3
      HA1=D8*H3
      HA2=1.125D0*H3
      HA3=2.375D0*H3
      HA4=0.625D0*H3
      HA5=0.125D0*H3
      HA6=H3/240.D0
      HA7=H3/120.D0
      HA8=H3/80.D0
      HA9=H3/60.D0
*
* JRI at MT sphere determination
*
      READ(UNIT1,103) TXT1,RMT
      IF(RMT.NE.0) THEN
         ILI=DLOG(RN/RMT)/H
         XMT=DBLE(ILI)*H
         RN=RMT*DEXP(XMT)
         JRI=N-ILI
      ENDIF
*
      READ(UNIT1,103) TXT1,PHI
      IF(PHI.EQ.DO) PHI=0.3D0
      PH=PHI
*
      READ(UNIT1,108) TXT1,NC1
  108 FORMAT(A30,I10)
      IF(NC1.EQ.0) NC1=30
*
      READ(UNIT1,109) TXT1,TXT3
      IF(TXT3.EQ.TXTD) THEN
         IDIRC=1
      ELSE
         IDIRC=0
      ENDIF
*
      READ(UNIT1,102) TXT1,TXT4
*
      READ(UNIT1,103) TXT1,ANUC
      IF(ANUC.GT.DO) THEN
         RNUC=(2.08D-05)*ANUC**(D1/D3)
         R=RN*DEXP(-(N-3)*H)
         IF(RNUC.LT.R) THEN
              WRITE(UNIT,112) RNUC
              RNUC=DO
         ENDIF
      ELSE
         RNUC=DO
      ENDIF
  112 FORMAT(5X,'* WARNING : RNUC=',D11.4,' TOO SMALL FOR GRID *',
     *   /,5X,'* IN CALCULATION WILL BE USED RNUC=0.0 *',/)
*
      READ(UNIT1,109) TXT1,TXT5
      IF(TXT5.EQ.TXTY) THEN
         XLATT=D1
      ELSE
         XLATT=DO
      ENDIF
*
      READ(UNIT1,114) TXT1,XALPH,IEXC
  114 FORMAT(A30,F10.5,I5)
      IF(XALPH.EQ.DO) XALPH=0.666666667
      IEXC=IEXC+1
*
      READ(UNIT1,103) TXT1,ZN
      DELRV=SMALL1*ZN
*
      READ(UNIT1,109) TXT1,TXTA6
      IF(INDSCF.EQ.1) TXTA6=TXTN
*
      READ(UNIT1,109) TXT1,TXT6
      IF(TXT6.EQ.TXTY.OR.INDSCF.EQ.1) THEN
         INEW=1
      ELSE
         INEW=0
      ENDIF
*
      READ(UNIT1,102) TXT1,TXT7
      IF(TXT7.EQ.TXTS) THEN
*
*  the INERT CORE function determination
*
         IZN=IDNINT(ZN)
         IFLOW=0
         DO 1 I=1,6
         IF(IZN.GE.INERT(I)) IFLOW=I
    1    CONTINUE
*
*   in the HYDROGEN case
*
         IF(IFLOW.EQ.0) THEN
            IFLOW=1
            INERT(1)=0
            INTXT(1)='HYDROGEN'
            INERTD(1)=0
            INERTS(1)=0
            I1=0
         ENDIF
*
         IF(IDIRC.EQ.1) THEN
              DO 2 I=1,INERTD(IFLOW)
              XN(I)=DBLE(IXN1(I))
              XL(I)=DBLE(IXL1(I))
              XJ(I)=XJ1(I)
              XZ(I)=DBLE(IXZ1(I))
    2         JPR(I)=0
              I1=INERTD(IFLOW)
         ELSE
              DO 3 I=1,INERTS(IFLOW)
              XN(I)=DBLE(IXN2(I))
              XL(I)=DBLE(IXL2(I))
              XJ(I)=XL(I)
              XZ(I)=DBLE(IXZ2(I))
    3         JPR(I)=0
              I1=INERTS(IFLOW)
         ENDIF
      ELSE
         WRITE(UNIT,121)
         STOP
      ENDIF
  121 FORMAT(5X,'<AUTOSTOP - NOT STANDART CORE DOSN''T WORK>')
*
*   the SEMICORE reading
*
      READ(UNIT1,109) TXT1,TXT8
      IF(TXT8.NE.TXTY) THEN
         I2=I1
      ELSE
         CALL READA1(I1,I2)
      ENDIF
      ICORE=I2
*
*   the VALENT electrons' reading
*
      READ(UNIT1,109) TXT1,TXT9
      IF(TXT9.NE.TXTY) THEN
         J=I2
      ELSE
         CALL READA1(I2,I3)
      ENDIF
      J=I3
*
      DO 4 I=1,J
      XE(I)=-DBLE((J+1-I)**2)/XN(I)**2
    4 XDELL(I)=DO
*
      SUMCOR=DO
      DO 5 I=1,ICORE
    5 SUMCOR=SUMCOR+XZ(I)
*
      SUMVAL=SUMCOR
      DO 6 I=(ICORE+1),J
    6 SUMVAL=SUMVAL+XZ(I)
*
      XION=ZN-SUMVAL
      Q=SUMCOR
*      < Q > USED IN  C O R E  V A R I A N T
      IF(TXT4.EQ.TXTS.AND.INDSCF.EQ.1) XION=DO
*
      IF(INDSCF.EQ.1) J=ICORE
*
      D=DEXP(H)
      R=RN*DEXP(-H*N)
      RSTART=-H*(N-1)+DLOG(RN)
      DO 7 K=1,N
      R=R*D
      RAD(K)=R
      A(K)=DO
      B(K)=DO
    7 VR(K)=DO
*
      WRITE(UNIT,201) NXT,XNAME,IZN,TXT1A
  201 FORMAT(5X,/,3X,61('-'),/,3X,'*',15X,'ELEMENT NO',I2,
     *   2X,A10,'(',I3,')',15X,'*',/,3X,'*',15X,
     *   'USED ENERGY UNIT - ',A12,13X,'*',/,3X,'*',15X,
     *   'MAIN INPUT/OUTPUT INFORMATION -',13X,'*')
*
      IF(TXTA6.NE.TXTY.AND.ITERX.NE.0) THEN
         IUNIT=17+NXT
         IF(INDSCF.NE.1) THEN
            TXTFL1=TXTFAL(NXT)
         ELSE
            TXTFL1=TXTFCO(NXT)
         ENDIF
*
         OPEN (IUNIT,FILE=TXTFL1,FORM='UNFORMATTED')
*
         READ (IUNIT) JNAME,NNAME
         IF(JNAME.NE.J.OR.NNAME.NE.N) THEN
            CLOSE(IUNIT,STATUS='KEEP')
            WRITE(UNIT,125) JNAME,J,NNAME,N
            CALL TOMAS
            WRITE(UNIT,204)
         ELSE
            READ (IUNIT) (XE(K),K=1,JNAME)
            IF(INDSCF.NE.1) THEN
               READ (IUNIT) (RH(K),K=1,NNAME)
               WRITE(UNIT,202) NXT
            ELSE
               READ (IUNIT) (CORE(K),K=1,NNAME)
               WRITE(UNIT,203) NXT
            ENDIF
            CLOSE (IUNIT,STATUS='KEEP')
         ENDIF
      ELSE
         CALL TOMAS
         WRITE(UNIT,204)
      ENDIF
  125    FORMAT(5X,'<< WRONG DATA FOR READING IN <ATOM> >>',
     *   /,5X,'<< JNAME=',I2,' J=',I2,' NNAME=',I3,' N=',I3,' >>',/)
  202 FORMAT(3X,'*',8X,'<ATOM> START DENSITY IS READ (FILE',
     *  ' outatm',I1,'.dat)',4X,'*')
  203 FORMAT(3X,'*',8X,'<CORE> START DENSITY IS READ (FILE',
     *  ' outatm',I1,'.dat)',4X,'*')
  204 FORMAT(3X,'*',15X,'NEW START DENSITY IS GENERATED',14X,'*')
*
      IF(INEW.EQ.1) THEN
         WRITE(UNIT,205) (NXT+17),NXT
      ELSE
         IF(INDSCF.EQ.0) THEN
            WRITE(UNIT,217)
         ELSE
            WRITE(UNIT,206)
         ENDIF
      ENDIF
  205 FORMAT(3X,'*',5X,'<OUT> DENSITY WILL BE WRITEN IN FILE',I3,
     * ' outatm',I1,'.dat',3X,'*')
  217 FORMAT(3X,'*',2X,'<OUT.INF.FOR <STR>> WILL BE WRITEN IN FILE 4'
     * ,' outatm.str',2X,'*')
  206 FORMAT(3X,'*',5X,'<OUTPUT> DENSITY WILL <NOT> BE WRITE',10X,'*')
*
      WRITE(UNIT,207) NXT,RMT,NXT,JRI
  207 FORMAT(3X,'*',6X,'BAND PARAMETERS - RMT(',I2,')=',F8.5,
     *    ' JRIS(',I2,')=',I3,6X,'*')
*
      WRITE(UNIT,208) XLATT,RNUC
  208 FORMAT(3X,'*',4X,'LATTER CORRECTION -',F5.2,
     *  ' AT.NUCLEAR RADIUS ',F8.6,4X,'*')
*
      IF(IEXC.EQ.2) THEN
         WRITE(UNIT,230) TEXCH(IEXC),XALPH
      ELSE
         WRITE(UNIT,231) TEXCH(IEXC)
      ENDIF
  230 FORMAT(3X,'*',10X,'Exchange by ',A17,F10.7,10X,'*')
  231 FORMAT(3X,'*',15X,'Exchange by ',A17,15X,'*')
*
      WRITE(UNIT,209) XION
  209 FORMAT(3X,'*',20X,'IONICITY -',F5.2,24X,'*')
*
      WRITE(UNIT,210) INTXT(IFLOW)
  210 FORMAT(3X,61('-'),/,3X,'*',17X,'INNER SHELLS FROM ',A8,16X,'*')
*
      DO 8 I=(I1+1),I3
      IF(I.GT.I2) THEN
         TFLOW=TVAL
      ELSE
         TFLOW=TCOR
      ENDIF
      N1=IDNINT(XN(I))
      L1=IDNINT(XL(I))+1
      J1=IDNINT(XJ(I)*D2)
      IF(INDSCF.EQ.0.OR.I.LE.I2) THEN
         IF(XZ(I).EQ.DO) THEN
            WRITE(UNIT,211) TFLOW,N1,TY(L1),J1
         ELSE
            WRITE(UNIT,212) TFLOW,N1,TY(L1),J1,XZ(I)
         ENDIF
      ENDIF
    8 CONTINUE
  211 FORMAT(3X,'*',10X,A10,' SHELL - ',I1,A1,I1,'/2',
     *       ' - <VIRTUAL>',13X,'*')
  212 FORMAT(3X,'*',10X,A10,' SHELL - ',I1,A1,I1,'/2',
     *       ' - OCCUP.-',F5.2,10X,'*')
*
      WRITE(UNIT,216) SUMCOR,(SUMVAL-SUMCOR)
  216 FORMAT(3X,'*',59X,'*',/,3X,'*',2X,'NUMBER OF ELECTRONS ',
     * '<CORE+SEMICORE>=',F5.2,' <VALENT>=',F5.2,1X,'*')
*
      IF(IDIRC.EQ.1) THEN
         WRITE(UNIT,214)
      ELSE
         WRITE(UNIT,215)
      ENDIF
  214 FORMAT(3X,61('-'),//,9X,'<RELATIVISTIC CALCULATION START>',/)
  215 FORMAT(3X,61('-'),//,8X,'< NOT RELATIV. CALCULATION START>',/)
*
      RETURN
      END
*
      SUBROUTINE WRITA(INEW,NFLOW,IENUN,IEXC,EY1,V1,SUMCOR)
*-----------------------------------------------------------*
*                                                           *
*     This Program Writes OUTPUT Information                *
*                                                           *
*-----------------------------------------------------------*
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      CHARACTER*10 XNAME,IDENT1,IDENT2
      CHARACTER*6 ANLJ,NAME2,TXT1A*13,TEXCH*17
      CHARACTER*11 TXTFAL,TXTFCO,TXTFL1
      INTEGER UNIT
* ---------------------------------------------------------
      COMMON  / B L A T /
     *  VR(300),RH(300),A(300),RAD(300),Y(300),B(300),
     *  A0(50),B0(50),XN(50),XL(50),
     *  XJ(50),XE(50),XZ(50),XDELL(50),DA(5),DB(5),VOC(5),
     *  RN,H,ZEFF,ZN,XION,PHI,EPS,DEL,DELRV,XALPH,XLATT,RNUC,
     *  CONVR,FN,FL,FJ,E,Q,EV,EZ,EY,V0,ET3,SUMMA,ANUC,FDLL,RBAR,
     *  VBAR,H3,FK,CS,G,Q11,Q22,HA1,HA2,HA3,HA4,HA5,HA6,HA7,HA8,
     *  HA9,TCS,ITSCF,INDSCF,NPRIN,N,J,NC1,IDIRC,NXT,NDEBU
      COMMON /ADD/ AMTC(300),CORE(300),FINT(50,5),
     *             RMT,STEP,RSTART,JPR(50),JRI,UNIT
      COMMON /TEXT1/ XNAME
      COMMON /TIP/ CORDNS1(50),EPSCOR1(50),TOP(300),EPEDE1,EPEDE2,
     * V_old(300),V_real(300),ICORE
      DIMENSION NIST(50),ANLJ(80),TXT1A(2),TEXCH(4)
      DIMENSION TXTFAL(5),TXTFCO(5)
* --------------------------------------------------------
      DATA ANLJ/' 1S1/2',' 2S1/2',' 2P1/2',' 2P3/2',
     * ' 3S1/2',' 3P1/2',' 3P3/2',' 3D3/2',' 3D5/2',' 4S1/2',
     * ' 4P1/2',' 4P3/2',' 4D3/2',' 4D5/2',' 4F5/2',' 4F7/2',
     * ' 5S1/2',' 5P1/2',' 5P3/2',' 5D3/2',' 5D5/2',' 5F5/2',
     * ' 5F7/2',' 5G7/2',' 5G9/2',' 6S1/2',' 6P1/2',' 6P3/2',
     * ' 6D3/2',' 6D5/2',' 6F5/2',' 6F7/2',' 6G7/2',' 6G9/2',
     * ' 6H9/2','6H11/2',' 7S1/2',' 7P1/2',' 7P3/2',' 7D3/2',
     * ' 7D5/2',' 7F5/2',' 7F7/2',' 7G7/2',' 7G9/2',' 7H9/2',
     * '7H11/2','7I11/2','7I13/2',' 8S1/2',
     * '1S','2S','2P','3S','3P','3D','4S','4P','4D','4F',
     * '5S','5P','5D','5F','5G','6S','6P','6D','6F','6G','6H',
     * '7S','7P','7D','7F','7G','7H','7I','8S','8P'/
* ----------------------------------------------------------
      DATA P4/0.0795774715D0/,D2/2.D0/,D1/1.D0/,EVOLTS/13.6058/,
     * D3/3.D0/,DO/0.0D0/,D12/0.5D0/,DZER/1.D-2/,PRO2/2.25D0/
      DATA IDENT1/'CORE'/,IDENT2/'AMTC'/,NAME2/' WRITE'/
      DATA TXT1A/'RYDBERGS','ELECTRONVOLTS'/
      DATA TEXCH /'Barth-Hedin','Slater X-alpha:','Barth-Hedin-Janak',
     *            'Vosko-Wilk-Nusair'/
      DATA TXTFAL /'outatm1.all','outatm2.all','outatm3.all',
     *             'outatm4.all','outatm5.all'/
      DATA TXTFCO /'outatm1.dat','outatm2.dat','outatm3.dat',
     *             'outatm4.dat','outatm5.dat'/
* ---------------------------------------------------------
      PROBA=DLOG(RH(2)/RH(1))/H-D2
      IF(DABS(PROBA).LT.DZER) PROBA=DO
      S3=D3+PRO2*PROBA
*
      WRITE (  3,101) XNAME,TXT1A(IENUN+1),TEXCH(IEXC)
      WRITE (UNIT,101) XNAME,TXT1A(IENUN+1),TEXCH(IEXC)
  101 FORMAT(2X,60('-'),//,25X,A10,//,2X,60('-'),/,
     * 12X,'THE USED ENERGY UNITS ARE ',A13,/,12X,
     * 'EXCHANGE and CORRELATION by ',A17,/)
*
      IF(INEW.EQ.1) THEN
         IUNIT=17+NXT
         IF(INDSCF.NE.1) THEN
            TXTFL1=TXTFAL(NXT)
         ELSE
            TXTFL1=TXTFCO(NXT)
         ENDIF
*
         OPEN (IUNIT,FILE=TXTFL1,FORM='UNFORMATTED')
*
         WRITE (IUNIT) J,NFLOW
         WRITE (IUNIT) (XE(K),K=1,J)
         IF(INDSCF.NE.1) THEN
            WRITE (IUNIT) (RH(K),K=1,NFLOW)
         ELSE
            WRITE (IUNIT) (CORE(K),K=1,NFLOW)
         ENDIF
         CLOSE(IUNIT,STATUS='KEEP')
      ELSE
         IF(INDSCF.NE.1) THEN
            WRITE(4) RSTART,STEP,RMT,JRI,N,ICORE
            WRITE(4) (CORE(I),I=1,N)
            WRITE(4) (AMTC(I),I=1,N)
            WRITE(4) (RAD(I),I=1,N)
*
            WRITE(4) (CORDNS1(IC),IC=1,ICORE)
            WRITE(4) (EPSCOR1(IC),IC=1,ICORE)
	 ENDIF
      ENDIF
*
*     Calculate Energy of Valence Electrons
*
      DO 8 I=1,N
    8 A(I)=AMTC(I)-CORE(I)
*
      QTOTV=ATINT(A,RAD,N,H,3,S3)
      QVAL=ATINT(A,RAD,JRI,H,3,S3)
      QOUT=QTOTV-QVAL
*
      CALL POISAT(A,N,B,QTOTV,RSTART,H)
*
      DO 9 I=1,N
      B(I)=B(I)*(RAD(I))**0.5D0*D12
      FLOW=D12*(A(I)/(RAD(I)*RAD(I)))
*
      CALL XCPOT(FLOW,FLOW,(FLOW+FLOW),VXC1,VXC2,EXC)
*
      Y(I)=(D12*RAD(I))*EXC
    9 CONTINUE
*
      EPOTV=ATINT(A,B,JRI,H,3,S3)*D12
*     EXVAL=ATINT(A,Y,N,H,3,S3)
      EXVAL=ATINT(A,Y,JRI,H,3,S3)
      EKINV=ATINT(A,VR,JRI,H,3,S3)
      SUMVAL=(SUMMA-SUMCOR)
      EKINVT=SUMVAL-EKINV
*
      ETOTVA=EKINVT+EPOTV+EXVAL
*
      IOPEN=0
      DO 1 I=1,J
      IF(JPR(I).EQ.1) IOPEN=1
      IF(IENUN.EQ.0) THEN
         XE(I)=XE(I)*D2
         XDELL(I)=XDELL(I)*D2
      ELSE
         XE(I)=XE(I)*D2*EVOLTS
         XDELL(I)=XDELL(I)*D2*EVOLTS
      ENDIF
    1 CONTINUE
*
      WRITE (UNIT,105)
      IF(ITSCF.EQ.10.AND.IOPEN.EQ.1) WRITE (3,105)
  105 FORMAT(1X,/,2X,'ORBITAL',4X,'N',5X,'L',5X,'J',4X,
     * 'ELECTRONS',5X,'EIGENVALUE',9X,'DELTA OF EIGENV.'/)
*
      DO 2 I=1,J
      IF(IDIRC.EQ.0) THEN
         FLOW=D12*XN(I)*(XN(I)-D1)+XL(I)+51.001D0
      ELSE
         FLOW=XN(I)*(XN(I)-D2)+XL(I)+XJ(I)+1.501D0
      ENDIF
      NI=IDNINT(FLOW)
      NIST(I)=NI
      IF(ITSCF.EQ.10.AND.JPR(I).EQ.1) THEN
        WRITE (UNIT,106)
     *  ANLJ(NI),XN(I),XL(I),XJ(I),XZ(I),XE(I),XDELL(I)
        WRITE (3,106)
     *  ANLJ(NI),XN(I),XL(I),XJ(I),XZ(I),XE(I),XDELL(I)
      ENDIF
  106 FORMAT(3X,A6,3F6.1,F8.2,3X,G18.5,2X,D16.3)
    2 CONTINUE
*
      QTOTAL=ATINT(RH,RAD,N,H,3,S3)
      QTOTA1=ATINT(RH,RAD,JRI,H,3,S3)
      QOUT=QTOTAL-QTOTA1
*
      EV=ATINT(RH,VR,N,H,3,S3)
      EV1=ATINT(RH,VR,JRI,H,3,S3)
      EPEDE1=ATINT(RH,TOP,N,H,3,S3)
      EPEDE2=ATINT(RH,TOP,JRI,H,3,S3)
*
* ----- OK: ERROR FOUND - DIVIDE BY ZERO, RNUC < MACHINE EPSILON -----
* ----> D0 is not initialized! USE DO !!!
*      IF(RNUC.GT.D0) THEN
      IF(RNUC.GT.DO) THEN
         NC=1
         DO 11 I=1,N
         IF(RAD(I).LT.RNUC) NC=I
   11    CONTINUE
         DO 12 I=1,NC
         X=RAD(I)/RNUC
   12    Y(I)=ZN*X*(1.5D0-0.5D0*X**2)
         DO 10 I=(NC+1),N
   10    Y(I)=ZN
         EZ1=-ATINT(RH,Y,N,H,3,S3)
         EZ=D12*(EZ+EZ1)
      ENDIF
*
      EKIN=SUMMA-EV
*
      ET1=SUMMA-D12*EPEDE1+(V0-EY)
      ET3=SUMMA-D12*(EV-EZ-EY)+(V0-EY)
      ET4=SUMMA-D12*(EV1-EZ-EY1)+(V1-EY1)
*
      IF(INDSCF.NE.1) THEN
        COH1=D2*(ET1-ET4)
*
        WRITE(4) (D2*SUMMA),(D2*SUMCOR),(D2*SUMVAL)
        WRITE(4) (D2*ET1),(D2*EKINVT),(D2*ETOTVA)
        WRITE(4) (D2*EPOTV),(D2*EXVAL),QOUT
*
*   transition from ATOMIC ENERGY units to RYDBERGS
*
         IF(IENUN.EQ.0) THEN
            EVOLT=D1
         ELSE
            EVOLT=EVOLTS
         ENDIF
*
*   transition to ELECTRONVOLTS
*
         SUMMA = D2 * SUMMA * EVOLT
         ET1   = D2 * ET1 * EVOLT
         ET3   = D2 * ET3 * EVOLT
         ET4   = D2 * ET4 * EVOLT
         EKIN  = D2 * EKIN * EVOLT
         EZ    = D2 * EZ * EVOLT
         EV    = D2 * (EV-EY) * EVOLT
         EV1   = D2 * (EV1-EY1) * EVOLT
*        V0    = D2 * EY * 0.75D0 * EVOLT
         V0    = D2 * V0 * EVOLT
         V1    = D2 * V1 * EVOLT
*
         EPOTV = D2 * EPOTV * EVOLT
         EXVAL = D2 * EXVAL * EVOLT
         SUMVAL= D2 * SUMVAL* EVOLT
         EKINV = D2 * EKINV * EVOLT
         EKINVT= D2 * EKINVT* EVOLT
         ETOTVA= D2 * ETOTVA* EVOLT
*
         ELEL=(EV-EZ)*D12
         EPOTAL=ET3-EKIN
         RATIO=(-EPOTAL)/EKIN
         WRITE (UNIT,107) ZN,QTOTAL,SUMMA,EKIN,V0,EV,EZ,ELEL,
     *                    EPOTAL,RATIO,ET1
         WRITE (3,107) ZN,QTOTAL,SUMMA,EKIN,V0,EV,EZ,ELEL,
     *                    EPOTAL,RATIO,ET1
*
         WRITE (3,115) SUMVAL,EKINV,EKINVT,EPOTV,EXVAL,ETOTVA
         WRITE (UNIT,115) SUMVAL,EKINV,EKINVT,EPOTV,EXVAL,ETOTVA
      ENDIF
  115 FORMAT(8X,'  INFORMATION ON VALENCE ELECTRON ENERGY',/,
     * 8X,30('-'),/,15X,'SUM Val.El.EigenVal=',F14.4,/,5X,
     * '(Tot.Pot.*Val.Dens) Kin.Contr=',F14.4,/,
     * 17X,'Val.El.KIN.Energy=',F14.4,
     * /,17X,'Val.El.POT.Energy=',F14.4,
     * /,17X,'Val.El.EXC.Energy=',F14.4,/,8X,30('-'),/,
     * 22X,'TOTAL ENERGY=',F14.4,///)
  107 FORMAT (1X,/,19X,' NUCLEAR CHARGE=',F13.6,/,
     * 8X,'INTEGRAL OF CHARGE DENSITY=',F13.6,/,8X,30('-'),
     * /,5X,'SUM OF THE ENERGY EIGENVALUES=',F14.4,
     * /,20X,'KINETIC ENERGY=',F14.4,
     * /,19X,'EXCHANGE ENERGY=',F14.4,
     * /,10X,'COULOMB POTENTIAL ENERGY=',F14.4,/,8X,30('-'),
     * /,7X,'NUCLEAR-EL.POTENTIAL ENERGY=',F14.4,
     * /,4X,'EL-EL COULOMB POTENTIAL ENERGY=',F14.4,
     * /,18X,'POTENTIAL ENERGY=',F14.4,
     * /,12X,'VIRIAL RATIO (Pot/Kin)=',F14.4,
     * /,8X,30('-'),/,22X,'TOTAL ENERGY=',F14.4,//)
*
      WRITE (UNIT,112)
      IF(ITSCF.EQ.10.AND.IOPEN.EQ.1) WRITE (3,112)
  112 FORMAT(2X,/,2X,'ORBITAL',6X,'<R**-3>',7X,'<R**-1>',9X,
     *  '<R>',9X,'<R**2>',8X,'<R**4>',/)
*
      DO 7 K=1,J
      NI=NIST(K)
      IF(JPR(K).GE.1) THEN
         WRITE (UNIT,114) ANLJ(NI),(FINT(K,IM),IM=1,5)
         IF(ITSCF.EQ.10)
     *     WRITE (3,114) ANLJ(NI),(FINT(K,IM),IM=1,5)
  114    FORMAT(1X,A8,5(3X,F11.5))
      ENDIF
    7 CONTINUE
*
      RETURN
      END

