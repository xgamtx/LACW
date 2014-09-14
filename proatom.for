*-------------------------------------------------------------------*
*              A T O M     P A C K                                  *
*     This is Atomic main Programm for Self-Consistent-Field        *
*                  Calculations.                                    *
*---------------------------------------------------------*---------*
*       L A P W    M E T H O D                            *
*       ----------------------                            *
*       GESPACK  VERSION -  3.10.91                       *
*       Source - ELLIS's programm for DV Method           *
*       ---------------------------------------           *
*       It was rewritten by NIKOLAEV A.V.(Moscow,USSR)    *
*       --------------------------------------------      *
*       TERMINAL UNIT = 1                                 *
*---------------------------------------------------------*
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
      COMMON /ADD/ AMTC(300),CORE(300),FINT(50,5),
     *             RMT,STEP,RSTART,JPR(50),JRI,UNIT
      COMMON /TEXT1/ XNAME
      COMMON /HUBB/ DDENS(300,4)
      COMMON /TIP/ CORDNS1(50),EPSCOR1(50),TOP(300),EPEDE1,EPEDE2,
     * V_old(300),V_real(300),ICORE
*
      COMMON /CXCPOT/ XCCP,XCCF,XCRP,XCRF,XALPHA,OTH,FTH,AA,BB,IXCH1
      COMMON /EXCH/ BLX,BLXM,ALPHA,ICORR,ICHNG1,ICHNG2,ICHNG3,ICHNG4
*
      DIMENSION  VSTOR(300),USTOR(300)
* -----------------------------------------------------------
      DATA D12/0.5D0/,DO/0.0D0/,D2/2.0D0/,D1/1.0D0/,
     *  POW/2.1D0/,IOPEN/0/,CONTOP/0.004D0/,SMALL/1.D-6/,
     *  PI/3.1415926535897932384D0/,D3/3.D0/,D4/4.D0/
*
*               < A T O M I C     V A R I A N T >
*
      ITERX=1
      INDSCF=0
      UNIT=6
*
*     Terminal Unit is defined by <UNIT>
*
      OPEN ( 4,FILE='outatm.str',FORM='UNFORMATTED')
      OPEN (16,FILE='inatom',FORM='FORMATTED')
      OPEN ( 3,FILE='outatm.see',FORM='FORMATTED')
      VOC(2)=DO
*
    1 LAST=0
*
      CALL READA(PH,ICORE,INEW,ITERX,IENUN,IEXC)
*
      ALPHA=XALPH
      ICORR=IEXC-1
*
      CALL SXCPOT
*
      IF(NXT.LE.0) THEN
         CLOSE(16,STATUS='KEEP')
         CLOSE( 3,STATUS='KEEP')
	 IF(INEW.NE.1) THEN
            CLOSE( 4,STATUS='KEEP')
         ELSE
            CLOSE( 4,STATUS='DELETE')
         ENDIF
         STOP
      ENDIF
*
      DO 20 L1=1,4
      DO 20 K=1,N
   20 DDENS(K,L1)=DO
*
      CONVR=ZN*D2
      QTOTAL=ZN-XION
      IOPEN=0
      PSI=D1-PHI
*
*     ------ Start of New Iteration -----
*
      DO 2 NCY = 0 , (NC1+1)
*
*     ------ Density Mixing --------------
*
      IF(CONVR.GT.CONTOP.AND.IOPEN.EQ.0) THEN
        IF(NCY.GT.0) THEN
          DO 3 K=1,N
          FLOW=PSI*AMTC(K)+PHI*RH(K)
          RH(K)=FLOW
    3     AMTC(K)=FLOW
        ELSE
          DO 4 K=1,N
    4     AMTC(K)=RH(K)
        ENDIF
*
        DO 5 K=1,N
        AMTCK=D12*(RH(K)/(RAD(K)*RAD(K)))
*
        CALL XCPOT(AMTCK,AMTCK,(AMTCK+AMTCK),VXC1,VXC2,EXC)
*
        CORE(K)=(D12*RAD(K))*EXC
    5   Y(K)=(D12*RAD(K))*VXC1
*
        CALL POISAT(RH,N,A,QTOTAL,RSTART,H)
        CALL POTAT(NCY,UNIT,NERRO,N)
*
      ELSE
*
*        ----- Potential Function Mixing ------
*
        IF(IOPEN.EQ.0) WRITE(UNIT,102)
  102   FORMAT(5X,/,60('-'),
     *     /,5X,'----> PRATT''S SCHEME IS TURNING ON <-----')
*
        CALL POISAT(RH,N,A,QTOTAL,RSTART,H)
*
        NTEST=MOD(NCY+1,2)
*
        DO 6 K=1,N
        AMTCK=D12*(RH(K)/(RAD(K)*RAD(K)))
*
        CALL XCPOT(AMTCK,AMTCK,(AMTCK+AMTCK),VXC1,VXC2,EXC)
*
        CORE(K)=(D12*RAD(K))*EXC
    6   Y(K)=D12*RAD(K)*VXC1
*
        IF(NTEST.EQ.0.AND.IOPEN.EQ.1) THEN
*
*      ------ PRATT Average Scheme ------
*
          DO 7 K=1,N
    7     B(K)=VR(K)
          CONVR1=CONVR/POW
*
          CALL POTAT(NCY,UNIT,NERRO,N)
*
          IF(CONVR.LT.CONVR1) THEN
              PHF=PH
              AFLOW=PHI
          ELSE
              PHF=PHI
              AFLOW=CONTOP
          ENDIF
          PSF=D1-PHF
          DO 8 K=1,N
          DETERM=USTOR(K)+VR(K)-B(K)-VSTOR(K)
          IF(DABS(DETERM).LT.SMALL) THEN
              ALPH=PSF
          ELSE
              ALPH=((USTOR(K)*VR(K)-B(K)*VSTOR(K))/DETERM-
     *             VR(K))/(B(K)-VR(K))
*             ALPH=DMAX1(ALPH,AFLOW)
              IF(ALPH.LT.AFLOW)   ALPH=AFLOW
*             ALPH=DMIN1(ALPH,PSF)
              IF(ALPH.GT.PSF)     ALPH=PSF
          ENDIF
          VR(K)=ALPH*B(K)+(D1-ALPH)*VR(K)
    8     CONTINUE
        ELSE
*
*   ------ Arithmetic Average Scheme ------
*
          IOPEN=1
          DO 9 K=1,N
    9     USTOR(K)=VR(K)
          CONVR1=CONVR/POW
*
          CALL POTAT(NCY,UNIT,NERRO,N)
*
          IF(CONVR.LT.CONVR1) THEN
              PHF=PH
          ELSE
              PHF=PHI
          ENDIF
          PSF=D1-PHF
          DO 10 K=1,N
          VRFLOW=VR(K)
          VSTOR(K)=VRFLOW
   10     VR(K)=PSF*USTOR(K)+PHF*VRFLOW
        ENDIF
      ENDIF
      IF(NCY.GT.0) THEN
         WRITE (UNIT,105) NCY
         WRITE (UNIT,103) CONVR,NERRO
      ENDIF
  105 FORMAT(/,2X,'***************************************',/,
     * 5X,'* * * Atomic S-C-F Cycle Number = ',I2,'   * * *')
  103 FORMAT(' Maximum Error in R.V(R) IS',1PD11.3,
     * ' at the',I4,'th Point')
*
      IF(CONVR.LE.DELRV.OR.LAST.EQ.1) THEN
         LAST=1
         EY=ATINT(RH,Y,N,H,2,XP)
         EY1=ATINT(RH,Y,JRI,H,2,XP)
*        EZ=-ZN*A(1)/DSQRT(RAD(1))*D12
         V0=ATINT(RH,CORE,N,H,2,XP)
         V1=ATINT(RH,CORE,JRI,H,2,XP)
         IF(RNUC.EQ.0) THEN
            EZ=-ZN*A(1)/(RAD(1)**0.5D0)*D12
         ELSE
            DO 21 KNUC=1,N
            IF(RAD(KNUC).LT.RNUC) NC=KNUC
   21       CONTINUE
*
            DO 22 KNUC=1,NC
            RH(KNUC)=RAD(KNUC)**2
   22       Y(KNUC)=RAD(KNUC)**2.5D0
            VREN=ATINT(RH,RAD,NC,H,2,XP)
            EZ=-ATINT(A,Y,NC,H,2,XP)*ZN*D12/VREN
         ENDIF
      ENDIF
*
*      ----- End of Mixing -------
*
      SUMMA=DO
      SUMCOR=DO
*
      DO 11 K=1,N
          CORE(K)=DO
          RH(K)=DO
   11 CONTINUE
      DO 12 I=1,J
      FN=XN(I)
      FL=XL(I)
      FJ=XJ(I)
      E=XE(I)
      FDLL=XDELL(I)
      FZ=XZ(I)
*
      IF(IDIRC.NE.0) THEN
*
         CALL DIRAT(UNIT,XNORM,FZ)
*
         DO 14 K=1,N
   14    RH(K)=RH(K)+XZ(I)*(A(K)*A(K)+B(K)*B(K))
         IF(LAST.EQ.1) THEN
            IF(I.LE.ICORE) THEN
               DO 15 K=1,N
   15          CORE(K)=CORE(K)+XZ(I)*(A(K)*A(K)+B(K)*B(K))
            ELSE IF(XZ(I).GT.D12) THEN
               L1=IDNINT(FL)+1
               DO 18 K=1,N
   18          DDENS(K,L1)=XZ(I)*(A(K)*A(K)+B(K)*B(K))
            ENDIF
         ENDIF
      ELSE
*
         CALL HSDAT(UNIT,XNORM)
*
         DO 16 K=1,N
   16    RH(K)=RH(K)+XZ(I)*A(K)*A(K)
         IF(LAST.EQ.1) THEN
            IF(I.LE.ICORE) THEN
               DO 17 K=1,N
   17          CORE(K)=CORE(K)+XZ(I)*A(K)*A(K)
*
               DO 19 K=1,N
   19          DDENS(K,L1)=XZ(I)*A(K)*A(K)
            ENDIF
         ENDIF
      ENDIF
      IF(I.LE.ICORE) THEN
         SUMCOR=SUMCOR+XZ(I)*E
      ENDIF
      SUMMA=SUMMA+XZ(I)*E
      XE(I)=E
      XDELL(I)=FDLL
*
      IF(LAST.EQ.1) THEN
         IF(JPR(I).GE.1) CALL AVER(XNORM,I,H)
*
         IF(I.LE.ICORE) THEN
            EPSCOR1(I)=-D2*E
	    IF(IDIRC.NE.0) THEN
               CORDNS1(I)=XZ(I)*(A(JRI)*A(JRI)+B(JRI)*B(JRI))
	    ELSE
               CORDNS1(I)=XZ(I)*A(JRI)*A(JRI)
            ENDIF
	 ENDIF
      ENDIF
*
   12 CONTINUE
      IF(LAST.EQ.0) THEN
         IF(NCY.EQ.NC1) LAST=1
      ELSE
*
         CALL WRITA(INEW,N,IENUN,IEXC,EY1,V1,SUMCOR)
         CALL HABAT(NXT,N,RAD,IENUN)
*
         GO TO 1
      ENDIF
    2 CONTINUE
      END
*
      SUBROUTINE SXCPOT
*-----------------------------------------------------------*
*                                                           *
*        Set constants for <XCPOT>                          *
*                                                           *
*-----------------------------------------------------------*
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON /CXCPOT/ XCCP,XCCF,XCRP,XCRF,XALPHA,OTH,FTH,AA,BB,IXCH1
*
      COMMON /EXCH/ BLX,BLXM,ALPHA,ICORR,ICHNG1,ICHNG2,ICHNG3,ICHNG4
*
      IXCH1=ICORR+1
*
      GO TO ( 1, 2, 3, 4) , IXCH1
*
*     BARTH-HEDIN J. PHYS. C5,1629 (1972)
*
    1 XCCP=0.0504D0
      XCCF=0.0254D0
      XCRP=30.D0
      XCRF=75.D0
      OTH=1.D0/3.D0
      FTH=4.D0/3.D0
      AA=1.D0/2.D0**OTH
      BB=1.D0-AA
      RETURN
*
*     SLATER X-ALPHA
*
    2 OTH=1.D0/3.D0
      XALPHA=6.D0*ALPHA*(3.D0/16.D0/3.141592654**2)**OTH
      RETURN
*
*     BARTH-HEDIN-JANAK PHYS. REV. B12,1257 (1975)
*
    3 XCCP=0.045D0
      XCCF=0.0225D0
      XCRP=21.D0
      XCRF=53.D0
      OTH=1.D0/3.D0
      FTH=4.D0/3.D0
      AA=1.D0/2.D0**OTH
      BB=1.D0-AA
      RETURN
*
*     VOSKO-WILK-NUSAIR CAN. J. PHYS. 58,1200 (1980)
*
    4 OTH=1.D0/3.D0
      FTH=4.D0/3.D0
      AA=2.D0**FTH-2.D0
      RETURN
      END
*
      SUBROUTINE XCPOT(RHO1,RHO2,RHO,V1,V2,EXC)
*-------------------------------------------------------------------*
*                                                                   *
*        Calculate Exchange Correlation Potential                   *
*                                                                   *
*-------------------------------------------------------------------*
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
*
      COMMON /CXCPOT/ XCCP,XCCF,XCRP,XCRF,XALPHA,OTH,FTH,AA,BB,IXCH1
*
      DATA AP,XPO,BP,CP,QP,CP1,CP2,CP3
     *  /0.0621814D0,-0.10498D0,3.72744D0,12.9352D0,6.1519908D0,
     *   1.2117833D0, 1.1435257D0,-0.031167608D0/
      DATA AF,XFO,BF,CF,QF,CF1,CF2,CF3
     *  /0.0310907D0,-0.32500D0,7.060428D0,18.0578D0,4.7309269D0,
     *   2.9847935D0,2.7100059D0,-0.1446006D0/
      DATA DMIN /1.0D-11/,DO/0.0D0/
*
      IF(RHO.LT.DMIN) THEN
         V1=0.D0
         V2=0.D0
         EXC=0.D0
         RETURN
      ENDIF
      IF(RHO1.LT.DO.OR.RHO2.LT.DO) THEN
         WRITE(1,100) RHO1,RHO2,RHO
         WRITE(3,100) RHO1,RHO2,RHO
      ENDIF
  100 FORMAT(5X,'* AUTOSTOP /XCPOT/ Density < 0',3F10.5)
      IF(DABS((RHO1+RHO2-RHO)/RHO).GT.DMIN) THEN
         WRITE(1,101) RHO1,RHO2,RHO
         WRITE(3,101) RHO1,RHO2,RHO
         STOP
      ENDIF
  101 FORMAT(5X,'* AUTOSTOP /XCPOT/',F10.6,'+',F10.6,'.NE.',F10.6)
*
      RS1=(RHO/3.D0)**OTH
      RS=1.D0/RS1
*
      GO TO (21,22,21,23),IXCH1
*
*     BARTH-HEDIN EXCHANGE CORRELATION
*     J. PHYS. C5,1629(1972)
*
   21 RSF=RS/XCRF
      RSF2=RSF*RSF
      RSF3=RSF2*RSF
      RSP=RS/XCRP
      RSP2=RSP*RSP
      RSP3=RSP2*RSP
      FCF=(1.D0+RSF3)*DLOG(1.D0+1.D0/RSF)+0.5D0*RSF-RSF2-OTH
      FCP=(1.D0+RSP3)*DLOG(1.D0+1.D0/RSP)+0.5D0*RSP-RSP2-OTH
      EPSCP=-XCCP*FCP
      EPSCF=-XCCF*FCF
      EPSXP=-0.91633059D0/RS
      CNY=5.1297628D0*(EPSCF-EPSCP)
      X=RHO1/RHO
      FX=(X**FTH+(1.D0-X)**FTH-AA)/BB
      EXC=EPSXP+EPSCP+FX*(CNY+FTH*EPSXP)/5.1297628D0
      ARS=-1.22177412D0/RS+CNY
      BRS=-XCCP*DLOG(1.D0+XCRP/RS)-CNY
      TRX1=(2.D0*X)**OTH
      V1=ARS*TRX1+BRS
      TRX2=(2.D0*RHO2/RHO)**OTH
      V2=ARS*TRX2+BRS
      RETURN
*
*     SLATER EXCHANGE POTENTIAL
*
   22 EXC=-0.75D0*XALPHA*(0.5D0*RHO)**OTH
      V1=-XALPHA*(RHO1)**OTH
      V2=-XALPHA*(RHO2)**OTH
      RETURN
*
*     VOSKO-WILK-NUSAIR EXCHANGE CORRELATION
*     CAN. J. PHYS. 58,1200(1980)
*
   23 X=DSQRT(RS)
      XPX=X*X+BP*X+CP
      XFX=X*X+BF*X+CF
      S=(RHO2-RHO1)/RHO
      SP=1.D0+S
      SM=1.D0-S
      S4=S**4-1.D0
      FS=(SP**FTH+SM**FTH-2.D0)/AA
      BETA=1.D0/(2.74208D0+3.182D0*X+0.09873D0*X*X+0.18268D0*X**3)
      DFS=FTH*(SP**OTH-SM**OTH)/AA
      DBETA=-(0.27402D0*X+0.09873D0+1.591D0/X)*BETA**2
      ATNP=DTAN(QP/(2.D0*X+BP))
      ATNF=DTAN(QF/(2.D0*X+BF))
      ECP=AP*(DLOG(X*X/XPX)+CP1*ATNP-CP3*(DLOG((X-XPO)**2/XPX)+
     *    CP2*ATNP))
      ECF=AF*(DLOG(X*X/XFX)+CP1*ATNF-CF3*(DLOG((X-XFO)**2/XFX)+
     *    CF2*ATNF))
      EC=ECP+FS*(ECF-ECP)*(1.D0+S4*BETA)
      TP1=(X*X+BP*X)/XPX
      TF1=(X*X+BF*X)/XFX
      UCP=ECP-AP/3.D0*(1.D0-TP1-CP3*(X/(X-XPO)-TP1-XPO*X/XPX))
      UCF=ECF-AF/3.D0*(1.D0-TF1-CF3*(X/(X-XFO)-TF1-XFO*X/XFX))
      UCO=UCP+(UCF-UCP)*FS
      UC2O=UCO+(ECF-ECP)*SM*DFS
      UC1O=UCO-(ECF-ECP)*SP*DFS
      DUC=(UCF-UCP)*BETA*S4*FS+(ECF-ECP)*(-RS/3.D0)*DBETA*S4*FS
      DUC2=DUC+(ECF-ECP)*BETA*SM*(4.D0*S**3*FS+S4*DFS)
      DUC1=DUC-(ECF-ECP)*BETA*SP*(4.D0*S**3*FS+S4*DFS)
      UC1=UC1O+DUC1
      UC2=UC2O+DUC2
      EPX=-0.91633059D0/RS*(1.D0+FTH*FS/5.1297628D0)
      AMYX2=-1.22177412D0/RS*SP**OTH
      AMYX1=-1.22177412D0/RS*SM**OTH
      EXC=EC+EPX
      V1=UC1+AMYX1
      V2=UC2+AMYX2
      RETURN
      END

