*
*     Old part begins here.
*
      SUBROUTINE BASIS1(E,K,IT)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)

	include 'atomdata.fh'

	COMMON /BAS/  P(LCUT,NCUT),PR(LCUT,NCUT),
     *	PE(LCUT,NCUT),PER(LCUT,NCUT),XE(LCUT,NCUT),
     *	QN(LCUT,NCUT),ELMIN(5,NCUT),ELMAX(5,NCUT),
     *	DELTA(5,NCUT),ITS(5,NCUT),NORLW,NOWAV,ISPIN,ISPIN1,ISPIN2
	COMMON / OU / P1(NGRID),P2(NGRID),Q1(NGRID),Q2(NGRID),
     *	BG(NGRID),BGD(NGRID),HMD(NGRID),RA(NGRID),IANS
*
	COMMON /ROS/ CA(NCUT),CC1(NCUT),CC2(NCUT),
     *		SKX,SKY,SKZ,PI,AREV1,N1,NL,NBAS
	COMMON /READER/ VMAG(3),BSTART,CONVRG,DESTOR,IEVS,IMAG,
     *	NBASM,ITRMAX,KFIX,MQ,ITURN,IRESTB,MOVSTR,MOV,MOVI,LPOINT
*
	C=274.0746d0
	CIN=1.d0/(C*C)
*
      CALL DIRAC(E,K,IT,RG)
      CALL DIRE (E,K,IT,RG)
*
      CALL XSUM (RP1, 1,1,1,IT)
      CALL XSUM (RPP, 1,2,1,IT)
      CALL XSUM (RP2, 2,2,1,IT)
*
      IF(NORLW.EQ.0) THEN
         CALL XSUM (RQ1, 1,1,2,IT)
         CALL XSUM (RQQ, 1,2,2,IT)
         CALL XSUM (RQ2, 2,2,2,IT)
      ELSE
         RQ1=0.d0
         RQQ=0.d0
      ENDIF
*
      RNORM=DSQRT(RP1+RQ1*CIN)
*
      XM=1.d0+CIN*(E-BGX(JRIS(IT),IT))
      RKOE=(-RPP-RQQ*CIN)/(RP1+RQ1*CIN)

	IF (K.EQ.1) THEN
	  DO I = 1, JRIS(IT)
		WRITE(121, *) P1(I) / (RNORM)
	    WRITE(122, *) (P2(I)+P1(I)*RKOE)/(RNORM)
		WRITE(123, *) XM*Q1(I)/(RNORM)
   	    WRITE(124, *) (XM*(Q1(I)*RKOE+Q2(I))+Q1(I)*CIN)/(RNORM)
	  END DO
	END IF
 
*
      PE(K,IT)=(P2(JRIS(IT))+P1(JRIS(IT))*RKOE)/(RNORM*RMTS(IT))
*
      PER(K,IT)=(XM*(Q1(JRIS(IT))*RKOE+Q2(JRIS(IT)))+
     .                Q1(JRIS(IT))*CIN)/(RNORM*RMTS(IT))
*
      P(K,IT)=P1(JRIS(IT))/(RNORM*RMTS(IT))
*
      PR(K,IT)=XM*Q1(JRIS(IT))/(RNORM*RMTS(IT))
*
      XE(K,IT)=(RP2+2.D0*RKOE*RPP+RKOE*RKOE*RP1+
     .     CIN*(RQ2+2.D0*RKOE*RQQ+RKOE*RKOE*RQ1))/RNORM/RNORM
      RETURN
      END

      SUBROUTINE DIRAC(E,K,IT,RG)
*
*     D.D.Koelling & B.N.Harmon equations
*     Journal of Physics C v.10 N.16 1977
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
	include 'atomdata.fh'

	COMMON /BAS/  P(LCUT,NCUT),PR(LCUT,NCUT),
     *	PE(LCUT,NCUT),PER(LCUT,NCUT),XE(LCUT,NCUT),
     *	QN(LCUT,NCUT),ELMIN(5,NCUT),ELMAX(5,NCUT),
     *	DELTA(5,NCUT),ITS(5,NCUT),NORLW,NOWAV,ISPIN,ISPIN1,ISPIN2
	COMMON / OU / P1(NGRID),P2(NGRID),Q1(NGRID),Q2(NGRID),
     *	BG(NGRID),BGD(NGRID),HMD(NGRID),RA(NGRID),IANS
*
      DIMENSION BGXF(4),RADF(4),COEF(4),DA(4),DB(4),XM(4),XK(4)
*
      DATA DO/0.0D0/,D1/1.0D0/,D2/2.0D0/,D3/3.0D0/,D12/0.5D0/,
     *     D6/6.0D0/,D8/8.0D0/,DMIN/1.D-14/
      DATA COEF / 0.0D0, 0.5D0, 0.5D0, 1.0D0 /
*
	C=274.0746d0
	CIN=1.d0/(C*C)
*
	JRI=JRIS(IT)
      H=STEP(IT)
      XO=RSTART(IT)
*
      H3=H/D3
      HA1=D8*H3
      HA2=1.125D0*H3
      HA3=2.375D0*H3
      HA4=0.625D0*H3
      HA5=0.125D0*H3
*
      FL=DBLE(K-1)
      FLL=DBLE(K*(K-1))
      CS=C
      CINV=D1/CS
*
      TEST=DABS(CIN-(CINV*CINV))
      IF(TEST.GT.DMIN) STOP 100
      CIN=CINV*CINV
*
      VOC1=(-BG(1)*RA(1)*CINV)
      G=(FLL+D1-VOC1*VOC1)**0.5D0
*
      Q11=  D1-G
      Q22=(-D1-G)
      P1(1)=VOC1
      Q1(1)=(G-D1)*CS
*
      DO 2 L=1,4
      FLOW1=XO+H*(DBLE(L-1)+D12)
      FLOW=DEXP(FLOW1)
      POT=(BG(L)*RA(L)+BG(L+1)*RA(L+1))*D12/FLOW
*
      RADF(1)=RA(L)
      RADF(2)=FLOW
      RADF(3)=FLOW
      RADF(4)=RA(L+1)
*
      BGXF(1)=BG(L)
      BGXF(2)=POT
      BGXF(3)=POT
      BGXF(4)=BG(L+1)
*
      XKK=DO
      XMM=DO
      DO 3 J=1,4
      PF=P1(L)+COEF(J)*XKK*H
      QF=Q1(L)+COEF(J)*XMM*H
*
      EMV=(E-BGXF(J))*RADF(J)
      X2M=RADF(J)+EMV*CIN
*
      XKK=Q11*PF+X2M*QF
      XMM=Q22*QF+(FLL/X2M-EMV)*PF
      XK(J)=XKK
      XM(J)=XMM
    3 CONTINUE
*
      P1(L+1)=P1(L)+(XK(1)+D2*(XK(2)+XK(3))+XK(4))*H/D6
      Q1(L+1)=Q1(L)+(XM(1)+D2*(XM(2)+XM(3))+XM(4))*H/D6
*
      DA(L)=Q11*P1(L+1)+X2M*Q1(L+1)
      DB(L)=Q22*Q1(L+1)+(FLL/X2M-EMV)*P1(L+1)
    2 CONTINUE
*
      SUM=DO
      DG1=DEXP(H*(G+G+D1))
      RG1=D1
*
      DO 4 L=5,JRI
      EMV=(E-BG(L))*RA(L)
      X2M=RA(L)+EMV*CIN
*
      AKK=P1(L-4)+HA1*(DA(3)-D12*DA(2)+DA(1))
      BKK=Q1(L-4)+HA1*(DB(3)-D12*DB(2)+DB(1))
*
      DA(4)=Q11*AKK+X2M*BKK
      DB(4)=Q22*BKK+(FLL/X2M-EMV)*AKK
*
      P1(L)=P1(L-1)+HA2*DA(4)+HA3*DA(3)-HA4*DA(2)+HA5*DA(1)
      Q1(L)=Q1(L-1)+HA2*DB(4)+HA3*DB(3)-HA4*DB(2)+HA5*DB(1)
*
      RG1=RG1*DG1
      SUM=SUM+P1(L)*P1(L)*RG1
*
      DA(1)=DA(2)
      DB(1)=DB(2)
      DA(2)=DA(3)
      DB(2)=DB(3)
      DA(3)=Q11*P1(L)+X2M*Q1(L)
      DB(3)=Q22*Q1(L)+(FLL/X2M-EMV)*P1(L)
    4 CONTINUE
*
      RG=D1/(SUM*H*RA(6))**0.5D0
*
      RETURN
      END
*----------------------------------------------------------
      SUBROUTINE DIRE(E,K,IT,RG)
*
*     Energy Derivatives of D.D.Koelling & B.N.Harmon eq.
*                      (Journal of Physics C v.10 N.16 1977)
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
	include 'atomdata.fh'

	COMMON /BAS/  P(LCUT,NCUT),PR(LCUT,NCUT),
     *	PE(LCUT,NCUT),PER(LCUT,NCUT),XE(LCUT,NCUT),
     *	QN(LCUT,NCUT),ELMIN(5,NCUT),ELMAX(5,NCUT),
     *	DELTA(5,NCUT),ITS(5,NCUT),NORLW,NOWAV,ISPIN,ISPIN1,ISPIN2
	COMMON / OU / P1(NGRID),P2(NGRID),Q1(NGRID),Q2(NGRID),
     *	BG(NGRID),BGD(NGRID),HMD(NGRID),RA(NGRID),IANS
*
      DIMENSION BGXF(4),RADF(4),COEF(4),DA(4),DB(4),XM(4),XK(4)
*
      DATA DO/0.0D0/,D1/1.0D0/,D2/2.0D0/,D3/3.0D0/,D12/0.5D0/,
     *     D6/6.0D0/,D8/8.0D0/,DMIN/1.D-14/
      DATA COEF / 0.0D0, 0.5D0, 0.5D0, 1.0D0 /
*
	C=274.0746d0
	CIN=1.d0/(C*C)
*
      JRI=JRIS(IT)
      H=STEP(IT)
      XO=RSTART(IT)
*
      H3=H/D3
      HA1=D8*H3
      HA2=1.125D0*H3
      HA3=2.375D0*H3
      HA4=0.625D0*H3
      HA5=0.125D0*H3
*
      FL=DBLE(K-1)
      FLL=DBLE(K*(K-1))
      CS=C
      CINV=D1/CS
*
      TEST=DABS(CIN-(CINV*CINV))
      IF(TEST.GT.DMIN) STOP 100
      CIN=CINV*CINV
*
      VOC1=(-BG(1)*RA(1)*CINV)
      G=(FLL+D1-VOC1*VOC1)**0.5D0
*
      Q11=  D1-G
      Q22=(-D1-G)
*
      P2(1)=DO
      Q2(1)=DO
*
      DO 2 L=1,4
      FLOW1=XO+H*(DBLE(L-1)+D12)
      FLOW=DEXP(FLOW1)
      POT=(BG(L)*RA(L)+BG(L+1)*RA(L+1))*D12/FLOW
*
      RADF(1)=RA(L)
      RADF(2)=FLOW
      RADF(3)=FLOW
      RADF(4)=RA(L+1)
*
      BGXF(1)=BG(L)
      BGXF(2)=POT
      BGXF(3)=POT
      BGXF(4)=BG(L+1)
*
      XKK=DO
      XMM=DO
      DO 3 J=1,4
      PF=P2(L)+COEF(J)*XKK*H
      QF=Q2(L)+COEF(J)*XMM*H
      PF1=P1(L)+COEF(J)*(P1(L+1)-P1(L))
      QF1=Q1(L)+COEF(J)*(Q1(L+1)-Q1(L))
*
      EMV=(E-BGXF(J))*RADF(J)
      X2M=RADF(J)+EMV*CIN
      FLOW=FLL/X2M
*
      ADD1=RADF(J)*QF1*CIN
      ADD2=((FLOW/X2M)*CIN+D1)*RADF(J)*PF1
*
      XKK=Q11*PF+X2M*QF+ADD1
      XMM=Q22*QF+(FLOW-EMV)*PF-ADD2
      XK(J)=XKK
      XM(J)=XMM
    3 CONTINUE
*
      P2(L+1)=P2(L)+(XK(1)+D2*(XK(2)+XK(3))+XK(4))*H/D6
      Q2(L+1)=Q2(L)+(XM(1)+D2*(XM(2)+XM(3))+XM(4))*H/D6
*
      DA(L)=Q11*P2(L+1)+X2M*Q2(L+1)+ADD1
      DB(L)=Q22*Q2(L+1)+(FLOW-EMV)*P2(L+1)-ADD2
    2 CONTINUE
*
      DO 4 L=5,JRI
      EMV=(E-BG(L))*RA(L)
      X2M=RA(L)+EMV*CIN
      FLOW=FLL/X2M
*
      ADD1=RA(L)*Q1(L)*CIN
      ADD2=((FLOW/X2M)*CIN+D1)*RA(L)*P1(L)
*
      AKK=P2(L-4)+HA1*(DA(3)-D12*DA(2)+DA(1))
      BKK=Q2(L-4)+HA1*(DB(3)-D12*DB(2)+DB(1))
*
      DA(4)=Q11*AKK+X2M*BKK+ADD1
      DB(4)=Q22*BKK+(FLOW-EMV)*AKK-ADD2
*
      P2(L)=P2(L-1)+HA2*DA(4)+HA3*DA(3)-HA4*DA(2)+HA5*DA(1)
      Q2(L)=Q2(L-1)+HA2*DB(4)+HA3*DB(3)-HA4*DB(2)+HA5*DB(1)
*
      DA(1)=DA(2)
      DB(1)=DB(2)
      DA(2)=DA(3)
      DB(2)=DB(3)
      DA(3)=Q11*P2(L)+X2M*Q2(L)+ADD1
      DB(3)=Q22*Q2(L)+(FLOW-EMV)*P2(L)-ADD2
    4 CONTINUE
*
      DG=DEXP(H*G)
*
      DO 5 L=1,JRI
      RG=RG*DG
*
      P1(L)=P1(L)*RG
      Q1(L)=Q1(L)*RG
*
      P2(L)=P2(L)*RG
      Q2(L)=Q2(L)*RG
    5 CONTINUE
*
      RETURN
      END

      SUBROUTINE XSUM(SUM, INDA1,INDA2,INDB,INDEX)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
*
	include 'atomdata.fh'

	COMMON / OU / P(NGRID,2,2),	! concat of P1(NGRID),P2(NGRID),Q1(NGRID),Q2(NGRID),
     *	BG(NGRID),BGD(NGRID),HMD(NGRID),RA(NGRID),IANS
*
      DIMENSION F(300)
*
      JRI=JRIS(INDEX)
      H=STEP(INDEX)
      DO 5 I=1,JRI
*   5 F(I)=P(I,INDA1,INDB)*P(I,INDA2,INDB)*RAD(I,INDEX)
    5 F(I)=P(I,INDA1,INDB)*P(I,INDA2,INDB)*RA(I)
      SUM=DIADL2(F,JRI,H,1,H)
      RETURN
      END
*
