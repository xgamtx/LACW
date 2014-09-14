*-------------------------------------------------------------------------
*      DOUBLE PRECISION FUNCTION DIADL2(X,N,H,NMX,DX)
*      DOUBLE PRECISION FUNCTION DIADL1_(X,Y,N,M,H,M1,M2,NMX,DX)
*      DOUBLE PRECISION FUNCTION ATINT_(X,Y,N,M,H,M1,M2,NMX,DX)
*      SUBROUTINE POISON(PSQ,M,M1,M2,J,W,Z,XO,H,E,F)
*      SUBROUTINE SXCPOT(A)
*      SUBROUTINE XCPOT(RHO1,RHO2,RHO,V1,V2,EXC,A)
*-------------------------------------------------------------------------
      DOUBLE PRECISION FUNCTION DIADL2(X,N,H,NMX,DX)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION X(N)
      DATA  FF1/0.2D1/,FF2/0.4D1/,FF3/0.0D0/,FF4/0.3D1/
      IF(NMX.EQ.1) THEN
         BEG=FF3
      ELSE
         BEG=X(1)/DX
      ENDIF
      NTEST=MOD(N,2)
      IF(NTEST.EQ.1) THEN
         ISTART=1
      ELSE
         ISTART=2
         BEG=BEG+(X(1)+X(2))*H/FF1
      ENDIF
      SA=FF3
      SB=FF3
      DO 1 K=ISTART+2,N,2
      SA=SA+X(K)
    1 SB=SB+X(K-1)
      DIADL2=H*(FF1*SA+FF2*SB-X(N)+X(ISTART))/FF4+BEG
      RETURN
      END

      DOUBLE PRECISION FUNCTION DIADL1_(X,Y,N,M,H,M1,M2,NMX,DX)
      INTEGER N,M,K,M1,M2,NMX,ISTART
      DOUBLE PRECISION X(M1,M2),Y(M1,M2),H,FF1,FF2,FF3,FF4,DX,BEG,SA,SB
      DATA FF1/.2D1/,FF2/.4D1/,FF3/.0D0/,FF4/.3D1/
      IF(NMX.EQ.1) THEN
         BEG=FF3
       ELSE
         BEG=X(1,M)*Y(1,M)/DX
      ENDIF
      IF(MOD(N,2).EQ.1) THEN
         ISTART=1
       ELSE
         ISTART=2
         BEG=BEG+(X(1,M)*Y(1,M)+X(2,M)*Y(2,M))*H/FF1
      ENDIF
      SA=FF3
      SB=FF3
      DO K=ISTART+2,N,2
        SA=SA+X(K,M)*Y(K,M)
        SB=SB+X(K-1,M)*Y(K-1,M)
      END DO
      DIADL1_=H*(FF1*SA+FF2*SB-X(N,M)*Y(N,M)+X(ISTART,M)*Y(ISTART,M))/
     .        FF4+BEG
      RETURN
      END
*
*
      DOUBLE PRECISION FUNCTION ATINT_(X,Y,N,M,H,M1,M2,NMX,DX)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER N,M,K,M1,M2,NMX
      DOUBLE PRECISION X(M1,M2),Y(M1,M2),H,DX
      DATA  D2/2.D0/,D4/4.D0/,D3/3.D0/
*
      START=0.0d0
      IF(NMX.EQ.2) DX=DLOG(X(2,M)*Y(2,M)/(X(1,M)*Y(1,M)))/H
      IF(NMX.GT.1) START=X(1,M)*Y(1,M)/DX
      SUM1=0.0d0
      SUM2=0.0d0
      DO 1 K=2,N,2
      SUM1=SUM1+X(K-1,M)*Y(K-1,M)
    1 SUM2=SUM2+X(K,M)*Y(K,M)
      ATINT_=H*(D2*SUM1+D4*SUM2-X(1,M)*Y(1,M))/D3+START
      RETURN
      END
*
*
      SUBROUTINE POISON(PSQ,M,M1,M2,J,W,Z,XO,H,E,F)
      INTEGER M1,M2,J,I,ITOP,JV
      DOUBLE PRECISION PSQ(M1,M2),W(*),E(*),F(*),Z,XO,H,
     +                 A,B,EDL,C,C2,G,X
      A=1.D0-H*H/48.D0
      B=-2.D0-H*H*10.D0/48.D0
      EDL=DEXP(H*.5D0)
      C=H*H/(2.D0*3.D0)
      C2=-B/A
      E(1)=0.D0
      F(1)=DEXP(H*.5D0)
      X=XO+H
      ITOP=J-1
*
      DO 5 I=2,ITOP
      G=C*DEXP(X*.5D0)*(EDL*PSQ(I+1,M)+10.D0*PSQ(I,M)+PSQ(I-1,M)/EDL)
      F(I)=C2-1.D0/F(I-1)
      E(I)=(G/A+E(I-1))/F(I)
 5    X=X+H
      W(J)=2.D0*Z*DEXP(-X*.5D0)
*
      DO 6 I=1,ITOP
      JV=J-I
  6   W(JV)=E(JV)+W(JV+1)/F(JV)
      RETURN
      END
*
      SUBROUTINE SXCPOT(A)
*-------------------------------------------------------------*
*     Set Constants for <XCPOT>                               *
*-------------------------------------------------------------*
      DOUBLE PRECISION PI
      PARAMETER( PI=3.14159265358979323846264338328D0 )
      STRUCTURE /cxcpot/
        DOUBLE PRECISION XCCP,XCCF,XCRP,XCRF,XALPHA,OTH,FTH,AA,BB,ALPHA
        INTEGER IXCH
      END STRUCTURE
      RECORD /cxcpot/ A
*
      GO TO (31,32,33,34),A.IXCH+1
*
*     BARTH-HEDIN J. PHYS. C5,1629 (1972)
*
   31 A.XCCP=0.0504D0
      A.XCCF=0.0254D0
      A.XCRP=30.D0
      A.XCRF=75.D0
      A.OTH=1.D0/3.D0
      A.FTH=4.D0/3.D0
      A.AA=1.D0/2.D0**A.OTH
      A.BB=1.D0-A.AA
      RETURN
*
*     SLATER X-ALPHA
*
   32 A.OTH=1.D0/3.D0
      A.XALPHA=6.D0*A.ALPHA*(3.D0/(16.D0*PI**2.D0))**A.OTH
      RETURN
*
*     BARTH-HEDIN-JANAK PHYS. REV. B12,1257 (1975)
*
   33 A.XCCP=0.045D0
      A.XCCF=0.0225D0
      A.XCRP=21.D0
      A.XCRF=53.D0
      A.OTH=1.D0/3.D0
      A.FTH=4.D0/3.D0
      A.AA=1.D0/2.D0**A.OTH
      A.BB=1.D0-A.AA
      RETURN
*
*     VOSKO-WILK-NUSAIR CAN. J. PHYS. 58,1200 (1980)
*
   34 A.OTH=1.D0/3.D0
      A.FTH=4.D0/3.D0
      A.AA=2.D0**A.FTH-2.D0
      RETURN
      END
*
      SUBROUTINE XCPOT(RHO1,RHO2,RHO,V1,V2,EXC,A)
*-------------------------------------------------------------*
*     Calculate Exchange Correlation Potential                *
*-------------------------------------------------------------*
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      STRUCTURE /cxcpot/
        DOUBLE PRECISION XCCP,XCCF,XCRP,XCRF,XALPHA,OTH,FTH,AA,BB,ALPHA
        INTEGER IXCH
      END STRUCTURE
      RECORD /cxcpot/ A
      DATA AP,XPO,BP,CP,QP,CP1,CP2,CP3
     1  /0.0621814D0,-0.10498D0,3.72744D0,12.9352D0,6.1519908D0,
     2   1.2117833D0, 1.1435257D0,-0.031167608D0/
      DATA AF,XFO,BF,CF,QF,CF2,CF3
     1  /0.0310907D0,-0.32500D0,7.060428D0,18.0578D0,4.7309269D0,
     2   2.7100059D0,-0.1446006D0/
*
      RS1=(RHO/3.D0)**A.OTH
      RS=1.D0/RS1
      GO TO (21,22,21,23),A.IXCH+1
*
*     BARTH-HEDIN EXCHANGE CORRELATION
*     J. PHYS. C5,1629(1972)
*
   21 RSF=RS/A.XCRF
      RSF2=RSF*RSF
      RSF3=RSF2*RSF
      RSP=RS/A.XCRP
      RSP2=RSP*RSP
      RSP3=RSP2*RSP
      FCF=(1.D0+RSF3)*DLOG(1.D0+1.D0/RSF)+0.5D0*RSF-RSF2-A.OTH
      FCP=(1.D0+RSP3)*DLOG(1.D0+1.D0/RSP)+0.5D0*RSP-RSP2-A.OTH
      EPSCP=-A.XCCP*FCP
      EPSCF=-A.XCCF*FCF
      EPSXP=-0.91633059D0/RS
      CNY=5.1297628D0*(EPSCF-EPSCP)
      X=RHO1/RHO
      FX=(X**A.FTH+(1.D0-X)**A.FTH-A.AA)/A.BB
      EXC=EPSXP+EPSCP+FX*(CNY+A.FTH*EPSXP)/5.1297628D0
      ARS=-1.22177412D0/RS+CNY
      BRS=-A.XCCP*DLOG(1.D0+A.XCRP/RS)-CNY
      TRX1=(2.D0*X)**A.OTH
      V1=ARS*TRX1+BRS
      TRX2=(2.D0*RHO2/RHO)**A.OTH
      V2=ARS*TRX2+BRS
      RETURN
*
*     SLATER EXCHANGE POTENTIAL
*
   22 EXC=-0.75D0*A.XALPHA*((0.5D0*RHO)**A.OTH)
      V1=-A.XALPHA*((RHO1)**A.OTH)
      V2=-A.XALPHA*((RHO2)**A.OTH)
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
      FS=(SP**A.FTH+SM**A.FTH-2.D0)/A.AA
      BETA=1.D0/(2.74208D0+3.182D0*X+0.09873D0*X*X+0.18268D0*X**3)
      DFS=A.FTH*(SP**A.OTH-SM**A.OTH)/A.AA
      DBETA=-(0.27402D0*X+0.09873D0+1.591D0/X)*(BETA**2)
      ATNP=DTAN(QP/(2.D0*X+BP))
      ATNF=DTAN(QF/(2.D0*X+BF))
      ECP=AP*(DLOG(X*X/XPX)+CP1*ATNP-CP3*(DLOG((X-XPO)**2/XPX)+
     1    CP2*ATNP))
      ECF=AF*(DLOG(X*X/XFX)+CP1*ATNF-CF3*(DLOG((X-XFO)**2/XFX)+
     1    CF2*ATNF))
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
      EPX=-0.91633059D0/RS*(1.D0+A.FTH*FS/5.1297628D0)
      AMYX2=-1.22177412D0/RS*(SP**A.OTH)
      AMYX1=-1.22177412D0/RS*(SM**A.OTH)
      EXC=EC+EPX
      V1=UC1+AMYX1
      V2=UC2+AMYX2
      RETURN
      END

	integer function IsInternalEx(x,y,z,dists_OUT,
     +	nt,nts,rmts,rsx,rsy,rsz,NQG)

	implicit double precision(a-h,o-z)

	integer ret
	double precision dists_OUT(NQG,*)

	double precision rsx(NQG,*),rsy(NQG,*),rsz(NQG,*),rmts(*)
	integer nts(*)

	ret=0
	do it=1,nt
	do is=1,nts(it)
	  dist=dist3d(x,y,z,rsx(is,it),rsy(is,it),rsz(is,it))
	  dists_OUT(is,it)=dist
	  if (dist.lt.rmts(it)) ret=ret+1
	enddo
	enddo
	IsInternalEx=ret
	return
	end


! Find min & max coords of the system
	subroutine StructBounds(nt,nts,rsx,rsy,rsz,NQG,
     +	xmin,ymin,zmin,xmax,ymax,zmax)
	implicit double precision(a-h,o-z)
	integer nt,NQG,nts(*)
	double precision rsx(NQG,*),rsy(NQG,*),rsz(NQG,*)
        double precision xmin,ymin,zmin,xmax,ymax,zmax

	xmin=1.0d+12
	ymin=1.0d+12
	zmin=1.0d+12
	xmax=-1.0d+12
	ymax=-1.0d+12
	zmax=-1.0d+12

	do it=1,nt
	do is=1,nts(it)
	   if (rsx(is,it).lt.xmin) xmin=rsx(is,it)
	   if (rsy(is,it).lt.ymin) ymin=rsy(is,it)
	   if (rsz(is,it).lt.zmin) zmin=rsz(is,it)
	   if (rsx(is,it).gt.xmax) xmax=rsx(is,it)
	   if (rsy(is,it).gt.ymax) ymax=rsy(is,it)
	   if (rsz(is,it).gt.zmax) zmax=rsz(is,it)
	enddo
	enddo

	return
	end


	double precision function
     +	CountElectrons(rad,amtc,NGRID,mgrid,nt,nts,step)

	implicit double precision(a-h,o-z)
	integer NGRID,mgrid,nt,nts(*)
	double precision rad(NGRID,*),amtc(NGRID,*),step(*)

	s=0.0d0
	do it=1,nt
	  s=s+dble(nts(it))*
     +		ATINT_(amtc,rad,mgrid,it,step(it),NGRID,it,1,1.d0)
	enddo
        CountElectrons=s
c Atint_() takes the second dimension of rad&amtc as the 7th parameter
c but does not use it, so I pass ''it'', which will only mean that
c there's surely not less than ''it'' lines in the arrays.

	return
	end


      double precision function FindYatom(xs,ys,n,x,it,index_RETURN)
*     ----------------------------------------------------------
*     Finds y(x), where the function is passed as SORTED arrays
*     of points (x,y) -- xs[n,?] and ys[n,?]
*     ----------------------------------------------------------

      implicit double precision(a-h,o-z)
      dimension xs(n,*),ys(n,*)

      data dmin/1.0d-12/

      index_RETURN=1
      FindYatom=x
      if (n.le.0) return        !then
cd         write(6,*) 'FindYatom() emergency return, n=',n
c         return
c      endif

      m=1
      do i=1,n-1
         if (x.ge.xs(i,it).and.x.lt.xs(i+1,it)) m=i
      enddo
c      if (x.lt.xs(1,it)) m=1
      if (n.gt.1.and.x.ge.xs(n,it)) m=n-1

      ret=ys(m,it)
      if (dabs(xs(m+1,it)-xs(m,it)).ge.dmin) ret=ys(m,it)+
     . (ys(m+1,it)-ys(m,it))*(x-xs(m,it))/(xs(m+1,it)-xs(m,it))

d      write(26,*) 'm,x,ys(m),ret:',m,x,ys(m,it),ret

      index_RETURN=m
      FindYatom=ret
      return
      end     !//FindYatom()
