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
*
      DOUBLE PRECISION FUNCTION F_NJM(RO,KP,M)
*        Calc. sub-integral value for integr. NJm
      DOUBLE PRECISION RO,KP,BJM,DUMMY
      INTEGER M
      CALL BESS(M,KP*RO,BJM,DUMMY)
      F_NJM=BJM*BJM*RO
      RETURN
      END

      DOUBLE PRECISION FUNCTION SIMP_NJM(KP,M,A,Eps)
*     Integr. I1 - calc. by Simpson's method. Uses F_NJM
      EXTERNAL F_NJM
      DOUBLE PRECISION KP,Eps,A,s1,s2,s3,h,h2,x,F_NJM
      INTEGER I,M,JMAX
      PARAMETER (JMAX=500)
      I=0
      M=IABS(M)
      h=A/2.d0  ! low integr. limit is 0
      SIMP_NJM=2.d0*F_NJM(h,KP,M)
      s1=F_NJM(0.d0,KP,M)+F_NJM(A,KP,M)+SIMP_NJM
      SIMP_NJM=s1+SIMP_NJM
      s3=SIMP_NJM
      DO WHILE(DABS(1.D0-SIMP_NJM/(2.d0*s3)).gt.Eps)
          s2=0.D0
          s3=SIMP_NJM
          h2=h
          h=h/2.d0
          DO x=0.D0+h,A,h2
            s2=s2+F_NJM(x,KP,M)
          END DO
          s2=2.d0*s2
          s1=s1+s2
          SIMP_NJM=s1+s2
          IF(I.gt.JMAX) THEN
            WRITE(*,'(1x,a32,i5)') 'SIMP_NJM---> limit of iterations',
     .            JMAX
            STOP
          ENDIF
          I=I+1
      END DO
      SIMP_NJM=1.D0/DSQRT(SIMP_NJM*h/3.d0)
      RETURN
      END
*
      DOUBLE PRECISION FUNCTION S_NJM(VKP1,VKP2,M1,M2,A)
      DOUBLE PRECISION VKP1,VKP2,A,BSJ1,BSJ2,dummy
      INTEGER M1,M2
      CALL BESS(M1,VKP1*A,dummy,BSJ1)
      CALL BESS(M2,VKP2*A,dummy,BSJ2)
      S_NJM=2.0D0/(A*A*DABS(BSJ1*BSJ2))
d      write(16,*) 's_njm:',m1,vkp1,bsj1,m2,vkp2,bsj2
      RETURN
      END
*
*! begin fragment SP 3.1
      DOUBLE PRECISION FUNCTION S_NYM(VKP1,VKP2,M1,M2,B,Eps)
      DOUBLE PRECISION VKP1,VKP2,B,BS1,BS2,Eps,dummy
      INTEGER M1,M2
      CALL BESSYPROC(M1,VKP1*B,Eps,dummy,BS1)
      CALL BESSYPROC(M2,VKP2*B,Eps,dummy,BS2)
      S_NYM=2.0D0/(B*B*DABS(BS1*BS2))
d      write(16,*) 's_nym:',m1,b*bs1,m2,b*bs2
      RETURN
      END
*! end fragment SP 3.1
*
*
*$ begin fragment
*
* The S_CJY() function calculates the coefficient to matrix elements
* in cored LACW model (analogue of S_NJM() of cylinder LACW).
*
* abXX arrays are filled by fkappa.for::FindBessAtAB(),
* pdKappas are obtained from FindKappas() from the same source file.
*
      double precision function S_CJY(a,b,m,n,m1,n1,
     . Mmax,Nmax,abJ,abJD,abY,abYD)

      implicit double precision(a-h,o-z)
c      dimension pdKappas(Mmax,Nmax)
      dimension abJ(Mmax,Nmax,2),abJD(Mmax,Nmax,2),
     . abY(Mmax,Nmax,2),abYD(Mmax,Nmax,2)

*
* Here tmp  = Jm(xa)/Ym(xa) = Jm(xb)/Ym(xb), x=kappa(m,n)
* I use ChooseBessAtAB() to retain the (-1)^m <--- if m<0
* Another way:
*
*      tmp=abJ(iabs(m)+1,n,1)/abY(iabs(m)+1,n,1)
*
* Normally, we'd have to apply this multiplication (by (-1)^m)
* to both abJ and abY here.
*
*      if (m.lt.0.and.mod(m,2).ne.0) ...
*
* Here we can write ....n,2, since Jm(xa)/Ym(xa) = Jm(xb)/Ym(xb)
*
      tmp=abJ(iabs(m)+1,n,1)/abY(iabs(m)+1,n,1)
      tmp1=abJ(iabs(m1)+1,n1,1)/abY(iabs(m1)+1,n1,1)

      m_sav=m
      m1_sav=m1
      m=iabs(m)+1
      m1=iabs(m1)+1
c      tmp=ChooseBessAtAB(m_sav,n,MMAX,NMAX,abJ,1)/abY(m,n,1)
c      tmp1=ChooseBessAtAB(m1_sav,n1,MMAX,NMAX,abJ,1)/abY(m1,n1,1)
      b1=a*(abJD(m,n,1)-abYD(m,n,1)*tmp)
      b2=b*(abJD(m,n,2)-abYD(m,n,2)*tmp)
      b11=a*(abJD(m1,n1,1)-abYD(m1,n1,1)*tmp1)
      b21=b*(abJD(m1,n1,2)-abYD(m1,n1,2)*tmp1)
      S_CJY=2.0d0/dsqrt(dabs((b1*b1-b2*b2)*(b11*b11-b21*b21))) 
      m=m_sav
      m1=m1_sav
d      write(23,*) 's_cjy:',m,n,tmp,b1,b2,m1,n1,tmp1,b11,b21

      return
      end
*$ end fragment
*$ begin fragment
*
* The S_CJYCR() function calculates the coefficient to matrix elements
* in cored LACW model (analogue of S_NJM() of cylinder LACW).
*
* abXX arrays are filled by fkappa.for::FindBessAtAB(),
* pdKappas are obtained from FindKappas() from the same source file.
*
      double precision function S_CJYCR(a,b,m,n,m1,n1,
     . Mmax,Nmax,abJ,abJD,abY,abYD,pdKappas,v,intCR,intCRC)

      implicit double precision(a-h,o-z)
      dimension pdKappas(Mmax,Nmax)
      dimension abJ(Mmax,Nmax,2),abJD(Mmax,Nmax,2),
     . abY(Mmax,Nmax,2),abYD(Mmax,Nmax,2)
	double precision trap_int_inf, tii1, tii2, v
      double precision intCR(MMAX, NMAX)
      integer          intCRC(MMAX, NMAX)


      tmp=abJ(iabs(m)+1,n,2)/abY(iabs(m)+1,n,2)
      tmp1=abJ(iabs(m1)+1,n1,2)/abY(iabs(m1)+1,n1,2)

      m_sav=m
      m1_sav=m1
      m=iabs(m)+1
      m1=iabs(m1)+1

	if(intCRC(m,n).eq.1)then
	tii1 = intCR(m,n)
	else
      tii1 = trap_int_inf(a, dsqrt(v - 
     *pdKappas(m,n)*pdKappas(m,n)), m - 1)
	intCR(m,n) = tii1
	intCRC(m,n) = 1
	endif

	if(intCRC(m1,n1).eq.1)then
	tii2 = intCR(m1,n1)
	else
      tii2 = trap_int_inf(a, dsqrt(v - 
     *pdKappas(m1,n1)*pdKappas(m1,n1)), m1 - 1)
	intCR(m1,n1) = tii2
	intCRC(m1,n1) = 1
	endif


      ba=abJ(m,n,1)-abY(m,n,1)*tmp
      b1=a*(abJD(m,n,1)-abYD(m,n,1)*tmp)
      b2=b*(abJD(m,n,2)-abYD(m,n,2)*tmp)
      ba1=abJ(m1,n1,1)-abY(m1,n1,1)*tmp1
      b11=a*(abJD(m1,n1,1)-abYD(m1,n1,1)*tmp1)
      b21=b*(abJD(m1,n1,2)-abYD(m1,n1,2)*tmp1)

      S_CJYCR=2.0d0/dsqrt(dabs((b1*b1-b2*b2)*
     *tii1*(b11*b11-b21*b21)*tii2)) 

      S_CJYCR=2.0d0/dsqrt(dabs((b1*b1-b2*b2)*
     *(b11*b11-b21*b21))) 

	S_CJYCR=1.0d0/dsqrt(
     *dabs((b1*b1 - b2*b2)/2.0d0 + ba*ba*tii1)*
     *dabs((b11*b11 - b21*b21)/2.0d0 + ba1*ba1*tii2))

      m=m_sav
      m1=m1_sav
d      write(23,*) 's_cjy:',m,n,tmp,b1,b2,m1,n1,tmp1,b11,b21

      return
      end
*$ end fragment
*$ begin fragment
*
* The S_CJYCRI() function calculates the coefficient to matrix elements
* in cored LACW model (analogue of S_NJM() of cylinder LACW).
*
* abXX arrays are filled by fkappa.for::FindBessAtAB(),
* pdKappas are obtained from FindKappas() from the same source file.
*
      double precision function S_CJYCRI(a,b,m,n,m1,n1,
     . Mmax,Nmax,abJ,abJD,abY,abYD,pdKappas,v,intCR,intCRC)

      implicit double precision(a-h,o-z)
      dimension pdKappas(Mmax,Nmax)
      dimension abJ(Mmax,Nmax,2),abJD(Mmax,Nmax,2),
     . abY(Mmax,Nmax,2),abYD(Mmax,Nmax,2)
	double precision trap_int_infI, tii1, tii2, v
      double precision intCR(MMAX, NMAX)
      integer          intCRC(MMAX, NMAX)


      tmp=abJ(iabs(m)+1,n,1)/abY(iabs(m)+1,n,1)
      tmp1=abJ(iabs(m1)+1,n1,1)/abY(iabs(m1)+1,n1,1)

      m_sav=m
      m1_sav=m1
      m=iabs(m)+1
      m1=iabs(m1)+1

	if(intCRC(m,n).eq.1)then
	tii1 = intCR(m,n)
	else
      tii1 = trap_int_infI(b, dsqrt(v - 
     *pdKappas(m,n)*pdKappas(m,n)), m - 1)
	intCR(m,n) = tii1
	intCRC(m,n) = 1
	endif

	if(intCRC(m1,n1).eq.1)then
	tii2 = intCR(m1,n1)
	else
      tii2 = trap_int_infI(b, dsqrt(v - 
     *pdKappas(m1,n1)*pdKappas(m1,n1)), m1 - 1)
	intCR(m1,n1) = tii2
	intCRC(m1,n1) = 1
	endif


      ba=abJ(m,n,2)-abY(m,n,2)*tmp
      b1=a*(abJD(m,n,1)-abYD(m,n,1)*tmp)
      b2=b*(abJD(m,n,2)-abYD(m,n,2)*tmp)
      ba1=abJ(m1,n1,2)-abY(m1,n1,2)*tmp1
      b11=a*(abJD(m1,n1,1)-abYD(m1,n1,1)*tmp1)
      b21=b*(abJD(m1,n1,2)-abYD(m1,n1,2)*tmp1)

c      S_CJYCRI=2.0d0/dsqrt(dabs((b1*b1-b2*b2)*
c     *tii1*(b11*b11-b21*b21)*tii2)) 
c
c      S_CJYCRI=2.0d0/dsqrt(dabs((b1*b1-b2*b2)*
c     *(b11*b11-b21*b21))) 

	S_CJYCRI=1.0d0/dsqrt(
     *dabs((b1*b1 - b2*b2)/2.0d0 + ba*ba*tii1)*
     *dabs((b11*b11 - b21*b21)/2.0d0 + ba1*ba1*tii2))

      m=m_sav
      m1=m1_sav
d      write(23,*) 's_cjy:',m,n,tmp,b1,b2,m1,n1,tmp1,b11,b21

      return
      end
*$ end fragment
*
      DOUBLE PRECISION FUNCTION SIMP_NJM1(KP,M,A)
      DOUBLE PRECISION KP,A,tmp,dummy
      INTEGER M
      CALL BESS(M,KP*A,dummy,TMP)
      SIMP_NJM1=DSQRT(2.D0)/(A*DABS(TMP))
      RETURN
      END
