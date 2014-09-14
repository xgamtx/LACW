* --------------------------------------------------------------------- *
*                  CYLINDER LAPW Pre-Calculation MODULE                 *
*            (c) 1998 Kepp O.M., Nikolaev A.V., D'yachkov P.N.          *
* --------------------------------------------------------------------- *
*   This module uses routines, imported from A.V.Nikolaev's             *
*   STRC program for bar LAPW calculations.                             *
*   These routines are:                                                 *
*                      DIADL1, POISON, SXCPOT, XCPOT;                   *
*                      VECREL, SUMAX - adapted for this program release *
* --------------------------------------------------------------------- *
*									*
*	Tubular Full-potential LACW - Service Pack 5.10			*
*	(c) 1999-2001-? by Dmitry V. Kirin & Pavel N. D'yachkov		*
*									*
* --------------------------------------------------------------------- *
*	`The Universe is run by the complex interweaving of three	*
*	 elements: energy, matter, and enlightened self-interest.'	*
*		--G'Kar, `Babylon 5' episode 1.11: `Survivors'		*
* --------------------------------------------------------------------- *
*
*! begin fragment
	implicit double precision(a-h,o-z)
*! end
	INTEGER I,J,M,N,NP,IT,JT	!,K,NT,NQG,NCUT,NGRID,MAXMNP,MUL

*! begin fragment SP5.0
	include 'atomdata.fh'
	include 'coordsys.fh'
	include 'projstrc.fh'
	include 'strresult.fh'
	record /TBasis/ Basis

	include 'fkappa.fh'
	include 'stdio.fh'
*! end fragment SP5.0

      DOUBLE PRECISION A,C,GMAX,
     +	WPOT1(NGRID),
     +	AMTC0(NGRID,NCUT),
     +	RHO,DUMMY,TMP,CORX,Ecache(NGRID),Fcache(NGRID),
     +	FLOWX,FLOWY,FLOWZ,DR(20),RMAX,ADN(10,NCUT,NCUT,NQG),
     +	XNAN(10,NCUT,NCUT,NQG),BSX(3),BSY(3),BSZ(3),
     +	ENER(MUL*MAXMNP),RES(MAXMNP),KPNM,VKP
      INTEGER NSR(20),NSH,
     +	JQJ(NCUT,NCUT,NQG),INDX(3,MUL*MAXMNP),COUNT,COUNTX,
     +	INDXX(4,MUL*MAXMNP),IND(MUL*MAXMNP)

*! begin fragment SP6.0
      DOUBLE PRECISION
     +   ENKP(MAXMNP,MAXMNP),
     +  CalcCJ, CalcCY, CalcCCR, 
     +  resCJ, resCY, resCCR, resKmn
*! end fragment SP6.0


c	external DIADL1_
c	double precision DIADL1_

*! begin fragment SP5.0
	target INDXX
*! end fragment SP5.0

*! begin fragment
	parameter (fk_DX=0.125d0, fk_XSTART=0.1d0) !, fk_XEND=120.0d0)
	! fk_XXXXXX parameters are used by FindKappa() to find zeros of
	! zerofunc(). The choice of fk_XSTART is restricted by the fact
	! that Ym(x) functions are -INFINITY at x=0; fk_XEND may be big;
	! fk_DX will be decreased automatically (but don't make its initial
	! value too big, for the search routine may fail.

	dimension EnerSort(MUL*MAXMNP)

	double precision B	  ! the tube internal radius
	double precision ChiralityAngle
	integer iUseCoredTube	
!     0 (Cylindrical LACW), 1 (Tubular LACW), 2 (Inverted LACW), 
!     or 3 (Tubular LACW in crystal)
*! end
	DATA BSX/1.d0,0.d0,0.d0/,BSY/0.d0,1.d0,0.d0/,BSZ/0.d0,0.d0,1.d0/
	DATA RMAX/1.5D0/
	DATA COUNTX/1/

c* legacy input file format parsing code requires some additional var defs:
c	include 'oldread_strcy.fh'
*
********
	A=-1.0d0
	B=-1.0d0
	C=-1.0d0
	V=-1.0d0
	ChiralityAngle=0.0d0
	iUseCoredTube=0
        Project.bAddPCompanion=.FALSE.
	zf_params.eps=1.0d-12
	B_SMALL=1.0d-5        ! B<B_SMALL are treated as zero
	PHI_SMALL=1.0d-8

	Co.IXCH=-1	! invalid exchange mode (was ==5 before SP 4.0)

	ierr=16
	DEBUG=8

* open the logfile:
	open(16,file='strcyout.log',form='formatted')	!,err=62)

	write(*,*)'**** STRCY **** LACW Service Pack 6.1 (based on 5.07)'
	write(16,*)'**** STRCY **** LACW Service Pack 6.1 (based on 5.07)'

c* legacy input file format parsing code
c	include 'oldread_strcy.for'

*********
	if (iReadProject('inlacw',16).ne.0) then	!``````````

	write(*,*) 'Job Name: "',
     +		Project.sJobTitle(1:lentrim(Project.sJobTitle)),'"'
	write(16,*) 'Job Name: "',
     +		Project.sJobTitle(1:lentrim(Project.sJobTitle)),'"'
	write(*,*) 'Method clone: ',TCalcVar(Project.iUseCoredTube+1)
	write(16,*) 'Method clone: ',TCalcVar(Project.iUseCoredTube+1)

c	TJNM=Project.sJobTitle		! CHARACTER*80 TJNM
	A=CellStruct.A
	B=CellStruct.B
	C=CellStruct.C
	V=CellStruct.V
	ChiralityAngle=CellStruct.ChiralityAngle
	GMax=Project.EnergyCutOff
	Ncell=FP.iNcell
        iUseCoredTube=Project.iUseCoredTube

* re-normalise atom location coords:
	if (Project.iUseCoredTube.eq.2) A=B	! dummy: A=0 in that case
	do it=1,nt
	  do i=1,nts(it)
	    RSX(i,it)=RSX(i,it)/A
	    RSY(i,it)=RSY(i,it)/A
	    RSZ(i,it)=RSZ(i,it)/A
	  enddo
	enddo
	BSZ(3)=C/A                                        ! ???
	if (Project.iUseCoredTube.eq.2) A=0.0d0	! restore zero A

	endif						!``````````

	if (A.lt.0.0d0.or.B.lt.0.0d0.or.C.lt.0.0d0) then
	  stop '*ERR* One or more of A,B,C are negative'
	endif

	zf_params.A=A
	zf_params.B=B
	zf_params.V=V

	write(16,'(8x,''A = '',f8.4,'' a.u.'')') A
	write(16,'(8x,''B = '',f8.4,'' a.u.'')') B
	write(16,'(8x,''C = '',f8.4,'' a.u.'')') C
	write(16,'(8x,''V = '',f8.4,'' a.u.'')') V
	write(16,'(8x,''ä = '',f8.4,'' rad'')') ChiralityAngle
	write(16,*) 'EnergyCutOff:',GMax*evs,' eV'
	write(16,*) 'AddPCompanion:',Project.bAddPCompanion

	CALL SXCPOT(Co)

* The atomic electron density, AMTC, will be changed due to
* the overlapping of charge density (OCD action).
* But we will need the original AMTC to build full potential. Make a copy:
	DO IT=1,NT
	DO I=1,NGRID
	  AMTC0(I,IT)=AMTC(I,IT)
	ENDDO
	ENDDO
*
d	do it=1,nt
d	write(16,*) 'Atom type',it
d	do i=2,NGRID
d	write(16,*) rad(i,it),
d     +	ATINT_(amtc,rad,i,it,step(it),NGRID,NCUT,0,1.0d0)
d	enddo
d	enddo
d	stop

*! added dummy line in SP 3.1:
      if (iUseCoredTube.eq.2) A=B
*
*----------------
c      DO IT=1,NT
c         TMP=DIADL1_(AMTC,RAD,NGRID,IT,STEP(IT),NGRID,NCUT,1,1.D0)
c         CALL POISON(AMTC,IT,NGRID,NCUT,NGRID,WPOT1,TMP,RSTART(IT),
c     +               STEP(IT),Ecache,Fcache)
c         DO I=1,JRIS(IT)
c           RHO=AMTC(I,IT)/RAD(I,IT)/RAD(I,IT)
c           CALL XCPOT(.5D0*RHO,.5D0*RHO,RHO,CORX,TMP,DUMMY,Co)
cc       -- Coulomb part + Exchange-correlation part of potential: --
c           WPOT(I,IT)=-2.D0*ZMAIN(IT)/RAD(I,IT)+WPOT1(I)
c     .                 /DSQRT(RAD(I,IT)) + CORX
c         END DO
c      END DO
*----------------
*		----------------------------------
*		Overlap Charge Density model (OCD)
*		----------------------------------
* Charge density of atoms in neighbour cells overlaps with
* that of atoms in the current cell. We add the extra density
* and take its spherically symmetrical average. AMTC is being updated;
* ATINT_(AMTC,RAD,...) is becoming no longer equal to the nuclear charge.
*
* NOTE: AMTC is Rho(r)*r*r, with ``Rho'' standing for charge density.
*
      DO 9 IT=1,NT
        IF(NTS(IT).GT.NQG) THEN
           WRITE(16,145) IT,TANAME(IT),NTS(IT)
           WRITE(*,145) IT,TANAME(IT),NTS(IT)
           STOP
        ENDIF
        DO IQ=1,NTS(IT)
          FLOWX=RSX(IQ,IT)-RSX(1,IT)
          FLOWY=RSY(IQ,IT)-RSY(1,IT)
          FLOWZ=RSZ(IQ,IT)-RSZ(1,IT)
          CALL VECREL(RMAX,FLOWX,FLOWY,FLOWZ,BSX,BSY,BSZ,DR,NSR,NSH)
          JQJ(IT,IT,IQ)=MIN0(NSH,10)
          DO J=1,NSH
            XNAN(J,IT,IT,IQ)=DBLE(NSR(J))
            ADN(J,IT,IT,IQ)=A*DR(J)
          END DO
        END DO
        DO 9 JT=1,NT
          IF(IT.NE.JT) THEN
            DO IQ=1,NTS(JT)
              FLOWX=RSX(IQ,JT)-RSX(1,IT)
              FLOWY=RSY(IQ,JT)-RSY(1,IT)
              FLOWZ=RSZ(IQ,JT)-RSZ(1,IT)
              CALL VECREL(RMAX,FLOWX,FLOWY,FLOWZ,BSX,BSY,BSZ,DR,NSR,NSH)
              JQJ(JT,IT,IQ)=MIN0(NSH,10)
              DO J=1,JQJ(JT,IT,IQ)
                XNAN(J,JT,IT,IQ)=DBLE(NSR(J))
                ADN(J,JT,IT,IQ)=A*DR(J)
              END DO
              WRITE(16,147) IT,TANAME(IT),JT,TANAME(JT),IQ,
     +                           FLOWX,FLOWY,FLOWZ,NSH
              WRITE(*,147) IT,TANAME(IT),JT,TANAME(JT),IQ,
     +                           FLOWX,FLOWY,FLOWZ,NSH
            END DO
  147       FORMAT(5X,60('-'),/,5X,'Atom',I3,1X,A10,
     +        ' /OCD Action from Atom',I3,1X,A10,' Site No',I3,/,
     +      5X,'Shift Vector',3f7.3,' Shells for OCD',I4,/,5X,60('-'))
          ENDIF
    9 CONTINUE

      DO IT=1,NT
        DO IQ=1,NTS(IT)
          CALL SUMAX(ADN,XNAN,AMTC,RAD,JRIS,STRJ,RSTART,STEP,IT,IT,IQ,
     .               NATOM(IT),JQJ(IT,IT,IQ),NCUT)
        END DO
        DO JT=1,NT
         IF(JT.NE.IT) THEN
           DO IQ=1,NTS(JT)
             CALL SUMAX(ADN,XNAN,AMTC,RAD,JRIS,STRJ,RSTART,STEP,JT,IT,
     .                  IQ,NATOM(JT),JQJ(JT,IT,IQ),NCUT)
           END DO
         ENDIF
        ENDDO	! JT
      ENDDO	! IT

   10 CONTINUE
*
* End of the OCD model. AMTC has been adjusted to reflect
* the interference of neighbour cells' atoms.
*-------------------------------------------------------------------------
*
*	Obtaining spherically symmetrical atomic potential (Poisson eqn.)
*	-----------------------------------------------------------------
*
* Given atomic charge density, AMTC(r), solve the Poisson equation
*		Laplacian[ W(r) ] = -4pi*Rho(r)
* to obtain atomic potential WPOT(r).
*
* NOTES:
*    +	DIADL1_ and ATINT_(y,x,...) return the integral of y(x)dx over
*	the grid x. ATINT_ seems to be more precise (it gives more reasonable
*	value of total electronic charge of neutral atom).
*
*    +	POISON needs to know the normalisation coeff
*	for the integrand function. The same number of grid points
*	should be used as the 3rd argument of DIADL1_/ATINT_ and
*	as the 5th argument of POISON.
*
	WPOT(:,:)=0.0d0
	DO IT=1,NT
	  WPOT1(:)=0.0d0

cc	  TMP=DIADL1_(AMTC,RAD,NGRID,IT,STEP(IT),NGRID,NCUT,1,1.D0)
c	  TMP=ATINT_(AMTC,RAD,NGRID,IT,STEP(IT),NGRID,NCUT,0,1.D0)
c	  CALL POISON(AMTC,IT,NGRID,NCUT,NGRID,WPOT1,TMP,RSTART(IT),
c     +			STEP(IT),Ecache,Fcache)

c	  TMP=DIADL1_(AMTC,RAD,JRIS(IT),IT,STEP(IT),NGRID,NCUT,1,1.D0)
	  TMP=ATINT_(AMTC,RAD,JRIS(IT),IT,STEP(IT),NGRID,NCUT,0,1.D0)
	  CALL POISON(AMTC,IT,NGRID,NCUT,JRIS(IT),WPOT1,TMP,RSTART(IT),
     +			STEP(IT),Ecache,Fcache)

d	  write(24,*) it,TMP
	  DO I=1,NGRID	!JRIS(IT)
	    RHO=AMTC(I,IT)/RAD(I,IT)/RAD(I,IT)
	    CALL XCPOT(.5D0*RHO,.5D0*RHO,RHO,CORX,TMP,DUMMY,Co)
c	-- Coulomb part + Exchange-correlation part of potential: --
	    WPOT(I,IT)=-2.D0*ZMAIN(IT)/RAD(I,IT)+
     .			WPOT1(I)/DSQRT(RAD(I,IT)) + CORX
	  END DO
	END DO
	open(121, file='WPOT.txt')
	open (432,file='buff.dat',form='unformatted')
	write(432) WPOT
	write(121,*) WPOT
	close(432)
	close(121)

d	do it=1,nt
d	open(26,file='wpot'//char(48+it)//'.lst',form='formatted')
d	do i=1,NGRID
d	   write(26,*) rad(i,it),amtc(i,it),wpot(i,it)
d	enddo
d	close(26)
d	enddo
d	stop
*
* Done with the Poisson equation; got atomic potential, WPOT(r).
*-------------------------------------------------------------------------

	GO TO 99          ! When error in reading of outatm.str...
  201	WRITE(*,105) IT
  105	FORMAT(/5X,'*ERR** in reading atom data for atom type ',I2)
	GO TO 98
  202	WRITE(*,106) IT
  106	FORMAT(5X,'*ERR* END in reading atom data for atom type ',I3)
   98	STOP '---> Reading of OUTATM.STR.'
*	-- Calculating of indexes M,N,P --
   99	I=1
*
*! added dummy line in SP 3.1:
	if (iUseCoredTube.eq.2) A=0.0d0
*
*! begin fragment SP1.3
* Try to read an explicit basis from SelectB.dat.
* ADD FILE TO PROJECT: SelectB.for
	countx=iReadExplicitBasis(INDX,3,MUL*MAXMNP,m,n,np)
	if (countx.le.1) then
	  countx=0
	  m=MMAX-1
	  n=NMAX-1
	else

	  write(16,*) '*WRN* BASIS OVERRIDE! The basis has been set',
     +	' explicitly from the SelectB.dat file. Rename or remove it',
     +	' to use the standard energy-cutoff method'
	  write(16,*) countx,' functions will be used'

	  GMAX=1000.d0/evs       ! or another big value;
* we already know the basis function set without this energy-cutoff parameter
	endif
*! end fragment

*! begin fragment
	if (iUseCoredTube.eq.1.or.iUseCoredTube.eq.3.or.
	*iUseCoredTube.eq.4) then

	  print *,'Searching for Kappas in range',fk_XSTART,dsqrt(GMAX)
	  call FindKappa(AllKappas,MMAX,NMAX,fk_DX,fk_XSTART,dsqrt(GMAX),
     .	m+1,n+1,NAct)
*	^^^^^^^changed 'MMAX,NMAX'-->'m,n' in SP1.3
*! begin fragment: override kappas with those of SP0.0
d	  tmp=(B/A)/0.03
d	  if (tmp.lt.1.0d0) then
d	    write(16,*) 'overriding kappas:',tmp
d	    do m=1,MMAX
d	    do n=1,Nact(m)
d	      call bess_nul(n,m-1,RES,MMAX)
d	      tmp1=RES(n)/A       !AllKappas(m,n)*tmp+(1.0d0-tmp)*RES(n)/A
d	      write(16,*) m,n,AllKappas(m,n),tmp1
d	      AllKappas(m,n)=tmp1
d	    enddo
d	    enddo
d	  endif
*! end fragment

	call FindBessAtAB(A,B,zf_params.eps,
     .		AllKappas,MMAX,NMAX,NAct,abJ,abJD,abY,abYD)

*! commented fragment SP 5.07: Kappas are saved into Outstrcy.dat
c	! save for future use:
c	open(33,file='allkappa.dat',form='unformatted')
c	write(33) A,B,C
c	write(33) (NAct(m),m=1,MMAX)
c	write(33) ((AllKappas(m,n),n=1,NAct(m)),m=1,MMAX)
c	write(33) (((abJ(m,n,i),abJD(m,n,i),abY(m,n,i),abYD(m,n,i),
c     +			i=1,2),n=1,NAct(m)),m=1,MMAX)
c	close(33)
*! end fragment SP 5.07

*! begin old fragment SP <5.0
d	open(33,file='allkappa.dat.old',form='unformatted')
d	write(33) a,b,c
d	call SaveArray2D_1(33,AllKappas,MMAX,NMAX,NAct)
d	close(33)
*! end old fragment SP <5.0

d	do m=1,MMAX
d	do n=1,NAct(m)
d	  write(17,*) m,n,AllKappas(m,n)
d	  write(17,*) (abJ(m,n,i),i=1,2)
d	  write(17,*) (abJD(m,n,i),i=1,2)
d	  write(17,*) (abY(m,n,i),i=1,2)
d	  write(17,*) (abYD(m,n,i),i=1,2)
d	enddo
d	enddo

d	stop '*TEST* Stopped after FindKappa()'
	endif
*! end fragment

*! begin fragment SP1.3
	if (countx.ne.0) then     ! basis override (with an explicit one)

	do i=1,countx
	  M=INDX(1,I)
	  N=INDX(2,I)
	  NP=INDX(3,I)
*! old fragment SP 3.01
c	  if (iUseCoredTube.ne.0) then
c	    KPNM=AllKappas(iabs(M)+1,N)
c	  else
c	    CALL BESS_NUL(N,IABS(M),RES,MAXMNP)
c	    KPNM=RES(N)/A
c	  endif
*! new fragment SP 3.1
c Bess_nul(N,M,dest,dest_size,bess_family)
c searches for N roots of a Bessel function J (bess_family==1)
c or Y (bess_family==2) of order M and stores those roots into
c the dest array.
	  if (iUseCoredTube.eq.0) then
	    CALL BESS_NUL(N,IABS(M),RES,MAXMNP,1)
	    KPNM=RES(N)/A
	  else if (iUseCoredTube.eq.1.or.iUseCoredTube.eq.3.or.
	*    iUseCoredTube.eq.4) then
	    KPNM=AllKappas(iabs(M)+1,N)
	  else if (iUseCoredTube.eq.2) then
	    CALL BESS_NUL(N,IABS(M),RES,MAXMNP,2)
	    KPNM=RES(N)/B
	  endif
*! end fragment
*! old fragment SP1x
c	  VKP=2.D0*PI*DBLE(NP)/C
*! new fragment SP2.0
	  VKP=(2.D0*PI*DBLE(NP)-dble(M)*ChiralityAngle)/C
c	  VKP=2.D0*PI*(DBLE(NP)-dble(M)*ChiralityAngle)/C
*! end fragment SP2.0
	  ENER(I)=KPNM*KPNM+VKP*VKP
	  ENKP(iabs(M)+1,N)=KPNM
	enddo

	! sort by energy and prepare for saving into Outstrcy.da0:
	CALL INDEXX(COUNTx,MUL*MAXMNP,ENER,IND)    ! quick sort by ENER()
	do i=1,countx
	  INDXX(1,i)=INDX(1,IND(I))
	  INDXX(2,i)=INDX(2,IND(I))
	  INDXX(3,i)=INDX(3,IND(I))
	  INDXX(4,i)=IND(I)  ! Save index for using in output of energies
	enddo

	else      ! <-------- countx.eq.0 (no basis override)
*! end fragment

	COUNT=1
	I=1
	DO M=0, MAXMNP
	DO N=1, MAXMNP

*! old fragment
c	  CALL BESS_NUL(N,IABS(M),RES,MAXMNP)
c	  KPNM=RES(N)/A
*! new fragment SP 3.01
c	  if (iUseCoredTube.ne.0) then
c	  KPNM=AllKappas(M+1,N)
cd	  print *,'sel m,n,x:',m,n,kpnm,kpnm*kpnm
c	  else
c	    CALL BESS_NUL(N,IABS(M),RES,MAXMNP)
c	    KPNM=RES(N)/A
c	  endif
*! new fragment SP 3.1
	  if (iUseCoredTube.eq.0) then
	    CALL BESS_NUL(N,IABS(M),RES,MAXMNP,1)
	    KPNM=RES(N)/A
	  else if (iUseCoredTube.eq.1.or.iUseCoredTube.eq.3.or.
	*           iUseCoredTube.eq.4) then
	    KPNM=AllKappas(iabs(M)+1,N)
d	    print *,'sel m,n,x:',m,n,kpnm,kpnm*kpnm
	  else if (iUseCoredTube.eq.2) then
	    CALL BESS_NUL(N,IABS(M),RES,MAXMNP,2)
	    KPNM=RES(N)/B
	  endif
*! end fragment

	  DO NP=0, MAXMNP
*! old fragment SP1x
c	    VKP=2.D0*PI*DBLE(NP)/C
*! new fragment SP2.0
	    VKP=(2.D0*PI*DBLE(NP)-dble(M)*ChiralityAngle)/C
c	    VKP=2.D0*PI*(DBLE(NP)-dble(M)*ChiralityAngle)/C
*! end fragment SP2.0
	    ENER(I)=KPNM*KPNM+VKP*VKP
          ENKP(iabs(M)+1,N)=KPNM

c         take energy < cut off with all P
	    IF(ENER(I).LT.GMAX.OR.VKP.GE.0) THEN
	      IF(ENER(I).GT.GMAX) EXIT
	      INDX(1,I)=M
	      INDX(2,I)=N
	      INDX(3,I)=NP
d	      write(16,'(1x,a,4i5,f20.10)') '>',I,M,N,NP,ENER(I)
	      I=I+1
	      IF(I.GT.MUL*MAXMNP) THEN
	        WRITE(*,112) MUL*MAXMNP, '-->INDEX.'
	        STOP
	      ENDIF
          ENDIF
	  END DO
c   -- Here we test part of ENER(i), exclude Kp^2, for M and N only. --
	  IF(KPNM*KPNM.GT.GMAX) EXIT
	END DO
c   -- Here we test part of ENER(i), exclude Kp^2 and for N=1 for M. --
	IF(KPNM*KPNM.GT.GMAX.AND.N.EQ.1) THEN
	  COUNT=I-1
	  EXIT
	ENDIF
	END DO
*-------------------------------------------------------------------------
	CALL INDEXX(COUNT,MUL*MAXMNP,ENER,IND)    ! quick sort by ENER()
*-------------------------------------------------------------------------
d	do i=1,COUNT
d	  write(16,'(1x,a,4i5,f20.10)') '?',
d     +	  IND(I),(INDX(J,IND(I)),J=1,3),ENER(IND(I))
d	enddo

	COUNTX=1
	DO I=1,COUNT         !  here we add '-' - values to basis set.
	  INDXX(1,COUNTX)=INDX(1,IND(I))
	  INDXX(2,COUNTX)=INDX(2,IND(I))
	  INDXX(3,COUNTX)=INDX(3,IND(I))
	  INDXX(4,COUNTX)=IND(I)  ! Save index for using in output of energies
d	  write(16,*) '+=',countx,(indxx(j,countx),j=1,4),ener(ind(i))
	  COUNTX=COUNTX+1
	  IF(COUNTX.GT.MUL*MAXMNP) THEN
	    WRITE(*,112) MUL*MAXMNP,'-->INDEXX{1}.'
	  STOP
        ENDIF
        IF(INDX(1,IND(I)).GT.0) THEN
          INDXX(1,COUNTX)=-INDX(1,IND(I))       ! for M < 0
          INDXX(2,COUNTX)=INDX(2,IND(I))
          INDXX(3,COUNTX)=INDX(3,IND(I))
          INDXX(4,COUNTX)=IND(I)
          IF(COUNTX.GT.MUL*MAXMNP) THEN
             WRITE(*,112) MUL*MAXMNP,'-->INDEXX{2}.'
             STOP
          ENDIF
          COUNTX=COUNTX+1
        ENDIF
        IF(INDX(3,IND(I)).GT.0) THEN
cr      if ((ener(ind(i))+4.0d0*
cr     . dble(1+2*INDX(3,IND(I)))*PI/C*PI/C).le.GMAX) then
          INDXX(1,COUNTX)=INDX(1,IND(I))        ! for NP < 0
          INDXX(2,COUNTX)=INDX(2,IND(I))
          INDXX(3,COUNTX)=-INDX(3,IND(I))
          INDXX(4,COUNTX)=IND(I)
          COUNTX=COUNTX+1
          IF(COUNTX.GT.MUL*MAXMNP) THEN
             WRITE(*,112) MUL*MAXMNP,'-->INDEXX{3}.'
             STOP
          ENDIF
cr      endif
        ENDIF
        IF(INDX(1,IND(I)).GT.0 .AND. INDX(3,IND(I)).GT.0) THEN
cr      if ((ener(ind(i))+4.0d0*
cr     . dble(1+2*INDX(3,IND(I)))*PI/C*PI/C).le.GMAX) then
          INDXX(1,COUNTX)=-INDX(1,IND(I))       ! for M < 0
          INDXX(2,COUNTX)=INDX(2,IND(I))
          INDXX(3,COUNTX)=-INDX(3,IND(I))       ! ...and for NP < 0
          INDXX(4,COUNTX)=IND(I)
          COUNTX=COUNTX+1
          IF(COUNTX.GT.MUL*MAXMNP) THEN
             WRITE(*,112) MUL*MAXMNP,'-->INDEXX{4}.'
             STOP
          ENDIF
cr      endif
        ENDIF
      END DO
      COUNTX=COUNTX-1

*! begin fragment SP1.3
	endif     ! <-------- countx.eq.0 (no basis override)
*! end fragment

*! begin fragment SP2.10 form SelectB.da~ file
	open(33,FILE='SelectB.da~',FORM='FORMATTED')
	write(33,212)
  212	format('# Strcy.exe has generated this file for you.',/,
     +	'# It contains the basis used in the calculation.',/,
     +	'# Rename it to SelectB.dat to manually define basis.',/,
     +	'# If (m,n,p+1) companions are added, they are not shown here.')
*	Start StrCy, then rename SelectB.da~ to SelectB.dat
*	to modify the real basis manually
	do i=1,countx
	  write(33,'(3i8)') (indxx(j,i),j=1,3)
	enddo
	close(33)
*! end fragment

*! begin fragment SP1.3 add (m,n,p+1) companions
	! countx will rise
	count=countx

	! save energies of basis functions:
	do i=1,count
	  indx(1,i)=indxx(1,i)
	  indx(2,i)=indxx(2,i)
	  indx(3,i)=indxx(3,i)
	  EnerSort(i)=ener(indxx(4,i))
d	  write(16,*) 'bwf',i,(indx(j,i),j=1,3),EnerSort(i)
	enddo

	! check for (m,n,p+1) presence in the basis:
	do i=1,count
	  m=indx(1,i)
	  n=indx(2,i)
	  np=indx(3,i)
	  ener(i)=EnerSort(i)
	  it=0
	  do j=1,count
	    if (i.ne.j.and.
     .	m.eq.indx(1,j).and.n.eq.indx(2,j).and.np+1.eq.indx(3,j)) it=1
	  enddo
	  if (it.eq.0.and.Project.bAddPCompanion) then	! add (m,n,p+1)
	    COUNTX=COUNTX+1
	    IF(COUNTX.GT.MUL*MAXMNP) THEN
	      WRITE(*,112) MUL*MAXMNP,'-->INDEXX{4}.'
	      STOP
	    ENDIF
	    indx(1,countx)=m
	    indx(2,countx)=n
	    indx(3,countx)=np+1
	    VKP=2.D0*PI/C
*! old fragment SP1x
c	    ener(countx)=ener(i)+VKP*VKP*dble(2*NP+1)
*! new fragment SP2.0
	    ener(countx)=
     .	ener(i)+VKP*VKP*(dble(2*NP+1)-dble(M)*ChiralityAngle/PI)
c     .	ener(i)+VKP*VKP*(dble(2*NP+1)-dble(M)*ChiralityAngle*2)
*! end fragment SP2.0
d	    countx=countx-1  !!!CHEATING!!!
d	    write(16,*) 'Added ',m,n,np+1,': E = ',ener(countx)
	  endif
	enddo

	! re-sort the basis:
	CALL INDEXX(COUNTx,MUL*MAXMNP,ENER,IND)    ! quick sort by ENER()

	do i=1,countx
	  indxx(1,i)=indx(1,ind(i))
	  indxx(2,i)=indx(2,ind(i))
	  indxx(3,i)=indx(3,ind(i))
	  indxx(4,i)=ind(i)
d	  write(16,*) 'se',i,(indxx(j,i),j=1,3),ener(indxx(4,i)),ind(i)
	enddo
*! end fragment


	WRITE (*,110) COUNTX
	WRITE (16,110) COUNTX
  110	FORMAT(/5X,'***** Basis has ',I5,' function(s). *****')
	N=0
	if(Project.iUseCoredTube.eq.3)then
  901	  FORMAT(3X,A34,3X,A34)
  902	  FORMAT(1X,I3,3X,2(I3,3X),I3,3X,F8.3,6X,
	*    F8.4,1X,F8.4,1X,F8.4,1X,F8.4)
        WRITE(16,901) 'No    M     N     P    Energy, eV.',
	*    '   Knm      Cj       Cy     Ccr*Km'
	else if(Project.iUseCoredTube.eq.4)then
        WRITE(16,901) 'No    M     N     P    Energy, eV.',
	*    '   Knm      Cj       Cy     Ccri*Km'
	else
	  WRITE(16,'(//,3X,A34,/)') 'No    M     N     P    Energy, eV.'
      end if
	DO I=1, COUNTX
c	  IF(INDXX(4,I).NE.N) THEN
	  IF(ENER(INDXX(4,I)).NE.ENER(N)) THEN
	    if(Project.iUseCoredTube.eq.3.or.
	*       Project.iUseCoredTube.eq.4)then
            resKmn = ENKP(iabs(INDXX(1,I))+1,INDXX(2,I))
d	      if(I.eq.184)then
d              resCJ = CalcCJ(INDXX(1,I),A,B,resKmn,V)
d	      else
q              resCJ = CalcCJ(INDXX(1,I),A,B,resKmn,V)
d	      endif

q            resCY = CalcCY(INDXX(1,I),resCJ,B,resKmn)
q            resCCR = CalcCCR(INDXX(1,I),resCJ,A,B,resKmn,V)
            WRITE(16,902) I,(INDXX(J,I),J=1,3),ENER(INDXX(4,I))*EVS,
	*        resKmn, resCJ, resCY, resCCR
	    else
	      WRITE(16,211) I,(INDXX(J,I),J=1,3),ENER(INDXX(4,I))*EVS
	    end if
	    N=INDXX(4,I)
	  ELSE
	    WRITE(16,210) I,(INDXX(J,I),J=1,3)
	  ENDIF
	END DO
	WRITE(16,'(/5x,A16)') '<<< Complete >>>'
	CLOSE(16)                         ! End of log-file.


*
*-------------------------------------------------------------------------
*! begin fragment SP5.0: save all Strcy results to file:
	Basis.INDXX=>INDXX
	Basis.COUNTX=COUNTX
	open(33,file='outstrcy.dat',form='unformatted')
	if (iWriteStrResult(33,Basis).eq.0)
     +		stop 'Error saving Strcy results to outstrcy.dat'
	close(33)
*! end fragment SP5.0
*-------------------------------------------------------------------------
*
	GO TO 50
*  29	FORMAT(/,<K>X,'<<< ',A<J-I>,' >>>',//)
*! old fragment:   30	FORMAT(5X,A18,':',F7.3,','/5X,A14,':',F7.3,'.')
*! new fragment:
   30	FORMAT(
     +	5x,a18,',',1x,'A:',f7.3,' a.u.,',/,
     +	5x,a15,',',4x,'B:',f7.3,' a.u.,',/,
     +	5x,a11,',',8x,'C:',f7.3,' a.u.,',/,
     +	5x,a15,',',4x,'ä:',f14.10,' rad')
*! end fragment
   32	FORMAT(/1X,'Atom type ',I2,': ',A12,\5X,'Nuclear charge ',I3/)
   33	FORMAT(1X,I3,' atom(s) found.')
   43	FORMAT(7X,A14,':',F10.5)
   44	FORMAT(1X,A2,'-',A21,' of potential inside MT''s selected.')
   45	STOP '*ERR* Can''t create Outstrcy.dat.'
   47	STOP '*ERR* Can''t write to the data interchange file.'
  112	FORMAT(/'*ERR* Count of elements > MAX.=',I6,1X,A14)
  145	FORMAT(5X,'*ERR* Atom Type',I3,1X,A12,' is too LARGE',' Nts',I5)
  210	FORMAT(1X,I3,3X,2(I3,3X),I3)
  211	FORMAT(1X,I3,3X,2(I3,3X),I3,3X,F8.3)
   50	WRITE(*,'(A6//,10X,A32)') ' done.','<<< STRCY module complete >>>'

	END
*
*
	SUBROUTINE VECREL(RMAX,FLOWX,FLOWY,FLOWZ,BSX,BSY,BSZ,DR,NSR,NSH)
*-------------------------------------------------------*
*		Real Space Vectors Generation		*
*-------------------------------------------------------*
	INTEGER NSR(*),NSH,NSHL,NUMR,N,N1,NR
	DOUBLE PRECISION RMAX,FLOWX,FLOWY,FLOWZ,BSX(*),BSY(*),BSZ(*),
     +		DR(*),CSX(300),CSY(300),CSZ(300),D(300),RA,DX,DA,DB,
     +		SX,SY,SZ
*
	WRITE(*,100)
  100	FORMAT(5X,/,5X,60('-'),
     *	/,5X,'Result from VECREL for Real Space Vectors',
     *	/,13X,'No',5X,'Sx',8X,'Sy',8X,'Sz',8X,'D',/)
	RA=DSQRT((BSX(1)+BSX(2)+BSX(3))**2+(BSY(1)+BSY(2)+BSY(3))**2+
     .		(BSZ(1)+BSZ(2)+BSZ(3))**2)/2.D0			! ???
	RA=RA+RMAX+DSQRT(FLOWX**2+FLOWY**2+FLOWZ**2)
	DX=0.D0
	DO N=1,3
	  DX=DX+(BSX(N)**2+BSY(N)**2+BSZ(N)**2)
	END DO
	NUMR=2*(IDINT(RA/DSQRT(DX/3.D0))+1)+1
	NR=0
	DO 1 N=1,NUMR
	  SX=FLOWX
	  SY=FLOWY
	  SZ=DBLE(N-NUMR/2-1)*BSZ(3)+FLOWZ
	  DX=DSQRT(SX*SX+SY*SY+SZ*SZ)
	  IF(DX.GT.RMAX) GO TO 1
	  NR=NR+1
	  IF(NR.GT.300) THEN
	    WRITE(*,101) NR
	    STOP
	  ENDIF
	  D(NR)=DX
	  CSX(NR)=SX
	  CSY(NR)=SY
	  CSZ(NR)=SZ
    1	CONTINUE
*
*	Sort Vectors in Order of Increasing D
*
	N1=1
	DA=D(1)
	DO N=2,NR
	  IF(D(N).LT.DA) DA=D(N)
	END DO
	DR(1)=DA
*
	NSH=0
	NSHL=0
	DO 3 K=1,NR
	  DX=D(1)
	  N1=1
	  DO N=2,NR
	   IF(D(N).LT.DX) THEN
	     DX=D(N)
	     N1=N
	   ENDIF
	  END DO
*
	  ASX=CSX(N1)
	  ASY=CSY(N1)
	  ASZ=CSZ(N1)
	  DB=D(N1)
*
	  IF(DB.LT.DA+2.5D-05) THEN		! ??? Eps?
	    NSHL=NSHL+1
	    WRITE(*,102) K,ASX,ASY,ASZ,DB
	  ELSE
	    NSH=NSH+1
	    IF(NSH.LE.10) THEN
	      WRITE(*,103) NSH,NSHL,DR(NSH)
	      IF(NSH.LT.10) DR(NSH+1)=DB
	      WRITE(*,102) K,ASX,ASY,ASZ,DB
	      NSR(NSH)=NSHL
	  ENDIF
  102	  FORMAT(10X,I5,4F10.6)
  103	  FORMAT(2X,'Shell Number',I5,' with',I5,' Points, Dist',F10.6/)
	  NSHL=1
	  DA=DB
	ENDIF
	D(N1)=1.D+10
    3	CONTINUE
	IF(NSH.GT.10) THEN
	   WRITE(*,104) NSH,NSHL,DA
	   NSH=10
	ENDIF
  101	FORMAT(1X,'*ERR* Internal error: out of range of dimensions.',I5)
  104	FORMAT(5X,'Max.Shell Number',I5,' with',I5,' Points, Dist',F10.6)
	RETURN
	END
*
	SUBROUTINE SUMAX(ADN,XNAN,AMTC,RAD,JRIS,STRJRI,RSTART,STEP,
     .			 JT,IT,IQ,JTOP,JSUM,NCUT)
c--------------------------------------------------------------------c
c            Make Out OCD model for Electron Density                 c
c--------------------------------------------------------------------c
	INTEGER JRIS(*),JT,IT,IQ,JTOP,JSUM,JRI,JRIJ,NCUT,I,J,JA,
     +		JTOP1,JBOT
	DOUBLE PRECISION FLOW(300),STOR(300),XTOP,XTOP1,XBOT,XBOT1,
     +		RAD(300,NCUT),ADN(10,NCUT,NCUT,*),AMTC(300,NCUT),
     +		XNAN(10,NCUT,NCUT,*),RSTART(*),STEP(*),STRJRI(*),
     +		ADFLOW,EpsLoc,XINT,XINT1,XINT2,XSTOP
	DATA EpsLoc/1.D-3/
*
	JRI=JRIS(IT)
	JRIJ=JRIS(JT)
	DO I=JRIJ+1,JTOP
	  FLOW(I)=AMTC(I,JT)
	END DO
	FLOW(JRIJ)=STRJRI(JT)
	DO I=1,JRI
	  STOR(I)=0.D0
	END DO
	XSTOP=DLOG(RAD(JTOP,JT))
*
	DO JA=1,JSUM
	  ADFLOW=ADN(JA,JT,IT,IQ)
	  IF(ADFLOW.GT.EpsLoc) THEN
	    DO 5 I=1,JRI
	      XTOP=DLOG(ADFLOW+RAD(I,IT))
	      XBOT=DLOG(ADFLOW-RAD(I,IT))
	      IF( XBOT.LT.XSTOP ) THEN
		JTOP1=IDINT((XTOP-RSTART(JT))/STEP(JT))+2
		XTOP1=RSTART(JT)+STEP(JT)*DBLE(JTOP1-1)
		JBOT=IDINT((XBOT-RSTART(JT))/STEP(JT))+1
		XBOT1=RSTART(JT)+STEP(JT)*DBLE(JBOT-1)
		IF(JTOP1.LE.JTOP) THEN
		  XINT1=.5d0*(XTOP1-XTOP)*(2.d0*FLOW(JTOP1)+
     .		    (FLOW(JTOP1-1)-FLOW(JTOP1))*(XTOP1-XTOP)/STEP(JT))
		ELSE
		  JTOP1=JTOP
		  XINT1=0.D0
		ENDIF
		IF(JBOT.GE.JRIJ) THEN
		  XINT2=.5d0*(XBOT-XBOT1)*(2.d0*FLOW(JBOT)+
     .		    (XBOT-XBOT1)/STEP(JT)*(FLOW(JBOT+1)-FLOW(JBOT)))
		ELSE
		  WRITE(*,7) JBOT,JRIJ,JT,IT
   7		  FORMAT(5X,'*** WARNING SUMAX JBOT=',I3,' JRIJ=',I3,
     +			' JT=',I1,' IT=',I1)
		  XINT2=0.D0
		ENDIF
		IF( JTOP1-JBOT.LT.0 ) THEN
		  WRITE(*,8) JBOT,JTOP1,JT,IT
   8		  FORMAT(5X,'*** SUMAX AUTOSTOP JBOT=',I3,' JTOP1=',I3,
     +			' JT=',I1,' IT=',I1)
		  STOP
		ELSE
		  XINT=-(XINT1+XINT2)/STEP(JT)*2.d0
		  DO J=JBOT,JTOP1-1
		    XINT=XINT+FLOW(J)+FLOW(J+1)
		  END DO
		  XINT=XINT*.5d0*STEP(JT)
		ENDIF
		STOR(I)=STOR(I)+XINT*RAD(I,IT)*.5d0*XNAN(JA,JT,IT,IQ)/ADFLOW
	      ENDIF
   5	    CONTINUE
	  ENDIF
	END DO
	DO I=1,JRI
	  AMTC(I,IT)=AMTC(I,IT)+STOR(I)
	END DO
	RETURN
	END
