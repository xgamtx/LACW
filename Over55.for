* --------------------------------------------------------------------- *
*		CYLINDER LAPW MAIN PROGRAM MODULE			                          *
*		(c) 1998 Kepp O.M., D'yachkov P.N.			                        *
* --------------------------------------------------------------------- *
*	Tubular Full-potential LACW - Service Pack 5.10  	                *
*	(c) 1999-2001-? by Dmitry V. Kirin & Pavel N. D'yachkov		        *
* --------------------------------------------------------------------- *
*	Tubular Full-potential LACW - Service Pack 6.0                   	*
*	(c) 2003-2004-? by Dmitry V. Makaev & Pavel N. D'yachkov		        *
* --------------------------------------------------------------------- *
* Version history:							*
*	Service Pack 0		Original Cylindrical LACW		*
*	Service Pack 2		Tubular LACW				*
*	Service Pack 3		Inverted/Nanopore LACW (untested)	*
*	Service Pack 4		Full-Potential Extension		*
*	Service Pack 5		Unified project data format		*
* --------------------------------------------------------------------- *
*	`Don't be too proud of this technological terror		*
*	you've constructed. The ability to destroy a planet		*
*	is insignificant next to the power of the Force.'		*
*			--Darth Vader, `Star Wars' episode IV		*
* --------------------------------------------------------------------- *

	IMPLICIT DOUBLE PRECISION(A-H,O-Z)
c	KPNM(1) - Kappa(N,M), KPNM(2) - Kappa(N1,M1)
c	VKP(1) - Kp(P), VKP(2) - Kp(P'), VK - k
c	Values of integrals I1 & I2 are:
c	FSI1(1,1) - I1(P,N,M).Re, FSI1(2,1) - I1(P,N,M).Img,
c	FSI1(1,2) - I1(P',N',M').Re, FSI1(2,2) - I1(P',N',M').Img,
c	FSI2(1,1) - I2(P,N,M).Re, FSI2(2,1) - I2(P,N,M).Img,
c	FSI2(1,2) - I2(P',N',M').Re, FSI2(2,2) - I2(P',N',M').Img,
*
	include 'atomdata.fh'
	include 'coordsys.fh'
	include 'projstrc.fh'
	include 'varsizes.fh'
	include 'stdio.fh'
*
	include 'strresult.fh'
	record /TBasis/ Basis
*
	INTEGER NProgPNT,NEF,Complete
	INTEGER ISH,JSH,CALC,NP,NP1,IIM,MOV,NPNT,
     +		ITRX,MAXBND,MAXPNT,MAXnEigv,NPROGRES,NCXPOINT,
     +		IBNDS,NTPNTS,NPNTS,NCPNT,NextPnt,NLmax,NLsav,
     +		IWORK(5*MUL*MAXMNP),IFAIL,LWORK,INFO,IDOSPTS,
     +		ISmode,I,J,II,IQ,IT,K
	PARAMETER (MAXSUM=50, MAXBND=MUL*MAXMNP, MAXPNT=200)
	PARAMETER (NPROGRES=66, Complete=1)
	INTEGER COUNTX,LI12CALC,IRO(NCUT,NQG)

	integer NL
	common NL

	COMMON /BAS/  P(LCUT,NCUT),PR(LCUT,NCUT),
     *	PE(LCUT,NCUT),PER(LCUT,NCUT),XE(LCUT,NCUT),
     *	QN(LCUT,NCUT),ELMIN(5,NCUT),ELMAX(5,NCUT),
     *	DELTA(5,NCUT),ITS(5,NCUT),NORLW,NOWAV,ISPIN,ISPIN1,ISPIN2
	COMMON / OU / P1(NGRID),P2(NGRID),Q1(NGRID),Q2(NGRID),
     *	BG(NGRID),BGD(NGRID),HMD(NGRID),RA(NGRID),IANS

	DOUBLE PRECISION CZ,A,B,Eps,RO,KPNM(2),VKP(2),TMP,TM1,TM2,TIME
	DOUBLE PRECISION EpsG,RES(MAXSUM),PHI,TMP1,TMPH,TMPIS,
     +		SH_RE,SH_IMG,SSI12(LCUT),HSI12(LCUT),S_RE,S_IMG,
     +		BDCALC(LCUT*2),BECALC(LCUT*2),BFCALC(LCUT*2),
     +		BJCALC(2,LCUT*2),BYCALC(2,LCUT*2),
     +		V_out,VAtShft,VL,VU,VK,TIME_of_all,
     +		BSmCALC(3+NCUT*LCUT,NQG), k_start


	double precision, allocatable:: EBND(:,:),D(:),Diag(:),RWORK(:)
	double complex, allocatable:: WORK(:)
c In SP 0..4, the arrays^ were static:
c	DOUBLE PRECISION EBND(MAXBND,MAXPNT),
c     +		D(MUL*MAXMNP),Diag(MUL*MAXMNP),RWORK(7*MUL*MAXMNP)
c	PARAMETER (LWORK=64*MUL*MAXMNP)
c	DOUBLE COMPLEX WORK(LWORK)


	DOUBLE COMPLEX, ALLOCATABLE:: H(:,:), S(:,:), PoolTmp(:,:,:)
	DOUBLE PRECISION, ALLOCATABLE:: PoolBuild(:,:), PoolRes(:,:,:)
	INTEGER, ALLOCATABLE:: IndxD(:)

	DOUBLE COMPLEX SI1(LCUT,LCUT),SI2(LCUT,LCUT),TMPc,
     +	SI1_(LCUT,LCUT),SI2_(LCUT,LCUT),TMPcc

	LOGICAL LBCALC(LCUT*2),LCALC(LCUT*2),
     +	LSmCALC(NQG),Create_I3M,TBZERO,Add_to_Saved_I3M,SaveI3M,
     +	BeforeExit,LSmode

	CHARACTER TBASE*80,TATIME*14,
     +	TCTIME*19,TJNM*80,JOB*1,UPLO*1,   !RANG*1,
     +	TNMBR*1(10),outfile*11,TDAT*4,TDOS*3,TIS*3,
     +	OutfL*10(NCUT*LCUT+2)

c	.. External Functions ..
	DOUBLE PRECISION FCT, FD, FE, FF, S_NJM, S_NYM
	DOUBLE COMPLEX SIMP_I12
	EXTERNAL FCT, FD, FE, FF, S_NJM, S_NYM, SIMP_I12
*
*! begin fragment
	double precision ChiralityAngle
	integer iUseCoredTube	! 0 (Cylindrical LACW), 1 (Tubular LACW), or 2 (Inverted LACW)

! exchange structure for zerofunc():
	structure /zerofunc_params_struct/
	  double precision A,B,V,eps
	  integer m
	end structure
	record /zerofunc_params_struct/ zf_params
	common /param__zerofunc/ zf_params

*! end
*! begin fragment SP 4.1a: FP
	integer iFP_NumPoints1,iFP_NumPoints2,iFP_NumPoints3,
     +	iFP_NumPoints123
	double precision dFP_Start1,dFP_Start2,dFP_Start3,
     +	dFP_d1,dFP_d2,dFP_d3

	integer iVfpMatrixRequest
* iVfpMatrixRequest shows if the FP matrix elements should be
* (re-)calculated during the calculation of the nearest k-point.
* After the (re-)calc, the matrix <i|Vfp|j> will be stored in fpmatrix.dat.

* (x,y,z)set are coords of real-space points where VfpMap is calculated
	double precision, allocatable:: xset(:),yset(:),zset(:)
	double precision, allocatable:: VfpMap(:)
	double complex, allocatable:: cBwfMap(:,:),VfpMEs(:)
	double complex VfpMEInt
	character, allocatable:: int_map(:)
*! end fragment SP 4.1a: FP
*
*! begin fragment SP 6.0
      double precision IntCR(MMAX, NMAX)
      integer          IntCRC(MMAX, NMAX)
*! end fragment SP 6.0


	character*4 tsee
	DATA TATIME/'One-point time'/,TCTIME/'Calculation''s time '/,
     .		JOB/'V'/,UPLO/'L'/,   !RANG/'I'/,
     .		TDOS/'dos'/,TDAT/'.dat'/,
     .		Create_I3M/.TRUE./,LSmode/.FALSE./,
     .		TNMBR/'0','1','2','3','4','5','6','7','8','9'/,
     .		ISmode/0/,TIS/'_is'/
	DATA MAXnEigv/0/,Add_to_Saved_I3M,TBZERO/2*.FALSE./
	data tsee/'SEE_'/
*
	ilog=17
	ierr=6
	iout=6
	DEBUG=7
d	open(ilog,file='lacw.log',form='formatted')

	write(*,*)'**** OVER55 **** LACW Service Pack 6.1 (based on 5.07)'
d	write(16,*)'**** OVER55 **** LACW Service Pack 6.0 (based on 5.07)'

* Read the project --- projstrc.for::iReadProject(fn,logf)
	if (iReadProject('inlacw',ierr).le.0)
     +	stop 'Error reading the project from Inlacw'
* Read the basis --- strresult.for::iReadStrResult(f,Basis)
	open(33,file='outstrcy.dat',form='unformatted')
	if (iReadStrResult(33,Basis).le.0)
     +	stop 'Outstrcy.dat is missing or corrupt. Restart Strcy.exe.'
	close(33)

* translate to old-style variables:
	COUNTX=Basis.COUNTX
	TJNM=Project.sJobTitle(1:80)	! TJNM*80
	VL=Plots.EnergyBottom
	VU=Plots.EnergyTop
	IBNDS=Project.iBands
	Eps=Project.EpsInt
	EpsG=Project.EpsCoords
	NTPNTS=Project.iPointsTotal-1
	NPNTS=Project.iPointsNow
	NL=Project.NL
	NLmax=Project.NLmax
	NLsav=Project.NLsav	!!! will change and should be updated in the Project structure
	NelCell=Project.NelCell
	TBZERO=(Project.WvStart.eq.0.0d0)
	SaveI3M=Project.SaveI3M
	ISmode=Project.bIntersphericalDos

	A=CellStruct.A
	B=CellStruct.B
	V=CellStruct.V
	CZ=CellStruct.C
	open(121, file='countx.txt')
	write(121, *) COUNTX
	write(121, *) A, B, CZ, Eps, NTPNTS, NLMAX
	do I = 1, COUNTX
	  write(121, *) Basis.INDXX(1, I), 
	+    Basis.INDXX(2, I), 
	+	Basis.INDXX(3, I)

	end do
	close(121)

	ChiralityAngle=CellStruct.ChiralityAngle
	iUseCoredTube=Project.iUseCoredTube

	zf_params.A=A
	zf_params.B=B
	zf_params.V=V
	zf_params.eps=Project.EpsBess

	HfWdth=Plots.Halfwidth
	HfWdt_EF=Plots.HalfwidthFermi
	IDOSPTS=Plots.iRenderPitch
	LSmode=.not.Plots.SeparateBandPlotFiles


	VAtShft=0.0d0
	V_out=0.0d0

*--------------------------------------------------------------------
c     clear intCRC
      do i_c = 1, MMAX
	  do j_c = 1, NMAX
	    IntCRC(i_c, j_c) = 0
	  end do
	end do
*--------------------------------------------------------------------
c	print *,'Memory alloc request [EBND]: ',
c     +		COUNTX*Project.iPointsTotal*IDBL,' bytes'
c	allocate(EBND(COUNTX,Project.iPointsTotal),stat=istat)
c	if (istat.ne.0) stop '	...failed!'

	LWORK=64*COUNTX
	print *,'Memory alloc request [D,RWORK]: ',
     +		(2+7)*COUNTX*IDBL+LWORK*ICMPLX,' bytes'
	allocate(D(COUNTX),Diag(COUNTX),RWORK(7*COUNTX),WORK(LWORK),
     +		stat=istat)
	if (istat.ne.0) stop '	...failed!'
*--------------------------------------------------------------------
* Reset the time counter:
	TIME_of_all=0.D0
	CALL TTIME(TIME_of_all,TCTIME,0)
*
*------------- Full Potential Scheme: Loading FP Map ----------------
	iVfpMatrixRequest=0
	if (FP.iUseFP.ne.0) then	!````

	! load FP spatial map:
	open(26,file='fpmaps.dat',form='unformatted',err=142,status='old')

	! read the file header:
	read(26,err=145,end=145)	! icNone
	read(26,err=145,end=145) nt_chk,
     +	iFP_NumPoints1,iFP_NumPoints2,iFP_NumPoints3,iFP_NumPoints123,
     +	dFP_Start1,dFP_Start2,dFP_Start3,
     +	dFP_d1,dFP_d2,dFP_d3,iFP_CoordSystem,	!must be ==FP.iCoordSystem
     +	tmp,VAtShft,tmp1
	if (FP.iVShiftTo.eq.1) VAtShft=tmp1
	if (FP.iVShiftTo.eq.2) VAtShft=tmp
* Three double-precision values after CoordSystem are the average,
* maximum (highest), and minimum (lowest) potentials outside MT spheres.
* Have one of them become the shift of the energy scale.

	if (nt_chk.ne.nt) stop
     +	'FP: Error: Invalid or corrupt fpmaps.dat. Restart Strcy.'

	if (iFP_NumPoints1*iFP_NumPoints2*iFP_NumPoints3
     +		.ne.iFP_NumPoints123) stop
     +	'FP: Error: The fpmaps.dat is incomplete. Restart Fpxcy2.'

c	FP.iUseFP=1

  145	close(26)

  142	continue	! skip reading FP maps
*
*------------- Full Potential Scheme: done Loading FP Map -----------
*------------- Full Potential Scheme: Load FP Matrix Elements -------
*
	iVfpMatrixRequest=1
	allocate(VfpMEs(COUNTX*(COUNTX+1)/2),
     +	stat=istat)
	if (istat.eq.0.and.DEBUG.ge.7.and.ierr.ne.0) write(ierr,*)
     +	'FP: MatrElem: Allocated',COUNTX*(COUNTX+1)/2*ICMPLX,' bytes'
	if (istat.ne.0) stop '*ERR* Low memory [FPMEs]'
*
	open(17,file='fpmatrix.dat',form='unformatted',
     +	status='old',err=143)
	read(17) i
	if (i.ne.COUNTX) goto 144	! re-calc MEs and overwrite the file
	read(17) (VfpMEs(i),i=1,COUNTX*(COUNTX+1)/2)
c	read(17) VAtShft	! fpmaps.dat may be deleted; we also pass VAtShft through fpmatrix.dat
	iVfpMatrixRequest=0	! have read the data; no re-calc request
	write(ierr,*) 'FP:Matrix elements read from fpmatrix.dat'
  144	close(17)
  143	continue

	if (iVfpMatrixRequest.ne.0) then	! we will calculate FP matrix elements, so we need spatial map arrays

	write(ierr,*) 'FP: Request for FPMEs'

	! allocate memory:
	allocate(VfpMap(iFP_NumPoints123),cBwfMap(iFP_NumPoints123,2),
     +	xset(iFP_NumPoints123),
     +	yset(iFP_NumPoints123),
     +	zset(iFP_NumPoints123),
     +	int_map(iFP_NumPoints123),
     +	stat=istat)
	if (istat.eq.0.and.DEBUG.ge.7.and.ierr.ne.0) write(ierr,*)
     +	'FP: Mapping: Allocated',
     +	6*iFP_NumPoints123*IDBL+iFP_NumPoints123*ICHAR,' bytes'
	if (istat.ne.0) stop '*ERR* Low memory [spatialmaps]'

	istat=0
	open(26,file='fpmaps.dat',form='unformatted',err=142,status='old')
	do while(istat.lt.3)
	   read(26,err=145,end=145) i	! icXXXXX, info class
	   if (i.eq.icSpatialGrid) then
	      read(26,err=145,end=145)
     +	(xset(i),yset(i),zset(i),i=1,iFP_NumPoints123)
	      print *,'FP: read spatial grid...'
	      istat=istat+1
	   elseif (i.eq.icIntMap) then
	      read(26,err=145,end=145) (int_map(i),i=1,iFP_NumPoints123)
	      print *,'FP: read int/ext map indices...'
	      istat=istat+1
	   elseif (i.eq.icVfp) then
	      read(26,err=145,end=145) (VfpMap(i),i=1,iFP_NumPoints123)
	      print *,'FP: read spatial map of FP...'
	      istat=istat+1
	   endif
	enddo
	close(26)

	endif	! iVfpMatrixRequest.ne.0
	endif	!````FP.iUseFP.ne.0
*
*------------- Full Potential Scheme: done Load FP Matrix Elements --
*
c   -- Searching for min. value of potential on MT-spheres --
	if (FP.iUseFP.eq.0) then	! standard MT approximation
	  VAtShft=BGX(JRIS(1),1)
	  DO IT=1,NT
	    IF(VAtShft.LT.BGX(JRIS(IT),IT)) VAtShft=BGX(JRIS(IT),IT)
	  END DO
	else	! FP is on
c VAtShft is already set to be either average, or min, or max
c of Vfp outside MT spheres. (See the reading of fpmaps.dat.)
	endif

c !!! Shifting the potential !!!
	V_out=0.0d0
	DO IT=1,NT
	  DO I=1,JRIS(IT)
	    BGX(I,IT)=BGX(I,IT)-VAtShft
	  END DO
	END DO
	if (FP.iUseFP.ne.0.and.iVfpMatrixRequest.ne.0) then	! shift FP
	  do i=1,iFP_NumPoints123
	    VfpMap(i)=VfpMap(i)-VAtShft
	  enddo
	endif
c !!! Done potential shifting !!!

	WRITE(*,44) VAtShft*EVS

	CLOSE(9)

c	-- attempting to allocate the matrices for main sets of matr. elements --
	ALLOCATE(H(COUNTX,COUNTX),S(COUNTX,COUNTX),STAT=IFAIL)
	IF(IFAIL.NE.0) THEN
	  PRINT '(/A46/7X,A51/7X,A14)',
     .		'*ERR* Can''t allocate the core set of matrices.',
     .		'Basis set is too big or the total amount of memory ',
     .		'for this application is too small?'
	  STOP
	ENDIF
	print *,'Allocated H,S; ',2*COUNTX*COUNTX*ICMPLX,' bytes'
	PRINT '(1X,A28\)','Allocating of main matrices '
c    -- Attempting to allocate the matrices for partial charges --
*! SP1.1: the line below deleted:
c	IF(IQChrg.NE.0) THEN
	II=NLsav*NT+1  ! All partial charges in range.
	ALLOCATE(PoolTmp(COUNTX,COUNTX,II),STAT=IFAIL)
	IF(IFAIL.NE.0) THEN
	  DO WHILE(IFAIL.NE.0 .AND. II.GT.1)
	    II=II-1
	    ALLOCATE(PoolTmp(COUNTX,COUNTX,II),STAT=IFAIL)
	  END DO
	  IF(IFAIL.NE.0) THEN
	    PRINT '(//1X,A41/6X,A42)',
     .		'*ERR* Can''t allocate temporary matrices. ',
     .		'Try to reduce the dimension of basis set.'
	  ELSE
	    IF(II.NE.IQChrg) PRINT '(//1X,A52,I3,A37)',
     .		'*ERR* Matrices for partial charges are too big. Only',
     .		II,' layer(s) allowed for this basis set.'
	  ENDIF
	  STOP
	ENDIF
	print *,'Allocated PoolTmp; ',COUNTX*COUNTX*II*ICMPLX,' bytes'
	II=NT*NLsav+1
	ALLOCATE(PoolRes(COUNTX,NTPNTS+1,II),STAT=IFAIL)
	IF(IFAIL.NE.0) THEN
	  DO WHILE(IFAIL.NE.0 .AND. II.GT.1)
	    II=II-1
	    ALLOCATE(PoolRes(COUNTX,NTPNTS+1,II),STAT=IFAIL)
	  END DO
	  IF(IFAIL.NE.0) THEN
	    PRINT '(//1X,A22,A38)',
     .		'*ERR* Can''t allocate matrices of results. ',
     .		'Try to reduce the dimension of basis set.'
	  ELSE
	    IF(II.NE.IQChrg+1) PRINT '(//1X,A36,I3,A35)',
     .		'*ERR* Too many data for memory. Only',
     .		II,' matrices is allowed in this time.'
	  ENDIF
	  STOP
	ENDIF
	print *,'Allocated PoolRes; ',COUNTX*(NTPNTS+1)*II*IDBL,' bytes'
*! SP1.1: the line below commented:
c	ENDIF
	PRINT *, 'and big matrix block successful.'
*
	C=274.0746d0
	CIN=1.d0/(C*C)
	NORLW=0
* ------- Calc. case: P(),PR(),PE(),PER(),XE() calc. -------
	DO IT=1,NT
	  DO I=1,NGRID
	    RAD(I,IT)=DEXP(RSTART(IT)+STEP(IT)*DBLE(I-1))
	  END DO
	  RMTS(IT)=RAD(JRIS(IT),IT)
	END DO
c
c --- If NLmax > NL, EE(i,IT)=EE(NL,IT), NL < i <= NLmax ---
c
	OPEN(121, FILE='PSIALL.txt')
	  OPEN(122, FILE='PSIEALL.txt')
	  OPEN(123, FILE='PSIRALL.txt')
	  OPEN(124, FILE='PSIERALL.txt')
 
      DO IT=1, NT
	  DO i=1, JRIS(IT)
	    RA(i)=RAD(i,IT)
	    IF(bUseVpot(it)) THEN
	      BG(i)=VPot(it)	! For virtual potential of atom type IT
	    ELSE
	      BG(i)=BGX(i,IT)	! This is a normal case.
	    END IF
	  END DO
	  DO K=1, NLmax
	    IF(K.GT.NL) EE(K,IT)=EE(NL,IT)
	    CALL BASIS1(EE(K,IT),K,IT)
	  END DO
	END DO

	CLOSE(121)
      CLOSE(122)
     	CLOSE(123)
     	CLOSE(124)
*	  open (235,file='pout.txt',form='formatted')
*	    write(235,*) 'P,PE,PR,PER'
*	    do i=1,8
*	      write(235,*) P(i,1),PE(i,1),PR(i,1),PER(i,1)
*	      write(235,*) '---------------------------'
*	    enddo
*	  close(235)

c   -- Testing of the file BND.STR --
	OPEN(9,FILE='bnd.str',FORM='UNFORMATTED',STATUS='OLD',ERR=30)
	READ(9,ERR=35) TBASE,NPNT,NCPNT,I,MAXnEigv
	IF(TBASE.NE.TJNM) THEN    ! Title; if error ==> NEW file
	PRINT '(//1x,A38,A35)', '*ERR* Yor file of bands data is not',
     .		' a previously stored file. Correct it.'
	STOP
	ENDIF

	IF(NPNT.NE.NTPNTS) THEN
	  WRITE(*,29) NPNT,NTPNTS
	  STOP '----> Error in saved data.'
	ENDIF

	IF(I.NE.NLmax) THEN
	  WRITE(*,48) I
	  STOP '----> Correct this and continue.'
	ENDIF

*-------------------------------------------------------------------------
* Allocate EBND:
	print *,'Memory alloc request [EBND]: ',
     +		COUNTX*Project.iPointsTotal*IDBL,' bytes'
	allocate(EBND(MAXnEigv,NPNT+1),stat=istat)
c	allocate(EBND(COUNTX,Project.iPointsTotal),stat=istat)
	if (istat.ne.0) stop '	...failed!'
*-------------------------------------------------------------------------
	READ(9) ((EBND(I,J),J=1,NCPNT),I=1,MAXnEigv)    ! Reading of eigenvalues

c   -- Choosing for work mode: proceed with calculations or to extract data --
	IF(NCPNT-1.EQ.NPNT) THEN
c   -- here we try to extract calculated data --
	PRINT *,'All points in range are calculated.'
	WRITE(*,'(1X,A18,F8.3,A4)') 'Bands starts from ',VL,' eV,'
	PRINT '(1X,A24,F8.3,A4)', 'Upper bound of bands is ',VU,' eV.'
c     -- Reallocating of memory for operathing with results --
	DEALLOCATE(PoolTmp,PoolRes,H,S)
	ALLOCATE(PoolBuild(2+NT*NLsav,MAXnEigv*(NPNT+1)),STAT=IFAIL)
	IF(IFAIL.NE.0) STOP '*ERR* Can''t allocate matrix block.'
	ALLOCATE(IndxD(MAXnEigv*(NPNT+1)),STAT=IFAIL)
	IF(IFAIL.NE.0) STOP '*ERR* Can''t allocate matrix of indexes.'
	DO I=1, MAXnEigv    ! Convert energies from rydbergs to electronvolts.
	  DO J=1, NCPNT
	    EBND(I,J)=EBND(I,J)*EVS
	  END DO
	END DO

	IF(IQChrg.GT.0) THEN
	  PRINT '(/1X,A50\)',
     .		'Partial density of states will be found for shell:'
	  I=1
	  IT=1
	  DO WHILE(I.LE.IQChrg)
	    DO WHILE(I.LE.IQChrg .AND. IQChrgFlg(1,I).NE.IT)
	      I=I+1
	    END DO
	    IF(I.GT.IQChrg) EXIT
	    IF(IQChrgFlg(1,I).EQ.IT)
     .		PRINT '(/1X,A10,I2,A2\)','Atom type ',IT,': '
	    DO WHILE(I.LE.IQChrg .AND. IQChrgFlg(1,I).EQ.IT)
	      PRINT '(5X,A1\)',QchText(IQChrgFlg(2,I)+1)
	      I=I+1
	    END DO
	    IT=IT+1
	    IF(IT.GT.NT) EXIT
	  END DO
	ENDIF
	IF(ISmode) THEN
	  PRINT '(A1/1X,A53,A14)', ',',
     .		'also partial density of states for intersphera space ',
     .		'will be found.'
	ELSE
	  PRINT '(1X)'
	ENDIF
	DO K=2, NT*NLsav+2
	  DO I=1, MAXnEigv
	    READ(9) (PoolBuild(K,(I-1)*NCPNT+J),J=1,NCPNT)
	  END DO
	END DO
*
*-------------------------------------------------------------------------
*! changed fragment SP 5.04: args 4 and 5 are allocated dimensions
* of the EBND array (arg.1); EBND became dynamic since SP 5.04.
	CALL ExtrCalcData(EBND,MAXnEigv,NTPNTS+1,MAXnEigv,NTPNTS+1,CZ,
     +		VU,VL,11,LSMode,TBZERO,IWORK)
*! old fragment:
c	CALL ExtrCalcData(EBND,MAXnEigv,NTPNTS+1,MAXBND,MAXPNT,CZ,VU,VL,
c     +		11,LSMode,TBZERO,IWORK)
*! end fragment
*-------------------------------------------------------------------------
*
* DEVELOPER'S REMARK [Mith, 28 Mar 2000]
* ------------------
*
* In LACW SP 3.11 and younger, two intuitively ambiguous parameters
* are used to set up the number of k points in which E(k) is calculated.
* Here's the line from the Inlpcw file (without initial commenting '*'):
*
*  -10.0    40.0   161    5.D-4     1.D-4     10    11    5    8
*
* "10" means that the k axis is splitted onto 10 *intervals*;
* "11" means that 11 points should be calculated this time.
* We may want to do less, and re-start the process when convenient;
* but the maximum value is the previous number PLUS ONE.
*
* In other words, 10 is more important, but it defines that
* there will be 11 solutions to the secular equation until
* the calculation is completed and E(k) & DoS graphs can be plotted.
*
* And it is 10 what is passed into ExtrCalcData() as NTPNTS!
* (and worked there wrong until SP 3.11---funny, isn't it?..)
*
* Resume: We must pass NTPNTS+1 to ExtrCalcData() (not NTPNTS!).
*
*-------------------------------------------------------------------------
* One more remark [14 Jan 2001]:
*
* Since SP 5, the project file has new format; the number of points
* can be specified using either of the two parameters in the [Project] section:
*	CalcIntervalsTotal=20	; total count of points inside Brillouine zone,
*	CalcPointsTotal=21	; ^+1
* One may also use the new variable in the Project structure;
* the relation to the legacy variable is: Project.iPointsTotal == NTPNTS+1
*
*-------------------------------------------------------------------------
* And more remark [16 Feb 2001]:
*
* To visualise the calculated part of an _incomplete_ bnd.str,
* use a stand-alone program, EbndConv.
*
*-------------------------------------------------------------------------
*
c	-- Fermi Level calculations --
	DO I=1,MAXnEigv
	  DO J=1,NPNT+1
	    PoolBuild(1,(I-1)*(NPNT+1)+J)=EBND(I,J)
	  END DO
	END DO
	CALL indexx2(2+NT*NLsav,MAXnEigv*(NPNT+1),MAXnEigv*(NPNT+1),1,
     .			PoolBuild,IndxD)
c	-- rearranging of the matrix block PoolBuild by indexes IndxD --
	DO I=1, MAXnEigv*(NPNT+1)
	  CALL swap(I,IndxD,PoolBuild,1+NT*NLsav,MAXnEigv*(NPNT+1))
	END DO
	WRITE(*,81) PoolBuild(1,1), PoolBuild(1,(NPNT+1)*MAXnEigv)
  81	FORMAT(/1X,'Bottom of the bands ',F9.4,' eV,',/
     .		1X,'Top of the bands ',F9.4,' eV.')
c	-- calc. of full DOS and recalculating of Ql(alpha): Ql(alpha)/Qtot --
	DO I=1, MAXnEigv*(NPNT+1)
	  TEMP=0.d0
	  DO J=2, NT*NLsav+2
	    TEMP=TEMP+PoolBuild(J,I)    ! Qtot for current value of I.
	  END DO
	  DO J=2, NT*NLsav+2
	    PoolBuild(J,I)=PoolBuild(J,I)/TEMP
	  END DO
	END DO
! Count valent electrons:
	NelTot=0
	DO IT=1,NT
	IF(NelCell.GT.0) THEN
	  NelTot=NelTot+NTS(IT)*NelCell
	ELSE
	  NelTot=NelTot+NTS(IT)*NEL(IT)
	ENDIF
	ENDDO
! No more than 2e on each level; find the number of [half-]occupied
! levels, determine the offset and select EFermi from EBND:
	NEF=NelTot*(NPNT+1)/2
	TEMP=PoolBuild(1,NEF)

	OPEN(7,FILE='ef.dat',FORM='FORMATTED')
	 TMP = -PI/CZ
	 IF(TBZERO) TMP=0
	 WRITE(7,17) TMP,TEMP
	 WRITE(7,17) PI/CZ,TEMP
	CLOSE(7)
	WRITE(*,49) TEMP

c	-- Preparing of file names for output of DOS's --
	DO I=1,IQChrg+1
	  IF(I.GT.1) THEN
	    IF(IQChrgFlg(1,I-1).LT.10) THEN
	      OutfL(I)=TDOS//TNMBR(IQChrgFlg(1,I-1)+1)//
     .		QchText(IQChrgFlg(2,I-1)+1)//TDAT
	    ELSE
	      OutfL(I)=TDOS//TNMBR(IQChrgFlg(1,I-1)/10+1)//
     .		TNMBR(IQChrgFlg(1,I-1)-10*(IQChrgFlg(1,I-1)/10)+1)
     .                //QchText(IQChrgFlg(2,I-1)+1)//TDAT
	    ENDIF
	  ELSE
	    OutfL(i)=TDOS//TDAT  ! for full DOS.
	  ENDIF
	END DO
	IF(ISmode.GT.0) THEN
	  IQChrg=IQChrg+1
	  OutfL(i)=TDOS//TIS//TDAT
	  IQChrgFlg(1,IQChrg)=-1	! negative value here is flag of Qis.
	ENDIF
c	-- Building of DOS's by GAUSS broaderning --
	CALL GAUSS(10,PoolBuild,IQChrg,IQChrgFlg,MAXnEigv*(NPNT+1),
     +		VL,VU,HfWdth,HfWdt_EF,TEMP,IDOSPTS,
     +		2.D0/DBLE(NPNT+1),RWORK,D,NLsav,2+NT*NLsav,OutfL)
	OPEN(10,FILE='Dos.log',FORM='FORMATTED',ERR=89)
	DO I=1,IQChrg+1
	  IF (I.EQ.1) THEN
	    WRITE(10,87) ' Full ',D(i),RWORK(i)
c	     WRITE(10,*) '*** Full DOS: ***'
	  ELSEIF (IQChrgFlg(1,I-1).GT.0) THEN
	    WRITE(10,87) '  '//TNMBR(IQChrgFlg(1,I-1)+1)//
     .		QchText(IQChrgFlg(2,I-1)+1)//'  ',D(i),RWORK(i)
c	    WRITE(10,88) IQChrgFlg(1,I-1),IQChrgFlg(2,I-1),
c     .		QchText(IQChrgFlg(2,I-1)+1)
	  ELSE
	    WRITE(10,87) '  IS  ',D(i),RWORK(i)
c	     WRITE(10,'(/1X,A18)') 'Intersphera state:'
	  ENDIF
c	   WRITE(10,87) HfWdth, D(i)
c	   WRITE(10,87) HfWdt_EF, RWORK(i)
	END DO
	CLOSE(10)
	STOP '  <<< The End >>>'
   87	FORMAT(1X,A6,2(2X,F5.3))

   89	PRINT '(//1x,A38)','*ERR* Can''t open log-file for dos''s.'
c   87	FORMAT(5X,'For Halfwidth ',F6.3,'eV, DOS on Fermi Level is ',
c     .		F9.3,'.')
c   88	FORMAT(/1X,'Atom type ',I2,',',' state ',I1,' (',A1,'):')
c   89	PRINT '(//1x,A34)','*ERR* Can''t open dos''s log-file.'
	STOP

	ELSE
	  DO K=1, NT*NLsav+1
	    DO I=1,MAXnEigv
	      READ(9) (PoolRes(I,J,K),J=1,NCPNT)
	    END DO
	  END DO
	  ENDIF
	CLOSE(9)
*
	IF(NCPNT+NPNTS.GT.NTPNTS) THEN
	  WRITE(*,27) NPNTS,NTPNTS-NCPNT
	  NPNTS=NTPNTS-NCPNT+1
	ENDIF
	WRITE(*,26) NCPNT,NPNTS
	GO TO 31

   35	CLOSE(9)
	PRINT *,
     .	'Data not found in BND.STR ==> file replaced. NEW calculation.'
	GOTO 31

   30	PRINT *,'Previously saved data not found ==> NEW calculation.'
	NCPNT=0     ! counter of calculated points inside Brillouine zone

   31	IF(NPNT.GT.1) THEN
	  OPEN(16,FILE='lacw.log',FORM='FORMATTED',STATUS='OLD',
     .		ACCESS='APPEND',ERR=45)   ! the file exists; append to its end
	  WRITE(16,'(//10X,A33)') '*** Continuing of calculation ***'
	  GO TO 46
	ENDIF

   45	OPEN(16,FILE='lacw.log',FORM='FORMATTED') ! ...otherwise, assume a new file.

*-------------------------------------------------------------------------
* Allocate EBND (2):
	if (.not.allocated(EBND)) then
	print *,'Memory alloc request [EBND]: ',
     +		COUNTX*Project.iPointsTotal*IDBL,' bytes'
c	allocate(EBND(MAXnEigv,NCPNT),stat=istat)
	allocate(EBND(COUNTX,Project.iPointsTotal),stat=istat)
	if (istat.ne.0) stop '	...failed!'
	endif
*-------------------------------------------------------------------------
	write(16,*)'*** OVER55 *** LACW Service Pack 6.0 (based on 5.07)'
	WRITE(16,*) TJNM

   46	PRINT '(/1X,A)', 'Check your input data:'
	DO IT=1,NT
	  PRINT '(10(2x,F7.2))', (EE(I,IT)*EVS,I=1,NL)
	END DO
	PRINT '(/)'
	WRITE(16,44) VAtShft*EVS
	IF(SaveI3M) THEN
	  OPEN(15,FILE='_intgrls.da0',ACCESS='DIRECT',ACTION='READWRITE',
     .		FORM='UNFORMATTED',RECL=2*INTL+3*(2*NLmax*IDBL),
     .		STATUS='OLD',ERR=40)
	  READ(15,REC=1) I,TBASE
	  IF(TBASE.NE.TJNM) THEN    ! Title; if error ==> NEW file
	    PRINT *,'*ERR* Can''t recognize this swap file ',
     .		'for current input data. Delete it and re-start Over55.'
	    STOP '----> Reading of double integrals.'
	  ENDIF
	  IF(I.NE.Complete) THEN
	    PRINT '(1x,A54/6X,A35/)',
     .	'*WRN* File of swapping of double integrals is incomplete.',
     .	'Skipping. File will be formed now.'
	    CLOSE(15,STATUS='DELETE')
	    GO TO 40
	  ENDIF
	  PRINT '(1X,A48/)',
     .		'File with stored double integrals has been read.'
	  Create_I3M=.FALSE.      ! use saved data, but not a new.
	  GO TO 41
   40	  OPEN(15,FILE='_intgrls.da0',ACCESS='DIRECT',ACTION='READWRITE',
     .		FORM='UNFORMATTED',RECL=2*INTL+3*(2*NLmax*IDBL),
     .		STATUS='NEW')
	  IF(NCPNT.GT.0) PRINT *,
     .		'WARNING: old file with double integrals formed now!'
*   -- Writing of 0 as marker of incomplete file and job name to file header --
	  WRITE(15,REC=1) 0,TJNM
	ENDIF

   41	NCXPOINT=COUNTX*(COUNTX-1)/2+COUNTX	! number of elements to be calculated in triangular matrices (S,H)

* --- finding of atoms with same RO ---
	DO IT=1, NT
	  DO I=1, NTS(IT)
	    IRO(IT,I)=0
	  END DO
	END DO
	DO IT=1, NT
	  K=0
	  DO I=1, NTS(IT)
	    IF(IRO(IT,I).EQ.0) THEN
	      K=K+1
	      IWORK(I)=0
	      DO J=1, NTS(IT)
		IF(I.EQ.1)
     .		RWORK(J)=DSQRT(RSX(J,IT)*RSX(J,IT)+RSY(J,IT)*RSY(J,IT))
		IF(DABS(RWORK(J)-RWORK(I)).LT.EpsG) THEN
		  IRO(IT,J)=I
		  IWORK(I)=IWORK(I)+1
		ENDIF
	      END DO
	    ENDIF
	  END DO
	IF(NPNT.EQ.0) THEN	! show geometry parameters only once
	  WRITE(16,'(/1X,A17)') 'System geometry:'
	  IF(K.GT.1) THEN
	    WRITE(*,53) IT
	    WRITE(*,56) K,NTS(IT)
	  ENDIF
	  WRITE(16,'(1X,A14,I2,A1)') 'For atom type ',IT,':'
	  WRITE(16,'(1X,A38)') '  No    Distance   Unique ID Frequency'
	  WRITE(16,58)
	  DO I=1, NTS(IT)
	    IF(IWORK(I).NE.0) THEN
	      WRITE(16,54) I,RWORK(I),IRO(IT,I),IWORK(I)
	    ELSE
	      WRITE(16,57) I,RWORK(I),IRO(IT,I)
	    ENDIF
	  END DO
	  WRITE(16,58)
	  WRITE(16,56) K,NTS(IT)
	ENDIF
	END DO
*
*------------- Full Potential Scheme: BWF Mapping -------------------
*
	if (FP.iUseFP.ne.0.and.iVfpMatrixRequest.ne.0) then	! FP on

        write(*,*) 'Full potential is on!'
        write(16,*) 'Full potential is on!'

	write(*,*) 'Looking for BWF maps (bwfmaps.dat)...'
* Uncomment these lines to use previously calculated BWF maps:
	open(14,file='bwfmaps.dat',access='direct',action='readwrite',
     +		form='unformatted',recl=(iFP_NumPoints123*ICMPLX),
     +		status='old',err=140)
	goto 141
  140	write(*,*) 'Building BWF spatial maps...'
	write(16,*) 'Building BWF spatial maps...'
	open(14,file='bwfmaps.dat',access='direct',action='readwrite',
     +		form='unformatted',recl=(iFP_NumPoints123*ICMPLX))

	do i=1,COUNTX

	M = Basis.INDXX(1,I)
	N = Basis.INDXX(2,I)
	NP= Basis.INDXX(3,I)

	VKP(1)=(2.D0*PI*DBLE(NP)-DBLE(M)*ChiralityAngle)/CZ

	if (iUseCoredTube.eq.0) then
	  CALL BESS_NUL(N,IABS(M),RES,MAXSUM,1)
	  KPNM(1)=RES(N)/A
	else if (iUseCoredTube.eq.1.or.iUseCoredTube.eq.3.or.
	*iUseCoredTube.eq.4) then
	  KPNM(1)=AllKappas(iabs(M)+1,N)
	else if (iUseCoredTube.eq.2) then
	  CALL BESS_NUL(N,IABS(M),RES,MAXSUM,2)
	  KPNM(1)=RES(N)/B
	endif

	do ii=1,iFP_NumPoints123
	Rho=xset(ii)	! !!!!!!!!WARNING: Ambiguous coord conversion
	Phi=yset(ii)
	Z=zset(ii)

	if (iUseCoredTube.eq.0) then
	  call Bess(M,KPNM(1)*Rho,ret_MOD,temp)
	elseif (iUseCoredTube.eq.2) then
	  call BessYproc(M,KPNM(1)*Rho,zf_params.eps,ret_MOD,temp)
	elseif (iUseCoredTube.eq.1.or.iUseCoredTube.eq.4) then
	  call Bess(M,KPNM(1)*Rho,dJ,temp)
	  call BessYproc(M,KPNM(1)*Rho,zf_params.eps,dY,temp)
	  tmp=abJ(iabs(M)+1,N,1)/abY(iabs(M)+1,N,1)
	  ret_MOD=dJ-tmp*dY
	elseif (iUseCoredTube.eq.3) then
	  call Bess(M,KPNM(1)*Rho,dJ,temp)
	  call BessYproc(M,KPNM(1)*Rho,zf_params.eps,dY,temp)
	  tmp=abJ(iabs(M)+1,N,2)/abY(iabs(M)+1,N,2)
	  ret_MOD=dJ-tmp*dY
	endif
	ret_ARG=VKP(1)*Z+dble(M)*Phi

	cBwfMap(ii,1)=dcmplx(ret_MOD*dcos(ret_ARG),ret_MOD*dsin(ret_ARG))
c	cBwfMap(ii,1)=dcmplx(ret_MOD,ret_ARG)	!!!!! not a real Psi(MNP)

	enddo	!ii (walking over (x..z)set)

	write(14,rec=i) (cBwfMap(j,1),j=1,iFP_NumPoints123)
	write(*,'(a1,i3,''%''\)') char(13),100*i/countx

	enddo	! cycle over BWFs
        write(*,*) ' done.'

  141	continue	! skipped BWF mapping

	endif	! FP is on and iVfpMatrixRequest.ne.0

d	stop '*TEST* FP: BWF mapping done.'
*
*------------- Full Potential Scheme: BWF Mapping (done) ------------
*
*
*
d	write(16,*) 'FP.iUseFP,iVfpMatrixRequest,VAtShft,V_out:',
d     +	FP.iUseFP,iVfpMatrixRequest,VAtShft,V_out
d	stop '*TEST* On the main cycle entry'
*
	ITRX=1
	WRITE(*,18) ITRX
	WRITE(16,18) ITRX
	WRITE(*,*)

	IF(TBZERO) THEN
	  k_start = 0.0d0
	ELSE
	  k_start = -PI/CZ
	END IF
	open(121, file='kstart_kstop.txt')
	write(121, *) k_start, PI/CZ
	close(121)
	OPEN(121, FILE='NATOMS.txt')
	DO IT = 1, NT
		WRITE(121, *) NTS(IT)
	END DO 
	CLOSE(121)
	OPEN(121, FILE='COORD.txt')
	DO IT = 1, NT
	  DO IQ = 1, NTS(IT)
	    RO=DSQRT(RSX(IQ,IT)*RSX(IQ,IT)+RSY(IQ,IT)*RSY(IQ,IT))

           PHI=0.D0
           IF(RO.GT.Eps) PHI=DASIN(DABS(RSY(IQ,IT))/RO)
           IF(RSY(IQ,IT).LT.0.D0) THEN
             IF(RSX(IQ,IT).LT.0.D0) THEN
               PHI=PI+PHI
             ELSE
               PHI=-PHI
             ENDIF
           ELSE
             IF(RSX(IQ,IT).LT.0.D0) PHI=PI-PHI
           ENDIF
		 WRITE(121, *) Ro, Phi, RSZ(IQ, IT)
 
	  END DO
	END DO
	close(121)

* ----- Main calc. begins here -----
	DO 77 NPNT=NCPNT, NCPNT+NPNTS-1
	NL=NLsav
	IF(TBZERO) THEN
	  VK = DBLE(NPNT)*PI/DBLE(NTPNTS)/CZ
*	  VK = DBLE(NPNT)*PI/DBLE(NTPNTS-1)/CZ
	ELSE
	  VK = -PI/CZ+DBLE(2*NPNT)*PI/DBLE(NTPNTS)/CZ
*	  VK = -PI/CZ+DBLE(2*NPNT)*PI/DBLE(NTPNTS-1)/CZ
	ENDIF
	NProgPNT=0
	NextPnt=NCXPOINT/NPROGRES	! First step on progres bar.
	WRITE(*,19) NPNT+1,NPNT-NCPNT+1,NPNTS,NTPNTS
	TIME=0.D0
	CALL TTIME(TIME,TATIME,0)
	CALL PROGRES(0) ! draw progres bar for first time...
	IF(NPNT.eq.0) THEN
	  outfile=TSEE//TNMBR(NPNT+1)//TDAT
	  OPEN(11,FILE=outfile, FORM='FORMATTED')	! temporary log-file.
	ENDIF
	MOV=2	! start record number for read/write of double integrals

*! begin debug
d	write(24,*) '***** J/Y coeffs *****'
d	write(24,*) '    i     m   n   p',
d     .	'        Ja        Ya       Ja''       Ya''',
d     .	'        Jb        Yb       Jb''       Yb'''
d	do ish=1,countx
d	  m=Basis.INDXX(1,ish)
d	  n=Basis.INDXX(2,ish)
d	  np=Basis.INDXX(3,ish)
d	  write(24,'(i4,2x,3i4,9f10.4)') ish,m,n,np,
d     .	AllKappas(iabs(m)+1,n),
d     .	abJ(iabs(m)+1,n,1),abY(iabs(m)+1,n,1),
d     .	abJD(iabs(m)+1,n,1),abYD(iabs(m)+1,n,1),
d     .	abJ(iabs(m)+1,n,2),abY(iabs(m)+1,n,2),
d     .	abJD(iabs(m)+1,n,2),abYD(iabs(m)+1,n,2)
d	enddo
d	stop '*TEST* AllKappas dumped to stream #24'
*! end debug

c ----- Calc. on points begins here -----

	iVfpMEs=0	! offset in the VfpMEs array

	DO 7 ISH=1, COUNTX
	  M1= Basis.INDXX(1,ISH)
	  N1= Basis.INDXX(2,ISH)
	  NP1=Basis.INDXX(3,ISH)

	  IF(NPNT.EQ.0) WRITE(11,*) 'ROW ',ISH  ! for test
c -- Due to Hermitian condition H(ISH,JSH)=H(JSH,ISH)* only --
c  -- the upper triangle of matrix H is to be calculated. --

	  if (FP.iUseFP.ne.0.and.iVfpMatrixRequest.ne.0) then	! cache BWF[ish] map for FPMEs calc
	    read(14,rec=ish) (cBwfMap(j,2),j=1,iFP_NumPoints123)
c and, later,	    read(14,rec=jsh) (cBwfMap(j,1),j=1,iFP_NumPoints123)
	  endif

	DO 7 JSH=ISH, COUNTX
	  M = Basis.INDXX(1,JSH)
	  N = Basis.INDXX(2,JSH)
	  NP= Basis.INDXX(3,JSH)
*! Developer's comment [Mith, 20.01.2001]:
* ISH=1..COUNTX and JSH=ISH..COUNTX,
* --> S(ISH,JSH), H(ISH,JSH),
* and the exponent coeff is exp( (pi*(NP-NP1)/Cz) *Rho + (M-M1)*phi);
* thus, the multiplication is dconjg(Psi(ISH))*Psi(JSH),

	DO I=1, NT*NLsav+1
	  PoolTmp(ISH,JSH,I)=DCMPLX(0.D0,0.D0)
	END DO
	IF(M.GT.MAXSUM.OR.M1.GT.MAXSUM) THEN
	  WRITE(*,20) 'Number of solutions is greater than ',MAXSUM
	  STOP '--- OVERLAP{1} --->'
	ENDIF

*! begin fragment SP 3.1
c Bess_nul(N,M,dest,dest_size,bess_family)
c searches for N roots of a Bessel function J (bess_family==1)
c or Y (bess_family==2) of order M and stores those roots into
c the dest array.
	if (iUseCoredTube.eq.0) then
	  CALL BESS_NUL(N,IABS(M),RES,MAXSUM,1)
	  KPNM(1)=RES(N)/A
	  CALL BESS_NUL(N1,IABS(M1),RES,MAXSUM,1)
	  KPNM(2)=RES(N1)/A
	else if (iUseCoredTube.eq.1.or.iUseCoredTube.eq.3.or.
	*iUseCoredTube.eq.4) then
	  KPNM(1)=AllKappas(iabs(M)+1,N)
	  KPNM(2)=AllKappas(iabs(M1)+1,N1)
	else if (iUseCoredTube.eq.2) then
	  CALL BESS_NUL(N,IABS(M),RES,MAXSUM,2)
	  KPNM(1)=RES(N)/B
	  CALL BESS_NUL(N1,IABS(M1),RES,MAXSUM,2)
	  KPNM(2)=RES(N1)/B
	endif
c*! end fragment

* FP: Read BWF maps and prepare FP matrix element, <i|Vfp|j>:
	if (FP.iUseFP.ne.0) then	! FP on

	iVfpMEs=iVfpMEs+1

	if (iVfpMatrixRequest.ne.0) then	! FP ME re-calc request

c	  read(14,rec=ish) (cBwfMap(j,2),j=1,iFP_NumPoints123)
c	  read(14,rec=jsh) (cBwfMap(j,1),j=1,iFP_NumPoints123)
	  if (jsh.eq.ish) then
	    do j=1,iFP_NumPoints123
	      cBwfMap(j,1)=cBwfMap(j,2)
	    enddo
	  else
	    read(14,rec=jsh) (cBwfMap(j,1),j=1,iFP_NumPoints123)
	  endif

	! integration using the spatial grid over elementary cell:
	VfpMEInt=dcmplx(0.0d0,0.0d0)
	do ii=1,iFP_NumPoints123
	if (int_map(ii).ne.char(0)) cycle	! skip MT spheres
	x=xset(ii)
	y=yset(ii)
	z=zset(ii)
	domega=dVolumeCoeff(FP.iCoordSystem,x,y,z)

c	! cBwfMap(ii,m)==(mod,arg), not (real,imag)!!!
c	tm1=dreal(cBwfMap(ii,1))*dreal(cBwfMap(ii,2))
c	tm2=dimag(cBwfMap(ii,1))-dimag(cBwfMap(ii,2))

	VfpMEInt=VfpMEInt+domega*VfpMap(ii)*
c     +	tm1*dcmplx(dcos(tm2),dsin(tm2))
     +	dconjg(cBwfMap(ii,2))*cBwfMap(ii,1)
c     +	dconjg(cBwfMap(ii,1))*cBwfMap(ii,2)	!!!--WRONG!
	! ^and *(dx*dy*dz)*CjMN*CjM1N1/(2*Pi*Cz) later!
d	print *,ii,VfpMap(ii),cBwfMap(ii,2),cBwfMap(ii,1)

	enddo	!ii

	VfpMEs(iVfpMEs)=VfpMEInt*dFP_d1*dFP_d2*dFP_d3/(2.0d0*Pi)	! multiply by domega and 1/(2pi)


d	if (iUseCoredTube.eq.0) then
d	TMP=S_NJM(KPNM(1),KPNM(2),M,M1,A)/CZ
*! new fragment SP 3.1
d	else if (iUseCoredTube.eq.2) then
d	   YEps=zf_params.eps
d	   TMP=S_NYM(KPNM(1),KPNM(2),M,M1,B,YEps)/CZ
*! end fragment SP 3.1
d	else if (iUseCoredTube.eq.1) then
d	   tmp=S_CJY(A,B,M,N,M1,N1,MMAX,NMAX,abJ,abJD,abY,abYD)/CZ
d	endif
d	write(16,'(8i3,2f10.4)') ish,jsh,M,N,NP,M1,N1,NP1,
d     +	dreal(VfpMEs(iVfpMEs))*TMP,dimag(VfpMEs(iVfpMEs))*TMP
d	write(*,'(8i3,2f10.4)') ish,jsh,M,N,NP,M1,N1,NP1,
d     +	dreal(VfpMEs(iVfpMEs))*TMP,dimag(VfpMEs(iVfpMEs))*TMP


	endif	! FP ME re-calc request

	endif	! FP on

d	cycle	! *TEST* FP matrix elements


*! old fragment SP 1.x
c	VKP(1)=2.D0*PI*DBLE(NP)/CZ+VK
c	VKP(2)=2.D0*PI*DBLE(NP1)/CZ+VK
*! new fragment SP 2.02

c	IF(TBZERO) THEN
c	  VK = DBLE(NPNT)*PI/DBLE(NTPNTS)/CZ
c	  VK = (DBLE(NPNT)*PI/DBLE(NTPNTS) - DBLE(M)*ChiralityAngle)/CZ
*	  VK = DBLE(NPNT)*PI/DBLE(NTPNTS-1)/CZ
c	ELSE
c	  VK = -PI/CZ+DBLE(2*NPNT)*PI/DBLE(NTPNTS)/CZ
*	  VK = -PI/CZ+DBLE(2*NPNT)*PI/DBLE(NTPNTS-1)/CZ
c	ENDIF

	VKP(1)=(2.D0*PI*DBLE(NP)-DBLE(M)*ChiralityAngle)/CZ + VK
	VKP(2)=(2.D0*PI*DBLE(NP1)-DBLE(M1)*ChiralityAngle)/CZ + VK
*! end fragment
	S_RE = 0.D0
	S_IMG= 0.D0
	SH_RE= 0.D0
	SH_IMG=0.D0

	IF(NLmax.GT.LCUT) THEN
	  WRITE(*,20) 'NLmax from INSUPER exceed max. value', LCUT
	  STOP '--- 8th parameter in your INLPWCY incorrect --->'
	ENDIF
	IF(NL.GT.LCUT) THEN
	  WRITE(*,20) 'NL from INSUPER exceed max. value', LCUT
	  STOP '----> 7th parameter in your INLPWCY is incorrect.'
	ENDIF
         DO 1 IT=1, NT
c         .. for any type of atoms ..
          DO J=1, 2*NLmax
            LCALC(J)=.FALSE.
            BDCALC(J)=0.D0
            BECALC(J)=0.D0
            BFCALC(J)=0.D0
          END DO
          LI12CALC=0
          DO J=1,NLmax
            DO K=1,NLmax
              SI1(J,K) =DCMPLX(0.D0,0.D0)
              SI1_(J,K)=DCMPLX(0.D0,0.D0)
              SI2(J,K) =DCMPLX(0.D0,0.D0)
              SI2_(J,K)=DCMPLX(0.D0,0.D0)
            END DO
          END DO
          IF(.NOT.Create_I3M .AND. SaveI3M) THEN
c         .. here we are read previously saved values of double integrals ..
            READ(15,REC=MOV,ERR=82) II,K,(BDCALC(I),BECALC(I),BFCALC(I),
     .                       I=K,K+II-1)
            DO I=K, K+II-1
              LCALC(I)=.TRUE.    ! ...and set flags.
            END DO
          ENDIF
          GO TO 84
   82     WRITE(*,83) MOV,II,K+II,2*NLmax,'read.'
          STOP '---> Internal error.'
c         .. for all atoms of type IT ..
   84     DO 2 IQ=1, NTS(IT)
           LSmCALC(IQ)=.FALSE.
           DO I=1, 3+NT*NLsav
             BSmCALC(I,IQ)=0.d0
           END DO
           DO I=1, 2*NLmax
             BJCALC(1,I)=0.D0
             BJCALC(2,I)=0.D0
*! begin fragment
             BYCALC(1,I)=0.D0
             BYCALC(2,I)=0.D0
*! end fragment
             LBCALC(I)=.FALSE. ! For calc. of BESSEL J functions.
           END DO
           RO=DSQRT(RSX(IQ,IT)*RSX(IQ,IT)+RSY(IQ,IT)*RSY(IQ,IT))
*
           PHI=0.D0
           IF(RO.GT.Eps) PHI=DASIN(DABS(RSY(IQ,IT))/RO)
           IF(RSY(IQ,IT).LT.0.D0) THEN
             IF(RSX(IQ,IT).LT.0.D0) THEN
               PHI=PI+PHI
             ELSE
               PHI=-PHI
             ENDIF
           ELSE
             IF(RSX(IQ,IT).LT.0.D0) PHI=PI-PHI
           ENDIF
*
           IF(.NOT.LSmCALC(IRO(IT,IQ))) THEN
            IM=0
            IIM=0
            CALC=1 ! forward calc.
            DO WHILE(CALC.NE.0)
             IF(.NOT.LBCALC(IM+NLmax)) THEN
               CALL BESS(IM-M,KPNM(1)*RO,BJCALC(1,IM+NLmax),TMP)
               CALL BESS(IM-M1,KPNM(2)*RO,BJCALC(2,IM+NLmax),TMP)
*! begin fragment
      if (iUseCoredTube.ne.0) then	! ".ne.0" is OK (for 1 and 2)
           YEps=zf_params.eps
d           print *,'Ymx:',im,m,m1,RO,kpnm(1),kpnm(2)
           CALL BESSYPROC(IM-M,KPNM(1)*RO,YEps,BYCALC(1,IM+NLmax),TMP)
           CALL BESSYPROC(IM-M1,KPNM(2)*RO,YEps,BYCALC(2,IM+NLmax),TMP)
d           write(*,*) 'im,m,m1:',im,m,m1
d           write(24,*)
d     . IM-M,BJCALC(1,IM+NLmax),BYCALC(1,IM+NLmax),
d     . IM-M1,BJCALC(2,IM+NLmax),BYCALC(2,IM+NLmax)
      endif
*! end fragment
               LBCALC(IM+NLmax)=.TRUE.
             ENDIF
             IF(RO.GT.Eps .OR. DABS(BJCALC(1,IM+NLmax)).GT.Eps
     .                      .AND. DABS(BJCALC(2,IM+NLmax)).GT.Eps) THEN
c -- we don't calculate this part, when any of Bessel J-functions are 0 --
             IF(.NOT.LCALC(IM+NLmax)) THEN
*          -- new calc. of integrals, if need (next line only) --
               IF(.NOT.Create_I3M.AND.SaveI3M) Add_to_Saved_I3M=.TRUE.
               CALL D01DAFM(0.D0,RMTS(IT),0.D0,PI/2.D0,FD,Eps,
     .                BDCALC(IM+NLmax),NPTS,IM,2.D0*PI*DBLE(NP-NP1)/CZ,
     .                KPNM(2),KPNM(1))
               BDCALC(IM+NLmax)=2.D0*BDCALC(IM+NLmax)
               CALL D01DAFM(0.D0,RMTS(IT),0.D0,PI/2.D0,FE,Eps,
     .                BECALC(IM+NLmax),NPTS,IM,2.D0*PI*DBLE(NP-NP1)/CZ,
     .                KPNM(2),KPNM(1))
               BECALC(IM+NLmax)=2.D0*BECALC(IM+NLmax)
               CALL D01DAFM(0.D0,RMTS(IT),0.D0,PI/2.D0,FF,Eps,
     .                BFCALC(IM+NLmax),NPTS,IM,2.D0*PI*DBLE(NP-NP1)/CZ,
     .                KPNM(2),KPNM(1))
               BFCALC(IM+NLmax)=2.D0*BFCALC(IM+NLmax)
               LCALC(IM+NLmax)=.TRUE.
             ENDIF
             SSI12(IIM+1)=0.D0  ! here we calc. I1 and I2 integrals
             HSI12(IIM+1)=0.D0
             BeforeExit=.FALSE.
             DO 4 L=IIM, NL
c        -- this part calculates for atom type IT --
              IF(LI12CALC.LT.IIM+1) THEN
                SI1(IIM+1,L+1) = SIMP_I12(Eps,RMTS(IT),VKP(1),KPNM(1),
     .                                    IIM,L,1)
                SI1_(IIM+1,L+1)= SIMP_I12(Eps,RMTS(IT),VKP(2),KPNM(2),
     .                                    IIM,L,1)
                SI2(IIM+1,L+1) = SIMP_I12(Eps,RMTS(IT),VKP(1),KPNM(1),
     .                                    IIM,L,2)
                SI2_(IIM+1,L+1)= SIMP_I12(Eps,RMTS(IT),VKP(2),KPNM(2),
     .                                    IIM,L,2)
              ENDIF
              IF(MOD(IIM+L,2).EQ.0) THEN
               TMP=(DBLE(SI2_(IIM+1,L+1))*PE(L+1,IT)-DBLE(SI1_(IIM+1,
     .          L+1))*PER(L+1,IT))*(DBLE(SI2(IIM+1,L+1))*PE(L+1,IT)-
     .          DBLE(SI1(IIM+1,L+1))*PER(L+1,IT))+XE(L+1,IT)*
     .          (DBLE(SI1_(IIM+1,L+1))*PR(L+1,IT)-DBLE(SI2_(IIM+1,L+1))*
     .          P(L+1,IT))*(DBLE(SI1(IIM+1,L+1))*PR(L+1,IT)-
     .          DBLE(SI2(IIM+1,L+1))*P(L+1,IT))
               TM1=(DBLE(SI2_(IIM+1,L+1))*DBLE(SI1(IIM+1,L+1))+
     .          DBLE(SI1_(IIM+1,L+1))*DBLE(SI2(IIM+1,L+1)))*PE(L+1,IT)*
     .          PR(L+1,IT)-DBLE(SI2_(IIM+1,L+1))*DBLE(SI2(IIM+1,L+1))*
     .          PE(L+1,IT)*P(L+1,IT)-DBLE(SI1_(IIM+1,L+1))*
     .          DBLE(SI1(IIM+1,L+1))*PER(L+1,IT)*PR(L+1,IT)
               ELSE
               TMP=(DIMAG(SI2_(IIM+1,L+1))*PE(L+1,IT)-DIMAG(SI1_(IIM+1,
     .          L+1))*PER(L+1,IT))*(DIMAG(SI2(IIM+1,L+1))*PE(L+1,IT)-
     .          DIMAG(SI1(IIM+1,L+1))*PER(L+1,IT))+XE(L+1,IT)*
     .          (DIMAG(SI1_(IIM+1,L+1))*PR(L+1,IT)-
     .          DIMAG(SI2_(IIM+1,L+1))*P(L+1,IT))*
     .          (DIMAG(SI1(IIM+1,L+1))*PR(L+1,IT)-
     .          DIMAG(SI2(IIM+1,L+1))*P(L+1,IT))
               TM1=(DIMAG(SI2_(IIM+1,L+1))*DIMAG(SI1(IIM+1,L+1))+
     .          DIMAG(SI1_(IIM+1,L+1))*DIMAG(SI2(IIM+1,L+1)))*
     .          PE(L+1,IT)*PR(L+1,IT)-DIMAG(SI2_(IIM+1,L+1))*
     .          DIMAG(SI2(IIM+1,L+1))*PE(L+1,IT)*P(L+1,IT)-
     .          DIMAG(SI1_(IIM+1,L+1))*DIMAG(SI1(IIM+1,L+1))*
     .          PER(L+1,IT)*PR(L+1,IT)
              ENDIF
              TMP1=DBLE(2*L+1)/FCT(L,IIM)
              TM1=(TMP*EE(L+1,IT)+TM1)*TMP1
              TM2=TMP*TMP1
              HSI12(IIM+1)=HSI12(IIM+1)+TM1
              SSI12(IIM+1)=SSI12(IIM+1)+TM2
*           -- calculation of elements of sums for partial charges --
              IF( IIM.LE.L .AND. L.LE.NLsav ) THEN
                 I=(IT-1)*NLsav+L+1
*! old fragment:
c                 BSmCALC(2+I,IRO(IT,IQ))=BSmCALC(2+I,IRO(IT,IQ))
c     .                      +TM2*BJCALC(1,IM+NLmax)*BJCALC(2,IM+NLmax)
*! new fragment
              if (iUseCoredTube.eq.0) then
                 BSmCALC(2+I,IRO(IT,IQ))=BSmCALC(2+I,IRO(IT,IQ))+TM2*
     . BJCALC(1,IM+NLmax)*BJCALC(2,IM+NLmax)
*! new fragment SP 3.1
              else if (iUseCoredTube.eq.2) then
                 BSmCALC(2+I,IRO(IT,IQ))=BSmCALC(2+I,IRO(IT,IQ))+TM2*
     . BYCALC(1,IM+NLmax)*BYCALC(2,IM+NLmax)
*! end fragment SP 3.1
              else if (iUseCoredTube.eq.1.or.iUseCoredTube.eq.4) then
* ChooseBessAtAB() returns abXX with the sign correction for negative m's,
* but here we don't need to do it, since (-1)/(-1)==(+1)
                 tmp=abJ(iabs(m)+1,n,1)/abY(iabs(m)+1,n,1)
                 tmp1=abJ(iabs(m1)+1,n1,1)/abY(iabs(m1)+1,n1,1)
c                 tmp=ChooseBessAtAB(m,n,MMAX,NMAX,abJ,1)
c     .              /ChooseBessAtAB(m,n,MMAX,NMAX,abY,1)
c                 tmp1=ChooseBessAtAB(m1,n1,MMAX,NMAX,abJ,1)
c     .               /ChooseBessAtAB(m1,n1,MMAX,NMAX,abY,1)
                 BSmCALC(2+I,IRO(IT,IQ))=BSmCALC(2+I,IRO(IT,IQ))+TM2*
     . (BJCALC(1,IM+NLmax)-BYCALC(1,IM+NLmax)*tmp)*
     . (BJCALC(2,IM+NLmax)-BYCALC(2,IM+NLmax)*tmp1)

d      write(24,'(2(2i4,3f10.4))')
d     . im,m,BJCALC(1,IM+NLmax),BYCALC(1,IM+NLmax),tmp,
d     . im,m1,BJCALC(2,IM+NLmax),BYCALC(2,IM+NLmax),tmp1
d      print *,'tm12~:',im,m,m1,n,n1,abY(m,n,1),abY(m,n,2)
d      if (dabs(abY(iabs(m)+1)).lt.zf_params.eps) stop 'game over'
              else if (iUseCoredTube.eq.3) then
                 tmp=abJ(iabs(m)+1,n,2)/abY(iabs(m)+1,n,2)
                 tmp1=abJ(iabs(m1)+1,n1,2)/abY(iabs(m1)+1,n1,2)
                 BSmCALC(2+I,IRO(IT,IQ))=BSmCALC(2+I,IRO(IT,IQ))+TM2*
     . (BJCALC(1,IM+NLmax)-BYCALC(1,IM+NLmax)*tmp)*
     . (BJCALC(2,IM+NLmax)-BYCALC(2,IM+NLmax)*tmp1)
                 endif
*! end fragment
              ENDIF
              IF( DABS(TM1).LT.Eps .AND. DABS(TM2).LT.Eps ) THEN
                IF(BeforeExit) THEN
                  EXIT
                ELSE
                  BeforeExit=.TRUE.
                ENDIF
              ELSE
                BeforeExit=.FALSE.
              ENDIF
    4        CONTINUE
            IF(LI12CALC.LT.IIM+1) LI12CALC=IIM+1
*! SP1.1: TMP divided by 2:
            TMP=RMTS(IT)**4.D0 /2.d0
            SSI12(IIM+1)=SSI12(IIM+1)*TMP
            HSI12(IIM+1)=HSI12(IIM+1)*TMP
*! begin old fragment
c            TMP=BJCALC(1,IM+NLmax)*BJCALC(2,IM+NLmax)
c     .          *(BDCALC(IM+NLmax)-SSI12(IIM+1))
c            TMPH=BJCALC(1,IM+NLmax)*BJCALC(2,IM+NLmax)
c     .          *((VKP(1)*VKP(2)+V_out)*BDCALC(IM+NLmax)+KPNM(1)*KPNM(2)
c*     .         *((VKP(1)*VKP(2))*BDCALC(IM+NLmax)+KPNM(1)*KPNM(2)
c     .          *BECALC(IM+NLmax)+DBLE(IM*IM)*BFCALC(IM+NLmax)
c     .          -HSI12(IIM+1))
c            TMPIS=BJCALC(1,IM+NLmax)*BJCALC(2,IM+NLmax)*BDCALC(IM+NLmax)
*! end old fragment
*! begin fragment
            if (iUseCoredTube.eq.0) then
               tm1=BJCALC(1,IM+NLmax)
               tm2=BJCALC(2,IM+NLmax)
*! new fragment SP 3.1
            else if (iUseCoredTube.eq.2) then
               tm1=BYCALC(1,IM+NLmax)
               tm2=BYCALC(2,IM+NLmax)
*! end fragment SP 3.1
            else if (iUseCoredTube.eq.1.or.iUseCoredTube.eq.4) then
               tmp=abJ(iabs(m)+1,n,1)/abY(iabs(m)+1,n,1)
               tmp1=abJ(iabs(m1)+1,n1,1)/abY(iabs(m1)+1,n1,1)
c               tmp=ChooseBessAtAB(m,n,MMAX,NMAX,abJ,1)
c     .            /ChooseBessAtAB(m,n,MMAX,NMAX,abY,1)
c               tmp1=ChooseBessAtAB(m1,n1,MMAX,NMAX,abJ,1)
c     .             /ChooseBessAtAB(m1,n1,MMAX,NMAX,abY,1)
               tm1=(BJCALC(1,IM+NLmax)-BYCALC(1,IM+NLmax)*tmp)
               tm2=(BJCALC(2,IM+NLmax)-BYCALC(2,IM+NLmax)*tmp1)
d      write(24,'(2(2i4,3f10.4))')
d     . im,m,BJCALC(1,IM+NLmax),BYCALC(1,IM+NLmax),tmp,
d     . im,m1,BJCALC(2,IM+NLmax),BYCALC(2,IM+NLmax),tmp1
            else if (iUseCoredTube.eq.3) then
               tmp=abJ(iabs(m)+1,n,2)/abY(iabs(m)+1,n,2)
               tmp1=abJ(iabs(m1)+1,n1,2)/abY(iabs(m1)+1,n1,2)
               tm1=(BJCALC(1,IM+NLmax)-BYCALC(1,IM+NLmax)*tmp)
               tm2=(BJCALC(2,IM+NLmax)-BYCALC(2,IM+NLmax)*tmp1)
            endif
            TMP=tm1*tm2
     .          *(BDCALC(IM+NLmax)-SSI12(IIM+1))
            TMPH=tm1*tm2
     .          *((VKP(1)*VKP(2)+V_out)*BDCALC(IM+NLmax)+KPNM(1)*KPNM(2)
*     .         *((VKP(1)*VKP(2))*BDCALC(IM+NLmax)+KPNM(1)*KPNM(2)
     .          *BECALC(IM+NLmax)+DBLE(IM*IM)*BFCALC(IM+NLmax)
     .          -HSI12(IIM+1))
            TMPIS=tm1*tm2*BDCALC(IM+NLmax)
d      write(23,*) 'tmbd',m,n,np,m1,n1,np1,im,
d     . BDCALC(IM+NLmax),SSI12(IIM+1),HSI12(IIM+1),
d     . (VKP(1)*VKP(2)+V_out)*BDCALC(IM+NLmax),
d     . KPNM(1)*KPNM(2)*BECALC(IM+NLmax),
d     . DBLE(IM*IM)*BFCALC(IM+NLmax)
*! end fragment
	   ELSE
	    TMP=0.D0	! For case, when any of Bessel J-functions are 0
	    TMPH=0.D0
	    TMPIS=0.D0
	   ENDIF

	    IF(RO.GT.Eps) THEN
	      IF(IM.NE.0.AND.IM.NE.-1.AND.
     .			DABS(TMP).LT.Eps.AND.DABS(TMPH).LT.Eps) THEN
		IF(CALC.EQ.1) THEN
		  CALC=2	! backward calc.
		  IM=-1
		ELSE
		  CALC=0	! calc. finished.
		ENDIF
	      ELSE
		IF(CALC.EQ.1) THEN
		  IM=IM+1
		ELSE
		  IM=IM-1
		ENDIF
	      ENDIF
	    ELSE
	      IF(IIM.EQ.MAX0(IABS(M),IABS(M1)).OR.IIM.EQ.NLmax-1) THEN
		IF(CALC.EQ.1) THEN
		  CALC=2
		  IM=-1
		ELSE
		  CALC=0
		ENDIF
	      ELSEIF(CALC.EQ.1) THEN
		IM=IM+1
	      ELSE
		IM=IM-1
	      ENDIF
	    ENDIF
	    IIM=IABS(IM)
	    IF(IIM.EQ.NL) THEN
	     NL=NL+1                        ! calc. of wave function
	     IF(NL.GT.NLmax) THEN
	       NL=NLmax
*! begin fragment: commented and new lines:
d	WRITE(*,47) IT,IQ,NL,ISH,JSH,M,N,NP,M1,N1,NP1,IM,TMP,TMPH,Eps
c	WRITE(16,47) IT,IQ,NL,ISH,JSH,M,N,NP,M1,N1,NP1,IM,TMP,TMPH,Eps
d	write(*,*)'*WRN* Accuracy loss; point, row, col:',NPNT+1,ISH,JSH
	write(16,'(A44,2i6)')
     +	'*WRN* Accuracy loss in matrix elements, i,j:',ish,jsh
d	CALL PROGRES(NProgPNT)	! redraw progres bar and points on it
*! end fragment
	       IF(CALC.EQ.1) THEN	! IM > NLmax not allowed, ==>,
		CALC=2			! change calc. direction or exit.
		IM=-1
               ELSE
		CALC=0
	       ENDIF
	     ELSE
	       DO I=0, IIM	! additional I1 and I2 for new NL
		SI1(I+1,NL+1) =
     .			SIMP_I12(Eps,RMTS(IT),VKP(1),KPNM(1), I,NL,1)
		SI1_(I+1,NL+1)=
     .			SIMP_I12(Eps,RMTS(IT),VKP(2),KPNM(2), I,NL,1)
		SI2(I+1,NL+1) =
     .			SIMP_I12(Eps,RMTS(IT),VKP(1),KPNM(1), I,NL,2)
		SI2_(I+1,NL+1)=
     .			SIMP_I12(Eps,RMTS(IT),VKP(2),KPNM(2), I,NL,2)
	       END DO
	     ENDIF
	    ENDIF
	    BSmCALC(1,IRO(IT,IQ))=BSmCALC(1,IRO(IT,IQ))+TMP
	    BSmCALC(2,IRO(IT,IQ))=BSmCALC(2,IRO(IT,IQ))+TMPH
	    BSmCALC(3+NT*NLsav,IRO(IT,IQ))=
     .				BSmCALC(3+NT*NLsav,IRO(IT,IQ))+TMPIS
	   END DO
*
c	-- Here we try to calculate (-1)**(M+M') --
	   IF(IABS(MOD(M+M1,2)).EQ.1) THEN
	     DO I=1, 3+NT*NLsav
	       BSmCALC(I,IRO(IT,IQ))=-BSmCALC(I,IRO(IT,IQ))
	     END DO
	   ENDIF
	   LSmCALC(IRO(IT,IQ))=.TRUE.
	  ENDIF
*
*! old fragment:
c	  TMP=2.D0*PI*DBLE(NP-NP1)/CZ*RSZ(IQ,IT)+DBLE(M-M1)*PHI
*! new fragment SP 3.01:
	  TMP=(VKP(1)-VKP(2))*RSZ(IQ,IT)+DBLE(M-M1)*PHI
*! end fragment
	  S_RE=S_RE+DCOS(TMP)*BSmCALC(1,IRO(IT,IQ))
	  S_IMG=S_IMG+DSIN(TMP)*BSmCALC(1,IRO(IT,IQ))
	  SH_RE=SH_RE+DCOS(TMP)*BSmCALC(2,IRO(IT,IQ))
	  SH_IMG=SH_IMG+DSIN(TMP)*BSmCALC(2,IRO(IT,IQ))
*
	  DO I=(IT-1)*NLsav+1, IT*NLsav	! in range of this value of atom type
	    PoolTmp(ISH,JSH,I)=PoolTmp(ISH,JSH,I)
     .		+DCMPLX(DCOS(TMP)*BSmCALC(2+I,IRO(IT,IQ)),
     .			DSIN(TMP)*BSmCALC(2+I,IRO(IT,IQ)))
	  END DO		! "1+NT*NLsav" for calc. of Qis { }-part
	  PoolTmp(ISH,JSH,1+NT*NLsav)=PoolTmp(ISH,JSH,1+NT*NLsav)
     .		+DCMPLX(DCOS(TMP)*BSmCALC(3+NT*NLsav,IRO(IT,IQ)),
     .			DSIN(TMP)*BSmCALC(3+NT*NLsav,IRO(IT,IQ)))
    2	CONTINUE

*-------------------------------------------------------------------------

c -- On first point - save the integrals, otherwise - can use saved, if any --
	IF(SaveI3M.AND.Create_I3M .OR. Add_to_Saved_I3M) THEN
	  II=0			! Saving of integrals here
	  K=0
	  DO I=1, 2*NLmax
	    IF(LCALC(I)) THEN
	      IF(K.EQ.0) K=I     ! first element in array
	      II=II+1            ! count of elements
	    ENDIF
	  END DO
	  IF(II.GT.0) THEN
	    WRITE(15,REC=MOV,ERR=85) II,K,
     .		(BDCALC(I),BECALC(I),BFCALC(I),I=K,K+II-1)
	  ELSE
	    WRITE(15,REC=MOV,ERR=85) 0,0,
     .		(BDCALC(I),BECALC(I),BFCALC(I),I=1,2*NLmax)
	  ENDIF
	  GO TO 86
   85	  WRITE(*,83) MOV,II,K+II,2*NLmax,'write.'
	  STOP '*ERR* An internal error occured!'
   86	  IF(Add_to_Saved_I3M) THEN
c	    WRITE (*,14) NPNT+1,IM,M,N,NP,M1,N1,NP1
cc	    WRITE(16,19) NPNT+1,NPNT-NCPNT+1,NPNTS,NTPNTS
c	    WRITE(16,14) NPNT+1,IM,M,N,NP,M1,N1,NP1
c	    CALL PROGRES(NProgPNT)     ! redraw progres bar and points on it
	    Add_to_Saved_I3M=.FALSE.
	  ENDIF
	ENDIF

	DO I=(IT-1)*NLsav+1, IT*NLsav         ! only for Ql(alpha) {}-parts.
	PoolTmp(ISH,JSH,I)=PoolTmp(ISH,JSH,I)*DCMPLX(RMTS(IT)**4.D0,0.D0)
	END DO

	MOV=MOV+1   ! next record for read/write to file of double integrals
    1	CONTINUE

*! old fragment:	TMP=S_NJM(KPNM(1),KPNM(2),M,M1,A)/CZ
*! begin fragment
	if (iUseCoredTube.eq.0) then
	  TMP=S_NJM(KPNM(1),KPNM(2),M,M1,A)/CZ
*! new fragment SP 3.1
	else if (iUseCoredTube.eq.2) then
	  YEps=zf_params.eps
	  TMP=S_NYM(KPNM(1),KPNM(2),M,M1,B,YEps)/CZ
*! end fragment SP 3.1
	else if (iUseCoredTube.eq.1) then
	  TMP=S_CJY(A,B,M,N,M1,N1,MMAX,NMAX,abJ,abJD,abY,abYD)/CZ
*! new fragment SP 6.0
	else if (iUseCoredTube.eq.3) then
	  TMP=S_CJYCR(A,B,M,N,M1,N1,MMAX,NMAX,abJ,abJD,abY,abYD,
     *  AllKappas,V,IntCR,IntCRC)/CZ
	else if (iUseCoredTube.eq.4) then
	  TMP=S_CJYCRI(A,B,M,N,M1,N1,MMAX,NMAX,abJ,abJD,abY,abYD,
     *  AllKappas,V,IntCR,IntCRC)/CZ
*! end fragment SP 6.0
	endif
d	write(16,*) 'S_XXX>> ',tmp,kpnm(1),kpnm(2)
*! end fragment

	TM1=-TMP*S_RE		! real part of overlap integral matrix
	TM2=-TMP*SH_RE		! real part of hamiltonian matrix
	PoolTmp(ISH,JSH,1+NT*NLsav)=DCMPLX(-TMP,0.d0)
     .					*PoolTmp(ISH,JSH,1+NT*NLsav)
c	-- Here we add a delta-function --
	IF(NP.EQ.NP1.AND.M.EQ.M1.AND.N.EQ.N1) THEN
	  TM1=1.D0+TM1
	  TM2=(VKP(1)*VKP(1)+KPNM(1)*KPNM(1)+V_out)+TM2
*	  TM2=(VKP(1)*VKP(1)+KPNM(1)*KPNM(1))+TM2
	  PoolTmp(ISH,JSH,1+NT*NLsav)=1.D0+PoolTmp(ISH,JSH,1+NT*NLsav)
	ENDIF

	S(ISH,JSH)=DCMPLX(TM1,-TMP*S_IMG)
	H(ISH,JSH)=DCMPLX(TM2,-TMP*SH_IMG)

*! SP 4a: add the matrix element of full potential:
	if (FP.iUseFP.ne.0) H(ISH,JSH)=H(ISH,JSH)+VfpMEs(iVfpMEs)*TMP

	IF(ISH.NE.JSH) THEN		! Mirroring of the elements
	  S(JSH,ISH)=DCONJG(S(ISH,JSH))
	  H(JSH,ISH)=DCONJG(H(ISH,JSH))
	ENDIF
	DO I=1, NT*NLsav
	  PoolTmp(ISH,JSH,I)=2.D0*TMP*PoolTmp(ISH,JSH,I)
	  IF(ISH.NE.JSH) THEN
	    PoolTmp(JSH,ISH,I)=DCONJG(PoolTmp(ISH,JSH,I))
	  ENDIF
	END DO
	PoolTmp(ISH,JSH,1+NT*NLsav)=2.D0*PoolTmp(ISH,JSH,1+NT*NLsav)
	IF(ISH.NE.JSH) THEN
	  PoolTmp(JSH,ISH,1+NT*NLsav)=DCONJG(PoolTmp(ISH,JSH,1+NT*NLsav))
	ENDIF
*! begin fragment SP 4a:
	if (FP.iUseFP.ne.0) then
	IF(NPNT.EQ.0) WRITE(11,111) JSH,M,N,NP,M1,N1,NP1,
     .		DBLE(S(ISH,JSH)),DIMAG(S(ISH,JSH)),
     .		DBLE(H(ISH,JSH)),DIMAG(H(ISH,JSH)),
     .		DBLE(VfpMEs(iVfpMEs)*TMP),DIMAG(VfpMEs(iVfpMEs)*TMP)
	else
*! old fragment:
	IF(NPNT.EQ.0) WRITE(11,11) JSH,M,N,NP,M1,N1,NP1,
     .		DBLE(S(ISH,JSH)),DIMAG(S(ISH,JSH)),
     .		DBLE(H(ISH,JSH)),DIMAG(H(ISH,JSH)) ! for test.
*! end fragment
	endif
*! end fragment SP 4a
	IF(ITRX.EQ.1) THEN
c	-- ATTENTION: this equation is correct only for upper-triang. matrix ! --
	  I=NCXPOINT-(COUNTX-ISH)*(COUNTX-ISH-1)/2-2*COUNTX+ISH+JSH ! cur. point
	  IF( NProgPNT.LT.NPROGRES .AND. I.EQ.NextPnt ) THEN
	    CALL PROGRES(-1)   ! plot one point
	    NextPnt=NextPnt+(NCXPOINT-I+1)/(NPROGRES-NProgPNT+1)
	    IF(NProgPNT.LT.NPROGRES) NProgPNT=NProgPNT+1
	  ENDIF
	ENDIF
    7	CONTINUE

	if (NPNT.eq.0) close(11)         ! reserved for test - log-file.
*
	IF(ITRX.EQ.1) WRITE(*,'(A1)') '<'         ! progress bar complete.
	IF(NL.GT.NLsav) WRITE(*,50) NL,NLmax
*
*! begin fragment SP 4a:
* Save the FP matrix elements if they were re-calculated with this k-point:
	if (FP.iUseFP.ne.0.and.iVfpMatrixRequest.ne.0) then
	  close(14)	! bwfmaps.dat
	  deallocate(VfpMap,cBwfMap,xset,yset,zset,int_map)
	  iVfpMatrixRequest=0
	  open(17,file='fpmatrix.dat',form='unformatted')
	  write(17) COUNTX
	  write(17) (VfpMEs(i),i=1,COUNTX*(COUNTX+1)/2)
	  write(17) VAtShft
	  close(17)
	  write(16,*) 'FP: Request for fpmatrix.dat is processed'
	  write(*,*) 'FP: Request for fpmatrix.dat is processed'
	endif
*! end fragment SP 4a
*
      CALL CHLSK1C(COUNTX,COUNTX,S,Diag)
      CALL CHLSK2C(COUNTX,COUNTX,H,S,Diag)
      CALL CHLSK3C(COUNTX,COUNTX,H,S,Diag)
      WRITE(*,'(/1X,A29)') 'Cholesky parts done.'
*
      NEigv=COUNTX
      info=0
c --- Here we call NAG/LAPACK routine for eigenvalues and eigenvectors ---
      CALL F02HAF(JOB,UPLO,COUNTX,H,COUNTX,D,RWORK,WORK,LWORK,INFO)
c    -- after this routine, D contain eigenvalues,         --
c    -- H - corresponding eigenvectors in orthogonal basis set --
*
*      CALL ZHEEVX (JOB, RANG, UPLO, COUNTX, H, COUNTX, VL, VU, 1, IBNDS,
*     +            0.D0, NEigv, D, Vec, COUNTX, WORK, 64*MUL*MAXMNP,
*     +            RWORK, IWORK, IFAIL, INFO)
*
      IF(INFO) THEN
        WRITE(*,36)
        STOP '---> DIAG: matrix elements are bad!'
      ENDIF
      print 37, NEigv, 'eigenvalues has been found,'
c   -- reducing of eigenvectors (stored in H) to source basis set --
      CALL REIGNR_c(COUNTX,COUNTX,NEigv,H,S,Diag)
! 12222
*   -- Test output. H contain eigenvectors in source basis set. --
d      if(npnt.eq.0) then
d        write(12) ((H(i,j),j=1,countx),i=1,countx)
d        write(12) (D(i),i=1,countx)
d        close(12)
d      endif
*   -- end of debug part --
*
	IF(NPNT.EQ.0) THEN 
		OPEN(121, FILE='LENH_LEND.txt')
		INQUIRE(IOLENGTH = LENH) H
		INQUIRE(IOLENGTH = LEND) D
		WRITE(121, *) LENH, LEND
		CLOSE(121)
	END IF

	OPEN(121, FILE = 'H_AMNP.TXT', ACCESS = 'DIRECT',
	+ RECL = LENH+1000)
C	ÇÀÏÈÑÛÂÀÅÌ ÇÍÀ×ÅÍÈß H Ñ ÏÎÐßÄÊÎÂÛÌ ÍÎÌÅÐÎÌ NPNT Â ÔÀÉËÅ, ÃÄÅ NPNT - ÝÒÎ ÍÎÌÅÐ ÒÎ×ÊÈ K(ÂÎËÍÎÂÎÃÎ ÂÅÊÒÎÐÀ)
	WRITE(121, REC = NPNT+1) ((H(I,J),J=1,COUNTX),I=1,COUNTX)
	CLOSE(121)

	OPEN(121, FILE = 'D_ENERGIES.TXT', 
	+ACCESS = 'DIRECT', RECL = LEND+1000)
	WRITE(121, REC = NPNT+1) D
	CLOSE(121)

	IF(NPNT+1.GT.MAXPNT) THEN
	  WRITE(*,'(/1X,A37,I3)') '*ERR* Too many points. Max. value is ',
     .		MAXPNT
	  STOP
	ENDIF
	DO I=1, NEigv
c	  EBND(I,NPNT+1)=D(I)+VAtShft     ! in Rydbergs.
	  EBND(I,NPNT+1)=D(I)     ! in Rydbergs.
	END DO
	IF(NEigv.GT.MAXnEigv) MAXnEigv=NEigv

c   -- Calc. of partial charges, based on PoolTmp and H - eigenvectors --
      IQ=0	! here: a simple flag.
      PRINT '(1X,A\)','Building partial charges...'
      DO II=1, NT*NLsav+1
	TMP=0.d0
	DO I=1,MAXnEigv
	  TMPc=DCMPLX(0.d0,0.d0)
	  DO J=1,MAXnEigv
	    TMPcc=DCONJG(H(I,J))
	    DO K=1,MAXnEigv
	      TMPc=TMPc+TMPcc*H(K,J)*PoolTmp(I,K,II)
	    END DO
	  END DO
	  PoolRes(I,NPNT+1,II)=DBLE(TMPc)
	  TMP=TMP+DABS(DIMAG(TMPc)) ! Checking of the imag. part of partial charges
	END DO
*! begin ? fragment
	IF(DABS(TMP).GT.Eps) THEN
	  IF(II/NLsav.LT.10) THEN
	    IF(IQ.eq.0) PRINT '(/)'
*! begin fragment SP 5.08. Bug out or bug in?
	    if (II.eq.NT*NLsav+1) then
              PRINT 90, 'X','x',DABS(TMP)
	      WRITE(16,90) 'X','x',DABS(TMP)
	    else
	      IT=(II-1)/NLsav+1
	      K=MOD(II-1,NLsav)+1
	      WRITE(*,90) CHAR(IT+48),QchText(K),DABS(TMP)
	      WRITE(16,90) CHAR(IT+48),QchText(K),DABS(TMP)
	    endif
*! old fragment:
c            PRINT 90, TNMBR((II-1)/NLsav+2),
c     .		QchText(II-(II/NLsav)*NLsav),DABS(TMP)
c	    WRITE(16,90) TNMBR((II-1)/NLSav+2),
c     .		QchText(II-(II/NLsav)*NLsav),DABS(TMP)
*! end fragment SP 5.08
	  ELSE
            PRINT 90, TNMBR(((II-1)/NLsav+1)/10+2)//TNMBR(II/NLsav+2),
     .		QchText(II-(II/NLsav)*NLsav),DABS(TMP)
	    WRITE(16,90) TNMBR(((II-1)/NLsav+1)/10+2)//TNMBR(II/NLsav+2),
     .		QchText(II-(II/NLsav)*NLsav),DABS(TMP)
	  ENDIF
	  IF(IQ.eq.0) THEN
d??!!!	    PRINT '(/)'
	    IQ=1
	  ENDIF
	ENDIF
*! end ? fragment
d        print *,ii,' imaginary part of partial charges',tmp
c!!!        IF(DABS(TMP).GT.Eps) THEN
c          PRINT '(//1X,A16,i2,A37,F8.5)', '*ERR* Atom type ',
c     . (ii-1)/NLsav+1,'L =',mod(ii-1,NLsav),
c     . ': Imaginary part of partial charges:',tmp
c         STOP
c       ENDIF
      END DO
      PRINT 70
      IF(NPNT.EQ.0) THEN
        WRITE(16,'(/1X,A18\)') 'Eigenvalues for k='
        IF(TBZERO) THEN
          WRITE(16,'(A2)') '0:'
         ELSE
          WRITE(16,'(A6)') '-PI/C:'
        ENDIF
        DO i=1,NEigv    ! write eigenvalues to log-file
          WRITE(16,'(1X,I3,3X,2(G14.4,2X))') i,D(i),D(i)*EVS
        END DO
      ENDIF
*      if(npnt.eq.ncpnt+npnts-1) then
*       write(16,*) 'On point',npnt+1
*       write(16,*) 'The values of partial charges are:'
*       DO II=1, NT*NLsav+1
*         write(16,*) 'Charge series',II
*         DO I=1,NEigv
*          WRITE(16,'(1X,i3,1x,f14.4)') i,PoolRes(I,npnt+1,II)
*         END DO
*         write(16,'(1x/)')
*       END DO
*      endif
      IF(Create_I3M .AND. SaveI3M)
     .   WRITE(15,REC=1) Complete,TJNM
      Create_I3M=.FALSE.

      PRINT '(1X,A\)', 'Writing bands data to file...'
c   -- Save calculated data to file BND.STR --
      OPEN(9,FILE='bnd.str',FORM='UNFORMATTED',STATUS='OLD',ERR=75)
      CLOSE(9,STATUS='DELETE')  ! Replace existing data, stored in EBND()
   75 OPEN(9,FILE='bnd.str',FORM='UNFORMATTED',STATUS='NEW')    ! for first point
      WRITE(9) TJNM,NTPNTS,NPNT+1,NLmax,MAXnEigv
      WRITE(9) ((EBND(I,J),J=1,NPNT+1),I=1,MAXnEigv)
      DO K=1, NT*NLsav+1
        DO I=1, MAXnEigv
          WRITE(9) (PoolRes(I,J,K),J=1,NPNT+1)
        END DO
      END DO
      CLOSE(9)
      WRITE(*,70)
*! old line:      CALL TTIME(TIME,TATIME,0)
*! new line:
      CALL TTIME(TIME,TATIME,16)
      PRINT '(/)'
   77 CONTINUE          ! end of cycle on points.
*
*! added de-allocation code for FP arrays in SP 4a:
c	deallocate(VfpMap,VfpMEs,cBwfMap)
	if (FP.iUseFP.ne.0) deallocate(VfpMEs)
*
*! added de-allocation code for the arrays which were static before SP 5.04:
	deallocate(EBND)
	deallocate(D,Diag,RWORK,WORK)


	DEALLOCATE (H,S,PoolTmp,PoolRes)

	NCPNT=NCPNT+NPNTS
	IF(NCPNT.LT.NTPNTS) WRITE(*,32) NTPNTS-NCPNT

	CALL TTIME(TIME_of_all,TCTIME,16)

	CLOSE(16)         ! End of log-file.

	IF(NCPNT.LT.NTPNTS) WRITE(*,34)
*
c   11 FORMAT(1X,I3,1X,2(3I3,1X),2x,2(G12.4,1x),1x,2(G12.4,1x))
   11 FORMAT(1X,I3,1X,2(3I3,1X),2(2x,2G12.4))
  111 FORMAT(1X,I3,1X,2(3I3,1X),3(2x,2G12.4))
   12 FORMAT(/,1X,'On this point, #',I3,', additional data has been
     + saved for summs, based on I1 and I2.')
*! begin new fragment
   14 FORMAT(1X,'On this point, #',I3,', additional data has been',
     +' saved for integrals I3, I3'' and Im for L=',I3,
     +' M,N,P x M'',N'',P'' ',3I3,' x',3I3,'.')
*! end new fragment
*! begin old fragment
c   14 FORMAT(/1X,'On this point, #',I3,', additional data has been',\,
c     +1X,'saved for integrals I3,',/,' I3'' and Im for L=',I3,
c     +' M,N,P x M'',N'',P'' ',3I3,' x',3I3,'.')
*! end old fragment
   16 FORMAT(1X,'---> Your basis,',I5,' functions, is too large. Limit
     + is', i5,'.')
   17 FORMAT(F5.3,1X,F9.3)
   18 FORMAT(28X,'< Iteration No ',I3,' >')
   19 FORMAT(1X,'Point ',I3,20X,'[',2(I3,'/'),I3'+1 ]')
   20 FORMAT(/2X,A50,I4)
*   21 FORMAT(/,<K>X,'<<< ',A<J-I>,' >>>',//)
   22 FORMAT(2F9.3,I5,2G10.4,5I5)
   23 FORMAT(10F10.3)
   25 FORMAT(/1X,'Your count of band parameters is less than as your',
     +' L (found as ',I1,'),',/,' or string of band parameters is bad.')
   26 FORMAT(1X,'Continuing of calculation. Old data has been read,',/,
     +1X,I3,' point(s) calculated, ',I3,' - calculating at this time.')
   27 FORMAT(1X,'WARNING: your count of points, ',I3,' truncated to',
     +I3)
   29 FORMAT(/1X,
     +'Your total count of points in your bands file and in main',/,
     + 1X,'input file are different, ',I3,' and ',I3)
   32 FORMAT(1X,'End in this time, ',I3,
     +       ' points left to be calculating.',/)
   33 FORMAT(/1X,'Your basis set, ',I3,', functions, is too small.',/,
     +       'Must be ',I3,' or greather.')
   34 FORMAT(/1X,'To be continued...')
   36 FORMAT(/1X,'An error occurs in eigenvalues!')
   37 FORMAT(1X,I4,A28)
   38 FORMAT(/1X,'Your count of bands ',I3,', is too small,',/,1X,
     +       'must be ',I3,' or greather.')
   39 FORMAT(1X,'Basis has ',I3,' functions.')
   42 FORMAT(/1X,'Your count of bands, ',I3,', is greather than',/,1X,
     +       'your basis set, ',I3,' functions.')
   44 FORMAT(//1X,'Shift of potential is ',F10.5,' eV')
   47 FORMAT(/1X,'RECOMMENDATION - increase your count of band ',
     +     'parameters: on atom type ',I2,' No ',I2,/,
     +      16X,I2,' is too small, your results may be incorrect.'/,
     +      7X,'LOCATION: Row ',I3, ' Column ',I3,
     +     ': M,N,P x M'',N'',P'' ',3I3,' x',3I3,', m=',I3,/6X,
     +     'Last elements of sum are:',G12.4,',',G12.4'; Eps=',G12.4,/)
   48 FORMAT(1X,'ATTENTION! You have changed max. bound of L. Can''t ',
     +       'continue with ',\,'new value.',/,' Set to ',I2,' now and',
     +     \,' restart or begin from first point with your new value.')
   49 FORMAT(/1X,'Fermi Level:',F8.3,' eV.')
   50 FORMAT(6X,'Max. value for L on this point is ',I2,', limit: ',
     +      I2,\,'.',/)
   53 FORMAT(5X,'For atom type ',I2,' - '\)
   54 FORMAT(3X,I2,1X,F11.6,6X,I2,8X,I2)
   56 FORMAT(I2,' atom(s) of ',I2,' with different distances to Z-axis',
     +      /)
   57 FORMAT(3X,I2,1X,F11.6,6X,I2)
   58 FORMAT(1X, '+----+------------+---------+---------+')
   61 FORMAT(1X,'No value in field of precision for calc. of ',\
     +     'integrals, ',F10.5,' by default.')
   62 FORMAT(1X,'No value in field of precision of coordinates',\
     +      F10.5,' by default.')
   64 FORMAT(I3)
   70 FORMAT(' done.')
   83 FORMAT(//1X,'*ERR* Random access overflow in operation with ',\
     +      'saved integrals -'/,7X,'Record:',I7/,
     +      7X,'Start element:',I2/,7X,'Last element:',I2/,
     +      7X,'Allowed value for number of element:',I2/,
     +      7X,'Action: ',A6)
   90 FORMAT(1X,'*WRN* Atom type ',A2,', state ',A1,': ',!'-',7X,
     +  'Sum of modules of imaginary parts of partial charges:',F10.7)
	END
*
	DOUBLE PRECISION FUNCTION FCT(L,M)
c	-- this function calc. (L+M)!/(L-M)! --
	INTEGER I,M,L
	FCT=1.D0
	IF(M.NE.0) THEN
	  DO I=L-M+1, L+M
	    FCT=FCT*DBLE(I)
	  END DO
	ENDIF
	RETURN
	END
*
	SUBROUTINE PROGRES(NPT)
c ---- Draw progress bar on screen ----
c	NPT = -1 ?  If yes - draw one point, no - redraw progres bar.
	INTEGER I,NPT
	IF(NPT.NE.-1) THEN
	  WRITE(*,1)
	  DO I=1,5
	    WRITE(*,2) '+'
	  END DO
	  WRITE(*,'(/,6X,A1,\)') '>'
	  IF(NPT.NE.0) THEN
	    DO I=1,NPT
	      WRITE(*,4) '.'
	    END DO
	  ENDIF
	ELSE
	  WRITE(*,4) '.'
	ENDIF
    1	FORMAT(7X,'0%',10X,'20%',10X,'40%',10X,'60%',10X,'80%',9X,
     +		'100%',/7X,'+'\)
    2	FORMAT(12('-'),A1\)
    3	FORMAT(A1,/7X)
    4	FORMAT(A1\)

	RETURN
	END
*
	SUBROUTINE GAUSS(IFL,Pool,N,Indx,NDIM,EL,ER,HFWDT,HFWDTEF,
     +			EF,IDNST,XNORM,WRK,DosSave,NLs,NPool,outf)
c    -- EF: on entry, contain the energy on Fermi Level,
c    -- WRK and DosSave are temporary matrices for building of DOS's,
c    -- on exit - for both vectors, first elem. contained full DOS
c    -- for HFWDT in DosSave and for HFWDTEF - in WRK.
c    -- Simultaneously, next N elem. containing DOS on Fermi level
c    -- for simultaneous partial DOS's in order of Pool(NPool,*).
	INTEGER Indx(2,*),NLs,IFL,IDNST,NPTS,NDIM,ISAVE,I,J,K,II,Ix,
     +		NPool,TPOINTS
	DOUBLE PRECISION PI,Pool(NPool,*),DDE,EL,ER,ED,EF,HFWDT,HFWDTEF,
     +		XNORM,WRK(*),DosSave(*),TMP
	CHARACTER outf*(*)(*)
	PARAMETER( TPOINTS=70 )
	PI=4.D0*DATAN(1.d0)
*
	NPTS=IDNINT((ER-EL)*DBLE(IDNST))	! Number of points
	DDE=(ER-EL)/DBLE(NPTS)
	ISAVE=IDINT((EF-EL)/DDE)
	EL=EF-DBLE(ISAVE)*DDE	! correct EL for DOS(EF) to be calculated
	IF(ER.GT.Pool(1,NDIM)-5.d0*HFWDT) THEN
	  WRITE(*,2) ER
	  ER=Pool(1,NDIM)-5.d0*HFWDT
	  WRITE(*,3) ER
	ENDIF
	DO I=0, N
	  OPEN(IFL+i,FILE=outf(i+1),FORM='FORMATTED',ERR=5)
	END DO
	J=0
	II=0
	DO K=1, NPTS+1
	  DO I=1, N+1
	    WRK(i)=0.D0
	  END DO
	  ENER=EL+DDE*(K-1)
	  DO L=1, NDIM
	    ED=ENER-Pool(1,L)
	    TMP=DEXP(-ED*ED/(HFWDT*HFWDT*2.D0))
	    WRK(1)=WRK(1)+TMP	! Full DOS is calculated here.
	    DO i=2, N+1
	      IF(Indx(1,i-1).GT.0) THEN
		Ix=(Indx(1,i-1)-1)*NLs+Indx(2,i-1)+1	! same as (IT-1)*NLsav+L+1
	      ELSE
		Ix=NPool-1	! for case, when intersphera DOS selected.
	      ENDIF
	      WRK(i)=WRK(i)+Pool(1+Ix,L)*TMP	! partial DOS's.
	    END DO
	  END DO
	  DO i=1, N+1
	    WRK(i)=WRK(i)/HFWDT/DSQRT(2.0D0*PI)*XNORM
	    IF(K.EQ.ISAVE+1) DosSave(i)=WRK(i)
	  END DO
	  IF(II.EQ.(NPTS+1)/TPOINTS) THEN
	    II=0
	    PRINT '(A1\)', '.'	! Draw simple progres indicator, when Gauss procedure works.
	  ELSE
	    II=II+1
	  ENDIF
	  IF(ENER.GE.EL+5.d0*HFWDT .AND. ENER.LE.ER-5.d0*HFWDT) THEN
	    DO i=1, N+1		! Write out to opened file(s)
	      WRITE (IFL+i-1,700) ENER, WRK(i)
	    END DO
	    J=J+1
	  ENDIF
	END DO
	WRITE(*,701) J
c	-- calculating of full DOS and partial DOS's on Fermi level --
	DO i=1, N+1
	  WRK(i)=0.d0
	END DO
	DO L=1, NDIM
	  ED=EL+DDE*DBLE(ISAVE)-Pool(1,L)
	  TMP=DEXP(-ED*ED/(HFWDTEF*HFWDTEF*2.D0))
	  WRK(1)=WRK(1)+TMP
	  DO i=2, N+1
	    IF(Indx(1,i-1).GT.0) THEN
	      Ix=(Indx(1,i-1)-1)*NLs+Indx(2,i-1)+1
	    ELSE
	      Ix=NPool-1
	    ENDIF
	    WRK(i)=WRK(i)+Pool(1+Ix,L)*TMP
	  END DO
	END DO
	DO i=1, N+1
	  WRK(i)=WRK(i)/HFWDTEF/DSQRT(2.0D0*PI)*XNORM
	END DO
	DO i=0, N	! Close all files after writing.
	  CLOSE(IFL+i)
	END DO
	GO TO 7
    5	PRINT '(//1X,A33,i2)','*ERR* Can''t open file for point #',i+1
	STOP
    7	RETURN
    1	FORMAT(1X,'*WRN* Your left bound of energy interval,',F7.2,'eV,'\,
     .		' increased '\)
    2	FORMAT(1X,'*WRN* Your right bound of energy interval,',F7.2,'eV,',
     .		\' reduced '\)
    3	FORMAT('to ',F7.2,'eV,',/7X,'by using the procedure of Gaussian'\,
     .		' broadening.'/)
 700	FORMAT(F8.3,1X,F9.3)
 701	FORMAT(//1X,I4,' point(s) written.'/)
 702	FORMAT(/' *ERR* Count of points in your density of points file,'\,
     .	1X,I4,',',/7X,'is greather than limit, ',I4,'.')
	END
*
	SUBROUTINE swap(II,Indx,Pool,N,IEnd)
c	-- rearranging of Pool by indexes Indx --
	INTEGER i,II,Indx(*),N,IEnd
	DOUBLE PRECISION Pool(N+1,*),TEMP
	DO i=1,N+1
	  TEMP=Pool(i,II)
	  Pool(i,II)=Pool(i,Indx(II))
	  Pool(i,Indx(II))=TEMP
	END DO
	i=II
	DO WHILE(Indx(i).NE.II)
	  IF(i.EQ.IEnd) THEN
	    PRINT '(//1x,A37,I5)',
     .		'*ERR* Can''t find index for element #',II
	    STOP
	  ENDIF
	  i=i+1
	END DO
	Indx(i)=Indx(II)
	Indx(II)=II
	RETURN
	END
