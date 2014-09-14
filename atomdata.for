	integer function iReadAtomic(fn)
	implicit double precision(a-h,o-z)

	include 'atomdata.fh'

	character*(*) fn

	OPEN(4,FILE=fn,FORM='UNFORMATTED',STATUS='OLD',ERR=63)
	print *,'nt=',nt
	DO IT=1,NT
	READ(4,ERR=201,END=202) RSTART(IT),STEP(IT),TMP,JRIS(IT),J,K
	READ(4,ERR=201,END=202) (CORE(I,IT),I=1,J)
	READ(4,ERR=201,END=202) (AMTC(I,IT),I=1,J)
	READ(4,ERR=201,END=202) (RAD(I,IT),I=1,J)
	NATOM(IT)=J
	STRJ(IT)=AMTC(JRIS(IT),IT)
c	-- Skipping this part - not need more. --
	READ(4)
	READ(4)
	READ(4)
	READ(4)
	READ(4)
	END DO
	CLOSE(4)

* Calculate effective charges for atom types:
	if (logf.ne.0) write(logf,*) 'Effective charges:'
	open(121, file='JRIS_RAD.txt')
	do it=1,nt
	rmts(it)=rad(jris(it),it)	!*0.9d0
	write(121, *) jris(it), rmts(it)
	write(121, *) rad(:, it)
	qeff(it)=dble(zmain(it))
     +	-ATINT_(amtc,rad,jris(it),it,step(it),NGRID,NCUT,1,1.d0)
	if (logf.ne.0) write(logf,*) it,zmain(it),qeff(it)
	enddo
	close(121)

	if (logf.ne.0) then	! verify total charges
	rho=0.0d0
	write(logf,*) 'Total charges:'
	do it=1,nt
	  tmp=ATINT_(amtc,rad,NGRID,it,step(it),NGRID,NCUT,1,1.d0)
	  write(logf,*) it,zmain(it),tmp
	  rho=rho+tmp*nts(it)
	enddo
	write(logf,*) 'Total charge in cell:',rho
	endif

	iReadAtomic=1
	return

! failure stops:
   63	continue	!stop 'RDF: Failed to open OutAtm.str. Re-run Atom'
  201	continue	!stop 'RDF: Error reading OutAtm.str'
  202	continue	!stop 'RDF: Unexpected end of OutAtm.str'
	iReadAtomic=0
	return

	end
