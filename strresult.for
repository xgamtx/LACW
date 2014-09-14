	integer function iReadStrResult(f,Basis)	!,Nact,AllKappas,abJ,abJD,abY,abYD)
	implicit double precision(a-h,o-z)

	include 'atomdata.fh'
	include 'strresult.fh'
	record /TBasis/ Basis
	integer f
	integer n,istat

	include 'stdio.fh'

	read(f) n
        Basis.COUNTX=n
        if (associated(Basis.INDXX)) deallocate(Basis.INDXX)
	allocate(Basis.INDXX(4,n),stat=istat)
	read(f) ((Basis.INDXX(j,i),j=1,3),i=1,n)
	do it=1,nt
c	  read(f) (Wpot(i,it),i=1,jris(it))
	  read(f) (Wpot(i,it),i=1,NGRID)
	enddo

	if (DEBUG.ge.9.and.ierr.ne.0) then	!`````
	write(ierr,*) 'iReadStrResult: Basis: (M,N,P):'
	do i=1,n
	  write(ierr,*) (Basis.INDXX(j,i),j=1,3)
	enddo
	write(ierr,*) 'iReadStrResult: Wpot: (it,(i,rad,Wpot)):'
	do it=1,nt
	  write(ierr,*) 'Atom type',it
	  do i=1,NGRID
	    write(ierr,*) i,rad(i,it),Wpot(i,it)
	  enddo
	enddo
	endif	!`````

	read(f) (NAct(m),m=1,MMAX)
	read(f) ((AllKappas(m,n),n=1,NAct(m)),m=1,MMAX)
	read(f) (((abJ(m,n,i),abJD(m,n,i),abY(m,n,i),abYD(m,n,i),
     +			i=1,2),n=1,NAct(m)),m=1,MMAX)

d	do m=1,MMAX
d	do n=1,NAct(m)
d	  write(17,'(2i4,9f20.10)') m,n,AllKappas(m,n),
d     +		(abJ(m,n,i),abJD(m,n,i),abY(m,n,i),abYD(m,n,i),i=1,2)
d	enddo
d	enddo

	iReadStrResult=1
	return
	end
*-------------------------------------------------------------------------
	integer function iWriteStrResult(f,Basis)	!,Nact,AllKappas,abJ,abJD,abY,abYD)
	implicit double precision(a-h,o-z)

	include 'atomdata.fh'
	include 'strresult.fh'
	record /TBasis/ Basis
	integer f

	if (.not.associated(Basis.INDXX)) then
	  iWriteStrResult=0
	  return
	endif

	write(f) Basis.COUNTX
	write(f) ((Basis.INDXX(j,i),j=1,3),i=1,Basis.COUNTX)
	do it=1,nt
	  write(f) (Wpot(i,it),i=1,NGRID)	!i=1,jris(it))
	enddo

	write(f) (NAct(m),m=1,MMAX)
	write(f) ((AllKappas(m,n),n=1,NAct(m)),m=1,MMAX)
	write(f) (((abJ(m,n,i),abJD(m,n,i),abY(m,n,i),abYD(m,n,i),
     +			i=1,2),n=1,NAct(m)),m=1,MMAX)

	iWriteStrResult=1
	return
	end
