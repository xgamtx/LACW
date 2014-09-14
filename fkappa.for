*-------------------------------------------------------------------------
* This module contains two additional functions for cored LACW package.
*
* Exports:
*    double zerofunc(double x);
*    void FindKappa(double **Dest, int Mmax, int Nmax, double dx,x0,x1);
*
* Imports:
*    void Bess(IN int m, IN double x, OUT double& Jx, OUT double& JDx);
*    void BessYproc(IN int m, IN double x, IN double eps,
*                   OUT double& Yx, OUT double& YDx);
*    struct zerofunc_params_struct zf_params;
*
*-------------------------------------------------------------------------

*! begin fragment
	double precision function zerofunc(x)

	implicit double precision(a-h,o-z)
	DIMENSION BI(0:250),DI(0:250),BK(0:250),DK(0:250)

	include 'fkappa.fh'
	include 'projstrc.fh'

	parameter (BIG=1.0d36,SMALL=1.0d-36)

	m=zf_params.m
	a=zf_params.a
	b=zf_params.b
	v=zf_params.v
	eps=zf_params.eps

d	print *,'zfi: m,x: ',m,x

	zerofunc=BIG

	if (dabs(x*a).gt.750.0d0) return	!!! abXX may be unfilled !!!

	if(Project.iUseCoredTube.eq.4)then

	xcr = dsqrt(v - x*x)

	call Bess(m,x*a,bja,bjap)
	call Bess(m,x*b,bjb,bjbp)
	call BessYproc(m,x*a,eps,bya,byap)
	call BessYproc(m,x*b,eps,byb,bybp)

      call IKNA(m,xcr*b,NM,BI,DI,BK,DK)

      zerofunc=x*(bjbp*bya-bja*bybp) -
	*xcr*DI(m)*(bjb*bya-bja*byb)/BI(m)

	else if(Project.iUseCoredTube.eq.3)then
c	if(v.ge.0.0d0)then

	if(v.gt.x*x)then

	xcr = dsqrt(v - x*x)

	call Bess(m,x*a,bja,bjap)
	call Bess(m,x*b,bjb,tmp)
	call BessYproc(m,x*a,eps,bya,byap)
	call BessYproc(m,x*b,eps,byb,tmp)

      call IKNA(m,xcr*a,NM,BI,DI,BK,DK)

      zerofunc=x*(bjap*byb-bjb*byap) -
	*xcr*DK(m)*(bja*byb-bjb*bya)/BK(m)

	else if(x*x.gt.v)then

	xcr = dsqrt(x*x - v)

	call Bess(m,x*a,bja,bjap)
	call Bess(m,x*b,bjb,tmp)
	call BessYproc(m,x*a,eps,bya,byap)
	call BessYproc(m,x*b,eps,byb,tmp)

	call BessYproc(m,xcr*a,eps,by2a,by2ap)

      zerofunc=x*(bjap*byb-bjb*byap) -
	*xcr*by2ap*(bja*byb-bjb*bya)/by2a

	else
	  stop '*ERR*  V = X*X'
	end if

	else

	if (b.lt.SMALL) then	! Original Cylindrical LACW, 0 == B <= r <= A
	  call Bess(m,x*a,bja,tmp)
	  zerofunc=bja
	else if (a.lt.SMALL) then	! Inverted LACW, B <= r <= +INF
	  call BessYproc(m,x*b,eps,byb,tmp)
	  zerofunc=byb
	else			! Tubular LACW, 0 < B <= r <= A
	  call Bess(m,x*a,bja,tmp)
	  call Bess(m,x*b,bjb,tmp)
	  call BessYproc(m,x*a,eps,bya,tmp)
	  call BessYproc(m,x*b,eps,byb,tmp)
	  zerofunc=bja*byb-bjb*bya
	endif

	end if

	return
	end

*-------------------------------------------------------------------------
* FindKappa(Dest,Mdim,Ndim,dx,x0,x1,Mmax,Nmax,Nact)
* scans the [x0,x1] interval for zeros of zerofunc(x) function.
* The step (dx) is divided by  dxcoeff  until the number of zeros
* in this interval changes no more.
*
* Mdim,Ndim     -- allocated dimensions of the Dest array;
* Mmax,Nmax     -- requested m,n maximum values;
* Nact[Mmax]    -- array to return numbers of roots actually found
*
	subroutine FindKappa(Dest,Mdim,Ndim,dx_,x0,x1,Mmax,Nmax,Nact)

	implicit double precision(a-h,o-z)
	dimension Dest(Mdim,Ndim)
	dimension Nact(Mmax)

	parameter (SMALL=0.0d0, BIG=1.0d36)
	parameter (dxcoeff=4.0d0)

	double precision zerofunc
	external zerofunc

	include 'fkappa.fh'

* Fill Dest with big-energy (big-kappa) values for the case
* if there's too big number of zeros is requested:
	do m=1,Mdim
	  do n=1,Ndim
	    Dest(m,n)=BIG
	  enddo
	  Nact(m)=0
	enddo

	eps=zf_params.eps
	dx=dx_	! <--MSFPS does not allow changing parameter variables (?!)

d	print *,'FindKappa(): Mmax,Nmax,dx,x0,x1:',Mmax,Nmax,dx,x0,x1
d	print *,'FK(): A,B,eps:',zf_params.A,zf_params.B,zf_params.eps

	x0_new=x0
* We'll use x0 or bigger value to stay out of the x=0 peculiarity of
* Y(x) functions, for the m-cycle is started for all m=1,Mdim, and
* for small starting point x0 we obtain an out-of-range value of Y
* on the tail of m's (m=33, 34, ...), where there's surely no roots.

	! write to log file:
	write(16,*)
	write(16,*) 'Kappas found:'

* Start search: for each  m, decrease  dx  until the number of zeros
* in the [x0,x1] interval does not change:
	dx_sav=dx
	do m=1,Mmax

	  zf_params.m=m-1

	  nr_old=-1	! was 0, when used "loop--until"s
	  nr_new=0
	  ndxr=-1
	  dx=dx_sav*dxcoeff

!!	  loop
	  do while (nr_old.ne.nr_new)

	  nr_old=nr_new
	  nr_new=0
	  dx=dx/dxcoeff
	  ndxr=ndxr+1

	  if (x0_new.ge.x1) exit
	  x=x0_new-dx
!	  x=x0-dx  ! x0 is the unmodified starting point of search
!!	loop
	  y=BIG
	  do while (dabs(y).ge.BIG.and.x.lt.x1)
	    x=x+dx
!mdv		if(x.gt.x1)x=x1
		if(x.gt.x1)x=x1
	    y=zerofunc(x)
	  enddo
!!	until (dabs(y).lt.BIG.or.x.ge.x1)
d	  print *,'x0_old,x0_new:',x0_new,x
	  x0_new=x

!!	  while (x.lt.x1) do
	  do while (x.lt.x1)
	    y0=y
	    x=x+dx
!mdv		if(x.gt.x1)x=x1
		if(x.gt.x1)x=x1
	    y=zerofunc(x)
	    if (dabs(y).ge.BIG) then
	      print *,'*ERR* ZeroFunction fail; m,x:',m-1,x
	      return	!!! WARNING: 'return' or 'stop'?
	    endif
d	    print *,'xs,xe,ys,ye:',x-dx,x,y0,y
	    if ((y0.lt.SMALL.and.y.gt.SMALL).or.
     .		(y0.gt.SMALL.and.y.lt.SMALL)) then
	! function sign has changed; find the routine more precisely:
	      z=FindZero(zerofunc,x-dx,x,eps)
	      if (z.ge.BIG) continue !stop 'ERROR: Function routine search error'
d	      print *,'Found zero at x =',z
	      if (nr_new.gt.Nmax) exit
	      nr_new=nr_new+1
	      Dest(m,nr_new)=z

d	      CALL BESS_NUL(N,IABS(M-1),temp,50)
d	      print *,'cmp n,o ->',z,temp(nr_new)/zf_params.a
c	      if (nr_new.ge.Nmax.and.x.le.x1) x1=x-dx    !!!
	    endif
	  enddo
!!	  endwhile

d	  print *,m-1,dx,', root numbers (o->n):',nr_old,nr_new
	  enddo
!!	  until (nr_old.eq.nr_new)

	  ! write to log file:
	  if (nr_new.ne.0) then
	  write(16,*) 'm =',m-1
	  do i=1,nr_new
	    write(16,*) i,Dest(m,i)
d	    write(*,*) i,Dest(m,i)
	  enddo
	  endif

	  Nact(m)=nr_new
d	  print *,'dx reduction steps needed:',ndxr

	enddo

	return
	end

*-------------------------------------------------------------------------
* FindBessAtAB(IN a,b,pdKappas,Mmax,Nmax,NAct, OUT abJ,abJD,abY,abYD)
* Calculates J(Kappa*r),Y,J',Y' for all KappaMN, r=a,b.
* NOTE: pdKappas(1,*) corresponds to "physical" m=0.
*
      subroutine FindBessAtAB(a,b,YEps,
     . pdKappas,Mmax,Nmax,NAct, abJ,abJD,abY,abYD)

      implicit double precision(a-h,o-z)
      dimension pdKappas(Mmax,Nmax),NAct(Mmax)
      dimension abJ(Mmax,Nmax,2),abJD(Mmax,Nmax,2),
     . abY(Mmax,Nmax,2),abYD(Mmax,Nmax,2)

      parameter (BIG=1.0d12)

*! begin
c      parameter (NBESS=2048)
c      dimension bj(NBESS),by(NBESS)
*! end

      do m=1,Mmax
      do n=1,NAct(m)	!Nmax

      x=pdKappas(m,n)
      if (x.ge.BIG) exit

*! begin
c      nm=max0(2*m,2*idnint(x*a),20)
c      nm=nm*2     ! Y() calc requires 2 times larger set
c      if (nm.gt.NBESS) stop '*ERR* Accuracy loss [FindBessAB]'
c      if (dabs(x*a).gt.750) exit       !!! abXX may be unfilled !!!

c      call CfBesselJY(nm,x*a,bj,by)

c      abJ(m,n,1)=bj(m+1)
c      abY(m,n,1)=by(m+1)
c      abJD(m,n,1)=(bj(m)-bj(m+2))*0.5d0
c      abYD(m,n,1)=(by(m)-by(m+2))*0.5d0

c      call CfBesselJY(nm,x*b,bj,by)

c      abJ(m,n,2)=bj(m+1)
c      abY(m,n,2)=by(m+1)
c      abJD(m,n,2)=(bj(m)-bj(m+2))*0.5d0
c      abYD(m,n,2)=(by(m)-by(m+2))*0.5d0
*! else
d      print *,'FBAAB: m,x:',m-1,x
      call Bess(m-1,x*a,abJ(m,n,1),abJD(m,n,1))
      call BessYproc(m-1,x*a,YEps,abY(m,n,1),abYD(m,n,1))
      call Bess(m-1,x*b,abJ(m,n,2),abJD(m,n,2))
      call BessYproc(m-1,x*b,YEps,abY(m,n,2),abYD(m,n,2))
d      write(16,'(a3,2i4,5f10.4,2x,5f10.4,4x,2f10.4)') 'fb',m-1,n,
d     . x*a,abJ(m,n,1),abY(m,n,1),abJD(m,n,1),abYD(m,n,1),
d     . x*b,abJ(m,n,2),abY(m,n,2),abJD(m,n,2),abYD(m,n,2)
d     . ,abJ(m,n,1)/abY(m,n,1),abJ(m,n,2)/abY(m,n,2)
d      if (m-1.eq.0) then
d         call BessYproc(m-1,x*b,YEps,tmp,tmp1)
d         write(16,*) 'sp: m,n,x,xb,Y,Y'':',m-1,n,x,x*b,tmp,tmp1
d      endif
*! end

      enddo
      enddo

      return
      end

*--------------------------------------------------------------------
* double ChooseBessAtAB(m,n,Mmax,Nmax,abXX,i)
* Returns Bessel function (stored at abXX, filled by FindBessAtAB())
* value at r=a (i=1) or r=b (i=2).
* m may be negative (Bessel will be *(-1)^|m|;
* m -- "physical" (from 0).
*
      double precision function ChooseBessAtAB(m,n,Mmax,Nmax,abXX,i)

      implicit double precision(a-h,o-z)
      dimension abXX(Mmax,Nmax,2)

      ret=abXX(iabs(m)+1,n,i)
      if (m.lt.0.and.mod(m,2).ne.0) ret=-ret

d      print *,'CBAAB:',m,n,ret
      ChooseBessAtAB=ret
      return
      end

*! end fragment
