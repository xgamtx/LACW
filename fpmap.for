	subroutine FPMakeMaps(file_name,nt,
     +	x0,y0,z0,dx,dy,dz,nx,ny,nz,iCoordSys,
     +	rad,amtc,wpot,step,rmts,qeff,NGRID,
     +	nts,rsx,rsy,rsz,NQG,Ncell,C,Co,logf)

	implicit double precision(a-h,o-z)

	parameter (EVS=13.6058D0)

	include 'coordsys.fh'
	include 'varsizes.fh'

	! return info class:
	parameter ( icNone=0, icDensity=1, icVspcharg=2, icVnuclei=3,
     .	icVxc=4, icVfp=5, !icVmnp=6, icVmnpIS=7,
     .	icIntMap=8, icSpatialGrid=256 )

	include 'xcpot.fh'
	record /cxcpot/ Co

	integer start_index

	double precision rad(NGRID,*),amtc(NGRID,*),wpot(NGRID,*),
     +	step(*),rmts(*),qeff(*),rsx(NQG,*),rsy(NQG,*),rsz(NQG,*)
	integer nts(*)	!,jris(*)

	double precision, allocatable:: xset(:),yset(:),zset(:),
     +	RhoMap(:),V2Map(:),V1Map(:),VxcMap(:),VfpMap(:)!,dists(:,:)
	character, allocatable:: int_map(:)

	character file_name*(*)

	integer logf

	common /debug/ DEBUG,ierr
	data dmin/1.0d-10/

	!allocate memory for temp array:
c	allocate(dists(NQG,nt),stat=istat)
c	if (istat.eq.0.and.DEBUG.ge.7.and.ierr.ne.0) write(ierr,*)
c     +	'FP: allocated the "dists" array,',NQG*nt*8,' bytes'
c	if (istat.ne.0) stop '*ERR* Low memory [FPdists]'

	pi=4.0d0*datan(1.0d0)
c	max_index=(nx-1)*(ny-1)*(nz-1)
	max_index=nx*ny*nz
	domega=dx*dy*dz	! elementary volume factor

d	write(logf,*) x0,y0,z0,dx,dy,dz,nx,ny,nz,Ncell,C
d	do it=1,nt
d	write(logf,*) it,'--->',jris(it),rmts(it),qeff(it),nts(it)
d	do i=1,NGRID
d	write(logf,*) rad(i,it),amtc(i,it),wpot(i,it)
d	enddo
d	enddo

d	do it=1,nt
d	write(logf,*) 'Atom type',it
d	do i=2,NGRID
d	write(logf,*) rad(i,it),
d     +	ATINT_(amtc,rad,i,it,step(it),NGRID,it,0,1.0d0)
d	enddo
d	enddo
d	stop

	! try to allocate memory for spatial maps:
	if (logf.ne.0) write(logf,*)'FP: Memory alloc request [FPSPM],',
     +	max_index*(IDBL*8+ICHAR),' bytes...'
	write(*,*) 'FP: Spatial maps of size',
     +	max_index*(IDBL*8+ICHAR),' bytes...'

	allocate(RhoMap(max_index),V2Map(max_index),V1Map(max_index),
     +	VxcMap(max_index),VfpMap(max_index),
     +	xset(max_index),yset(max_index),zset(max_index),
     +	int_map(max_index),stat=istat)

	if (istat.eq.0) write(*,*) '...allocated successfully'
	if (istat.ne.0) stop '*ERR* Low memory [FPSM]'

	! initialize structures:
c	if (start_index.le.1) then	! only if the calculation is beginning
	do index=1,max_index
	   RhoMap(index)=0.0d0
	   V2Map(index)=0.0d0
	   V1Map(index)=0.0d0
	   VxcMap(index)=0.0d0
	   VfpMap(index)=0.0d0
	   int_map(index)=char(0)
	enddo
c	endif

*--------------------------------------------------------------------
	start_index=0

* Try to read fpmaps.dat. If incomplete, continue the calculation.
	open(26,file='fpmaps.dat',form='unformatted',err=142,status='old')
	do while(.true.)
	   read(26,err=145,end=143) i	! icXXXXX, info class
	   if (i.eq.icNone) then
	      read(26) nt_chk,ix,iy,iz,index
	      if (nt_chk.ne.nt.or.ix.ne.nx.or.iy.ne.ny.or.iz.ne.nz) stop
     +	'FP: Error: Invalid or corrupt fpmaps.dat. Restart Fpxcy2.'
           elseif (i.eq.icSpatialGrid) then
		read(26,err=145,end=145)
     +			(xset(i),yset(i),zset(i),i=1,max_index)
           elseif (i.eq.icVfp) then
		read(26,err=145,end=145) (VfpMap(i),i=1,index)
           elseif (i.eq.icVspcharg) then
		read(26,err=145,end=145) (V1Map(i),i=1,max_index)	! the charge in each dV affects the V1 potential everywhere
           elseif (i.eq.icVnuclei) then
		read(26,err=145,end=145) (V2Map(i),i=1,index)
           elseif (i.eq.icVxc) then
		read(26,err=145,end=145) (VxcMap(i),i=1,index)
           elseif (i.eq.icDensity) then
		read(26,err=145,end=145) (RhoMap(i),i=1,index)
           elseif (i.eq.icIntMap) then
		read(26,err=145,end=145) (int_map(i),i=1,index)
	   endif
	enddo

  143	start_index=index

  145	close(26)
  142	continue

*--------------------------------------------------------------------
d	write(logf,*) x0,y0,z0,dx,dy,dz,nx,ny,nz
	x=(x0+dx*dble(nx))
	Omega=pi*(x*x-x0*x0)*(dz*dble(nz))
	if (logf.ne.0) write(logf,*)
     +		'Omega by geometry (==pi*(a*a-b*b)*c):',Omega

*--------------------------------------------------------------------
* Prepare spatial map:
	s=0.0d0
	index=1
	   z=z0+dz/2.0d0
	do iz=1,nz
	   y=y0+dy/2.0d0
d	   write(*,'(a1\)') '.'
	do iy=1,ny
	   x=x0+dx/2.0d0
	do ix=1,nx
	   xset(index)=x
	   yset(index)=y
	   zset(index)=z
	   s=s+domega*dVolumeCoeff(iCoordSys,x,y,z)	! [coordsys.for]
c	   s=s+x*dx*dy*dz
c	   call ToCartesian3D(x,y,z,
c     +		xset(index),yset(index),zset(index),iCoordSys)
	   index=index+1	    ! offset++
	   x=x+dx
	enddo     !ix
	   y=y+dy
	enddo     !iy
	   z=z+dz
	enddo     !iz
	if (logf.ne.0) write(logf,*)
     +		'Omega by grid (==integral(1*dOmega)):',s
	if (logf.ne.0) write(logf,'(''(difference is'',i3,''%)'')')
     +		idnint(100.0d0*dabs(s-Omega)/Omega)

c	! save the grid to file:
c	open(26,file='mapgrid.dat',form='unformatted')
c	write(26) max_index
c	write(26) (xset(i),yset(i),zset(i),i=1,max_index)
c	close(26)

	! setup XCPOT constants:
	call sXCPOT(Co)

*--------------------------------------------------------------------

* Using AtomFunc.for::CountElectrons(rad,amtc,NGRID,jri,nt,step)
* to get the most exact value of electronic charge of the cell's atoms:
	sum_el=CountElectrons(rad,amtc,NGRID,NGRID,nt,nts,step)
	if (logf.ne.0) write(logf,*) 'Total nuclear charge:  ',sum_el
	write(*,*) 'Total nuclear charge:  ',sum_el

* Now let's calculate total atomic electron density that is inside
* our potential tube PER ONE CELL. We integrate it up to Ncell cells
* in each direction of Z. We would want to obtain the number of
* electrons, with its abs value equal to the sum of nuclear charges
* of the atoms. If it's not equal to that, it's due to the roughness
* of the spatial grid, and we cheat by re-normalising the atomic density
* or effective charges.

	sum_rho=0.0d0
	do icell=-Ncell,Ncell
	rho=0.0d0

	do it=1,nt
	do is=1,nts(it)
c Modify the above cycles to sum the density of a single atom
c (we would then obtain its full electron charge, or something).
d	write(logf,*) 'aaa',it,is,rsz(is,it)+dble(icell)*C

	do index=1,max_index
	x=xset(index)
	y=yset(index)
	z=zset(index)
	call ToCartesian3D(x,y,z,xx,yy,zz,iCoordSys)

	   dist=dist3d(xx,yy,zz,
     +			rsx(is,it),rsy(is,it),rsz(is,it)+dble(icell)*C)
	   if (dist.lt.dmin) cycle
	   drho=FindYatom(rad,amtc,NGRID,dist,it,idummy)/(dist*dist)
	   rho=rho+drho*dVolumeCoeff(iCoordSys,x,y,z)*domega
	enddo	! index

	enddo
	enddo

	rho=rho/(4.0d0*pi)
d	write(logf,*) icell,rho
	sum_rho=sum_rho+rho
	enddo	! icell
	if (logf.ne.0) write(logf,*) 'Total electron density:',sum_rho
	write(*,*) 'Total electron density:',sum_rho

* Shift the effective charges by the value of the difference:
	tmp=sum_rho/sum_el
c	tmp=sum_rho-sum_el
	do it=1,nt
c	  qeff(it)=qeff(it)-tmp
	  qeff(it)=qeff(it)*tmp
	enddo
	if (logf.ne.0) write(logf,*) 'Effective charges adjusted; *=',tmp
	write(*,*) 'Effective charges adjusted; *=',tmp
d	stop '*STOP* After Qeff adjustment'


*--------------------------------------------------------------------

	print *,'FP: Preparing spatial map...'
	iSavePercent=100*start_index/max_index

	sum_rho=0.0d0
	do index=start_index+1,max_index
	x=xset(index)
	y=yset(index)
	z=zset(index)
	call ToCartesian3D(x,y,z,xx,yy,zz,iCoordSys)

	rho=0.0d0
	v2=0.0d0
	is_internal=0

	do icell=-Ncell,Ncell	! sum over this and neighbouring elementary cells
c	do icell=-Ncell-1,Ncell+1	! sum over this and neighbouring elementary cells

d	rho=0.0d0
	do it=1,nt	! sum over atom types
	do is=1,nts(it)	! sum over atoms of a selected type

c	if (dabs(-zz+rsz(is,it)+dble(icell)*C)
c     +		.gt.C*dble(Ncell*2+1)*0.5d0) cycle

	   dist=dist3d(xx,yy,zz,
     +			rsx(is,it),rsy(is,it),rsz(is,it)+dble(icell)*C)
	   if (dist.lt.rmts(it)) is_internal=is_internal+1
	   if (dist.lt.dmin) cycle
	   drho=FindYatom(rad,amtc,NGRID,dist,it,idummy)/(dist*dist)
	   rho=rho+drho

	   if (dist.lt.rmts(it)) VfpMap(index)=wpot(idummy,it)

d	   write(logf,*) is,icell,dist
	   if (dist.lt.rmts(it)) dist=rmts(it)
	   v2=v2-2.0d0*qeff(it)/dist

d	   if (iz.eq.1) write(logf,'(i2,f8.4,f10.6,'' |'',6f8.4)')
d     +	iq,dist,drho,xx,yy,zz,rsx(is,it),rsy(is,it),rsz(is,it)+dble(iq)*C

	enddo
	enddo
d	write(logf,*) icell,rho*x*domega/(4.0d0*PI),sum_rho
	enddo
	int_map(index)=char(is_internal)
	RhoMap(index)=rho/(4.0d0*PI)
	V2Map(index)=v2

	! exchange-correlation potential:
	rho_serv1=RhoMap(index)*2.0d0*PI
	rho_serv2=rho_serv1
	rho_serv=rho_serv1+rho_serv2
	call XCPOT(rho_serv1,rho_serv2,rho_serv,Vxc1,Vxc2,dummy,Co)
	VxcMap(index)=Vxc1

c	if (iCoordSys.eq.iCoords3dCartesian) then
c	   domega1=domega
c	elseif (iCoordSys.eq.iCoords3dSpherical) then
c	   domega1=x*x*dsin(y)*domega
c	elseif (iCoordSys.eq.iCoords3dCylindrical) then
c	   domega1=x*domega
c	else
c	   stop 'Strcy: Internal error [FP/CoordSys]'
c	endif
	domega1=dVolumeCoeff(iCoordSys,x,y,z)*domega
	t=RhoMap(index)*domega1	! charge in the elementary volume
	sum_rho=sum_rho+t	!+domega1 to check---should be ==Omega. (It is.)
c	if (t.ge.dmin) then
	if (t.ge.dmin.and.is_internal.eq.0) then

	! update the V1 potential map:
	do indexo=1,max_index
	xo=xset(indexo)
	yo=yset(indexo)
	zo=zset(indexo)
	call ToCartesian3D(xo,yo,zo,xxo,yyo,zzo,iCoordSys)

	do icello=-Ncell,Ncell
c	if (icello.eq.0) cycle	! skip the rho element in the current cell (we've evidently allowed for its contribution to V1)
	dist=dist3d(xx,yy,zz,xxo,yyo,zzo+dble(icello)*C)
	if (dist.ge.dmin) !.and.ichar(int_map(indexo)).eq.0)
     .	V1Map(indexo)=V1Map(indexo)+2.0d0*t/dist
	enddo	! icello

	enddo	! indexo

	endif	! need to sum (the point is outside MTS)


	i=100*index/max_index
	write(*,'(a1,i3,''%''\)') char(13),i
	if (index.eq.max_index.or.
     +		(i.gt.iSavePercent.and.mod(i,10).eq.0)) then
*=========================================================================
*! begin fragment SP 5.09: moved from below
	! save results each 10% of the calculation:
d	write(*,*) ' Emergency save at ',i
	iSavePercent=i

	if (index.eq.max_index) then	! finalise the calculation
*-------------------------------------------------------------------------
	do i=1,max_index
c Comment the condition line to skip update of VfpMap() for internal cells.
c Vfp will be left equal to wpot(), the pure atomic spherical potential;
c however, Vfp will not be continuous on MT-borders.
c	   if (int_map(i).ne.char(0)) cycle
	   VfpMap(i)=V1Map(i)+V2Map(i)+VxcMap(i)
	enddo

	! calculate average, minimum, and maximum potential outside MT spheres:
	VfpAvgOutside=0.0d0
	VfpMinOutside=+1.0d+6
	VfpMaxOutside=-1.0d+6
	j=0
	do i=1,max_index
d	  if (int_map(i).eq.char(0)) write(25,*) i,'ext'
d	  if (int_map(i).ne.char(0)) write(25,*) i,'int'
	  if (int_map(i).ne.char(0)) cycle
	  j=j+1
	  t=VfpMap(i)
	  VfpAvgOutside=VfpAvgOutside+t
d	  write(24,*) i,t,VfpAvgOutside
	  if (t.gt.VfpMaxOutside) VfpMaxOutside=t
	  if (t.lt.VfpMinOutside) VfpMinOutside=t
d	  x=xset(i)
d	  y=yset(i)
d	  z=zset(i)
d	  call ToCartesian3D(x,y,z,xx,yy,zz,iCoordSys)
d	  write(logf,'(a4,4f10.4)') 'aver',xx,yy,zz,VfpMap(i)
	enddo
	if (logf.ne.0) then
d	write(logf,*) 'sum_rho:',sum_rho
	write(logf,*) idnint(100.d0*dble(j)/dble(max_index)),
     +		'% of points outside MTS'
	VfpAvgOutside=VfpAvgOutside/dble(j)
	write(logf,*) 'Potential in the interspherical region:'
        write(logf,*) '	maximum:',VfpMaxOutside*evs,' eV'
        write(logf,*) '	average:',VfpAvgOutside*evs,' eV'
	write(logf,*) '	minimum:',VfpMinOutside*evs,' eV'
	endif	! logf.ne.0
d	write(*,*) 'sum_rho:',sum_rho
	write(*,*) 'Potential in the interspherical region:'
        write(*,*) '	maximum:',VfpMaxOutside*evs,' eV'
        write(*,*) '	average:',VfpAvgOutside*evs,' eV'
	write(*,*) '	minimum:',VfpMinOutside*evs,' eV'
*-------------------------------------------------------------------------
	endif

	! save maps to binary file:
	open(26,file=file_name,form='unformatted')
	write(26) icNone
	write(26) nt,nx,ny,nz,index,x0,y0,z0,dx,dy,dz,iCoordSys,
     +	VfpAvgOutside,VfpMaxOutside,VfpMinOutside
	write(26) icSpatialGrid
	write(26) (xset(i),yset(i),zset(i),i=1,max_index)
	write(26) icVfp
	write(26) (VfpMap(i),i=1,index)
	write(26) icVspcharg
	write(26) (V1Map(i),i=1,max_index)	! the charge in each dV affects the V1 potential everywhere
	write(26) icVnuclei
	write(26) (V2Map(i),i=1,index)
	write(26) icVxc
	write(26) (VxcMap(i),i=1,index)
	write(26) icDensity
	write(26) (RhoMap(i),i=1,index)
	write(26) icIntMap
	write(26) (int_map(i),i=1,index)
	close(26)

	endif	! if the last point, or each 10%
*! end fragment
*=========================================================================

	enddo	! main cycle over spatial map points

c	deallocate(dists)


*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	open(26,file='fpmaps.txt',form='formatted')
	open(27,file='v1map.txt',form='formatted')
	open(28,file='v2map.txt',form='formatted')
	open(29,file='vxcmap.txt',form='formatted')
	open(25,file='vfpmap.txt',form='formatted')
	open(30,file='rhomap.txt',form='formatted')
	! select one constant coordinate to plot potential surfaces:
	z_sel=1.0d+6
	phi_sel=1.0d+6
	do index=1,max_index
	   x=xset(index)
	   y=yset(index)
	   z=zset(index)
	   call ToCartesian3D(x,y,z,xx,yy,zz,iCoordSys)
	   if (zz.lt.z_sel) z_sel=zz
	   if (y.lt.phi_sel) phi_sel=y
	enddo
	! plot potential and density surfaces
	do index=1,max_index
	   x=xset(index)
	   y=yset(index)
	   z=zset(index)
	   call ToCartesian3D(x,y,z,xx,yy,zz,iCoordSys)
	   if (zz.eq.z_sel) then	!.and.int_map(index).ne.char(1))	!(nz/2+1))
		write(26,'(6f10.5)') xx,yy,
     +	RhoMap(index),V1Map(index),V2Map(index),VxcMap(index)
		write(27,'(3f10.5)') xx,yy,V1Map(index)
		write(28,'(3f10.5)') xx,yy,V2Map(index)
		write(29,'(3f10.5)') xx,yy,VxcMap(index)
		write(25,'(3f10.5)') xx,yy,VfpMap(index)
		write(30,'(3f10.5)') xx,yy,RhoMap(index)
	   endif
c	   if (dabs(phi_sel+pi-y).lt.dy) x=-x
c	   if (y.eq.phi_sel.or.dabs(phi_sel+pi-y).lt.dy) then	!.and.int_map(index).ne.char(1))	!(nz/2+1))
c		write(26,'(6f10.5)') x,z,
c     +	RhoMap(index),V1Map(index),V2Map(index),VxcMap(index)
c		write(27,'(3f10.5)') x,z,V1Map(index)
c		write(28,'(3f10.5)') x,z,V2Map(index)
c		write(29,'(3f10.5)') x,z,VxcMap(index)
c		write(25,'(3f10.5)') x,z,VfpMap(index)
c		write(30,'(3f10.5)') x,z,RhoMap(index)
c	   endif
	enddo
	close(26)
	close(27)
	close(28)
	close(29)
	close(25)
	close(30)

*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	deallocate(RhoMap,V2Map,V1Map,VxcMap,VfpMap,
     +	xset,yset,zset,int_map,stat=istat)
	if (istat.eq.0.and.DEBUG.ge.7.and.ierr.ne.0) write(ierr,*)
     +	'FP: memory de-allocation successful'

	return
	end
