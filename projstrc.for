	integer function iReadProject(fn,logf)
	implicit double precision(a-h,o-z)

	include 'projstrc.fh'
	include 'coordsys.fh'
	include 'atomdata.fh'

	character*(*) fn
	integer logf

	external GetLogical
	logical GetLogical

	character buf*256,str*256
	double precision rs(3),eetmp(LCUT)

* test the input file:
	open(25,file=fn,form='formatted',err=215)
	close(25)

*****
	B_SMALL=1.0d-5
	PHI_SMALL=1.0d-8
*****
	A=-1.0d0
	B=-1.0d0
	C=-1.0d0
	V=-1.0d0
	VI=-1.0d0
	ChiralityAngle=0.0d0
	Project.bAddPCompanion=.FALSE.
	Project.EnergyCutOff=0.0d0	! no basis

	ifs=iIniFileOpenSection(fn,'Structure',25,0)	!..,1 !!!
	dowhile (iIniFileReadEntryString(ifs,buf).gt.0)
	m=index(buf,'=')
	if (m.le.1) continue
	str=buf(1:m-1)
	buf=buf(m+1:)
        call BTrimEx(buf,256)
	if (str.eq.'A') read(buf,*) A
	if (str.eq.'B') read(buf,*) B
	if (str.eq.'C') read(buf,*) C
	if (str.eq.'F') read(buf,*) ChiralityAngle
	if (str.eq.'V') then
        read(buf,*) V
	  V=V/evs	! to electronvolts
	end if
	if (str.eq.'VI') then
        read(buf,*) VI
	  VI=VI/evs	! to electronvolts
	end if
	if (str.eq.'ExternalRadius') read(buf,*) A
	if (str.eq.'InternalRadius') read(buf,*) B
	if (str.eq.'CellLength') read(buf,*) C
	if (str.eq.'ChiralityAngle') read(buf,*) ChiralityAngle
	if (str.eq.'ChiralityStepping') then
	  read(buf,*) i
	  if (i.gt.0) ChiralityAngle=2.0d0*pi/i
	endif
	if (str.eq.'AddPCompanion')
     +		Project.bAddPCompanion=GetLogical(buf,.FALSE.)
	if (str.eq.'EnergyCutOff') then
	  read(buf,*) Project.EnergyCutOff
	  Project.EnergyCutOff=Project.EnergyCutOff/evs	! to electronvolts
	endif
	enddo

*****

	if (dabs(ChiralityAngle).lt.PHI_SMALL) ChiralityAngle=0.0d0
	if (B.lt.B_SMALL) B=0.0d0
	if (A.lt.B_SMALL) A=0.0d0
	if (A.eq.B) stop
     +	'*ERR* Internal and external radii (A and B) can''t be equal'
*
      if (V.ge.0.0d0.and.VI.ge.0.0d0) then
c The case present V and VI[nternal] parameter is mean nanotube in crystal
c method may be well realized in future
	  Project.iUseCoredTube=5
      else if (VI.ge.0.0d0) then
c The case present VI[nternal] parameter is mean nanotube inside the crystal
c method based on iUseCoredTube = 1
	  Project.iUseCoredTube=4
      else if (V.ge.0.0d0) then
c The case present V parameter is mean nanotube in crystal
c method based on iUseCoredTube = 1
	  Project.iUseCoredTube=3
	else if (A.lt.B) then
c	  write(*,'(''*'',/,''*'',a78,/,''*'')')
c     +		'* Inverted LACW: 3D medium with empty cylindrical channel'
	  Project.iUseCoredTube=2
	  A=B	!=0.0d0	!!! DUMMY! Physically, A==+INFINITY
c The case A<B, A!=0 is not interesting: physically, it would mean that
c electrons are allowed to travel inside the internal cylinder (of radius A)
c and outside the external one (of radius B). Potential barriers are considered
c to be impenetrable, so those two potential regions are quite independent
c and can be investigated separately.
	else if (B.gt.B_SMALL) then
c The case 0 < B < A is the tubular LACW case
c	  write(*,'(''*'',/,''*'',a78,/,''*'')')
c     +		'Tubular LACW: Cylinder with internal hole'
	  Project.iUseCoredTube=1
	else
c The case 0 = B < A is the original LACW case
c	  write(*,'(''*'',/,''*'',a78,/,''*'')')
c     +		'Original LACW: Cylinder without internal hole'
	  Project.iUseCoredTube=0
	endif
*
	if (Project.iUseCoredTube.lt.0.or.Project.iUseCoredTube.gt.4)
     +	stop 'Strcy: Internal error: Unknown program variant requested'
*
d	if (logf.ne.0) write(logf,*)
d     +	'***** METHOD CLONE: ',TCalcVar(Project.iUseCoredTube+1)
*
*****
	it=0
	IQChrg=0
	bQCF(:,:)=.FALSE.
	do
	  write(unit=buf,fmt='(i5)') it+1
	  call BTrimEx(buf,256)
	  buf='Structure/AtomTypes/'//buf
	  ifs=iIniFileOpenSection(fn,buf,25,0)
	  if (ifs.eq.0) exit
	  it=it+1
	  nts(it)=0
	  nel(it)=0
	  zmain(it)=0

	  bUseVpot(it)=.FALSE.
	  Vpot(it)=0.0d0

	dowhile (iIniFileReadEntryString(ifs,buf).gt.0)

	if (index(buf,'`').eq.1.and.nts(it).lt.NQG) then
	  buf=buf(2:)
	  if (iStringToDoubleArray(buf,rs,3,ierr).ne.0) then
	  i=nts(it)+1
	  rsx(i,it)=rs(1)
	  rsy(i,it)=rs(2)
	  rsz(i,it)=rs(3)
	  nts(it)=i
	  endif
	endif

	m=index(buf,'=')
	if (m.le.1) cycle
	str=buf(1:m-1)
	buf=buf(m+1:)
        call BTrimEx(buf,256)
	if (str.eq.'@') read(buf,'(a12)') tAName(it)
	if (str.eq.'Z'.or.str.eq.'NuclearCharge') read(buf,*) zmain(it)
	if (str.eq.'Q'.or.str.eq.'ValentElectrons') read(buf,*) nel(it)
	if (str.eq.'VirtualPotential') then
	  bUseVpot(it)=.TRUE.	! WARNING: These two params should be independent on the level of input file
	  read(buf,*) Vpot(it)
	endif
	if (str.eq.'OrbitalFlags') then
	  do i=1,5	! WARNING: it should be "1,LCUT"
	    if (index(buf,QchText(i)//'[*]').gt.0) then
	      IQChrg=IQChrg+1
	      IQChrgFlg(1,IQChrg)=IT	! Save atom type
	      IQChrgFlg(2,IQChrg)=I-1	! and it's shell number.
	      bQCF(it,i)=.TRUE.
	    endif
	  enddo
	endif
	if (str.eq.'OrbitalEnergies') then
	  m=iStringToDoubleArray(buf,eetmp,LCUT,ierr)
	  do i=1,m
	    EE(i,it)=eetmp(i)/evs
	  enddo
	endif
	enddo	! dowhile -- reading the IniFile

	if (it.ge.NCUT) exit
	enddo
	if (it.ne.0) nt=it	!!!!! re-write for new format

	if (logf.ne.0) then
	write(logf,*) nt,' atom types:'
	do it=1,nt
	  write(logf,'(3a,i4,a,i3,a,i3)') '@=',tAName(it),
     +		': ',nts(it),' atoms; NuclearCharge = ',zmain(it),
     +		'; ValentElectrons = ',nel(it)
d	  do i=1,nts(it)
d	    write(logf,'(a1,3f20.10)') '`',rsx(i,it),rsy(i,it),rsz(i,it)
d	  enddo
	enddo
	endif
*****

	!
	!	Exchange-correlation parameters
	!
	Co.IXCH=-1
	ifs=iIniFileOpenSection(fn,'Exchange',25,0)
	dowhile (iIniFileReadEntryString(ifs,buf).gt.0)
	m=index(buf,'=')
	if (m.le.1) continue
	str=buf(1:m-1)
	buf=buf(m+1:)
        call BTrimEx(buf,256)

	if (str.eq.'Mode') then
	  do i=1,4
	    if(index(buf,TMTHDS(i)).gt.0) then
	      Co.IXCH=I-1
	      exit
	    endif
	  enddo
	endif
	if (str.eq.'SlaterAlpha') read(buf,*) Co.ALPHA
	enddo	! instrcy/Exchange/

	if (Co.IXCH.lt.0) stop 'FPX: Unknown exchange mode'
d	if (logf.ne.0) write(logf,*) 'FPX: Exchange mode: ',
d     +	TMTHDS(Co.IXCH+1),'; SlaterAlpha:',Co.ALPHA

	CALL SXCPOT(Co)

*-------------------------------------------------------------------------
	if (iReadAtomic('outatm.str').eq.0)
     +	stop 'RDF: OutAtm.str is absent or corrupt. Re-run Atom'
*-------------------------------------------------------------------------

c	Project.iUseCoredTube=0
c	Project.bAddPCompanion=.FALSE.
c	Project.EnergyCutOff=?????
	Project.sJobTitle='untitled'
	Project.iBands=0	! must be >=(int)(Nelectrons/2)
	Project.iPointsTotal=1	! number of k points
	Project.iPointsNow=1
	Project.WvStart=0.0d0
	Project.EpsInt=0.0005d0
	Project.EpsCoords=0.0001d0
	Project.EpsBess=1.0d-12
	Project.NL=5
	Project.NLmax=LCUT-1
	Project.NLsav=Project.NL	! NL may vary up to NLmax depending on the convergency of integrals
	Project.NelCell=0
	Project.SaveI3M=.TRUE.
	Project.bIntersphericalDos=.FALSE.

	ifs=iIniFileOpenSection(fn,'Project',25,0)
	dowhile (iIniFileReadEntryString(ifs,buf).gt.0)
	m=index(buf,'=')
	if (m.le.1) continue
	str=buf(1:m-1)
	buf=buf(m+1:)
        call BTrimEx(buf,256)
	if (str.eq.'JobTitle') read(buf,*) Project.sJobTitle
	if (str.eq.'NumberOfBands') read(buf,*) Project.iBands

	if (str.eq.'IntegralAccuracy') read(buf,*) Project.EpsInt
	if (str.eq.'BesselAccuracy') read(buf,*) Project.EpsBess
	if (str.eq.'InputCoordsAccuracy')
     +		read(buf,*) Project.EpsCoords

	if (str.eq.'LValueInitial') read(buf,*) Project.NL
	if (str.eq.'LValueMaximum') read(buf,*) Project.NLmax

	if (str.eq.'CacheDoubleIntegrals')
     +		Project.SaveI3M=GetLogical(buf,.TRUE.)	!!! "... =<space>0" fails!
	if (str.eq.'MakeIntersphericalDos')
     +		Project.bIntersphericalDos=GetLogical(buf,.FALSE.)	!!! "... =<space>0" fails!

	if (str.eq.'AddPCompanion')
     +		Project.bAddPCompanion=GetLogical(buf,.FALSE.)
	if (str.eq.'EnergyCutOff') then
	  read(buf,*) Project.EnergyCutOff
	  Project.EnergyCutOff=Project.EnergyCutOff/evs	! to electronvolts
	endif

	if (str.eq.'kStart') then	! allow settings like "kStart=-pi/C"
	  if (index('-p',buf).gt.0.or.index('-P',buf).gt.0) then
	    Project.WvStart=-PI/C
	  else
	    read(buf,*) Project.WvStart
	  endif
	endif
	if (str.eq.'CalcIntervalsTotal') then
		read(buf,*) i
                Project.iPointsTotal=i+1
	endif
	if (str.eq.'CalcPointsTotal')
     +		read(buf,*) Project.iPointsTotal
	if (str.eq.'CalcPointsNow')
     +		read(buf,*) Project.iPointsNow

	enddo

*-------------------------------------------------------------------------
	Plots.EnergyBottom=0.0d0
	Plots.EnergyTop=0.0d0
	Plots.Halfwidth=0.25d0
	Plots.HalfwidthFermi=0.025d0
	Plots.iRenderPitch=15
	Plots.SeparateBandPlotFiles=.TRUE.

	ifs=iIniFileOpenSection(fn,'Project/Graphs',25,0)
	dowhile (iIniFileReadEntryString(ifs,buf).gt.0)
	m=index(buf,'=')
	if (m.le.1) continue
	str=buf(1:m-1)
	buf=buf(m+1:)
        call BTrimEx(buf,256)
	if (str.eq.'EnergyBottom') read(buf,*) Plots.EnergyBottom
	if (str.eq.'EnergyTop') read(buf,*) Plots.EnergyTop
	if (str.eq.'Halfwidth') read(buf,*) Plots.Halfwidth
	if (str.eq.'HalfwidthFermi')
     +		read(buf,*) Plots.HalfwidthFermi
	if (str.eq.'RenderPitch') read(buf,*) Plots.iRenderPitch
	if (str.eq.'SeparateBandPlotFiles')
     +		Plots.SeparateBandPlotFiles=GetLogical(buf,.FALSE.)	!!! "... = 0" fails!
	enddo

*-------------------------------------------------------------------------

	FP.iUseFP=0
	FP.iVShiftTo=0

	ifs=iIniFileOpenSection(fn,'FullPotential',25,0)
	dowhile (iIniFileReadEntryString(ifs,buf).gt.0)
	m=index(buf,'=')
	if (m.le.1) continue
	str=buf(1:m-1)
	buf=buf(m+1:)
        call BTrimEx(buf,256)
	if (str.eq.'UseFullPotential') read(buf,*) FP.iUseFP
	if (str.eq.'PotentialShiftTo') then
	  i=index('0Hh1Ll2Aa',buf(1:1))
	  if (i.eq.0.and.logf.ne.0) write(logf,*)
     +	  'Warning: Unknown type of FullPotential.PotentialShiftTo'
	  if (i.gt.0) FP.iVShiftTo=(i-1)/3
d	  print *,'FP.iVShiftTo:',FP.iVShiftTo
	endif
	enddo

c	if (FP.iUseFP.ne.0) then	!-----

	FP.iCoordSystem=2
	FP.iNumPointsRho=10
	FP.iNumPointsTheta=10
	FP.iNumPointsPhi=10
	FP.iNumPointsX=10
	FP.iNumPointsY=10
	FP.iNumPointsZ=10
	FP.iNcell=1

	ifs=iIniFileOpenSection(fn,'FullPotential/SpatialMap',25,0)
	dowhile (iIniFileReadEntryString(ifs,buf).gt.0)
	m=index(buf,'=')
	if (m.le.1) continue
	str=buf(1:m-1)
	buf=buf(m+1:)
        call BTrimEx(buf,256)

	if (str.eq.'CoordinateSystem')
     +	read(buf,*) FP.iCoordSystem
	if (str.eq.'NumPointsX') read(buf,*) FP.iNumPointsX
	if (str.eq.'NumPointsY') read(buf,*) FP.iNumPointsY
	if (str.eq.'NumPointsRho') read(buf,*) FP.iNumPointsRho
	if (str.eq.'NumPointsTheta')read(buf,*)FP.iNumPointsTheta
	if (str.eq.'NumPointsPhi') read(buf,*) FP.iNumPointsPhi
	if (str.eq.'NumPointsZ') read(buf,*) FP.iNumPointsZ
	if (str.eq.'NeighbourCells') read(buf,*) FP.iNcell
	enddo

* Write debug info:
	if (logf.ne.0) then
d	write(logf,*) 'Numbers of the FP spatial map points:'
	if (FP.iCoordSystem.eq.iCoords3dCartesian) then
	   write(logf,*) 'NumPoints(X,Y,Z):',
     +		FP.iNumPointsX,FP.iNumPointsY,FP.iNumPointsZ
	elseif (FP.iCoordSystem.eq.iCoords3dSpherical) then
	   write(logf,*) 'NumPoints(Rho,Theta,Phi):',
     +		FP.iNumPointsRho,FP.iNumPointsTheta,FP.iNumPointsPhi
	elseif (FP.iCoordSystem.eq.iCoords3dCylindrical) then
	   write(logf,*) 'NumPoints(Rho,Phi,Z):',
     +		FP.iNumPointsRho,FP.iNumPointsPhi,FP.iNumPointsZ
	endif
	endif

* Calculate the spatial map steps (in a.u.):
	if (Project.iUseCoredTube.ne.2) then
           FP.dX=2.0d0*(A-B)/FP.iNumPointsX
           FP.dY=2.0d0*(A-B)/FP.iNumPointsY
	   FP.dRho=(A-B)/FP.iNumPointsRho
	   FP.dTheta=pi/FP.iNumPointsTheta
	   FP.dPhi=2.0d0*pi/FP.iNumPointsPhi
	   if (FP.iNumPointsZ.eq.1) then
	      FP.dZ=0
	   else
	      FP.dZ=C/FP.iNumPointsZ
	   endif
! When transferring the line above^ to Over55.for, be careful:
! the length of the elementary cell along the z axis ("C")
! is called "CZ" there.
	else
	   write(*,*) 'FP mapping is not ready yet for Inverted LACW'
c	   stop 'FPX: Internal Error [FP/Map/Invert]'
	   iReadProject=0
	   return
	endif

c	endif	!----- if (FP.iUseFP.ne.0)

! count valent electrons:
c	Project.NelTot=0
c	do it=1,nt
c	  if (Project.NelCell.gt.0) then
c	    Project.NelTot=Project.NelTot+nts(it)*Project.NelCell
c	  else
c	    Project.NelTot=Project.NelTot+nts(it)*nel(it)
c	  endif
c	enddo

! copy to records:
	CellStruct.A=A
      CellStruct.B=B
      CellStruct.C=C
      CellStruct.V=V
	if(Project.iUseCoredTube.eq.4)then
        CellStruct.V=VI
	end if
      CellStruct.ChiralityAngle=ChiralityAngle

	iReadProject=1
	return

! failure return:
  215	iReadProject=0
	return

	end
