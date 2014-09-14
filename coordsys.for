*-------------------------------------------------------------------------
	double precision function dist3d(x1,y1,z1,x2,y2,z2)
	implicit double precision(a-h,o-z)

	dx=x2-x1
	dy=y2-y1
	dz=z2-z1
	dist3d=dsqrt(dx*dx+dy*dy+dz*dz)

	return
	end

*-------------------------------------------------------------------------
	double precision function dist3dx(a1,a2,a3,b1,b2,b3,iCoordSys)
	implicit double precision(a-h,o-z)

	include 'coordsys.fh'

	call ToCartesian3D(a1,a2,a3,x1,y1,z1,iCoordSys)
	call ToCartesian3D(b1,b2,b3,x2,y2,z2,iCoordSys)
	dx=x2-x1
	dy=y2-y1
	dz=z2-z1
	dist3dx=dsqrt(dx*dx+dy*dy+dz*dz)

	return
	end

*-------------------------------------------------------------------------
	subroutine ToCartesian3D(a1,a2,a3,x,y,z,iCoordSys)
	implicit double precision(a-h,o-z)

	include 'coordsys.fh'
c	parameter (iCoords3dCartesian=0,iCoords3dSpherical=1,
c     +	iCoords3dCylindrical=2)

	if (iCoordSys.eq.iCoords3dCartesian) then	! a123=x,y,z
	   x=a1
	   y=a2
	   z=a3
	elseif (iCoordSys.eq.iCoords3dSpherical) then	! a123=rho,theta,phi
	   x=a1*dcos(a2)*dcos(a3)
	   y=a1*dcos(a2)*dsin(a3)
	   z=a1*dsin(a2)
	elseif (iCoordSys.eq.iCoords3dCylindrical) then	! a123=rho,phi,z
	   x=a1*dcos(a2)
	   y=a1*dsin(a2)
	   z=a3
	endif

	return
	end
*-------------------------------------------------------------------------
	double precision function dVolumeCoeff(iCoordSys,a1,a2,a3)
* Given the coordinates of a point in 3D space, dVolumeCoeff returns
* the coefficient in dV=(?)*da1*da2*da3, where (?) is 1 for Cartesian,
* r*r*sin(theta) for spherical, and r for cylidrical coordinates.
	implicit double precision (a-h,o-z)
	integer iCoordSys
	double precision a1,a2,a3
	include 'coordsys.fh'

	if (iCoordSys.eq.iCoords3dCartesian) then	! a123=x,y,z
	  dVolumeCoeff=1.0d0
	elseif (iCoordSys.eq.iCoords3dSpherical) then	! a123=rho,theta,phi
	  dVolumeCoeff=a1*a1*dsin(a2)
	elseif (iCoordSys.eq.iCoords3dCylindrical) then	! a123=rho,phi,z
	  dVolumeCoeff=a1
	else
	  stop 'Fatal Error [CoordSys.for]: Unknown coordinate system'
	endif

	return
	end
