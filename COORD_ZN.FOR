	INTEGER N,NC
	DOUBLE PRECISION A,alpha,PI,r,z
*
	PI=4.D0*DATAN(1.d0)
	READ(*,*) N,A,NC
	alpha=2.*PI/n
!?	r=a*dcos(PI/6.)/(alpha/2.D0)
	r=a*dsqrt(3.D0)/4.D0/dsin(alpha/4.D0)
	do i=1,NC
	  z=a*dble(3*(i-1))
	  CALL Outp(0.D0,z+0.D0,r,alpha,n)
	  CALL Outp(0.D0,z+a,r,alpha,n)
	  CALL Outp(alpha/2.D0,z+1.5d0*a,r,alpha,n)
	  CALL Outp(alpha/2.D0,z+2.5d0*a,r,alpha,n)
	enddo
	END

	SUBROUTINE Outp(angle,mtin,rad,alpha,n)
	INTEGER N,I
	DOUBLE PRECISION a,alpha,angle,mtin,rad
	a=angle
	DO I=1,n
	  print '(3(5X,F15.10))',rad*dsin(a),rad*dcos(a),mtin
	  a=a+alpha
	END DO
	RETURN
	END
