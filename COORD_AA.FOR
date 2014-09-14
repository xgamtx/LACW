	INTEGER N,NC
	DOUBLE PRECISION A,alpha,PI,r,z
*
	PI=4.D0*DATAN(1.D0)
	READ(*,*) N,A,NC
	alpha=2.*PI/3./n
	n=n*2
c	r=a/alpha
	r=a/(2.0d0*dsin(alpha/2.0d0))
	z=a*dcos(PI/6.0d0)
	do i=1,NC
	  CALL Outp(0.D0,z*dble(2*i-2),r,alpha,n,0)
	  CALL Outp(-alpha/2.D0,z*dble(2*i-1),r,alpha,n,1)
	enddo
	END

	SUBROUTINE Outp(angle,mtin,rad,alpha,n,ILFlag)
	INTEGER N,I
	DOUBLE PRECISION a,alpha,angle,mtin,rad
	a=angle
	DO I=1,n
	  print '(3(5X,F15.10))',rad*dsin(a),rad*dcos(a),mtin
	  IF(MOD(I+ILFlag,2)) THEN
	  a=a+alpha
	   ELSE
	    A=A+2.D0*alpha
	  ENDif
	END DO
	RETURN
	END
