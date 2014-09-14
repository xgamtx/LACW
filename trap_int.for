********** K *************
	double precision function funBKP(m, a, x, k)
	implicit double precision(a-h,o-z)
	double precision x, k
	integer m
	DIMENSION BI(0:250),DI(0:250),BK(0:250),DK(0:250)

	funBKP = 0.0D0

      if(x*k < 700.0D0)then

      call IKNA(m,x*k,NM,BI,DI,BK,DK)

	funBKP = BK(m)**2*x      

      call IKNA(m,a*k,NM,BI,DI,BK,DK)

	funBKP = funBKP / BK(m)**2

	end if

	return
	end
      
      double precision function trap_int(a, b, k, m)
      double precision a,b,k,eps,s1,s2,delta,funBKP
	integer i, n, m

	eps = 1.0d-1

	n = 10
	s1 = 2*eps
	s2 = 0
	dowhile(dabs(s1-s2)>eps)
	  delta = (b-a)/n
	  s2 = s1
	  s1 = 0
	  do i=0,n-1
	    s1 = s1 + delta*(
     *      funBKP(m, a, a+dble(i)*delta, k)+
	*      funBKP(m, a, a+dble(i+1)*delta, k))/dble(2)
	  enddo
        n=n*2
	enddo

	trap_int = s1

	return
	end

      double precision function trap_int_inf(a, k, m)
      double precision a,k,eps,s1,s2,trap_int
	integer n, m

	eps = 1.0d-0

	n = 0
	s1 = 2*eps
	s2 = 0
	dowhile(dabs(s1-s2)>eps)
	  s2 = s1
        s1 = trap_int(a, a*4.0d0+dble(n), k, m)
        n=n + a*2
	enddo

	trap_int_inf = s1

	return
	end


********** Y *************

	double precision function funBYP(m, a, x, k, eps)
	implicit double precision(a-h,o-z)
	double precision x, k, eps, bya, byap
	integer m
      
	call BessYproc(m,x*k,eps,bya,byap)

	funBYP = bya**2*x      

	call BessYproc(m,a*k,eps,bya,byap)

	funBYP = funBYP / bya**2

	return
	end
      
      double precision function trap_int_y(a, b, k, m)
      double precision a,b,k,eps,s1,s2,delta,funBYP
	integer i, n, m

	eps = 1.0d-1

	n = 10
	s1 = 2*eps
	s2 = 0
	dowhile(dabs(s1-s2)>eps)
	  delta = (b-a)/n
	  s2 = s1
	  s1 = 0
	  do i=0,n-1
	    s1 = s1 + delta*(
     *      funBYP(m, a, a+dble(i)*delta, k, eps)+
	*      funBYP(m, a, a+dble(i+1)*delta, k, eps))/dble(2)
	  enddo
        n=n*2
	enddo

	trap_int_y = s1

	return
	end

      double precision function trap_int_inf_y(a, k, m)
      double precision a,k,eps,s1,s2,trap_int_y
	integer n, m

	eps = 1.0d-0

	n = 0
	s1 = 2*eps
	s2 = 0
	dowhile(dabs(s1-s2)>eps)
	  s2 = s1
        s1 = trap_int_y(a, a*4.0d0+dble(n), k, m)
        n=n + a*2
	enddo

	trap_int_inf_y = s1

	return
	end

********** I *************
	double precision function funBIP(m, b, x, k)
	implicit double precision(a-h,o-z)
	double precision x, k
	integer m
	DIMENSION BI(0:250),DI(0:250),BK(0:250),DK(0:250)
      
      call IKNA(m,x*k,NM,BI,DI,BK,DK)

	funBIP = BI(m)**2*x      

      call IKNA(m,b*k,NM,BI,DI,BK,DK)

	funBIP = funBIP / BI(m)**2

	return
	end
      
      double precision function trap_intI(a, b, k, m)
      double precision a,b,k,eps,s1,s2,delta,funBIP
	integer i, n, m

	eps = 1.0d-1

	n = 10
	s1 = 2*eps
	s2 = 0
	dowhile(dabs(s1-s2)>eps)
	  delta = (b-a)/n
	  s2 = s1
	  s1 = 0
	  do i=0,n-1
	    s1 = s1 + delta*(
     *      funBIP(m, b, a+dble(i)*delta, k)+
	*      funBIP(m, b, a+dble(i+1)*delta, k))/dble(2)
	  enddo
        n=n*2
	enddo

	trap_intI = s1

	return
	end

      double precision function trap_int_infI(b, k, m)
      double precision b,k,trap_intI
	integer m

      trap_int_infI = trap_intI(0.0d0, b, k, m)

	return
	end
