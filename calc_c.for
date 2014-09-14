	double precision function CalcCJ(m, a, b, k, v)
	implicit double precision(a-h,o-z)
	double precision a, b, k, v, eps, kcr
	integer m, s_m
	double precision trap_int_inf, tii

	eps=1.0d-12

      s_m = iabs(m)
	
      if(v.gt.k*k)then

	kcr = dsqrt(v - k*k)

	call Bess(s_m,k*a,bja,bjap)
	call Bess(s_m,k*b,bjb,bjbp)
	call BessYproc(s_m,k*a,eps,bya,byap)
	call BessYproc(s_m,k*b,eps,byb,bybp)

      tii = trap_int_inf(a, kcr, s_m)

      CalcCJ = (dabs(a**2/2.0d0*(bjap-bjb*byap/byb)**2-
	*b**2/2.0d0*(bjbp-bjb*bybp/byb)**2+
	*(bja-bjb*bya/byb)**2*tii))**(-0.5d0)

	else if(k*k.gt.v)then

	kcr = dsqrt(k*k - v)

	call Bess(s_m,k*a,bja,bjap)
	call Bess(s_m,k*b,bjb,bjbp)
	call BessYproc(s_m,k*a,eps,bya,byap)
	call BessYproc(s_m,k*b,eps,byb,bybp)

      tii = trap_int_inf_y(a, kcr, s_m)

      CalcCJ = (dabs(a**2/2.0d0*(bjap-bjb*byap/byb)**2-
	*b**2/2.0d0*(bjbp-bjb*bybp/byb)**2+
	*(bja-bjb*bya/byb)**2*tii))**(-0.5d0)

	else
	  stop '*ERR*  V = X*X'
	end if

	return
	end


	double precision function CalcCY(m, cj, b, k)
	implicit double precision(a-h,o-z)
	double precision b, k, eps
	integer m, s_m

	eps=1.0d-12

      s_m = iabs(m)

	call Bess(s_m,k*b,bjb,tmp)
	call BessYproc(s_m,k*b,eps,byb,tmp)

      CalcCY = -cj*bjb/byb

	return
	end

	double precision function CalcCCR(m, cj, a, b, k, v)
	implicit double precision(a-h,o-z)
	double precision a, b, k, v, eps
	integer m, s_m
	double precision trap_int_inf, tii
	DIMENSION BI(0:250),DI(0:250),BK(0:250),DK(0:250)

	eps=1.0d-12

	!xcr = dsqrt(v - k*k)

      s_m = iabs(m)

	call Bess(s_m,k*a,bja,bjap)
	call Bess(s_m,k*b,bjb,tmp)
	call BessYproc(s_m,k*a,eps,bya,byap)
	call BessYproc(s_m,k*b,eps,byb,tmp)

      !call IKNA(s_m,xcr*a,NM,BI,DI,BK,DK)

      CalcCCR = cj*(bja-(bjb*bya/byb))!/BK(s_m)

	return
	end
