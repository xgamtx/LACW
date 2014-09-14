      FUNCTION plgndr(l,m,x)
*
*     by IMSL
*
      INTEGER*2 l,m
      double precision plgndr,x
      INTEGER*2 i,ll
      double precision fact,pll,pmm,pmmp1,somx2
      if(m.lt.0.or.m.gt.l.or.dabs(x).gt.1.) then 
	   STOP 'Bad arguments in function PLGNDR'
	endif
      pmm=1.
      if(m.gt.0) then
        somx2=dsqrt((1.-x)*(1.+x))
        fact=1.
        do 11 i=1,m
          pmm=-pmm*fact*somx2
          fact=fact+2.d0
11      continue
      endif
      if(l.eq.m) then
        plgndr=pmm
      else
        pmmp1=x*(2*m+1)*pmm
        if(l.eq.m+1) then
          plgndr=pmmp1
        else
          do 12 ll=m+2,l
            pll=(x*(2*ll-1)*pmmp1-(ll+m-1)*pmm)/(ll-m)
            pmm=pmmp1
            pmmp1=pll
12        continue
          plgndr=pll
        endif
      endif
      return
      END











