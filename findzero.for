*====================================================================
* FindZero of user function F(x) in interval [a,b] with accuracy eps;
* set eps=0 to calculate with maximum accuracy ("machine epsilon").
* Return value: an x in [a,b], F(x)=0; or +INF if not found.
*....................................................................

      double precision function FindZero(F,a,b,eps)

      double precision a,b,eps,ya,yb,x,y,small,big
      double precision F
      external F
      data SMALL/1.0d-36/, BIG/1.0d36/

d      print *,'In function FindZero()...'
      do while (dabs(b-a).gt.eps)
         ya=f(a)
         yb=f(b)
         x=(a+b)/2.0d0
         if (x.eq.a) exit   ! for the case when eps==0
d         if (dabs(x-a).lt.SMALL) exit
         y=f(x)
d         print *,'axb,ya,y,yb:',a,x,b,ya,y,yb
         if ((y.lt.SMALL.and.ya.gt.SMALL).or.
     .       (y.gt.SMALL.and.ya.lt.SMALL)) then
            b=x
         elseif ((y.lt.SMALL.and.yb.gt.SMALL).or.
     .           (y.gt.SMALL.and.yb.lt.SMALL)) then
            a=x
         else
            FindZero=BIG
            return
         endif
d         print *,'eps: |a-b| = ',dabs(a-b)
      enddo
      FindZero=x
      return
      end



*-------------------------------------------------------------------------
      double precision function FindY(xs,ys,n,x)
*     ----------------------------------------------------------
*     Finds y(x), where the function is passed as SORTED arrays
*     of points (x,y) -- xs[n] and ys[n]
*     ----------------------------------------------------------

      implicit double precision(a-h,o-z)
      dimension xs(*),ys(*)

      data dmin/1.0d-6/

      FindY=x
      if (n.le.0) return

      m=1
      do i=1,n-1
         if (x.ge.xs(i).and.x.lt.xs(i+1)) m=i
      enddo
c      if (x.lt.xs(1)) m=1
      if (n.gt.1.and.x.ge.xs(n)) m=n-1

      ret=ys(m)
      if (dabs(xs(m+1)-xs(m)).ge.dmin)
     . ret=ys(m)+(ys(m+1)-ys(m))*(x-xs(m))/(xs(m+1)-xs(m))

      FindY=ret
      return
      end     !//FindY()
