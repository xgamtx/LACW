* TEST routine:
c      implicit integer(a-z)
c      parameter (MAXMNP=40,MUL=8)
c      dimension INDXX(3,MUL*MAXMNP)
c      countx=iReadExplicitBasis(INDXX,3,MUL*MAXMNP,Mmax,Nmax,Pmax)
c      print *,'Total functions:',countx
c      print *,'Max abs values of m,n,p:',Mmax,Nmax,Pmax
c      do i=1,countx
c         print *,i,':',INDXX(1,i),INDXX(2,i),INDXX(3,i)
c      enddo
c      stop
c      end

*--------------------------------------------------------------------
* int iReadExplicitBasis(OUT int INDXX[][], IN int INDXX_dim1,MAX,
*			 OUT int Mmax,Nmax,Pmax);
* Reads the (m,n,p) numbers (maximum = MAX) from the 'selectb.dat' file
* into INDXX; sets Mmax,Nmax,Pmax = max|m|,...
* INDXX is an array [1..INDXX_dim1,MAX] (INDXX_dim1>=3 to store (m,n,p)
*
      integer function iReadExplicitBasis(INDXX,INDXX_dim1,MAX,
     . Mmax,Nmax,Pmax)
      implicit integer(a-z)
      character buf*(80)
      dimension INDXX(INDXX_dim1,MAX)

      ret=0
      Mmax=0
      Nmax=0	! =1?
      Pmax=0
      open(33,file='selectb.dat',status='old',
     . form='formatted',iostat=res)
      if (res.ne.0) then
         close(33,status='delete')
         goto 213	! failed to open the file
      endif
      do while (.true.)
         buf=' '
         do while (lentrim(buf).eq.0.or.index(buf,'#').gt.0)
            read(33,'(a80)',iostat=res) buf
            if (res.ne.0) exit
         enddo
         if (res.ne.0) exit	! positive on error, negative on EOF
d         print *,'"',buf(1:lentrim(buf)),'"'
         read(unit=buf,fmt=210,iostat=res) m,n,p
  210    format(3i8)
         if (res.ne.0) continue
         is_duplicate=0
         do i=1,ret
            if (INDXX(1,i).eq.m.and.INDXX(2,i).eq.n.and.
     .          INDXX(3,i).eq.p) is_duplicate=1
         enddo
         if (is_duplicate.ne.0) then	! already know it
            print *,'*WRN* Duplicate function',m,n,p
            goto 211
         endif
         ret=ret+1
         INDXX(1,ret)=m
         INDXX(2,ret)=n
         INDXX(3,ret)=p
         if (iabs(m).gt.Mmax) Mmax=iabs(m)
         if (iabs(n).gt.Nmax) Nmax=iabs(n)
         if (iabs(p).gt.Pmax) Pmax=iabs(p)
  211    if (ret.ge.MAX) exit
      enddo
  212 close(33)

  213 iReadExplicitBasis=ret
      return
      end
