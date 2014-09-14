*====================================================================
* Module origname:        Lentrim
* Module specs:           String trimming function
* Project:                Units/LATW
* Functions:              1
* Style warnings:         0
*====================================================================
* TODO:                   0
*====================================================================

* LINK WITH: -


*--------------------------------------------------------------------
*       int iMyLentrimEx(char *s, int maxlen);
*
* Returns the length of a character string of allocated length maxlen
* with trailing space characters (0..32) removed.
*--------------------------------------------------------------------
      integer function iMyLentrimEx(s,nmax)

      character s*(*)
      integer nmax,i

d      print *,'s="',s(1:nmax),'"',nmax
      do i=nmax,1,-1
         if (s(i:i+1).gt.char(32)) exit
      enddo

d      print *,'returning',i
      iMyLentrimEx=i
      return
      end

*====================================================================
* The lentrim() function is implemented in Watcom Fortran 77, but is
* absent in MS Fortran Powerstation.
*--------------------------------------------------------------------
c$ifdef MSFPS
      integer function lentrim(s)
      character*80 s

      n=iMyLentrimEx(s,80)
d      print *,'s="',s(1:n),'"',n

      lentrim=n
      return
      end
c$endif

	subroutine BTrimEx(s,nmax)
	character s*(*)
	integer nmax
	do i=1,iMyLentrimEx(s,nmax)
	  if (s(i:i).le.char(32)) s=s(1:i-1)//char(32)//s(i+1:)
	enddo
	dowhile (iMyLentrimEx(s,nmax).gt.0.and.s(1:1).le.char(32))
	  s=s(2:)
	enddo
	return
	end
