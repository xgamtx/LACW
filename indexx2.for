      SUBROUTINE indexx2(nn,mn,n,nc,arr,indx)
*
c      Sorting of indexes of an array, stored in ARR.
c      By IMSL.
c      mn here is "effective" size of arrays arr() and indx(),
c      nc - number of string in arr for sorting process.
*
      INTEGER nn,mn,n,nc,indx(*),M,NSTACK
      DOUBLE PRECISION arr(nn,n)
      PARAMETER (M=7,NSTACK=50)
      INTEGER i,indxt,ir,itemp,j,jstack,k,l,istack(NSTACK)
      DOUBLE PRECISION a
      do 11 j=1,mn
        indx(j)=j
11    continue
      jstack=0
      l=1
      ir=mn
1     if(ir-l.lt.M)then
        do 13 j=l+1,ir
          indxt=indx(j)
          a=arr(nc,indxt)
          do 12 i=j-1,l,-1
            if(arr(nc,indx(i)).le.a)goto 2
            indx(i+1)=indx(i)
12        continue
          i=l-1
2         indx(i+1)=indxt
13      continue
        if(jstack.eq.0)return
        ir=istack(jstack)
        l=istack(jstack-1)
        jstack=jstack-2
      else
        k=(l+ir)/2
        itemp=indx(k)
        indx(k)=indx(l+1)
        indx(l+1)=itemp
        if(arr(nc,indx(l)).gt.arr(nc,indx(ir)))then
          itemp=indx(l)
          indx(l)=indx(ir)
          indx(ir)=itemp
        endif
        if(arr(nc,indx(l+1)).gt.arr(nc,indx(ir)))then
          itemp=indx(l+1)
          indx(l+1)=indx(ir)
          indx(ir)=itemp
        endif
        if(arr(nc,indx(l)).gt.arr(nc,indx(l+1)))then
          itemp=indx(l)
          indx(l)=indx(l+1)
          indx(l+1)=itemp
        endif
        i=l+1
        j=ir
        indxt=indx(l+1)
        a=arr(nc,indxt)
3       continue
          i=i+1
        if(arr(nc,indx(i)).lt.a)goto 3
4       continue
          j=j-1
        if(arr(nc,indx(j)).gt.a)goto 4
        if(j.lt.i)goto 5
        itemp=indx(i)
        indx(i)=indx(j)
        indx(j)=itemp
        goto 3
5       indx(l+1)=indx(j)
        indx(j)=indxt
        jstack=jstack+2
        if(jstack.gt.NSTACK) stop 'NSTACK too small in indexx'
        if(ir-i+1.ge.j-l)then
          istack(jstack)=ir
          istack(jstack-1)=i
          ir=j-1
        else
          istack(jstack)=j-1
          istack(jstack-1)=l
          l=i
        endif
      endif
      goto 1
      END
