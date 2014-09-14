      SUBROUTINE CHLSK1C(NN,N,S,DIAG)
c     Here: NN - dimension of H() and S(), N - effective count of elements --
      INTEGER I,J,K,N,NN
      DOUBLE COMPLEX S(NN,NN), FLOWA
      DOUBLE PRECISION FLOWB, DIAG(NN)
*
      DO 10 J=1,N
      FLOWB=DBLE(S(J,J))
      DO 20 K=1,(J-1)
*
*     lower triangle
*
  20  FLOWB=FLOWB-DBLE(S(J,K)*DCONJG(S(J,K)))
*
      if(flowb.le.0.D0) then
         WRITE(*,501) j,n,dble(flowb)
 501     format(5x,'*** J',i5,' N',i5,' flowb',f20.10)
	 stop 501
      endif
      DIAG(J)=1.D0/DSQRT(FLOWB)
      S(J,J) = DCMPLX(DSQRT(FLOWB),0.D0)
*
      DO 10 I=(J+1),N
      FLOWA=dconjg(S(J,I))
      DO 2 K=1,(J-1)
*                     - s1 x dconj( s2)
      FLOWA=FLOWA-S(I,K)*DCONJG(S(J,K))
    2 CONTINUE
*
*     lower triangle
*
      S(I,J) = FLOWA*DCMPLX(DIAG(J),0.D0)
   10 CONTINUE
      RETURN
      END
*
      SUBROUTINE CHLSK2C(NN,N,H,S,DIAG)
      INTEGER I,J,K,N,NN
      DOUBLE COMPLEX H(NN,NN),S(NN,NN),FLOWA
      DOUBLE PRECISION DIAG(NN)
*
      DO 1 J=1,N
      DO 1 I=J,N
*
*     lower triangle
*
      FLOWA=H(I,J)
      DO 2 K=1,J-1
*                - h x dconj( s )
      FLOWA=FLOWA-H(K,I)*DCONJG(S(J,K))
    2 CONTINUE
*
*     upper triangle
*
      H(J,I) =FLOWA*DCMPLX(DIAG(J),0.D0)
    1 CONTINUE
      RETURN
      END
*
      SUBROUTINE CHLSK3C(NN,N,H,S,DIAG)
      INTEGER I,J,K,N,NN
      DOUBLE COMPLEX H(NN,NN),S(NN,NN),FLOWA
      DOUBLE PRECISION DIAG(NN)
*
      DO 2 J=1,N
      FLOWA=H(J,J)
      DO 3 K=1,J-1
*                - s x dconj( h )
      FLOWA=FLOWA-S(J,K)*DCONJG(H(J,K))
    3 CONTINUE
      H(J,J) =FLOWA*DCMPLX(DIAG(J),0.D0)
*
      DO 2 I=(J+1),N
*
*     upper triangle
*
         FLOWA=H(J,I)
         DO 4 K=1,J-1
*                   - s x dconj( h )
         FLOWA=FLOWA-S(I,K)*DCONJG(H(J,K))
    4    CONTINUE
         DO 5 K=J,I-1
*                   - s x h
         FLOWA=FLOWA-S(I,K)*H(K,J)
    5    CONTINUE
*
*     lower triangle
*
         H(I,J) =FLOWA*DCMPLX(DIAG(I),0.D0)
    2 CONTINUE
*
      RETURN
      END
*
      SUBROUTINE REIGNR_c(NN,N,NLMQ,VEC,S,DIAG)
      INTEGER I,J,K,N,NN,NLMQ
      DOUBLE COMPLEX VEC(NN,*),S(NN,*),FLOWA
      DOUBLE PRECISION DIAG(*)
*
*     invert L by solving   L * L(-1) = 1
*
      DO 1 J=1,N-1
      DO 1 I=J+1,N
      FLOWA=-S(I,J)*DCMPLX(DIAG(J),0.d0)
      DO 2 K=J+1,I-1
*                - s x pl
      FLOWA=FLOWA-S(I,K)*S(J,K)
    2 CONTINUE
*
*     upper triangle
*
      S(J,I) =FLOWA*DCMPLX(DIAG(I),0.d0)
    1 CONTINUE
*
*     X = L(-H) * Y
*
      DO 3 J=1,NLMQ
      DO 3 I=1,N
      FLOWA=DCMPLX(DIAG(I),0.d0)*VEC(I,J)
      DO 4 K=I+1,N
*                + dconj( pl ) x vec
      FLOWA=FLOWA+DCONJG(S(I,K))*VEC(K,J)
    4 CONTINUE
      VEC(I,J) = FLOWA
    3 CONTINUE
      RETURN
      END
