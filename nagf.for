      SUBROUTINE F06AAZ ( SRNAME, INFO )
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     MARK 15 REVISED. IER-915 (APR 1991).
C     .. Scalar Arguments ..
      INTEGER            INFO
      CHARACTER*13       SRNAME
C     ..
C
C  Purpose
C  =======
C
C  F06AAZ  is an error handler for the Level 2 BLAS routines.
C
C  It is called by the Level 2 BLAS routines if an input parameter is
C  invalid.
C
C  Parameters
C  ==========
C
C  SRNAME - CHARACTER*13.
C           On entry, SRNAME specifies the name of the routine which
C           called F06AAZ.
C
C  INFO   - INTEGER.
C           On entry, INFO specifies the position of the invalid
C           parameter in the parameter-list of the calling routine.
C
C
C  Auxiliary routine for Level 2 Blas.
C
C  Written on 20-July-1986.
C
C     .. Local Scalars ..
      INTEGER            IERR, IFAIL
      CHARACTER*4        VARBNM
C     .. Local Arrays ..
      CHARACTER*80       REC (1)
C     .. External Functions ..
      INTEGER            P01ACF
      EXTERNAL           P01ACF
C     ..
C     .. Executable Statements ..
      WRITE (REC (1),99999) SRNAME, INFO
      IF (SRNAME(1:3).EQ.'F06') THEN
         IERR = -1
         VARBNM = '    '
      ELSE
         IERR = -INFO
         VARBNM = 'INFO'
      END IF
      IFAIL = 0
      IFAIL = P01ACF (IFAIL, IERR, SRNAME(1:6), VARBNM, 1, REC)
C
      RETURN
C
99999 FORMAT ( ' ** On entry to ', A13, ' parameter number ', I2,
     $         ' had an illegal value' )
C
C     End of F06AAZ.
C
      END
*
      DOUBLE PRECISION FUNCTION F06BNF( A, B )
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     .. Scalar Arguments ..
      DOUBLE PRECISION                  A, B
C     ..
C
C  F06BNF returns the value
C
C     p = sqrt( a*a + b*b )
C
C  via the function name.
C
C
C  Nag Fortran 77 O( 1 ) basic linear algebra routine.
C
C  -- Written on 17-January-1985.
C     Sven Hammarling, Nag Central Office.
C
C
C     .. Parameters ..
      DOUBLE PRECISION      ZERO
      PARAMETER           ( ZERO = 0.0D+0 )
C     .. Local Scalars ..
      DOUBLE PRECISION      P
C     .. Intrinsic Functions ..
      INTRINSIC             ABS, SQRT
C     ..
C     .. Executable Statements ..
      IF( A.EQ.ZERO )THEN
         P = ABS( B )
      ELSE IF( B.EQ.ZERO )THEN
         P = ABS( A )
      ELSE IF( ABS( A ).GE.ABS( B ) )THEN
         P = ABS( A )*SQRT( 1 + ( B/A )**2 )
      ELSE
         P = ABS( B )*SQRT( 1 + ( A/B )**2 )
      END IF
C
      F06BNF = P
      RETURN
C
C     End of F06BNF. ( SPYTH )
C
      END
*
      COMPLEX*16       FUNCTION F06CLF( A, B, FAIL )
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     MARK 13 REVISED. IER-605 (MAR 1988).
C     .. Scalar Arguments ..
      COMPLEX*16                        A, B
      LOGICAL                           FAIL
C     ..
C
C  F06CLF returns the value div given by
C
C     div = ( a/b      if a/b does not overflow,
C           (
C           ( 0.0      if a .eq. 0.0,
C           (
C           ( cflmax   if a .ne. 0.0 and a/b would overflow,
C
C  where
C
C     cflmax = ( flmax*sign( re( a/b ) ), flmax*sign( im( a/b ) ) )
C
C  and flmax is a large value, via the function name. In addition if a/b
C  would  overflow then  fail  is returned as  true, otherwise  fail  is
C  returned as false.
C
C  Note that when a and b are both zero, fail is returned as .true., but
C  div  is returned as  0.0. In all other cases of overflow  div is such
C  that abs( re( div ) ) and abs( im( div ) ) are flmax.
C
C  For  real  x and y,  if  y = 0,  sign( x/y )  is taken as  sign( x ).
C
C  Nag Fortran 77 O( 1 ) basic linear algebra routine.
C
C  -- Written on 27-April-1983.
C     Sven Hammarling, Nag Central Office.
C  -- Amended on 4-December-1987.
C     Sven Hammarling and Mick Pont, Nag Central Office.
C        To avoid extremely unlikely division by zero.
C
C
C     .. Parameters ..
      DOUBLE PRECISION         ONE
      PARAMETER              ( ONE  = 1.0D+0 )
      COMPLEX*16               ZERO
      PARAMETER              ( ZERO = ( 0.0D+0, 0.0D+0 ) )
C     .. Local Scalars ..
      COMPLEX*16               VALUE
      DOUBLE PRECISION         AI, AR, BI, BIG, BR, DIV, FLMAX, FLMIN,
     $                         NUMI, NUMR, TEMP
      LOGICAL                  FIRST
C     .. External Functions ..
      DOUBLE PRECISION         X02AMF
      EXTERNAL                 X02AMF
C     .. Intrinsic Functions ..
      INTRINSIC                ABS, DBLE, DCMPLX, DIMAG, MAX, SIGN
C     .. Save statement ..
      SAVE                     BIG, FIRST, FLMAX
C     .. Data statements ..
      DATA                     FIRST/ .TRUE. /
C     ..
C     .. Executable Statements ..
      IF( A.EQ.ZERO )THEN
         VALUE = ZERO
         IF( B.EQ.ZERO )THEN
            FAIL = .TRUE.
         ELSE
            FAIL = .FALSE.
         END IF
      ELSE
C
         IF( FIRST )THEN
            FIRST = .FALSE.
            FLMIN =  X02AMF( )
            FLMAX =  1/FLMIN
            BIG   =  FLMAX/2
         END IF
C
         AR    =  DBLE ( A )
         AI    =  DIMAG( A )
         BR    =  DBLE ( B )
         BI    =  DIMAG( B )
         TEMP  =  MAX( ABS( AR ), ABS( AI ), ABS( BR ), ABS( BI ) )
         IF( TEMP.GE.BIG )THEN
            AR = AR/2
            AI = AI/2
            BR = BR/2
            BI = BI/2
         END IF
         IF( DCMPLX( BR, BI ).EQ.ZERO ) THEN
            VALUE =  DCMPLX( SIGN( FLMAX, DBLE ( A ) ),
     $                       SIGN( FLMAX, DIMAG( A ) )  )
            FAIL  = .TRUE.
         ELSE
            IF( ABS( BR ).GE.ABS( BI ) )THEN
               TEMP = BI/BR
               DIV  = BR     + TEMP*BI
               NUMR = AR     + TEMP*AI
               NUMI = AI     - TEMP*AR
            ELSE
               TEMP = BR/BI
               DIV  = BI      + TEMP*BR
               NUMR = AI      + TEMP*AR
               NUMI = TEMP*AI - AR
            END IF
            IF( ABS( DIV ).GE.ONE )THEN
               VALUE =  DCMPLX( NUMR/DIV, NUMI/DIV )
               FAIL  = .FALSE.
            ELSE
               TEMP  =  ABS( DIV )*FLMAX
               IF( ( ABS( NUMR ).LE.TEMP ).AND.
     $             ( ABS( NUMI ).LE.TEMP )      )THEN
                  VALUE =  DCMPLX( NUMR/DIV, NUMI/DIV )
                  FAIL  = .FALSE.
               ELSE
                  IF( DIV.GE.DBLE( ZERO ) )THEN
                     VALUE = DCMPLX( SIGN( FLMAX,  NUMR ),
     $                               SIGN( FLMAX,  NUMI )  )
                  ELSE
                     VALUE = DCMPLX( SIGN( FLMAX, -NUMR ),
     $                               SIGN( FLMAX, -NUMI )  )
                  END IF
                  FAIL = .TRUE.
               END IF
            END IF
         END IF
      END IF
C
      F06CLF = VALUE
      RETURN
C
C     End of F06CLF. ( CDIV )
C
      END
*
      SUBROUTINE F06EDF( N, ALPHA, X, INCX )
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     .. Entry Points ..
      ENTRY      DSCAL ( N, ALPHA, X, INCX )
C     .. Scalar Arguments ..
      DOUBLE PRECISION   ALPHA
      INTEGER            INCX, N
C     .. Array Arguments ..
      DOUBLE PRECISION   X( * )
C     ..
C
C  F06EDF performs the operation
C
C     x := alpha*x
C
C
C  Nag Fortran 77 version of the Blas routine DSCAL.
C  Nag Fortran 77 O( n ) basic linear algebra routine.
C
C  -- Written on 26-November-1982.
C     Sven Hammarling, Nag Central Office.
C
C
C     .. Parameters ..
      DOUBLE PRECISION   ONE         , ZERO
      PARAMETER        ( ONE = 1.0D+0, ZERO = 0.0D+0 )
C     .. Local Scalars ..
      INTEGER            IX
C     ..
C     .. Executable Statements ..
      IF( N.GT.0 )THEN
         IF( ALPHA.EQ.ZERO )THEN
            DO 10, IX = 1, 1 + ( N - 1 )*INCX, INCX
               X( IX ) = ZERO
   10       CONTINUE
         ELSE IF( ALPHA.EQ.( -ONE ) )THEN
            DO 20, IX = 1, 1 + ( N - 1 )*INCX, INCX
               X( IX ) = -X( IX )
   20       CONTINUE
         ELSE IF( ALPHA.NE.ONE )THEN
            DO 30, IX = 1, 1 + ( N - 1 )*INCX, INCX
               X( IX ) = ALPHA*X( IX )
   30       CONTINUE
         END IF
      END IF
C
      RETURN
C
C     End of F06EDF. ( DSCAL )
C
      END
*
      COMPLEX*16       FUNCTION F06GBF( N, X, INCX, Y, INCY )
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     .. Entry Points ..
      COMPLEX*16                ZDOTC
      ENTRY                     ZDOTC ( N, X, INCX, Y, INCY )
C     .. Scalar Arguments ..
      INTEGER                           INCX, INCY, N
C     .. Array Arguments ..
      COMPLEX*16                        X( * ), Y( * )
C     ..
C
C  F06GBF returns the value
C
C     F06GBF = conjg( x' )*y
C
C
C  Nag Fortran 77 version of the Blas routine ZDOTC.
C  Nag Fortran 77 O( n ) basic linear algebra routine.
C
C  -- Written on 21-September-1982.
C     Sven Hammarling, Nag Central Office.
C
C
C     .. Parameters ..
      COMPLEX*16               ZERO
      PARAMETER              ( ZERO = ( 0.0D+0, 0.0D+0 ) )
C     .. Local Scalars ..
      COMPLEX*16               SUM
      INTEGER                  I, IX, IY
C     .. Intrinsic Functions ..
      INTRINSIC                DCONJG
C     ..
C     .. Executable Statements ..
      SUM = ZERO
      IF( N.GT.0 )THEN
         IF( ( INCX.EQ.INCY ).AND.( INCX.GT.0 ) )THEN
            DO 10, IX = 1, 1 + ( N - 1 )*INCX, INCX
               SUM = SUM + DCONJG( X( IX ) )*Y( IX )
   10       CONTINUE
         ELSE
            IF( INCY.GE.0 )THEN
               IY = 1
            ELSE
               IY = 1 - ( N - 1 )*INCY
            END IF
            IF( INCX.GT.0 )THEN
               DO 20, IX = 1, 1 + ( N - 1 )*INCX, INCX
                  SUM = SUM + DCONJG( X( IX ) )*Y( IY )
                  IY  = IY  + INCY
   20          CONTINUE
            ELSE
               IX = 1 - ( N - 1 )*INCX
               DO 30, I = 1, N
                  SUM = SUM + DCONJG( X( IX ) )*Y( IY )
                  IX  = IX  + INCX
                  IY  = IY  + INCY
   30          CONTINUE
            END IF
         END IF
      END IF
C
      F06GBF = SUM
      RETURN
C
C     End of F06GBF. ( ZDOTC )
C
      END
*
      SUBROUTINE F06GCF( N, ALPHA, X, INCX, Y, INCY )
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     .. Entry Points ..
      ENTRY      ZAXPY ( N, ALPHA, X, INCX, Y, INCY )
C     .. Scalar Arguments ..
      COMPLEX*16         ALPHA
      INTEGER            INCX, INCY, N
C     .. Array Arguments ..
      COMPLEX*16         X( * ), Y( * )
C     ..
C
C  F06GCF performs the operation
C
C     y := alpha*x + y
C
C
C  Nag Fortran 77 version of the Blas routine ZAXPY.
C  Nag Fortran 77 O( n ) basic linear algebra routine.
C
C  -- written on 28-April-1983.
C     Sven Hammarling, Nag Central Office.
C
C
C     .. Parameters ..
      COMPLEX*16         ZERO
      PARAMETER        ( ZERO = ( 0.0D+0, 0.0D+0 ) )
C     .. Local Scalars ..
      INTEGER            I, IX, IY
C     .. Executable Statements ..
      IF( N.GT.0 )THEN
         IF( ALPHA.NE.ZERO )THEN
            IF( ( INCX.EQ.INCY ).AND.( INCX.GT.0 ) )THEN
               DO 10, IX = 1, 1 + ( N - 1 )*INCX, INCX
                  Y( IX ) = ALPHA*X( IX ) + Y( IX )
   10          CONTINUE
            ELSE
               IF( INCY.GE.0 )THEN
                  IY = 1
               ELSE
                  IY = 1 - ( N - 1 )*INCY
               END IF
               IF( INCX.GT.0 )THEN
                  DO 20, IX = 1, 1 + ( N - 1 )*INCX, INCX
                     Y( IY ) = ALPHA*X( IX ) + Y( IY )
                     IY      = IY            + INCY
   20             CONTINUE
               ELSE
                  IX = 1 - ( N - 1 )*INCX
                  DO 30, I = 1, N
                     Y( IY ) = ALPHA*X( IX ) + Y( IY )
                     IX      = IX            + INCX
                     IY      = IY            + INCY
   30             CONTINUE
               END IF
            END IF
         END IF
      END IF
C
      RETURN
C
C     End of F06GCF. ( ZAXPY )
C
      END
*
      SUBROUTINE F06GDF( N, ALPHA, X, INCX )
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     .. Entry Points ..
      ENTRY      ZSCAL ( N, ALPHA, X, INCX )
C     .. Scalar Arguments ..
      COMPLEX*16         ALPHA
      INTEGER            INCX, N
C     .. Array Arguments ..
      COMPLEX*16         X( * )
C     ..
C
C  F06GDF performs the operation
C
C     x := alpha*x
C
C
C  Nag Fortran 77 version of the Blas routine ZSCAL.
C  Nag Fortran 77 O( n ) basic linear algebra routine.
C
C  -- Written on 26-November-1982.
C     Sven Hammarling, Nag Central Office.
C
C
C     .. Parameters ..
      COMPLEX*16         CZERO
      PARAMETER        ( CZERO = ( 0.0D+0, 0.0D+0 ) )
C     .. Local Scalars ..
      INTEGER            IX
C     ..
C     .. Executable Statements ..
      IF( N.GT.0 )THEN
         IF( ALPHA.EQ.CZERO )THEN
            DO 10, IX = 1, 1 + ( N - 1 )*INCX, INCX
               X( IX ) = CZERO
   10       CONTINUE
         ELSE
            DO 20, IX = 1, 1 + ( N - 1 )*INCX, INCX
               X( IX ) = ALPHA*X( IX )
   20       CONTINUE
         END IF
      END IF
C
      RETURN
C
C     End of F06GDF. ( ZSCAL )
C
      END
*
      SUBROUTINE F06JDF( N, ALPHA, X, INCX )
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     .. Entry Points ..
      ENTRY      ZDSCAL( N, ALPHA, X, INCX )
C     .. Scalar Arguments ..
      DOUBLE PRECISION   ALPHA
      INTEGER            INCX, N
C     .. Array Arguments ..
      COMPLEX*16         X( * )
C     ..
C
C  F06JDF performs the operation
C
C     x := alpha*x
C
C
C  Nag Fortran 77 version of the Blas routine ZDSCAL.
C  Nag Fortran 77 O( n ) basic linear algebra routine.
C
C  -- Written on 26-November-1982.
C     Sven Hammarling, Nag Central Office.
C
C
C     .. Parameters ..
      DOUBLE PRECISION   ONE
      PARAMETER        ( ONE  = 1.0D+0 )
      COMPLEX*16         ZERO
      PARAMETER        ( ZERO = ( 0.0D+0, 0.0D+0 ) )
C     .. Local Scalars ..
      INTEGER            IX
C     .. Intrinsic Functions ..
      INTRINSIC          DBLE
C     ..
C     .. Executable Statements ..
      IF( N.GT.0 )THEN
         IF( ALPHA.EQ.DBLE( ZERO ) )THEN
            DO 10, IX = 1, 1 + ( N - 1 )*INCX, INCX
               X( IX ) = ZERO
   10       CONTINUE
         ELSE IF( ALPHA.EQ.( -ONE ) )THEN
            DO 20, IX = 1, 1 + ( N - 1 )*INCX, INCX
               X( IX ) = -X( IX )
   20       CONTINUE
         ELSE IF( ALPHA.NE.ONE )THEN
            DO 30, IX = 1, 1 + ( N - 1 )*INCX, INCX
               X( IX ) = ALPHA*X( IX )
   30       CONTINUE
         END IF
      END IF
C
      RETURN
C
C     End of F06JDF. ( ZDSCAL )
C
      END
*
      SUBROUTINE F06KPF( N, X, INCX, Y, INCY, C, S )
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     .. Scalar Arguments ..
      DOUBLE PRECISION   C, S
      INTEGER            INCX, INCY, N
C     .. Array Arguments ..
      COMPLEX*16         X( * ), Y( * )
C     ..
C
C  F06KPF performs the plane rotation
C
C     ( x  y ) = ( x  y )*( c  -s ).
C                         ( s   c )
C
C
C
C  Nag Fortran 77 O( n ) basic linear algebra routine.
C
C  -- Written on 13-February-1986.
C     Sven Hammarling, Nag Central Office.
C
C
C     .. Parameters ..
      DOUBLE PRECISION   ONE         , ZERO
      PARAMETER        ( ONE = 1.0D+0, ZERO = 0.0D+0 )
C     .. Local Scalars ..
      COMPLEX*16         TEMP
      INTEGER            I, IX, IY
C     ..
C     .. Executable Statements ..
      IF( N.GT.0 )THEN
         IF( ( S.NE.ZERO ).OR.( C.NE.ONE ) )THEN
            IF( ( C.EQ.ZERO ).AND.( S.EQ.ONE ) )THEN
               IF( ( INCX.EQ.INCY ).AND.( INCX.GT.0 ) )THEN
                  DO 10, IX = 1, 1 + ( N - 1 )*INCX, INCX
                     TEMP    = -X( IX )
                     X( IX ) =  Y( IX )
                     Y( IX ) =  TEMP
   10             CONTINUE
               ELSE
                  IF( INCY.GE.0 )THEN
                     IY = 1
                  ELSE
                     IY = 1 - ( N - 1 )*INCY
                  END IF
                  IF( INCX.GT.0 )THEN
                     DO 20, IX = 1, 1 + ( N - 1 )*INCX, INCX
                        TEMP    = -X( IX )
                        X( IX ) =  Y( IY )
                        Y( IY ) =  TEMP
                        IY      =  IY       + INCY
   20                CONTINUE
                  ELSE
                     IX = 1 - ( N - 1 )*INCX
                     DO 30, I = 1, N
                        TEMP    = -X( IX )
                        X( IX ) =  Y( IY )
                        Y( IY ) =  TEMP
                        IX      =  IX       + INCX
                        IY      =  IY       + INCY
   30                CONTINUE
                  END IF
               END IF
            ELSE IF( ( C.EQ.ZERO ).AND.( S.EQ.( -ONE ) ) )THEN
               IF( ( INCX.EQ.INCY ).AND.( INCX.GT.0 ) )THEN
                  DO 40, IX = 1, 1 + ( N - 1 )*INCX, INCX
                     TEMP    =  X( IX )
                     X( IX ) = -Y( IX )
                     Y( IX ) =  TEMP
   40             CONTINUE
               ELSE
                  IF( INCY.GE.0 )THEN
                     IY = 1
                  ELSE
                     IY = 1 - ( N - 1 )*INCY
                  END IF
                  IF( INCX.GT.0 )THEN
                     DO 50, IX = 1, 1 + ( N - 1 )*INCX, INCX
                        TEMP    =  X( IX )
                        X( IX ) = -Y( IY )
                        Y( IY ) =  TEMP
                        IY      =  IY       + INCY
   50                CONTINUE
                  ELSE
                     IX = 1 - ( N - 1 )*INCX
                     DO 60, I = 1, N
                        TEMP    =  X( IX )
                        X( IX ) = -Y( IY )
                        Y( IY ) =  TEMP
                        IX      =  IX       + INCX
                        IY      =  IY       + INCY
   60                CONTINUE
                  END IF
               END IF
            ELSE
               IF( ( INCX.EQ.INCY ).AND.( INCX.GT.0 ) )THEN
                  DO 70, IX = 1, 1 + ( N - 1 )*INCX, INCX
                     TEMP    = X( IX )
                     X( IX ) = S*Y( IX ) + C*TEMP
                     Y( IX ) = C*Y( IX ) - S*TEMP
   70             CONTINUE
               ELSE
                  IF( INCY.GE.0 )THEN
                     IY = 1
                  ELSE
                     IY = 1 - ( N - 1 )*INCY
                  END IF
                  IF( INCX.GT.0 )THEN
                     DO 80, IX = 1, 1 + ( N - 1 )*INCX, INCX
                        TEMP    = X( IX )
                        X( IX ) = S*Y( IY ) + C*TEMP
                        Y( IY ) = C*Y( IY ) - S*TEMP
                        IY      = IY        + INCY
   80                CONTINUE
                  ELSE
                     IX = 1 - ( N - 1 )*INCX
                     DO 90, I = 1, N
                        TEMP    = X( IX )
                        X( IX ) = S*Y( IY ) + C*TEMP
                        Y( IY ) = C*Y( IY ) - S*TEMP
                        IX      = IX        + INCX
                        IY      = IY        + INCY
   90                CONTINUE
                  END IF
               END IF
            END IF
         END IF
      END IF
C
      RETURN
C
C     End of F06KPF. ( CSROT )
C
      END
*
      SUBROUTINE F06SAF( TRANS, M, N, ALPHA, A, LDA, X, INCX,
     $                   BETA, Y, INCY )
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     DOT VERSION FOR VECTOR MACHINES
C     .. Entry Points ..
      ENTRY      ZGEMV ( TRANS, M, N, ALPHA, A, LDA, X, INCX,
     $                   BETA, Y, INCY )
C     .. Scalar Arguments ..
      COMPLEX*16         ALPHA, BETA
      INTEGER            INCX, INCY, LDA, M, N
      CHARACTER*1        TRANS
C     .. Array Arguments ..
      COMPLEX*16         A( LDA, * ), X( * ), Y( * )
C     ..
C
C  Purpose
C  =======
C
C  ZGEMV  performs one of the matrix-vector operations
C
C     y := alpha*A*x + beta*y,   or   y := alpha*A'*x + beta*y,   or
C
C     y := alpha*conjg( A' )*x + beta*y,
C
C  where alpha and beta are scalars, x and y are vectors and A is an
C  m by n matrix.
C
C  Parameters
C  ==========
C
C  TRANS  - CHARACTER*1.
C           On entry, TRANS specifies the operation to be performed as
C           follows:
C
C              TRANS = 'N' or 'n'   y := alpha*A*x + beta*y.
C
C              TRANS = 'T' or 't'   y := alpha*A'*x + beta*y.
C
C              TRANS = 'C' or 'c'   y := alpha*conjg( A' )*x + beta*y.
C
C           Unchanged on exit.
C
C  M      - INTEGER.
C           On entry, M specifies the number of rows of the matrix A.
C           M must be at least zero.
C           Unchanged on exit.
C
C  N      - INTEGER.
C           On entry, N specifies the number of columns of the matrix A.
C           N must be at least zero.
C           Unchanged on exit.
C
C  ALPHA  - COMPLEX*16      .
C           On entry, ALPHA specifies the scalar alpha.
C           Unchanged on exit.
C
C  A      - COMPLEX*16       array of DIMENSION ( LDA, n ).
C           Before entry, the leading m by n part of the array A must
C           contain the matrix of coefficients.
C           Unchanged on exit.
C
C  LDA    - INTEGER.
C           On entry, LDA specifies the first dimension of A as declared
C           in the calling (sub) program. LDA must be at least
C           max( 1, m ).
C           Unchanged on exit.
C
C  X      - COMPLEX*16       array of DIMENSION at least
C           ( 1 + ( n - 1 )*abs( INCX ) ) when TRANS = 'N' or 'n'
C           and at least
C           ( 1 + ( m - 1 )*abs( INCX ) ) otherwise.
C           Before entry, the incremented array X must contain the
C           vector x.
C           Unchanged on exit.
C
C  INCX   - INTEGER.
C           On entry, INCX specifies the increment for the elements of
C           X. INCX must not be zero.
C           Unchanged on exit.
C
C  BETA   - COMPLEX*16      .
C           On entry, BETA specifies the scalar beta. When BETA is
C           supplied as zero then Y need not be set on input.
C           Unchanged on exit.
C
C  Y      - COMPLEX*16       array of DIMENSION at least
C           ( 1 + ( m - 1 )*abs( INCY ) ) when TRANS = 'N' or 'n'
C           and at least
C           ( 1 + ( n - 1 )*abs( INCY ) ) otherwise.
C           Before entry with BETA non-zero, the incremented array Y
C           must contain the vector y. On exit, Y is overwritten by the
C           updated vector y.
C
C  INCY   - INTEGER.
C           On entry, INCY specifies the increment for the elements of
C           Y. INCY must not be zero.
C           Unchanged on exit.
C
C
C  Level 2 Blas routine.
C
C  -- Written on 22-October-1986.
C     Jack Dongarra, Argonne National Lab.
C     Jeremy Du Croz, Nag Central Office.
C     Sven Hammarling, Nag Central Office.
C     Richard Hanson, Sandia National Labs.
C
C
C     .. Parameters ..
      COMPLEX*16         ONE
      PARAMETER        ( ONE  = ( 1.0D+0, 0.0D+0 ) )
      COMPLEX*16         ZERO
      PARAMETER        ( ZERO = ( 0.0D+0, 0.0D+0 ) )
C     .. Local Scalars ..
      COMPLEX*16         TEMP
      INTEGER            I, INFO, IX, IY, J, JY, KX, KY, LENX, LENY
      LOGICAL            NOCONJ
C     .. External Subroutines ..
      EXTERNAL           F06AAZ
C     .. Intrinsic Functions ..
      INTRINSIC          DCONJG, MAX
C     ..
C     .. Executable Statements ..
C
C     Test the input parameters.
C
      INFO = 0
      IF     ( .NOT.(TRANS.EQ.'N' .OR. TRANS.EQ.'n').AND.
     $         .NOT.(TRANS.EQ.'T' .OR. TRANS.EQ.'t').AND.
     $         .NOT.(TRANS.EQ.'C' .OR. TRANS.EQ.'c')      )THEN
         INFO = 1
      ELSE IF( M.LT.0 )THEN
         INFO = 2
      ELSE IF( N.LT.0 )THEN
         INFO = 3
      ELSE IF( LDA.LT.MAX( 1, M ) )THEN
         INFO = 6
      ELSE IF( INCX.EQ.0 )THEN
         INFO = 8
      ELSE IF( INCY.EQ.0 )THEN
         INFO = 11
      END IF
      IF( INFO.NE.0 )THEN
         CALL F06AAZ( 'F06SAF/ZGEMV ', INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF( ( M.EQ.0 ).OR.( N.EQ.0 ).OR.
     $    ( ( ALPHA.EQ.ZERO ).AND.( BETA.EQ.ONE ) ) )
     $   RETURN
C
      NOCONJ = (TRANS.EQ.'T' .OR. TRANS.EQ.'t')
C
C     Set  LENX  and  LENY, the lengths of the vectors x and y, and set
C     up the start points in  X  and  Y.
C
      IF( (TRANS.EQ.'N' .OR. TRANS.EQ.'n') )THEN
         LENX = N
         LENY = M
      ELSE
         LENX = M
         LENY = N
      END IF
      IF( INCX.GT.0 )THEN
         KX = 1
      ELSE
         KX = 1 - ( LENX - 1 )*INCX
      END IF
      IF( INCY.GT.0 )THEN
         KY = 1
      ELSE
         KY = 1 - ( LENY - 1 )*INCY
      END IF
C
C     Start the operations. In this version the inner loops are all
C     equivalent to dot-product operations.
C
C     First form  y := beta*y.
C
      IF( BETA.NE.ONE )THEN
         IF( INCY.EQ.1 )THEN
            IF( BETA.EQ.ZERO )THEN
               DO 10, I = 1, LENY
                  Y( I ) = ZERO
   10          CONTINUE
            ELSE
               DO 20, I = 1, LENY
                  Y( I ) = BETA*Y( I )
   20          CONTINUE
            END IF
         ELSE
            IY = KY
            IF( BETA.EQ.ZERO )THEN
               DO 30, I = 1, LENY
                  Y( IY ) = ZERO
                  IY      = IY   + INCY
   30          CONTINUE
            ELSE
               DO 40, I = 1, LENY
                  Y( IY ) = BETA*Y( IY )
                  IY      = IY           + INCY
   40          CONTINUE
            END IF
         END IF
      END IF
      IF( ALPHA.EQ.ZERO )
     $   RETURN
      JY = KY
      IF( (TRANS.EQ.'N' .OR. TRANS.EQ.'n') )THEN
C
C        Form  y := alpha*A*x + y.
C
         IF( INCX.EQ.1 )THEN
            DO 60, J = 1, M
               TEMP = ZERO
               DO 50, I = 1, N
                  TEMP = TEMP + A( J, I )*X( I )
   50          CONTINUE
               Y( JY ) = Y( JY ) + ALPHA*TEMP
               JY      = JY      + INCY
   60       CONTINUE
         ELSE
            DO 80, J = 1, M
               TEMP = ZERO
               IX   = KX
               DO 70, I = 1, N
                  TEMP = TEMP + A( J, I )*X( IX )
                  IX   = IX   + INCX
   70          CONTINUE
               Y( JY ) = Y( JY ) + ALPHA*TEMP
               JY      = JY      + INCY
   80       CONTINUE
         END IF
      ELSE
C
C        Form  y := alpha*A'*x + y  or  y := alpha*conjg( A' )*x + y.
C
         IF( INCX.EQ.1 )THEN
            DO 110, J = 1, N
               TEMP = ZERO
               IF( NOCONJ )THEN
                  DO 90, I = 1, M
                     TEMP = TEMP + A( I, J )*X( I )
   90             CONTINUE
               ELSE
                  DO 100, I = 1, M
                     TEMP = TEMP + DCONJG( A( I, J ) )*X( I )
  100             CONTINUE
               END IF
               Y( JY ) = Y( JY ) + ALPHA*TEMP
               JY      = JY      + INCY
  110       CONTINUE
         ELSE
            DO 140, J = 1, N
               TEMP = ZERO
               IX   = KX
               IF( NOCONJ )THEN
                  DO 120, I = 1, M
                     TEMP = TEMP + A( I, J )*X( IX )
                     IX   = IX   + INCX
  120             CONTINUE
               ELSE
                  DO 130, I = 1, M
                     TEMP = TEMP + DCONJG( A( I, J ) )*X( IX )
                     IX   = IX   + INCX
  130             CONTINUE
               END IF
               Y( JY ) = Y( JY ) + ALPHA*TEMP
               JY      = JY      + INCY
  140       CONTINUE
         END IF
      END IF
C
      RETURN
C
C     End of F06SAF (ZGEMV ).
C
      END
*
      SUBROUTINE F06SCF( UPLO, N, ALPHA, A, LDA, X, INCX,
     $                   BETA, Y, INCY )
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     DOT VERSION FOR VECTOR MACHINES
C     .. Entry Points ..
      ENTRY      ZHEMV ( UPLO, N, ALPHA, A, LDA, X, INCX,
     $                   BETA, Y, INCY )
C     .. Scalar Arguments ..
      COMPLEX*16         ALPHA, BETA
      INTEGER            INCX, INCY, LDA, N
      CHARACTER*1        UPLO
C     .. Array Arguments ..
      COMPLEX*16         A( LDA, * ), X( * ), Y( * )
C     ..
C
C  Purpose
C  =======
C
C  ZHEMV  performs the matrix-vector  operation
C
C     y := alpha*A*x + beta*y,
C
C  where alpha and beta are scalars, x and y are n element vectors and
C  A is an n by n hermitian matrix.
C
C  Parameters
C  ==========
C
C  UPLO   - CHARACTER*1.
C           On entry, UPLO specifies whether the upper or lower
C           triangular part of the array A is to be referenced as
C           follows:
C
C              UPLO = 'U' or 'u'   Only the upper triangular part of A
C                                  is to be referenced.
C
C              UPLO = 'L' or 'l'   Only the lower triangular part of A
C                                  is to be referenced.
C
C           Unchanged on exit.
C
C  N      - INTEGER.
C           On entry, N specifies the order of the matrix A.
C           N must be at least zero.
C           Unchanged on exit.
C
C  ALPHA  - COMPLEX*16      .
C           On entry, ALPHA specifies the scalar alpha.
C           Unchanged on exit.
C
C  A      - COMPLEX*16       array of DIMENSION ( LDA, n ).
C           Before entry with  UPLO = 'U' or 'u', the leading n by n
C           upper triangular part of the array A must contain the upper
C           triangular part of the hermitian matrix and the strictly
C           lower triangular part of A is not referenced.
C           Before entry with UPLO = 'L' or 'l', the leading n by n
C           lower triangular part of the array A must contain the lower
C           triangular part of the hermitian matrix and the strictly
C           upper triangular part of A is not referenced.
C           Note that the imaginary parts of the diagonal elements need
C           not be set and are assumed to be zero.
C           Unchanged on exit.
C
C  LDA    - INTEGER.
C           On entry, LDA specifies the first dimension of A as declared
C           in the calling (sub) program. LDA must be at least
C           max( 1, n ).
C           Unchanged on exit.
C
C  X      - COMPLEX*16       array of dimension at least
C           ( 1 + ( n - 1 )*abs( INCX ) ).
C           Before entry, the incremented array X must contain the n
C           element vector x.
C           Unchanged on exit.
C
C  INCX   - INTEGER.
C           On entry, INCX specifies the increment for the elements of
C           X. INCX must not be zero.
C           Unchanged on exit.
C
C  BETA   - COMPLEX*16      .
C           On entry, BETA specifies the scalar beta. When BETA is
C           supplied as zero then Y need not be set on input.
C           Unchanged on exit.
C
C  Y      - COMPLEX*16       array of dimension at least
C           ( 1 + ( n - 1 )*abs( INCY ) ).
C           Before entry, the incremented array Y must contain the n
C           element vector y. On exit, Y is overwritten by the updated
C           vector y.
C
C  INCY   - INTEGER.
C           On entry, INCY specifies the increment for the elements of
C           Y. INCY must not be zero.
C           Unchanged on exit.
C
C
C  Level 2 Blas routine.
C
C  -- Written on 22-October-1986.
C     Jack Dongarra, Argonne National Lab.
C     Jeremy Du Croz, Nag Central Office.
C     Sven Hammarling, Nag Central Office.
C     Richard Hanson, Sandia National Labs.
C
C
C     .. Parameters ..
      COMPLEX*16         ONE
      PARAMETER        ( ONE  = ( 1.0D+0, 0.0D+0 ) )
      COMPLEX*16         ZERO
      PARAMETER        ( ZERO = ( 0.0D+0, 0.0D+0 ) )
C     .. Local Scalars ..
      COMPLEX*16         TEMP
      INTEGER            I, INFO, IX, IY, J, JX, JY, KX, KY
C     .. External Subroutines ..
      EXTERNAL           F06AAZ
C     .. Intrinsic Functions ..
      INTRINSIC          DCONJG, MAX, DBLE
C     ..
C     .. Executable Statements ..
C
C     Test the input parameters.
C
      INFO = 0
      IF     ( .NOT.(UPLO.EQ.'U' .OR. UPLO.EQ.'u').AND.
     $         .NOT.(UPLO.EQ.'L' .OR. UPLO.EQ.'l')      )THEN
         INFO = 1
      ELSE IF( N.LT.0 )THEN
         INFO = 2
      ELSE IF( LDA.LT.MAX( 1, N ) )THEN
         INFO = 5
      ELSE IF( INCX.EQ.0 )THEN
         INFO = 7
      ELSE IF( INCY.EQ.0 )THEN
         INFO = 10
      END IF
      IF( INFO.NE.0 )THEN
         CALL F06AAZ( 'F06SCF/ZHEMV ', INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF( ( N.EQ.0 ).OR.( ( ALPHA.EQ.ZERO ).AND.( BETA.EQ.ONE ) ) )
     $   RETURN
C
C     Set up the start points in  X  and  Y.
C
      IF( INCX.GT.0 )THEN
         KX = 1
      ELSE
         KX = 1 - ( N - 1 )*INCX
      END IF
      IF( INCY.GT.0 )THEN
         KY = 1
      ELSE
         KY = 1 - ( N - 1 )*INCY
      END IF
C
C     Start the operations. In this version the inner loops are all
C     equivalent to dot-product operations.
C
C     First form  y := beta*y.
C
      IF( BETA.NE.ONE )THEN
         IF( INCY.EQ.1 )THEN
            IF( BETA.EQ.ZERO )THEN
               DO 10, I = 1, N
                  Y( I ) = ZERO
   10          CONTINUE
            ELSE
               DO 20, I = 1, N
                  Y( I ) = BETA*Y( I )
   20          CONTINUE
            END IF
         ELSE
            IY = KY
            IF( BETA.EQ.ZERO )THEN
               DO 30, I = 1, N
                  Y( IY ) = ZERO
                  IY      = IY   + INCY
   30          CONTINUE
            ELSE
               DO 40, I = 1, N
                  Y( IY ) = BETA*Y( IY )
                  IY      = IY           + INCY
   40          CONTINUE
            END IF
         END IF
      END IF
      IF( ALPHA.EQ.ZERO )
     $   RETURN
      JY = KY
      IF( (UPLO.EQ.'U' .OR. UPLO.EQ.'u') )THEN
C
C        Form  y  when A is stored in upper triangle.
C
         IF( INCX.EQ.1 )THEN
            DO 60, J = 1, N
               TEMP = ZERO
               DO 50, I = 1, J - 1
                  TEMP = TEMP + DCONJG( A( I, J ) )*X( I )
   50          CONTINUE
               TEMP = TEMP + DBLE( A( J, J ) )*X( J )
               DO 55, I = J + 1, N
                  TEMP = TEMP + A( J, I )*X( I )
   55          CONTINUE
               Y(JY) = Y( JY ) + ALPHA*TEMP
               JY    = JY      + INCY
   60       CONTINUE
         ELSE
            JX = KX
            DO 80, J = 1, N
               TEMP = ZERO
               IX   = KX
               DO 70, I = 1, J - 1
                  TEMP = TEMP + DCONJG( A( I, J ) )*X( IX )
                  IX   = IX   + INCX
   70          CONTINUE
               TEMP = TEMP + DBLE( A( J, J ) )*X( JX )
               IX   = JX   + INCX
               DO 75, I = J + 1, N
                  TEMP = TEMP + A( J, I )*X( IX )
                  IX   = IX   + INCX
   75          CONTINUE
               Y( JY ) = Y( JY ) + ALPHA*TEMP
               JX      = JX      + INCX
               JY      = JY      + INCY
   80       CONTINUE
         END IF
      ELSE
C
C        Form  y  when A is stored in lower triangle.
C
         IF( INCX.EQ.1 )THEN
            DO 100, J = 1, N
               TEMP = ZERO
               DO 85, I = 1, J - 1
                  TEMP = TEMP + A( J, I )*X( I )
   85          CONTINUE
               TEMP = TEMP + DBLE( A( J, J ) )*X( J )
               DO 90, I = J + 1, N
                  TEMP = TEMP + DCONJG( A( I, J ) )*X( I )
   90          CONTINUE
               Y( JY ) = Y( JY ) + ALPHA*TEMP
               JY      = JY      + INCY
  100       CONTINUE
         ELSE
            JX = KX
            DO 120, J = 1, N
               TEMP = ZERO
               IX   = KX
               DO 105, I = 1, J - 1
                  TEMP = TEMP + A( J, I )*X( IX )
                  IX   = IX   + INCX
  105          CONTINUE
               TEMP = TEMP + DBLE( A( J, J ) )*X( JX )
               IX   = JX   + INCX
               DO 110, I = J + 1, N
                  TEMP = TEMP + DCONJG( A( I, J ) )*X( IX )
                  IX   = IX   + INCX
  110          CONTINUE
               Y( JY ) = Y( JY ) + ALPHA*TEMP
               JX      = JX      + INCX
               JY      = JY      + INCY
  120       CONTINUE
         END IF
      END IF
C
      RETURN
C
C     End of F06SCF (ZHEMV ).
C
      END
*
      SUBROUTINE F06SRF( UPLO, N, ALPHA, X, INCX, Y, INCY, A, LDA )
C     MARK 13 RE-ISSUE. NAG COPYRIGHT 1988.
C     .. Entry Points ..
      ENTRY      ZHER2 ( UPLO, N, ALPHA, X, INCX, Y, INCY, A, LDA )
C     .. Scalar Arguments ..
      COMPLEX*16         ALPHA
      INTEGER            INCX, INCY, LDA, N
      CHARACTER*1        UPLO
C     .. Array Arguments ..
      COMPLEX*16         A( LDA, * ), X( * ), Y( * )
C     ..
C
C  Purpose
C  =======
C
C  ZHER2  performs the hermitian rank 2 operation
C
C     A := alpha*x*conjg( y' ) + conjg( alpha )*y*conjg( x' ) + A,
C
C  where alpha is a scalar, x and y are n element vectors and A is an n
C  by n hermitian matrix.
C
C  Parameters
C  ==========
C
C  UPLO   - CHARACTER*1.
C           On entry, UPLO specifies whether the upper or lower
C           triangular part of the array A is to be referenced as
C           follows:
C
C              UPLO = 'U' or 'u'   Only the upper triangular part of A
C                                  is to be referenced.
C
C              UPLO = 'L' or 'l'   Only the lower triangular part of A
C                                  is to be referenced.
C
C           Unchanged on exit.
C
C  N      - INTEGER.
C           On entry, N specifies the order of the matrix A.
C           N must be at least zero.
C           Unchanged on exit.
C
C  ALPHA  - COMPLEX*16      .
C           On entry, ALPHA specifies the scalar alpha.
C           Unchanged on exit.
C
C  X      - COMPLEX*16       array of dimension at least
C           ( 1 + ( n - 1 )*abs( INCX ) ).
C           Before entry, the incremented array X must contain the n
C           element vector x.
C           Unchanged on exit.
C
C  INCX   - INTEGER.
C           On entry, INCX specifies the increment for the elements of
C           X. INCX must not be zero.
C           Unchanged on exit.
C
C  Y      - COMPLEX*16       array of dimension at least
C           ( 1 + ( n - 1 )*abs( INCY ) ).
C           Before entry, the incremented array Y must contain the n
C           element vector y.
C           Unchanged on exit.
C
C  INCY   - INTEGER.
C           On entry, INCY specifies the increment for the elements of
C           Y. INCY must not be zero.
C           Unchanged on exit.
C
C  A      - COMPLEX*16       array of DIMENSION ( LDA, n ).
C           Before entry with  UPLO = 'U' or 'u', the leading n by n
C           upper triangular part of the array A must contain the upper
C           triangular part of the hermitian matrix and the strictly
C           lower triangular part of A is not referenced. On exit, the
C           upper triangular part of the array A is overwritten by the
C           upper triangular part of the updated matrix.
C           Before entry with UPLO = 'L' or 'l', the leading n by n
C           lower triangular part of the array A must contain the lower
C           triangular part of the hermitian matrix and the strictly
C           upper triangular part of A is not referenced. On exit, the
C           lower triangular part of the array A is overwritten by the
C           lower triangular part of the updated matrix.
C           Note that the imaginary parts of the diagonal elements need
C           not be set, they are assumed to be zero, and on exit they
C           are set to zero.
C
C  LDA    - INTEGER.
C           On entry, LDA specifies the first dimension of A as declared
C           in the calling (sub) program. LDA must be at least
C           max( 1, n ).
C           Unchanged on exit.
C
C
C  Level 2 Blas routine.
C
C  -- Written on 22-October-1986.
C     Jack Dongarra, Argonne National Lab.
C     Jeremy Du Croz, Nag Central Office.
C     Sven Hammarling, Nag Central Office.
C     Richard Hanson, Sandia National Labs.
C
C
C     .. Parameters ..
      COMPLEX*16         ZERO
      PARAMETER        ( ZERO = ( 0.0D+0, 0.0D+0 ) )
C     .. Local Scalars ..
      COMPLEX*16         TEMP1, TEMP2
      INTEGER            I, INFO, IX, IY, J, JX, JY, KX, KY
C     .. External Subroutines ..
      EXTERNAL           F06AAZ
C     .. Intrinsic Functions ..
      INTRINSIC          DCONJG, MAX, DBLE
C     ..
C     .. Executable Statements ..
C
C     Test the input parameters.
C
      INFO = 0
      IF     ( .NOT.(UPLO.EQ.'U' .OR. UPLO.EQ.'u').AND.
     $         .NOT.(UPLO.EQ.'L' .OR. UPLO.EQ.'l')      )THEN
         INFO = 1
      ELSE IF( N.LT.0 )THEN
         INFO = 2
      ELSE IF( INCX.EQ.0 )THEN
         INFO = 5
      ELSE IF( INCY.EQ.0 )THEN
         INFO = 7
      ELSE IF( LDA.LT.MAX( 1, N ) )THEN
         INFO = 9
      END IF
      IF( INFO.NE.0 )THEN
         CALL F06AAZ( 'F06SRF/ZHER2 ', INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF( ( N.EQ.0 ).OR.( ALPHA.EQ.ZERO ) )
     $   RETURN
C
C     Set up the start points in X and Y if the increments are not both
C     unity.
C
      IF( ( INCX.NE.1 ).OR.( INCY.NE.1 ) )THEN
         IF( INCX.GT.0 )THEN
            KX = 1
         ELSE
            KX = 1 - ( N - 1 )*INCX
         END IF
         IF( INCY.GT.0 )THEN
            KY = 1
         ELSE
            KY = 1 - ( N - 1 )*INCY
         END IF
         JX = KX
         JY = KY
      END IF
C
C     Start the operations. In this version the elements of A are
C     accessed sequentially with one pass through the triangular part
C     of A.
C
      IF( (UPLO.EQ.'U' .OR. UPLO.EQ.'u') )THEN
C
C        Form  A  when A is stored in the upper triangle.
C
         IF( ( INCX.EQ.1 ).AND.( INCY.EQ.1 ) )THEN
            DO 20, J = 1, N
               IF( ( X( J ).NE.ZERO ).OR.( Y( J ).NE.ZERO ) )THEN
                  TEMP1 = ALPHA*DCONJG( Y( J ) )
                  TEMP2 = DCONJG( ALPHA*X( J ) )
                  DO 10, I = 1, J - 1
                     A( I, J ) = A( I, J ) + X( I )*TEMP1 + Y( I )*TEMP2
   10             CONTINUE
                  A( J, J ) = DBLE( A( J, J ) ) +
     $                        DBLE( X( J )*TEMP1 + Y( J )*TEMP2 )
               ELSE
                  A( J, J ) = DBLE( A( J, J ) )
               END IF
   20       CONTINUE
         ELSE
            DO 40, J = 1, N
               IF( ( X( JX ).NE.ZERO ).OR.( Y( JY ).NE.ZERO ) )THEN
                  TEMP1 = ALPHA*DCONJG( Y( JY ) )
                  TEMP2 = DCONJG( ALPHA*X( JX ) )
                  IX    = KX
                  IY    = KY
                  DO 30, I = 1, J - 1
                     A( I, J ) = A( I, J ) + X( IX )*TEMP1
     $                                     + Y( IY )*TEMP2
                     IX        = IX        + INCX
                     IY        = IY        + INCY
   30             CONTINUE
                  A( J, J ) = DBLE( A( J, J ) ) +
     $                        DBLE( X( JX )*TEMP1 + Y( JY )*TEMP2 )
               ELSE
                  A( J, J ) = DBLE( A( J, J ) )
               END IF
               JX = JX + INCX
               JY = JY + INCY
   40       CONTINUE
         END IF
      ELSE
C
C        Form  A  when A is stored in the lower triangle.
C
         IF( ( INCX.EQ.1 ).AND.( INCY.EQ.1 ) )THEN
            DO 60, J = 1, N
               IF( ( X( J ).NE.ZERO ).OR.( Y( J ).NE.ZERO ) )THEN
                  TEMP1     = ALPHA*DCONJG( Y( J ) )
                  TEMP2     = DCONJG( ALPHA*X( J ) )
                  A( J, J ) = DBLE( A( J, J ) ) +
     $                        DBLE( X( J )*TEMP1 + Y( J )*TEMP2 )
                  DO 50, I = J + 1, N
                     A( I, J ) = A( I, J ) + X( I )*TEMP1 + Y( I )*TEMP2
   50             CONTINUE
               ELSE
                  A( J, J ) = DBLE( A( J, J ) )
               END IF
   60       CONTINUE
         ELSE
            DO 80, J = 1, N
               IF( ( X( JX ).NE.ZERO ).OR.( Y( JY ).NE.ZERO ) )THEN
                  TEMP1     = ALPHA*DCONJG( Y( JY ) )
                  TEMP2     = DCONJG( ALPHA*X( JX ) )
                  A( J, J ) = DBLE( A( J, J ) ) +
     $                        DBLE( X( JX )*TEMP1 + Y( JY )*TEMP2 )
                  IX        = JX
                  IY        = JY
                  DO 70, I = J + 1, N
                     IX        = IX        + INCX
                     IY        = IY        + INCY
                     A( I, J ) = A( I, J ) + X( IX )*TEMP1
     $                                     + Y( IY )*TEMP2
   70             CONTINUE
               ELSE
                  A( J, J ) = DBLE( A( J, J ) )
               END IF
               JX = JX + INCX
               JY = JY + INCY
   80       CONTINUE
         END IF
      END IF
C
      RETURN
C
C     End of F06SRF (ZHER2 ).
C
      END
*
      SUBROUTINE F06SSF( UPLO, N, ALPHA, X, INCX, Y, INCY, AP )
C     MARK 13 RE-ISSUE. NAG COPYRIGHT 1988.
C     .. Entry Points ..
      ENTRY      ZHPR2 ( UPLO, N, ALPHA, X, INCX, Y, INCY, AP )
C     .. Scalar Arguments ..
      COMPLEX*16         ALPHA
      INTEGER            INCX, INCY, N
      CHARACTER*1        UPLO
C     .. Array Arguments ..
      COMPLEX*16         AP( * ), X( * ), Y( * )
C     ..
C
C  Purpose
C  =======
C
C  ZHPR2  performs the hermitian rank 2 operation
C
C     A := alpha*x*conjg( y' ) + conjg( alpha )*y*conjg( x' ) + A,
C
C  where alpha is a scalar, x and y are n element vectors and A is an
C  n by n hermitian matrix.
C
C  Parameters
C  ==========
C
C  UPLO   - CHARACTER*1.
C           On entry, UPLO specifies whether the upper or lower
C           triangular part of the matrix A is supplied in the packed
C           array AP as follows:
C
C              UPLO = 'U' or 'u'   The upper triangular part of A is
C                                  supplied in AP.
C
C              UPLO = 'L' or 'l'   The lower triangular part of A is
C                                  supplied in AP.
C
C           Unchanged on exit.
C
C  N      - INTEGER.
C           On entry, N specifies the order of the matrix A.
C           N must be at least zero.
C           Unchanged on exit.
C
C  ALPHA  - COMPLEX*16      .
C           On entry, ALPHA specifies the scalar alpha.
C           Unchanged on exit.
C
C  X      - COMPLEX*16       array of dimension at least
C           ( 1 + ( n - 1 )*abs( INCX ) ).
C           Before entry, the incremented array X must contain the n
C           element vector x.
C           Unchanged on exit.
C
C  INCX   - INTEGER.
C           On entry, INCX specifies the increment for the elements of
C           X. INCX must not be zero.
C           Unchanged on exit.
C
C  Y      - COMPLEX*16       array of dimension at least
C           ( 1 + ( n - 1 )*abs( INCY ) ).
C           Before entry, the incremented array Y must contain the n
C           element vector y.
C           Unchanged on exit.
C
C  INCY   - INTEGER.
C           On entry, INCY specifies the increment for the elements of
C           Y. INCY must not be zero.
C           Unchanged on exit.
C
C  AP     - COMPLEX*16       array of DIMENSION at least
C           ( ( n*( n + 1 ) )/2 ).
C           Before entry with  UPLO = 'U' or 'u', the array AP must
C           contain the upper triangular part of the hermitian matrix
C           packed sequentially, column by column, so that AP( 1 )
C           contains a( 1, 1 ), AP( 2 ) and AP( 3 ) contain a( 1, 2 )
C           and a( 2, 2 ) respectively, and so on. On exit, the array
C           AP is overwritten by the upper triangular part of the
C           updated matrix.
C           Before entry with UPLO = 'L' or 'l', the array AP must
C           contain the lower triangular part of the hermitian matrix
C           packed sequentially, column by column, so that AP( 1 )
C           contains a( 1, 1 ), AP( 2 ) and AP( 3 ) contain a( 2, 1 )
C           and a( 3, 1 ) respectively, and so on. On exit, the array
C           AP is overwritten by the lower triangular part of the
C           updated matrix.
C           Note that the imaginary parts of the diagonal elements need
C           not be set, they are assumed to be zero, and on exit they
C           are set to zero.
C
C
C  Level 2 Blas routine.
C
C  -- Written on 22-October-1986.
C     Jack Dongarra, Argonne National Lab.
C     Jeremy Du Croz, Nag Central Office.
C     Sven Hammarling, Nag Central Office.
C     Richard Hanson, Sandia National Labs.
C
C
C     .. Parameters ..
      COMPLEX*16         ZERO
      PARAMETER        ( ZERO = ( 0.0D+0, 0.0D+0 ) )
C     .. Local Scalars ..
      COMPLEX*16         TEMP1, TEMP2
      INTEGER            I, INFO, IX, IY, J, JX, JY, K, KK, KX, KY
C     .. External Subroutines ..
      EXTERNAL           F06AAZ
C     .. Intrinsic Functions ..
      INTRINSIC          DCONJG, DBLE
C     ..
C     .. Executable Statements ..
C
C     Test the input parameters.
C
      INFO = 0
      IF     ( .NOT.(UPLO.EQ.'U' .OR. UPLO.EQ.'u').AND.
     $         .NOT.(UPLO.EQ.'L' .OR. UPLO.EQ.'l')      )THEN
         INFO = 1
      ELSE IF( N.LT.0 )THEN
         INFO = 2
      ELSE IF( INCX.EQ.0 )THEN
         INFO = 5
      ELSE IF( INCY.EQ.0 )THEN
         INFO = 7
      END IF
      IF( INFO.NE.0 )THEN
         CALL F06AAZ( 'F06SSF/ZHPR2 ', INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF( ( N.EQ.0 ).OR.( ALPHA.EQ.ZERO ) )
     $   RETURN
C
C     Set up the start points in X and Y if the increments are not both
C     unity.
C
      IF( ( INCX.NE.1 ).OR.( INCY.NE.1 ) )THEN
         IF( INCX.GT.0 )THEN
            KX = 1
         ELSE
            KX = 1 - ( N - 1 )*INCX
         END IF
         IF( INCY.GT.0 )THEN
            KY = 1
         ELSE
            KY = 1 - ( N - 1 )*INCY
         END IF
         JX = KX
         JY = KY
      END IF
C
C     Start the operations. In this version the elements of the array AP
C     are accessed sequentially with one pass through AP.
C
      KK = 1
      IF( (UPLO.EQ.'U' .OR. UPLO.EQ.'u') )THEN
C
C        Form  A  when upper triangle is stored in AP.
C
         IF( ( INCX.EQ.1 ).AND.( INCY.EQ.1 ) )THEN
            DO 20, J = 1, N
               IF( ( X( J ).NE.ZERO ).OR.( Y( J ).NE.ZERO ) )THEN
                  TEMP1 = ALPHA*DCONJG( Y( J ) )
                  TEMP2 = DCONJG( ALPHA*X( J ) )
                  K     = KK
                  DO 10, I = 1, J - 1
                     AP( K ) = AP( K ) + X( I )*TEMP1 + Y( I )*TEMP2
                     K       = K       + 1
   10             CONTINUE
                  AP( KK + J - 1 ) = DBLE( AP( KK + J - 1 ) ) +
     $                               DBLE( X( J )*TEMP1 + Y( J )*TEMP2 )
               ELSE
                  AP( KK + J - 1 ) = DBLE( AP( KK + J - 1 ) )
               END IF
               KK = KK + J
   20       CONTINUE
         ELSE
            DO 40, J = 1, N
               IF( ( X( JX ).NE.ZERO ).OR.( Y( JY ).NE.ZERO ) )THEN
                  TEMP1 = ALPHA*DCONJG( Y( JY ) )
                  TEMP2 = DCONJG( ALPHA*X( JX ) )
                  IX    = KX
                  IY    = KY
                  DO 30, K = KK, KK + J - 2
                     AP( K ) = AP( K ) + X( IX )*TEMP1 + Y( IY )*TEMP2
                     IX      = IX      + INCX
                     IY      = IY      + INCY
   30             CONTINUE
                  AP( KK + J - 1 ) = DBLE( AP( KK + J - 1 ) ) +
     $                               DBLE( X( JX )*TEMP1 +
     $                                     Y( JY )*TEMP2 )
               ELSE
                  AP( KK + J - 1 ) = DBLE( AP( KK + J - 1 ) )
               END IF
               JX = JX + INCX
               JY = JY + INCY
               KK = KK + J
   40       CONTINUE
         END IF
      ELSE
C
C        Form  A  when lower triangle is stored in AP.
C
         IF( ( INCX.EQ.1 ).AND.( INCY.EQ.1 ) )THEN
            DO 60, J = 1, N
               IF( ( X( J ).NE.ZERO ).OR.( Y( J ).NE.ZERO ) )THEN
                  TEMP1   = ALPHA*DCONJG( Y( J ) )
                  TEMP2   = DCONJG( ALPHA*X( J ) )
                  AP( KK ) = DBLE( AP( KK ) ) +
     $                       DBLE( X( J )*TEMP1 + Y( J )*TEMP2 )
                  K        = KK               + 1
                  DO 50, I = J + 1, N
                     AP( K ) = AP( K ) + X( I )*TEMP1 + Y( I )*TEMP2
                     K       = K       + 1
   50             CONTINUE
               ELSE
                  AP( KK ) = DBLE( AP( KK ) )
               END IF
               KK = KK + N - J + 1
   60       CONTINUE
         ELSE
            DO 80, J = 1, N
               IF( ( X( JX ).NE.ZERO ).OR.( Y( JY ).NE.ZERO ) )THEN
                  TEMP1    = ALPHA*DCONJG( Y( JY ) )
                  TEMP2    = DCONJG( ALPHA*X( JX ) )
                  AP( KK ) = DBLE( AP( KK ) ) +
     $                       DBLE( X( JX )*TEMP1 + Y( JY )*TEMP2 )
                  IX       = JX
                  IY       = JY
                  DO 70, K = KK + 1, KK + N - J
                     IX      = IX      + INCX
                     IY      = IY      + INCY
                     AP( K ) = AP( K ) + X( IX )*TEMP1 + Y( IY )*TEMP2
   70             CONTINUE
               ELSE
                  AP( KK ) = DBLE( AP( KK ) )
               END IF
               JX = JX + INCX
               JY = JY + INCY
               KK = KK + N - J + 1
   80       CONTINUE
         END IF
      END IF
C
      RETURN
C
C     End of F06SSF (ZHPR2 ).
C
      END
*
      SUBROUTINE F06THF( MATRIX, M, N, CONST, DIAG, A, LDA )
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C     .. Scalar Arguments ..
      CHARACTER*1        MATRIX
      COMPLEX*16         CONST, DIAG
      INTEGER            LDA, M, N
C     .. Array Arguments ..
      COMPLEX*16         A( LDA, * )
C     ..
C
C  F06THF forms the m by n matrix A given by
C
C     a( i, j ) = (  diag  i.eq.j,
C                 (
C                 ( const  i.ne.j.
C
C  If   MATRIX = 'G' or 'g'   then  A  is regarded  as a general matrix,
C  if   MATRIX = 'U' or 'u'   then  A  is regarded  as upper triangular,
C                             and only  elements  for which  i.le.j  are
C                             referenced,
C  if   MATRIX = 'L' or 'l'   then  A  is regarded  as lower triangular,
C                             and only  elements  for which  i.ge.j  are
C                             referenced.
C
C
C  Nag Fortran 77 O( n**2 ) basic linear algebra routine.
C
C  -- Written on 21-November-1986.
C     Sven Hammarling, Nag Central Office.
C
C
C     .. Local Scalars ..
      INTEGER            I, J
C     .. Intrinsic Functions ..
      INTRINSIC          MIN
C     ..
C     .. Executable Statements ..
      IF( ( MATRIX.EQ.'G' ).OR.( MATRIX.EQ.'g' ) )THEN
         DO 20 J = 1, N
            DO 10 I = 1, M
               A( I, J ) = CONST
   10       CONTINUE
   20    CONTINUE
         IF( CONST.NE.DIAG )THEN
            DO 30 I = 1, MIN( M, N )
               A( I, I ) = DIAG
   30       CONTINUE
         END IF
      ELSE IF( ( MATRIX.EQ.'U' ).OR.( MATRIX.EQ.'u' ) )THEN
         DO 50 J = 1, N
            DO 40 I = 1, MIN( M, J )
               A( I, J ) = CONST
   40       CONTINUE
   50    CONTINUE
         IF( CONST.NE.DIAG )THEN
            DO 60 I = 1, MIN( M, N )
               A( I, I ) = DIAG
   60       CONTINUE
         END IF
      ELSE IF( ( MATRIX.EQ.'L' ).OR.( MATRIX.EQ.'l' ) )THEN
         DO 80 J = 1, MIN( M, N )
            DO 70 I = J, M
               A( I, J ) = CONST
   70       CONTINUE
   80    CONTINUE
         IF( CONST.NE.DIAG )THEN
            DO 90 I = 1, MIN( M, N )
               A( I, I ) = DIAG
   90       CONTINUE
         END IF
      END IF
C
      RETURN
C
C     End of F06THF. ( CMLOAD )
C
      END
*
      SUBROUTINE F06VXF( SIDE, PIVOT, DIRECT, M, N, K1, K2, C, S, A,
     $                   LDA )
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C     .. Scalar Arguments ..
      INTEGER            K1, K2, LDA, M, N
      CHARACTER*1        DIRECT, PIVOT, SIDE
C     .. Array Arguments ..
      COMPLEX*16         A( LDA, * )
      DOUBLE PRECISION   C( * ), S( * )
C     ..
C
C  F06VXF  performs the transformation
C
C     A := P*A,    when   SIDE = 'L' or 'l'  (  Left-hand side )
C
C     A := A*P',   when   SIDE = 'R' or 'r'  ( Right-hand side )
C
C  where A is an m by n matrix and P is an orthogonal matrix, consisting
C  of  a  sequence  of  plane  rotations,  applied  in planes  k1 to k2,
C  determined by the parameters PIVOT and DIRECT as follows:
C
C     When  PIVOT  = 'V' or 'v'  ( Variable pivot )
C     and   DIRECT = 'F' or 'f'  ( Forward sequence ) then
C
C        P is given as a sequence of plane rotation matrices
C
C           P = P( k2 - 1 )*...*P( k1 + 1 )*P( k1 ),
C
C        where  P( k )  is a plane rotation matrix for the  ( k, k + 1 )
C        plane.
C
C     When  PIVOT  = 'V' or 'v'  ( Variable pivot )
C     and   DIRECT = 'B' or 'b'  ( Backward sequence ) then
C
C        P is given as a sequence of plane rotation matrices
C
C           P = P( k1 )*P( k1 + 1 )*...*P( k2 - 1 ),
C
C        where  P( k )  is a plane rotation matrix for the  ( k, k + 1 )
C        plane.
C
C     When  PIVOT  = 'T' or 't'  ( Top pivot )
C     and   DIRECT = 'F' or 'f'  ( Forward sequence ) then
C
C        P is given as a sequence of plane rotation matrices
C
C           P = P( k2 - 1 )*P( k2 - 2 )*...*P( k1 ),
C
C        where  P( k )  is a plane rotation matrix for the ( k1, k + 1 )
C        plane.
C
C     When  PIVOT  = 'T' or 't'  ( Top pivot )
C     and   DIRECT = 'B' or 'b'  ( Backward sequence ) then
C
C        P is given as a sequence of plane rotation matrices
C
C           P = P( k1 )*P( k1 + 1 )*...*P( k2 - 1 ),
C
C        where  P( k )  is a plane rotation matrix for the ( k1, k + 1 )
C        plane.
C
C     When  PIVOT  = 'B' or 'b'  ( Bottom pivot )
C     and   DIRECT = 'F' or 'f'  ( Forward sequence ) then
C
C        P is given as a sequence of plane rotation matrices
C
C           P = P( k2 - 1 )*P( k2 - 2 )*...*P( k1 ),
C
C        where  P( k )  is a  plane rotation  matrix  for the  ( k, k2 )
C        plane.
C
C     When  PIVOT  = 'B' or 'b'  ( Bottom pivot )
C     and   DIRECT = 'B' or 'b'  ( Backward sequence ) then
C
C        P is given as a sequence of plane rotation matrices
C
C           P = P( k1 )*P( k1 + 1 )*...*P( k2 - 1 ),
C
C        where  P( k )  is a  plane rotation  matrix  for the  ( k, k2 )
C        plane.
C
C  c( k ) and s( k )  must contain the  cosine and sine  that define the
C  matrix  P( k ).  The  two by two  plane rotation  part of the  matrix
C  P( k ), R( k ), is assumed to be of the form
C
C     R( k ) = (  c( k )  s( k ) ),   with  c( k ) and s( k )  real.
C              ( -s( k )  c( k ) )
C
C  If m, n or k1 are less than unity,  or k2 is not greater than k1,  or
C  SIDE = 'L' or 'l'  and  k2  is greater than  m, or  SIDE = 'R' or 'r'
C  and  k2  is greater than  n,  then an  immediate return  is effected.
C
C
C  Nag Fortran 77 O( n**2 ) basic linear algebra routine.
C
C  -- Written on 20-November-1986.
C     Sven Hammarling and Mick Pont, Nag Central Office.
C
C
C     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
C     .. Local Scalars ..
      COMPLEX*16         AIJ, TEMP
      DOUBLE PRECISION   CTEMP, STEMP
      INTEGER            I, J
      LOGICAL            LEFT, RIGHT
C     .. Intrinsic Functions ..
      INTRINSIC          MIN
C     ..
C     .. Executable Statements ..
      LEFT = ( SIDE.EQ.'L' ).OR.( SIDE.EQ.'l' )
      RIGHT = ( SIDE.EQ.'R' ).OR.( SIDE.EQ.'r' )
      IF( ( MIN( M, N, K1 ).LT.1 ).OR.( K2.LE.K1 ).OR.
     $    ( ( LEFT ).AND.( K2.GT.M ) ).OR.
     $    ( ( RIGHT ).AND.( K2.GT.N ) ) )RETURN
      IF( LEFT )THEN
         IF( ( PIVOT.EQ.'V' ).OR.( PIVOT.EQ.'v' ) )THEN
            IF( ( DIRECT.EQ.'F' ).OR.( DIRECT.EQ.'f' ) )THEN
               DO 20 J = 1, N
                  AIJ = A( K1, J )
                  DO 10 I = K1, K2 - 1
                     TEMP = A( I + 1, J )
                     A( I, J ) = S( I )*TEMP + C( I )*AIJ
                     AIJ = C( I )*TEMP - S( I )*AIJ
   10             CONTINUE
                  A( K2, J ) = AIJ
   20          CONTINUE
            ELSE IF( ( DIRECT.EQ.'B' ).OR.( DIRECT.EQ.'b' ) )THEN
               DO 40 J = 1, N
                  AIJ = A( K2, J )
                  DO 30 I = K2 - 1, K1, -1
                     TEMP = A( I, J )
                     A( I + 1, J ) = C( I )*AIJ - S( I )*TEMP
                     AIJ = S( I )*AIJ + C( I )*TEMP
   30             CONTINUE
                  A( K1, J ) = AIJ
   40          CONTINUE
            END IF
         ELSE IF( ( PIVOT.EQ.'T' ).OR.( PIVOT.EQ.'t' ) )THEN
            IF( ( DIRECT.EQ.'F' ).OR.( DIRECT.EQ.'f' ) )THEN
               DO 60 J = 1, N
                  TEMP = A( K1, J )
                  DO 50 I = K1, K2 - 1
                     AIJ = A( I + 1, J )
                     A( I + 1, J ) = C( I )*AIJ - S( I )*TEMP
                     TEMP = S( I )*AIJ + C( I )*TEMP
   50             CONTINUE
                  A( K1, J ) = TEMP
   60          CONTINUE
            ELSE IF( ( DIRECT.EQ.'B' ).OR.( DIRECT.EQ.'b' ) )THEN
               DO 80 J = 1, N
                  TEMP = A( K1, J )
                  DO 70 I = K2 - 1, K1, -1
                     AIJ = A( I + 1, J )
                     A( I + 1, J ) = C( I )*AIJ - S( I )*TEMP
                     TEMP = S( I )*AIJ + C( I )*TEMP
   70             CONTINUE
                  A( K1, J ) = TEMP
   80          CONTINUE
            END IF
         ELSE IF( ( PIVOT.EQ.'B' ).OR.( PIVOT.EQ.'b' ) )THEN
            IF( ( DIRECT.EQ.'F' ).OR.( DIRECT.EQ.'f' ) )THEN
               DO 100 J = 1, N
                  TEMP = A( K2, J )
                  DO 90 I = K1, K2 - 1
                     AIJ = A( I, J )
                     A( I, J ) = S( I )*TEMP + C( I )*AIJ
                     TEMP = C( I )*TEMP - S( I )*AIJ
   90             CONTINUE
                  A( K2, J ) = TEMP
  100          CONTINUE
            ELSE IF( ( DIRECT.EQ.'B' ).OR.( DIRECT.EQ.'b' ) )THEN
               DO 120 J = 1, N
                  TEMP = A( K2, J )
                  DO 110 I = K2 - 1, K1, -1
                     AIJ = A( I, J )
                     A( I, J ) = S( I )*TEMP + C( I )*AIJ
                     TEMP = C( I )*TEMP - S( I )*AIJ
  110             CONTINUE
                  A( K2, J ) = TEMP
  120          CONTINUE
            END IF
         END IF
      ELSE IF( RIGHT )THEN
         IF( ( PIVOT.EQ.'V' ).OR.( PIVOT.EQ.'v' ) )THEN
            IF( ( DIRECT.EQ.'F' ).OR.( DIRECT.EQ.'f' ) )THEN
               DO 140 J = K1, K2 - 1
                  IF( ( C( J ).NE.ONE ).OR.( S( J ).NE.ZERO ) )THEN
                     CTEMP = C( J )
                     STEMP = S( J )
                     DO 130 I = 1, M
                        TEMP = A( I, J + 1 )
                        A( I, J + 1 ) = CTEMP*TEMP - STEMP*A( I, J )
                        A( I, J ) = STEMP*TEMP + CTEMP*A( I, J )
  130                CONTINUE
                  END IF
  140          CONTINUE
            ELSE IF( ( DIRECT.EQ.'B' ).OR.( DIRECT.EQ.'b' ) )THEN
               DO 160 J = K2 - 1, K1, -1
                  IF( ( C( J ).NE.ONE ).OR.( S( J ).NE.ZERO ) )THEN
                     CTEMP = C( J )
                     STEMP = S( J )
                     DO 150 I = M, 1, -1
                        TEMP = A( I, J + 1 )
                        A( I, J + 1 ) = CTEMP*TEMP - STEMP*A( I, J )
                        A( I, J ) = STEMP*TEMP + CTEMP*A( I, J )
  150                CONTINUE
                  END IF
  160          CONTINUE
            END IF
         ELSE IF( ( PIVOT.EQ.'T' ).OR.( PIVOT.EQ.'t' ) )THEN
            IF( ( DIRECT.EQ.'F' ).OR.( DIRECT.EQ.'f' ) )THEN
               DO 180 J = K1 + 1, K2
                  CTEMP = C( J - 1 )
                  STEMP = S( J - 1 )
                  IF( ( CTEMP.NE.ONE ).OR.( STEMP.NE.ZERO ) )THEN
                     DO 170 I = 1, M
                        TEMP = A( I, J )
                        A( I, J ) = CTEMP*TEMP - STEMP*A( I, K1 )
                        A( I, K1 ) = STEMP*TEMP + CTEMP*A( I, K1 )
  170                CONTINUE
                  END IF
  180          CONTINUE
            ELSE IF( ( DIRECT.EQ.'B' ).OR.( DIRECT.EQ.'b' ) )THEN
               DO 200 J = K2, K1 + 1, -1
                  CTEMP = C( J - 1 )
                  STEMP = S( J - 1 )
                  IF( ( CTEMP.NE.ONE ).OR.( STEMP.NE.ZERO ) )THEN
                     DO 190 I = M, 1, -1
                        TEMP = A( I, J )
                        A( I, J ) = CTEMP*TEMP - STEMP*A( I, K1 )
                        A( I, K1 ) = STEMP*TEMP + CTEMP*A( I, K1 )
  190                CONTINUE
                  END IF
  200          CONTINUE
            END IF
         ELSE IF( ( PIVOT.EQ.'B' ).OR.( PIVOT.EQ.'b' ) )THEN
            IF( ( DIRECT.EQ.'F' ).OR.( DIRECT.EQ.'f' ) )THEN
               DO 220 J = K1, K2 - 1
                  IF( ( C( J ).NE.ONE ).OR.( S( J ).NE.ZERO ) )THEN
                     CTEMP = C( J )
                     STEMP = S( J )
                     DO 210 I = 1, M
                        TEMP = A( I, J )
                        A( I, J ) = STEMP*A( I, K2 ) + CTEMP*TEMP
                        A( I, K2 ) = CTEMP*A( I, K2 ) - STEMP*TEMP
  210                CONTINUE
                  END IF
  220          CONTINUE
            ELSE IF( ( DIRECT.EQ.'B' ).OR.( DIRECT.EQ.'b' ) )THEN
               DO 240 J = K2 - 1, K1, -1
                  IF( ( C( J ).NE.ONE ).OR.( S( J ).NE.ZERO ) )THEN
                     CTEMP = C( J )
                     STEMP = S( J )
                     DO 230 I = M, 1, -1
                        TEMP = A( I, J )
                        A( I, J ) = STEMP*A( I, K2 ) + CTEMP*TEMP
                        A( I, K2 ) = CTEMP*A( I, K2 ) - STEMP*TEMP
  230                CONTINUE
                  END IF
  240          CONTINUE
            END IF
         END IF
      END IF
C
      RETURN
C
C     End of F06VXF. ( CSGESR  )
C
      END
*
      SUBROUTINE F06ZRF(UPLO,TRANS,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C
C  Purpose
C  =======
C
C  ZHER2K  performs one of the hermitian rank 2k operations
C
C     C := alpha*A*conjg( B' ) + conjg( alpha )*B*conjg( A' ) + beta*C,
C
C  or
C
C     C := alpha*conjg( A' )*B + conjg( alpha )*conjg( B' )*A + beta*C,
C
C  where  alpha and beta  are scalars with  beta  real,  C is an  n by n
C  hermitian matrix and  A and B  are  n by k matrices in the first case
C  and  k by n  matrices in the second case.
C
C  Parameters
C  ==========
C
C  UPLO   - CHARACTER*1.
C           On  entry,   UPLO  specifies  whether  the  upper  or  lower
C           triangular  part  of the  array  C  is to be  referenced  as
C           follows:
C
C              UPLO = 'U' or 'u'   Only the  upper triangular part of  C
C                                  is to be referenced.
C
C              UPLO = 'L' or 'l'   Only the  lower triangular part of  C
C                                  is to be referenced.
C
C           Unchanged on exit.
C
C  TRANS  - CHARACTER*1.
C           On entry,  TRANS  specifies the operation to be performed as
C           follows:
C
C              TRANS = 'N' or 'n'    C := alpha*A*conjg( B' )          +
C                                         conjg( alpha )*B*conjg( A' ) +
C                                         beta*C.
C
C              TRANS = 'C' or 'c'    C := alpha*conjg( A' )*B          +
C                                         conjg( alpha )*conjg( B' )*A +
C                                         beta*C.
C
C           Unchanged on exit.
C
C  N      - INTEGER.
C           On entry,  N specifies the order of the matrix C.  N must be
C           at least zero.
C           Unchanged on exit.
C
C  K      - INTEGER.
C           On entry with  TRANS = 'N' or 'n',  K  specifies  the number
C           of  columns  of the  matrices  A and B,  and on  entry  with
C           TRANS = 'C' or 'c',  K  specifies  the number of rows of the
C           matrices  A and B.  K must be at least zero.
C           Unchanged on exit.
C
C  ALPHA  - COMPLEX         .
C           On entry, ALPHA specifies the scalar alpha.
C           Unchanged on exit.
C
C  A      - COMPLEX          array of DIMENSION ( LDA, ka ), where ka is
C           k  when  TRANS = 'N' or 'n',  and is  n  otherwise.
C           Before entry with  TRANS = 'N' or 'n',  the  leading  n by k
C           part of the array  A  must contain the matrix  A,  otherwise
C           the leading  k by n  part of the array  A  must contain  the
C           matrix A.
C           Unchanged on exit.
C
C  LDA    - INTEGER.
C           On entry, LDA specifies the first dimension of A as declared
C           in  the  calling  (sub)  program.   When  TRANS = 'N' or 'n'
C           then  LDA must be at least  max( 1, n ), otherwise  LDA must
C           be at least  max( 1, k ).
C           Unchanged on exit.
C
C  B      - COMPLEX          array of DIMENSION ( LDB, kb ), where kb is
C           k  when  TRANS = 'N' or 'n',  and is  n  otherwise.
C           Before entry with  TRANS = 'N' or 'n',  the  leading  n by k
C           part of the array  B  must contain the matrix  B,  otherwise
C           the leading  k by n  part of the array  B  must contain  the
C           matrix B.
C           Unchanged on exit.
C
C  LDB    - INTEGER.
C           On entry, LDB specifies the first dimension of B as declared
C           in  the  calling  (sub)  program.   When  TRANS = 'N' or 'n'
C           then  LDB must be at least  max( 1, n ), otherwise  LDB must
C           be at least  max( 1, k ).
C           Unchanged on exit.
C
C  BETA   - REAL            .
C           On entry, BETA specifies the scalar beta.
C           Unchanged on exit.
C
C  C      - COMPLEX          array of DIMENSION ( LDC, n ).
C           Before entry  with  UPLO = 'U' or 'u',  the leading  n by n
C           upper triangular part of the array C must contain the upper
C           triangular part  of the  hermitian matrix  and the strictly
C           lower triangular part of C is not referenced.  On exit, the
C           upper triangular part of the array  C is overwritten by the
C           upper triangular part of the updated matrix.
C           Before entry  with  UPLO = 'L' or 'l',  the leading  n by n
C           lower triangular part of the array C must contain the lower
C           triangular part  of the  hermitian matrix  and the strictly
C           upper triangular part of C is not referenced.  On exit, the
C           lower triangular part of the array  C is overwritten by the
C           lower triangular part of the updated matrix.
C           Note that the imaginary parts of the diagonal elements need
C           not be set,  they are assumed to be zero,  and on exit they
C           are set to zero.
C
C  LDC    - INTEGER.
C           On entry, LDC specifies the first dimension of C as declared
C           in  the  calling  (sub)  program.   LDC  must  be  at  least
C           max( 1, n ).
C           Unchanged on exit.
C
C
C  Level 3 Blas routine.
C
C  -- Written on 8-February-1989.
C     Jack Dongarra, Argonne National Laboratory.
C     Iain Duff, AERE Harwell.
C     Jeremy Du Croz, Numerical Algorithms Group Ltd.
C     Sven Hammarling, Numerical Algorithms Group Ltd.
C
C
C     .. Entry Points ..
      ENTRY             ZHER2K(UPLO,TRANS,N,K,ALPHA,A,LDA,B,LDB,BETA,C,
     *                  LDC)
C     .. Parameters ..
      DOUBLE PRECISION  ONE
      PARAMETER         (ONE=1.0D+0)
      COMPLEX*16        ZERO
      PARAMETER         (ZERO=(0.0D+0,0.0D+0))
C     .. Scalar Arguments ..
      COMPLEX*16        ALPHA
      DOUBLE PRECISION  BETA
      INTEGER           K, LDA, LDB, LDC, N
      CHARACTER*1       TRANS, UPLO
C     .. Array Arguments ..
      COMPLEX*16        A(LDA,*), B(LDB,*), C(LDC,*)
C     .. Local Scalars ..
      COMPLEX*16        TEMP1, TEMP2
      INTEGER           I, INFO, J, L, NROWA
      LOGICAL           UPPER
C     .. External Subroutines ..
      EXTERNAL          F06AAZ
C     .. Intrinsic Functions ..
      INTRINSIC         DBLE, DCONJG, MAX
C     .. Executable Statements ..
C
C     Test the input parameters.
C
      IF ((TRANS.EQ.'N' .OR. TRANS.EQ.'n')) THEN
         NROWA = N
      ELSE
         NROWA = K
      END IF
      UPPER = (UPLO.EQ.'U' .OR. UPLO.EQ.'u')
C
      INFO = 0
      IF (( .NOT. UPPER) .AND. ( .NOT. (UPLO.EQ.'L' .OR. UPLO.EQ.'l')))
     *    THEN
         INFO = 1
      ELSE IF (( .NOT. (TRANS.EQ.'N' .OR. TRANS.EQ.'n'))
     *         .AND. ( .NOT. (TRANS.EQ.'C' .OR. TRANS.EQ.'c'))) THEN
         INFO = 2
      ELSE IF (N.LT.0) THEN
         INFO = 3
      ELSE IF (K.LT.0) THEN
         INFO = 4
      ELSE IF (LDA.LT.MAX(1,NROWA)) THEN
         INFO = 7
      ELSE IF (LDB.LT.MAX(1,NROWA)) THEN
         INFO = 9
      ELSE IF (LDC.LT.MAX(1,N)) THEN
         INFO = 12
      END IF
      IF (INFO.NE.0) THEN
         CALL F06AAZ('F06ZRF/ZHER2K',INFO)
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF ((N.EQ.0) .OR. (((ALPHA.EQ.ZERO) .OR. (K.EQ.0))
     *    .AND. (BETA.EQ.ONE))) RETURN
C
C     And when  alpha.eq.zero.
C
      IF (ALPHA.EQ.ZERO) THEN
         IF (UPPER) THEN
            IF (BETA.EQ.DBLE(ZERO)) THEN
               DO 40 J = 1, N
                  DO 20 I = 1, J
                     C(I,J) = ZERO
   20             CONTINUE
   40          CONTINUE
            ELSE
               DO 80 J = 1, N
                  DO 60 I = 1, J - 1
                     C(I,J) = BETA*C(I,J)
   60             CONTINUE
                  C(J,J) = BETA*DBLE(C(J,J))
   80          CONTINUE
            END IF
         ELSE
            IF (BETA.EQ.DBLE(ZERO)) THEN
               DO 120 J = 1, N
                  DO 100 I = J, N
                     C(I,J) = ZERO
  100             CONTINUE
  120          CONTINUE
            ELSE
               DO 160 J = 1, N
                  C(J,J) = BETA*DBLE(C(J,J))
                  DO 140 I = J + 1, N
                     C(I,J) = BETA*C(I,J)
  140             CONTINUE
  160          CONTINUE
            END IF
         END IF
         RETURN
      END IF
C
C     Start the operations.
C
      IF ((TRANS.EQ.'N' .OR. TRANS.EQ.'n')) THEN
C
C        Form  C := alpha*A*conjg( B' ) + conjg( alpha )*B*conjg( A' ) +
C                   C.
C
         IF (UPPER) THEN
            DO 260 J = 1, N
               IF (BETA.EQ.DBLE(ZERO)) THEN
                  DO 180 I = 1, J
                     C(I,J) = ZERO
  180             CONTINUE
               ELSE IF (BETA.NE.ONE) THEN
                  DO 200 I = 1, J - 1
                     C(I,J) = BETA*C(I,J)
  200             CONTINUE
                  C(J,J) = BETA*DBLE(C(J,J))
               END IF
               DO 240 L = 1, K
                  IF ((A(J,L).NE.ZERO) .OR. (B(J,L).NE.ZERO)) THEN
                     TEMP1 = ALPHA*DCONJG(B(J,L))
                     TEMP2 = DCONJG(ALPHA*A(J,L))
                     DO 220 I = 1, J - 1
                        C(I,J) = C(I,J) + A(I,L)*TEMP1 + B(I,L)*TEMP2
  220                CONTINUE
                     C(J,J) = DBLE(C(J,J)) + DBLE(A(J,L)*TEMP1+B(J,L)
     *                        *TEMP2)
                  END IF
  240          CONTINUE
  260       CONTINUE
         ELSE
            DO 360 J = 1, N
               IF (BETA.EQ.DBLE(ZERO)) THEN
                  DO 280 I = J, N
                     C(I,J) = ZERO
  280             CONTINUE
               ELSE IF (BETA.NE.ONE) THEN
                  DO 300 I = J + 1, N
                     C(I,J) = BETA*C(I,J)
  300             CONTINUE
                  C(J,J) = BETA*DBLE(C(J,J))
               END IF
               DO 340 L = 1, K
                  IF ((A(J,L).NE.ZERO) .OR. (B(J,L).NE.ZERO)) THEN
                     TEMP1 = ALPHA*DCONJG(B(J,L))
                     TEMP2 = DCONJG(ALPHA*A(J,L))
                     DO 320 I = J + 1, N
                        C(I,J) = C(I,J) + A(I,L)*TEMP1 + B(I,L)*TEMP2
  320                CONTINUE
                     C(J,J) = DBLE(C(J,J)) + DBLE(A(J,L)*TEMP1+B(J,L)
     *                        *TEMP2)
                  END IF
  340          CONTINUE
  360       CONTINUE
         END IF
      ELSE
C
C        Form  C := alpha*conjg( A' )*B + conjg( alpha )*conjg( B' )*A +
C                   C.
C
         IF (UPPER) THEN
            DO 420 J = 1, N
               DO 400 I = 1, J
                  TEMP1 = ZERO
                  TEMP2 = ZERO
                  DO 380 L = 1, K
                     TEMP1 = TEMP1 + DCONJG(A(L,I))*B(L,J)
                     TEMP2 = TEMP2 + DCONJG(B(L,I))*A(L,J)
  380             CONTINUE
                  IF (I.EQ.J) THEN
                     IF (BETA.EQ.DBLE(ZERO)) THEN
                        C(J,J) = DBLE(ALPHA*TEMP1+DCONJG(ALPHA)*TEMP2)
                     ELSE
                        C(J,J) = BETA*DBLE(C(J,J)) +
     *                           DBLE(ALPHA*TEMP1+DCONJG(ALPHA)*TEMP2)
                     END IF
                  ELSE
                     IF (BETA.EQ.DBLE(ZERO)) THEN
                        C(I,J) = ALPHA*TEMP1 + DCONJG(ALPHA)*TEMP2
                     ELSE
                        C(I,J) = BETA*C(I,J) + ALPHA*TEMP1 +
     *                           DCONJG(ALPHA)*TEMP2
                     END IF
                  END IF
  400          CONTINUE
  420       CONTINUE
         ELSE
            DO 480 J = 1, N
               DO 460 I = J, N
                  TEMP1 = ZERO
                  TEMP2 = ZERO
                  DO 440 L = 1, K
                     TEMP1 = TEMP1 + DCONJG(A(L,I))*B(L,J)
                     TEMP2 = TEMP2 + DCONJG(B(L,I))*A(L,J)
  440             CONTINUE
                  IF (I.EQ.J) THEN
                     IF (BETA.EQ.DBLE(ZERO)) THEN
                        C(J,J) = DBLE(ALPHA*TEMP1+DCONJG(ALPHA)*TEMP2)
                     ELSE
                        C(J,J) = BETA*DBLE(C(J,J)) +
     *                           DBLE(ALPHA*TEMP1+DCONJG(ALPHA)*TEMP2)
                     END IF
                  ELSE
                     IF (BETA.EQ.DBLE(ZERO)) THEN
                        C(I,J) = ALPHA*TEMP1 + DCONJG(ALPHA)*TEMP2
                     ELSE
                        C(I,J) = BETA*C(I,J) + ALPHA*TEMP1 +
     *                           DCONJG(ALPHA)*TEMP2
                     END IF
                  END IF
  460          CONTINUE
  480       CONTINUE
         END IF
      END IF
C
      RETURN
C
C     End of F06ZRF (ZHER2K).
C
      END
*
      SUBROUTINE F07FRY(N,X,INCX)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C     ENTRY             ZLACGV(N,X,INCX)
C
C  Purpose
C  =======
C
C     ZLACGV conjugates a complex vector of length N.
C
C  Arguments
C  =========
C
C  N      - INTEGER
C           On entry, N specifies the length of the vector X.
C           Unchanged on exit.
C
C  X      - COMPLEX array, dimension( N )
C           On entry, X contains the vector to be conjugated.
C           On return, X contains conjg(X).
C
C  INCX   - INTEGER
C           On entry, INCX specifies the spacing betweeen successive
C           elements of X.
C           Unchanged on exit.
C
C  -- LAPACK auxiliary routine (adapted for NAG Library)
C     Univ. of Tennessee, Oak Ridge National Lab, Argonne National Lab,
C     Courant Institute, NAG Ltd., and Rice University
C
C     .. Scalar Arguments ..
      INTEGER           INCX, N
C     .. Array Arguments ..
      COMPLEX*16        X(*)
C     .. Local Scalars ..
      INTEGER           I, IOFF
C     .. Intrinsic Functions ..
      INTRINSIC         DCONJG
C     .. Executable Statements ..
C
      IF (INCX.EQ.1) THEN
         DO 20 I = 1, N
            X(I) = DCONJG(X(I))
   20    CONTINUE
      ELSE
         IOFF = 1
         IF (INCX.LT.0) IOFF = 1 - (N-1)*INCX
         DO 40 I = 1, N
            X(IOFF) = DCONJG(X(IOFF))
            IOFF = IOFF + INCX
   40    CONTINUE
      END IF
      RETURN
C
C     End of F07FRY (ZLACGV)
C
      END
*
      SUBROUTINE F07MDX(A,B,C,RT1,RT2,CS1,SN1)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C     ENTRY             DLAEV2(A,B,C,RT1,RT2,CS1,SN1)
C
C  Purpose
C  =======
C
C  DLAEV2 computes the eigendecomposition of a 2 by 2 symmetric matrix:
C     ((A,B);(B,C))
C
C  RT1 is the eigenvalue of larger absolute value, RT2 is the eigenvalue
C  of smaller absolute value, and (CS1,SN1) is the unit right
C  eigenvector for RT1, giving the decomposition
C
C     [ CS1  SN1 ]  .  [ A B ]  .  [ CS1 -SN1 ]  =  [ RT1  0  ]
C     [-SN1  CS1 ]     [ B C ]     [ SN1  CS1 ]     [  0  RT2 ]
C
C  RT1 is accurate to a few ulps barring over/underflow.
C  RT2 may be inaccurate if there is massive cancellation in the
C     determinant A*C-B*B; higher precision or correctly rounded or
C     correctly truncated arithmetic would be needed to compute RT2
C     accurately in all cases.
C  CS1 and SN1 are accurate to a few ulps barring over/underflow.
C  Overflow is possible only if RT1 is within a factor of 5 of overflow.
C  Underflow is harmless if the input data is 0 or exceeds
C     underflow_threshold / macheps.
C
C  Arguments
C  =========
C
C  A      (input) REAL
C         The (1,1) entry of the input matrix.
C
C  B      (input) REAL
C         The (1,2) entry and the conjugate of the (2,1) entry of the
C         input matrix.
C
C  C      (input) REAL
C         The (2,2) entry of the input matrix.
C
C  RT1    (output) REAL
C         The eigenvalue of larger absolute value.
C
C  RT2    (output) REAL
C         The eigenvalue of smaller absolute value.
C
C  CS1    (output) REAL
C  SN1    (output) REAL
C         The vector (CS1, SN1) is a unit right eigenvector for RT1.
C
C
C  -- LAPACK auxiliary routine (adapted for NAG Library)
C     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
C     Courant Institute, Argonne National Lab, and Rice University
C
C     .. Parameters ..
      DOUBLE PRECISION  ONE
      PARAMETER         (ONE=1.0D0)
      DOUBLE PRECISION  TWO
      PARAMETER         (TWO=2.0D0)
      DOUBLE PRECISION  ZERO
      PARAMETER         (ZERO=0.0D0)
      DOUBLE PRECISION  HALF
      PARAMETER         (HALF=0.5D0)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  A, B, C, CS1, RT1, RT2, SN1
C     .. Local Scalars ..
      DOUBLE PRECISION  AB, ACMN, ACMX, ACS, ADF, CS, CT, DF, RT, SM,
     *                  TB, TN
      INTEGER           SGN1, SGN2
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, SQRT
C     .. Executable Statements ..
C
C     Compute the eigenvalues
C
      SM = A + C
      DF = A - C
      ADF = ABS(DF)
      TB = B + B
      AB = ABS(TB)
      IF (ABS(A).GT.ABS(C)) THEN
         ACMX = A
         ACMN = C
      ELSE
         ACMX = C
         ACMN = A
      END IF
      IF (ADF.GT.AB) THEN
         RT = ADF*SQRT(ONE+(AB/ADF)**2)
      ELSE IF (ADF.LT.AB) THEN
         RT = AB*SQRT(ONE+(ADF/AB)**2)
      ELSE
C
C        Includes case AB=ADF=0
C
         RT = AB*SQRT(TWO)
      END IF
      IF (SM.LT.ZERO) THEN
         RT1 = HALF*(SM-RT)
         SGN1 = -1
C
C        Order of execution important.
C        To get fully accurate smaller eigenvalue,
C        next line needs to be executed in higher precision.
C
         RT2 = (ACMX/RT1)*ACMN - (B/RT1)*B
      ELSE IF (SM.GT.ZERO) THEN
         RT1 = HALF*(SM+RT)
         SGN1 = 1
C
C        Order of execution important.
C        To get fully accurate smaller eigenvalue,
C        next line needs to be executed in higher precision.
C
         RT2 = (ACMX/RT1)*ACMN - (B/RT1)*B
      ELSE
C
C        Includes case RT1 = RT2 = 0
C
         RT1 = HALF*RT
         RT2 = -HALF*RT
         SGN1 = 1
      END IF
C
C     Compute the eigenvector
C
      IF (DF.GE.ZERO) THEN
         CS = DF + RT
         SGN2 = 1
      ELSE
         CS = DF - RT
         SGN2 = -1
      END IF
      ACS = ABS(CS)
      IF (ACS.GT.AB) THEN
         CT = -TB/CS
         SN1 = ONE/SQRT(ONE+CT*CT)
         CS1 = CT*SN1
      ELSE
         IF (AB.EQ.ZERO) THEN
            CS1 = ONE
            SN1 = ZERO
         ELSE
            TN = -CS/TB
            CS1 = ONE/SQRT(ONE+TN*TN)
            SN1 = TN*CS1
         END IF
      END IF
      IF (SGN1.EQ.SGN2) THEN
         TN = CS1
         CS1 = -SN1
         SN1 = TN
      END IF
      RETURN
C
C     End of F07MDX (DLAEV2)
C
      END
*
      INTEGER FUNCTION F07ZAY(NAME)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C     MARK 16 RE-ISSUE. NAG COPYRIGHT 1993.
C
C     F07ZAY returns a unique positive integer code
C     corresponding to a six-letter NAG routine name
C     given in NAME. If NAME is not recognised, 0 is
C     returned.
C
C     Modified at Mark 16 to allow calls from F08 routines.
C
C     .. Scalar Arguments ..
      CHARACTER*6             NAME
C     .. Local Scalars ..
      INTEGER                 J, K
      CHARACTER               NAME4, NAME5
      CHARACTER*3             NAME13
C     .. Executable Statements ..
C
      NAME13 = NAME(1:3)
      NAME4 = NAME(4:4)
      NAME5 = NAME(5:5)
C
      IF (NAME13.EQ.'F07') THEN
C
         IF (NAME4.EQ.'A') THEN
            J = 0
         ELSE IF (NAME4.EQ.'B') THEN
            J = 1
         ELSE IF (NAME4.EQ.'F') THEN
            J = 2
         ELSE IF (NAME4.EQ.'H') THEN
            J = 3
         ELSE IF (NAME4.EQ.'M') THEN
            J = 4
         ELSE IF (NAME4.EQ.'N') THEN
            J = 5
         ELSE IF (NAME4.EQ.'T') THEN
            J = 6
         ELSE
            J = -1
         END IF
C
         IF (NAME5.EQ.'D') THEN
            K = 0
         ELSE IF (NAME5.EQ.'J') THEN
            K = 1
         ELSE IF (NAME5.EQ.'R') THEN
            K = 2
         ELSE IF (NAME5.EQ.'W') THEN
            K = 3
         ELSE
            K = -1
         END IF
C
         IF (J.LT.0 .OR. K.LT.0) THEN
            F07ZAY = 0
         ELSE
C           F07ZAY is in the range 1-28 for F07 routines.
            F07ZAY = 1 + 4*J + K
         END IF
C
      ELSE IF (NAME13.EQ.'F08') THEN
C
         IF (NAME4.EQ.'A') THEN
            J = 0
         ELSE IF (NAME4.EQ.'C') THEN
            J = 1
         ELSE IF (NAME4.EQ.'F') THEN
            J = 2
         ELSE IF (NAME4.EQ.'J') THEN
            J = 3
         ELSE IF (NAME4.EQ.'K') THEN
            J = 4
         ELSE IF (NAME4.EQ.'N') THEN
            J = 5
         ELSE IF (NAME4.EQ.'P') THEN
            J = 6
         ELSE IF (NAME4.EQ.'S') THEN
            J = 7
         ELSE
            J = -1
         END IF
C
         IF (NAME5.EQ.'E') THEN
            K = 0
         ELSE IF (NAME5.EQ.'F') THEN
            K = 1
         ELSE IF (NAME5.EQ.'G') THEN
            K = 2
         ELSE IF (NAME5.EQ.'H') THEN
            K = 3
         ELSE IF (NAME5.EQ.'J') THEN
            K = 4
         ELSE IF (NAME5.EQ.'K') THEN
            K = 5
         ELSE IF (NAME5.EQ.'S') THEN
            K = 6
         ELSE IF (NAME5.EQ.'T') THEN
            K = 7
         ELSE IF (NAME5.EQ.'U') THEN
            K = 8
         ELSE IF (NAME5.EQ.'V') THEN
            K = 9
         ELSE IF (NAME5.EQ.'W') THEN
            K = 10
         ELSE IF (NAME5.EQ.'X') THEN
            K = 11
         ELSE
            K = -1
         END IF
C
         IF (J.LT.0 .OR. K.LT.0) THEN
            F07ZAY = 0
         ELSE
C           F07ZAY is in the range 29-124 for F08 routines.
            F07ZAY = 29 + 12*J + K
         END IF
C
      ELSE
         F07ZAY = 0
      END IF
C
      RETURN
C
      END
*
      SUBROUTINE F07ZAZ(ISPEC,NAME,IVAL,RWFLAG)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C     MARK 16 RE-ISSUE. NAG COPYRIGHT 1993.
C
C     Modified at Mark 16 to allow calls from F08
C     routines.
C
C  -- NAG version of LAPACK auxiliary routine ILAENV
C
C  This is an unimplemented but functional version
C  of F07ZAZ, which returns the value 1 in IVAL for
C  all calls with RWFLAG = 0 and ISPEC/NAME
C  combinations which will be used by the Fortran
C  Library.
C
C  The parameter MAXIC is set to the largest value
C  that can be returned by F07ZAY. NCODES is the
C  number of Fortran Library routines that are
C  expected to call F07ZAZ. NSPECS is the number
C  of values that ISPEC can take.
C
C  Purpose
C  =======
C
C  F07ZAZ sets or returns problem-dependent
C  parameters for the local environment. See
C  ISPEC for a description of the parameters.
C
C  The problem-dependent parameters are contained
C  in the integer array IPARMS, and the value with
C  index ISPEC is set or copied to IVAL.
C
C  Arguments
C  =========
C
C  ISPEC (input) INTEGER
C     Specifies the parameter to be set or
C     returned by F07ZAZ.
C     = 1: the optimal blocksize; if this value
C          is 1, an unblocked algorithm will give
C          the best performance.
C     = 2: the minimum block size for which the
C          block routine should be used; if the
C          usable block size is less than this
C          value, an unblocked routine should be
C          used.
C     = 3: the crossover point (in a block
C          routine, for N less than this value,
C          an unblocked routine should be used).
C     = 4: the number of shifts, used in the
C          nonsymmetric eigenvalue routines.
C     = 5: the minimum column dimension for
C          blocking to be used; rectangular
C          blocks must have dimension at least
C          k by m, where k is given by
C          F07ZAZ(2,...) and m by F07ZAZ(5,...).
C     = 6: the crossover point for the SVD (when
C          reducing an m by n matrix to bidiagonal
C          form, if max(m,n)/min(m,n) exceeds this
C          value, a QR factorization is used first
C          to reduce the matrix to a triangular
C          form).
C     = 7: the number of processors.
C     = 8: the crossover point for the multishift
C          QR and QZ methods for nonsymmetric
C          eigenvalue problems.
C
C  NAME  (input) CHARACTER*(*)
C     The name of the calling subroutine.
C
C  IVAL  (input/output) INTEGER
C     the value of the parameter set or returned.
C
C  FLAG  (input) INTEGER
C     = 0: F07ZAZ returns in IVAL the value of
C          the parameter specified by ISPEC.
C     = 1: F07ZAZ sets the parameter specified
C          by ISPEC to the value in IVAL.
C
C  ==============================================
C
C     .. Parameters ..
      INTEGER           NSPECS, NCODES, MAXIC
      PARAMETER         (NSPECS=8,NCODES=58,MAXIC=124)
C     .. Scalar Arguments ..
      INTEGER           ISPEC, IVAL, RWFLAG
      CHARACTER*(*)     NAME
C     .. Local Scalars ..
      INTEGER           I, ICODE
C     .. Local Arrays ..
      INTEGER           IPARMS(NSPECS,NCODES), POINT(0:MAXIC)
C     .. External Functions ..
      INTEGER           F07ZAY
      EXTERNAL          F07ZAY
C     .. Save statement ..
      SAVE              IPARMS, POINT
C     .. Data statements ..
      DATA              (IPARMS(I,1),I=1,NSPECS)/1, 0, 0, 0, 0, 0, 0, 0/
      DATA              (IPARMS(I,2),I=1,NSPECS)/1, 1, 0, 0, 0, 0, 0, 0/
      DATA              (IPARMS(I,3),I=1,NSPECS)/1, 0, 0, 0, 0, 0, 0, 0/
      DATA              (IPARMS(I,4),I=1,NSPECS)/1, 1, 0, 0, 0, 0, 0, 0/
      DATA              (IPARMS(I,5),I=1,NSPECS)/1, 0, 0, 0, 0, 0, 0, 0/
      DATA              (IPARMS(I,6),I=1,NSPECS)/1, 0, 0, 0, 0, 0, 0, 0/
      DATA              (IPARMS(I,7),I=1,NSPECS)/1, 0, 0, 0, 0, 0, 0, 0/
      DATA              (IPARMS(I,8),I=1,NSPECS)/1, 0, 0, 0, 0, 0, 0, 0/
      DATA              (IPARMS(I,9),I=1,NSPECS)/1, 0, 0, 0, 0, 0, 0, 0/
      DATA              (IPARMS(I,10),I=1,NSPECS)/1, 0, 0, 0, 0, 0, 0,
     *                  0/
      DATA              (IPARMS(I,11),I=1,NSPECS)/1, 0, 0, 0, 0, 0, 0,
     *                  0/
      DATA              (IPARMS(I,12),I=1,NSPECS)/1, 0, 0, 0, 0, 0, 0,
     *                  0/
      DATA              (IPARMS(I,13),I=1,NSPECS)/1, 1, 0, 0, 0, 0, 0,
     *                  0/
      DATA              (IPARMS(I,14),I=1,NSPECS)/1, 1, 0, 0, 0, 0, 0,
     *                  0/
      DATA              (IPARMS(I,15),I=1,NSPECS)/1, 1, 0, 0, 0, 0, 0,
     *                  0/
      DATA              (IPARMS(I,16),I=1,NSPECS)/1, 0, 0, 0, 0, 0, 0,
     *                  0/
      DATA              (IPARMS(I,17),I=1,NSPECS)/1, 0, 0, 0, 0, 0, 0,
     *                  0/
      DATA              (IPARMS(I,18),I=1,NSPECS)/1, 1, 1, 0, 0, 0, 0,
     *                  0/
      DATA              (IPARMS(I,19),I=1,NSPECS)/1, 1, 1, 0, 0, 0, 0,
     *                  0/
      DATA              (IPARMS(I,20),I=1,NSPECS)/1, 1, 0, 0, 0, 0, 0,
     *                  0/
      DATA              (IPARMS(I,21),I=1,NSPECS)/1, 1, 1, 0, 0, 0, 0,
     *                  0/
      DATA              (IPARMS(I,22),I=1,NSPECS)/1, 1, 1, 0, 0, 0, 0,
     *                  0/
      DATA              (IPARMS(I,23),I=1,NSPECS)/1, 1, 0, 0, 0, 0, 0,
     *                  0/
      DATA              (IPARMS(I,24),I=1,NSPECS)/1, 1, 1, 0, 0, 0, 0,
     *                  0/
      DATA              (IPARMS(I,25),I=1,NSPECS)/1, 1, 1, 0, 0, 0, 0,
     *                  0/
      DATA              (IPARMS(I,26),I=1,NSPECS)/1, 1, 0, 0, 0, 0, 0,
     *                  0/
      DATA              (IPARMS(I,27),I=1,NSPECS)/1, 1, 1, 0, 0, 0, 0,
     *                  0/
      DATA              (IPARMS(I,28),I=1,NSPECS)/1, 1, 1, 0, 0, 0, 0,
     *                  0/
      DATA              (IPARMS(I,29),I=1,NSPECS)/1, 1, 0, 0, 0, 0, 0,
     *                  0/
      DATA              (IPARMS(I,30),I=1,NSPECS)/1, 1, 1, 0, 0, 0, 0,
     *                  0/
      DATA              (IPARMS(I,31),I=1,NSPECS)/1, 1, 1, 0, 0, 0, 0,
     *                  0/
      DATA              (IPARMS(I,32),I=1,NSPECS)/1, 1, 0, 0, 0, 0, 0,
     *                  0/
      DATA              (IPARMS(I,33),I=1,NSPECS)/1, 1, 1, 0, 0, 0, 0,
     *                  0/
      DATA              (IPARMS(I,34),I=1,NSPECS)/1, 1, 1, 0, 0, 0, 0,
     *                  0/
      DATA              (IPARMS(I,35),I=1,NSPECS)/1, 1, 0, 0, 0, 0, 0,
     *                  0/
      DATA              (IPARMS(I,36),I=1,NSPECS)/1, 1, 1, 0, 0, 0, 0,
     *                  0/
      DATA              (IPARMS(I,37),I=1,NSPECS)/1, 1, 1, 0, 0, 0, 0,
     *                  0/
      DATA              (IPARMS(I,38),I=1,NSPECS)/1, 1, 0, 0, 0, 0, 0,
     *                  0/
      DATA              (IPARMS(I,39),I=1,NSPECS)/1, 1, 1, 0, 0, 0, 0,
     *                  0/
      DATA              (IPARMS(I,40),I=1,NSPECS)/1, 1, 1, 0, 0, 0, 0,
     *                  0/
      DATA              (IPARMS(I,41),I=1,NSPECS)/1, 1, 0, 0, 0, 0, 0,
     *                  0/
      DATA              (IPARMS(I,42),I=1,NSPECS)/1, 1, 1, 0, 0, 0, 0,
     *                  0/
      DATA              (IPARMS(I,43),I=1,NSPECS)/1, 1, 1, 0, 0, 0, 0,
     *                  0/
      DATA              (IPARMS(I,44),I=1,NSPECS)/1, 0, 0, 0, 0, 0, 0,
     *                  0/
      DATA              (IPARMS(I,45),I=1,NSPECS)/1, 1, 1, 0, 0, 0, 0,
     *                  0/
      DATA              (IPARMS(I,46),I=1,NSPECS)/1, 0, 0, 0, 0, 0, 0,
     *                  0/
      DATA              (IPARMS(I,47),I=1,NSPECS)/1, 0, 0, 0, 0, 0, 0,
     *                  0/
      DATA              (IPARMS(I,48),I=1,NSPECS)/1, 1, 1, 0, 0, 0, 0,
     *                  0/
      DATA              (IPARMS(I,49),I=1,NSPECS)/1, 0, 0, 0, 0, 0, 0,
     *                  0/
      DATA              (IPARMS(I,50),I=1,NSPECS)/1, 0, 0, 0, 0, 0, 0,
     *                  0/
      DATA              (IPARMS(I,51),I=1,NSPECS)/1, 1, 1, 0, 0, 0, 0,
     *                  0/
      DATA              (IPARMS(I,52),I=1,NSPECS)/1, 0, 0, 0, 0, 0, 0,
     *                  0/
      DATA              (IPARMS(I,53),I=1,NSPECS)/1, 1, 1, 0, 0, 0, 0,
     *                  0/
      DATA              (IPARMS(I,54),I=1,NSPECS)/1, 0, 0, 0, 0, 0, 0,
     *                  0/
      DATA              (IPARMS(I,55),I=1,NSPECS)/0, 0, 0, 1, 0, 0, 0,
     *                  1/
      DATA              (IPARMS(I,56),I=1,NSPECS)/0, 0, 0, 1, 0, 0, 0,
     *                  1/
      DATA              (IPARMS(I,57),I=1,NSPECS)/1, 0, 0, 0, 0, 0, 0,
     *                  0/
      DATA              (IPARMS(I,58),I=1,NSPECS)/1, 0, 0, 0, 0, 0, 0,
     *                  0/
      DATA              POINT/0, 1, 2, 3, 4, 5, 0, 6, 0, 7, 8, 9, 10,
     *                  11, 0, 12, 0, 13, 0, 14, 0, 0, 0, 15, 0, 0, 16,
     *                  0, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27,
     *                  28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39,
     *                  40, 41, 42, 0, 0, 0, 0, 0, 43, 0, 0, 0, 0, 0, 0,
     *                  0, 0, 0, 44, 0, 0, 0, 0, 0, 0, 0, 45, 46, 47, 0,
     *                  0, 0, 48, 49, 50, 0, 0, 0, 51, 52, 0, 0, 0, 0,
     *                  53, 54, 0, 0, 0, 0, 55, 0, 0, 0, 0, 0, 56, 0, 0,
     *                  0, 0, 0, 57, 0, 0, 0, 0, 0, 58, 0, 0, 0, 0, 0/
C     .. Executable Statements ..
C
C     Convert the NAG name to an integer code.
      ICODE = POINT(F07ZAY(NAME))
C
      IF (ISPEC.LT.1 .OR. ISPEC.GT.NSPECS) THEN
C        Invalid value for ISPEC
         IVAL = -1
      ELSE IF (ICODE.EQ.0) THEN
C        Invalid value for NAME
         IVAL = -2
      ELSE IF (RWFLAG.EQ.0) THEN
C        Read the value of a parameter
         IVAL = IPARMS(ISPEC,ICODE)
      ELSE
C        Set the value of a parameter
         IPARMS(ISPEC,ICODE) = IVAL
      END IF
C
      RETURN
C
C     End of F07ZAZ
C
      END
*
      DOUBLE PRECISION FUNCTION F08ASU(X,Y,Z)
C     MARK 16 RELEASE. NAG COPYRIGHT 1992.
C     DOUBLE PRECISION                 DLAPY3
C     ENTRY                            DLAPY3(X,Y,Z)
C
C  Purpose
C  =======
C
C  DLAPY3 returns sqrt(x**2+y**2+z**2), taking care not to cause
C  unnecessary overflow.
C
C  Arguments
C  =========
C
C  X       (input) DOUBLE PRECISION
C  Y       (input) DOUBLE PRECISION
C  Z       (input) DOUBLE PRECISION
C          X, Y and Z specify the values x, y and z.
C
C  -- LAPACK auxiliary routine (adapted for NAG Library)
C     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
C     Courant Institute, Argonne National Lab, and Rice University
C
C     .. Parameters ..
      DOUBLE PRECISION                 ZERO
      PARAMETER                        (ZERO=0.0D0)
C     .. Scalar Arguments ..
      DOUBLE PRECISION                 X, Y, Z
C     .. Local Scalars ..
      DOUBLE PRECISION                 W, XABS, YABS, ZABS
C     .. Intrinsic Functions ..
      INTRINSIC                        ABS, MAX, SQRT
C     .. Executable Statements ..
C
      XABS = ABS(X)
      YABS = ABS(Y)
      ZABS = ABS(Z)
      W = MAX(XABS,YABS,ZABS)
      IF (W.EQ.ZERO) THEN
         F08ASU = ZERO
      ELSE
         F08ASU = W*SQRT((XABS/W)**2+(YABS/W)**2+(ZABS/W)**2)
      END IF
      RETURN
C
C     End of F08ASU (DLAPY3)
C
      END
*
      SUBROUTINE F08ASV(N,ALPHA,X,INCX,TAU)
C     MARK 16 RELEASE. NAG COPYRIGHT 1992.
C     ENTRY             ZLARFG(N,ALPHA,X,INCX,TAU)
C
C  Purpose
C  =======
C
C  ZLARFG generates a complex elementary reflector H of order n, such
C  that
C
C        H' * ( alpha ) = ( beta ),   H' * H = I.
C             (   x   )   (   0  )
C
C  where alpha and beta are scalars, with beta real, and x is an
C  (n-1)-element complex vector. H is represented in the form
C
C        H = I - tau * ( 1 ) * ( 1 v' ) ,
C                      ( v )
C
C  where tau is a complex scalar and v is a complex (n-1)-element
C  vector. Note that H is not hermitian.
C
C  If the elements of x are all zero and alpha is real, then tau = 0
C  and H is taken to be the unit matrix.
C
C  Otherwise  1 <= real(tau) <= 2  and  abs(tau-1) <= 1 .
C
C  Arguments
C  =========
C
C  N       (input) INTEGER
C          The order of the elementary reflector.
C
C  ALPHA   (input/output) COMPLEX*16
C          On entry, the value alpha.
C          On exit, it is overwritten with the value beta.
C
C  X       (input/output) COMPLEX*16 array, dimension
C                         (1+(N-2)*abs(INCX))
C          On entry, the vector x.
C          On exit, it is overwritten with the vector v.
C
C  INCX    (input) INTEGER
C          The increment between elements of X. INCX <> 0.
C
C  TAU     (output) COMPLEX*16
C          The value tau.
C
C  -- LAPACK auxiliary routine (adapted for NAG Library)
C     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
C     Courant Institute, Argonne National Lab, and Rice University
C
C  =====================================================================
C
C     .. Parameters ..
      DOUBLE PRECISION  ONE, ZERO
      PARAMETER         (ONE=1.0D+0,ZERO=0.0D+0)
C     .. Scalar Arguments ..
      COMPLEX*16        ALPHA, TAU
      INTEGER           INCX, N
C     .. Array Arguments ..
      COMPLEX*16        X(*)
C     .. Local Scalars ..
      DOUBLE PRECISION  ALPHI, ALPHR, BETA, RSAFMN, SAFMIN, XNORM
      INTEGER           J, KNT
      LOGICAL           DIVFLG
C     .. External Functions ..
      COMPLEX*16        F06CLF
      DOUBLE PRECISION  DZNRM2, F08ASU, X02AMF
      EXTERNAL          F06CLF, DZNRM2, F08ASU, X02AMF
C     .. External Subroutines ..
      EXTERNAL          ZDSCAL, ZSCAL
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, DBLE, DCMPLX, DIMAG, SIGN
C     .. Executable Statements ..
C
      IF (N.LE.0) THEN
         TAU = ZERO
         RETURN
      END IF
C
      XNORM = DZNRM2(N-1,X,INCX)
      ALPHR = DBLE(ALPHA)
      ALPHI = DIMAG(ALPHA)
C
      IF (XNORM.EQ.ZERO .AND. ALPHI.EQ.ZERO) THEN
C
C        H  =  I
C
         TAU = ZERO
      ELSE
C
C        general case
C
         BETA = -SIGN(F08ASU(ALPHR,ALPHI,XNORM),ALPHR)
         SAFMIN = X02AMF()
         RSAFMN = ONE/SAFMIN
C
         IF (ABS(BETA).LT.SAFMIN) THEN
C
C           XNORM, BETA may be inaccurate; scale X and recompute them
C
            KNT = 0
   20       CONTINUE
            KNT = KNT + 1
            CALL ZDSCAL(N-1,RSAFMN,X,INCX)
            BETA = BETA*RSAFMN
            ALPHI = ALPHI*RSAFMN
            ALPHR = ALPHR*RSAFMN
            IF (ABS(BETA).LT.SAFMIN) GO TO 20
C
C           New BETA is at most 1, at least SAFMIN
C
            XNORM = DZNRM2(N-1,X,INCX)
            ALPHA = DCMPLX(ALPHR,ALPHI)
            BETA = -SIGN(F08ASU(ALPHR,ALPHI,XNORM),ALPHR)
            TAU = DCMPLX((BETA-ALPHR)/BETA,-ALPHI/BETA)
            ALPHA = F06CLF(DCMPLX(ONE),ALPHA-BETA,DIVFLG)
            CALL ZSCAL(N-1,ALPHA,X,INCX)
C
C           If ALPHA is subnormal, it may lose relative accuracy
C
            ALPHA = BETA
            DO 40 J = 1, KNT
               ALPHA = ALPHA*SAFMIN
   40       CONTINUE
         ELSE
            TAU = DCMPLX((BETA-ALPHR)/BETA,-ALPHI/BETA)
            ALPHA = F06CLF(DCMPLX(ONE),ALPHA-BETA,DIVFLG)
            CALL ZSCAL(N-1,ALPHA,X,INCX)
            ALPHA = BETA
         END IF
      END IF
C
      RETURN
C
C     End of F08ASV (ZLARFG)
C
      END
*
      SUBROUTINE F08FSF(UPLO,N,A,LDA,D,E,TAU,WORK,LWORK,INFO)
C     MARK 16 RELEASE. NAG COPYRIGHT 1992.
C     .. Entry Points ..
      ENTRY             ZHETRD(UPLO,N,A,LDA,D,E,TAU,WORK,LWORK,INFO)
C
C  Purpose
C  =======
C
C  ZHETRD reduces a complex Hermitian matrix A to real symmetric
C  tridiagonal form T by a unitary similarity transformation:
C  Q' * A * Q = T.
C
C  Arguments
C  =========
C
C  UPLO    (input) CHARACTER*1
C          Specifies whether the upper or lower triangular part of the
C          Hermitian matrix A is stored:
C          = 'U':  Upper triangular
C          = 'L':  Lower triangular
C
C  N       (input) INTEGER
C          The order of the matrix A.  N >= 0.
C
C  A       (input/output) COMPLEX*16 array, dimension (LDA,N)
C          On entry, the Hermitian matrix A.  If UPLO = 'U', the leading
C          n-by-n upper triangular part of A contains the upper
C          triangular part of the matrix A, and the strictly lower
C          triangular part of A is not referenced.  If UPLO = 'L', the
C          leading n-by-n lower triangular part of A contains the lower
C          triangular part of the matrix A, and the strictly upper
C          triangular part of A is not referenced.
C          On exit, if UPLO = 'U', the diagonal and first superdiagonal
C          of A are overwritten by the corresponding elements of the
C          tridiagonal matrix T, and the elements above the first
C          superdiagonal, with the array TAU, represent the unitary
C          matrix Q as a product of elementary reflectors; if UPLO
C          = 'L', the diagonal and first subdiagonal of A are over-
C          written by the corresponding elements of the tridiagonal
C          matrix T, and the elements below the first subdiagonal, with
C          the array TAU, represent the unitary matrix Q as a product
C          of elementary reflectors. See Further Details.
C
C  LDA     (input) INTEGER
C          The leading dimension of the array A.  LDA >= max(1,N).
C
C  D       (output) DOUBLE PRECISION array, dimension (N)
C          The diagonal elements of the tridiagonal matrix T:
C          D(i) = A(i,i).
C
C  E       (output) DOUBLE PRECISION array, dimension (N-1)
C          The off-diagonal elements of the tridiagonal matrix T:
C          E(i) = A(i,i+1) if UPLO = 'U', E(i) = A(i+1,i) if UPLO = 'L'.
C
C  TAU     (output) COMPLEX*16 array, dimension (N-1)
C          The scalar factors of the elementary reflectors (see Further
C          Details).
C
C  WORK    (workspace) COMPLEX*16 array, dimension (LWORK)
C          On exit, if INFO = 0, WORK(1) returns the minimum value of
C          LWORK required to use the optimal blocksize.
C
C  LWORK   (input) INTEGER
C          The dimension of the array WORK.  LWORK >= 1.
C          For optimum performance LWORK should be at least N*NB,
C          where NB is the optimal blocksize.
C
C  INFO    (output) INTEGER
C          = 0:  successful exit.
C          < 0:  if INFO = -i, the i-th argument had an illegal value.
C
C  Further Details
C  ===============
C
C  If UPLO = 'U', the matrix Q is represented as a product of elementary
C  reflectors
C
C     Q = H(n-1) . . . H(2) H(1).
C
C  Each H(i) has the form
C
C     H(i) = I - tau * v * v'
C
C  where tau is a complex scalar, and v is a complex vector with
C  v(i+1:n) = 0 and v(i) = 1; v(1:i-1) is stored on exit in
C  A(1:i-1,i+1), and tau in TAU(i).
C
C  If UPLO = 'L', the matrix Q is represented as a product of elementary
C  reflectors
C
C     Q = H(1) H(2) . . . H(n-1).
C
C  Each H(i) has the form
C
C     H(i) = I - tau * v * v'
C
C  where tau is a complex scalar, and v is a complex vector with
C  v(1:i) = 0 and v(i+1) = 1; v(i+2:n) is stored on exit in A(i+2:n,i),
C  and tau in TAU(i).
C
C  The contents of A on exit are illustrated by the following examples
C  with n = 5:
C
C  if UPLO = 'U':                       if UPLO = 'L':
C
C    (  d   e   v2  v3  v4 )              (  d                  )
C    (      d   e   v3  v4 )              (  e   d              )
C    (          d   e   v4 )              (  v1  e   d          )
C    (              d   e  )              (  v1  v2  e   d      )
C    (                  d  )              (  v1  v2  v3  e   d  )
C
C  where d and e denote diagonal and off-diagonal elements of T, and vi
C  denotes an element of the vector defining H(i).
C
C  -- LAPACK routine (adapted for NAG Library)
C     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
C     Courant Institute, Argonne National Lab, and Rice University
C
C  =====================================================================
C
C     .. Parameters ..
      DOUBLE PRECISION  ONE
      PARAMETER         (ONE=1.0D0)
      COMPLEX*16        CONE
      PARAMETER         (CONE=1.0D0)
C     .. Scalar Arguments ..
      INTEGER           INFO, LDA, LWORK, N
      CHARACTER         UPLO
C     .. Array Arguments ..
      COMPLEX*16        A(LDA,*), TAU(*), WORK(LWORK)
      DOUBLE PRECISION  D(*), E(*)
C     .. Local Scalars ..
      INTEGER           I, IINFO, IWS, J, KK, LDWORK, NB, NBMIN, NX
      LOGICAL           UPPER
C     .. External Subroutines ..
      EXTERNAL          F06AAZ, F07ZAZ, F08FSY, F08FSZ, ZHER2K
C     .. Intrinsic Functions ..
      INTRINSIC         MAX
C     .. Executable Statements ..
C
C     Test the input parameters
C
      INFO = 0
      UPPER = (UPLO.EQ.'U' .OR. UPLO.EQ.'u')
      IF ( .NOT. UPPER .AND. .NOT. (UPLO.EQ.'L' .OR. UPLO.EQ.'l')) THEN
         INFO = -1
      ELSE IF (N.LT.0) THEN
         INFO = -2
      ELSE IF (LDA.LT.MAX(1,N)) THEN
         INFO = -4
      ELSE IF (LWORK.LT.1) THEN
         INFO = -9
      END IF
      IF (INFO.NE.0) THEN
         CALL F06AAZ('F08FSF/ZHETRD',-INFO)
         RETURN
      END IF
C
C     Quick return if possible
C
      IF (N.EQ.0) THEN
         WORK(1) = 1
         RETURN
      END IF
C
C     Determine the block size. Use LDA as the leading dimension of
C     WORK if there is room, otherwise use N.
C
      CALL F07ZAZ(1,'F08FSF',NB,0)
      NX = N
      IWS = 1
      IF (NB.GT.1 .AND. NB.LT.N) THEN
C
C        Determine when to cross over from blocked to unblocked code
C        (last block is always handled by unblocked code).
C
         CALL F07ZAZ(3,'F08FSF',NX,0)
         NX = MAX(NB,NX)
         IF (NX.LT.N) THEN
C
C           Determine if workspace is large enough for blocked code.
C
            LDWORK = N
            IWS = LDWORK*NB
            IF (LWORK.LT.IWS) THEN
C
C              Not enough workspace to use optimal NB:  determine the
C              minimum value of NB, and reduce NB or force use of
C              unblocked code by setting NX = N.
C
               NB = LWORK/LDWORK
               CALL F07ZAZ(2,'F08FSF',NBMIN,0)
               IF (NB.LT.NBMIN) NX = N
            END IF
         ELSE
            NX = N
         END IF
      ELSE
         NB = 1
      END IF
C
      IF (UPPER) THEN
C
C        Reduce the upper triangle of A.
C        Columns 1:kk are handled by the unblocked method.
C
         KK = N - ((N-NX+NB-1)/NB)*NB
         DO 40 I = N - NB + 1, KK + 1, -NB
C
C           Reduce columns i:i+nb-1 to tridiagonal form and form the
C           matrix W which is needed to update the unreduced part of
C           the matrix
C
            CALL F08FSY(UPLO,I+NB-1,NB,A,LDA,E,TAU,WORK,LDWORK)
C
C           Update the unreduced submatrix A(1:i-1,1:i-1), using an
C           update of the form:  A := A - V*W' - W*V'
C
            CALL ZHER2K(UPLO,'No transpose',I-1,NB,-CONE,A(1,I),LDA,
     *                  WORK,LDWORK,ONE,A,LDA)
C
C           Copy superdiagonal elements back into A, and diagonal
C           elements into D
C
            DO 20 J = I, I + NB - 1
               A(J-1,J) = E(J-1)
               D(J) = A(J,J)
   20       CONTINUE
   40    CONTINUE
C
C        Use unblocked code to reduce the last or only block
C
         CALL F08FSZ(UPLO,KK,A,LDA,D,E,TAU,IINFO)
      ELSE
C
C        Reduce the lower triangle of A
C
         DO 80 I = 1, N - NX, NB
C
C           Reduce columns i:i+nb-1 to tridiagonal form and form the
C           matrix W which is needed to update the unreduced part of
C           the matrix
C
            CALL F08FSY(UPLO,N-I+1,NB,A(I,I),LDA,E(I),TAU(I),WORK,
     *                  LDWORK)
C
C           Update the unreduced submatrix A(i+nb:n,i+nb:n), using
C           an update of the form:  A := A - V*W' - W*V'
C
            CALL ZHER2K(UPLO,'No transpose',N-I-NB+1,NB,-CONE,A(I+NB,I),
     *                  LDA,WORK(NB+1),LDWORK,ONE,A(I+NB,I+NB),LDA)
C
C           Copy subdiagonal elements back into A, and diagonal
C           elements into D
C
            DO 60 J = I, I + NB - 1
               A(J+1,J) = E(J)
               D(J) = A(J,J)
   60       CONTINUE
   80    CONTINUE
C
C        Use unblocked code to reduce the last or only block
C
         CALL F08FSZ(UPLO,N-I+1,A(I,I),LDA,D(I),E(I),TAU(I),IINFO)
      END IF
C
      WORK(1) = IWS
      RETURN
C
C     End of F08FSF (ZHETRD)
C
      END
*
      SUBROUTINE F08FSY(UPLO,N,NB,A,LDA,E,TAU,W,LDW)
C     MARK 16 RELEASE. NAG COPYRIGHT 1992.
C     ENTRY             ZLATRD(UPLO,N,NB,A,LDA,E,TAU,W,LDW)
C
C  Purpose
C  =======
C
C  ZLATRD reduces NB rows and columns of a complex Hermitian matrix A to
C  Hermitian tridiagonal form by a unitary similarity
C  transformation Q' * A * Q, and returns the matrices V and W which are
C  needed to apply the transformation to the unreduced part of A.
C
C  If UPLO = 'U', ZLATRD reduces the last NB rows and columns of a
C  matrix, of which the upper triangle is supplied;
C  if UPLO = 'L', ZLATRD reduces the first NB rows and columns of a
C  matrix, of which the lower triangle is supplied.
C
C  This is an auxiliary routine called by ZHETRD.
C
C  Arguments
C  =========
C
C  UPLO    (input) CHARACTER
C          Specifies whether the upper or lower triangular part of the
C          Hermitian matrix A is stored:
C          = 'U': Upper triangular
C          = 'L': Lower triangular
C
C  N       (input) INTEGER
C          The order of the matrix A.
C
C  NB      (input) INTEGER
C          The number of rows and columns to be reduced.
C
C  A       (input/output) COMPLEX*16 array, dimension (LDA,N)
C          On entry, the Hermitian matrix A.  If UPLO = 'U', the leading
C          n-by-n upper triangular part of A contains the upper
C          triangular part of the matrix A, and the strictly lower
C          triangular part of A is not referenced.  If UPLO = 'L', the
C          leading n-by-n lower triangular part of A contains the lower
C          triangular part of the matrix A, and the strictly upper
C          triangular part of A is not referenced.
C          On exit:
C          if UPLO = 'U', the last NB columns have been reduced to
C            tridiagonal form, with the diagonal elements overwriting
C            the diagonal elements of A; the elements above the diagonal
C            with the array TAU, represent the unitary matrix Q as a
C            product of elementary reflectors;
C          if UPLO = 'L', the first NB columns have been reduced to
C            tridiagonal form, with the diagonal elements overwriting
C            the diagonal elements of A; the elements below the diagonal
C            with the array TAU, represent the  unitary matrix Q as a
C            product of elementary reflectors.
C          See Further Details.
C
C  LDA     (input) INTEGER
C          The leading dimension of the array A.  LDA >= max(1,N).
C
C  E       (output) DOUBLE PRECISION array, dimension (N-1)
C          If UPLO = 'U', E(n-nb:n-1) contains the superdiagonal
C          elements of the last NB columns of the reduced matrix;
C          if UPLO = 'L', E(1:nb) contains the subdiagonal elements of
C          the first NB columns of the reduced matrix.
C
C  TAU     (output) COMPLEX*16 array, dimension (N-1)
C          The scalar factors of the elementary reflectors, stored in
C          TAU(n-nb:n-1) if UPLO = 'U', and in TAU(1:nb) if UPLO = 'L'.
C          See Further Details.
C
C  W       (output) COMPLEX*16 array, dimension (LDW,NB)
C          The n-by-nb matrix W required to update the unreduced part
C          of A.
C
C  LDW     (input) INTEGER
C          The leading dimension of the array W. LDW >= max(1,N).
C
C  Further Details
C  ===============
C
C  If UPLO = 'U', the matrix Q is represented as a product of elementary
C  reflectors
C
C     Q = H(n) H(n-1) . . . H(n-nb+1).
C
C  Each H(i) has the form
C
C     H(i) = I - tau * v * v'
C
C  where tau is a complex scalar, and v is a complex vector with
C  v(i:n) = 0 and v(i-1) = 1; v(1:i-1) is stored on exit in A(1:i-1,i),
C  and tau in TAU(i-1).
C
C  If UPLO = 'L', the matrix Q is represented as a product of elementary
C  reflectors
C
C     Q = H(1) H(2) . . . H(nb).
C
C  Each H(i) has the form
C
C     H(i) = I - tau * v * v'
C
C  where tau is a complex scalar, and v is a complex vector with
C  v(1:i) = 0 and v(i+1) = 1; v(i+1:n) is stored on exit in A(i+1:n,i),
C  and tau in TAU(i).
C
C  The elements of the vectors v together form the n-by-nb matrix V
C  which is needed, with W, to apply the transformation to the unreduced
C  part of the matrix, using a Hermitian rank-2k update of the form:
C  A := A - V*W' - W*V'.
C
C  The contents of A on exit are illustrated by the following examples
C  with n = 5 and nb = 2:
C
C  if UPLO = 'U':                       if UPLO = 'L':
C
C    (  a   a   a   v4  v5 )              (  d                  )
C    (      a   a   v4  v5 )              (  1   d              )
C    (          a   1   v5 )              (  v1  1   a          )
C    (              d   1  )              (  v1  v2  a   a      )
C    (                  d  )              (  v1  v2  a   a   a  )
C
C  where d denotes a diagonal element of the reduced matrix, a denotes
C  an element of the original matrix that is unchanged, and vi denotes
C  an element of the vector defining H(i).
C
C  -- LAPACK auxiliary routine (adapted for NAG Library)
C     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
C     Courant Institute, Argonne National Lab, and Rice University
C
C  =====================================================================
C
C     .. Parameters ..
      COMPLEX*16        ZERO, ONE, HALF
      PARAMETER         (ZERO=0.0D+0,ONE=1.0D+0,HALF=0.5D+0)
C     .. Scalar Arguments ..
      INTEGER           LDA, LDW, N, NB
      CHARACTER         UPLO
C     .. Array Arguments ..
      COMPLEX*16        A(LDA,*), TAU(*), W(LDW,*)
      DOUBLE PRECISION  E(*)
C     .. Local Scalars ..
      COMPLEX*16        ALPHA
      INTEGER           I, IW
C     .. External Functions ..
      COMPLEX*16        ZDOTC
      EXTERNAL          ZDOTC
C     .. External Subroutines ..
      EXTERNAL          F07FRY, F08ASV, ZAXPY, ZGEMV, ZHEMV, ZSCAL
C     .. Intrinsic Functions ..
      INTRINSIC         DBLE, MIN
C     .. Executable Statements ..
C
C     Quick return if possible
C
      IF (N.LE.0) RETURN
C
      IF ((UPLO.EQ.'U' .OR. UPLO.EQ.'u')) THEN
C
C        Reduce last NB columns of upper triangle
C
         DO 20 I = N, N - NB + 1, -1
            IW = I - N + NB
            IF (I.LT.N) THEN
C
C              Update A(1:i,i)
C
               A(I,I) = DBLE(A(I,I))
               CALL F07FRY(N-I,W(I,IW+1),LDW)
               CALL ZGEMV('No transpose',I,N-I,-ONE,A(1,I+1),LDA,
     *                    W(I,IW+1),LDW,ONE,A(1,I),1)
               CALL F07FRY(N-I,W(I,IW+1),LDW)
               CALL F07FRY(N-I,A(I,I+1),LDA)
               CALL ZGEMV('No transpose',I,N-I,-ONE,W(1,IW+1),LDW,
     *                    A(I,I+1),LDA,ONE,A(1,I),1)
               CALL F07FRY(N-I,A(I,I+1),LDA)
               A(I,I) = DBLE(A(I,I))
            END IF
            IF (I.GT.1) THEN
C
C              Generate elementary reflector H(i) to annihilate
C              A(1:i-2,i)
C
               ALPHA = A(I-1,I)
               CALL F08ASV(I-1,ALPHA,A(1,I),1,TAU(I-1))
               E(I-1) = ALPHA
               A(I-1,I) = ONE
C
C              Compute W(1:i-1,i)
C
               CALL ZHEMV('Upper',I-1,ONE,A,LDA,A(1,I),1,ZERO,W(1,IW),1)
               IF (I.LT.N) THEN
                  CALL ZGEMV('Conjugate transpose',I-1,N-I,ONE,W(1,IW+1)
     *                       ,LDW,A(1,I),1,ZERO,W(I+1,IW),1)
                  CALL ZGEMV('No transpose',I-1,N-I,-ONE,A(1,I+1),LDA,
     *                       W(I+1,IW),1,ONE,W(1,IW),1)
                  CALL ZGEMV('Conjugate transpose',I-1,N-I,ONE,A(1,I+1),
     *                       LDA,A(1,I),1,ZERO,W(I+1,IW),1)
                  CALL ZGEMV('No transpose',I-1,N-I,-ONE,W(1,IW+1),LDW,
     *                       W(I+1,IW),1,ONE,W(1,IW),1)
               END IF
               CALL ZSCAL(I-1,TAU(I-1),W(1,IW),1)
               ALPHA = -HALF*TAU(I-1)*ZDOTC(I-1,W(1,IW),1,A(1,I),1)
               CALL ZAXPY(I-1,ALPHA,A(1,I),1,W(1,IW),1)
            END IF
C
   20    CONTINUE
      ELSE
C
C        Reduce first NB columns of lower triangle
C
         DO 40 I = 1, NB
C
C           Update A(i:n,i)
C
            A(I,I) = DBLE(A(I,I))
            CALL F07FRY(I-1,W(I,1),LDW)
            CALL ZGEMV('No transpose',N-I+1,I-1,-ONE,A(I,1),LDA,W(I,1),
     *                 LDW,ONE,A(I,I),1)
            CALL F07FRY(I-1,W(I,1),LDW)
            CALL F07FRY(I-1,A(I,1),LDA)
            CALL ZGEMV('No transpose',N-I+1,I-1,-ONE,W(I,1),LDW,A(I,1),
     *                 LDA,ONE,A(I,I),1)
            CALL F07FRY(I-1,A(I,1),LDA)
            A(I,I) = DBLE(A(I,I))
            IF (I.LT.N) THEN
C
C              Generate elementary reflector H(i) to annihilate
C              A(i+2:n,i)
C
               ALPHA = A(I+1,I)
               CALL F08ASV(N-I,ALPHA,A(MIN(I+2,N),I),1,TAU(I))
               E(I) = ALPHA
               A(I+1,I) = ONE
C
C              Compute W(i+1:n,i)
C
               CALL ZHEMV('Lower',N-I,ONE,A(I+1,I+1),LDA,A(I+1,I),1,
     *                    ZERO,W(I+1,I),1)
               CALL ZGEMV('Conjugate transpose',N-I,I-1,ONE,W(I+1,1),
     *                    LDW,A(I+1,I),1,ZERO,W(1,I),1)
               CALL ZGEMV('No transpose',N-I,I-1,-ONE,A(I+1,1),LDA,
     *                    W(1,I),1,ONE,W(I+1,I),1)
               CALL ZGEMV('Conjugate transpose',N-I,I-1,ONE,A(I+1,1),
     *                    LDA,A(I+1,I),1,ZERO,W(1,I),1)
               CALL ZGEMV('No transpose',N-I,I-1,-ONE,W(I+1,1),LDW,
     *                    W(1,I),1,ONE,W(I+1,I),1)
               CALL ZSCAL(N-I,TAU(I),W(I+1,I),1)
               ALPHA = -HALF*TAU(I)*ZDOTC(N-I,W(I+1,I),1,A(I+1,I),1)
               CALL ZAXPY(N-I,ALPHA,A(I+1,I),1,W(I+1,I),1)
            END IF
C
   40    CONTINUE
      END IF
C
      RETURN
C
C     End of F08FSY (ZLATRD)
C
      END
*
      SUBROUTINE F08FSZ(UPLO,N,A,LDA,D,E,TAU,INFO)
C     MARK 16 RELEASE. NAG COPYRIGHT 1992.
C     ENTRY             ZHETD2(UPLO,N,A,LDA,D,E,TAU,INFO)
C
C  Purpose
C  =======
C
C  ZHETD2 reduces a complex Hermitian matrix A to real symmetric
C  tridiagonal form T by a unitary similarity transformation:
C  Q' * A * Q = T.
C
C  Arguments
C  =========
C
C  UPLO    (input) CHARACTER*1
C          Specifies whether the upper or lower triangular part of the
C          Hermitian matrix A is stored:
C          = 'U':  Upper triangular
C          = 'L':  Lower triangular
C
C  N       (input) INTEGER
C          The order of the matrix A.  N >= 0.
C
C  A       (input/output) COMPLEX*16 array, dimension (LDA,N)
C          On entry, the Hermitian matrix A.  If UPLO = 'U', the leading
C          n-by-n upper triangular part of A contains the upper
C          triangular part of the matrix A, and the strictly lower
C          triangular part of A is not referenced.  If UPLO = 'L', the
C          leading n-by-n lower triangular part of A contains the lower
C          triangular part of the matrix A, and the strictly upper
C          triangular part of A is not referenced.
C          On exit, if UPLO = 'U', the diagonal and first superdiagonal
C          of A are overwritten by the corresponding elements of the
C          tridiagonal matrix T, and the elements above the first
C          superdiagonal, with the array TAU, represent the unitary
C          matrix Q as a product of elementary reflectors; if UPLO
C          = 'L', the diagonal and first subdiagonal of A are over-
C          written by the corresponding elements of the tridiagonal
C          matrix T, and the elements below the first subdiagonal, with
C          the array TAU, represent the unitary matrix Q as a product
C          of elementary reflectors. See Further Details.
C
C  LDA     (input) INTEGER
C          The leading dimension of the array A.  LDA >= max(1,N).
C
C  D       (output) DOUBLE PRECISION array, dimension (N)
C          The diagonal elements of the tridiagonal matrix T:
C          D(i) = A(i,i).
C
C  E       (output) DOUBLE PRECISION array, dimension (N-1)
C          The off-diagonal elements of the tridiagonal matrix T:
C          E(i) = A(i,i+1) if UPLO = 'U', E(i) = A(i+1,i) if UPLO = 'L'.
C
C  TAU     (output) COMPLEX*16 array, dimension (N-1)
C          The scalar factors of the elementary reflectors (see Further
C          Details).
C
C  INFO    (output) INTEGER
C          = 0:  successful exit
C          < 0:  if INFO = -i, the i-th argument had an illegal value.
C
C  Further Details
C  ===============
C
C  If UPLO = 'U', the matrix Q is represented as a product of elementary
C  reflectors
C
C     Q = H(n-1) . . . H(2) H(1).
C
C  Each H(i) has the form
C
C     H(i) = I - tau * v * v'
C
C  where tau is a complex scalar, and v is a complex vector with
C  v(i+1:n) = 0 and v(i) = 1; v(1:i-1) is stored on exit in
C  A(1:i-1,i+1), and tau in TAU(i).
C
C  If UPLO = 'L', the matrix Q is represented as a product of elementary
C  reflectors
C
C     Q = H(1) H(2) . . . H(n-1).
C
C  Each H(i) has the form
C
C     H(i) = I - tau * v * v'
C
C  where tau is a complex scalar, and v is a complex vector with
C  v(1:i) = 0 and v(i+1) = 1; v(i+2:n) is stored on exit in A(i+2:n,i),
C  and tau in TAU(i).
C
C  The contents of A on exit are illustrated by the following examples
C  with n = 5:
C
C  if UPLO = 'U':                       if UPLO = 'L':
C
C    (  d   e   v2  v3  v4 )              (  d                  )
C    (      d   e   v3  v4 )              (  e   d              )
C    (          d   e   v4 )              (  v1  e   d          )
C    (              d   e  )              (  v1  v2  e   d      )
C    (                  d  )              (  v1  v2  v3  e   d  )
C
C  where d and e denote diagonal and off-diagonal elements of T, and vi
C  denotes an element of the vector defining H(i).
C
C  -- LAPACK routine (adapted for NAG Library)
C     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
C     Courant Institute, Argonne National Lab, and Rice University
C
C  =====================================================================
C
C     .. Parameters ..
      COMPLEX*16        ONE, ZERO, HALF
      PARAMETER         (ONE=1.0D0,ZERO=0.0D0,HALF=1.0D0/2.0D0)
C     .. Scalar Arguments ..
      INTEGER           INFO, LDA, N
      CHARACTER         UPLO
C     .. Array Arguments ..
      COMPLEX*16        A(LDA,*), TAU(*)
      DOUBLE PRECISION  D(*), E(*)
C     .. Local Scalars ..
      COMPLEX*16        ALPHA, TAUI
      INTEGER           I
      LOGICAL           UPPER
C     .. External Functions ..
      COMPLEX*16        ZDOTC
      EXTERNAL          ZDOTC
C     .. External Subroutines ..
      EXTERNAL          F06AAZ, F08ASV, ZAXPY, ZHEMV, ZHER2
C     .. Intrinsic Functions ..
      INTRINSIC         DBLE, MAX, MIN
C     .. Executable Statements ..
C
C     Test the input parameters
C
      INFO = 0
      UPPER = (UPLO.EQ.'U' .OR. UPLO.EQ.'u')
      IF ( .NOT. UPPER .AND. .NOT. (UPLO.EQ.'L' .OR. UPLO.EQ.'l')) THEN
         INFO = -1
      ELSE IF (N.LT.0) THEN
         INFO = -2
      ELSE IF (LDA.LT.MAX(1,N)) THEN
         INFO = -4
      END IF
      IF (INFO.NE.0) THEN
         CALL F06AAZ('F08FSZ/ZHETD2',-INFO)
         RETURN
      END IF
C
C     Quick return if possible
C
      IF (N.LE.0) RETURN
C
      IF (UPPER) THEN
C
C        Reduce the upper triangle of A
C
         A(N,N) = DBLE(A(N,N))
         DO 20 I = N - 1, 1, -1
C
C           Generate elementary reflector H(i) = I - tau * v * v'
C           to annihilate A(1:i-1,i+1)
C
            ALPHA = A(I,I+1)
            CALL F08ASV(I,ALPHA,A(1,I+1),1,TAUI)
            E(I) = ALPHA
C
            IF (TAUI.NE.ZERO) THEN
C
C              Apply H(i) from both sides to A(1:i,1:i)
C
               A(I,I+1) = ONE
C
C              Compute  x := tau * A * v  storing x in TAU(1:i)
C
               CALL ZHEMV(UPLO,I,TAUI,A,LDA,A(1,I+1),1,ZERO,TAU,1)
C
C              Compute  w := x - 1/2 * tau * (x'*v) * v
C
               ALPHA = -HALF*TAUI*ZDOTC(I,TAU,1,A(1,I+1),1)
               CALL ZAXPY(I,ALPHA,A(1,I+1),1,TAU,1)
C
C              Apply the transformation as a rank-2 update:
C                 A := A - v * w' - w * v'
C
               CALL ZHER2(UPLO,I,-ONE,A(1,I+1),1,TAU,1,A,LDA)
C
            END IF
            A(I,I+1) = E(I)
            D(I+1) = A(I+1,I+1)
            TAU(I) = TAUI
   20    CONTINUE
         D(1) = A(1,1)
      ELSE
C
C        Reduce the lower triangle of A
C
         A(1,1) = DBLE(A(1,1))
         DO 40 I = 1, N - 1
C
C           Generate elementary reflector H(i) = I - tau * v * v'
C           to annihilate A(i+2:n,i)
C
            ALPHA = A(I+1,I)
            CALL F08ASV(N-I,ALPHA,A(MIN(I+2,N),I),1,TAUI)
            E(I) = ALPHA
C
            IF (TAUI.NE.ZERO) THEN
C
C              Apply H(i) from both sides to A(i+1:n,i+1:n)
C
               A(I+1,I) = ONE
C
C              Compute  x := tau * A * v  storing y in TAU(i:n-1)
C
               CALL ZHEMV(UPLO,N-I,TAUI,A(I+1,I+1),LDA,A(I+1,I),1,ZERO,
     *                    TAU(I),1)
C
C              Compute  w := x - 1/2 * tau * (x'*v) * v
C
               ALPHA = -HALF*TAUI*ZDOTC(N-I,TAU(I),1,A(I+1,I),1)
               CALL ZAXPY(N-I,ALPHA,A(I+1,I),1,TAU(I),1)
C
C              Apply the transformation as a rank-2 update:
C                 A := A - v * w' - w * v'
C
               CALL ZHER2(UPLO,N-I,-ONE,A(I+1,I),1,TAU(I),1,A(I+1,I+1),
     *                    LDA)
C
            END IF
            A(I+1,I) = E(I)
            D(I) = A(I,I)
            TAU(I) = TAUI
   40    CONTINUE
         D(N) = A(N,N)
      END IF
C
      RETURN
C
C     End of F08FSZ (ZHETD2)
C
      END
*
      SUBROUTINE F08FTF(UPLO,N,A,LDA,TAU,WORK,LWORK,INFO)
C     MARK 16 RELEASE. NAG COPYRIGHT 1992.
C     .. Entry Points ..
      ENTRY             ZUNGTR(UPLO,N,A,LDA,TAU,WORK,LWORK,INFO)
C
C  Purpose
C  =======
C
C  ZUNGTR generates a complex unitary matrix Q which is defined as the
C  product of n-1 elementary reflectors of order n, as returned by
C  ZHETRD:
C
C  if UPLO = 'U', Q = H(n-1) . . . H(2) H(1),
C
C  if UPLO = 'L', Q = H(1) H(2) . . . H(n-1).
C
C  Arguments
C  =========
C
C  UPLO    (input) CHARACTER*1
C          Specifies whether the upper or lower triangle of the array A
C          holds details of the elementary reflectors, as returned by
C          ZHETRD:
C          = 'U': Upper triangle;
C          = 'L': Lower triangle.
C
C  N       (input) INTEGER
C          The order of the matrix Q. N >= 0.
C
C  A       (input/output) COMPLEX*16 array, dimension (LDA,N)
C          On entry, the vectors which define the elementary reflectors,
C          as returned by ZHETRD.
C          On exit, the n by n unitary matrix Q.
C
C  LDA     (input) INTEGER
C          The leading dimension of the array A. LDA >= N.
C
C  TAU     (input) COMPLEX*16 array, dimension (N-1)
C          TAU(i) must contain the scalar factor of the elementary
C          reflector H(i), as returned by ZHETRD.
C
C  WORK    (workspace) COMPLEX*16 array, dimension (LWORK)
C          On exit, if INFO = 0, WORK(1) returns the minimum value of
C          LWORK required to use the optimal blocksize.
C
C  LWORK   (input) INTEGER
C          The dimension of the array WORK. LWORK >= N-1.
C          For optimum performance LWORK should be at least (N-1)*NB,
C          where NB is the optimal blocksize.
C
C  INFO    (output) INTEGER
C          = 0: successful exit
C          < 0: if INFO = -i, the i-th argument had an illegal value
C
C  -- LAPACK routine (adapted for NAG Library)
C     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
C     Courant Institute, Argonne National Lab, and Rice University
C
C  =====================================================================
C
C     .. Parameters ..
      COMPLEX*16        ZERO, ONE
      PARAMETER         (ZERO=0.0D+0,ONE=1.0D+0)
C     .. Scalar Arguments ..
      INTEGER           INFO, LDA, LWORK, N
      CHARACTER         UPLO
C     .. Array Arguments ..
      COMPLEX*16        A(LDA,*), TAU(*), WORK(LWORK)
C     .. Local Scalars ..
      INTEGER           I, IINFO, J
      LOGICAL           UPPER
C     .. External Subroutines ..
      EXTERNAL          F06AAZ, F08CTY, ZUNGQR
C     .. Intrinsic Functions ..
      INTRINSIC         MAX
C     .. Executable Statements ..
C
C     Test the input arguments
C
      INFO = 0
      UPPER = (UPLO.EQ.'U' .OR. UPLO.EQ.'u')
      IF ( .NOT. UPPER .AND. .NOT. (UPLO.EQ.'L' .OR. UPLO.EQ.'l')) THEN
         INFO = -1
      ELSE IF (N.LT.0) THEN
         INFO = -2
      ELSE IF (LDA.LT.MAX(1,N)) THEN
         INFO = -4
      ELSE IF (LWORK.LT.MAX(1,N-1)) THEN
         INFO = -7
      END IF
      IF (INFO.NE.0) THEN
         CALL F06AAZ('F08FTF/ZUNGTR',-INFO)
         RETURN
      END IF
C
C     Quick return if possible
C
      IF (N.EQ.0) THEN
         WORK(1) = 1
         RETURN
      END IF
C
      IF (UPPER) THEN
C
C        Q was determined by a call to ZHETRD with UPLO = 'U'
C
C        Shift the vectors which define the elementary reflectors one
C        column to the left, and set the last row and column of Q to
C        those of the unit matrix
C
         DO 40 J = 1, N - 1
            DO 20 I = 1, J - 1
               A(I,J) = A(I,J+1)
   20       CONTINUE
            A(N,J) = ZERO
   40    CONTINUE
         DO 60 I = 1, N - 1
            A(I,N) = ZERO
   60    CONTINUE
         A(N,N) = ONE
C
C        Generate Q(1:n-1,1:n-1)
C
         CALL F08CTY(N-1,N-1,N-1,A,LDA,TAU,WORK,LWORK,IINFO)
C
      ELSE
C
C        Q was determined by a call to ZHETRD with UPLO = 'L'.
C
C        Shift the vectors which define the elementary reflectors one
C        column to the right, and set the first row and column of Q to
C        those of the unit matrix
C
         DO 100 J = N, 2, -1
            A(1,J) = ZERO
            DO 80 I = J + 1, N
               A(I,J) = A(I,J-1)
   80       CONTINUE
  100    CONTINUE
         A(1,1) = ONE
         DO 120 I = 2, N
            A(I,1) = ZERO
  120    CONTINUE
         IF (N.GT.1) THEN
C
C           Generate Q(2:n,2:n)
C
            CALL ZUNGQR(N-1,N-1,N-1,A(2,2),LDA,TAU,WORK,LWORK,IINFO)
         ELSE
            WORK(1) = 1
         END IF
      END IF
      RETURN
C
C     End of F08FTF (ZUNGTR)
C
      END
*
      SUBROUTINE F08HEW(F,G,CS,SN,R)
C     MARK 16 RELEASE. NAG COPYRIGHT 1992.
C     ENTRY             DLARTG(F,G,CS,SN,R)
C
C  Purpose
C  =======
C
C  DLARTG generate a plane rotation so that
C
C     [  CS  SN  ]  .  [ F ]  =  [ R ]   where CS**2 + SN**2 = 1.
C     [ -SN  CS  ]     [ G ]     [ 0 ]
C
C  This is a faster version of the BLAS1 routine DROTG, except for
C  the following differences:
C     F and G are unchanged on return.
C     If G=0, then CS=1 and SN=0.
C     If F=0 and (G .ne. 0), then CS=0 and SN=1 without doing any
C        floating point operations (saves work in DBDSQR when
C        there are zeros on the diagonal).
C
C  Arguments
C  =========
C
C  F       (input) DOUBLE PRECISION
C          The first component of vector to be rotated.
C
C  G       (input) DOUBLE PRECISION
C          The second component of vector to be rotated.
C
C  CS      (output) DOUBLE PRECISION
C          The cosine of the rotation.
C
C  SN      (output) DOUBLE PRECISION
C          The sine of the rotation.
C
C  R       (output) DOUBLE PRECISION
C          The nonzero component of the rotated vector.
C
C  -- LAPACK auxiliary routine (adapted for NAG Library)
C     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
C     Courant Institute, Argonne National Lab, and Rice University
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO
      PARAMETER         (ZERO=0.0D0)
      DOUBLE PRECISION  ONE
      PARAMETER         (ONE=1.0D0)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  CS, F, G, R, SN
C     .. Local Scalars ..
      DOUBLE PRECISION  T, TT
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, SQRT
C     .. Executable Statements ..
C
      IF (G.EQ.ZERO) THEN
         CS = ONE
         SN = ZERO
         R = F
      ELSE IF (F.EQ.ZERO) THEN
         CS = ZERO
         SN = ONE
         R = G
      ELSE
         IF (ABS(F).GT.ABS(G)) THEN
            T = G/F
            TT = SQRT(ONE+T*T)
            CS = ONE/TT
            SN = T*CS
            R = F*TT
         ELSE
            T = F/G
            TT = SQRT(ONE+T*T)
            SN = ONE/TT
            CS = T*SN
            R = G*TT
         END IF
      END IF
      RETURN
C
C     End of F08HEW (DLARTG)
C
      END
*
      SUBROUTINE F08JEX(A,B,C,RT1,RT2)
C     MARK 16 RELEASE. NAG COPYRIGHT 1992.
C     ENTRY             DLAE2(A,B,C,RT1,RT2)
C
C  Purpose
C  =======
C
C  DLAE2  computes the eigenvalues of a 2-by-2 symmetric matrix
C     [  A   B  ]
C     [  B   C  ].
C  On return, RT1 is the eigenvalue of larger absolute value, and RT2
C  is the eigenvalue of smaller absolute value.
C
C  Arguments
C  =========
C
C  A       (input) DOUBLE PRECISION
C          The (1,1) entry of the 2-by-2 matrix.
C
C  B       (input) DOUBLE PRECISION
C          The (1,2) and (2,1) entries of the 2-by-2 matrix.
C
C  C       (input) DOUBLE PRECISION
C          The (2,2) entry of the 2-by-2 matrix.
C
C  RT1     (output) DOUBLE PRECISION
C          The eigenvalue of larger absolute value.
C
C  RT2     (output) DOUBLE PRECISION
C          The eigenvalue of smaller absolute value.
C
C  Further Details
C  ===============
C
C  RT1 is accurate to a few ulps barring over/underflow.
C
C  RT2 may be inaccurate if there is massive cancellation in the
C  determinant A*C-B*B; higher precision or correctly rounded or
C  correctly truncated arithmetic would be needed to compute RT2
C  accurately in all cases.
C
C  Overflow is possible only if RT1 is within a factor of 5 of overflow.
C  Underflow is harmless if the input data is 0 or exceeds
C     underflow_threshold / macheps.
C
C
C  -- LAPACK auxiliary routine (adapted for NAG Library)
C     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
C     Courant Institute, Argonne National Lab, and Rice University
C
C     .. Parameters ..
      DOUBLE PRECISION  ONE
      PARAMETER         (ONE=1.0D0)
      DOUBLE PRECISION  TWO
      PARAMETER         (TWO=2.0D0)
      DOUBLE PRECISION  ZERO
      PARAMETER         (ZERO=0.0D0)
      DOUBLE PRECISION  HALF
      PARAMETER         (HALF=0.5D0)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  A, B, C, RT1, RT2
C     .. Local Scalars ..
      DOUBLE PRECISION  AB, ACMN, ACMX, ADF, DF, RT, SM, TB
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, SQRT
C     .. Executable Statements ..
C
C     Compute the eigenvalues
C
      SM = A + C
      DF = A - C
      ADF = ABS(DF)
      TB = B + B
      AB = ABS(TB)
      IF (ABS(A).GT.ABS(C)) THEN
         ACMX = A
         ACMN = C
      ELSE
         ACMX = C
         ACMN = A
      END IF
      IF (ADF.GT.AB) THEN
         RT = ADF*SQRT(ONE+(AB/ADF)**2)
      ELSE IF (ADF.LT.AB) THEN
         RT = AB*SQRT(ONE+(ADF/AB)**2)
      ELSE
C
C        Includes case AB=ADF=0
C
         RT = AB*SQRT(TWO)
      END IF
      IF (SM.LT.ZERO) THEN
         RT1 = HALF*(SM-RT)
C
C        Order of execution important.
C        To get fully accurate smaller eigenvalue,
C        next line needs to be executed in higher precision.
C
         RT2 = (ACMX/RT1)*ACMN - (B/RT1)*B
      ELSE IF (SM.GT.ZERO) THEN
         RT1 = HALF*(SM+RT)
C
C        Order of execution important.
C        To get fully accurate smaller eigenvalue,
C        next line needs to be executed in higher precision.
C
         RT2 = (ACMX/RT1)*ACMN - (B/RT1)*B
      ELSE
C
C        Includes case RT1 = RT2 = 0
C
         RT1 = HALF*RT
         RT2 = -HALF*RT
      END IF
      RETURN
C
C     End of F08JEX (DLAE2)
C
      END
*
      SUBROUTINE F08JFF(N,D,E,INFO)
C     MARK 16 RELEASE. NAG COPYRIGHT 1992.
C     .. Entry Points ..
      ENTRY             DSTERF(N,D,E,INFO)
C
C  Purpose
C  =======
C
C  DSTERF computes all the eigenvalues of a real symmetric tridiagonal
C  matrix T, using the Pal-Walker-Kahan root-free variants of the QL and
C  QR algorithms.
C
C  Arguments
C  =========
C
C  N       (input) INTEGER
C          The order of the matrix T.  N >= 0.
C
C  D       (input/output) DOUBLE PRECISION array, dimension (N)
C          On entry, the n diagonal elements of the tridiagonal matrix.
C          On exit, if INFO = 0, the eigenvalues in ascending order.
C
C  E       (input/output) DOUBLE PRECISION array, dimension (N-1)
C          On entry, the n-1 subdiagonal elements of the tridiagonal
C          matrix.
C          On exit, E has been overwritten.
C
C  INFO    (output) INTEGER
C          = 0: successful exit.
C          < 0: if INFO = -i, the i-th argument had an illegal value.
C          > 0: the algorithm has failed to find all the eigenvalues in
C               a total of 30*N iterations; if INFO = i, then i elements
C               of E have not converged to zero.
C
C  -- LAPACK routine (adapted for NAG Library)
C     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
C     Courant Institute, Argonne National Lab, and Rice University
C
C  =====================================================================
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE, TWO
      PARAMETER         (ZERO=0.0D0,ONE=1.0D0,TWO=2.0D0)
      INTEGER           MAXIT
      PARAMETER         (MAXIT=30)
C     .. Scalar Arguments ..
      INTEGER           INFO, N
C     .. Array Arguments ..
      DOUBLE PRECISION  D(*), E(*)
C     .. Local Scalars ..
      DOUBLE PRECISION  ALPHA, BB, C, EPS, GAMMA, OLDC, OLDGAM, P, R,
     *                  RT1, RT2, RTE, S, SIGMA
      INTEGER           I, ITN, J, K, L, LL, M, MM
C     .. External Functions ..
      DOUBLE PRECISION  F06BNF, X02AJF
      EXTERNAL          F06BNF, X02AJF
C     .. External Subroutines ..
      EXTERNAL          F06AAZ, F08JEX
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, SIGN, SQRT
C     .. Executable Statements ..
C
C     Test the input parameters.
C
      INFO = 0
      IF (N.LT.0) THEN
         INFO = -1
         CALL F06AAZ('F08JFF/DSTERF',-INFO)
         RETURN
      END IF
C
C     Quick return if possible
C
      IF (N.EQ.0) RETURN
C
C     Determine the unit roundoff for this environment.
C
      EPS = X02AJF()
C
C     Square the elements of E.
C
      DO 20 I = 1, N - 1
         E(I) = E(I)*E(I)
   20 CONTINUE
C
      ITN = N*MAXIT
      MM = 0
C
   40 CONTINUE
C
C     Determine where the matrix splits and choose QL or QR iteration
C     for each block, according to whether top or bottom diagonal
C     element is smaller.
C
      LL = MM + 1
      IF (LL.GT.N) GO TO 300
      IF (LL.GT.1) E(LL-1) = ZERO
      DO 60 MM = LL, N - 1
         IF (SQRT(E(MM)).LE.EPS*(ABS(D(MM))+ABS(D(MM+1)))) GO TO 80
   60 CONTINUE
   80 CONTINUE
C
      IF (ABS(D(LL)).LE.ABS(D(MM))) THEN
C
C        Perform QL iterations on rows and columns LL to MM;
C        unconverged eigenvalues are in rows and columns L to MM.
C
         L = LL
  100    CONTINUE
         IF (L.GT.MM) GO TO 40
C
C        Look for small offdiagonal element.
C
         DO 120 M = L, MM - 1
            IF (SQRT(E(M)).LE.EPS*(ABS(D(M))+ABS(D(M+1)))) GO TO 140
  120    CONTINUE
  140    CONTINUE
         IF (M.NE.N) E(M) = ZERO
C
         IF (M.GT.L+1) THEN
C
C           Perform QL iteration on rows and columns L to M.
C
            ITN = ITN - 1
            IF (ITN.LT.0) GO TO 260
C
C           Form shift.
C
            P = D(L)
            RTE = SQRT(E(L))
            SIGMA = (D(L+1)-P)/(TWO*RTE)
            R = F06BNF(SIGMA,ONE)
            SIGMA = P - (RTE/(SIGMA+SIGN(R,SIGMA)))
C
C           Inner loop.
C
            C = ONE
            S = ZERO
            GAMMA = D(M) - SIGMA
            P = GAMMA*GAMMA
            DO 160 I = M - 1, L, -1
               BB = E(I)
               R = P + BB
               IF (I.NE.M-1) E(I+1) = S*R
               OLDC = C
               C = P/R
               S = BB/R
               OLDGAM = GAMMA
               ALPHA = D(I)
               GAMMA = C*(ALPHA-SIGMA) - S*OLDGAM
               D(I+1) = OLDGAM + (ALPHA-GAMMA)
               IF (C.NE.ZERO) THEN
                  P = (GAMMA*GAMMA)/C
               ELSE
                  P = OLDC*BB
               END IF
  160       CONTINUE
            E(L) = S*P
            D(L) = SIGMA + GAMMA
C
         ELSE
            IF (M.EQ.L+1) THEN
C
C              If remaining matrix is 2 by 2, use F08JEX to compute its
C              eigensystem.
C
               CALL F08JEX(D(L),SQRT(E(L)),D(L+1),RT1,RT2)
               D(L) = RT1
               D(L+1) = RT2
               E(L) = ZERO
            END IF
            L = M + 1
         END IF
         GO TO 100
C
      ELSE
C
C        Perform QR iterations on rows and columns LL to MM;
C        unconverged eigenvalues are in rows and columns LL to M.
C
         M = MM
  180    CONTINUE
         IF (M.LT.LL) GO TO 40
C
C        Look for small offdiagonal element.
C
         DO 200 L = M, LL + 1, -1
            IF (SQRT(E(L-1)).LE.EPS*(ABS(D(L))+ABS(D(L-1)))) GO TO 220
  200    CONTINUE
  220    CONTINUE
         IF (L.NE.1) E(L-1) = ZERO
C
         IF (L.LT.M-1) THEN
C
C           Perform QR iteration on rows and columns L to M.
C
            ITN = ITN - 1
            IF (ITN.LT.0) GO TO 260
C
C           Form shift.
C
            P = D(M)
            RTE = SQRT(E(M-1))
            SIGMA = (D(M-1)-P)/(TWO*RTE)
            R = F06BNF(SIGMA,ONE)
            SIGMA = P - (RTE/(SIGMA+SIGN(R,SIGMA)))
C
C           Inner loop.
C
            C = ONE
            S = ZERO
            GAMMA = D(L) - SIGMA
            P = GAMMA*GAMMA
            DO 240 I = L, M - 1
               BB = E(I)
               R = P + BB
               IF (I.NE.L) E(I-1) = S*R
               OLDC = C
               C = P/R
               S = BB/R
               OLDGAM = GAMMA
               ALPHA = D(I+1)
               GAMMA = C*(ALPHA-SIGMA) - S*OLDGAM
               D(I) = OLDGAM + (ALPHA-GAMMA)
               IF (C.NE.ZERO) THEN
                  P = (GAMMA*GAMMA)/C
               ELSE
                  P = OLDC*BB
               END IF
  240       CONTINUE
            E(M-1) = S*P
            D(M) = SIGMA + GAMMA
C
         ELSE
            IF (L.EQ.M-1) THEN
C
C              If remaining matrix is 2 by 2, use F08JEX to compute its
C              eigenvalues.
C
               CALL F08JEX(D(M-1),SQRT(E(M-1)),D(M),RT1,RT2)
               D(M-1) = RT1
               D(M) = RT2
               E(M-1) = ZERO
            END IF
            M = L - 1
         END IF
         GO TO 180
C
      END IF
C
  260 CONTINUE
C
C     Failure to converge.
C
      DO 280 I = 1, N - 1
         IF (E(I).NE.ZERO) INFO = INFO + 1
  280 CONTINUE
      RETURN
C
  300 CONTINUE
C
C     Order eigenvalues.
C
      DO 340 I = 1, N - 1
         K = I
         P = D(I)
         DO 320 J = I + 1, N
            IF (D(J).LT.P) THEN
               K = J
               P = D(J)
            END IF
  320    CONTINUE
         IF (K.NE.I) THEN
            D(K) = D(I)
            D(I) = P
         END IF
  340 CONTINUE
C
      RETURN
C
C     End of F08JFF (DSTERF)
C
      END
*
      SUBROUTINE F08JSF(COMPZ,N,D,E,Z,LDZ,WORK,INFO)
C     MARK 16 RELEASE. NAG COPYRIGHT 1992.
C     .. Entry Points ..
      ENTRY             ZSTEQR(COMPZ,N,D,E,Z,LDZ,WORK,INFO)
C
C  Purpose
C  =======
C
C  ZSTEQR computes all eigenvalues and, optionally, eigenvectors of a
C  real symmetric tridiagonal matrix T, using a combination of the QL
C  and QR algorithms, with an implicit shift.
C
C  ZSTEQR can also compute the eigenvectors of a complex Hermitian
C  matrix A, which has been reduced to tridiagonal form by ZHETRD,
C  ZHPTRD or ZHBTRD.
C
C  Arguments
C  =========
C
C  COMPZ   (input) CHARACTER*1
C          = 'N': Compute eigenvalues only, no eigenvectors;
C          = 'I': Compute eigenvectors of the tridiagonal matrix T; Z is
C                 initialized by the routine;
C          = 'V': Compute eigenvectors of a symmetric matrix A which has
C                 been reduced to tridiagonal form T = Q'*A*Q; Z must
C                 contain Q on entry.
C
C  N       (input) INTEGER
C          The order of the matrix T.  N >= 0.
C
C  D       (input/output) DOUBLE PRECISION array, dimension (N)
C          On entry, the n diagonal elements of the tridiagonal matrix.
C          On exit, if INFO = 0, the eigenvalues in ascending order.
C
C  E       (input/output) DOUBLE PRECISION array, dimension (N-1)
C          On entry, the n-1 offdiagonal elements of the tridiagonal
C          matrix.
C          On exit, E has been overwritten.
C
C  Z       (input/output) COMPLEX*16 array, dimension (LDZ, N)
C          On entry, if COMPZ = 'V', Z must contain the unitary
C          matrix used in the reduction to tridiagonal form; if
C          COMPZ = 'N' or 'I', Z need not be set.
C          On exit, if INFO = 0, if COMPZ = 'I' or 'V', Z contains the
C          requested eigenvectors.
C          If COMPZ = 'N', Z is not referenced.
C
C  LDZ     (input) INTEGER
C          The leading dimension of the array Z.
C          If COMPZ = 'I' or 'V', LDZ >= max(1,N); otherwise LDZ >= 1.
C
C  WORK    (workspace) DOUBLE PRECISION array, dimension (max(1,2*N-2))
C          If COMPZ = 'N', WORK is not referenced.
C
C  INFO    (output) INTEGER
C          = 0: successful exit.
C          < 0: if INFO = -i, the i-th argument had an illegal value.
C          > 0: the algorithm has failed to find all the eigenvalues in
C               a total of 30*N iterations; if INFO = i, then i elements
C               of E have not converged to zero; on exit, D and E
C               represent a tridiagonal matrix which is unitarily
C               similar to the original matrix.
C
C  -- LAPACK routine (adapted for NAG Library)
C     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
C     Courant Institute, Argonne National Lab, and Rice University
C
C  =====================================================================
C
C     .. Parameters ..
      COMPLEX*16        CZERO, CONE
      PARAMETER         (CZERO=0.0D0,CONE=1.0D0)
      DOUBLE PRECISION  ZERO, ONE, TWO
      PARAMETER         (ZERO=0.0D0,ONE=1.0D0,TWO=2.0D0)
      INTEGER           MAXIT
      PARAMETER         (MAXIT=30)
C     .. Scalar Arguments ..
      INTEGER           INFO, LDZ, N
      CHARACTER         COMPZ
C     .. Array Arguments ..
      COMPLEX*16        Z(LDZ,*)
      DOUBLE PRECISION  D(*), E(*), WORK(*)
C     .. Local Scalars ..
      COMPLEX*16        TEMP
      DOUBLE PRECISION  B, C, EPS, F, G, P, R, RT1, RT2, S
      INTEGER           I, II, ITN, J, K, L, LL, M, MM
      LOGICAL           INITZ, WANTZ
C     .. External Functions ..
      DOUBLE PRECISION  F06BNF, X02AJF
      INTEGER           IDAMAX
      EXTERNAL          F06BNF, X02AJF, IDAMAX
C     .. External Subroutines ..
      EXTERNAL          F06AAZ, F06KPF, F06THF, F06VXF, F07MDX, F08HEW,
     *                  F08JEX, ZSCAL, ZSWAP
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, DCONJG, MAX, SIGN
C     .. Executable Statements ..
C
C     Decode and test the input parameters.
C
      INITZ = (COMPZ.EQ.'I' .OR. COMPZ.EQ.'i')
      WANTZ = INITZ .OR. (COMPZ.EQ.'V' .OR. COMPZ.EQ.'v')
C
      INFO = 0
      IF ( .NOT. (COMPZ.EQ.'N' .OR. COMPZ.EQ.'n') .AND. .NOT. WANTZ)
     *    THEN
         INFO = -1
      ELSE IF (N.LT.0) THEN
         INFO = -2
      ELSE IF ((LDZ.LT.1) .OR. (WANTZ .AND. LDZ.LT.MAX(1,N))) THEN
         INFO = -6
      END IF
      IF (INFO.NE.0) THEN
         CALL F06AAZ('F08JSF/ZSTEQR',-INFO)
         RETURN
      END IF
C
C     Quick return if possible
C
      IF (N.EQ.0) RETURN
C
C     Determine the unit roundoff for this environment.
C
      EPS = X02AJF()
C
C     If COMPZ = 'I', initialize Z to the unit matrix.
C
      IF (INITZ) CALL F06THF('General',N,N,CZERO,CONE,Z,LDZ)
C
      ITN = N*MAXIT
      MM = 0
C
   20 CONTINUE
      LL = MM + 1
      IF (LL.GT.N) GO TO 280
      IF (LL.GT.1) E(LL-1) = ZERO
C
C     Determine where the matrix splits and choose QL or QR iteration
C     for each block, according to whether top or bottom diagonal
C     element is smaller.
C
      DO 40 MM = LL, N - 1
         IF (ABS(E(MM)).LE.EPS*(ABS(D(MM))+ABS(D(MM+1)))) GO TO 60
   40 CONTINUE
   60 CONTINUE
C
      IF (ABS(D(LL)).LE.ABS(D(MM))) THEN
C
C        Perform QL iterations on rows and columns LL to MM;
C        unconverged eigenvalues are in rows and columns L to MM.
C
         L = LL
   80    CONTINUE
         IF (L.GT.MM) GO TO 20
C
C        Look for small offdiagonal element.
C
         DO 100 M = L, MM - 1
            IF (ABS(E(M)).LE.EPS*(ABS(D(M))+ABS(D(M+1)))) GO TO 120
  100    CONTINUE
  120    CONTINUE
         IF (M.NE.N) E(M) = ZERO
C
         IF (M.GT.L+1) THEN
C
C           Perform QL iteration on rows and columns L to M.
C
            ITN = ITN - 1
            IF (ITN.LT.0) GO TO 240
C
C           Form shift.
C
            P = D(L)
            G = (D(L+1)-P)/(TWO*E(L))
            R = F06BNF(G,ONE)
            G = D(M) - P + (E(L)/(G+SIGN(R,G)))
C
C           Inner loop.
C
            C = ONE
            S = ONE
            P = ZERO
            DO 140 I = M - 1, L, -1
               F = S*E(I)
               B = C*E(I)
               CALL F08HEW(G,F,C,S,R)
               IF (I.NE.M-1) E(I+1) = R
               G = D(I+1) - P
               R = (D(I)-G)*S + TWO*C*B
               P = S*R
               D(I+1) = G + P
               G = C*R - B
               IF (WANTZ) THEN
C
C                 Save parameters of rotation.
C
                  WORK(I) = C
                  WORK(N-1+I) = -S
               END IF
  140       CONTINUE
            E(L) = G
            D(L) = D(L) - P
            IF (WANTZ) THEN
C
C              Apply saved rotations.
C
               CALL F06VXF('R','V','B',N,M-L+1,1,M-L+1,WORK(L),
     *                     WORK(N-1+L),Z(1,L),LDZ)
            END IF
         ELSE
            IF (M.EQ.L+1) THEN
C
C              If remaining matrix is 2 by 2, use F08JEX or F07MDX
C              to compute its eigensystem.
C
               IF (WANTZ) THEN
                  CALL F07MDX(D(L),E(L),D(L+1),RT1,RT2,C,S)
                  CALL F06KPF(N,Z(1,L),1,Z(1,L+1),1,C,S)
               ELSE
                  CALL F08JEX(D(L),E(L),D(L+1),RT1,RT2)
               END IF
               D(L) = RT1
               D(L+1) = RT2
               E(L) = ZERO
            END IF
            L = M + 1
         END IF
         GO TO 80
C
      ELSE
C
C        Perform QR iterations on rows and columns LL to MM;
C        unconverged eigenvalues are in rows and columns LL to M.
C
         M = MM
  160    CONTINUE
         IF (M.LT.LL) GO TO 20
C
C        Look for small offdiagonal element.
C
         DO 180 L = M, LL + 1, -1
            IF (ABS(E(L-1)).LE.EPS*(ABS(D(L))+ABS(D(L-1)))) GO TO 200
  180    CONTINUE
  200    CONTINUE
         IF (L.NE.1) E(L-1) = ZERO
C
         IF (L.LT.M-1) THEN
C
C           Perform QR iteration on rows and columns L to M.
C
            ITN = ITN - 1
            IF (ITN.LT.0) GO TO 240
C
C           Form shift.
C
            P = D(M)
            G = (D(M-1)-P)/(TWO*E(M-1))
            R = F06BNF(G,ONE)
            G = D(L) - P + (E(M-1)/(G+SIGN(R,G)))
C
C           Inner loop.
C
            C = ONE
            S = ONE
            P = ZERO
            DO 220 I = L, M - 1
               F = S*E(I)
               B = C*E(I)
               CALL F08HEW(G,F,C,S,R)
               IF (I.NE.L) E(I-1) = R
               G = D(I) - P
               R = (D(I+1)-G)*S + TWO*C*B
               P = S*R
               D(I) = G + P
               G = C*R - B
               IF (WANTZ) THEN
C
C                 Save parameters of rotation.
C
                  WORK(I) = C
                  WORK(N-1+I) = S
               END IF
  220       CONTINUE
            E(M-1) = G
            D(M) = D(M) - P
            IF (WANTZ) THEN
C
C              Apply saved rotations.
C
               CALL F06VXF('R','V','F',N,M-L+1,1,M-L+1,WORK(L),
     *                     WORK(N-1+L),Z(1,L),LDZ)
            END IF
         ELSE
            IF (L.EQ.M-1) THEN
C
C              If remaining matrix is 2 by 2, use F08JEX or F07MDX
C              to compute its eigensystem.
C
               IF (WANTZ) THEN
                  CALL F07MDX(D(M-1),E(M-1),D(M),RT1,RT2,C,S)
                  CALL F06KPF(N,Z(1,M-1),1,Z(1,M),1,C,S)
               ELSE
                  CALL F08JEX(D(M-1),E(M-1),D(M),RT1,RT2)
               END IF
               D(M-1) = RT1
               D(M) = RT2
               E(M-1) = ZERO
            END IF
            M = L - 1
         END IF
         GO TO 160
C
      END IF
C
  240 CONTINUE
C
C     Failure to converge.
C
      DO 260 I = 1, N - 1
         IF (E(I).NE.ZERO) INFO = INFO + 1
  260 CONTINUE
      RETURN
C
  280 CONTINUE
C
C     Order eigenvalues and eigenvectors.
C
      DO 320 I = 1, N - 1
         K = I
         P = D(I)
         DO 300 J = I + 1, N
            IF (D(J).LT.P) THEN
               K = J
               P = D(J)
            END IF
  300    CONTINUE
         IF (K.NE.I) THEN
            D(K) = D(I)
            D(I) = P
            IF (WANTZ) CALL ZSWAP(N,Z(1,I),1,Z(1,K),1)
         END IF
  320 CONTINUE
C
      IF (WANTZ) THEN
C
C        Normalize eigenvectors so that element of largest absolute
C        value is real and positive
C
         DO 360 J = 1, N
            DO 340 I = 1, N
               WORK(I) = ABS(Z(I,J))
  340       CONTINUE
            II = IDAMAX(N,WORK,1)
            TEMP = Z(II,J)/WORK(II)
            CALL ZSCAL(N,DCONJG(TEMP),Z(1,J),1)
            Z(II,J) = WORK(II)
  360    CONTINUE
      END IF
      RETURN
C
C     End of F08JSF (ZSTEQR)
C
      END
*
      DOUBLE PRECISION FUNCTION F06BMF( SCALE, SSQ )
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     .. Scalar Arguments ..
      DOUBLE PRECISION                  SCALE, SSQ
C     ..
C
C  F06BMF returns the value norm given by
C
C     norm = ( scale*sqrt( ssq ), scale*sqrt( ssq ) .lt. flmax
C            (
C            ( flmax,             scale*sqrt( ssq ) .ge. flmax
C
C  via the function name.
C
C
C  Nag Fortran 77 O( 1 ) basic linear algebra routine.
C
C  -- Written on 22-October-1982.
C     Sven Hammarling, Nag Central Office.
C
C
C     .. Local Scalars ..
      DOUBLE PRECISION      FLMAX, FLMIN, NORM, SQT
      LOGICAL               FIRST
C     .. External Functions ..
      DOUBLE PRECISION      X02AMF
      EXTERNAL              X02AMF
C     .. Intrinsic Functions ..
      INTRINSIC             SQRT
C     .. Save statement ..
      SAVE                  FIRST, FLMAX
C     .. Data statements ..
      DATA                  FIRST/ .TRUE. /
C     ..
C     .. Executable Statements ..
      IF( FIRST )THEN
         FIRST = .FALSE.
         FLMIN =  X02AMF( )
         FLMAX =  1/FLMIN
      END IF
C
      SQT = SQRT( SSQ )
      IF( SCALE.LT.FLMAX/SQT )THEN
         NORM = SCALE*SQT
      ELSE
         NORM = FLMAX
      END IF
C
      F06BMF = NORM
      RETURN
C
C     End of F06BMF. ( SNORM )
C
      END
*
      SUBROUTINE F06GGF( N, X, INCX, Y, INCY )
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     .. Entry Points ..
      ENTRY      ZSWAP ( N, X, INCX, Y, INCY )
C     .. Scalar Arguments ..
      INTEGER            INCX, INCY, N
C     .. Array Arguments ..
      COMPLEX*16         X( * ), Y( * )
C     ..
C
C  F06GGF performs the operations
C
C     temp := x,   x := y,   y := temp.
C
C
C  Nag Fortran 77 version of the Blas routine ZSWAP.
C  Nag Fortran 77 O( n ) basic linear algebra routine.
C
C  -- Written on 26-November-1982.
C     Sven Hammarling, Nag Central Office.
C
C
C     .. Local Scalars ..
      COMPLEX*16         TEMP
      INTEGER            I, IX, IY
C     ..
C     .. Executable Statements ..
      IF( N.GT.0 )THEN
         IF( ( INCX.EQ.INCY ).AND.( INCY.GT.0 ) )THEN
            DO 10, IY = 1, 1 + ( N - 1 )*INCY, INCY
               TEMP    = X( IY )
               X( IY ) = Y( IY )
               Y( IY ) = TEMP
   10       CONTINUE
         ELSE
            IF( INCX.GE.0 )THEN
               IX = 1
            ELSE
               IX = 1 - ( N - 1 )*INCX
            END IF
            IF( INCY.GT.0 )THEN
               DO 20, IY = 1, 1 + ( N - 1 )*INCY, INCY
                  TEMP    = X( IX )
                  X( IX ) = Y( IY )
                  Y( IY ) = TEMP
                  IX      = IX      + INCX
   20          CONTINUE
            ELSE
               IY = 1 - ( N - 1 )*INCY
               DO 30, I = 1, N
                  TEMP    = X( IX )
                  X( IX ) = Y( IY )
                  Y( IY ) = TEMP
                  IY      = IY      + INCY
                  IX      = IX      + INCX
   30          CONTINUE
            END IF
         END IF
      END IF
C
      RETURN
C
C     End of F06GGF. ( ZSWAP )
C
      END
* 
      DOUBLE PRECISION FUNCTION F06JJF( N, X, INCX )
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     .. Entry Points ..
      DOUBLE PRECISION          DZNRM2
      ENTRY                     DZNRM2( N, X, INCX )
C     .. Scalar Arguments ..
      INTEGER                           INCX, N
C     .. Array Arguments ..
      COMPLEX*16                        X( * )
C     ..
C
C  F06JJF returns the euclidean norm of a vector via the function
C  name, so that
C
C     F06JJF := sqrt( conjg( x' )*x )
C
C
C  Nag Fortran 77 version of the Blas routine DZNRM2.
C  Nag Fortran 77 O( n ) basic linear algebra routine.
C
C  -- Written on 25-October-1982.
C     Sven Hammarling, Nag Central Office.
C
C
C     .. Parameters ..
      DOUBLE PRECISION      ONE         , ZERO
      PARAMETER           ( ONE = 1.0D+0, ZERO = 0.0D+0 )
C     .. Local Scalars ..
      DOUBLE PRECISION      NORM, SCALE, SSQ
C     .. External Functions ..
      DOUBLE PRECISION      F06BMF
      EXTERNAL              F06BMF
C     .. External Subroutines ..
      EXTERNAL              F06KJF
C     ..
C     .. Executable Statements ..
      IF( N.LT.1 )THEN
         NORM  = ZERO
      ELSE
         SCALE = ZERO
         SSQ   = ONE
         CALL F06KJF( N, X, INCX, SCALE, SSQ )
         NORM  = F06BMF( SCALE, SSQ )
      END IF
C
      F06JJF = NORM
      RETURN
C
C     End of F06JJF. ( DZNRM2 )
C
      END
*
      SUBROUTINE F08ASX(DIRECT,STOREV,N,K,V,LDV,TAU,T,LDT)
C     MARK 16 RELEASE. NAG COPYRIGHT 1992.
C     ENTRY             ZLARFT(DIRECT,STOREV,N,K,V,LDV,TAU,T,LDT)
C
C  Purpose
C  =======
C
C  ZLARFT forms the triangular factor T of a complex block reflector H
C  of order n, which is defined as a product of k elementary reflectors.
C
C  If DIRECT = 'F', H = H(1) H(2) . . . H(k) and T is upper triangular;
C
C  If DIRECT = 'B', H = H(k) . . . H(2) H(1) and T is lower triangular.
C
C  If STOREV = 'C', the vector which defines the elementary reflector
C  H(i) is stored in the i-th column of the array V, and
C
C     H  =  I - V * T * V'
C
C  If STOREV = 'R', the vector which defines the elementary reflector
C  H(i) is stored in the i-th row of the array V, and
C
C     H  =  I - V' * T * V
C
C  Arguments
C  =========
C
C  DIRECT  (input) CHARACTER*1
C          Specifies the order in which the elementary reflectors are
C          multiplied to form the block reflector:
C          = 'F': H = H(1) H(2) . . . H(k) (Forward)
C          = 'B': H = H(k) . . . H(2) H(1) (Backward)
C
C  STOREV  (input) CHARACTER*1
C          Specifies how the vectors which define the elementary
C          reflectors are stored (see also Further Details):
C          = 'C': columnwise
C          = 'R': rowwise
C
C  N       (input) INTEGER
C          The order of the block reflector H. N >= 0.
C
C  K       (input) INTEGER
C          The order of the triangular factor T (= the number of
C          elementary reflectors). K >= 1.
C
C  V       (input/output) COMPLEX*16 array, dimension
C                               (LDV,K) if STOREV = 'C'
C                               (LDV,N) if STOREV = 'R'
C          The matrix V. See further details.
C
C  LDV     (input) INTEGER
C          The leading dimension of the array V.
C          If STOREV = 'C', LDV >= max(1,N); if STOREV = 'R', LDV >= K.
C
C  TAU     (input) COMPLEX*16 array, dimension (K)
C          TAU(i) must contain the scalar factor of the elementary
C          reflector H(i).
C
C  T       (output) COMPLEX*16 array, dimension (LDT,K)
C          The k by k triangular factor T of the block reflector.
C          If DIRECT = 'F', T is upper triangular; if DIRECT = 'B', T is
C          lower triangular. The rest of the array is not used.
C
C  LDT     (input) INTEGER
C          The leading dimension of the array T. LDT >= K.
C
C  Further Details
C  ===============
C
C  The shape of the matrix V and the storage of the vectors which define
C  the H(i) is best illustrated by the following example with n = 5 and
C  k = 3. The elements equal to 1 are not stored; the corresponding
C  array elements are modified but restored on exit. The rest of the
C  array is not used.
C
C  DIRECT = 'F' and STOREV = 'C':         DIRECT = 'F' and STOREV = 'R':
C
C               V = (  1       )                 V = (  1 v1 v1 v1 v1 )
C                   ( v1  1    )                     (     1 v2 v2 v2 )
C                   ( v1 v2  1 )                     (        1 v3 v3 )
C                   ( v1 v2 v3 )
C                   ( v1 v2 v3 )
C
C  DIRECT = 'B' and STOREV = 'C':         DIRECT = 'B' and STOREV = 'R':
C
C               V = ( v1 v2 v3 )                 V = ( v1 v1  1       )
C                   ( v1 v2 v3 )                     ( v2 v2 v2  1    )
C                   (  1 v2 v3 )                     ( v3 v3 v3 v3  1 )
C                   (     1 v3 )
C                   (        1 )
C
C  -- LAPACK auxiliary routine (adapted for NAG Library)
C     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
C     Courant Institute, Argonne National Lab, and Rice University
C
C  =====================================================================
C
C     .. Parameters ..
      COMPLEX*16        ONE, ZERO
      PARAMETER         (ONE=1.0D+0,ZERO=0.0D+0)
C     .. Scalar Arguments ..
      INTEGER           K, LDT, LDV, N
      CHARACTER         DIRECT, STOREV
C     .. Array Arguments ..
      COMPLEX*16        T(LDT,*), TAU(*), V(LDV,*)
C     .. Local Scalars ..
      COMPLEX*16        VII
      INTEGER           I, J
C     .. External Subroutines ..
      EXTERNAL          F07FRY, ZGEMV, ZTRMV
C     .. Executable Statements ..
C
C     Quick return if possible
C
      IF (N.EQ.0) RETURN
C
      IF ((DIRECT.EQ.'F' .OR. DIRECT.EQ.'f')) THEN
         DO 40 I = 1, K
            IF (TAU(I).EQ.ZERO) THEN
C
C              H(i)  =  I
C
               DO 20 J = 1, I
                  T(J,I) = ZERO
   20          CONTINUE
            ELSE
C
C              general case
C
               VII = V(I,I)
               V(I,I) = ONE
               IF ((STOREV.EQ.'C' .OR. STOREV.EQ.'c')) THEN
C
C                 T(1:i-1,i) := - tau(i) * V(i:n,1:i-1)' * V(i:n,i)
C
                  CALL ZGEMV('Conjugate transpose',N-I+1,I-1,-TAU(I),
     *                       V(I,1),LDV,V(I,I),1,ZERO,T(1,I),1)
               ELSE
C
C                 T(1:i-1,i) := - tau(i) * V(1:i-1,i:n) * V(i,i:n)'
C
                  IF (I.LT.N) CALL F07FRY(N-I,V(I,I+1),LDV)
                  CALL ZGEMV('No transpose',I-1,N-I+1,-TAU(I),V(1,I),
     *                       LDV,V(I,I),LDV,ZERO,T(1,I),1)
                  IF (I.LT.N) CALL F07FRY(N-I,V(I,I+1),LDV)
               END IF
               V(I,I) = VII
C
C              T(1:i-1,i) := T(1:i-1,1:i-1) * T(1:i-1,i)
C
               CALL ZTRMV('Upper','No transpose','Non-unit',I-1,T,LDT,
     *                    T(1,I),1)
               T(I,I) = TAU(I)
            END IF
   40    CONTINUE
      ELSE
         DO 80 I = K, 1, -1
            IF (TAU(I).EQ.ZERO) THEN
C
C              H(i)  =  I
C
               DO 60 J = I, K
                  T(J,I) = ZERO
   60          CONTINUE
            ELSE
C
C              general case
C
               IF (I.LT.K) THEN
                  IF ((STOREV.EQ.'C' .OR. STOREV.EQ.'c')) THEN
                     VII = V(N-K+I,I)
                     V(N-K+I,I) = ONE
C
C                    T(i+1:k,i) :=
C                            - tau(i) * V(1:n-k+i,i+1:k)' * V(1:n-k+i,i)
C
                     CALL ZGEMV('Conjugate transpose',N-K+I,K-I,-TAU(I),
     *                          V(1,I+1),LDV,V(1,I),1,ZERO,T(I+1,I),1)
                     V(N-K+I,I) = VII
                  ELSE
                     VII = V(I,N-K+I)
                     V(I,N-K+I) = ONE
C
C                    T(i+1:k,i) :=
C                            - tau(i) * V(i+1:k,1:n-k+i) * V(i,1:n-k+i)'
C
                     CALL F07FRY(N-K+I-1,V(I,1),LDV)
                     CALL ZGEMV('No transpose',K-I,N-K+I,-TAU(I),
     *                          V(I+1,1),LDV,V(I,1),LDV,ZERO,T(I+1,I),1)
                     CALL F07FRY(N-K+I-1,V(I,1),LDV)
                     V(I,N-K+I) = VII
                  END IF
C
C                 T(i+1:k,i) := T(i+1:k,i+1:k) * T(i+1:k,i)
C
                  CALL ZTRMV('Lower','No transpose','Non-unit',K-I,
     *                       T(I+1,I+1),LDT,T(I+1,I),1)
               END IF
               T(I,I) = TAU(I)
            END IF
   80    CONTINUE
      END IF
      RETURN
C
C     End of F08ASX (ZLARFT)
C
      END
*
      SUBROUTINE F08ASY(SIDE,TRANS,DIRECT,STOREV,M,N,K,V,LDV,T,LDT,C,
     *                  LDC,WORK,LDWORK)
C     MARK 16 RELEASE. NAG COPYRIGHT 1992.
C     ENTRY             ZLARFB(SIDE,TRANS,DIRECT,STOREV,M,N,K,V,LDV,T,
C    *                  LDT,C,LDC,WORK,LDWORK)
C
C  Purpose
C  =======
C
C  ZLARFB applies a complex block reflector H or its transpose H' to a
C  complex m by n matrix C, from either the left or the right.
C
C  Arguments
C  =========
C
C  SIDE    (input) CHARACTER*1
C          = 'L': apply H or H' from the Left
C          = 'R': apply H or H' from the Right
C
C  TRANS   (input) CHARACTER*1
C          = 'N': apply H (No transpose)
C          = 'C': apply H' (Conjugate transpose)
C
C  DIRECT  (input) CHARACTER*1
C          Indicates how H is formed from a product of elementary
C          reflectors
C          = 'F': H = H(1) H(2) . . . H(k) (Forward)
C          = 'B': H = H(k) . . . H(2) H(1) (Backward)
C
C  STOREV  (input) CHARACTER*1
C          Indicates how the vectors which define the elementary
C          reflectors are stored:
C          = 'C': Columnwise
C          = 'R': Rowwise
C
C  M       (input) INTEGER
C          The number of rows of the matrix C.
C
C  N       (input) INTEGER
C          The number of columns of the matrix C.
C
C  K       (input) INTEGER
C          The order of the matrix T (= the number of elementary
C          reflectors whose product defines the block reflector).
C
C  V       (input) COMPLEX*16 array, dimension
C                                (LDV,K) if STOREV = 'C'
C                                (LDV,M) if STOREV = 'R' and SIDE = 'L'
C                                (LDV,N) if STOREV = 'R' and SIDE = 'R'
C          The matrix V. See further details.
C
C  LDV     (input) INTEGER
C          The leading dimension of the array V.
C          If STOREV = 'C' and SIDE = 'L', LDV >= max(1,M);
C          if STOREV = 'C' and SIDE = 'R', LDV >= max(1,N);
C          if STOREV = 'R', LDV >= K.
C
C  T       (input) COMPLEX*16 array, dimension (LDT,K)
C          The triangular k-by-k matrix T in the representation of the
C          block reflector.
C
C  LDT     (input) INTEGER
C          The leading dimension of the array T. LDT >= K.
C
C  C       (input/output) COMPLEX*16 array, dimension (LDC,N)
C          On entry, the m-by-n matrix C.
C          On exit, C is overwritten by H*C or H'*C or C*H or C*H'.
C
C  LDC     (input) INTEGER
C          The leading dimension of the array C. LDA >= max(1,M).
C
C  WORK    (workspace) COMPLEX*16 array, dimension (LDWORK,K)
C
C  LDWORK  (input) INTEGER
C          The leading dimension of the array WORK.
C          If SIDE = 'L', LDWORK >= max(1,N);
C          if SIDE = 'R', LDWORK >= max(1,M).
C
C  -- LAPACK auxiliary routine (adapted for NAG Library)
C     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
C     Courant Institute, Argonne National Lab, and Rice University
C
C  =====================================================================
C
C     .. Parameters ..
      COMPLEX*16        ONE
      PARAMETER         (ONE=1.0D+0)
C     .. Scalar Arguments ..
      INTEGER           K, LDC, LDT, LDV, LDWORK, M, N
      CHARACTER         DIRECT, SIDE, STOREV, TRANS
C     .. Array Arguments ..
      COMPLEX*16        C(LDC,*), T(LDT,*), V(LDV,*), WORK(LDWORK,*)
C     .. Local Scalars ..
      INTEGER           I, J
      CHARACTER         TRANST
C     .. External Subroutines ..
      EXTERNAL          F07FRY, ZCOPY, ZGEMM, ZTRMM
C     .. Intrinsic Functions ..
      INTRINSIC         DCONJG
C     .. Executable Statements ..
C
C     Quick return if possible
C
      IF (M.LE.0 .OR. N.LE.0) RETURN
C
      IF ((TRANS.EQ.'N' .OR. TRANS.EQ.'n')) THEN
         TRANST = 'C'
      ELSE
         TRANST = 'N'
      END IF
C
      IF ((STOREV.EQ.'C' .OR. STOREV.EQ.'c')) THEN
C
         IF ((DIRECT.EQ.'F' .OR. DIRECT.EQ.'f')) THEN
C
C           Let  V =  ( V1 )    (first K rows)
C                     ( V2 )
C           where  V1  is unit lower triangular.
C
            IF ((SIDE.EQ.'L' .OR. SIDE.EQ.'l')) THEN
C
C              Form  H * C  or  H' * C  where  C = ( C1 )
C                                                  ( C2 )
C
C              W := C' * V  =  (C1'*V1 + C2'*V2)  (stored in WORK)
C
C              W := C1'
C
               DO 20 J = 1, K
                  CALL ZCOPY(N,C(J,1),LDC,WORK(1,J),1)
                  CALL F07FRY(N,WORK(1,J),1)
   20          CONTINUE
C
C              W := W * V1
C
               CALL ZTRMM('Right','Lower','No transpose','Unit',N,K,ONE,
     *                    V,LDV,WORK,LDWORK)
               IF (M.GT.K) THEN
C
C                 W := W + C2'*V2
C
                  CALL ZGEMM('Conjugate transpose','No transpose',N,K,
     *                       M-K,ONE,C(K+1,1),LDC,V(K+1,1),LDV,ONE,WORK,
     *                       LDWORK)
               END IF
C
C              W := W * T'  or  W * T
C
               CALL ZTRMM('Right','Upper',TRANST,'Non-unit',N,K,ONE,T,
     *                    LDT,WORK,LDWORK)
C
C              C := C - V * W'
C
               IF (M.GT.K) THEN
C
C                 C2 := C2 - V2 * W'
C
                  CALL ZGEMM('No transpose','Conjugate transpose',M-K,N,
     *                       K,-ONE,V(K+1,1),LDV,WORK,LDWORK,ONE,
     *                       C(K+1,1),LDC)
               END IF
C
C              W := W * V1'
C
               CALL ZTRMM('Right','Lower','Conjugate transpose','Unit',
     *                    N,K,ONE,V,LDV,WORK,LDWORK)
C
C              C1 := C1 - W'
C
               DO 60 J = 1, K
                  DO 40 I = 1, N
                     C(J,I) = C(J,I) - DCONJG(WORK(I,J))
   40             CONTINUE
   60          CONTINUE
C
            ELSE IF ((SIDE.EQ.'R' .OR. SIDE.EQ.'r')) THEN
C
C              Form  C * H  or  C * H'  where  C = ( C1  C2 )
C
C              W := C * V  =  (C1*V1 + C2*V2)  (stored in WORK)
C
C              W := C1
C
               DO 80 J = 1, K
                  CALL ZCOPY(M,C(1,J),1,WORK(1,J),1)
   80          CONTINUE
C
C              W := W * V1
C
               CALL ZTRMM('Right','Lower','No transpose','Unit',M,K,ONE,
     *                    V,LDV,WORK,LDWORK)
               IF (N.GT.K) THEN
C
C                 W := W + C2 * V2
C
                  CALL ZGEMM('No transpose','No transpose',M,K,N-K,ONE,
     *                       C(1,K+1),LDC,V(K+1,1),LDV,ONE,WORK,LDWORK)
               END IF
C
C              W := W * T  or  W * T'
C
               CALL ZTRMM('Right','Upper',TRANS,'Non-unit',M,K,ONE,T,
     *                    LDT,WORK,LDWORK)
C
C              C := C - W * V'
C
               IF (N.GT.K) THEN
C
C                 C2 := C2 - W * V2'
C
                  CALL ZGEMM('No transpose','Conjugate transpose',M,N-K,
     *                       K,-ONE,WORK,LDWORK,V(K+1,1),LDV,ONE,
     *                       C(1,K+1),LDC)
               END IF
C
C              W := W * V1'
C
               CALL ZTRMM('Right','Lower','Conjugate transpose','Unit',
     *                    M,K,ONE,V,LDV,WORK,LDWORK)
C
C              C1 := C1 - W
C
               DO 120 J = 1, K
                  DO 100 I = 1, M
                     C(I,J) = C(I,J) - WORK(I,J)
  100             CONTINUE
  120          CONTINUE
            END IF
C
         ELSE
C
C           Let  V =  ( V1 )
C                     ( V2 )    (last K rows)
C           where  V2  is unit upper triangular.
C
            IF ((SIDE.EQ.'L' .OR. SIDE.EQ.'l')) THEN
C
C              Form  H * C  or  H' * C  where  C = ( C1 )
C                                                  ( C2 )
C
C              W := C' * V  =  (C1'*V1 + C2'*V2)  (stored in WORK)
C
C              W := C2'
C
               DO 140 J = 1, K
                  CALL ZCOPY(N,C(M-K+J,1),LDC,WORK(1,J),1)
                  CALL F07FRY(N,WORK(1,J),1)
  140          CONTINUE
C
C              W := W * V2
C
               CALL ZTRMM('Right','Upper','No transpose','Unit',N,K,ONE,
     *                    V(M-K+1,1),LDV,WORK,LDWORK)
               IF (M.GT.K) THEN
C
C                 W := W + C1'*V1
C
                  CALL ZGEMM('Conjugate transpose','No transpose',N,K,
     *                       M-K,ONE,C,LDC,V,LDV,ONE,WORK,LDWORK)
               END IF
C
C              W := W * T'  or  W * T
C
               CALL ZTRMM('Right','Lower',TRANST,'Non-unit',N,K,ONE,T,
     *                    LDT,WORK,LDWORK)
C
C              C := C - V * W'
C
               IF (M.GT.K) THEN
C
C                 C1 := C1 - V1 * W'
C
                  CALL ZGEMM('No transpose','Conjugate transpose',M-K,N,
     *                       K,-ONE,V,LDV,WORK,LDWORK,ONE,C,LDC)
               END IF
C
C              W := W * V2'
C
               CALL ZTRMM('Right','Upper','Conjugate transpose','Unit',
     *                    N,K,ONE,V(M-K+1,1),LDV,WORK,LDWORK)
C
C              C2 := C2 - W'
C
               DO 180 J = 1, K
                  DO 160 I = 1, N
                     C(M-K+J,I) = C(M-K+J,I) - DCONJG(WORK(I,J))
  160             CONTINUE
  180          CONTINUE
C
            ELSE IF ((SIDE.EQ.'R' .OR. SIDE.EQ.'r')) THEN
C
C              Form  C * H  or  C * H'  where  C = ( C1  C2 )
C
C              W := C * V  =  (C1*V1 + C2*V2)  (stored in WORK)
C
C              W := C2
C
               DO 200 J = 1, K
                  CALL ZCOPY(M,C(1,N-K+J),1,WORK(1,J),1)
  200          CONTINUE
C
C              W := W * V2
C
               CALL ZTRMM('Right','Upper','No transpose','Unit',M,K,ONE,
     *                    V(N-K+1,1),LDV,WORK,LDWORK)
               IF (N.GT.K) THEN
C
C                 W := W + C1 * V1
C
                  CALL ZGEMM('No transpose','No transpose',M,K,N-K,ONE,
     *                       C,LDC,V,LDV,ONE,WORK,LDWORK)
               END IF
C
C              W := W * T  or  W * T'
C
               CALL ZTRMM('Right','Lower',TRANS,'Non-unit',M,K,ONE,T,
     *                    LDT,WORK,LDWORK)
C
C              C := C - W * V'
C
               IF (N.GT.K) THEN
C
C                 C1 := C1 - W * V1'
C
                  CALL ZGEMM('No transpose','Conjugate transpose',M,N-K,
     *                       K,-ONE,WORK,LDWORK,V,LDV,ONE,C,LDC)
               END IF
C
C              W := W * V2'
C
               CALL ZTRMM('Right','Upper','Conjugate transpose','Unit',
     *                    M,K,ONE,V(N-K+1,1),LDV,WORK,LDWORK)
C
C              C2 := C2 - W
C
               DO 240 J = 1, K
                  DO 220 I = 1, M
                     C(I,N-K+J) = C(I,N-K+J) - WORK(I,J)
  220             CONTINUE
  240          CONTINUE
            END IF
         END IF
C
      ELSE IF ((STOREV.EQ.'R' .OR. STOREV.EQ.'r')) THEN
C
         IF ((DIRECT.EQ.'F' .OR. DIRECT.EQ.'f')) THEN
C
C           Let  V =  ( V1  V2 )    (V1: first K columns)
C           where  V1  is unit upper triangular.
C
            IF ((SIDE.EQ.'L' .OR. SIDE.EQ.'l')) THEN
C
C              Form  H * C  or  H' * C  where  C = ( C1 )
C                                                  ( C2 )
C
C              W := C' * V'  =  (C1'*V1' + C2'*V2') (stored in WORK)
C
C              W := C1'
C
               DO 260 J = 1, K
                  CALL ZCOPY(N,C(J,1),LDC,WORK(1,J),1)
                  CALL F07FRY(N,WORK(1,J),1)
  260          CONTINUE
C
C              W := W * V1'
C
               CALL ZTRMM('Right','Upper','Conjugate transpose','Unit',
     *                    N,K,ONE,V,LDV,WORK,LDWORK)
               IF (M.GT.K) THEN
C
C                 W := W + C2'*V2'
C
                  CALL ZGEMM('Conjugate transpose',
     *                       'Conjugate transpose',N,K,M-K,ONE,C(K+1,1),
     *                       LDC,V(1,K+1),LDV,ONE,WORK,LDWORK)
               END IF
C
C              W := W * T'  or  W * T
C
               CALL ZTRMM('Right','Upper',TRANST,'Non-unit',N,K,ONE,T,
     *                    LDT,WORK,LDWORK)
C
C              C := C - V' * W'
C
               IF (M.GT.K) THEN
C
C                 C2 := C2 - V2' * W'
C
                  CALL ZGEMM('Conjugate transpose',
     *                       'Conjugate transpose',M-K,N,K,-ONE,V(1,K+1)
     *                       ,LDV,WORK,LDWORK,ONE,C(K+1,1),LDC)
               END IF
C
C              W := W * V1
C
               CALL ZTRMM('Right','Upper','No transpose','Unit',N,K,ONE,
     *                    V,LDV,WORK,LDWORK)
C
C              C1 := C1 - W'
C
               DO 300 J = 1, K
                  DO 280 I = 1, N
                     C(J,I) = C(J,I) - DCONJG(WORK(I,J))
  280             CONTINUE
  300          CONTINUE
C
            ELSE IF ((SIDE.EQ.'R' .OR. SIDE.EQ.'r')) THEN
C
C              Form  C * H  or  C * H'  where  C = ( C1  C2 )
C
C              W := C * V'  =  (C1*V1' + C2*V2')  (stored in WORK)
C
C              W := C1
C
               DO 320 J = 1, K
                  CALL ZCOPY(M,C(1,J),1,WORK(1,J),1)
  320          CONTINUE
C
C              W := W * V1'
C
               CALL ZTRMM('Right','Upper','Conjugate transpose','Unit',
     *                    M,K,ONE,V,LDV,WORK,LDWORK)
               IF (N.GT.K) THEN
C
C                 W := W + C2 * V2'
C
                  CALL ZGEMM('No transpose','Conjugate transpose',M,K,
     *                       N-K,ONE,C(1,K+1),LDC,V(1,K+1),LDV,ONE,WORK,
     *                       LDWORK)
               END IF
C
C              W := W * T  or  W * T'
C
               CALL ZTRMM('Right','Upper',TRANS,'Non-unit',M,K,ONE,T,
     *                    LDT,WORK,LDWORK)
C
C              C := C - W * V
C
               IF (N.GT.K) THEN
C
C                 C2 := C2 - W * V2
C
                  CALL ZGEMM('No transpose','No transpose',M,N-K,K,-ONE,
     *                       WORK,LDWORK,V(1,K+1),LDV,ONE,C(1,K+1),LDC)
               END IF
C
C              W := W * V1
C
               CALL ZTRMM('Right','Upper','No transpose','Unit',M,K,ONE,
     *                    V,LDV,WORK,LDWORK)
C
C              C1 := C1 - W
C
               DO 360 J = 1, K
                  DO 340 I = 1, M
                     C(I,J) = C(I,J) - WORK(I,J)
  340             CONTINUE
  360          CONTINUE
C
            END IF
C
         ELSE
C
C           Let  V =  ( V1  V2 )    (V2: last K columns)
C           where  V2  is unit lower triangular.
C
            IF ((SIDE.EQ.'L' .OR. SIDE.EQ.'l')) THEN
C
C              Form  H * C  or  H' * C  where  C = ( C1 )
C                                                  ( C2 )
C
C              W := C' * V'  =  (C1'*V1' + C2'*V2') (stored in WORK)
C
C              W := C2'
C
               DO 380 J = 1, K
                  CALL ZCOPY(N,C(M-K+J,1),LDC,WORK(1,J),1)
                  CALL F07FRY(N,WORK(1,J),1)
  380          CONTINUE
C
C              W := W * V2'
C
               CALL ZTRMM('Right','Lower','Conjugate transpose','Unit',
     *                    N,K,ONE,V(1,M-K+1),LDV,WORK,LDWORK)
               IF (M.GT.K) THEN
C
C                 W := W + C1'*V1'
C
                  CALL ZGEMM('Conjugate transpose',
     *                       'Conjugate transpose',N,K,M-K,ONE,C,LDC,V,
     *                       LDV,ONE,WORK,LDWORK)
               END IF
C
C              W := W * T'  or  W * T
C
               CALL ZTRMM('Right','Lower',TRANST,'Non-unit',N,K,ONE,T,
     *                    LDT,WORK,LDWORK)
C
C              C := C - V' * W'
C
               IF (M.GT.K) THEN
C
C                 C1 := C1 - V1' * W'
C
                  CALL ZGEMM('Conjugate transpose',
     *                       'Conjugate transpose',M-K,N,K,-ONE,V,LDV,
     *                       WORK,LDWORK,ONE,C,LDC)
               END IF
C
C              W := W * V2
C
               CALL ZTRMM('Right','Lower','No transpose','Unit',N,K,ONE,
     *                    V(1,M-K+1),LDV,WORK,LDWORK)
C
C              C2 := C2 - W'
C
               DO 420 J = 1, K
                  DO 400 I = 1, N
                     C(M-K+J,I) = C(M-K+J,I) - DCONJG(WORK(I,J))
  400             CONTINUE
  420          CONTINUE
C
            ELSE IF ((SIDE.EQ.'R' .OR. SIDE.EQ.'r')) THEN
C
C              Form  C * H  or  C * H'  where  C = ( C1  C2 )
C
C              W := C * V'  =  (C1*V1' + C2*V2')  (stored in WORK)
C
C              W := C2
C
               DO 440 J = 1, K
                  CALL ZCOPY(M,C(1,N-K+J),1,WORK(1,J),1)
  440          CONTINUE
C
C              W := W * V2'
C
               CALL ZTRMM('Right','Lower','Conjugate transpose','Unit',
     *                    M,K,ONE,V(1,N-K+1),LDV,WORK,LDWORK)
               IF (N.GT.K) THEN
C
C                 W := W + C1 * V1'
C
                  CALL ZGEMM('No transpose','Conjugate transpose',M,K,
     *                       N-K,ONE,C,LDC,V,LDV,ONE,WORK,LDWORK)
               END IF
C
C              W := W * T  or  W * T'
C
               CALL ZTRMM('Right','Lower',TRANS,'Non-unit',M,K,ONE,T,
     *                    LDT,WORK,LDWORK)
C
C              C := C - W * V
C
               IF (N.GT.K) THEN
C
C                 C1 := C1 - W * V1
C
                  CALL ZGEMM('No transpose','No transpose',M,N-K,K,-ONE,
     *                       WORK,LDWORK,V,LDV,ONE,C,LDC)
               END IF
C
C              W := W * V2
C
               CALL ZTRMM('Right','Lower','No transpose','Unit',M,K,ONE,
     *                    V(1,N-K+1),LDV,WORK,LDWORK)
C
C              C1 := C1 - W
C
               DO 480 J = 1, K
                  DO 460 I = 1, M
                     C(I,N-K+J) = C(I,N-K+J) - WORK(I,J)
  460             CONTINUE
  480          CONTINUE
C
            END IF
C
         END IF
      END IF
C
      RETURN
C
C     End of F08ASY (ZLARFB)
C
      END
*
      SUBROUTINE F08ATF(M,N,K,A,LDA,TAU,WORK,LWORK,INFO)
C     MARK 16 RELEASE. NAG COPYRIGHT 1992.
C     .. Entry Points ..
      ENTRY             ZUNGQR(M,N,K,A,LDA,TAU,WORK,LWORK,INFO)
C
C  Purpose
C  =======
C
C  ZUNGQR generates an m-by-n complex matrix Q with orthonormal columns,
C  which is defined as the first n columns of a product of k elementary
C  reflectors of order m
C
C        Q  =  H(1) H(2) . . . H(k)
C
C  as returned by ZGEQRF.
C
C  Arguments
C  =========
C
C  M       (input) INTEGER
C          The number of rows of the matrix Q. M >= 0.
C
C  N       (input) INTEGER
C          The number of columns of the matrix Q. M >= N >= 0.
C
C  K       (input) INTEGER
C          The number of elementary reflectors whose product defines the
C          matrix Q. N >= K >= 0.
C
C  A       (input/output) COMPLEX*16 array, dimension (LDA,N)
C          On entry, the i-th column must contain the vector which
C          defines the elementary reflector H(i), for i = 1,2,...,k, as
C          returned by ZGEQRF in the first k columns of its array
C          argument A.
C          On exit, the m by n matrix Q.
C
C  LDA     (input) INTEGER
C          The first dimension of the array A. LDA >= max(1,M).
C
C  TAU     (input) COMPLEX*16 array, dimension (K)
C          TAU(i) must contain the scalar factor of the elementary
C          reflector H(i), as returned by ZGEQRF.
C
C  WORK    (workspace) COMPLEX*16 array, dimension (LWORK)
C          On exit, if INFO = 0, WORK(1) returns the minimum value of
C          LWORK required to use the optimal blocksize.
C
C  LWORK   (input) INTEGER
C          The dimension of the array WORK. LWORK >= max(1,N).
C          For optimum performance LWORK should be at least N*NB, where
C          NB is the optimal blocksize.
C
C  INFO    (output) INTEGER
C          = 0: successful exit;
C          < 0: if INFO = -i, the i-th argument has an illegal value
C
C  -- LAPACK routine (adapted for NAG Library)
C     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
C     Courant Institute, Argonne National Lab, and Rice University
C
C  =====================================================================
C
C     .. Parameters ..
      COMPLEX*16        ZERO
      PARAMETER         (ZERO=0.0D+0)
C     .. Scalar Arguments ..
      INTEGER           INFO, K, LDA, LWORK, M, N
C     .. Array Arguments ..
      COMPLEX*16        A(LDA,*), TAU(*), WORK(LWORK)
C     .. Local Scalars ..
      INTEGER           I, IB, IINFO, IWS, J, KI, KK, L, LDWORK, NB,
     *                  NBMIN, NX
C     .. External Subroutines ..
      EXTERNAL          F06AAZ, F07ZAZ, F08ASX, F08ASY, F08ATZ
C     .. Intrinsic Functions ..
      INTRINSIC         MAX, MIN
C     .. Executable Statements ..
C
C     Test the input arguments
C
      INFO = 0
      IF (M.LT.0) THEN
         INFO = -1
      ELSE IF (N.LT.0 .OR. N.GT.M) THEN
         INFO = -2
      ELSE IF (K.LT.0 .OR. K.GT.N) THEN
         INFO = -3
      ELSE IF (LDA.LT.MAX(1,M)) THEN
         INFO = -5
      ELSE IF (LWORK.LT.MAX(1,N)) THEN
         INFO = -8
      END IF
      IF (INFO.NE.0) THEN
         CALL F06AAZ('F08ATF/ZUNGQR',-INFO)
         RETURN
      END IF
C
C     Quick return if possible
C
      IF (N.LE.0) THEN
         WORK(1) = 1
         RETURN
      END IF
C
C     Determine the block size.
C
      CALL F07ZAZ(1,'F08ATF',NB,0)
      NBMIN = 2
      NX = 0
      IWS = N
      IF (NB.GT.1 .AND. NB.LT.K) THEN
C
C        Determine when to cross over from blocked to unblocked code.
C
         CALL F07ZAZ(3,'F08ATF',NX,0)
         NX = MAX(0,NX)
         IF (NX.LT.K) THEN
C
C           Determine if workspace is large enough for blocked code.
C
            LDWORK = N
            IWS = LDWORK*NB
            IF (LWORK.LT.IWS) THEN
C
C              Not enough workspace to use optimal NB:  reduce NB and
C              determine the minimum value of NB.
C
               NB = LWORK/LDWORK
               CALL F07ZAZ(2,'F08ATF',NBMIN,0)
               NBMIN = MAX(2,NBMIN)
            END IF
         END IF
      END IF
C
      IF (NB.GE.NBMIN .AND. NB.LT.K .AND. NX.LT.K) THEN
C
C        Use blocked code after the last block.
C        The first kk columns are handled by the block method.
C
         KI = ((K-NX-1)/NB)*NB
         KK = MIN(K,KI+NB)
C
C        Set A(1:kk,kk+1:n) to zero.
C
         DO 40 J = KK + 1, N
            DO 20 I = 1, KK
               A(I,J) = ZERO
   20       CONTINUE
   40    CONTINUE
      ELSE
         KK = 0
      END IF
C
C     Use unblocked code for the last or only block.
C
      IF (KK.LT.N) CALL F08ATZ(M-KK,N-KK,K-KK,A(KK+1,KK+1),LDA,TAU(KK+1)
     *                         ,WORK,IINFO)
C
      IF (KK.GT.0) THEN
C
C        Use blocked code
C
         DO 100 I = KI + 1, 1, -NB
            IB = MIN(NB,K-I+1)
            IF (I+IB.LE.N) THEN
C
C              Form the triangular factor of the block reflector
C              H = H(i) H(i+1) . . . H(i+ib-1)
C
               CALL F08ASX('Forward','Columnwise',M-I+1,IB,A(I,I),LDA,
     *                     TAU(I),WORK,LDWORK)
C
C              Apply H to A(i:m,i+ib:n) from the left
C
               CALL F08ASY('Left','No transpose','Forward','Columnwise',
     *                     M-I+1,N-I-IB+1,IB,A(I,I),LDA,WORK,LDWORK,
     *                     A(I,I+IB),LDA,WORK(IB+1),LDWORK)
            END IF
C
C           Apply H to rows i:m of current block
C
            CALL F08ATZ(M-I+1,IB,IB,A(I,I),LDA,TAU(I),WORK,IINFO)
C
C           Set rows 1:i-1 of current block to zero
C
            DO 80 J = I, I + IB - 1
               DO 60 L = 1, I - 1
                  A(L,J) = ZERO
   60          CONTINUE
   80       CONTINUE
  100    CONTINUE
      END IF
C
      WORK(1) = IWS
      RETURN
C
C     End of F08ATF (ZUNGQR)
C
      END
*
      SUBROUTINE F08CTZ(M,N,K,A,LDA,TAU,WORK,INFO)
C     MARK 16 RELEASE. NAG COPYRIGHT 1992.
C     ENTRY             ZUNG2L(M,N,K,A,LDA,TAU,WORK,INFO)
C
C  Purpose
C  =======
C
C  ZUNG2L generates an m by n complex matrix Q with orthonormal columns,
C  which is defined as the last n columns of a product of k elementary
C  reflectors of order m
C
C        Q  =  H(k) . . . H(2) H(1)
C
C  as returned by ZGEQLF.
C
C  Arguments
C  =========
C
C  M       (input) INTEGER
C          The number of rows of the matrix Q. M >= 0.
C
C  N       (input) INTEGER
C          The number of columns of the matrix Q. M >= N >= 0.
C
C  K       (input) INTEGER
C          The number of elementary reflectors whose product defines the
C          matrix Q. N >= K >= 0.
C
C  A       (input/output) COMPLEX*16 array, dimension (LDA,N)
C          On entry, the (n-k+i)-th column must contain the vector which
C          defines the elementary reflector H(i), for i = 1,2,...,k, as
C          returned by ZGEQLF in the last k columns of its array
C          argument A.
C          On exit, the m-by-n matrix Q.
C
C  LDA     (input) INTEGER
C          The first dimension of the array A. LDA >= max(1,M).
C
C  TAU     (input) COMPLEX*16 array, dimension (K)
C          TAU(i) must contain the scalar factor of the elementary
C          reflector H(i), as returned by ZGEQLF.
C
C  WORK    (workspace) COMPLEX*16 array, dimension (N)
C
C  INFO    (output) INTEGER
C          = 0: successful exit
C          < 0: if INFO = -i, the i-th argument has an illegal value
C
C  -- LAPACK routine (adapted for NAG Library)
C     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
C     Courant Institute, Argonne National Lab, and Rice University
C
C  =====================================================================
C
C     .. Parameters ..
      COMPLEX*16        ONE, ZERO
      PARAMETER         (ONE=1.0D+0,ZERO=0.0D+0)
C     .. Scalar Arguments ..
      INTEGER           INFO, K, LDA, M, N
C     .. Array Arguments ..
      COMPLEX*16        A(LDA,*), TAU(*), WORK(*)
C     .. Local Scalars ..
      INTEGER           I, II, J, L
C     .. External Subroutines ..
      EXTERNAL          F06AAZ, F08ASW, ZSCAL
C     .. Intrinsic Functions ..
      INTRINSIC         MAX
C     .. Executable Statements ..
C
C     Test the input arguments
C
      INFO = 0
      IF (M.LT.0) THEN
         INFO = -1
      ELSE IF (N.LT.0 .OR. N.GT.M) THEN
         INFO = -2
      ELSE IF (K.LT.0 .OR. K.GT.N) THEN
         INFO = -3
      ELSE IF (LDA.LT.MAX(1,M)) THEN
         INFO = -5
      END IF
      IF (INFO.NE.0) THEN
         CALL F06AAZ('F08CTZ/ZUNG2L',-INFO)
         RETURN
      END IF
C
C     Quick return if possible
C
      IF (N.LE.0) RETURN
C
C     Initialise columns 1:n-k to columns of the unit matrix
C
      DO 40 J = 1, N - K
         DO 20 L = 1, M
            A(L,J) = ZERO
   20    CONTINUE
         A(M-N+J,J) = ONE
   40 CONTINUE
C
      DO 80 I = 1, K
         II = N - K + I
C
C        Apply H(i) to A(1:m-k+i,1:n-k+i) from the left
C
         A(M-N+II,II) = ONE
         CALL F08ASW('Left',M-N+II,II-1,A(1,II),1,TAU(I),A,LDA,WORK)
         CALL ZSCAL(M-N+II-1,-TAU(I),A(1,II),1)
         A(M-N+II,II) = ONE - TAU(I)
C
C        Set A(m-k+i+1:m,n-k+i) to zero
C
         DO 60 L = M - N + II + 1, M
            A(L,II) = ZERO
   60    CONTINUE
   80 CONTINUE
      RETURN
C
C     End of F08CTZ (ZUNG2L)
C
      END
*
      SUBROUTINE F06GFF( N, X, INCX, Y, INCY )
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     .. Entry Points ..
      ENTRY      ZCOPY ( N, X, INCX, Y, INCY )
C     .. Scalar Arguments ..
      INTEGER            INCX, INCY, N
C     .. Array Arguments ..
      COMPLEX*16         X( * ), Y( * )
C     ..
C
C  F06GFF performs the operation
C
C     y := x
C
C
C  Nag Fortran 77 version of the Blas routine ZCOPY.
C  Nag Fortran 77 O( n ) basic linear algebra routine.
C
C  -- Written on 26-November-1982.
C     Sven Hammarling, Nag Central Office.
C
C
C     .. Local Scalars ..
      INTEGER            I, IX, IY
C     ..
C     .. Executable Statements ..
      IF( N.GT.0 )THEN
         IF( ( INCX.EQ.INCY ).AND.( INCY.GT.0 ) )THEN
            DO 10, IY = 1, 1 + ( N - 1 )*INCY, INCY
               Y( IY ) = X( IY )
   10       CONTINUE
         ELSE
            IF( INCX.GE.0 )THEN
               IX = 1
            ELSE
               IX = 1 - ( N - 1 )*INCX
            END IF
            IF( INCY.GT.0 )THEN
               DO 20, IY = 1, 1 + ( N - 1 )*INCY, INCY
                  Y( IY ) = X( IX )
                  IX      = IX      + INCX
   20          CONTINUE
            ELSE
               IY = 1 - ( N - 1 )*INCY
               DO 30, I = 1, N
                  Y( IY ) = X( IX )
                  IY      = IY + INCY
                  IX      = IX + INCX
   30          CONTINUE
            END IF
         END IF
      END IF
C
      RETURN
C
C     End of F06GFF. ( ZCOPY )
C
      END
*
      DOUBLE PRECISION FUNCTION F06JKF( N, X, INCX )
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     .. Entry Points ..
      DOUBLE PRECISION          DZASUM
      ENTRY                     DZASUM( N, X, INCX )
C     .. Scalar Arguments ..
      INTEGER                           INCX, N
C     .. Array Arguments ..
      COMPLEX*16                        X( * )
C     ..
C
C  F06JKF returns the value, asum, given by
C
C     asum = abs( real( x( 1 ) ) ) + abs( aimag( x( 1 ) ) ) + ... +
C            abs( real( x( n ) ) ) + abs( aimag( x( n ) ) ),
C
C  for the vector x. asum is returned via the function name F06JKF.
C
C
C  Nag Fortran 77 version of the Blas routine DZASUM.
C  Nag Fortran 77 O( n ) basic linear algebra routine.
C
C  -- Written on 12-February-1986.
C     Sven Hammarling, Nag Central Office.
C
C
C     .. Parameters ..
      DOUBLE PRECISION      ZERO
      PARAMETER           ( ZERO = 0.0D+0 )
C     .. Local Scalars ..
      DOUBLE PRECISION      SUM
      INTEGER               IX
C     .. Intrinsic Functions ..
      INTRINSIC             ABS, DIMAG, DBLE
C     ..
C     .. Executable Statements ..
      SUM = ZERO
      IF( N.GT.0 )THEN
         DO 10, IX = 1, 1 + ( N - 1 )*INCX, INCX
            SUM = SUM + ABS( DBLE ( X( IX ) ) )
     $                + ABS( DIMAG( X( IX ) ) )
   10    CONTINUE
      END IF
C
      F06JKF = SUM
      RETURN
C
C     End of F06JKF. ( DZASUM )
C
      END
*
      INTEGER FUNCTION F06JLF( N, X, INCX )
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     .. Entry Points ..
      INTEGER          IDAMAX
      ENTRY            IDAMAX( N, X, INCX )
C     .. Scalar Arguments ..
      INTEGER                  INCX, N
C     .. Array Arguments ..
      DOUBLE PRECISION         X( * )
C     ..
C
C  F06JLF returns the smallest value of i such that
C
C     abs( x( i ) ) = max( abs( x( j ) ) )
C                      j
C
C
C  Nag Fortran 77 version of the Blas routine IDAMAX.
C  Nag Fortran 77 O( n ) basic linear algebra routine.
C
C  -- Written on 31-May-1983.
C     Sven Hammarling, Nag Central Office.
C
C
C     .. Local Scalars ..
      DOUBLE PRECISION         XMAX
      INTEGER                  I, IMAX, IX
C     .. Intrinsic Functions ..
      INTRINSIC                ABS
C     ..
C     .. Executable Statements ..
      IF( N.GT.0 )THEN
         IMAX = 1
         IF( N.GT.1 )THEN
            XMAX = ABS( X( 1 ) )
            IX   = 1
            DO 10, I = 2, N
               IX = IX + INCX
               IF( XMAX.LT.ABS( X( IX ) ) )THEN
                  XMAX = ABS( X( IX ) )
                  IMAX = I
               END IF
   10       CONTINUE
         END IF
      ELSE
         IMAX = 0
      END IF
C
      F06JLF = IMAX
      RETURN
C
C     End of F06JLF. ( IDAMAX )
C
      END
*
      INTEGER FUNCTION F06JMF( N, X, INCX )
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     .. Entry Points ..
      INTEGER          IZAMAX
      ENTRY            IZAMAX( N, X, INCX )
C     .. Scalar Arguments ..
      INTEGER                  INCX, N
C     .. Array Arguments ..
      COMPLEX*16               X( * )
C     ..
C
C  F06JMF returns the smallest value of i such that
C
C     alpha( i ) = max( abs( real( x( j ) ) ) + abs( imag( x( j ) ) ) )
C                   j
C
C  Nag Fortran 77 version of the Blas routine IZAMAX.
C  Nag Fortran 77 O( n ) basic linear algebra routine.
C
C  -- Written on 31-May-1983.
C     Sven Hammarling, Nag Central Office.
C
C
C     .. Local Scalars ..
      DOUBLE PRECISION         TEMP, XMAX
      INTEGER                  I, IMAX, IX
C     .. Intrinsic Functions ..
      INTRINSIC                ABS, DIMAG, DBLE
C     ..
C     .. Executable Statements ..
      IF( N.GT.0 )THEN
         IMAX = 1
         IF( N.GT.1 )THEN
            XMAX = ABS( DBLE( X( 1 ) ) ) + ABS( DIMAG( X( 1 ) ) )
            IX   = 1
            DO 10, I = 2, N
               IX   = IX                     + INCX
               TEMP = ABS( DBLE( X( IX ) ) ) + ABS( DIMAG( X( IX ) ) )
               IF( XMAX.LT.TEMP )THEN
                  XMAX = TEMP
                  IMAX = I
               END IF
   10       CONTINUE
         END IF
      ELSE
         IMAX = 0
      END IF
C
      F06JMF = IMAX
      RETURN
C
C     End of F06JMF. ( IZAMAX )
C
      END
*
      SUBROUTINE F06SFF( UPLO, TRANS, DIAG, N, A, LDA, X, INCX )
C     MARK 13 RE-ISSUE. NAG COPYRIGHT 1988.
C     .. Entry Points ..
      ENTRY      ZTRMV ( UPLO, TRANS, DIAG, N, A, LDA, X, INCX )
C     .. Scalar Arguments ..
      INTEGER            INCX, LDA, N
      CHARACTER*1        DIAG, TRANS, UPLO
C     .. Array Arguments ..
      COMPLEX*16         A( LDA, * ), X( * )
C     ..
C
C  Purpose
C  =======
C
C  ZTRMV  performs one of the matrix-vector operations
C
C     x := A*x,   or   x := A'*x,   or   x := conjg( A' )*x,
C
C  where x is n element vector and A is an n by n unit, or non-unit,
C  upper or lower triangular matrix.
C
C  Parameters
C  ==========
C
C  UPLO   - CHARACTER*1.
C           On entry, UPLO specifies whether the matrix is an upper or
C           lower triangular matrix as follows:
C
C              UPLO = 'U' or 'u'   A is an upper triangular matrix.
C
C              UPLO = 'L' or 'l'   A is a lower triangular matrix.
C
C           Unchanged on exit.
C
C  TRANS  - CHARACTER*1.
C           On entry, TRANS specifies the operation to be performed as
C           follows:
C
C              TRANS = 'N' or 'n'   x := A*x.
C
C              TRANS = 'T' or 't'   x := A'*x.
C
C              TRANS = 'C' or 'c'   x := conjg( A' )*x.
C
C           Unchanged on exit.
C
C  DIAG   - CHARACTER*1.
C           On entry, DIAG specifies whether or not A is unit
C           triangular as follows:
C
C              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
C
C              DIAG = 'N' or 'n'   A is not assumed to be unit
C                                  triangular.
C
C           Unchanged on exit.
C
C  N      - INTEGER.
C           On entry, N specifies the order of the matrix A.
C           N must be at least zero.
C           Unchanged on exit.
C
C  A      - COMPLEX*16       array of DIMENSION ( LDA, n ).
C           Before entry with  UPLO = 'U' or 'u', the leading n by n
C           upper triangular part of the array A must contain the upper
C           triangular matrix and the strictly lower triangular part of
C           A is not referenced.
C           Before entry with UPLO = 'L' or 'l', the leading n by n
C           lower triangular part of the array A must contain the lower
C           triangular matrix and the strictly upper triangular part of
C           A is not referenced.
C           Note that when  DIAG = 'U' or 'u', the diagonal elements of
C           A are not referenced either, but are assumed to be unity.
C           Unchanged on exit.
C
C  LDA    - INTEGER.
C           On entry, LDA specifies the first dimension of A as declared
C           in the calling (sub) program. LDA must be at least
C           max( 1, n ).
C           Unchanged on exit.
C
C  X      - COMPLEX*16       array of dimension at least
C           ( 1 + ( n - 1 )*abs( INCX ) ).
C           Before entry, the incremented array X must contain the n
C           element vector x. On exit, X is overwritten with the
C           tranformed vector x.
C
C  INCX   - INTEGER.
C           On entry, INCX specifies the increment for the elements of
C           X. INCX must not be zero.
C           Unchanged on exit.
C
C
C  Level 2 Blas routine.
C
C  -- Written on 22-October-1986.
C     Jack Dongarra, Argonne National Lab.
C     Jeremy Du Croz, Nag Central Office.
C     Sven Hammarling, Nag Central Office.
C     Richard Hanson, Sandia National Labs.
C
C
C     .. Parameters ..
      COMPLEX*16         ZERO
      PARAMETER        ( ZERO = ( 0.0D+0, 0.0D+0 ) )
C     .. Local Scalars ..
      COMPLEX*16         TEMP
      INTEGER            I, INFO, IX, J, JX, KX
      LOGICAL            NOCONJ, NOUNIT
C     .. External Subroutines ..
      EXTERNAL           F06AAZ
C     .. Intrinsic Functions ..
      INTRINSIC          DCONJG, MAX
C     ..
C     .. Executable Statements ..
C
C     Test the input parameters.
C
      INFO = 0
      IF     ( .NOT.(UPLO .EQ.'U' .OR. UPLO .EQ.'u').AND.
     $         .NOT.(UPLO .EQ.'L' .OR. UPLO .EQ.'l')      )THEN
         INFO = 1
      ELSE IF( .NOT.(TRANS.EQ.'N' .OR. TRANS.EQ.'n').AND.
     $         .NOT.(TRANS.EQ.'T' .OR. TRANS.EQ.'t').AND.
     $         .NOT.(TRANS.EQ.'C' .OR. TRANS.EQ.'c')      )THEN
         INFO = 2
      ELSE IF( .NOT.(DIAG .EQ.'U' .OR. DIAG .EQ.'u').AND.
     $         .NOT.(DIAG .EQ.'N' .OR. DIAG .EQ.'n')      )THEN
         INFO = 3
      ELSE IF( N.LT.0 )THEN
         INFO = 4
      ELSE IF( LDA.LT.MAX( 1, N ) )THEN
         INFO = 6
      ELSE IF( INCX.EQ.0 )THEN
         INFO = 8
      END IF
      IF( INFO.NE.0 )THEN
         CALL F06AAZ( 'F06SFF/ZTRMV ', INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF( N.EQ.0 )
     $   RETURN
C
      NOCONJ = (TRANS.EQ.'T' .OR. TRANS.EQ.'t')
      NOUNIT = (DIAG .EQ.'N' .OR. DIAG .EQ.'n')
C
C     Set up the start point in X if the increment is not unity. This
C     will be  ( N - 1 )*INCX  too small for descending loops.
C
      IF( INCX.LE.0 )THEN
         KX = 1 - ( N - 1 )*INCX
      ELSE IF( INCX.NE.1 )THEN
         KX = 1
      END IF
C
C     Start the operations. In this version the elements of A are
C     accessed sequentially with one pass through A.
C
      IF( (TRANS.EQ.'N' .OR. TRANS.EQ.'n') )THEN
C
C        Form  x := A*x.
C
         IF( (UPLO.EQ.'U' .OR. UPLO.EQ.'u') )THEN
            IF( INCX.EQ.1 )THEN
               DO 20, J = 1, N
                  IF( X( J ).NE.ZERO )THEN
                     TEMP = X( J )
                     DO 10, I = 1, J - 1
                        X( I ) = X( I ) + TEMP*A( I, J )
   10                CONTINUE
                     IF( NOUNIT )
     $                  X( J ) = X( J )*A( J, J )
                  END IF
   20          CONTINUE
            ELSE
               JX = KX
               DO 40, J = 1, N
                  IF( X( JX ).NE.ZERO )THEN
                     TEMP = X( JX )
                     IX   = KX
                     DO 30, I = 1, J - 1
                        X( IX ) = X( IX ) + TEMP*A( I, J )
                        IX      = IX      + INCX
   30                CONTINUE
                     IF( NOUNIT )
     $                  X( JX ) = X( JX )*A( J, J )
                  END IF
                  JX = JX + INCX
   40          CONTINUE
            END IF
         ELSE
            IF( INCX.EQ.1 )THEN
               DO 60, J = N, 1, -1
                  IF( X( J ).NE.ZERO )THEN
                     TEMP = X( J )
                     DO 50, I = N, J + 1, -1
                        X( I ) = X( I ) + TEMP*A( I, J )
   50                CONTINUE
                     IF( NOUNIT )
     $                  X( J ) = X( J )*A( J, J )
                  END IF
   60          CONTINUE
            ELSE
               KX = KX + ( N - 1 )*INCX
               JX = KX
               DO 80, J = N, 1, -1
                  IF( X( JX ).NE.ZERO )THEN
                     TEMP = X( JX )
                     IX   = KX
                     DO 70, I = N, J + 1, -1
                        X( IX ) = X( IX ) + TEMP*A( I, J )
                        IX      = IX      - INCX
   70                CONTINUE
                     IF( NOUNIT )
     $                  X( JX ) = X( JX )*A( J, J )
                  END IF
                  JX = JX - INCX
   80          CONTINUE
            END IF
         END IF
      ELSE
C
C        Form  x := A'*x  or  x := conjg( A' )*x.
C
         IF( (UPLO.EQ.'U' .OR. UPLO.EQ.'u') )THEN
            IF( INCX.EQ.1 )THEN
               DO 110, J = N, 1, -1
                  TEMP = X( J )
                  IF( NOCONJ )THEN
                     IF( NOUNIT )
     $                  TEMP = TEMP*A( J, J )
                     DO 90, I = J - 1, 1, -1
                        TEMP = TEMP + A( I, J )*X( I )
   90                CONTINUE
                  ELSE
                     IF( NOUNIT )
     $                  TEMP = TEMP*DCONJG( A( J, J ) )
                     DO 100, I = J - 1, 1, -1
                        TEMP = TEMP + DCONJG( A( I, J ) )*X( I )
  100                CONTINUE
                  END IF
                  X( J ) = TEMP
  110          CONTINUE
            ELSE
               JX = KX + ( N - 1 )*INCX
               DO 140, J = N, 1, -1
                  TEMP = X( JX )
                  IX   = JX
                  IF( NOCONJ )THEN
                     IF( NOUNIT )
     $                  TEMP = TEMP*A( J, J )
                     DO 120, I = J - 1, 1, -1
                        IX   = IX   - INCX
                        TEMP = TEMP + A( I, J )*X( IX )
  120                CONTINUE
                  ELSE
                     IF( NOUNIT )
     $                  TEMP = TEMP*DCONJG( A( J, J ) )
                     DO 130, I = J - 1, 1, -1
                        IX   = IX   - INCX
                        TEMP = TEMP + DCONJG( A( I, J ) )*X( IX )
  130                CONTINUE
                  END IF
                  X( JX ) = TEMP
                  JX      = JX   - INCX
  140          CONTINUE
            END IF
         ELSE
            IF( INCX.EQ.1 )THEN
               DO 170, J = 1, N
                  TEMP = X( J )
                  IF( NOCONJ )THEN
                     IF( NOUNIT )
     $                  TEMP = TEMP*A( J, J )
                     DO 150, I = J + 1, N
                        TEMP = TEMP + A( I, J )*X( I )
  150                CONTINUE
                  ELSE
                     IF( NOUNIT )
     $                  TEMP = TEMP*DCONJG( A( J, J ) )
                     DO 160, I = J + 1, N
                        TEMP = TEMP + DCONJG( A( I, J ) )*X( I )
  160                CONTINUE
                  END IF
                  X( J ) = TEMP
  170          CONTINUE
            ELSE
               JX = KX
               DO 200, J = 1, N
                  TEMP = X( JX )
                  IX   = JX
                  IF( NOCONJ )THEN
                     IF( NOUNIT )
     $                  TEMP = TEMP*A( J, J )
                     DO 180, I = J + 1, N
                        IX   = IX   + INCX
                        TEMP = TEMP + A( I, J )*X( IX )
  180                CONTINUE
                  ELSE
                     IF( NOUNIT )
     $                  TEMP = TEMP*DCONJG( A( J, J ) )
                     DO 190, I = J + 1, N
                        IX   = IX   + INCX
                        TEMP = TEMP + DCONJG( A( I, J ) )*X( IX )
  190                CONTINUE
                  END IF
                  X( JX ) = TEMP
                  JX      = JX   + INCX
  200          CONTINUE
            END IF
         END IF
      END IF
C
      RETURN
C
C     End of F06SFF (ZTRMV ).
C
      END
*
      SUBROUTINE F06SNF( M, N, ALPHA, X, INCX, Y, INCY, A, LDA )
C     MARK 13 RE-ISSUE. NAG COPYRIGHT 1988.
C     .. Entry Points ..
      ENTRY      ZGERC ( M, N, ALPHA, X, INCX, Y, INCY, A, LDA )
C     .. Scalar Arguments ..
      COMPLEX*16         ALPHA
      INTEGER            INCX, INCY, LDA, M, N
C     .. Array Arguments ..
      COMPLEX*16         A( LDA, * ), X( * ), Y( * )
C     ..
C
C  Purpose
C  =======
C
C  ZGERC  performs the rank 1 operation
C
C     A := alpha*x*conjg( y' ) + A,
C
C  where alpha is a scalar, x is an m element vector, y is an n element
C  vector and A is an m by n matrix.
C
C  Parameters
C  ==========
C
C  M      - INTEGER.
C           On entry, M specifies the number of rows of the matrix A.
C           M must be at least zero.
C           Unchanged on exit.
C
C  N      - INTEGER.
C           On entry, N specifies the number of columns of the matrix A.
C           N must be at least zero.
C           Unchanged on exit.
C
C  ALPHA  - COMPLEX*16      .
C           On entry, ALPHA specifies the scalar alpha.
C           Unchanged on exit.
C
C  X      - COMPLEX*16       array of dimension at least
C           ( 1 + ( m - 1 )*abs( INCX ) ).
C           Before entry, the incremented array X must contain the m
C           element vector x.
C           Unchanged on exit.
C
C  INCX   - INTEGER.
C           On entry, INCX specifies the increment for the elements of
C           X. INCX must not be zero.
C           Unchanged on exit.
C
C  Y      - COMPLEX*16       array of dimension at least
C           ( 1 + ( n - 1 )*abs( INCY ) ).
C           Before entry, the incremented array Y must contain the n
C           element vector y.
C           Unchanged on exit.
C
C  INCY   - INTEGER.
C           On entry, INCY specifies the increment for the elements of
C           Y. INCY must not be zero.
C           Unchanged on exit.
C
C  A      - COMPLEX*16       array of DIMENSION ( LDA, n ).
C           Before entry, the leading m by n part of the array A must
C           contain the matrix of coefficients. On exit, A is
C           overwritten by the updated matrix.
C
C  LDA    - INTEGER.
C           On entry, LDA specifies the first dimension of A as declared
C           in the calling (sub) program. LDA must be at least
C           max( 1, m ).
C           Unchanged on exit.
C
C
C  Level 2 Blas routine.
C
C  -- Written on 22-October-1986.
C     Jack Dongarra, Argonne National Lab.
C     Jeremy Du Croz, Nag Central Office.
C     Sven Hammarling, Nag Central Office.
C     Richard Hanson, Sandia National Labs.
C
C
C     .. Parameters ..
      COMPLEX*16         ZERO
      PARAMETER        ( ZERO = ( 0.0D+0, 0.0D+0 ) )
C     .. Local Scalars ..
      COMPLEX*16         TEMP
      INTEGER            I, INFO, IX, J, JY, KX
C     .. External Subroutines ..
      EXTERNAL           F06AAZ
C     .. Intrinsic Functions ..
      INTRINSIC          DCONJG, MAX
C     ..
C     .. Executable Statements ..
C
C     Test the input parameters.
C
      INFO = 0
      IF     ( M.LT.0 )THEN
         INFO = 1
      ELSE IF( N.LT.0 )THEN
         INFO = 2
      ELSE IF( INCX.EQ.0 )THEN
         INFO = 5
      ELSE IF( INCY.EQ.0 )THEN
         INFO = 7
      ELSE IF( LDA.LT.MAX( 1, M ) )THEN
         INFO = 9
      END IF
      IF( INFO.NE.0 )THEN
         CALL F06AAZ( 'F06SNF/ZGERC ', INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF( ( M.EQ.0 ).OR.( N.EQ.0 ).OR.( ALPHA.EQ.ZERO ) )
     $   RETURN
C
C     Start the operations. In this version the elements of A are
C     accessed sequentially with one pass through A.
C
      IF( INCY.GT.0 )THEN
         JY = 1
      ELSE
         JY = 1 - ( N - 1 )*INCY
      END IF
      IF( INCX.EQ.1 )THEN
         DO 20, J = 1, N
            IF( Y( JY ).NE.ZERO )THEN
               TEMP = ALPHA*DCONJG( Y( JY ) )
               DO 10, I = 1, M
                  A( I, J ) = A( I, J ) + X( I )*TEMP
   10          CONTINUE
            END IF
            JY = JY + INCY
   20    CONTINUE
      ELSE
         IF( INCX.GT.0 )THEN
            KX = 1
         ELSE
            KX = 1 - ( M - 1 )*INCX
         END IF
         DO 40, J = 1, N
            IF( Y( JY ).NE.ZERO )THEN
               TEMP = ALPHA*DCONJG( Y( JY ) )
               IX   = KX
               DO 30, I = 1, M
                  A( I, J ) = A( I, J ) + X( IX )*TEMP
                  IX        = IX        + INCX
   30          CONTINUE
            END IF
            JY = JY + INCY
   40    CONTINUE
      END IF
C
      RETURN
C
C     End of F06SNF (ZGERC ).
C
      END
*
      SUBROUTINE F06ZAF(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,
     *                  LDC)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C
C  Purpose
C  =======
C
C  ZGEMM  performs one of the matrix-matrix operations
C
C     C := alpha*op( A )*op( B ) + beta*C,
C
C  where  op( X ) is one of
C
C     op( X ) = X   or   op( X ) = X'   or   op( X ) = conjg( X' ),
C
C  alpha and beta are scalars, and A, B and C are matrices, with op( A )
C  an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix.
C
C  Parameters
C  ==========
C
C  TRANSA - CHARACTER*1.
C           On entry, TRANSA specifies the form of op( A ) to be used in
C           the matrix multiplication as follows:
C
C              TRANSA = 'N' or 'n',  op( A ) = A.
C
C              TRANSA = 'T' or 't',  op( A ) = A'.
C
C              TRANSA = 'C' or 'c',  op( A ) = conjg( A' ).
C
C           Unchanged on exit.
C
C  TRANSB - CHARACTER*1.
C           On entry, TRANSB specifies the form of op( B ) to be used in
C           the matrix multiplication as follows:
C
C              TRANSB = 'N' or 'n',  op( B ) = B.
C
C              TRANSB = 'T' or 't',  op( B ) = B'.
C
C              TRANSB = 'C' or 'c',  op( B ) = conjg( B' ).
C
C           Unchanged on exit.
C
C  M      - INTEGER.
C           On entry,  M  specifies  the number  of rows  of the  matrix
C           op( A )  and of the  matrix  C.  M  must  be at least  zero.
C           Unchanged on exit.
C
C  N      - INTEGER.
C           On entry,  N  specifies the number  of columns of the matrix
C           op( B ) and the number of columns of the matrix C. N must be
C           at least zero.
C           Unchanged on exit.
C
C  K      - INTEGER.
C           On entry,  K  specifies  the number of columns of the matrix
C           op( A ) and the number of rows of the matrix op( B ). K must
C           be at least  zero.
C           Unchanged on exit.
C
C  ALPHA  - COMPLEX         .
C           On entry, ALPHA specifies the scalar alpha.
C           Unchanged on exit.
C
C  A      - COMPLEX          array of DIMENSION ( LDA, ka ), where ka is
C           k  when  TRANSA = 'N' or 'n',  and is  m  otherwise.
C           Before entry with  TRANSA = 'N' or 'n',  the leading  m by k
C           part of the array  A  must contain the matrix  A,  otherwise
C           the leading  k by m  part of the array  A  must contain  the
C           matrix A.
C           Unchanged on exit.
C
C  LDA    - INTEGER.
C           On entry, LDA specifies the first dimension of A as declared
C           in the calling (sub) program. When  TRANSA = 'N' or 'n' then
C           LDA must be at least  max( 1, m ), otherwise  LDA must be at
C           least  max( 1, k ).
C           Unchanged on exit.
C
C  B      - COMPLEX          array of DIMENSION ( LDB, kb ), where kb is
C           n  when  TRANSB = 'N' or 'n',  and is  k  otherwise.
C           Before entry with  TRANSB = 'N' or 'n',  the leading  k by n
C           part of the array  B  must contain the matrix  B,  otherwise
C           the leading  n by k  part of the array  B  must contain  the
C           matrix B.
C           Unchanged on exit.
C
C  LDB    - INTEGER.
C           On entry, LDB specifies the first dimension of B as declared
C           in the calling (sub) program. When  TRANSB = 'N' or 'n' then
C           LDB must be at least  max( 1, k ), otherwise  LDB must be at
C           least  max( 1, n ).
C           Unchanged on exit.
C
C  BETA   - COMPLEX         .
C           On entry,  BETA  specifies the scalar  beta.  When  BETA  is
C           supplied as zero then C need not be set on input.
C           Unchanged on exit.
C
C  C      - COMPLEX          array of DIMENSION ( LDC, n ).
C           Before entry, the leading  m by n  part of the array  C must
C           contain the matrix  C,  except when  beta  is zero, in which
C           case C need not be set on entry.
C           On exit, the array  C  is overwritten by the  m by n  matrix
C           ( alpha*op( A )*op( B ) + beta*C ).
C
C  LDC    - INTEGER.
C           On entry, LDC specifies the first dimension of C as declared
C           in  the  calling  (sub)  program.   LDC  must  be  at  least
C           max( 1, m ).
C           Unchanged on exit.
C
C
C  Level 3 Blas routine.
C
C  -- Written on 8-February-1989.
C     Jack Dongarra, Argonne National Laboratory.
C     Iain Duff, AERE Harwell.
C     Jeremy Du Croz, Numerical Algorithms Group Ltd.
C     Sven Hammarling, Numerical Algorithms Group Ltd.
C
C
C     .. Entry Points ..
      ENTRY             ZGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,
     *                  BETA,C,LDC)
C     .. Parameters ..
      COMPLEX*16        ONE
      PARAMETER         (ONE=(1.0D+0,0.0D+0))
      COMPLEX*16        ZERO
      PARAMETER         (ZERO=(0.0D+0,0.0D+0))
C     .. Scalar Arguments ..
      COMPLEX*16        ALPHA, BETA
      INTEGER           K, LDA, LDB, LDC, M, N
      CHARACTER*1       TRANSA, TRANSB
C     .. Array Arguments ..
      COMPLEX*16        A(LDA,*), B(LDB,*), C(LDC,*)
C     .. Local Scalars ..
      COMPLEX*16        TEMP
      INTEGER           I, INFO, J, L, NCOLA, NROWA, NROWB
      LOGICAL           CONJA, CONJB, NOTA, NOTB
C     .. External Subroutines ..
      EXTERNAL          F06AAZ
C     .. Intrinsic Functions ..
      INTRINSIC         DCONJG, MAX
C     .. Executable Statements ..
C
C     Set  NOTA  and  NOTB  as  true if  A  and  B  respectively are not
C     conjugated or transposed, set  CONJA and CONJB  as true if  A  and
C     B  respectively are to be  transposed but  not conjugated  and set
C     NROWA, NCOLA and  NROWB  as the number of rows and  columns  of  A
C     and the number of rows of  B  respectively.
C
      NOTA = (TRANSA.EQ.'N' .OR. TRANSA.EQ.'n')
      NOTB = (TRANSB.EQ.'N' .OR. TRANSB.EQ.'n')
      CONJA = (TRANSA.EQ.'C' .OR. TRANSA.EQ.'c')
      CONJB = (TRANSB.EQ.'C' .OR. TRANSB.EQ.'c')
      IF (NOTA) THEN
         NROWA = M
         NCOLA = K
      ELSE
         NROWA = K
         NCOLA = M
      END IF
      IF (NOTB) THEN
         NROWB = K
      ELSE
         NROWB = N
      END IF
C
C     Test the input parameters.
C
      INFO = 0
      IF (( .NOT. NOTA) .AND. ( .NOT. CONJA)
     *    .AND. ( .NOT. (TRANSA.EQ.'T' .OR. TRANSA.EQ.'t'))) THEN
         INFO = 1
      ELSE IF (( .NOT. NOTB) .AND. ( .NOT. CONJB)
     *         .AND. ( .NOT. (TRANSB.EQ.'T' .OR. TRANSB.EQ.'t'))) THEN
         INFO = 2
      ELSE IF (M.LT.0) THEN
         INFO = 3
      ELSE IF (N.LT.0) THEN
         INFO = 4
      ELSE IF (K.LT.0) THEN
         INFO = 5
      ELSE IF (LDA.LT.MAX(1,NROWA)) THEN
         INFO = 8
      ELSE IF (LDB.LT.MAX(1,NROWB)) THEN
         INFO = 10
      ELSE IF (LDC.LT.MAX(1,M)) THEN
         INFO = 13
      END IF
      IF (INFO.NE.0) THEN
         CALL F06AAZ('F06ZAF/ZGEMM ',INFO)
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF ((M.EQ.0) .OR. (N.EQ.0) .OR. (((ALPHA.EQ.ZERO) .OR. (K.EQ.0))
     *     .AND. (BETA.EQ.ONE))) RETURN
C
C     And when  alpha.eq.zero.
C
      IF (ALPHA.EQ.ZERO) THEN
         IF (BETA.EQ.ZERO) THEN
            DO 40 J = 1, N
               DO 20 I = 1, M
                  C(I,J) = ZERO
   20          CONTINUE
   40       CONTINUE
         ELSE
            DO 80 J = 1, N
               DO 60 I = 1, M
                  C(I,J) = BETA*C(I,J)
   60          CONTINUE
   80       CONTINUE
         END IF
         RETURN
      END IF
C
C     Start the operations.
C
      IF (NOTB) THEN
         IF (NOTA) THEN
C
C           Form  C := alpha*A*B + beta*C.
C
            DO 180 J = 1, N
               IF (BETA.EQ.ZERO) THEN
                  DO 100 I = 1, M
                     C(I,J) = ZERO
  100             CONTINUE
               ELSE IF (BETA.NE.ONE) THEN
                  DO 120 I = 1, M
                     C(I,J) = BETA*C(I,J)
  120             CONTINUE
               END IF
               DO 160 L = 1, K
                  IF (B(L,J).NE.ZERO) THEN
                     TEMP = ALPHA*B(L,J)
                     DO 140 I = 1, M
                        C(I,J) = C(I,J) + TEMP*A(I,L)
  140                CONTINUE
                  END IF
  160          CONTINUE
  180       CONTINUE
         ELSE IF (CONJA) THEN
C
C           Form  C := alpha*conjg( A' )*B + beta*C.
C
            DO 240 J = 1, N
               DO 220 I = 1, M
                  TEMP = ZERO
                  DO 200 L = 1, K
                     TEMP = TEMP + DCONJG(A(L,I))*B(L,J)
  200             CONTINUE
                  IF (BETA.EQ.ZERO) THEN
                     C(I,J) = ALPHA*TEMP
                  ELSE
                     C(I,J) = ALPHA*TEMP + BETA*C(I,J)
                  END IF
  220          CONTINUE
  240       CONTINUE
         ELSE
C
C           Form  C := alpha*A'*B + beta*C
C
            DO 300 J = 1, N
               DO 280 I = 1, M
                  TEMP = ZERO
                  DO 260 L = 1, K
                     TEMP = TEMP + A(L,I)*B(L,J)
  260             CONTINUE
                  IF (BETA.EQ.ZERO) THEN
                     C(I,J) = ALPHA*TEMP
                  ELSE
                     C(I,J) = ALPHA*TEMP + BETA*C(I,J)
                  END IF
  280          CONTINUE
  300       CONTINUE
         END IF
      ELSE IF (NOTA) THEN
         IF (CONJB) THEN
C
C           Form  C := alpha*A*conjg( B' ) + beta*C.
C
            DO 400 J = 1, N
               IF (BETA.EQ.ZERO) THEN
                  DO 320 I = 1, M
                     C(I,J) = ZERO
  320             CONTINUE
               ELSE IF (BETA.NE.ONE) THEN
                  DO 340 I = 1, M
                     C(I,J) = BETA*C(I,J)
  340             CONTINUE
               END IF
               DO 380 L = 1, K
                  IF (B(J,L).NE.ZERO) THEN
                     TEMP = ALPHA*DCONJG(B(J,L))
                     DO 360 I = 1, M
                        C(I,J) = C(I,J) + TEMP*A(I,L)
  360                CONTINUE
                  END IF
  380          CONTINUE
  400       CONTINUE
         ELSE
C
C           Form  C := alpha*A*B'          + beta*C
C
            DO 500 J = 1, N
               IF (BETA.EQ.ZERO) THEN
                  DO 420 I = 1, M
                     C(I,J) = ZERO
  420             CONTINUE
               ELSE IF (BETA.NE.ONE) THEN
                  DO 440 I = 1, M
                     C(I,J) = BETA*C(I,J)
  440             CONTINUE
               END IF
               DO 480 L = 1, K
                  IF (B(J,L).NE.ZERO) THEN
                     TEMP = ALPHA*B(J,L)
                     DO 460 I = 1, M
                        C(I,J) = C(I,J) + TEMP*A(I,L)
  460                CONTINUE
                  END IF
  480          CONTINUE
  500       CONTINUE
         END IF
      ELSE IF (CONJA) THEN
         IF (CONJB) THEN
C
C           Form  C := alpha*conjg( A' )*conjg( B' ) + beta*C.
C
            DO 560 J = 1, N
               DO 540 I = 1, M
                  TEMP = ZERO
                  DO 520 L = 1, K
                     TEMP = TEMP + DCONJG(A(L,I))*DCONJG(B(J,L))
  520             CONTINUE
                  IF (BETA.EQ.ZERO) THEN
                     C(I,J) = ALPHA*TEMP
                  ELSE
                     C(I,J) = ALPHA*TEMP + BETA*C(I,J)
                  END IF
  540          CONTINUE
  560       CONTINUE
         ELSE
C
C           Form  C := alpha*conjg( A' )*B' + beta*C
C
            DO 620 J = 1, N
               DO 600 I = 1, M
                  TEMP = ZERO
                  DO 580 L = 1, K
                     TEMP = TEMP + DCONJG(A(L,I))*B(J,L)
  580             CONTINUE
                  IF (BETA.EQ.ZERO) THEN
                     C(I,J) = ALPHA*TEMP
                  ELSE
                     C(I,J) = ALPHA*TEMP + BETA*C(I,J)
                  END IF
  600          CONTINUE
  620       CONTINUE
         END IF
      ELSE
         IF (CONJB) THEN
C
C           Form  C := alpha*A'*conjg( B' ) + beta*C
C
            DO 680 J = 1, N
               DO 660 I = 1, M
                  TEMP = ZERO
                  DO 640 L = 1, K
                     TEMP = TEMP + A(L,I)*DCONJG(B(J,L))
  640             CONTINUE
                  IF (BETA.EQ.ZERO) THEN
                     C(I,J) = ALPHA*TEMP
                  ELSE
                     C(I,J) = ALPHA*TEMP + BETA*C(I,J)
                  END IF
  660          CONTINUE
  680       CONTINUE
         ELSE
C
C           Form  C := alpha*A'*B' + beta*C
C
            DO 740 J = 1, N
               DO 720 I = 1, M
                  TEMP = ZERO
                  DO 700 L = 1, K
                     TEMP = TEMP + A(L,I)*B(J,L)
  700             CONTINUE
                  IF (BETA.EQ.ZERO) THEN
                     C(I,J) = ALPHA*TEMP
                  ELSE
                     C(I,J) = ALPHA*TEMP + BETA*C(I,J)
                  END IF
  720          CONTINUE
  740       CONTINUE
         END IF
      END IF
C
      RETURN
C
C     End of F06ZAF (ZGEMM ).
C
      END
*
      SUBROUTINE F06ZFF(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C
C  Purpose
C  =======
C
C  ZTRMM  performs one of the matrix-matrix operations
C
C     B := alpha*op( A )*B,   or   B := alpha*B*op( A )
C
C  where  alpha  is a scalar,  B  is an m by n matrix,  A  is a unit, or
C  non-unit,  upper or lower triangular matrix  and  op( A )  is one  of
C
C     op( A ) = A   or   op( A ) = A'   or   op( A ) = conjg( A' ).
C
C  Parameters
C  ==========
C
C  SIDE   - CHARACTER*1.
C           On entry,  SIDE specifies whether  op( A ) multiplies B from
C           the left or right as follows:
C
C              SIDE = 'L' or 'l'   B := alpha*op( A )*B.
C
C              SIDE = 'R' or 'r'   B := alpha*B*op( A ).
C
C           Unchanged on exit.
C
C  UPLO   - CHARACTER*1.
C           On entry, UPLO specifies whether the matrix A is an upper or
C           lower triangular matrix as follows:
C
C              UPLO = 'U' or 'u'   A is an upper triangular matrix.
C
C              UPLO = 'L' or 'l'   A is a lower triangular matrix.
C
C           Unchanged on exit.
C
C  TRANSA - CHARACTER*1.
C           On entry, TRANSA specifies the form of op( A ) to be used in
C           the matrix multiplication as follows:
C
C              TRANSA = 'N' or 'n'   op( A ) = A.
C
C              TRANSA = 'T' or 't'   op( A ) = A'.
C
C              TRANSA = 'C' or 'c'   op( A ) = conjg( A' ).
C
C           Unchanged on exit.
C
C  DIAG   - CHARACTER*1.
C           On entry, DIAG specifies whether or not A is unit triangular
C           as follows:
C
C              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
C
C              DIAG = 'N' or 'n'   A is not assumed to be unit
C                                  triangular.
C
C           Unchanged on exit.
C
C  M      - INTEGER.
C           On entry, M specifies the number of rows of B. M must be at
C           least zero.
C           Unchanged on exit.
C
C  N      - INTEGER.
C           On entry, N specifies the number of columns of B.  N must be
C           at least zero.
C           Unchanged on exit.
C
C  ALPHA  - COMPLEX         .
C           On entry,  ALPHA specifies the scalar  alpha. When  alpha is
C           zero then  A is not referenced and  B need not be set before
C           entry.
C           Unchanged on exit.
C
C  A      - COMPLEX          array of DIMENSION ( LDA, k ), where k is m
C           when  SIDE = 'L' or 'l'  and is  n  when  SIDE = 'R' or 'r'.
C           Before entry  with  UPLO = 'U' or 'u',  the  leading  k by k
C           upper triangular part of the array  A must contain the upper
C           triangular matrix  and the strictly lower triangular part of
C           A is not referenced.
C           Before entry  with  UPLO = 'L' or 'l',  the  leading  k by k
C           lower triangular part of the array  A must contain the lower
C           triangular matrix  and the strictly upper triangular part of
C           A is not referenced.
C           Note that when  DIAG = 'U' or 'u',  the diagonal elements of
C           A  are not referenced either,  but are assumed to be  unity.
C           Unchanged on exit.
C
C  LDA    - INTEGER.
C           On entry, LDA specifies the first dimension of A as declared
C           in the calling (sub) program.  When  SIDE = 'L' or 'l'  then
C           LDA  must be at least  max( 1, m ),  when  SIDE = 'R' or 'r'
C           then LDA must be at least max( 1, n ).
C           Unchanged on exit.
C
C  B      - COMPLEX          array of DIMENSION ( LDB, n ).
C           Before entry,  the leading  m by n part of the array  B must
C           contain the matrix  B,  and  on exit  is overwritten  by the
C           transformed matrix.
C
C  LDB    - INTEGER.
C           On entry, LDB specifies the first dimension of B as declared
C           in  the  calling  (sub)  program.   LDB  must  be  at  least
C           max( 1, m ).
C           Unchanged on exit.
C
C
C  Level 3 Blas routine.
C
C  -- Written on 8-February-1989.
C     Jack Dongarra, Argonne National Laboratory.
C     Iain Duff, AERE Harwell.
C     Jeremy Du Croz, Numerical Algorithms Group Ltd.
C     Sven Hammarling, Numerical Algorithms Group Ltd.
C
C
C     .. Entry Points ..
      ENTRY             ZTRMM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,
     *                  LDB)
C     .. Parameters ..
      COMPLEX*16        ONE
      PARAMETER         (ONE=(1.0D+0,0.0D+0))
      COMPLEX*16        ZERO
      PARAMETER         (ZERO=(0.0D+0,0.0D+0))
C     .. Scalar Arguments ..
      COMPLEX*16        ALPHA
      INTEGER           LDA, LDB, M, N
      CHARACTER*1       DIAG, SIDE, TRANSA, UPLO
C     .. Array Arguments ..
      COMPLEX*16        A(LDA,*), B(LDB,*)
C     .. Local Scalars ..
      COMPLEX*16        TEMP
      INTEGER           I, INFO, J, K, NROWA
      LOGICAL           LSIDE, NOCONJ, NOUNIT, UPPER
C     .. External Subroutines ..
      EXTERNAL          F06AAZ
C     .. Intrinsic Functions ..
      INTRINSIC         DCONJG, MAX
C     .. Executable Statements ..
C
C     Test the input parameters.
C
      LSIDE = (SIDE.EQ.'L' .OR. SIDE.EQ.'l')
      IF (LSIDE) THEN
         NROWA = M
      ELSE
         NROWA = N
      END IF
      NOCONJ = (TRANSA.EQ.'T' .OR. TRANSA.EQ.'t')
      NOUNIT = (DIAG.EQ.'N' .OR. DIAG.EQ.'n')
      UPPER = (UPLO.EQ.'U' .OR. UPLO.EQ.'u')
C
      INFO = 0
      IF (( .NOT. LSIDE) .AND. ( .NOT. (SIDE.EQ.'R' .OR. SIDE.EQ.'r')))
     *    THEN
         INFO = 1
      ELSE IF (( .NOT. UPPER) .AND. ( .NOT. (UPLO.EQ.'L' .OR. UPLO.EQ.
     *         'l'))) THEN
         INFO = 2
      ELSE IF (( .NOT. (TRANSA.EQ.'N' .OR. TRANSA.EQ.'n'))
     *         .AND. ( .NOT. (TRANSA.EQ.'T' .OR. TRANSA.EQ.'t'))
     *         .AND. ( .NOT. (TRANSA.EQ.'C' .OR. TRANSA.EQ.'c'))) THEN
         INFO = 3
      ELSE IF (( .NOT. (DIAG.EQ.'U' .OR. DIAG.EQ.'u'))
     *         .AND. ( .NOT. (DIAG.EQ.'N' .OR. DIAG.EQ.'n'))) THEN
         INFO = 4
      ELSE IF (M.LT.0) THEN
         INFO = 5
      ELSE IF (N.LT.0) THEN
         INFO = 6
      ELSE IF (LDA.LT.MAX(1,NROWA)) THEN
         INFO = 9
      ELSE IF (LDB.LT.MAX(1,M)) THEN
         INFO = 11
      END IF
      IF (INFO.NE.0) THEN
         CALL F06AAZ('F06ZFF/ZTRMM ',INFO)
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF (N.EQ.0) RETURN
C
C     And when  alpha.eq.zero.
C
      IF (ALPHA.EQ.ZERO) THEN
         DO 40 J = 1, N
            DO 20 I = 1, M
               B(I,J) = ZERO
   20       CONTINUE
   40    CONTINUE
         RETURN
      END IF
C
C     Start the operations.
C
      IF (LSIDE) THEN
         IF ((TRANSA.EQ.'N' .OR. TRANSA.EQ.'n')) THEN
C
C           Form  B := alpha*A*B.
C
            IF (UPPER) THEN
               DO 100 J = 1, N
                  DO 80 K = 1, M
                     IF (B(K,J).NE.ZERO) THEN
                        TEMP = ALPHA*B(K,J)
                        DO 60 I = 1, K - 1
                           B(I,J) = B(I,J) + TEMP*A(I,K)
   60                   CONTINUE
                        IF (NOUNIT) TEMP = TEMP*A(K,K)
                        B(K,J) = TEMP
                     END IF
   80             CONTINUE
  100          CONTINUE
            ELSE
               DO 160 J = 1, N
                  DO 140 K = M, 1, -1
                     IF (B(K,J).NE.ZERO) THEN
                        TEMP = ALPHA*B(K,J)
                        B(K,J) = TEMP
                        IF (NOUNIT) B(K,J) = B(K,J)*A(K,K)
                        DO 120 I = K + 1, M
                           B(I,J) = B(I,J) + TEMP*A(I,K)
  120                   CONTINUE
                     END IF
  140             CONTINUE
  160          CONTINUE
            END IF
         ELSE
C
C           Form  B := alpha*B*A'   or   B := alpha*B*conjg( A' ).
C
            IF (UPPER) THEN
               DO 240 J = 1, N
                  DO 220 I = M, 1, -1
                     TEMP = B(I,J)
                     IF (NOCONJ) THEN
                        IF (NOUNIT) TEMP = TEMP*A(I,I)
                        DO 180 K = 1, I - 1
                           TEMP = TEMP + A(K,I)*B(K,J)
  180                   CONTINUE
                     ELSE
                        IF (NOUNIT) TEMP = TEMP*DCONJG(A(I,I))
                        DO 200 K = 1, I - 1
                           TEMP = TEMP + DCONJG(A(K,I))*B(K,J)
  200                   CONTINUE
                     END IF
                     B(I,J) = ALPHA*TEMP
  220             CONTINUE
  240          CONTINUE
            ELSE
               DO 320 J = 1, N
                  DO 300 I = 1, M
                     TEMP = B(I,J)
                     IF (NOCONJ) THEN
                        IF (NOUNIT) TEMP = TEMP*A(I,I)
                        DO 260 K = I + 1, M
                           TEMP = TEMP + A(K,I)*B(K,J)
  260                   CONTINUE
                     ELSE
                        IF (NOUNIT) TEMP = TEMP*DCONJG(A(I,I))
                        DO 280 K = I + 1, M
                           TEMP = TEMP + DCONJG(A(K,I))*B(K,J)
  280                   CONTINUE
                     END IF
                     B(I,J) = ALPHA*TEMP
  300             CONTINUE
  320          CONTINUE
            END IF
         END IF
      ELSE
         IF ((TRANSA.EQ.'N' .OR. TRANSA.EQ.'n')) THEN
C
C           Form  B := alpha*B*A.
C
            IF (UPPER) THEN
               DO 400 J = N, 1, -1
                  TEMP = ALPHA
                  IF (NOUNIT) TEMP = TEMP*A(J,J)
                  DO 340 I = 1, M
                     B(I,J) = TEMP*B(I,J)
  340             CONTINUE
                  DO 380 K = 1, J - 1
                     IF (A(K,J).NE.ZERO) THEN
                        TEMP = ALPHA*A(K,J)
                        DO 360 I = 1, M
                           B(I,J) = B(I,J) + TEMP*B(I,K)
  360                   CONTINUE
                     END IF
  380             CONTINUE
  400          CONTINUE
            ELSE
               DO 480 J = 1, N
                  TEMP = ALPHA
                  IF (NOUNIT) TEMP = TEMP*A(J,J)
                  DO 420 I = 1, M
                     B(I,J) = TEMP*B(I,J)
  420             CONTINUE
                  DO 460 K = J + 1, N
                     IF (A(K,J).NE.ZERO) THEN
                        TEMP = ALPHA*A(K,J)
                        DO 440 I = 1, M
                           B(I,J) = B(I,J) + TEMP*B(I,K)
  440                   CONTINUE
                     END IF
  460             CONTINUE
  480          CONTINUE
            END IF
         ELSE
C
C           Form  B := alpha*B*A'   or   B := alpha*B*conjg( A' ).
C
            IF (UPPER) THEN
               DO 560 K = 1, N
                  DO 520 J = 1, K - 1
                     IF (A(J,K).NE.ZERO) THEN
                        IF (NOCONJ) THEN
                           TEMP = ALPHA*A(J,K)
                        ELSE
                           TEMP = ALPHA*DCONJG(A(J,K))
                        END IF
                        DO 500 I = 1, M
                           B(I,J) = B(I,J) + TEMP*B(I,K)
  500                   CONTINUE
                     END IF
  520             CONTINUE
                  TEMP = ALPHA
                  IF (NOUNIT) THEN
                     IF (NOCONJ) THEN
                        TEMP = TEMP*A(K,K)
                     ELSE
                        TEMP = TEMP*DCONJG(A(K,K))
                     END IF
                  END IF
                  IF (TEMP.NE.ONE) THEN
                     DO 540 I = 1, M
                        B(I,K) = TEMP*B(I,K)
  540                CONTINUE
                  END IF
  560          CONTINUE
            ELSE
               DO 640 K = N, 1, -1
                  DO 600 J = K + 1, N
                     IF (A(J,K).NE.ZERO) THEN
                        IF (NOCONJ) THEN
                           TEMP = ALPHA*A(J,K)
                        ELSE
                           TEMP = ALPHA*DCONJG(A(J,K))
                        END IF
                        DO 580 I = 1, M
                           B(I,J) = B(I,J) + TEMP*B(I,K)
  580                   CONTINUE
                     END IF
  600             CONTINUE
                  TEMP = ALPHA
                  IF (NOUNIT) THEN
                     IF (NOCONJ) THEN
                        TEMP = TEMP*A(K,K)
                     ELSE
                        TEMP = TEMP*DCONJG(A(K,K))
                     END IF
                  END IF
                  IF (TEMP.NE.ONE) THEN
                     DO 620 I = 1, M
                        B(I,K) = TEMP*B(I,K)
  620                CONTINUE
                  END IF
  640          CONTINUE
            END IF
         END IF
      END IF
C
      RETURN
C
C     End of F06ZFF (ZTRMM ).
C
      END
*
      SUBROUTINE F08ASW(SIDE,M,N,V,INCV,TAU,C,LDC,WORK)
C     MARK 16 RELEASE. NAG COPYRIGHT 1992.
C     ENTRY             ZLARF(SIDE,M,N,V,INCV,TAU,C,LDC,WORK)
C
C  Purpose
C  =======
C
C  ZLARF applies a complex elementary reflector H to a complex m by n
C  matrix C, from either the left or the right. H is represented in the
C  form
C
C        H = I - tau * v * v'
C
C  where tau is a complex scalar and v is a complex vector.
C
C  If tau = 0, then H is taken to be the unit matrix.
C
C  To apply H' (the conjugate transpose of H), supply conjg(tau) instead
C  tau.
C
C  Arguments
C  =========
C
C  SIDE    (input) CHARACTER*1
C          = 'L': form  H * C
C          = 'R': form  C * H
C
C  M       (input) INTEGER
C          The number of rows of the matrix C.
C
C  N       (input) INTEGER
C          The number of columns of the matrix C.
C
C  V       (input) COMPLEX*16 array, dimension
C                     (1 + (M-1)*abs(INCV)) if SIDE = 'L'
C                  or (1 + (N-1)*abs(INCV)) if SIDE = 'R'
C          The vector v in the representation of H. V is not used if
C          TAU = 0.
C
C  INCV    (input) INTEGER
C          The increment between elements of v. INCV <> 0.
C
C  TAU     (input) COMPLEX*16
C          The value tau in the representation of H.
C
C  C       (input/output) COMPLEX*16 array, dimension (LDC,N)
C          On entry, the m-by-n matrix C.
C          On exit, C is overwritten by the matrix H * C if SIDE = 'L',
C          or C * H if SIDE = 'R'.
C
C  LDC     (input) INTEGER
C          The leading dimension of the array C. LDA >= max(1,M).
C
C  WORK    (workspace) COMPLEX*16 array, dimension
C                         (N) if SIDE = 'L'
C                      or (M) if SIDE = 'R'
C
C  -- LAPACK auxiliary routine (adapted for NAG Library)
C     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
C     Courant Institute, Argonne National Lab, and Rice University
C
C  =====================================================================
C
C     .. Parameters ..
      COMPLEX*16        ONE, ZERO
      PARAMETER         (ONE=1.0D+0,ZERO=0.0D+0)
C     .. Scalar Arguments ..
      COMPLEX*16        TAU
      INTEGER           INCV, LDC, M, N
      CHARACTER         SIDE
C     .. Array Arguments ..
      COMPLEX*16        C(LDC,*), V(*), WORK(*)
C     .. External Subroutines ..
      EXTERNAL          ZGEMV, ZGERC
C     .. Executable Statements ..
C
      IF ((SIDE.EQ.'L' .OR. SIDE.EQ.'l')) THEN
C
C        Form  H * C
C
         IF (TAU.NE.ZERO) THEN
C
C           w := C' * v
C
            CALL ZGEMV('Conjugate transpose',M,N,ONE,C,LDC,V,INCV,ZERO,
     *                 WORK,1)
C
C           C := C - v * w'
C
            CALL ZGERC(M,N,-TAU,V,INCV,WORK,1,C,LDC)
         END IF
      ELSE
C
C        Form  C * H
C
         IF (TAU.NE.ZERO) THEN
C
C           w := C * v
C
            CALL ZGEMV('No transpose',M,N,ONE,C,LDC,V,INCV,ZERO,WORK,1)
C
C           C := C - w * v'
C
            CALL ZGERC(M,N,-TAU,WORK,1,V,INCV,C,LDC)
         END IF
      END IF
      RETURN
C
C     End of F08ASW (ZLARF)
C
      END
*
      SUBROUTINE F08ATZ(M,N,K,A,LDA,TAU,WORK,INFO)
C     MARK 16 RELEASE. NAG COPYRIGHT 1992.
C     ENTRY             ZUNG2R(M,N,K,A,LDA,TAU,WORK,INFO)
C
C  Purpose
C  =======
C
C  ZUNG2R generates an m by n complex matrix Q with orthonormal columns,
C  which is defined as the first n columns of a product of k elementary
C  reflectors of order m
C
C        Q  =  H(1) H(2) . . . H(k)
C
C  as returned by ZGEQRF.
C
C  Arguments
C  =========
C
C  M       (input) INTEGER
C          The number of rows of the matrix Q. M >= 0.
C
C  N       (input) INTEGER
C          The number of columns of the matrix Q. M >= N >= 0.
C
C  K       (input) INTEGER
C          The number of elementary reflectors whose product defines the
C          matrix Q. N >= K >= 0.
C
C  A       (input/output) COMPLEX*16 array, dimension (LDA,N)
C          On entry, the i-th column must contain the vector which
C          defines the elementary reflector H(i), for i = 1,2,...,k, as
C          returned by ZGEQRF in the first k columns of its array
C          argument A.
C          On exit, the m by n matrix Q.
C
C  LDA     (input) INTEGER
C          The first dimension of the array A. LDA >= max(1,M).
C
C  TAU     (input) COMPLEX*16 array, dimension (K)
C          TAU(i) must contain the scalar factor of the elementary
C          reflector H(i), as returned by ZGEQRF.
C
C  WORK    (workspace) COMPLEX*16 array, dimension (N)
C
C  INFO    (output) INTEGER
C          = 0: successful exit
C          < 0: if INFO = -i, the i-th argument has an illegal value
C
C  -- LAPACK routine (adapted for NAG Library)
C     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
C     Courant Institute, Argonne National Lab, and Rice University
C
C  =====================================================================
C
C     .. Parameters ..
      COMPLEX*16        ONE, ZERO
      PARAMETER         (ONE=1.0D+0,ZERO=0.0D+0)
C     .. Scalar Arguments ..
      INTEGER           INFO, K, LDA, M, N
C     .. Array Arguments ..
      COMPLEX*16        A(LDA,*), TAU(*), WORK(*)
C     .. Local Scalars ..
      INTEGER           I, J, L
C     .. External Subroutines ..
      EXTERNAL          F06AAZ, F08ASW, ZSCAL
C     .. Intrinsic Functions ..
      INTRINSIC         MAX
C     .. Executable Statements ..
C
C     Test the input arguments
C
      INFO = 0
      IF (M.LT.0) THEN
         INFO = -1
      ELSE IF (N.LT.0 .OR. N.GT.M) THEN
         INFO = -2
      ELSE IF (K.LT.0 .OR. K.GT.N) THEN
         INFO = -3
      ELSE IF (LDA.LT.MAX(1,M)) THEN
         INFO = -5
      END IF
      IF (INFO.NE.0) THEN
         CALL F06AAZ('F08ATZ/ZUNG2R',-INFO)
         RETURN
      END IF
C
C     Quick return if possible
C
      IF (N.LE.0) RETURN
C
C     Initialise columns k+1:n to columns of the unit matrix
C
      DO 40 J = K + 1, N
         DO 20 L = 1, M
            A(L,J) = ZERO
   20    CONTINUE
         A(J,J) = ONE
   40 CONTINUE
C
      DO 80 I = K, 1, -1
C
C        Apply H(i) to A(i:m,i:n) from the left
C
         IF (I.LT.N) THEN
            A(I,I) = ONE
            CALL F08ASW('Left',M-I+1,N-I,A(I,I),1,TAU(I),A(I,I+1),LDA,
     *                  WORK)
         END IF
         IF (I.LT.M) CALL ZSCAL(M-I,-TAU(I),A(I+1,I),1)
         A(I,I) = ONE - TAU(I)
C
C        Set A(1:i-1,i) to zero
C
         DO 60 L = 1, I - 1
            A(L,I) = ZERO
   60    CONTINUE
   80 CONTINUE
      RETURN
C
C     End of F08ATZ (ZUNG2R)
C
      END
*
      SUBROUTINE F08CTY(M,N,K,A,LDA,TAU,WORK,LWORK,INFO)
C     MARK 16 RELEASE. NAG COPYRIGHT 1992.
C     ENTRY             ZUNGQL(M,N,K,A,LDA,TAU,WORK,LWORK,INFO)
C
C  Purpose
C  =======
C
C  ZUNGQL generates an m-by-n complex matrix Q with orthonormal columns,
C  which is defined as the last n columns of a product of k elementary
C  reflectors of order m
C
C        Q  =  H(k) . . . H(2) H(1)
C
C  as returned by ZGEQLF.
C
C  Arguments
C  =========
C
C  M       (input) INTEGER
C          The number of rows of the matrix Q. M >= 0.
C
C  N       (input) INTEGER
C          The number of columns of the matrix Q. M >= N >= 0.
C
C  K       (input) INTEGER
C          The number of elementary reflectors whose product defines the
C          matrix Q. N >= K >= 0.
C
C  A       (input/output) COMPLEX*16 array, dimension (LDA,N)
C          On entry, the (n-k+i)-th column must contain the vector which
C          defines the elementary reflector H(i), for i = 1,2,...,k, as
C          returned by ZGEQLF in the last k columns of its array
C          argument A.
C          On exit, the m by n matrix Q.
C
C  LDA     (input) INTEGER
C          The first dimension of the array A. LDA >= max(1,M).
C
C  TAU     (input) COMPLEX*16 array, dimension (K)
C          TAU(i) must contain the scalar factor of the elementary
C          reflector H(i), as returned by ZGEQLF.
C
C  WORK    (workspace) COMPLEX*16 array, dimension (LWORK)
C          On exit, if INFO = 0, WORK(1) returns the minimum value of
C          LWORK required to use the optimal blocksize.
C
C  LWORK   (input) INTEGER
C          The dimension of the array WORK. LWORK >= max(1,N).
C          For optimum performance LWORK should be at least N*NB, where
C          NB is the optimal blocksize.
C
C  INFO    (output) INTEGER
C          = 0: successful exit;
C          < 0: if INFO = -i, the i-th argument has an illegal value
C
C  -- LAPACK routine (adapted for NAG Library)
C     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
C     Courant Institute, Argonne National Lab, and Rice University
C
C  =====================================================================
C
C     .. Parameters ..
      COMPLEX*16        ZERO
      PARAMETER         (ZERO=0.0D+0)
C     .. Scalar Arguments ..
      INTEGER           INFO, K, LDA, LWORK, M, N
C     .. Array Arguments ..
      COMPLEX*16        A(LDA,*), TAU(*), WORK(LWORK)
C     .. Local Scalars ..
      INTEGER           I, IB, IINFO, IWS, J, KK, L, LDWORK, NB, NBMIN,
     *                  NX
C     .. External Subroutines ..
      EXTERNAL          F06AAZ, F07ZAZ, F08ASX, F08ASY, F08CTZ
C     .. Intrinsic Functions ..
      INTRINSIC         MAX, MIN
C     .. Executable Statements ..
C
C     Test the input arguments
C
      INFO = 0
      IF (M.LT.0) THEN
         INFO = -1
      ELSE IF (N.LT.0 .OR. N.GT.M) THEN
         INFO = -2
      ELSE IF (K.LT.0 .OR. K.GT.N) THEN
         INFO = -3
      ELSE IF (LDA.LT.MAX(1,M)) THEN
         INFO = -5
      ELSE IF (LWORK.LT.MAX(1,N)) THEN
         INFO = -8
      END IF
      IF (INFO.NE.0) THEN
         CALL F06AAZ('F08CTY/ZUNGQL',-INFO)
         RETURN
      END IF
C
C     Quick return if possible
C
      IF (N.LE.0) THEN
         WORK(1) = 1
         RETURN
      END IF
C
C     Determine the block size.
C
      CALL F07ZAZ(1,'F08CTY',NB,0)
      NBMIN = 2
      NX = 0
      IWS = N
      IF (NB.GT.1 .AND. NB.LT.K) THEN
C
C        Determine when to cross over from blocked to unblocked code.
C
         CALL F07ZAZ(3,'F08CTY',NX,0)
         NX = MAX(0,NX)
         IF (NX.LT.K) THEN
C
C           Determine if workspace is large enough for blocked code.
C
            LDWORK = N
            IWS = LDWORK*NB
            IF (LWORK.LT.IWS) THEN
C
C              Not enough workspace to use optimal NB:  reduce NB and
C              determine the minimum value of NB.
C
               NB = LWORK/LDWORK
               CALL F07ZAZ(2,'F08CTY',NBMIN,0)
               NBMIN = MAX(2,NBMIN)
            END IF
         END IF
      END IF
C
      IF (NB.GE.NBMIN .AND. NB.LT.K .AND. NX.LT.K) THEN
C
C        Use blocked code after the first block.
C        The last kk columns are handled by the block method.
C
         KK = MIN(K,((K-NX+NB-1)/NB)*NB)
C
C        Set A(m-kk+1:m,1:n-kk) to zero.
C
         DO 40 J = 1, N - KK
            DO 20 I = M - KK + 1, M
               A(I,J) = ZERO
   20       CONTINUE
   40    CONTINUE
      ELSE
         KK = 0
      END IF
C
C     Use unblocked code for the first or only block.
C
      CALL F08CTZ(M-KK,N-KK,K-KK,A,LDA,TAU,WORK,IINFO)
C
      IF (KK.GT.0) THEN
C
C        Use blocked code
C
         DO 100 I = K - KK + 1, K, NB
            IB = MIN(NB,K-I+1)
            IF (N-K+I.GT.1) THEN
C
C              Form the triangular factor of the block reflector
C              H = H(i+ib-1) . . . H(i+1) H(i)
C
               CALL F08ASX('Backward','Columnwise',M-K+I+IB-1,IB,
     *                     A(1,N-K+I),LDA,TAU(I),WORK,LDWORK)
C
C              Apply H to A(1:m-k+i+ib-1,1:n-k+i-1) from the left
C
               CALL F08ASY('Left','No transpose','Backward',
     *                     'Columnwise',M-K+I+IB-1,N-K+I-1,IB,A(1,N-K+I)
     *                     ,LDA,WORK,LDWORK,A,LDA,WORK(IB+1),LDWORK)
            END IF
C
C           Apply H to rows 1:m-k+i+ib-1 of current block
C
            CALL F08CTZ(M-K+I+IB-1,IB,IB,A(1,N-K+I),LDA,TAU(I),WORK,
     *                  IINFO)
C
C           Set rows m-k+i+ib:m of current block to zero
C
            DO 80 J = N - K + I, N - K + I + IB - 1
               DO 60 L = M - K + I + IB, M
                  A(L,J) = ZERO
   60          CONTINUE
   80       CONTINUE
  100    CONTINUE
      END IF
C
      WORK(1) = IWS
      RETURN
C
C     End of F08CTY (ZUNGQL)
C
      END
*
      SUBROUTINE F06KJF( N, X, INCX, SCALE, SUMSQ )
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     .. Scalar Arguments ..
      DOUBLE PRECISION   SCALE, SUMSQ
      INTEGER            INCX, N
C     .. Array Arguments ..
      COMPLEX*16         X( * )
C     ..
C
C  F06KJF returns the values scl and ssq such that
C
C     ( scl**2 )*ssq = x( 1 )**2 +...+ x( n )**2 + ( scale**2 )*sumsq,
C
C  where x( i ) = abs( X( 1 + ( i - 1 )*INCX ) ). The value of sumsq is
C  assumed to be at least unity and the value of ssq will then satisfy
C
C     1.0 .le. ssq .le. ( sumsq + 2*n ).
C
C  scale is assumed to be non-negative and scl returns the value
C
C     scl = max( scale, abs( real( x( i ) ) ), abs( aimag( x( i ) ) ) ),
C            i
C
C  scale and sumsq must be supplied in SCALE and SUMSQ respectively.
C  SCALE and SUMSQ are overwritten by scl and ssq respectively.
C
C  The routine makes only one pass through the vector X.
C
C
C  Nag Fortran 77 basic linear algebra routine.
C  Nag Fortran 77 O( n ) basic linear algebra routine.
C
C  -- Written on 27-April-1983.
C     Sven Hammarling, Nag Central Office.
C
C
C     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER        ( ZERO = 0.0D+0 )
C     .. Local Scalars ..
      DOUBLE PRECISION   TEMP1
      INTEGER            IX
C     .. Intrinsic Functions ..
      INTRINSIC          ABS, DIMAG, DBLE
C     ..
C     .. Executable Statements ..
      IF( N.GT.0 )THEN
         DO 10, IX = 1, 1 + ( N - 1 )*INCX, INCX
            IF( DBLE( X( IX ) ).NE.ZERO )THEN
               TEMP1 = ABS( DBLE( X( IX ) ) )
               IF( SCALE.LT.TEMP1 )THEN
                  SUMSQ = 1     + SUMSQ*( SCALE/TEMP1 )**2
                  SCALE = TEMP1
               ELSE
                  SUMSQ = SUMSQ +       ( TEMP1/SCALE )**2
               END IF
            END IF
            IF( DIMAG( X( IX ) ).NE.ZERO )THEN
               TEMP1 = ABS( DIMAG( X( IX ) ) )
               IF( SCALE.LT.TEMP1 )THEN
                  SUMSQ = 1     + SUMSQ*( SCALE/TEMP1 )**2
                  SCALE = TEMP1
               ELSE
                  SUMSQ = SUMSQ +       ( TEMP1/SCALE )**2
               END IF
            END IF
   10    CONTINUE
      END IF
C
      RETURN
C
C     End of F06KJF. ( SCSSQ )
C
      END
