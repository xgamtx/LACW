      SUBROUTINE P01ABZ
C     MARK 11.5(F77) RELEASE. NAG COPYRIGHT 1986.
C
C     Terminates execution when a hard failure occurs.
C
C     ******************** IMPLEMENTATION NOTE ********************
C     The following STOP statement may be replaced by a call to an
C     implementation-dependent routine to display a message and/or
C     to abort the program.
C     *************************************************************
C     .. Executable Statements ..
      STOP
      END
*
      DOUBLE PRECISION FUNCTION X02AJF()
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C
C     RETURNS  (1/2)*B**(1-P)  IF ROUNDS IS .TRUE.
C     RETURNS  B**(1-P)  OTHERWISE
C
      DOUBLE PRECISION X02CON
      DATA X02CON /1.11022302462516D-16 /
C     .. Executable Statements ..
      X02AJF = X02CON
      RETURN
      END
*
      DOUBLE PRECISION FUNCTION X02AMF()
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C
C     RETURNS THE 'SAFE RANGE' PARAMETER
C     I.E. THE SMALLEST POSITIVE MODEL NUMBER Z SUCH THAT
C     FOR ANY X WHICH SATISFIES X.GE.Z AND X.LE.1/Z
C     THE FOLLOWING CAN BE COMPUTED WITHOUT OVERFLOW, UNDERFLOW OR OTHER
C     ERROR
C
C        -X
C        1.0/X
C        SQRT(X)
C        LOG(X)
C        EXP(LOG(X))
C        Y**(LOG(X)/LOG(Y)) FOR ANY Y
C
      DOUBLE PRECISION X02CON
      DATA X02CON /2.22507385850721D-308 /
C     .. Executable Statements ..
      X02AMF = X02CON
      RETURN
      END
*
      INTEGER FUNCTION X02BHF()
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C
C     RETURNS THE MODEL PARAMETER, B.
C
C     .. Executable Statements ..
      X02BHF =     2
      RETURN
      END
*
      INTEGER FUNCTION X02BJF()
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C
C     RETURNS THE MODEL PARAMETER, p.
C
C     .. Executable Statements ..
      X02BJF =    53
      RETURN
      END
      INTEGER FUNCTION X02BKF()
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C
C     RETURNS THE MODEL PARAMETER, EMIN.
C
C     .. Executable Statements ..
      X02BKF =  -1021
      RETURN
      END
*
      INTEGER FUNCTION X02BLF()
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C
C     RETURNS THE MODEL PARAMETER, EMAX.
C
C     .. Executable Statements ..
      X02BLF =  1024
      RETURN
      END
*
      SUBROUTINE X04AAF(I,NERR)
C     MARK 7 RELEASE. NAG COPYRIGHT 1978
C     MARK 7C REVISED IER-190 (MAY 1979)
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     MARK 14 REVISED. IER-829 (DEC 1989).
C     IF I = 0, SETS NERR TO CURRENT ERROR MESSAGE UNIT NUMBER
C     (STORED IN NERR1).
C     IF I = 1, CHANGES CURRENT ERROR MESSAGE UNIT NUMBER TO
C     VALUE SPECIFIED BY NERR.
C
C     .. Scalar Arguments ..
      INTEGER           I, NERR
C     .. Local Scalars ..
      INTEGER           NERR1
C     .. Save statement ..
      SAVE              NERR1
C     .. Data statements ..
      DATA              NERR1/0/
C     .. Executable Statements ..
      IF (I.EQ.0) NERR = NERR1
      IF (I.EQ.1) NERR1 = NERR
      RETURN
      END
*
      SUBROUTINE X04ABF(I,NADV)
C     MARK 7 RELEASE. NAG COPYRIGHT 1978
C     MARK 7C REVISED IER-190 (MAY 1979)
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     MARK 14 REVISED. IER-830 (DEC 1989).
C      IF I = 0, SETS NADV TO CURRENT ADVISORY MESSAGE UNIT NUMBER
C     (STORED IN NADV1).
C     IF I = 1, CHANGES CURRENT ADVISORY MESSAGE UNIT NUMBER TO
C     VALUE SPECIFIED BY NADV.
C
C     .. Scalar Arguments ..
      INTEGER           I, NADV
C     .. Local Scalars ..
      INTEGER           NADV1
C     .. Save statement ..
      SAVE              NADV1
C     .. Data statements ..
      DATA              NADV1/6/
C     .. Executable Statements ..
      IF (I.EQ.0) NADV = NADV1
      IF (I.EQ.1) NADV1 = NADV
      RETURN
      END
*
      SUBROUTINE X04BAF(NOUT,REC)
C     MARK 11.5(F77) RELEASE. NAG COPYRIGHT 1986.
C
C     X04BAF writes the contents of REC to the unit defined by NOUT.
C
C     Trailing blanks are not output, except that if REC is entirely
C     blank, a single blank character is output.
C     If NOUT.lt.0, i.e. if NOUT is not a valid Fortran unit identifier,
C     then no output occurs.
C
C     .. Scalar Arguments ..
      INTEGER           NOUT
      CHARACTER*(*)     REC
C     .. Local Scalars ..
      INTEGER           I
C     .. Intrinsic Functions ..
      INTRINSIC         LEN
C     .. Executable Statements ..
      IF (NOUT.GE.0) THEN
C        Remove trailing blanks
         DO 20 I = LEN(REC), 2, -1
            IF (REC(I:I).NE.' ') GO TO 40
   20    CONTINUE
C        Write record to external file
   40    WRITE (NOUT,FMT=99999) REC(1:I)
      END IF
      RETURN
C
99999 FORMAT (A)
      END
*
      SUBROUTINE X04CBZ(STRING,START,FINISH)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C     Returns START as first non blank position of STRING, and FINISH
C     as the last non-blank.
C     .. Scalar Arguments ..
      INTEGER           FINISH, START
      CHARACTER*(*)     STRING
C     .. Local Scalars ..
      INTEGER           L
C     .. Intrinsic Functions ..
      INTRINSIC         LEN
C     .. Executable Statements ..
      IF (STRING.EQ.' ') THEN
         START = 0
         FINISH = 0
      ELSE
         L = LEN(STRING)
         START = 1
   20    IF (STRING(START:START).EQ.' ' .AND. START.LT.L) THEN
            START = START + 1
            GO TO 20
         END IF
         FINISH = L
   40    IF (STRING(FINISH:FINISH).EQ.' ' .AND. FINISH.GT.1) THEN
            FINISH = FINISH - 1
            GO TO 40
         END IF
      END IF
      RETURN
      END
*
      SUBROUTINE X04DBF(MATRIX,DIAG,M,N,A,LDA,USEFRM,FORMAT,TITLE,
     *                  LABROW,RLABS,LABCOL,CLABS,NCOLS,INDENT,IFAIL)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C     MARK 16A REVISED. IER-1049 (JUN 1993).
C     Prints a general complex matrix.
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='X04DBF')
      COMPLEX*16        CZERO
      PARAMETER         (CZERO=(0.0D+0,0.0D+0))
      DOUBLE PRECISION  FRNGLW, FRNGUP, ONE, ZERO
      PARAMETER         (FRNGLW=1.0D-3,FRNGUP=9.99999991D+3,ONE=1.0D+0,
     *                  ZERO=0.0D+0)
C     .. Scalar Arguments ..
      INTEGER           IFAIL, INDENT, LDA, M, N, NCOLS
      CHARACTER         DIAG, LABCOL, LABROW, MATRIX, USEFRM
      CHARACTER*(*)     FORMAT, TITLE
C     .. Array Arguments ..
      COMPLEX*16        A(LDA,*)
      CHARACTER*(*)     CLABS(*), RLABS(*)
C     .. Local Scalars ..
      COMPLEX*16        DIAGEL
      DOUBLE PRECISION  AA, BASE, LBASE, MAXEL, MINEL, RNDIGS
      INTEGER           CL, CLBWID, CLEFT, CRIGHT, CSHIFT, FINIS2,
     *                  FINISH, I, IERR, INCOLS, INDNT, J, K, LWID, ND,
     *                  NDIGS1, NDIGS2, NELEMS, NELS, NOUT, NREC, NSETS,
     *                  NTITLE, NUMWID, OFFSET, RLBWID, START
      LOGICAL           ABOVE, BRACK, GENERL, LOWER, PRDIAG
      CHARACTER*10      ICFORM, IRFORM
      CHARACTER*81      FORMT
      CHARACTER*132     BLANKS, TILTER
      CHARACTER*133     INFBOT, INFILE
      CHARACTER*189     FORM
      CHARACTER*266     INFIL2
C     .. Local Arrays ..
      CHARACTER*80      REC(2)
C     .. External Functions ..
      INTEGER           P01ABF, X02BHF, X02BJF, X02BKF, X02BLF
      EXTERNAL          P01ABF, X02BHF, X02BJF, X02BKF, X02BLF
C     .. External Subroutines ..
      EXTERNAL          X04ABF, X04BAF, X04CBZ
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, DIMAG, LEN, LOG10, MAX, MIN, DBLE
C     .. Executable Statements ..
C
      IERR = 0
      GENERL = MATRIX .EQ. 'G' .OR. MATRIX .EQ. 'g'
      LOWER = MATRIX .EQ. 'L' .OR. MATRIX .EQ. 'l'
      ABOVE = USEFRM .EQ. 'A' .OR. USEFRM .EQ. 'a'
      BRACK = USEFRM .EQ. 'B' .OR. USEFRM .EQ. 'b'
      IF (NCOLS.LE.0 .OR. NCOLS.GT.132) THEN
         INCOLS = 80
      ELSE
         INCOLS = NCOLS
      END IF
      IF (INDENT.LT.0 .OR. INDENT.GE.INCOLS) THEN
         INDNT = 0
      ELSE
         INDNT = INDENT
      END IF
      INCOLS = INCOLS - INDNT
      BLANKS = ' '
C
C     Check for incorrect arguments.
      IF ( .NOT. GENERL .AND. .NOT. LOWER .AND. MATRIX.NE.'U' .AND.
     *    MATRIX.NE.'u') THEN
         IERR = 1
         NREC = 1
         WRITE (REC,FMT=99999) MATRIX
      ELSE IF ((LABROW.NE.'N' .AND. LABROW.NE.'n' .AND. LABROW.NE.
     *         'I' .AND. LABROW.NE.'i' .AND. LABROW.NE.'C' .AND.
     *         LABROW.NE.'c') .OR. (LABCOL.NE.'N' .AND. LABCOL.NE.
     *         'n' .AND. LABCOL.NE.'I' .AND. LABCOL.NE.'i' .AND.
     *         LABCOL.NE.'C' .AND. LABCOL.NE.'c')) THEN
         IERR = 7
         NREC = 2
         WRITE (REC,FMT=99996) LABROW, LABCOL
      ELSE IF ( .NOT. ABOVE .AND. .NOT. BRACK .AND. USEFRM.NE.'D' .AND.
     *         USEFRM.NE.'d') THEN
         IERR = 4
         NREC = 1
         WRITE (REC,FMT=99995) USEFRM
      ELSE IF (M.GT.LDA) THEN
         IERR = 3
         NREC = 1
         WRITE (REC,FMT=99998) M, LDA
      ELSE IF ( .NOT. GENERL) THEN
         IF (DIAG.NE.'U' .AND. DIAG.NE.'u' .AND. DIAG.NE.'N' .AND.
     *       DIAG.NE.'n' .AND. DIAG.NE.'B' .AND. DIAG.NE.'b') THEN
            IERR = 2
            NREC = 1
            WRITE (REC,FMT=99997) MATRIX, DIAG
         END IF
      END IF
      IF (IERR.NE.0) GO TO 560
C
C     Get the advisory message unit number.
      CALL X04ABF(0,NOUT)
C
      FORMT = FORMAT
      IF (FORMT.EQ.'*') THEN
C        Construct an E FORMAT that is wide enough to distinguish
C        between adjacent machine numbers.
         BASE = X02BHF()
         LBASE = LOG10(BASE)
         RNDIGS = X02BJF()*LBASE
         NDIGS1 = RNDIGS
         IF (NDIGS1.LT.RNDIGS) NDIGS1 = NDIGS1 + 1
C        NDIGS1 is the number of significant decimal digits required
C        for the mantissa.
         RNDIGS = LOG10(MAX(-X02BKF(),X02BLF())*LBASE) + 2
         NDIGS2 = RNDIGS
         IF (NDIGS2.LT.RNDIGS) NDIGS2 = NDIGS2 + 1
C        NDIGS2 is the number of decimal places required for the
C        exponent, including letter 'E' and sign.
         FORMT = '1P,E   .   '
         WRITE (FORMT(5:7),FMT='(I3)') NDIGS1 + NDIGS2 + 3
         WRITE (FORMT(9:11),FMT='(I3)') NDIGS1 - 1
      ELSE IF (FORMT.EQ.' ') THEN
C        Construct either a fixed point FORMAT, if the elements to be
C        printed are a reasonable size, i.e. all lie inside the range
C        0.001 to 9999.9999, or a floating-point FORMAT otherwise,
C        printing to 5 significant digits.
C        First find the largest and smallest elements to be printed,
C        ignoring zeros.
         MAXEL = ONE
         MINEL = ONE
         IF (GENERL) THEN
            DO 40 J = 1, N
               DO 20 I = 1, M
                  AA = MAX(ABS(DBLE(A(I,J))),ABS(DIMAG(A(I,J))))
                  IF (AA.GT.MAXEL) MAXEL = AA
                  IF (AA.NE.ZERO .AND. AA.LT.MINEL) MINEL = AA
   20          CONTINUE
   40       CONTINUE
         ELSE
            IF (DIAG.EQ.'N' .OR. DIAG.EQ.'n') THEN
               DO 50 J = 1, MIN(M,N)
                  AA = MAX(ABS(DBLE(A(J,J))),ABS(DIMAG(A(J,J))))
                  IF (AA.GT.MAXEL) MAXEL = AA
                  IF (AA.NE.ZERO .AND. AA.LT.MINEL) MINEL = AA
   50          CONTINUE
            END IF
            DO 100 J = 1, N
               IF (LOWER) THEN
                  DO 60 I = J + 1, M
                     AA = MAX(ABS(DBLE(A(I,J))),ABS(DIMAG(A(I,J))))
                     IF (AA.GT.MAXEL) MAXEL = AA
                     IF (AA.NE.ZERO .AND. AA.LT.MINEL) MINEL = AA
   60             CONTINUE
               ELSE
                  DO 80 I = 1, J - 1
                     AA = MAX(ABS(DBLE(A(I,J))),ABS(DIMAG(A(I,J))))
                     IF (AA.GT.MAXEL) MAXEL = AA
                     IF (AA.NE.ZERO .AND. AA.LT.MINEL) MINEL = AA
   80             CONTINUE
               END IF
  100       CONTINUE
         END IF
C
         IF (MINEL.GE.FRNGLW .AND. MAXEL.LT.FRNGUP) THEN
C           If all elements to be printed are moderately sized,
C           use a fixed point FORMAT ...
            IF (MAXEL.LE.ONE) THEN
               FORMT = 'F8.4'
            ELSE
               FORMT = 'F11.4'
            END IF
         ELSE
C           ... otherwise use a floating-point FORMAT.
            FORMT = '1PE13.4'
         END IF
      END IF
C
C     Construct the format statement to be used internally.
      CALL X04CBZ(FORMT,START,FINISH)
      IF (FINISH-START+1.GT.80) THEN
C        The length of FORMAT is too great.
         IERR = 5
         NREC = 2
         WRITE (REC,FMT=99994) FORMAT(START:START+74)
         GO TO 560
      END IF
C
      IF (BRACK) THEN
C        Insert brackets and comma into the edit descriptor.
         FORM = '(999(:,'' ('','//FORMT(START:FINISH)
     *          //','','','//FORMT(START:FINISH)//','')''))'
      ELSE
         FORM = '(999(:,'//FORMT(START:FINISH)//'))'
      END IF
C
C     Decide how wide each column of numbers is going to be,
C     by writing the number ZERO to an internal file and measuring
C     the width needed. Since the width may include trailing blanks,
C     we also write ZERO twice with the same format, and compute
C     the required field width from the two widths. Note that if
C     FORMAT has more than one edit descriptors in it, with different
C     field widths, for example FORMAT = 'E12.3,E13.4', then it may
C     not be possible to compute a sensible value, in which case
C     the columns of the output matrix will be a little skew.
      IF (ABOVE) THEN
         WRITE (INFILE,FMT=FORM,ERR=120) ZERO
         WRITE (INFIL2,FMT=FORM,ERR=120) ZERO, ZERO
      ELSE
         WRITE (INFILE,FMT=FORM,ERR=120) CZERO
         WRITE (INFIL2,FMT=FORM,ERR=120) CZERO, CZERO
      END IF
      CALL X04CBZ(INFILE,START,FINISH)
      CALL X04CBZ(INFIL2,START,FINIS2)
      NUMWID = FINIS2 - FINISH
C     NUMWID is the width of a number as printed using FORMT.
      GO TO 140
  120 CONTINUE
C     The format in FORM caused an error when used to print a number.
      IERR = 6
      NREC = 2
      WRITE (REC,FMT=99993) FORMAT(1:MIN(75,LEN(FORMAT)))
      GO TO 560
  140 CONTINUE
C
C     What kind of row labelling is required?
      IF (LABROW.EQ.'N' .OR. LABROW.EQ.'n') THEN
C        No row labelling.
         RLBWID = 1
      ELSE IF (LABROW.EQ.'I' .OR. LABROW.EQ.'i') THEN
C        Numeric row labelling.
         WRITE (INFILE,FMT='(I16)') M
         CALL X04CBZ(INFILE,START,FINISH)
         RLBWID = FINISH - START + 2
         IRFORM = '(I    )'
         WRITE (IRFORM(3:6),FMT='(I4)') RLBWID
      ELSE
C        User supplied row labelling.
         RLBWID = 1
         DO 160 I = 1, M
            CALL X04CBZ(RLABS(I),START,FINISH)
            RLBWID = MAX(RLBWID,FINISH-START+2)
  160    CONTINUE
      END IF
C
C     What kind of column labelling is required?
      IF (LABCOL.EQ.'I' .OR. LABCOL.EQ.'i') THEN
C        Numeric column labelling.
         WRITE (INFILE,FMT='(I16)') N
         CALL X04CBZ(INFILE,START,FINISH)
         CLBWID = FINISH - START + 2
         ICFORM = '(999I    )'
         WRITE (ICFORM(6:9),FMT='(I4)') NUMWID
      ELSE IF (LABCOL.NE.'N' .AND. LABCOL.NE.'n') THEN
C        User supplied column labelling.
         CLBWID = LEN(CLABS(1))
      END IF
C
      NELEMS = (INCOLS-1-RLBWID)/NUMWID
      IF (NELEMS.LT.1) THEN
         IERR = 8
         NREC = 2
         WRITE (REC,FMT=99992) INCOLS + INDNT, INDNT
         GO TO 560
      END IF
C     NELEMS is the number of elements that can fit into INCOLS columns.
C
      NSETS = (N-1)/NELEMS + 1
C     NSETS is the number of pieces that the matrix must be split into.
C
C     Print the title, splitting it up if more than INCOLS-1 characters.
      CALL X04CBZ(TITLE,START,FINISH)
      IF (FINISH.NE.0) THEN
         NTITLE = (FINISH-1)/(INCOLS-1) + 1
         DO 180 I = 1, NTITLE - 1
            TILTER = BLANKS(1:INDNT+1)//TITLE((I-1)*(INCOLS-1)
     *               +1:I*(INCOLS-1))
            CALL X04BAF(NOUT,TILTER)
  180    CONTINUE
         TILTER = BLANKS(1:INDNT+1)//TITLE((NTITLE-1)*(INCOLS-1)
     *            +1:FINISH)
         CALL X04BAF(NOUT,TILTER)
      END IF
C
C     Exit after printing the title if M or N is less than 1.
      IF (M.LT.1 .OR. N.LT.1) GO TO 560
C
C     Print the matrix, with row and column labels if requested.
      CSHIFT = 0
C     CSHIFT is the offset into the current set of columns, when
C     the matrix cannot be printed in one go but has to be split.
      DO 540 I = 1, NSETS
         IF (I.EQ.NSETS) THEN
            NELS = N - (NSETS-1)*NELEMS
         ELSE
            NELS = NELEMS
         END IF
         IF (LABCOL.EQ.'I' .OR. LABCOL.EQ.'i') THEN
C           Construct the numeric column labels.
            INFILE = ' '
            WRITE (INFILE(RLBWID+2+INDNT:),FMT=ICFORM) (CL,CL=CSHIFT+1,
     *        CSHIFT+NELS)
         ELSE IF (LABCOL.NE.'N' .AND. LABCOL.NE.'n') THEN
C           Process the user-supplied column labels.
            INFILE = ' '
            LWID = MIN(CLBWID,NUMWID)
            OFFSET = RLBWID + 1 + INDNT
            DO 200 K = CSHIFT + 1, CSHIFT + NELS
               CALL X04CBZ(CLABS(K)(1:LWID),START,FINISH)
               IF (START.EQ.0) THEN
                  START = 1
                  FINISH = 1
               END IF
               INFILE(OFFSET+NUMWID-FINISH+START:
     *           OFFSET+NUMWID+FINISH-START) = CLABS(K) (START:FINISH)
               OFFSET = OFFSET + NUMWID
  200       CONTINUE
         END IF
C        Output the column labels.
         IF (LABCOL.NE.'N' .AND. LABCOL.NE.'n') THEN
            CALL X04BAF(NOUT,INFILE(1:NCOLS))
         END IF
C
C        Now print each row in turn.
         DO 520 J = 1, M
            INFILE = ' '
C
C           Insert the row label.
            IF (LABROW.EQ.'I' .OR. LABROW.EQ.'i') THEN
               WRITE (INFILE(INDNT+1:INDNT+RLBWID),FMT=IRFORM) J
            ELSE IF (LABROW.NE.'N' .AND. LABROW.NE.'n') THEN
               CALL X04CBZ(RLABS(J),START,FINISH)
               IF (START.EQ.0) THEN
                  START = 1
                  FINISH = 1
               END IF
               INFILE(INDNT+RLBWID-FINISH+START:INDNT+RLBWID) = RLABS(J)
     *           (START:FINISH)
            END IF
C
            IF (GENERL) THEN
C              General rectangular matrix.
               IF (ABOVE) THEN
                  WRITE (INFILE(INDNT+RLBWID+2:),FMT=FORM,ERR=220)
     *              (DBLE(A(J,CL)),CL=CSHIFT+1,CSHIFT+NELS)
  220             CONTINUE
                  INFBOT(1:INDNT+RLBWID+1) = ' '
                  WRITE (INFBOT(INDNT+RLBWID+2:),FMT=FORM,ERR=240)
     *              (DIMAG(A(J,CL)),CL=CSHIFT+1,CSHIFT+NELS)
  240             CONTINUE
               ELSE
                  WRITE (INFILE(INDNT+RLBWID+2:),FMT=FORM,ERR=260)
     *              (A(J,CL),CL=CSHIFT+1,CSHIFT+NELS)
  260             CONTINUE
               END IF
            ELSE
C              Upper or lower triangular matrix.
               ND = MAX(0,J-(I-1)*NELEMS)
C              ND is the position of the Jth row diagonal element.
               IF (LOWER) THEN
                  CLEFT = CSHIFT + 1
                  CRIGHT = CSHIFT + MIN(ND-1,NELS)
               ELSE
                  CLEFT = CSHIFT + ND + 1
                  CRIGHT = CSHIFT + NELS
               END IF
C              CLEFT and CRIGHT are the leftmost and rightmost elements
C              of the current row to be printed, excluding the diagonal.
               PRDIAG = DIAG .NE. 'B' .AND. DIAG .NE. 'b' .AND. ND .GT.
     *                  0 .AND. ND .LE. NELS
C              PRDIAG is true if a diagonal element appears in the
C              current matrix row section, and it is to be printed.
               IF (PRDIAG) THEN
                  IF (DIAG.EQ.'U' .OR. DIAG.EQ.'u') THEN
                     DIAGEL = ONE
                  ELSE
                     DIAGEL = A(J,J)
                  END IF
               END IF
C
               IF (LOWER) THEN
C                 Lower triangular matrix.
                  IF (PRDIAG) THEN
                     IF (ABOVE) THEN
                        WRITE (INFILE(INDNT+RLBWID+2:),FMT=FORM,ERR=280)
     *                    (DBLE(A(J,CL)),CL=CLEFT,CRIGHT), DBLE(DIAGEL)
  280                   CONTINUE
                        INFBOT(1:INDNT+RLBWID+1) = ' '
                        WRITE (INFBOT(INDNT+RLBWID+2:),FMT=FORM,ERR=300)
     *                    (DIMAG(A(J,CL)),CL=CLEFT,CRIGHT),
     *                    DIMAG(DIAGEL)
  300                   CONTINUE
                     ELSE
                        WRITE (INFILE(INDNT+RLBWID+2:),FMT=FORM,ERR=320)
     *                    (A(J,CL),CL=CLEFT,CRIGHT), DIAGEL
  320                   CONTINUE
                     END IF
                  ELSE
                     IF (ABOVE) THEN
                        WRITE (INFILE(INDNT+RLBWID+2:),FMT=FORM,ERR=340)
     *                    (DBLE(A(J,CL)),CL=CLEFT,CRIGHT)
  340                   CONTINUE
                        INFBOT(1:INDNT+RLBWID+1) = ' '
                        WRITE (INFBOT(INDNT+RLBWID+2:),FMT=FORM,ERR=360)
     *                    (DIMAG(A(J,CL)),CL=CLEFT,CRIGHT)
  360                   CONTINUE
                     ELSE
                        WRITE (INFILE(INDNT+RLBWID+2:),FMT=FORM,ERR=380)
     *                    (A(J,CL),CL=CLEFT,CRIGHT)
  380                   CONTINUE
                     END IF
                  END IF
               ELSE
C                 Upper triangular matrix.
                  IF (PRDIAG) THEN
                     IF (ABOVE) THEN
                        WRITE (INFILE(INDNT+RLBWID+2+NUMWID*(ND-1):),
     *                    FMT=FORM,ERR=400) DBLE(DIAGEL),
     *                    (DBLE(A(J,CL)),CL=CLEFT,CRIGHT)
  400                   CONTINUE
                        INFBOT = ' '
                        WRITE (INFBOT(INDNT+RLBWID+2+NUMWID*(ND-1):),
     *                    FMT=FORM,ERR=420) DIMAG(DIAGEL),
     *                    (DIMAG(A(J,CL)),CL=CLEFT,CRIGHT)
  420                   CONTINUE
                     ELSE
                        WRITE (INFILE(INDNT+RLBWID+2+NUMWID*(ND-1):),
     *                    FMT=FORM,ERR=440) DIAGEL,
     *                    (A(J,CL),CL=CLEFT,CRIGHT)
  440                   CONTINUE
                     END IF
                  ELSE
                     IF (CLEFT.LE.CRIGHT) THEN
C                       Have to do the check on CLEFT and CRIGHT to
C                       avoid INDNT+RLBWID+2+NUMWID*ND possibly being
C                       out of range.
                        IF (ABOVE) THEN
                           WRITE (INFILE(INDNT+RLBWID+2+NUMWID*ND:),
     *                       FMT=FORM,ERR=460) (DBLE(A(J,CL)),CL=CLEFT,
     *                       CRIGHT)
  460                      CONTINUE
                           INFBOT = ' '
                           WRITE (INFBOT(INDNT+RLBWID+2+NUMWID*ND:),
     *                       FMT=FORM,ERR=480) (DIMAG(A(J,CL)),CL=CLEFT,
     *                       CRIGHT)
  480                      CONTINUE
                        ELSE
                           WRITE (INFILE(INDNT+RLBWID+2+NUMWID*ND:),
     *                       FMT=FORM,ERR=500) (A(J,CL),CL=CLEFT,CRIGHT)
  500                      CONTINUE
                        END IF
                     END IF
                  END IF
               END IF
            END IF
C
C           Output the (partial) matrix row.
            CALL X04BAF(NOUT,INFILE(1:NCOLS))
            IF (ABOVE) THEN
               CALL X04BAF(NOUT,INFBOT(1:NCOLS))
               IF (J.LT.M) CALL X04BAF(NOUT,' ')
            END IF
C
  520    CONTINUE
C
         CSHIFT = CSHIFT + NELEMS
         IF (I.NE.NSETS) THEN
            CALL X04BAF(NOUT,' ')
         END IF
  540 CONTINUE
C
  560 IFAIL = P01ABF(IFAIL,IERR,SRNAME,NREC,REC)
C
      RETURN
C
99999 FORMAT (1X,'** On entry, MATRIX is not valid : MATRIX = ''',A,
     *       '''')
99998 FORMAT (1X,'** On entry, M.gt.LDA : M = ',I16,', LDA = ',I16)
99997 FORMAT (1X,'** On entry, MATRIX = ''',A,''', but DIAG is not val',
     *       'id : DIAG = ''',A,'''')
99996 FORMAT (1X,'** On entry, either LABROW or LABCOL is not valid',
     *       /4X,'LABROW = ''',A,''', LABCOL = ''',A,'''.')
99995 FORMAT (1X,'** On entry, USEFRM is not valid : USEFRM = ''',A,
     *       '''         .')
99994 FORMAT (1X,'** On entry, FORMAT has more than 80 characters: the',
     *       ' first 75 are',/4X,A)
99993 FORMAT (1X,'** The format specifier in FORMAT cannot be used to ',
     *       'print a number: FORMAT =',/4X,A)
99992 FORMAT (1X,'** On entry, NCOLS-INDENT is not wide enough to hold',
     *       ' at least one matrix',/4X,'column: NCOLS = ',I16,', INDE',
     *       'NT =',I16)
      END
*
      SUBROUTINE P01ABW(N,NAME,INFORM,IERR,SRNAME)
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C
C     P01ABW increases the value of IERR by 1 and, if
C
C        ( mod( INFORM, 10 ).ne.1 ).or.( mod( INFORM/10, 10 ).ne.0 )
C
C     writes a message on the current error message channel giving the
C     value of N, a message to say that N is invalid and the strings
C     NAME and SRNAME.
C
C     NAME must be the name of the actual argument for N and SRNAME must
C     be the name of the calling routine.
C
C     This routine is intended for use when N is an invalid input
C     parameter to routine SRNAME. For example
C
C        IERR = 0
C        IF( N.NE.'Valid value' )
C     $     CALL P01ABW( N, 'N', IDIAG, IERR, SRNAME )
C
C  -- Written on 15-November-1984.
C     Sven Hammarling, Nag Central Office.
C
C     .. Scalar Arguments ..
      INTEGER           IERR, INFORM
      CHARACTER*(*)     N
      CHARACTER*(*)     NAME, SRNAME
C     .. Local Scalars ..
      INTEGER           NERR
C     .. Local Arrays ..
      CHARACTER*65      REC(3)
C     .. External Subroutines ..
      EXTERNAL          X04AAF, X04BAF
C     .. Intrinsic Functions ..
      INTRINSIC         MOD
C     .. Executable Statements ..
      IERR = IERR + 1
      IF ((MOD(INFORM,10).NE.1) .OR. (MOD(INFORM/10,10).NE.0)) THEN
         CALL X04AAF(0,NERR)
         WRITE (REC,FMT=99999) NAME, SRNAME, N
         CALL X04BAF(NERR,' ')
         CALL X04BAF(NERR,REC(1))
         CALL X04BAF(NERR,REC(2))
         CALL X04BAF(NERR,REC(3))
      END IF
      RETURN
C
C
C     End of P01ABW.
C
99999 FORMAT (' *****  Parameter  ',A,'  is invalid in routine  ',A,
     *  '  ***** ',/8X,'Value supplied is',/8X,A)
      END
*
      SUBROUTINE P01ABY(N,NAME,INFORM,IERR,SRNAME)
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C
C     P01ABY increases the value of IERR by 1 and, if
C
C        ( mod( INFORM, 10 ).ne.1 ).or.( mod( INFORM/10, 10 ).ne.0 )
C
C     writes a message on the current error message channel giving the
C     value of N, a message to say that N is invalid and the strings
C     NAME and SRNAME.
C
C     NAME must be the name of the actual argument for N and SRNAME must
C     be the name of the calling routine.
C
C     This routine is intended for use when N is an invalid input
C     parameter to routine SRNAME. For example
C
C        IERR = 0
C        IF( N.LT.1 )CALL P01ABY( N, 'N', IDIAG, IERR, SRNAME )
C
C  -- Written on 23-February-1984.  Sven.
C
C     .. Scalar Arguments ..
      INTEGER           IERR, INFORM, N
      CHARACTER*(*)     NAME, SRNAME
C     .. Local Scalars ..
      INTEGER           NERR
C     .. Local Arrays ..
      CHARACTER*65      REC(2)
C     .. External Subroutines ..
      EXTERNAL          X04AAF, X04BAF
C     .. Intrinsic Functions ..
      INTRINSIC         MOD
C     .. Executable Statements ..
      IERR = IERR + 1
      IF ((MOD(INFORM,10).NE.1) .OR. (MOD(INFORM/10,10).NE.0)) THEN
         CALL X04AAF(0,NERR)
         WRITE (REC,FMT=99999) NAME, SRNAME, N
         CALL X04BAF(NERR,' ')
         CALL X04BAF(NERR,REC(1))
         CALL X04BAF(NERR,REC(2))
      END IF
      RETURN
C
C
C     End of P01ABY.
C
99999 FORMAT (' *****  Parameter  ',A,'  is invalid in routine  ',A,
     *  '  ***** ',/8X,'Value supplied is ',I6)
      END
*
      SUBROUTINE F02EAZ(AMAX,RMIN,RMAX,SIGMA,SCALE)
C     MARK 16 RELEASE. NAG COPYRIGHT 1992.
C
C     F02EAZ determines a scaling factor SIGMA such that if AMAX
C     lies outside the range RMIN to RMAX, then AMAX*SIGMA lies within
C     this range (except that SIGMA must not overflow or underflow).
C     SIGMA is constrained to be a power of the base of
C     floating-point arithmetic.
C
C     SCALE is set to .TRUE. if scaling is required.
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE
      PARAMETER         (ZERO=0.0D+0,ONE=1.0D+0)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  AMAX, RMAX, RMIN, SIGMA
      LOGICAL           SCALE
C     .. Local Scalars ..
      DOUBLE PRECISION  BASE, FAC, SIGLOG
      INTEGER           IPSIG
C     .. External Functions ..
      INTEGER           X02BHF, X02BKF, X02BLF
      EXTERNAL          X02BHF, X02BKF, X02BLF
C     .. Intrinsic Functions ..
      INTRINSIC         INT, LOG, MIN, MOD
C     .. Executable Statements ..
C
      BASE = X02BHF()
      IF (AMAX.GT.RMAX) THEN
C
C        SIGMA should be the largest power of the base <= RMAX/AMAX
C
         SCALE = .TRUE.
         SIGLOG = LOG(AMAX) - LOG(RMAX)
         FAC = ONE/BASE
      ELSE IF (AMAX.LT.RMIN .AND. AMAX.GT.ZERO) THEN
C
C        SIGMA should be the smallest power of the base > RMIN/AMAX
C
         SCALE = .TRUE.
         SIGLOG = LOG(RMIN) - LOG(AMAX)
         FAC = BASE
      ELSE
         SCALE = .FALSE.
      END IF
      SIGMA = ONE
      IF (SCALE) THEN
         IPSIG = MIN(INT(SIGLOG/LOG(BASE))+1,-X02BKF(),X02BLF()-1)
C
C        SIGMA = FAC**IPSIG
C
         GO TO 40
   20    CONTINUE
         FAC = FAC*FAC
   40    IF (MOD(IPSIG,2).GT.0) SIGMA = SIGMA*FAC
         IPSIG = IPSIG/2
         IF (IPSIG.GT.0) GO TO 20
      END IF
      RETURN
      END
*
      INTEGER FUNCTION P01ABF(IFAIL,IERROR,SRNAME,NREC,REC)
C     MARK 11.5(F77) RELEASE. NAG COPYRIGHT 1986.
C     MARK 13 REVISED. IER-621 (APR 1988).
C     MARK 13B REVISED. IER-668 (AUG 1988).
C
C     P01ABF is the error-handling routine for the NAG Library.
C
C     P01ABF either returns the value of IERROR through the routine
C     name (soft failure), or terminates execution of the program
C     (hard failure). Diagnostic messages may be output.
C
C     If IERROR = 0 (successful exit from the calling routine),
C     the value 0 is returned through the routine name, and no
C     message is output
C
C     If IERROR is non-zero (abnormal exit from the calling routine),
C     the action taken depends on the value of IFAIL.
C
C     IFAIL =  1: soft failure, silent exit (i.e. no messages are
C                 output)
C     IFAIL = -1: soft failure, noisy exit (i.e. messages are output)
C     IFAIL =-13: soft failure, noisy exit but standard messages from
C                 P01ABF are suppressed
C     IFAIL =  0: hard failure, noisy exit
C
C     For compatibility with certain routines included before Mark 12
C     P01ABF also allows an alternative specification of IFAIL in which
C     it is regarded as a decimal integer with least significant digits
C     cba. Then
C
C     a = 0: hard failure  a = 1: soft failure
C     b = 0: silent exit   b = 1: noisy exit
C
C     except that hard failure now always implies a noisy exit.
C
C     S.Hammarling, M.P.Hooper and J.J.du Croz, NAG Central Office.
C
C     .. Scalar Arguments ..
      INTEGER                 IERROR, IFAIL, NREC
      CHARACTER*(*)           SRNAME
C     .. Array Arguments ..
      CHARACTER*(*)           REC(*)
C     .. Local Scalars ..
      INTEGER                 I, NERR
      CHARACTER*72            MESS
C     .. External Subroutines ..
      EXTERNAL                P01ABZ, X04AAF, X04BAF
C     .. Intrinsic Functions ..
      INTRINSIC               ABS, MOD
C     .. Executable Statements ..
      IF (IERROR.NE.0) THEN
C        Abnormal exit from calling routine
         IF (IFAIL.EQ.-1 .OR. IFAIL.EQ.0 .OR. IFAIL.EQ.-13 .OR.
     *       (IFAIL.GT.0 .AND. MOD(IFAIL/10,10).NE.0)) THEN
C           Noisy exit
            CALL X04AAF(0,NERR)
            DO 20 I = 1, NREC
               CALL X04BAF(NERR,REC(I))
   20       CONTINUE
            IF (IFAIL.NE.-13) THEN
               WRITE (MESS,FMT=99999) SRNAME, IERROR
               CALL X04BAF(NERR,MESS)
               IF (ABS(MOD(IFAIL,10)).NE.1) THEN
C                 Hard failure
                  CALL X04BAF(NERR,
     *                     ' ** NAG hard failure - execution terminated'
     *                        )
                  CALL P01ABZ
               ELSE
C                 Soft failure
                  CALL X04BAF(NERR,
     *                        ' ** NAG soft failure - control returned')
               END IF
            END IF
         END IF
      END IF
      P01ABF = IERROR
      RETURN
C
99999 FORMAT (' ** ABNORMAL EXIT from NAG Library routine ',A,': IFAIL',
     *  ' =',I6)
      END
*
      INTEGER FUNCTION P01ACF(IFAIL,IERROR,SRNAME,VARBNM,NREC,REC)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C
C     P01ACF is the error-handling routine for the F06 AND F07
C     Chapters of the NAG Fortran Library. It is a slightly modified
C     version of P01ABF.
C
C     P01ACF either returns the value of IERROR through the routine
C     name (soft failure), or terminates execution of the program
C     (hard failure). Diagnostic messages may be output.
C
C     If IERROR = 0 (successful exit from the calling routine),
C     the value 0 is returned through the routine name, and no
C     message is output
C
C     If IERROR is non-zero (abnormal exit from the calling routine),
C     the action taken depends on the value of IFAIL.
C
C     IFAIL =  1: soft failure, silent exit (i.e. no messages are
C                 output)
C     IFAIL = -1: soft failure, noisy exit (i.e. messages are output)
C     IFAIL =-13: soft failure, noisy exit but standard messages from
C                 P01ACF are suppressed
C     IFAIL =  0: hard failure, noisy exit
C
C     For compatibility with certain routines included before Mark 12
C     P01ACF also allows an alternative specification of IFAIL in which
C     it is regarded as a decimal integer with least significant digits
C     cba. Then
C
C     a = 0: hard failure  a = 1: soft failure
C     b = 0: silent exit   b = 1: noisy exit
C
C     except that hard failure now always implies a noisy exit.
C
C     S.Hammarling, M.P.Hooper and J.J.du Croz, NAG Central Office.
C
C     .. Scalar Arguments ..
      INTEGER                 IERROR, IFAIL, NREC
      CHARACTER*(*)           SRNAME, VARBNM
C     .. Array Arguments ..
      CHARACTER*(*)           REC(*)
C     .. Local Scalars ..
      INTEGER                 I, NERR, VARLEN
      CHARACTER*72            MESS
C     .. External Subroutines ..
      EXTERNAL                P01ABZ, X04AAF, X04BAF
C     .. Intrinsic Functions ..
      INTRINSIC               ABS, LEN, MOD
C     .. Executable Statements ..
      IF (IERROR.NE.0) THEN
         VARLEN = 0
         DO 20 I = LEN(VARBNM), 1, -1
            IF (VARBNM(I:I).NE.' ') THEN
               VARLEN = I
               GO TO 40
            END IF
   20    CONTINUE
   40    CONTINUE
C        Abnormal exit from calling routine
         IF (IFAIL.EQ.-1 .OR. IFAIL.EQ.0 .OR. IFAIL.EQ.-13 .OR.
     *       (IFAIL.GT.0 .AND. MOD(IFAIL/10,10).NE.0)) THEN
C           Noisy exit
            CALL X04AAF(0,NERR)
            DO 60 I = 1, NREC
               CALL X04BAF(NERR,REC(I))
   60       CONTINUE
            IF (IFAIL.NE.-13) THEN
               IF (VARLEN.NE.0) THEN
                  WRITE (MESS,FMT=99999) SRNAME, VARBNM(1:VARLEN),
     *              IERROR
               ELSE
                  WRITE (MESS,FMT=99998) SRNAME
               END IF
               CALL X04BAF(NERR,MESS)
               IF (ABS(MOD(IFAIL,10)).NE.1) THEN
C                 Hard failure
                  CALL X04BAF(NERR,
     *                     ' ** NAG hard failure - execution terminated'
     *                        )
                  CALL P01ABZ
               ELSE
C                 Soft failure
                  CALL X04BAF(NERR,
     *                        ' ** NAG soft failure - control returned')
               END IF
            END IF
         END IF
      END IF
      P01ACF = IERROR
      RETURN
C
99999 FORMAT (' ** ABNORMAL EXIT from NAG Library routine ',A,': ',A,
     *       ' =',I6)
99998 FORMAT (' ** ABNORMAL EXIT from NAG Library routine ',A)
      END
