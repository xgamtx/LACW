      SUBROUTINE ExtrCalcData(DAT,MCnt,NCnt,Msz,Nsz,C,Top,Bot,OCHNL,
     +                        Mode,IfZERO,IWRK)
c   -- This routine extracts calculated DATA from file BND.STR,
c   -- contained energy values for bands.
*
* See Developer's Remark near the call to ExtrCalcData().
* Briefly: It happened so that NCnt passed here was NTPNTS
* and it was equal to the number of *intervals* onto which the k axis
* is being splitted. But the DAT (Ebnd) array contains NTPNTS+1
* points, 0..NTPNTS. Since 3.11, NCnt equals to the number of points,
* i.e., the number of diagonalisations done.
*
      INTEGER MCnt,NCnt,Msz,Nsz,OCHNL,IWRK(*),I,J,K,II,Dir
      DOUBLE PRECISION DAT(Msz,Nsz),C,Top,Bot,PI,Kcur,Kprv,TMP,TMP1
      LOGICAL Mode,IfZERO,LPlot,CrossTop,CrossBot,IfOpened
      CHARACTER*1 TNMBR(10), outf*11, TDAT*4, TBND*5
      DATA Dir/1/, TDAT/'.dat'/, TBND/'bnd_0'/,K/0/,
     .     TNMBR/'0','1','2','3','4','5','6','7','8','9'/,
     .     CrossTop,CrossBot/.FALSE.,.FALSE./

	print *,'EBND: ',Msz,' eigenvalues at ',Nsz,' points'

      PI=4.D0*DATAN(1.d0)
      J=0
      DO 55, I=1, 99    ! Cleaning of dirty files from last sessions, if any.
        outf=TBND//TNMBR(I/10+1)//TNMBR(I-10*(I/10)+1)//TDAT
        OPEN(OCHNL,FILE=outf,FORM='FORMATTED',STATUS='OLD',ERR=55)
        J=J+1
        IF(J.EQ.1)
     .     WRITE(*,'(1X,A35,\)') 'Cleaning rests of old plot files...'
        CLOSE(OCHNL, STATUS='DELETE')
   55 CONTINUE
      IF(J.GT.0) THEN
         WRITE(*,35)
         WRITE(*,'(I2,A18/)') J,' file(s) deleted.'
      ENDIF
      DO I=1,MCnt               ! Finding of "solid" bands
        II=0
        DO J=1,NCnt
          IF(DAT(I,J).LE.Top .AND. DAT(I,J).GE.Bot) THEN
            IF(J.EQ.1) II=2
            IF(II.NE.2) II=1
          ELSE
            IF(II.EQ.2) II=1
          ENDIF
        END DO
        IF(Mode) THEN
          K=K+1
         ELSE
          IF(II.EQ.2) K=K+1
        ENDIF
        IWRK(I)=II
      END DO
      WRITE(*,'(/1X,A38,\)') 'Building of plot file of the bands...'
      I=1
      J=1
      OPEN(OCHNL,FILE='bnd.dat',FORM='FORMATTED')
      IF(Mode) THEN
        DO WHILE(I.LE.MCnt)
          DO WHILE(IWRK(I).EQ.0 .AND. I.LT.MCnt)
            I=I+1
          END DO
          IF(IfZERO) THEN
            Kcur=DBLE(J-1)*PI/DBLE(NCnt-1)/C
           ELSE
            Kcur=-PI/C+2.D0*DBLE(J-1)*PI/DBLE(NCnt-1)/C
          ENDIF
          IF(J.GT.1 .AND. J.LT.NCnt) THEN
             IF(IfZERO) THEN
               Kprv=DBLE(J-1-Dir)*PI/DBLE(NCnt-1)/C
              ELSE
               Kprv=-PI/C+2.D0*DBLE(J-1-Dir)*PI/DBLE(NCnt-1)/C
             ENDIF
             IF((DAT(I,J)-Top)*(DAT(I,J-Dir)-Top) .LT. 0.d0) THEN
               TMP1=Kprv+(Kprv-Kcur)/(DAT(I,J-Dir)-DAT(I,J))
     .              *(Top-DAT(I,J-Dir))
               IF(.NOT.CrossTop) THEN
                 IF(.NOT.LPlot) THEN
                   TMP=-PI/C
                   IF(IfZERO) TMP=0.D0
                   IF(Dir.EQ.1) THEN
                      WRITE(OCHNL,17) TMP, Top
                      WRITE(OCHNL,17) TMP1, Top
                    ELSE
                      WRITE(OCHNL,17) PI/C, Top
                      WRITE(OCHNL,17) TMP1, Top
                   ENDIF
                  ELSE
                   WRITE(OCHNL,17) TMP1, Top
                 ENDIF
                ELSE
                 WRITE(OCHNL,17) TMP1, Top
               ENDIF
               CrossTop=.TRUE.
             ENDIF
             IF((DAT(I,J)-Bot)*(DAT(I,J-Dir)-Bot) .LT. 0.d0) THEN
               TMP1=Kprv+(Kprv-Kcur)/(DAT(I,J-Dir)-DAT(I,J))
     .              *(Bot-DAT(I,J-Dir))
               IF(.NOT.CrossBot) THEN
                 IF(.NOT.LPlot) THEN
                   TMP=-PI/C
                   IF(IfZERO) TMP=0.D0
                   IF(Dir.EQ.1) THEN
                      WRITE(OCHNL,17) TMP, Bot
                      WRITE(OCHNL,17) TMP1, Bot
                    ELSE
                      WRITE(OCHNL,17) TMP1, Bot
                   ENDIF
                  ELSE
                   WRITE(OCHNL,17) TMP1, Bot
                 ENDIF
                ELSE
                 WRITE(OCHNL,17) TMP1, Bot
               ENDIF
               CrossBot=.TRUE.
             ENDIF
          ENDIF
          IF(DAT(I,J).LE.Top .AND. DAT(I,J).GE.Bot) THEN
             WRITE(OCHNL,17) Kcur, DAT(I,J)
             LPlot=.TRUE.       ! LPlot: status of current point,
           ELSE                 ! inside/outside of interval [Bot;Top].
             LPlot=.FALSE.
          ENDIF
          IF(Dir.EQ.1 .AND. J.EQ.Ncnt) THEN
             IF(.NOT.LPlot) THEN
               IF(CrossTop) WRITE(OCHNL,17) PI/C, Top
               IF(CrossBot) WRITE(OCHNL,17) PI/C, Bot
             ENDIF
             Dir=-1
             I=I+1
             LPlot=.FALSE.
             CrossTop=.FALSE.
             CrossBot=.FALSE.
           ELSEIF(Dir.EQ.-1 .AND. J.EQ.1) THEN
             IF(.NOT.LPlot) THEN
               TMP=-PI/C
               IF(IfZERO) TMP=0.D0
               IF(CrossTop) WRITE(OCHNL,17) TMP, Top
               IF(CrossBot) WRITE(OCHNL,17) TMP, Bot
             ENDIF
             Dir=1
             I=I+1
             LPlot=.FALSE.
             CrossTop=.FALSE.
             CrossBot=.FALSE.
           ELSE
             J=J+Dir
          ENDIF
        END DO
        CLOSE(OCHNL)
        WRITE(*,35)
      ELSE
        DO WHILE(I.LE.MCnt)
          DO WHILE(IWRK(I).NE.2 .AND. I.LT.MCnt)        ! plot only "solid" bands.
            I=I+1
          END DO
          IF(I.GT.MCnt) EXIT
          IF(IfZERO) THEN
            Kcur = DBLE(J-1)*PI/DBLE(NCnt-1)/C
           ELSE
            Kcur = -PI/C+2.D0*DBLE(J-1)*PI/DBLE(NCnt-1)/C
          ENDIF
          IF(DAT(I,J).LE.Top .AND. DAT(I,J).GE.Bot)
     .               WRITE(OCHNL,17) Kcur, DAT(I,J)
          IF(Dir.EQ.1 .AND. J.EQ.Ncnt) THEN
             Dir=-1
             I=I+1
           ELSEIF(Dir.EQ.-1 .AND. J.EQ.1) THEN
             Dir=1
             I=I+1
           ELSE
             J=J+Dir
          ENDIF
        END DO
        CLOSE(OCHNL)
        II=0
        IfOpened=.FALSE.
        I=1
        DO WHILE(I.LE.MCnt)
          LPlot=.FALSE.
          DO WHILE(IWRK(I).NE.1 .AND. I.LT.MCnt)        ! Proceeding with "cutted" bands.
            I=I+1
          END DO
          IF(I.GT.MCnt) EXIT
          DO J=1,NCnt
            IF(IfZERO) THEN
              Kcur=DBLE(J-1)*PI/DBLE(NCnt-1)/C
             ELSE
              Kcur=-PI/C+2.D0*DBLE(J-1)*PI/DBLE(NCnt-1)/C
            ENDIF
            IF(J.GT.1) THEN
              IF(IfZERO) THEN
                 Kprv=DBLE(J-2)*PI/DBLE(NCnt-1)/C
               ELSE
                 Kprv=-PI/C+2.D0*DBLE(J-2)*PI/DBLE(NCnt-1)/C
              ENDIF
              IF((DAT(I,J)-Top)*(DAT(I,J-1)-Top).LT.0.d0) THEN
                TMP=Kprv+(Kprv-Kcur)/(DAT(I,J-1)-DAT(I,J))
     .              *(Top-DAT(I,J-1))
                IF(IfOpened) THEN
                  WRITE(OCHNL,17) TMP, Top
                  CLOSE(OCHNL)
                  IfOpened=.FALSE.
                 ELSE
                  CALL MkFile(OCHNL)    ! For assigning a new file
                  IfOpened=.TRUE.
                  WRITE(OCHNL,17) TMP, Top
                ENDIF
              ENDIF
              IF((DAT(I,J)-Bot)*(DAT(I,J-1)-Bot).LT.0.d0) THEN
                TMP=Kprv+(Kprv-Kcur)/(DAT(I,J-1)-DAT(I,J))
     .              *(Bot-DAT(I,J-1))
                IF(IfOpened) THEN
                  WRITE(OCHNL,17) TMP, Bot
                  CLOSE(OCHNL)
                  IfOpened=.FALSE.
                 ELSE
                  CALL MkFile(OCHNL)    ! For assigning a new file
                  WRITE(OCHNL,17) TMP, Bot
                  IfOpened=.TRUE.
                ENDIF
              ENDIF
            ENDIF
            IF(DAT(I,J).LE.Top .AND. DAT(I,J).GE.Bot) THEN
              IF(.NOT.IfOpened) THEN
                CALL MkFile(OCHNL)      ! For assigning a new file
                IfOpened=.TRUE.
              ENDIF
              LPlot=.TRUE.
             ELSE
              CLOSE(OCHNL)
              LPlot=.FALSE.
              IfOpened=.FALSE.
            ENDIF
            IF(LPlot) THEN
              WRITE(OCHNL,17) Kcur, DAT(I,J)
              IF(J.EQ.NCnt) THEN
                CLOSE(OCHNL)
                IfOpened=.FALSE.
              ENDIF
            ENDIF
          END DO
          I=I+1
        END DO
        WRITE(*,35)
        II=-1
        CALL MkFile(II)         ! How many files we are created?
        IF(II.GT.0) WRITE(*,25) II
      ENDIF
      IF(K.GT.0) WRITE(*,20) K
   17 FORMAT(F5.3,1X,F9.3)
   20 FORMAT(I3,' band(s) saved to main bands file, BND.DAT.')
   25 FORMAT(I3,' file(s), contained cutted bands, formed '\,
     +       ' additionally.'/)
   35 FORMAT('done.'/)
      RETURN
      END
*
      SUBROUTINE MkFile(INDX)
      INTEGER INDX,NUMBER
c   -- INDX: if -1, then write out number of created files. --
      CHARACTER*1 TNMBR(10), outfile*11, TDAT*4, TBND*5
      SAVE NUMBER
      DATA NUMBER/0/, TDAT/'.dat'/, TBND/'bnd_0'/,
     .     TNMBR/'0','1','2','3','4','5','6','7','8','9'/
*
      IF(INDX.EQ.-1) THEN
        INDX=NUMBER
        RETURN
      ENDIF
      NUMBER=NUMBER+1
      IF(NUMBER.GT.99) THEN
         WRITE(*,1)
         STOP ' *ERR* Building of additional files of bands.'
      ENDIF
      outfile=TBND//TNMBR(NUMBER/10+1)//TNMBR(NUMBER-10*(NUMBER/10)+1)
     +        //TDAT
      OPEN(INDX, FILE=outfile, FORM='FORMATTED')
      RETURN
    1 FORMAT(/1X,'*ERR* I''m sorry, can''t create greather than 99',\
     .       'files!')
      END
