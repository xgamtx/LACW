        SUBROUTINE TTIME(TOTAL_TIME,TTEXT,IOUT)
c ---- This subprogram calculates time of calculations, for PC's only! ----
c      Here TOTAL_TIME - REAL*16 variable for storing current data/time
c           TTEXT - comment,
c           IOUT - if <> 0, output channel for LOG-file.
*
        INTEGER*4 IDAYS,I,IOUT
	INTEGER*2 IYR,IMON,IDAY,IHR,IMIN,ISEC,I100th,IMNT(12)
	DOUBLE PRECISION CUR_TIME ,TOTAL_TIME
	CHARACTER*(*) TTEXT
*
	DATA IMNT/30,28,31,30,31,30,31,31,30,31,30,31/
	CALL GETDAT(IYR,IMON,IDAY)
	CALL GETTIM(Ihr,Imin,Isec,I100th)
	IDAYS=365*(IYR-1981)+(IYR-1981)/4+1
	IF(IMON.GT.2.AND.MOD(IYR,4).EQ.0) IDAYS=IDAYS+1
	IF(IMON.GT.1) THEN
        DO 1, J=1, IMON-1
   1	   IDAYS=IDAYS+IMNT(J)
	ENDIF
	IDAYS=IDAYS+IDAY-1
	CUR_TIME=DBLE(24*3600*IDAYS+Ihr*3600+Imin*60+Isec)+
     .              dble(I100th)/100D0
	IF(TOTAL_TIME.NE.0D0) THEN
	   TOTAL_TIME=CUR_TIME-TOTAL_TIME
 	   I=LEN(TTEXT)
           WRITE(*,2) TTEXT
           IF(IOUT.GT.0) WRITE(IOUT,2) TTEXT
	   IHR=IDNINT(TOTAL_TIME)/3600
           IF(IHR.NE.0) THEN
              WRITE(*,3) IHR
              IF(IOUT.GT.0) WRITE(IOUT,3) IHR
           ENDIF
	   IMIN=(IDNINT(TOTAL_TIME)-IHR*3600)/60
           IF(IMIN.NE.0) THEN
             WRITE(*,4) IMIN
             IF(IOUT.GT.0) WRITE(IOUT,4) IMIN
           ENDIF
	   TOTAL_TIME=TOTAL_TIME-DBLE(3600*IHR+60*IMIN)
           WRITE(*,5) DBLE(IDNINT(TOTAL_TIME*100.0D0))/100.0D0
           IF(IOUT.GT.0)
     .        WRITE(IOUT,5) DBLE(IDNINT(TOTAL_TIME*100.0D0))/100.0D0
	ELSE
	   TOTAL_TIME=CUR_TIME
	ENDIF
   2	FORMAT(1X,A,1X,\)
   3	FORMAT(I4,' hour(s)',\)
   4	FORMAT(1X,I2,' min.',\)
   5	FORMAT(1X,F5.2,' sec.')
      RETURN
	END
