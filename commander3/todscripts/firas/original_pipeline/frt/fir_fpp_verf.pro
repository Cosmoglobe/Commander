	PRO  FIR_FPP_VERF, JSTART, JSTOP, PLTDEV, CHANNEL, REFSET
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
;+
;	PRO  FIR_FPP_VERF, JSTART, JSTOP, PLTDEV, CHANNEL, REFSET
;
;	AUTHOR:	Larry P. Rosen, STX, 15 April 1991
;
;	PURPOSE:  This facility makes a "timeline" plot to help verify the
;		results of running FPP firas facility.  The plots compare the
;		values FPP gets for gain or fakeit status with those values in
;		the reference data sets.  Also plotted is the "Bad Time Flag"
;		which indicates the flag set in FPP.
;	PARAMETERS:
;		JSTART: Start of data for plots.
;		JSTOP:  Stop time of data for plots.
;		PLTDEV: Plot device, ie. 'REGIS', 'PS', etc.
;		CHANNEL: (optional) Channel to plot, ie. 'RH'. Default = ALL
;		REFSET: (optional) 'GAIN' or 'FAKE' data to plot. Default = both
;	INPUT:
;		CSDR$FIRAS_REF:FEX_GAIN_R.DAT;
;		CSDR$FIRAS_REF:FEX_FAKEIT_R.DAT;
;		CSDR$FIRAS_REF:FEX_GAIN_L.DAT;
;		CSDR$FIRAS_REF:FEX_FAKEIT_L.DAT;
;		CSDR$FIRAS_IN:FPP_SDF_xx.ED_**;
;
;	OUTPUT:  Plot of data for time range.  If PLTDEV is not 'REGIS' then
;		a plot file is created called FIR_FPP.XXX, where XXX indicates
;		type of plot, ie. PS for postscript.
;
;	Modifications:  July 1991, L. Rosen, STX
;		FEX_GAIN and FEX_FAKEIT were split into left and right sides.
;		This program was modified to accomodate the change.
;		Also, added channel and reference data qualifiers so that
;		user can plot one channel or all, and user can pick gain, fakeit
;		or both.
;		R:  0 = Gain,  1 = Fakeit,  2 = both
;		C:  0 = RH, 1 = RL, 2 = LH, 3 = LL, 4 = ALL
;
;               Oct 4, 1991 - K.Jensen (STX) - Changed names of IDL structure
;               functions to be compatible with Pass2-Phase1 L1 structure
;               procedure names.
;
;               Oct 4, 1991 - K.Jensen (STX) - Corrected the "PLOT_TO" statement
;               for the PLTDEV = 'TEK' option.
;
;               Oct 4, 1991 - K.Jensen (STX) - Added statements to properly
;               close the Plot File.
;-
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
;
print,'FIR_FPP_VERF, JSTART, JSTOP, PLTDEV ( ,CHANNEL, REFSET )'
;
; Parameters
DATATYPE = 'FPP_SDF_'
REFSETGR = 'CSDR$FIRAS_REF:FEX_GAIN_R.DAT;'
REFSETFR = 'CSDR$FIRAS_REF:FEX_FAKEIT_R.DAT;'
REFSETGL = 'CSDR$FIRAS_REF:FEX_GAIN_L.DAT;'
REFSETFL = 'CSDR$FIRAS_REF:FEX_FAKEIT_L.DAT;'
ARC = 'CSDR$FIRAS_IN:'
FEXT = '.ED_'
C = SIZE (CHANNEL)
CHANS = ['RH','RL','LH','LL']
IF (C(1) EQ 0) THEN  CHAN = CHANS ELSE BEGIN
   CHAN = STRUPCASE (CHANNEL)
   IF (CHAN EQ 'ALL') THEN CHAN = CHANS
ENDELSE
NCHAN = N_ELEMENTS (CHAN)
IF (NCHAN EQ 4) THEN C = 4 ELSE C = WHERE (CHAN EQ CHANS)
R = SIZE (REFSET)
SETS = ['G','F']
IF (R(1) EQ 0) THEN SET = SETS  ELSE BEGIN
   SET = STRUPCASE (STRMID (REFSET,0,1) )
   IF (SET EQ 'B') THEN SET = SETS
ENDELSE
NSETS = N_ELEMENTS (SET)
IF (NSETS EQ 2) THEN R = 2 ELSE R = WHERE (SET EQ SETS)
;-------------------------------------------------------------------------------
PLTDEV = STRUPCASE(PLTDEV)
SET_PLOT, PLTDEV
IF (PLTDEV EQ 'PS') THEN DEVICE, FILENAME='FIR_FPP.PS'
IF (PLTDEV EQ 'TEK') THEN BEGIN
   OPENW,LUNP, 'FIR_FPP.TEK',/GET_LUN
   DEVICE,PLOT_TO = LUNP 
ENDIF
TSTART = '89001000000000'   &   TSTOP = '90365235959999'
STRPUT, TSTART, JSTART   &   STRPUT, TSTOP, JSTOP
DAY = STRMID (TSTART,0,5)
TREF1 = CONV_GMT_SEC (BYTE(TSTART))
TREF2 = CONV_GMT_SEC (BYTE(TSTOP))
OFFSET = 365*(LONG(STRMID(TSTART,0,2))-89) + LONG(STRMID(TSTART,2,3)) - 324
;
;------------------------Getting Reference Gain Data----------------------------
;
IF (R(0) EQ 0 OR R(0) EQ 2) THEN BEGIN
   GOFFSET = 10 * OFFSET
;--------------------------------Right Channels---------------------------------
   IF (C(0) EQ 0 OR C(0) EQ 1 OR C(0) EQ 4) THEN BEGIN
      PRINT,' GETTING GAIN RIGHT REFERENCE DATA'
      OPENR, LUN, REFSETGR, /GET_LUN, /SHARE
      FINFO = FSTAT (LUN)
      NUMREC = FINFO.SIZE / FINFO.REC_LEN
      DATAS = ASSOC ( LUN, FEX_GAIN_ST (NUMREC-GOFFSET), GOFFSET )
      DATA = DATAS(0)
      TIM = CONV_ADT_SEC (DATA.CT_HEAD.TIME)
      WHILE (TIM(0) GT TREF1) DO BEGIN
         PRINT, ' GAIN OFFSET OF ',GOFFSET,' IS TOO BIG'
         GOFFSET = GOFFSET - OFFSET
         DATAS = ASSOC ( LUN, FEX_GAIN_ST(NUMREC-GOFFSET), GOFFSET )
         DATA = DATAS(0)
         TIM = CONV_ADT_SEC (DATA.CT_HEAD.TIME)
      ENDWHILE
      TRANGE = WHERE (TIM GE TREF1 AND TIM LE TREF2, COUNT)
      IF (COUNT GT 0) THEN BEGIN
         TRANGE = [TRANGE(0)-1, TRANGE]
         NUM = COUNT + 1
      ENDIF ELSE BEGIN
         TRANGE = WHERE (TIM LT TREF1, COUNT)
         T = TRANGE (COUNT-1)
         TRANGE = [T,T+1]
         NUM = 2
      ENDELSE
      PRINT,' NUMBER OF GAIN RIGHT RECORDS IS ',NUM
      EVEN = 2*INDGEN (NUM)
      ODD = EVEN + 1
      NGR = NUM + NUM
      G_TIMER = DBLARR(NGR,/NOZERO)
      G_TIMER(EVEN) = CONV_ADT_SEC (DATA(TRANGE).START_TIME)
      G_TIMER(ODD)  = CONV_ADT_SEC (DATA(TRANGE).STOP_TIME)
      G_GAINR = INTARR(2,NGR,/NOZERO)
      G_GAINR(*,EVEN) = DATA(TRANGE).GAIN
      G_GAINR(*,ODD) = G_GAINR(*,EVEN)
      DATA = 0 & DATAS = 0
      CLOSE, LUN  &  FREE_LUN, LUN
   ENDIF
;--------------------------------Left Channels----------------------------------
   IF (C(0) EQ 2 OR C(0) EQ 3 OR C(0) EQ 4) THEN BEGIN
      PRINT,' GETTING GAIN LEFT '
      OPENR, LUN, REFSETGL, /GET_LUN, /SHARE
      FINFO = FSTAT (LUN)
      NUMREC = FINFO.SIZE / FINFO.REC_LEN
      DATAS = ASSOC ( LUN, FEX_GAIN_ST (NUMREC-GOFFSET), GOFFSET )
      DATA = DATAS(0)
      TIM = CONV_ADT_SEC (DATA.CT_HEAD.TIME)
      WHILE (TIM(0) GT TREF1) DO BEGIN
         PRINT, ' GAIN OFFSET OF ',GOFFSET,' IS TOO BIG'
         GOFFSET = GOFFSET - OFFSET
         DATAS = ASSOC ( LUN, FEX_GAIN_ST(NUMREC-GOFFSET), GOFFSET )
         DATA = DATAS(0)
         TIM = CONV_ADT_SEC (DATA.CT_HEAD.TIME)
      ENDWHILE
      TRANGE = WHERE (TIM GE TREF1 AND TIM LE TREF2, COUNT)
      IF (COUNT GT 0) THEN BEGIN
         TRANGE = [TRANGE(0)-1, TRANGE]
         NUM = COUNT + 1
      ENDIF ELSE BEGIN
         TRANGE = WHERE (TIM LT TREF1, COUNT)
         T = TRANGE (COUNT-1)
         TRANGE = [T,T+1]
         NUM = 2
      ENDELSE
      PRINT,' NUMBER OF GAIN LEFT RECORDS IS ',NUM
      EVEN = 2*INDGEN (NUM)
      ODD = EVEN + 1
      NGL = NUM + NUM
      G_TIMEL = DBLARR(NGL,/NOZERO)
      G_TIMEL(EVEN) = CONV_ADT_SEC (DATA(TRANGE).START_TIME)
      G_TIMEL(ODD)  = CONV_ADT_SEC (DATA(TRANGE).STOP_TIME)
      G_GAINL = INTARR(2,NGL,/NOZERO)
      G_GAINL(*,EVEN) = DATA(TRANGE).GAIN
      G_GAINL(*,ODD) = G_GAINL(*,EVEN)
      DATA = 0 & DATAS = 0
      CLOSE, LUN  &  FREE_LUN, LUN
   ENDIF
ENDIF
;
;-----------------Getting Reference Fakeit Data---------------------------------
;
IF (R(0) EQ 1 OR R(0) EQ 2) THEN BEGIN
   FOFFSET = 7 * OFFSET
;--------------------------------Right Channels---------------------------------
   IF (C(0) EQ 0 OR C(0) EQ 1 OR C(0) EQ 4) THEN BEGIN
      PRINT,' GETTING FAKEIT RIGHT DATA'
      OPENR, LUN, REFSETFR, /GET_LUN, /SHARE
      FINFO = FSTAT (LUN)
      NUMREC = FINFO.SIZE / FINFO.REC_LEN
      DATAS = ASSOC ( LUN, FEX_FAKEIT_ST (NUMREC-FOFFSET), FOFFSET )
      DATA = DATAS(0)
      TIM = CONV_ADT_SEC (DATA.CT_HEAD.TIME)
      WHILE (TIM(0) GT TREF1) DO BEGIN
         PRINT, ' FAKEIT OFFSET OF ',FOFFSET,' IS TOO BIG'
         FOFFSET = FOFFSET - OFFSET
         DATAS = ASSOC ( LUN, FEX_FAKEIT_ST (NUMREC-FOFFSET), FOFFSET )
         DATA = DATAS(0)
         TIM = CONV_ADT_SEC (DATA.CT_HEAD.TIME)
      ENDWHILE
      TRANGE = WHERE (TIM GE TREF1 AND TIM LE TREF2,COUNT)
      IF (COUNT GT 0) THEN BEGIN
         TRANGE = [TRANGE(0)-1, TRANGE]
         NUM = COUNT + 1
      ENDIF ELSE BEGIN
         TRANGE = WHERE (TIM LT TREF1, COUNT)
         T = TRANGE(COUNT-1)
         TRANGE = [T,T+1]
         NUM = 2
      ENDELSE
      PRINT,' NUMBER OF FAKEIT RECORDS IS ',NUM
      EVEN = 2*INDGEN (NUM)
      ODD = EVEN + 1
      NFR = NUM + NUM
      F_TIMER = DBLARR(NFR,/NOZERO)
      F_TIMER(EVEN) = CONV_ADT_SEC (DATA(TRANGE).START_TIME)      
      F_TIMER(ODD)  = CONV_ADT_SEC (DATA(TRANGE).STOP_TIME)
      F_FAKER = INTARR(2,NFR,/NOZERO)
      F_FAKER(*,EVEN) = DATA(TRANGE).FAKEIT
      F_FAKER(*,ODD) = F_FAKER(*,EVEN)
      DATA = 0 & DATAS = 0
      CLOSE, LUN  &  FREE_LUN, LUN
   ENDIF
;--------------------------------Left Channels----------------------------------
   IF (C(0) EQ 2 OR C(0) EQ 3 OR C(0) EQ 4) THEN BEGIN
      PRINT,' GETTING FAKEIT LEFT'
      OPENR, LUN, REFSETFL, /GET_LUN, /SHARE
      FINFO = FSTAT (LUN)
      NUMREC = FINFO.SIZE / FINFO.REC_LEN
      DATAS = ASSOC ( LUN, FEX_FAKEIT_ST (NUMREC-FOFFSET),FOFFSET )
      DATA = DATAS(0)
      TIM = CONV_ADT_SEC (DATA.CT_HEAD.TIME)
      WHILE (TIM(0) GT TREF1) DO BEGIN
         PRINT, ' FAKEIT OFFSET OF ',FOFFSET,' IS TOO BIG'
         FOFFSET = FOFFSET - OFFSET
         DATAS = ASSOC ( LUN, FEX_FAKEIT_ST (NUMREC-FOFFSET), FOFFSET )
         DATA = DATAS(0)
         TIM = CONV_ADT_SEC (DATA.CT_HEAD.TIME)
      ENDWHILE      
      TRANGE = WHERE (TIM GE TREF1 AND TIM LE TREF2, COUNT)
      IF (COUNT GT 0) THEN BEGIN
         TRANGE = [TRANGE(0)-1, TRANGE]
         NUM = COUNT + 1
      ENDIF ELSE BEGIN
         TRANGE = WHERE (TIM LT TREF1, COUNT)
         T = TRANGE(COUNT-1)
         TRANGE = [T,T+1]
         NUM = 2
      ENDELSE
      PRINT,' NUMBER OF FAKEIT RECORDS IS ',NUM
      EVEN = 2*INDGEN (NUM)
      ODD = EVEN + 1
      NFL = NUM + NUM
      F_TIMEL = DBLARR(NFL,/NOZERO)
      F_TIMEL(EVEN) = CONV_ADT_SEC (DATA(TRANGE).START_TIME)
      F_TIMEL(ODD)  = CONV_ADT_SEC (DATA(TRANGE).STOP_TIME)
      F_FAKEL = INTARR(2,NFL,/NOZERO)
      F_FAKEL(*,EVEN) = DATA(TRANGE).FAKEIT
      F_FAKEL(*,ODD) = F_FAKEL(*,EVEN)
      DATA = 0  &  DATAS = 0
      CLOSE, LUN  &  FREE_LUN, LUN
   ENDIF
ENDIF
;
;-------------------------PROCEESS SCIENCE CHANNELS-----------------------------
; (J = channel number)
;
FOR I = 0, 3 DO BEGIN
   IF ((C(0)-4)*I EQ 0) THEN BEGIN
      IF (C(0) EQ 4) THEN J = I ELSE J = C(0)
      DATASET = DATATYPE + CHANS(J)
      FILE = ARC + DATASET + FEXT + DAY + "*;"
      A = FINDFILE (FILE,COUNT=NUM)
      IF (NUM EQ 0) THEN $
         PRINT,' NO FILES FOUND OF TYPE : ',FILE $
      ELSE BEGIN
;-------------------------Get Science Data--------------------------------------
         PRINT,' GETTING SCIENCE DATA ',CHANS(J)
         OPENR, LUN, A(0), /GET_LUN, /SHARE
         FINFO = FSTAT (LUN)
         NUMREC = FINFO.SIZE / FINFO.REC_LEN
         NUMR1 = NUMREC - 1
         DATAS = ASSOC ( LUN, FPP_SDF_ST (NUMREC) )
         DATA = DATAS(0)
         TIM = CONV_ADT_SEC (DATA.COLLECT_TIME.MIDPOINT_TIME)
         TRANGE = WHERE (TIM GE TREF1 AND TIM LE TREF2, NUM)
         IF (NUM EQ 0) THEN PRINT,' NO SCIENCE RECORDS IN ',A(0) ELSE BEGIN
            TMID = CONV_ADT_SEC (DATA(TRANGE).COLLECT_TIME.MIDPOINT_TIME)
            S_GAIN = DATA(TRANGE).SCI_HEAD.GAIN
            S_FAKE = DATA(TRANGE).DQ_DATA.FAKE
            BADFLAG = FIX (DATA(TRANGE).COLLECT_TIME.BADTIME_FLAG)
            DATA = 0  &  DATAS = 0
            N1 = NUM - 1
            T1 = TMID(0)
            T2 = TMID(N1)
            IF (T2-T1 GT 14400.) THEN BEGIN
               CONV = 3600.
               UNIT = 'HOURS' & ENDIF $
            ELSE BEGIN
               CONV = 60.
               UNIT = 'MINUTES'
            ENDELSE
            S_TIME = (TMID - T1) / CONV
            X2 = S_TIME(N1)
            TIME1 = CONV_SEC_JDT (T1)
            TIME2 = CONV_SEC_JDT (T2)
            TIME1 = CONV_JDT_JDS(TIME1)
            TIME2 = CONV_JDT_JDS(TIME2)
            XTITL = 'Time in ' + UNIT + ' from ' + TIME1 + ' to ' + TIME2
            CLOSE, LUN  &  FREE_LUN, LUN
;-------------------------Make The Plots----------------------------------------
            !P.MULTI(2) = 3
            !P.CHARSIZE = 1.2
;-------------------------Plot Gain Data----------------------------------------
            IF (R(0) EQ 0 OR R(0) EQ 2) THEN BEGIN
               PLOT, S_TIME, S_GAIN, XRANGE=[0,X2], YRANGE=[-2,+7], PSYM=7, $
                  SYMSIZE=.4, TITLE=DATASET, XTITLE=XTITL, YTITLE='GAIN'
               IF (J EQ 0 OR J EQ 1) THEN BEGIN
                  G_TIME = (G_TIMER - T1) / CONV
                  GAIN = G_GAINR(J,*)
                  NG1 = NGR/2 - 1
               ENDIF ELSE BEGIN
                  G_TIME = (G_TIMEL - T1) / CONV
                  GAIN = G_GAINL(J-2,*)
                  NG1 = NGL/2 - 1
               ENDELSE
               PLOT, G_TIME, GAIN, XRANGE=[0,X2], YRANGE=[-2,+7], PSYM=7, $
                  SYMSIZE=.4,TITLE='REFERENCE DATA', XTITLE=XTITL, YTITLE='GAIN'
               K = 0
               FOR L=0,NG1 DO BEGIN
                  K1 = K + 1
                  OPLOT, G_TIME(K:K1), GAIN(K:K1), XRANGE=[0,X2], $
                     YRANGE=[-2,7],PSYM=-3
                  K = K + 2
               ENDFOR
;-------------------------Plot BadTime Flag-------------------------------------
               PLOT, S_TIME, BADFLAG, XRANGE=[0,X2], YRANGE=[-1,18], PSYM=7, $
                  SYMSIZE=.4, TITLE=DATASET,XTITLE=XTITL, YTITLE='BAD TIME FLAG'
               IF (PLTDEV EQ 'REGIS') THEN WAIT,3
            ENDIF
;-------------------------Plot Fakeit Data--------------------------------------
            IF (R(0) EQ 1 OR R(0) EQ 2) THEN BEGIN
               PLOT, S_TIME, S_FAKE, XRANGE=[0,X2], YRANGE=[-2,+2], PSYM=7, $
                  SYMSIZE=.4, TITLE=DATASET, XTITLE=XTITL, $
                  YTITLE='FAKEIT STATUS'
               IF (J EQ 0 OR J EQ 1) THEN BEGIN
                  F_TIME = (F_TIMER - T1) / CONV
                  FAKE = F_FAKER(J,*)
                  NF1 = NFR/2 - 1
               ENDIF ELSE BEGIN
                  F_TIME = (F_TIMEL - T1) / CONV
                  FAKE = F_FAKEL(J-2,*)
                  NF1 = NFL/2 - 1
               ENDELSE
               PLOT, F_TIME, FAKE, XRANGE=[0,X2], YRANGE=[-2,+2], PSYM=7, $
                  SYMSIZE=.4, TITLE='REFERENCE DATA', XTITLE=XTITL, $
                  YTITLE='FAKEIT STATUS'
               K = 0
               FOR L=0,NF1 DO BEGIN
                  K1 = K + 1
                  OPLOT, F_TIME(K:K1), FAKE(K:K1), XRANGE=[0,X2], $
                     YRANGE=[-2,+2], PSYM=-3
                  K = K + 2
               ENDFOR
;-------------------------Plot BadTime Flag-------------------------------------
               PLOT, S_TIME, BADFLAG, XRANGE=[0,X2], YRANGE=[-1,18], PSYM=7, $
                  SYMSIZE=.4, TITLE=DATASET,XTITLE=XTITL, $
                  YTITLE='BAD TIME FLAG'
               IF (PLTDEV EQ 'REGIS') THEN WAIT,3
            ENDIF
         ENDELSE
      ENDELSE
   ENDIF
ENDFOR	; channel do loop
;-----------------------If Interactive, Ask If User Wants Hard Copy------------
ANS = 'NO '
IF (PLTDEV EQ 'REGIS') THEN BEGIN
   PRINT,' DO YOU WANT A HARD PS COPY? (Y/N)  (only the last channel)'
   READ,ANS
   ANS = STRUPCASE (STRMID (ANS,0,1))
   IF (ANS EQ 'Y') THEN BEGIN
      SET_PLOT, 'PS'
      DEVICE, FILENAME='FIR_FPP.PS',/LANDSCAPE
      IF (R(0) EQ 0 OR R(0) EQ 2) THEN BEGIN
         PLOT, S_TIME, S_GAIN, XRANGE=[0,X2], YRANGE=[-2,+7], PSYM=7, $
            SYMSIZE=.4, TITLE=DATASET, XTITLE=XTITL, YTITLE='GAIN'
         IF (J EQ 0 OR J EQ 1) THEN BEGIN
            G_TIME = (G_TIMER - T1) / CONV
            GAIN = G_GAINR(J,*)
            NG1 = NGR/2 - 1
         ENDIF ELSE BEGIN
            G_TIME = (G_TIMEL - T1) / CONV
            GAIN = G_GAINL(J-2,*)
            NG1 = NGL/2 - 1
         ENDELSE
         PLOT, G_TIME, GAIN, XRANGE=[0,X2], YRANGE=[-2,+7], PSYM=7, $
            SYMSIZE=.4,TITLE='REFERENCE DATA', XTITLE=XTITL, YTITLE='GAIN'
         K = 0
         FOR L=0,NG1 DO BEGIN
            K1 = K + 1
            OPLOT, G_TIME(K:K1), GAIN(K:K1), XRANGE=[0,X2], $
               YRANGE=[-2,7],PSYM=-3
            K = K + 2
         ENDFOR
         PLOT, S_TIME, BADFLAG, XRANGE=[0,X2], YRANGE=[-1,18], PSYM=7, $
            SYMSIZE=.4, TITLE=DATASET,XTITLE=XTITL, YTITLE='BAD TIME FLAG'
      ENDIF
         IF (R(0) EQ 1 OR R(0) EQ 2) THEN BEGIN
         PLOT, S_TIME, S_FAKE, XRANGE=[0,X2], YRANGE=[-2,+2], PSYM=7, $
            SYMSIZE=.4, TITLE=DATASET, XTITLE=XTITL, YTITLE='FAKEIT STATUS'
         IF (J EQ 0 OR J EQ 1) THEN BEGIN
            F_TIME = (F_TIMER - T1) / CONV
            FAKE = F_FAKER(J,*)
            NF1 = NFR/2 - 1
         ENDIF ELSE BEGIN
            F_TIME = (F_TIMEL - T1) / CONV
            FAKE = F_FAKEL(J-2,*)
            NF1 = NFL/2 - 1
         ENDELSE
         PLOT, F_TIME, FAKE, XRANGE=[0,X2], YRANGE=[-2,+2], PSYM=7, $
            SYMSIZE=.4, TITLE='REFERENCE DATA', XTITLE=XTITL, $
            YTITLE='FAKEIT STATUS'
         K = 0
         FOR L=0,NF1 DO BEGIN
            K1 = K + 1
            OPLOT, F_TIME(K:K1), FAKE(K:K1), XRANGE=[0,X2], $
               YRANGE=[-2,+2], PSYM=-3
            K = K + 2
         ENDFOR
         PLOT, S_TIME, BADFLAG, XRANGE=[0,X2], YRANGE=[-1,18], PSYM=7, $
            SYMSIZE=.4, TITLE=DATASET,XTITLE=XTITL, YTITLE='BAD TIME FLAG'
      ENDIF
   ENDIF
ENDIF
;
IF(PLTDEV EQ 'TEK')THEN BEGIN
 DEVICE,PLOT_TO=0
 CLOSE,LUNP
ENDIF
IF(PLTDEV EQ 'PS')THEN DEVICE,/CLOSE
;
STOP
END
