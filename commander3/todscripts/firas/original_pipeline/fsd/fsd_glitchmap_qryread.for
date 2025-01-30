	Subroutine FSD_Glitchmap_QryRead(ct_nr, rse, sci_rec, scitype,
	2                                i_stop, s_flag)
C-----------------------------------------------------------------------------
C	This program uses query routines to allow the user to select a time
C	 range to open the selected data set.
C
C  CALLING ARGUMENTS:
C    INPUT:
C	 I*  4    ct_nr
C	Ch*128    rse(16)
C   Rec /NFS_SDF/ sci_rec
C	Ch*  3    scitype     ! NFS (default), FPP, or FDQ.  From cmd line.
C	 I*  4    chan_id     ! RH=1, RL=2, LH=3, LL=4
C	Ch* 30    time_range
C	 L*  4    begin
C    OUTPUT:
C	 L*  4    i_stop
C	 L*  4    s_flag
C	 L*  4    eof         !end of file flag
C	 L*  4    more
C
C-----------------------------------------------------------------------------
C    Changes:
C
C       VERSION 4.2.1, SER 2375, QUOC CHUNG STX, 12/29/88
C       USE QUERY BUILDER TO ACCESS RAW SCIENCE USING RSE.
C
C	SPR 3945, 4178, Allow glitchmap to work on the new FPP_SDF and FDQ_SDF
C	    files, as well as the raw science, NFS_SDF, so that it can run
C	    whether or not FPP or FDQ have.  FSD_Glitchmap, _Parse,
C	    _Init, _QryRead, and FSD.CLD.  Fred Shuman, STX / 1989 Aug 29.
C------------------------------------------------------------------------------

	Implicit None

	Include 'CT$Library:CTUser.inc'
	Include '(CCT_Status_Record)'
	Include '(FSD_Glitchmap)'
	Include '($SSDEF)'

	Record /CCT_Status/ struc

	Dictionary 'CDU_DAFS'
	Record /CDU_DAFS/ dafs_rec !output DAFS record

	Integer   *  4      ct_nr
	Character *128      rse(16)
	Character *  3      scitype
	Logical   *  4      i_stop
	Logical   *  4      s_flag

	Logical   *  4      i_chan
	Logical   *  1      oldrse /.False./

	Character * 64      filespec
	Character *  2      fac_channel_id(4)
	Character * 14      start_time, stop_time

	Integer   *  2      ios           !io status

	Integer   *  4      rstatus
	Integer   *  4      status
	Integer   *  4      iostat
	Integer   *  4      len
C
C       functions
C
	Integer   *  4      LIB$Get_Lun
	Integer   *  4      CT_Connect_Query_Nonraw

	Data   fac_channel_id /'RH','RL','LH','LL'/

	External            FSD_CTQGet
	External            FSD_CTQEOF
	External            FSD_CTOpen

	External            CT_Connect_Query_Nonraw

C
C  Get a logical unit number:
C
	If (begin) Then
	   status = LIB$Get_Lun(ct_nr)
	   If (status .Ne. SS$_Normal) Then
	      Call Exit (status)
	   End If
C
C  Open the selected archive file:  either NFS_SDF, FPP_SDF, or FDQ_SDF.
C
	   filespec = 'CSDR$FIRAS_ARCHIVE:' // scitype // '_SDF_'
	2             // fac_channel_id(chan_id) // '/' // time_range
	   Open (UNIT =     ct_nr,
	3        FILE =     filespec,
	4        STATUS =   'old',
	5        USEROPEN = CT_Connect_Query_Nonraw,
	6        IOSTAT =   rstatus )

	   If (rstatus .Ne. 0) Then
	      Call LIB$Signal( FSD_CTOpen, %Val(2), %Val(chan_id),
	2                      %Val(rstatus) )
	      more = .False.
	      i_stop = .True.
	      s_flag = .True.
	      Return

	   End If

	   Call CT_Query_Arcv(, ct_nr, rse, struc)
	   begin = .False.
	End If  !begin
C
C  Get a record:
C
	If (more) Then
	   Call CT_Query_Get(, ct_nr, sci_rec, struc)
	   If (struc.cterr .Eq. CTP_Normal) Then
	      more = .True.
	      eof = .False.
	   Else If (struc.cterr .Eq. CTP_EndOfFile) Then
	      more = .False.
	      eof = .True.
	      rstatus = struc.cterr
	      Call LIB$Signal(FSD_CTQEOF, %Val(1), %Val(rstatus))
	      Call CT_Close_Arcv(, ct_nr, struc)
	   Else
	      rstatus = struc.cterr
	      Call LIB$Signal(FSD_CTQGet, %Val(1), %Val(rstatus))
	      Call CT_Close_Arcv(, ct_nr, struc)
	      status = SS$_Abort
	      Call Exit(status)
	   End If
	End If  ! more
	Return
	End
