	Subroutine FSD_Glitchmap_Init(ct_nr, rse)
C------------------------------------------------------------------------------
C	This program uses the Query Builder to allow the user to select a time
C	 range for the selected data set.
C
C------------------------------------------------------------------------------
C   Changes:
C
C       VERSION 4.2.1, SER 2375, QUOC CHUNG STX, 12/22/88
C       USE QUERY BUILDER TO ACCESS RAW SCIENCE USING RSE.
C
C	SPR 3945, 4178, Allow glitchmap to work on the new FPP_SDF and FDQ_SDF
C	    files, as well as the raw science, NFS_SDF, so that it can run
C	    whether or not FPP or FDQ have.  FSD_Glitchmap, _Parse,
C	    _Init, _QryRead, and FSD.CLD.  Fred Shuman, STX / 1989 Aug 29.
C
C       SPR 9847, Update FSD_Glitchmap to use new FEX_SAMPRATE. Move initialize
C           cobetrieve statement to the main program due to cct_open_fig.
C           Nilo G. Gonzales/Hughes STX, 1992 Aug 5.
C------------------------------------------------------------------------------

	Implicit None

	Include 'CT$Library:CTUser.inc'
	Include '(CCT_Status_Record)'
	Include '(FSD_Glitchmap)'
	Include '($SSDef)'

	Record /CCT_Status/ struc

	Integer   *  4      ct_nr
	Character * 128     rse(16)

	Logical   *  4      i_chan
	Logical   *  1      oldrse /.False./

	Integer   *  2      dbid        ! database ID for CT_Query_Bld
	Integer   *  2      ios         ! io status

	Integer   *  4      rstatus
	Integer   *  4      status
	Integer   *  4      lun3
	Integer   *  4      unit

	Character * 14      start_time, stop_time

	Logical   *  4      i_stop

	External            FSD_CTQBld
	External            FSD_FOpen

	begin = .True.
	dbid = CTU_$FIRAS

	If (begin) Then
C
C *** Use the Query Builder to get a timerange and RSEs
C
	   Call CT_Query_Bld(dbid, time_range, rse, struc, oldrse)
	   If (struc.cterr .Ne. CTP_Normal) Then
	      rstatus = struc.cterr
	      Call LIB$Signal(FSD_CTQBld, %Val(2), %Val(chan_id),
	2                     %Val(rstatus))
	      status = SS$_Abort
	      Call Exit(status)
	   End If        ! (struc.cterr .Ne. CTP_Normal)
	End If           ! (begin)

	Return
	End
