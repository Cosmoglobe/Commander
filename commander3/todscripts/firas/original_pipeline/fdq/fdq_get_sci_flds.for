	INTEGER*4 FUNCTION FDQ_GET_SCI_FLDS (SCI_REC, NCHAN, 
	1	VALID_CHAN, IDX_REC)
C/
C/	PROGRAM NAME:
C/	  FDQ_GET_SCI_FLDS
C/
C/	PROGRAM DESCRIPTION:
C/	  This routine picks up those Index fields that are found in
C/	  the Science files.
C/
C/	AUTHOR:
C/	  Edwin H. Fung
C/	  GSFC
C/	  April 3, 1986
C/
C/	MODIFIED BY:
C/
C/	  Shirley M. Read
C/	  STX
C/	  January 14, 1988
C/	  REASON: 	Converted from subroutine to function for interface
C/			with Fut_Error condition handler. Added error checking
C/		        and calls to Lib$Signal.
C/
C/        Harte Y. Wang
C/        STX
C/        Dec. 12, 1990
C/        SPR 7855, Missing Gain data in FDQ_IDX record
C/
C/        Harte Y. Wang
C/        STX
C/        Jan. 29, 1991
C/        Reason:   New requirements for FDQ
C/                Get the Gain data from FDQ_SDF_XX        
C/  
C/	CALLING SEQUENCE:
C/	  Status = FDQ_GET_SCI_FLDS (SCI_REC, VALID_CHAN, IDX_REC)
C/
C/	INPUT PARAMETERS:
C/	  SCI_REC(4)	BYTE	Buffer containing the current records
C/					for the 4 Science files
C/	  NCHAN			I*2	Number of channels with valid data in
C/					this pass
C/	  VALID_CHAN(4)		I*2	Array listing the channels with
C/					valid data in this pass. 1 = RH,
C/					2 = RL, 3 = LH, 4 = LL.
C/
C/	OUTPUT PARAMETERS:
C/	  IDX_REC(512)		BYTE	The index buffer containing the fields
C/					from the Science files filled in when
C/					the routine exits.  However, this routine
C/					should not modify other fields in the
C/					Index record buffer which are not related
C/					to the Science files
C/
C/	INPUT/OUTPUT FILES:
C/	  NONE
C/
C/	SUBROUTINES CALLED:
C/	  TBD
C/
C/	INCLUDE FILES USED:
C/	  FIR_IDX_PARMS.INC
C/
C/	ERROR HANDLING:
C/	  TBD
C/
C/	METHOD USED:
C/	  TBD
C/

	IMPLICIT	NONE

        dictionary 'fdq_idx'
        record /fdq_idx/  idx_rec

        dictionary 'nfs_sdf'
        record /nfs_sdf/   sci_rec(4)

	EXTERNAL   	FDQ_NORMAL
	EXTERNAL	FDQ_ABERR

	integer*4 	RETSTAT		! Return status
	integer*4 	SUCCESS / 1 /, ERROR / 2 /  ! Values for status

	INTEGER*2	NCHAN, VALID_CHAN(4), Gain(4)

	INTEGER*2	CH, GAIN_VAL(0:7) 
	1		/ 1, 3, 10, 30, 100, 300, 1000, 3000/,
	1		I

!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!			   !
!     Code begins here     !
!			   !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!


C	Set return status to success.

	RETSTAT = SUCCESS


	DO	I = 1, NCHAN
!
!     NOTE:
!     Only update those science channels that are valid in this pass.
!     For those which are not valid, keep it the same as the previous
!     pass!
!
	  CH = VALID_CHAN(I)

!!!	  IDX_REC.chan(ch).up_sci_mode = SCI_REC(CH).sci_head.sc_head1a
!!!	  IDX_REC.chan(ch).up_adds_per_grp = SCI_REC(CH).sci_head.sc_head9
!!!	  IDX_REC.chan(ch).up_sweeps_per_ifg = SCI_REC(CH).sci_head.sc_head11
!!!	  IDX_REC.chan(ch).mtm_speed_at_xmit = SCI_REC(CH).sci_head.mtm_speed
!!!	  IDX_REC.chan(ch).mtm_length_at_xmit = SCI_REC(CH).sci_head.mtm_length

!     The above 5 lines have to be changed to the following because of
!     mismatching datatypes:  the fields in IDX_REC are byte-fields while 
!     those in SCI_REC are word-fields

	  call lib$movc3 (1, sci_rec(ch).sci_head.sc_head1a, 
	1	idx_rec.chan(ch).up_sci_mode)
	  call lib$movc3 (1, sci_rec(ch).sci_head.sc_head9,
	1	IDX_REC.chan(ch).up_adds_per_grp)
	  call lib$movc3 (1, SCI_REC(CH).sci_head.sc_head11,
	1	IDX_REC.chan(ch).up_sweeps_per_ifg)
	  call lib$movc3 (1, SCI_REC(CH).sci_head.mtm_speed,
	1	IDX_REC.chan(ch).mtm_speed_at_xmit)
	  call lib$movc3 (1, SCI_REC(CH).sci_head.mtm_length,
	1	IDX_REC.chan(ch).mtm_length_at_xmit)

!
!     Move converted gains into 2-byte field in Index buffer
!
          If (sci_rec(ch).sci_head.gain .ge. 0) then
            idx_rec.chan(ch).science_gain = gain_val(sci_rec(ch).sci_head.gain)
          else
            idx_rec.chan(ch).science_gain = -1
          Endif

	ENDDO

!	Set function to return status

	IF (RETSTAT.EQ.SUCCESS) THEN
	  FDQ_GET_SCI_FLDS = %loc(FDQ_NORMAL)
	ELSE
	  FDQ_GET_SCI_FLDS = %loc(FDQ_ABERR)
	ENDIF

	RETURN
	END
