	Integer * 4 Function  FFL_Get_Reference (ref_open_dir, ref_open_seq,
     &                                           fflc, ffli)

c-------------------------------------------------------------------------------
c
c	Function FFL_GET_REFERENCE
c
c	This function opens and gets the reference datasets for FFL via the
c	Get_Config routines.
c
c	Author:   Gene Eplee
c		  General Sciences Corp.
c		  513-7768
c		  22 September 1992
c		  SER 10763
c
c-------------------------------------------------------------------------------
c
c	Input:
c	    fflc.* \	These structures are declared in the
c	    ffli.* /	 FFL_Config and FFL_Invoc include files
c
c	Output:
c	    ref_open_dir	L * 1	open status flag
c	    ref_open_seq	L * 1	open status flag
c
c	Subroutines called:
c		cct_get_config_idx_tod
c		cct_get_config_tod
c		cct_open_config
c		ct_gmt_to_binary
c		LIB$Signal
c
c	Include files:
c		ffl_config.txt
c		ffl_invoc.txt
c		fut_params.txt
c
c-------------------------------------------------------------------------------
c
c   Changes:
c	Name of facility changed from FFI to FFL.  Fred Shuman, HSTX
c	1995 May 19.
c
c	In order to remove some of the obscurity arising from Include files that
c	   harbor hidden Common blocks, converted these Commons to Structures.
c	   This necessitates adding their Record names to the calling lists of
c	   functions that use their variables.
c	Fred Shuman, HSTX, 1995 June 14.
c
c	Removed include file ffl_spec.txt; it was needed only for defining a
c	   couple of variables that were passed through this routine between the
c	   extinct routine FFL_Compute_Constants and the main program; these
c	   variables are now defined in-line in the only routine that uses them,
c	   FFL_Produce_Spectra.
c	Fred Shuman, HSTX, 1995 June 30.
c
c	Subsumed the FFL_SPEC.TXT include file into FFL_CONFIG.TXT .
c	Fred Shuman, HSTX, 1995 July 21.
c-------------------------------------------------------------------------------

	Implicit None

	Include '(fut_params)'
c
c  Call arguments:
c
	Logical  * 1	ref_open_dir	!  reference file open flag
	Logical  * 1	ref_open_seq	!  reference file open flag
	Include '(ffl_config)'		! defines record structure fflc and
		!  certain other vbles; Includes CCT_Get_Config; Dictionary
		!  structures FEX_GRTCOAWT, FEX_GRTTRANS.
		!defines structure fishspec and parameters vrecsize and vspecsiz
	Include '(ffl_invoc)'		! defines record structure ffli
c
c  All other variables and functions:
c
	Integer  * 4	bin_time(2)	!  ADT time string
	Integer  * 4	cstatus		!  return status
	Integer  * 4	status		!  return status

	Integer  * 4	cct_get_config_idx_tod
	Integer  * 4	cct_get_config_tod
	Integer  * 4	cct_open_config
	Real	 * 8	uoe_adt2t68

	External	ffl_cfgdirget
	External	ffl_cfgdiropen
	External	ffl_cfgseqget
	External	ffl_cfgseqopen
	External	ffl_normal
	External	cct_normal

	status = %Loc(ffl_normal)
c
c  Get the aperture cover eject time and covert it to T68 seconds.
c
	Call ct_gmt_to_binary(fac_apco_gmt,bin_time)
	fflc.apco_eject_time = uoe_adt2t68(bin_time)
c
c  Open the configuration files for sequential and direct access.
c
	ref_open_dir = .False.
	ref_open_seq = .False.
	Call ct_gmt_to_binary(ref_gmt_start, fflc.ref_start)
	Call ct_gmt_to_binary(ref_gmt_stop, fflc.ref_stop)

	cstatus = cct_open_config (fflc.ref_start, fflc.ref_stop, ndset_tod,
     &                             dset_tod, size_tod, seq_access,
     &                             ncache, fflc.tod_lun, fflc.tod_index,
     &                             fflc.tod_stat, ref_count)
	If (cstatus .Ne. %Loc(cct_normal)) Then
	   status = %Loc(ffl_cfgseqopen)
	   Call LIB$Signal (ffl_cfgseqopen, %Val(1), %Val(cstatus))
	Else
	   ref_open_seq = .True.

	   cstatus = cct_open_config (fflc.ref_start, fflc.ref_stop, ndset_dir,
     &                                dset_dir, size_dir, dir_access,
     &                                ncache, fflc.dir_lun, fflc.dir_index,
     &                                fflc.dir_stat, ref_count)
	   If (cstatus .Ne. %Loc(cct_normal)) Then
	      status = %Loc(ffl_cfgdiropen)
	      Call LIB$Signal (ffl_cfgdiropen, %Val(1), %Val(cstatus))
	   Else
	      ref_open_dir = .True.
	   EndIf
	EndIf

	If (status .Eq. %Loc(ffl_normal)) Then
c
c  Get the reference datasets.
c
	   cstatus = cct_get_config_tod (ffli.jstart, ndset_tod, size_tod,
     &				fflc.tod_lun, fflc.tod_index, fflc.config,
     &				new_tod_segment, fflc.tod_stat)
	   If (cstatus .Ne. %Loc(cct_normal)) Then
	      status = %Loc(ffl_cfgseqget)
	      Call LIB$Signal (ffl_cfgseqget, %Val(1), %Val(cstatus))
	   Else

	      cstatus = cct_get_config_idx_tod (ffli.jstart, ndset_dir,
     &						fflc.dir_lun, fflc.dir_index,
     &						new_dir_segment, fflc.dir_stat)
	      If (cstatus .Ne. %Loc(cct_normal)) Then
	         status = %Loc(ffl_cfgdirget)
	         Call LIB$Signal (ffl_cfgdirget, %Val(1), %Val(cstatus))
	      EndIf
	   EndIf

	EndIf	!   (status from Open_Config


	FFL_Get_Reference = status

	Return
	End
