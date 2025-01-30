! Builds facility FEP_ENGPLOTS and FEP_DWELLPLOTS.
!
! Author: Rob Kummerer
!	  March 5, 1987
!	  STX
! Change Log:
!
!	  Shirley M. Read
!	  September 9, 1988
!	  STX
!	  Reason: Remove reference to obsolete Futilib.Olb.
!
!	  R. Kummerer
!	  February 12, 1990
!	  STX
!	  Reason: SPR 5795, Remove dependency on TEMPLATE.
!	          SPR 6008, VMS 5.2 options file changes.
!
!	  H. Wang
!	  February 22, 1990
!	  STX
!	  Reason: SPR 3921, REDESIGN VERSION, FEP READ the reference data
!                 from the reference archive. 
!
!	  R. Kummerer
!	  March 2, 1990
!	  STX
!	  Reason: SPR 6354, Prevent MMS abort on warning "Multiply defined
!		  symbol WAIT".
!
!	  H. Wang
!	  Mar 2, 1990
!	  STX
!	  Reason: SPR 3253, REDESIGN VERSION, FEP/Dwell, READ the reference data
!                 from the reference archive. 
!
.suffixes
.suffixes .exe .olb .obj .for .cld .msg .tlb .txt

fep : fepbld.tlb, fep_engplots.exe, fep_dwellplot.exe
  ! Facility FEP is up to date.

! Define macros.

text_libs = fepbld.tlb/lib + csdr$library:futlib.tlb/lib + -
			     csdr$library:csdrlib.tlb/lib

.for.obj
  $(fort)$(fflags)/extend $(mms$source) + $(text_libs)

inc_files =	fep_invoc.txt, -
		fep_menu.txt, -
		fep_firnames.txt, -
		fep_data.txt

! Build FEP_ENGPLOTS program.

engplots_obj =	fep_engplots.obj, -
		fep_init.obj, -
		fep_get_command.obj, -
		fep_display_menu.obj, -
		fep_display_main.obj, -
		fep_display_plot_menu.obj, -
		fep_match_fields.obj, -
		fep_convert_to_eng.obj, -
		fep_get_grt_attr.obj, -
		fep_get_grt_conv.obj, -
		fep_fix_counts_to_ohms.obj, -
		fep_counts_to_ohms.obj, -
		fep_get_calres.obj, -
		fep_decode_dbname.obj, -
		fep_grt_lookup.obj, -
		fep_invert.obj, -
		fep_get_curve.obj, -
		fep_plot.obj, -
		fep_convolve.obj, -
		fep_report.obj, -
		fep_exit.obj, -
		fep_msg.obj

! Source code dependencies.

fep_engplots.obj		: fep_engplots.for, -
				  fepbld.tlb(fep_invoc), -
				  fepbld.tlb(fep_menu), -
				  fepbld.tlb(fep_firnames), -
				  fepbld.tlb(fep_data), -
	                          csdr$library:futlib.tlb, -
				  ct$library:ctuser.inc

fep_init.obj			: fep_init.for, -
				  fepbld.tlb(fep_invoc), -
	                          csdr$library:futlib.tlb

fep_get_command.obj		: fep_get_command.for, -
				  fepbld.tlb(fep_invoc), -
	                          csdr$library:futlib.tlb, -
                                  ct$library:ctuser.inc


fep_display_menu.obj		: fep_display_menu.for, -
	                          csdr$library:futlib.tlb, -
				  fepbld.tlb(fep_menu)

fep_display_main.obj		: fep_display_main.for, -
	                          csdr$library:futlib.tlb, -
				  fepbld.tlb(fep_menu)

fep_display_plot_menu.obj	: fep_display_plot_menu.for, -
	                          csdr$library:futlib.tlb, -
				  fepbld.tlb(fep_menu)

fep_match_fields.obj		: fep_match_fields.for, -
	                          csdr$library:futlib.tlb, -
				  fepbld.tlb(fep_menu), -
				  fepbld.tlb(fep_invoc)

fep_convert_to_eng.obj		: fep_convert_to_eng.for, -
				  fepbld.tlb(fep_menu), -
				  fepbld.tlb(fep_firnames), -
				  fepbld.tlb(fep_data), -
	                          csdr$library:futlib.tlb

fep_get_grt_attr.obj		: fep_get_grt_attr.for, -
	                          csdr$library:futlib.tlb

fep_get_grt_conv.obj		: fep_get_grt_conv.for, -
	                          csdr$library:futlib.tlb

fep_fix_counts_to_ohms.obj	: fepbld.tlb(fep_data), -
	                          csdr$library:futlib.tlb

fep_counts_to_ohms.obj		: fep_counts_to_ohms.for

fep_get_calres.obj		: fep_get_calres.for

fep_decode_dbname.obj		: fep_decode_dbname.for

fep_grt_lookup.obj		: fep_grt_lookup.for

fep_invert.obj			: fep_invert.for

fep_get_curve.obj		: fep_get_curve.for

fep_plot.obj			: fep_plot.for, -
				  fepbld.tlb(fep_menu)

fep_convolve.obj		: fep_convolve.for

fep_report.obj			: fep_report.for, -
				  fepbld.tlb(fep_invoc), -
	                          csdr$library:futlib.tlb

fep_exit.obj			: fep_exit.for, -
				  fepbld.tlb(fep_invoc), -
	                          csdr$library:futlib.tlb

! Messages...

fep_msg.obj : fep_msg.msg
  MESSAGE fep_msg


! Build the libraries.

fepbld.tlb : $(inc_files)
  LIBRARY/CREATE/TEXT fepbld $(inc_files)
  ! FEPBLD.TLB has been built.

fep_engplots.exe :	fepbld.olb( $(engplots_obj) ),-
			csdr$library:futlib.olb,-
	                xanadu:[lib]vilib.olb,-
	                xanadu:[lib]xlib.olb,-
	                graphics:grpshr.exe,-
			imsl$dir:imsl.olb, -
			csdr$library:v5.opt
 - $(link) $(linkflags)	fepbld.olb/lib/inc=(fep_engplots), -
	  		csdr$library:futlib.olb/lib, -
			csdr$library:csdrlib.olb/lib,-
			csdr$library:csdrmsg.olb/inc=(upm_msg,cut_msg),-
	                xanadu:[lib]vilib.olb/lib, -
	                xanadu:[lib]xlib.olb/lib, -
	                graphics:grpshr/lib, -
			imsl$dir:imsl.olb/lib, -
			csdr$library:v5.opt/option
  !
  ! The "Multiply defined symbol WAIT" warning may be ignored.
  !
  ! FEP_ENGPLOTS has been built.


!Build FEP_DWELLPLOTS program.
!		cnts2ohms, -
!		read_firas_hkp, -
!cnts2ohms.obj		: cnts2ohms.for
!read_firas_hkp.obj	: read_firas_hkp.for, ct$library:ctuser.inc, nfs_hkp^

dwellplot_obj =	fep_dwellplot, -
		fep_convert_to_ohms, -
		fep_dmakeplot, -
		fep_dplotter, -
		fep_dtemperatures, -
		fep_get_grt_conv, -
		fep_dwell_address, -
		fep_grtcoeffs, -
		fep_grt_lookup, -
		fep_get_calres, -
		fep_read_firas_hkp

fep_dwellplot.obj	 : fep_dwellplot.for, nfs_hkp^
fep_convert_to_ohms.obj	 : fep_convert_to_ohms.for
fep_dmakeplot.obj	 : fep_dmakeplot.for
fep_dplotter.obj	 : fep_dplotter.for
fep_dtemperatures.obj	 : fep_dtemperatures.for, -
		           fepbld.tlb(fep_dwell_dbword)
fep_get_grt_conv.obj	 : fep_get_grt_conv.for
fep_dwell_address.obj	 : fep_dwell_address.for
fep_grtcoeffs.obj	 : fep_grtcoeffs.for
fep_grt_lookup.obj	 : fep_grt_lookup.for
fep_get_calres.obj       : fep_get_calres.for
fep_read_firas_hkp.obj	 : fep_read_firas_hkp.for

fep_dwellplot.exe :  fepbld.olb($(dwellplot_obj)), -             
		     csdr$library:futlib.olb,-
	             xanadu:[lib]vilib.olb,-
	             xanadu:[lib]xlib.olb,-
	             graphics:grpshr.exe,-
                     imsl$dir:imsl.olb,-
                     csdr$library:v5.opt
 - $(link)$(linkflags)	fepbld.olb/lib/inc=(fep_dwellplot), -
		        csdr$library:futlib.olb/lib,-
	                xanadu:[lib]vilib.olb/lib, -
	                xanadu:[lib]xlib.olb/lib, -
	                graphics:grpshr/lib, -
                        imsl$dir:imsl.olb/lib, -
                        csdr$library:csdrlib.olb/lib, -
                        csdr$library:v5.opt/option
 !
 ! The "Multiply defined symbol WAIT" warning may be ignored.
 !
 ! FEP_DWELLPLOT has been built.
