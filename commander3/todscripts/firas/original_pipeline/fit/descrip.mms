! Builds facility FIT
!
! Author: Fred Shuman
!	  STX			Last revision: 1987 Nov 20
!	  October 23, 1987		to remove fit.opt
!
! 1988 Sep  7:  Removed references to FUTI.  LO/SMR
! 1989 Oct 24:  Added references to new module FIT_Combine_Temps.  FGS
! 1990 Feb 12:  SPR 6012, VMS 5.2 options file changes.  RJK
! 1991 Sep 27:  Removed references to module FIT_Combine_Temps, which is now
!               FUT_Combine_HiLo.  FGS
!
.suffixes
.suffixes .exe .olb .obj .for .cld .msg .tlb .txt

fit : fitbld.tlb, fit.exe
  ! Facility FIT is up to date.

! Define macros.

text_libs = fitbld.tlb/lib + csdr$library:futlib.tlb/lib + -
			     csdr$library:csdrlib.tlb/lib

.FOR.OBJ
  $(fort)$(fflags)/extend_source $(mms$source) + $(text_libs)

inc_files =	fit_invoc.txt, -
		fit_menu.txt, -
		fit_data.txt

! Build FIT (INSTRUMENT_TRENDS) program.

trendplots_obj = fit_main.obj, -
		 fit_init.obj, -
		 fit_get_command.obj, -
		 fit_display_menu.obj, -
		 fit_display_main.obj, -
		 fit_display_plot_menu.obj, -
		 fit_match_fields.obj, -
		 fit_extract_eng_stats.obj, -
		 fit_merge_eng_stats.obj, -
		 fit_plot.obj, -
		 fit_report.obj, -
		 fit_exit.obj, -
		 fit_msg.obj

! Source code dependencies.

fit_main.obj            	: fit_main.for, -
				  fitbld.tlb(fit_invoc), -
				  fitbld.tlb(fit_menu), -
				  fitbld.tlb(fit_data), -
	                          csdr$library:futlib.tlb, -
				  ct$library:ctuser.inc, -
				  fdq_etr^

fit_init.obj			: fit_init.for, -
				  fitbld.tlb(fit_invoc), -
	                          csdr$library:futlib.tlb

fit_get_command.obj		: fit_get_command.for, -
				  fitbld.tlb(fit_invoc), -
	                          csdr$library:futlib.tlb, -
                                  ct$library:ctuser.inc


fit_display_menu.obj		: fit_display_menu.for, -
	                          csdr$library:futlib.tlb, -
				  fitbld.tlb(fit_menu)

fit_display_main.obj		: fit_display_main.for, -
	                          csdr$library:futlib.tlb, -
				  fitbld.tlb(fit_menu)

fit_display_plot_menu.obj	: fit_display_plot_menu.for, -
	                          csdr$library:futlib.tlb, -
				  fitbld.tlb(fit_menu)

fit_match_fields.obj		: fit_match_fields.for, -
	                          csdr$library:futlib.tlb, -
				  fitbld.tlb(fit_menu), -
				  fitbld.tlb(fit_invoc)

fit_extract_eng_stats.obj	: fit_extract_eng_stats.for, -
				  fitbld.tlb(fit_invoc), -
				  fitbld.tlb(fit_menu), -
				  fitbld.tlb(fit_data), -
	                          csdr$library:futlib.tlb

fit_merge_eng_stats.obj		: fit_merge_eng_stats.for, -
				  fitbld.tlb(fit_invoc), -
				  fitbld.tlb(fit_menu), -
	                          csdr$library:futlib.tlb

fit_plot.obj			: fit_plot.for, -
				  fitbld.tlb(fit_menu)

fit_report.obj			: fit_report.for, -
				  fitbld.tlb(fit_invoc), -
	                          csdr$library:futlib.tlb

fit_exit.obj			: fit_exit.for, -
				  fitbld.tlb(fit_invoc), -
	                          csdr$library:futlib.tlb

! Messages...

fit_msg.obj : fit_msg.msg
  MESSAGE fit_msg


! Build the libraries.

fitbld.tlb : $(inc_files)
  LIBRARY/CREATE/TEXT fitbld $(inc_files)
  ! FITBLD.TLB has been built.

fit.exe :               fitbld.olb( $(trendplots_obj) ),-
	                csdr$library:futlib.olb,-
	                csdr$library:futlib.tlb,-
	                xanadu:[lib]vilib.olb,-
	                xanadu:[lib]xlib.olb,-
	                graphics:grpshr.exe,-
	                csdr$library:v5.opt
  $(link) $(linkflags)  fitbld.olb/lib/inc=(fit_main), -
	                csdr$library:futlib.olb/lib, -
	                csdr$library:csdrlib.olb/lib, -
			csdr$library:csdrmsg.olb/inc=(upm_msg,cut_msg),-
	                csdr$library:igselib.olb/lib, -
	                xanadu:[lib]vilib.olb/lib, -
	                xanadu:[lib]xlib.olb/lib, -
	                graphics:grpshr/lib, -
	                csdr$library:v5.opt/option

  ! FIT_TRENDPLOTS has been built.
