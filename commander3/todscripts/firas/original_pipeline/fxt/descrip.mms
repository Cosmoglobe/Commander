! Builds facility FXT_LOG_EXTREMA.
!
! Author: Shirley M. Read
!	  January 10, 1989
!	  STX Incorporated
!
! Changes:
!	SPR 6023, February 12, 1990, R. Kummerer. VMS 5.2 options file changes.
!

fxt : fxtbld.tlb, fxt.exe
  ! Facility FXT is up to date.

inc_files = fxt_msg.txt

text_libs = fxtbld.tlb/lib + csdr$library:futlib.tlb/lib + -
			     csdr$library:csdrlib.tlb/lib

.for.obj
  $(fort)/extend_source $(fflags) $(mms$source) + $(text_libs)

fxt_object =	fxt_log_extrema.obj,-
		fxt_parse_command.obj,-
		fxt_init_report.obj,-
		fxt_cat_info.obj,-
		fxt_get_minmax.obj,-
		fxt_open_archive.obj,-
		fxt_get_flags.obj,-
		fxt_compare_eng.obj,-
		fxt_trigger_plot.obj,-
		fxt_write_extrema.obj,-
		fxt_write_engplot.obj,-
		fxt_close_xtrm.obj,-
		fxt_msg.obj


fxt_log_extrema.obj     : fxtbld.tlb(fxt_msg),-
                          ct$library:ctuser.inc,-
			  fdq_eng^,-
			  fxt_eng_xtrm^
fxt_parse_command.obj 	: fxtbld.tlb(fxt_msg)
fxt_init_report.obj 	: fxtbld.tlb(fxt_msg)
fxt_cat_info.obj 	: fxtbld.tlb(fxt_msg),-
			  ct$library:ctuser.inc
fxt_get_minmax.obj 	: fxtbld.tlb(fxt_msg),-
                          ct$library:ctuser.inc,-
			  fxt_eng_xtrm^
fxt_open_archive.obj 	: fxtbld.tlb(fxt_msg),-
                          ct$library:ctuser.inc
fxt_get_flags.obj	: fxtbld.tlb(fxt_msg),-
			  csdr$library:futlib.tlb,-
			  ct$library:ctuser.inc,-
			  fex_limflags^
fxt_compare_eng.obj 	: fxtbld.tlb(fxt_msg),-
			  fdq_eng^,-
			  fxt_eng_xtrm^
fxt_trigger_plot.obj 	: fxtbld.tlb(fxt_msg)
fxt_write_extrema.obj	: fxtbld.tlb(fxt_msg),-
			  csdr$library:futlib.tlb,-
			  ct$library:ctuser.inc,-
			  fxt_eng_xtrm^
fxt_write_engplot.obj	: fxtbld.tlb(fxt_msg),-
			  csdr$library:futlib.tlb
fxt_close_xtrm.obj 	: fxtbld.tlb(fxt_msg),-
			  ct$library:ctuser.inc

fxt_msg.obj : fxt_msg.msg
  MESSAGE fxt_msg

fxtbld.tlb : $(inc_files)
  LIBRARY/CREATE/TEXT fxtbld $(inc_files)
  ! FXTBLD.TLB has been built.

fxt.exe :	fxtbld.olb( $(fxt_object) ),-
		fxtbld.tlb,-
		csdr$library:futlib.olb,-
		csdr$library:v5.opt
  $(link) $(linkflags) fxtbld.olb/lib/inc=(fxt_log_extrema),-
                       csdr$library:futlib.olb/lib,-
	   	       csdr$library:csdrlib.olb/lib,-
	   	       csdr$library:csdrmsg.olb/inc=(upm_msg), -
	   	       csdr$library:v5.opt/option
  ! FXTBLD.OLB has been built.

