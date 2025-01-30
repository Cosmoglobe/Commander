! Builds facility FPP_FIRAS_PREPROCESSOR
!
! Author: QUOC C CHUNG
!	  MARCH 3, 1989
!	  STX Incorporated
!
! Changes:
!	SPR 6014, February 12, 1990, R. Kummerer, VMS 5.2 options file changes.
!	New design, February-April 1991, Larry P. Rosen, STX

fpp : fppbld.tlb, fpp.exe
  ! Facility Fpp is up to date.

inc_files = fpp_msg.txt

text_libs = fppbld.tlb/lib + csdr$library:futlib.tlb/lib + -
			     csdr$library:csdrlib.tlb/lib

.for.obj
  $(fort)/extend_source $(fflags) $(mms$source) + $(text_libs)

fpp_object =	fpp.obj,-
		fpp_parse_command.obj,-
		fpp_init_report.obj,-
		fpp_cat_info.obj,-
		fpp_open_archive.obj,-
		fpp_open_fex.obj,-
		fpp_check_time.obj,-
		fpp_midpoint_time.obj,-
		fpp_read2.obj,-
		fpp_collect_time.obj,-
		fpp_gain.obj,-
		fpp_fakeit.obj,-
		fpp_close_archive.obj,-
		fpp_check_segment_overlap.obj,-
		fpp_close_report.obj,-
                fpp_msg.obj

fpp.obj : 		      fppbld.tlb,-
		              csdr$library:futlib.tlb,-
			      ct$library:ctuser.inc,-
                              nfs_sdf^, -
                              nfs_anc^

fpp_parse_command.obj :       fppbld.tlb,-
		              csdr$library:futlib.tlb,-
			      ct$library:ctuser.inc

fpp_init_report.obj :         fppbld.tlb,-
		              csdr$library:futlib.tlb,-

fpp_cat_info.obj :            fppbld.tlb,-
		              csdr$library:futlib.tlb,-
			      ct$library:ctuser.inc,-
                              nfs_sdf^, -
                              nfs_anc^

fpp_open_archive.obj :        fppbld.tlb,-
			      ct$library:ctuser.inc,-
                              nfs_sdf^, -
                              nfs_anc^

fpp_open_fex.obj :            fppbld.tlb,-
		              csdr$library:futlib.tlb,-
			      ct$library:ctuser.inc,-
			      fex_fakeit^,-
			      fex_gain^

fpp_check_time.obj :          fppbld.tlb,-
			      ct$library:ctuser.inc,-
                              nfs_sdf^, -
                              nfs_anc^

fpp_midpoint_time.obj :       fppbld.tlb,-
			      fex_mtmsweep^

fpp_read2.obj :               fppbld.tlb,-
		              csdr$library:futlib.tlb,-
			      ct$library:ctuser.inc,-
                              nfs_anc^

fpp_collect_time.obj :        fppbld.tlb,-
		              csdr$library:futlib.tlb,-
			      ct$library:ctuser.inc,-
                              nfs_sdf^, -
                              nfs_anc^

fpp_gain.obj :                fppbld.tlb,-
			      ct$library:ctuser.inc,-
                              nfs_sdf^, -
			      fex_gain^

fpp_fakeit.obj :              fppbld.tlb,-
			      ct$library:ctuser.inc,-
                              nfs_sdf^, -
			      fex_fakeit^

fpp_close_archive.obj :       fppbld.tlb,-
			      ct$library:ctuser.inc

fpp_check_segment_overlap.obj  : fppbld.tlb,-
		                 csdr$library:futlib.tlb,-
				 ct$library:ctuser.inc

fpp_close_report.obj :        fppbld.tlb,-
		              csdr$library:futlib.tlb

fpp_msg.obj : fpp_msg.msg
  MESSAGE fpp_msg

fppbld.tlb : $(inc_files)
  LIBRARY/CREATE/TEXT fppbld $(inc_files)
  ! FPPBLD.TLB has been built.

fpp.exe :	fppbld.olb( $(fpp_object) ),-
		fppbld.tlb,-
		csdr$library:futlib.olb,-
		csdr$library:v5.opt
  $(link) $(linkflags) fppbld.olb/lib/inc=(fpp),-
                       csdr$library:futlib.olb/lib,-
	   	       csdr$library:csdrlib.olb/lib,-
	   	       csdr$library:csdrmsg.olb/inc=(upm_msg,cct_msg,cut_msg),-
	   	       csdr$library:v5.opt/option
  ! FPPBLD.OLB has been built.
