! Builds facility FNT which consists of FNT_NOISETEST.
!
! Author: Rob Kummerer
!	  November 15, 1986
!	  STI Incorporated
!
!	  Adapted from work done by W. Young.
!
! Changes:
!	SPR 6013, R. Kummerer, February 12, 1990. VMS 5.2 options file changes.
!	SPR 6355, R. Kummerer, March 2, 1990. Prevent MMS abort on warning
!		"Multiply defined symbol WAIT".

fnt : fntbld.tlb, fnt.exe
  !Facility FNT is up-to-date.

! Define macros

text_libs = fntbld.tlb/lib + csdr$library:futlib.tlb/lib + -
			     csdr$library:csdrlib.tlb/lib
inc_files = fnt_invoc.txt

.for.obj
  $(fort)/extend_source $(fflags) $(mms$source) + $(text_libs)

fnt_obj = fnt.obj, -
	  fnt_parse_noise.obj, -
	  fnt_read_noise.obj, -
	  fnt_convert_to_volts.obj, -
	  fnt_power_density.obj, -
	  fnt_display_noise.obj, -
	  fnt_write_noise.obj, -
	  fnt_get_transfcn.obj, -
	  fnt_msg.obj

! Dependencies.

fnt.obj 		: fnt.for,-
			  fntbld.tlb(fnt_invoc),-
			  csdr$library:futlib.tlb,-
			  nfs_sdf^,-
			  fnt_noise^

fnt_parse_noise.obj	: fnt_parse_noise.for, -
			  fntbld.tlb(fnt_invoc),-
			  csdr$library:futlib.tlb, -
			  ct$library:ctuser.inc

fnt_read_noise.obj	: fnt_read_noise.for, -
			  csdr$library:futlib.tlb, -
			  ct$library:ctuser.inc,-
			  nfs_sdf^,-
			  nfs_hkp^

fnt_convert_to_volts.obj : fnt_convert_to_volts.for, -
			   fntbld.tlb(fnt_invoc),-
			   csdr$library:futlib.tlb,-
			   fnt_noise^

fnt_power_density.obj	: fnt_power_density.for, -
			  csdr$library:futlib.tlb,-
			  nfs_sdf^,-
			  fnt_noise^

fnt_display_noise.obj	: fnt_display_noise.for, -
			  csdr$library:futlib.tlb,-
			  fnt_noise^

fnt_write_noise.obj	: fnt_write_noise.for, -
			  csdr$library:futlib.tlb, -
			  ct$library:ctuser.inc,-
			  fnt_noise^

fnt_get_transfcn.obj	: fnt_get_transfcn.for, -
			  csdr$library:futlib.tlb

fnt_msg.obj : fnt_msg.msg
  MESSAGE fnt_msg

fntbld.tlb : $(inc_files)
  LIBRARY/CREATE/TEXT fntbld $(inc_files)
  ! FNTBLD.TLB has been built.

fnt.exe : fntbld.olb( $(fnt_obj) ),-
	  csdr$library:futlib.olb,-
	  imsl$dir:imsl.olb,-
	  csdr$library:v5.opt
 - $(link) $(linkflags) fntbld.olb/lib/inc=(fnt),-
	  		imsl$dir:imsl.olb/lib,-
			csdr$library:futlib.olb/lib,-
			csdr$library:csdrlib.olb/lib,-
		        xanadu:[lib]vilib/lib, -
		        xanadu:[lib]xlib/lib, -
		        graphics:grpshr/lib, -
			csdr$library:v5.opt/option
  !
  ! The "Multiply defined symbol WAIT" warning may be ignored.
  !
  ! FNT.EXE has been built.
