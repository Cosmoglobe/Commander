!********************************************************************
! File DESCRIP.MMS for FSS facility
!
! Author: D. Bouler, STX, CDAC, Jan, 1991
!         Modified from FES descrip file by Reid Wilson
!
!********************************************************************
! Changes:
!    L. Rosen, HSTX, Feb 1995, SER 12245
!    Modified for new routines to include records from neighbor pixels.
!
!********************************************************************

FSS : fssbld.tlb, fss.exe
 !Facility FSS is up to date

.for.obj
  $(FORT) /extend_source $(FFLAGS) $(MMS$SOURCE) + $(TXTLIBS)

txtlibs=fssbld.tlb/lib+csdr$library:futlib.tlb/lib+csdr$library:csdrlib.tlb/lib

fss_inc_files = fss_include.txt

fssbld.tlb : $(fss_inc_files)
  lib/create/text fssbld.tlb $(fss_inc_files)
  !FSS text library has been built.

fss_obj =  fss.obj,                -
           fss_close_skymaps.obj,  -
           fss_form_groups.obj,    -
           fss_get_quals.obj,      -
           fss_get_rse.obj,        -
           fss_init_report,        -
           fss_sort.obj,           -
           fss_msg.obj,            -
           fss_open_skymaps.obj,   -
           fss_read_config.obj,    -
           fss_read_eng.obj,       -
           fss_read_science.obj,   -
           fss_read_skymap.obj,    -
           fss_write_report.obj,   -
           fss_write_skymap.obj,   -
           fss_get_sort_neighbors.obj, -
           fss_match_group.obj,    -
           fss_match_pixel.obj,    -
           fss_avg_gallat.obj

fss_msg.obj : fss_msg.msg
 MESSAGE/OBJECT=$(MMS$TARGET) $(MMS$SOURCE)

fss.exe : fssbld.olb($(fss_obj)), -
          csdr$library:csdrlib.olb, -
          csdr$library:futlib.olb, -
          csdr$library:csdrmsg.olb, -
          imsl$dir:imsl.olb, -
	  csdr$library:v5.opt
 $(link) $(linkflags) fssbld.olb/lib/include=(fss,fss_msg), -
                      csdr$library:csdrlib/lib, -
                      csdr$library:futlib/lib, -
                      csdr$library:csdrmsg/lib, -
                      imsl$dir:imsl/lib, -
                      csdr$library:v5.opt/option

fss.obj                 : fssbld.tlb(fss_include),  fss_sssky^, -
                          fex_grtrawwt^,            fex_grttrans^
fss_form_groups.obj     : fssbld.tlb(fss_include),  fss_sssky^
fss_sort.obj            : fssbld.tlb(fss_include),  fss_sssky^
fss_read_config.obj     : fssbld.tlb(fss_include),  fex_mincoadd^,   -
                          fex_grtrawwt^,            fex_grttrans^
fss_read_eng.obj        : fssbld.tlb(fss_include),  fdq_eng^,   -
                          fex_grtrawwt^,            fex_grttrans^
fss_read_science.obj    : fssbld.tlb(fss_include),  fdq_sdf^,   fss_sssky^
fss_read_skymap.obj     : fssbld.tlb(fss_include),  fss_sssky^
fss_write_skymap.obj    : fssbld.tlb(fss_include),  fss_sssky^
fss_get_sort_neighbors.obj   : fssbld.tlb(fss_include),  fss_sssky^
fss_match_group.obj     : fssbld.tlb(fss_include),  fss_sssky^
fss_match_pixel.obj     : fssbld.tlb(fss_include),  fss_sssky^
fss_avg_gallat.obj      : fssbld.tlb(fss_include),  fss_sssky^
