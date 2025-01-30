! Builds facility FEC_EXTRACT_CALIBRATION
!
! Author: Reid Wilson
!         STI,Inc.
!
! Changes:
!	SPR 6007, R. Kummerer, VMS 5.2 options file changes.
!       SER 6851, H. Wang, STX, FEC improvement.
!                 June 19, 1990.          
FEC : fec.exe
 !Facility FEC is up to date

! define macros
TXTLIBS = FECBLD.TLB/LIB + CSDR$LIBRARY:FUTLIB.TLB/LIB+CSDR$LIBRARY:CSDRLIB.TLB/LIB

.for.obj
  $(FORT) /extend_source $(FFLAGS) $(MMS$SOURCE) + $(TXTLIBS)

! macros for fec
fec_inc_files = -
 fec_inc.txt, -
 fec_msg.txt

fec_obj = -
 FEC.obj, -
 FEC_GET_QUALS.obj,-
 FEC_LOAD.obj,-
 FEC_MSG.obj, -
 FEC_ANALYZE.obj,-
 FEC_COMMAND_CHANGE.obj,-
 FEC_HOT_SPOT_COMMAND_CHANGE.obj,-
 FEC_HOT_FLAG_PLATEAU.obj,-
 FEC_OPEN.obj,-
 FEC_CLOSE.obj,-
 FEC_SORT.obj,-
 FEC_WRITE_RECORDS.obj,-
 FEC_REPORT_PLATEAU.obj
!
!   SYMBOLS FOR INCLUDE FILES

!
!       ** COMPILE FILE FOR PROGRAM FEC_EXTRACT_CALIBRATION_INDEX **
!
fec_load.obj : FEC_LOAD.FOR -
  CSDR$LIBRARY:FUTLIB.TLB FECBLD.TLB($(FEC_INC_FILES))
fec_msg.obj : fec_msg.msg
  MESSAGE/NOLIST/OBJECT=FEC_MSG.OBJ FEC_MSG.MSG
fec_analyze.obj : FEC_ANALYZE.FOR -
  CSDR$LIBRARY:FUTLIB.TLB FECBLD.TLB($(FEC_INC_FILES)) 
fec_open.obj : FEC_OPEN.FOR -
  CSDR$LIBRARY:FUTLIB.TLB FECBLD.TLB($(FEC_INC_FILES)) 
fec_command_change.obj : FEC_COMMAND_CHANGE.FOR -
  CSDR$LIBRARY:FUTLIB.TLB FECBLD.TLB($(FEC_INC_FILES)) 
fec_hot_spot_command_change.obj : FEC_HOT_SPOT_COMMAND_CHANGE.FOR -
  CSDR$LIBRARY:FUTLIB.TLB FECBLD.TLB($(FEC_INC_FILES)) 
fec_hot_flag_plateau.obj : FEC_HOT_FLAG_PLATEAU.FOR -
  CSDR$LIBRARY:FUTLIB.TLB FECBLD.TLB($(FEC_INC_FILES)) 
fec_write_records.obj : FEC_WRITE_RECORDS.FOR -
  CSDR$LIBRARY:FUTLIB.TLB FECBLD.TLB($(FEC_INC_FILES)) -
  FEC_SSCAL^
fec_report_plateau.obj : FEC_REPORT_PLATEAU.FOR -
  CSDR$LIBRARY:FUTLIB.TLB FECBLD.TLB($(FEC_INC_FILES)) 
fec.obj : FEC.FOR    CSDR$LIBRARY:FUTLIB.TLB FECBLD.TLB($(FEC_INC_FILES)) 
fec_get_quals.obj : FEC_GET_QUALS.FOR   CSDR$LIBRARY:FUTLIB.TLB FECBLD.TLB($(FEC_INC_FILES)) 
fec_close.obj : FEC_CLOSE.FOR   CSDR$LIBRARY:FUTLIB.TLB FECBLD.TLB($(FEC_INC_FILES)) $(CTLIB) 
fec_sort.obj : FEC_SORT.FOR   CSDR$LIBRARY:FUTLIB.TLB FECBLD.TLB($(FEC_INC_FILES)) 

! remember all target statements must begin in column one, all action 
! statements must begin in any column > one.



fec.exe : fecbld.olb($(fec_obj)) -
          csdr$library:futlib.olb -
          csdr$library:v5.opt
      $(link) $(LINKFLAGS) fecbld.olb/lib/include=fec, -
                           CSDR$LIBRARY:futlib/L/include=fut_msg, -
                           csdr$library:csdrlib/L, -
                           csdr$library:v5.opt/option
