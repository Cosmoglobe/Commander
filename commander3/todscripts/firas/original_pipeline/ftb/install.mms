! Installs facility FTB.
!
! Author: Rob Kummerer
!	  December 19, 1986
!	  STX
!
! Changes:
!	Add FTB_MAKE_RSE, FTB_SNAPSHOT, FTB_MSG references. R. Kummerer,
!	Feb 18, 1987.
!	Remove GETRSE reference. R. Kummerer, Mar 25, 1987.
!	Add FTB_TESTLOG reference. R. Kummerer, Apr 6, 1987.
!       Add JDATE reference, L. Olson, 5/21/87
!	Add FTB_WORD31 reference. R. Kummerer, Aug 31, 1988.
!	Add FTB_TEMPCONTROL reference. R. Kummerer, May 23, 1989.
!	Remove FTB_MAKE_RSE. SPR 3161. R. Kummerer, July 18, 1989.
!	Remove FTB_SNAPSHOT. SPR 3493. R. Kummerer, August 25, 1989.
!	Add FTB_GAPFLAG. SER 4023. R. Kummerer, September 6, 1989.
!	Add FTB_OVERLAP. SER 4210. R. Kummerer, November 6, 1989.
!	Add FTB_BINFLD.  SER 6852. N. Gonzales, June 1, 1990.
!	Add FTB_Checksum.  SER 6852. N. Gonzales, June 1, 1990.
!	Add FTB_VALID_TIME. SER 6852. N. Gonzales, June 1, 1990.
!       Add FTB_CHECK_DQ. SER 6852. N. Gonzales, June 14, 1990.
!       Remove FTB_GAPFLAG. SPR 9726. S. Alexander, June 29, 1992.
!	Add FTB_TRANS. SER 10764, 10772. S. Alexander, April 5, 1993.
!	Add FTB_LTRANS. SER 12285. S. Brodd, December 20, 1995.
!
.suffixes
.suffixes .olb .obj .msg .tlb .cld .exe .hlb .hlp 

ftb : csdr$cld:ftb.cld, -
      csdr$library:csdrmsg.olb(ftb_msg), -
      csdr$help:csdrhelp.hlb(ftb.hlp), -
      csdr$system:ftb_gettim.exe, -
      csdr$system:ftb_ctviewer.exe, -
      csdr$system:ftb_convtime.exe, -
      csdr$system:ftb_jdate.exe, -
      csdr$system:ftb_mtm_jitter.exe, -
      csdr$system:ftb_word31.exe, -
      csdr$system:ftb_testlog.exe, -
      csdr$system:ftb_histo.exe, -
      csdr$system:ftb_tempcontrol.exe, -
      csdr$system:ftb_overlap.exe, -
      csdr$system:ftb_binfld.exe, -
      csdr$system:ftb_checksum.exe, -
      csdr$system:ftb_valid_time.exe, -
      csdr$system:ftb_check_dq.exe, -
      csdr$system:ftb_trans.exe, -
      csdr$system:ftb_ltrans.exe
  ! Facility FTB has been installed.

csdr$cld:ftb.cld : ftb.cld
  COPY ftb.cld csdr$cld:ftb.cld;0
  ! FTB.CLD has been installed.

csdr$system:ftb_gettim.exe : ftb_gettim.exe
  COPY ftb_gettim.exe csdr$system:ftb_gettim.exe;0
  ! FTB_GETTIM.EXE has been installed.

csdr$system:ftb_jdate.exe : ftb_jdate.exe
  COPY ftb_jdate.exe csdr$system:ftb_jdate.exe;0
  ! FTB_JDATE.EXE has been installed.

csdr$system:ftb_ctviewer.exe : ftb_ctviewer.exe
  COPY ftb_ctviewer.exe csdr$system:ftb_ctviewer.exe;0
  ! FTB_CTVIEWER.EXE has been installed.

csdr$system:ftb_convtime.exe : ftb_convtime.exe
  COPY ftb_convtime.exe csdr$system:ftb_convtime.exe;0
  ! FTB_CONVTIME.EXE has been installed.

csdr$system:ftb_mtm_jitter.exe : ftb_mtm_jitter.exe
  COPY ftb_mtm_jitter.exe csdr$system:ftb_mtm_jitter.exe;0
  ! FTB_MTM_JITTER.EXE has been installed.

csdr$system:ftb_word31.exe : ftb_word31.exe
  COPY ftb_word31.exe csdr$system:ftb_word31.exe;0
  ! FTB_WORD31.EXE has been installed.

csdr$system:ftb_testlog.exe : ftb_testlog.exe
  COPY ftb_testlog.exe csdr$system:ftb_testlog.exe;0
  ! FTB_TESTLOG.EXE has been installed.

csdr$system:ftb_histo.exe : ftb_histo.exe
  COPY ftb_histo.exe csdr$system:ftb_histo.exe;0
  ! FTB_HISTO.EXE has been installed.

csdr$system:ftb_tempcontrol.exe : ftb_tempcontrol.exe
  COPY ftb_tempcontrol.exe csdr$system:ftb_tempcontrol.exe;0
  ! FTB_TEMPCONTROL.EXE has been installed.

csdr$system:ftb_overlap.exe : ftb_overlap.exe
  COPY ftb_overlap.exe csdr$system:ftb_overlap.exe;0
  ! FTB_OVERLAP.EXE has been installed.

csdr$system:ftb_binfld.exe : ftb_binfld.exe
  COPY ftb_binfld.exe csdr$system:ftb_binfld.exe;0
  ! FTB_BINFLD.EXE has been installed.

csdr$system:ftb_checksum.exe : ftb_checksum.exe
  COPY ftb_checksum.exe csdr$system:ftb_checksum.exe;0
  ! FTB_CHECKSUM.EXE has been installed.

csdr$system:ftb_valid_time.exe : ftb_valid_time.exe
  COPY ftb_valid_time.exe csdr$system:ftb_valid_time.exe;0
  ! FTB_VALID_TIME.EXE has been installed.

csdr$system:ftb_check_dq.exe : ftb_check_dq.exe
  COPY ftb_check_dq.exe csdr$system:ftb_check_dq.exe;0
  ! FTB_CHECK_DQ.EXE has been installed.

csdr$system:ftb_trans.exe : ftb_trans.exe
  COPY ftb_trans.exe csdr$system:ftb_trans.exe;0
  ! FTB_TRANS.EXE has been installed.

csdr$system:ftb_ltrans.exe : ftb_ltrans.exe
  COPY ftb_ltrans.exe csdr$system:ftb_ltrans.exe;0
  ! FTB_LTRANS.EXE has been installed.
