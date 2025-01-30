! Installs facility FUT.
!
! Author: Rob Kummerer
!	  October 14, 1986
!	  STI Incorporated


fut : csdr$library:futlib.olb, -
      csdr$library:futlib.tlb, -
      csdr$library:csdrmsg.olb(fut_msg),-
      csdr$help:csdrhelp.hlb(firas_utilities=fut.hlp)
  ! Facility FUT has been installed.

csdr$library:futlib.olb : futlib.olb
  COPY futlib.olb csdr$library:futlib.olb;0
  ! FUTLIB.OLB has been installed.

csdr$library:futlib.tlb : futlib.tlb
  COPY futlib.tlb csdr$library:futlib.tlb;0
  ! FUTLIB.TLB has been installed.
