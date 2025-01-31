!==============================================================================
!
!   FIRAS Regression Testing Facility (FRT) Installation
!
!
!  MMS Revision History:
!
!    Date   	Version   SPR#   Programmer   Comments
!    ----   	-------   ----   ----------   --------
!
!  10/08/91       8.8            A.C.Raugh    Creation
!
!==============================================================================

!  Clear suffix list:

.suffixes


! **********************  FRT Facility Installation  **************************

FRT : csdr$idl:frt_idl.tlb
   ! ***                                     ***
   ! ***  FRT_IDL.TLB installation complete  ***
   ! ***                                     ***

CSDR$IDL:FRT_IDL.TLB : frt_idl.tlb
   copy frt_idl.tlb csdr$idl:frt_idl.tlb;0
