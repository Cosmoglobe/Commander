#-------------------------------------------------------------------------------
#
#	RLSS0.AWK
#
#	This GAWK script is designed to edit the first line of FIRAS FISH
#	PEM.RLSS calibration model solution ASCII file.  GAWK changes the
#	calibration model solution label from 93hybrid to 93hybrid0 and
#	converts the Unix-format timetag to a VMS-format timetag.
#
#	Author:   Gene Eplee
#		  General Sciences Corp.
#		  513-7768
#		  9 December 1993
#
#-------------------------------------------------------------------------------
{
   if ($2 ~ /RLSS/) {
      sub(/93hybrid /, "93hybrid0", $0)
      sub(/Wed May 12/, "12-MAY-1993", $0)
      print $0
    }
    else {
      print $0
    }
}
