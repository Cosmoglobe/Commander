#-------------------------------------------------------------------------------
#
#	RLSF0.AWK
#
#	This GAWK script is designed to edit the first line of FIRAS FISH
#	PEM.RLSF calibration model solution ASCII file.  GAWK changes the
#	calibration model solution label from 93hybrid to 93hybrid0 and
#	converts the Unix-format timetag to a VMS-format timetag.
#
#	Author:   Gene Eplee
#		  General Sciences Corp.
#		  513-7768
#		  4 April 1994
#
#-------------------------------------------------------------------------------
{
   if ($2 ~ /RLSF/) {
      sub(/93hybrid /, "93hybrid0", $0)
      sub(/Tue Apr 20/, "20-APR-1993", $0)
      print $0
    }
    else {
      print $0
    }
}
