#-------------------------------------------------------------------------------
#
#	LLSF0.AWK
#
#	This GAWK script is designed to edit the first line of FIRAS FISH
#	PEM.LLSF calibration model solution ASCII file.  GAWK changes the
#	calibration model solution label from 93hybrid to 93hybrid0 and
#	converts the Unix-format timetag to a VMS-format timetag.
#
#	Author:   Gene Eplee
#		  General Sciences Corp.
#		  513-7768
#		  18 April 1994
#
#-------------------------------------------------------------------------------
{
   if ($2 ~ /LLSF/) {
      sub(/93hybrid /, "93hybrid0", $0)
      sub(/Fri Apr 15/, "15-APR-1994", $0)
      print $0
    }
    else {
      print $0
    }
}
