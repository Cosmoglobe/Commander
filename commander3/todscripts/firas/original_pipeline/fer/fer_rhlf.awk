#-------------------------------------------------------------------------------
#
#	RHLF.AWK
#
#	This GAWK script is designed to edit the first line of FIRAS FISH
#	PEMx.RHSF calibration model solution ASCII files, where x is [1-14] and
#	represents one of the JCJ bolometer parameter model solutions.  GAWK
#	changes the scan mode from RHSF to RHLF, changes the calibration model
#	solution label from 93hybrid to 93hybridx, and converts the Unix-format
#	timetag to a VMS-format timetag.
#
#	Author:   Gene Eplee
#		  General Sciences Corp.
#		  513-7768
#		  9 December 1993
#
#-------------------------------------------------------------------------------
{
   if ($2 ~ /RHSF/) {
      sub(/RHSF/, "RHLF", $0)
      result = $3$1
      if (length($1) == 1) {
         sub(/[1-9]/, "!", $0)
         sub(/93hybrid /, result, $0)
      }
      else {
         sub(/(1[0-4])/, " !", $0)
         sub(/93hybrid  /, result, $0)
      }
      sub(/Tue Apr 20/, "20-APR-1993", $0)
      print $0
    }
    else {
      print $0
    }
}
