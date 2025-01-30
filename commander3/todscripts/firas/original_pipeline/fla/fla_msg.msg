!-------------------------------------------------------------------------------
!
!	FLA_MSG.MSG
!
!	Message definition file for FIRAS facility FLA_LINES_ANALYSIS.
!
!	Author:   Gene Eplee
!		  General Sciences Corp.
!		  513-7768
!		  1 July 1993
!
!       Modifications:   S. Brodd, HSTX, 7/31/96, SPR 12333.
!                        Fix facility number.
!-------------------------------------------------------------------------------

.title		FLA_MSG Message Definition File
.ident		'Version 13.9'
.facility	FLA,722/Prefix=FLA_

.base 0
.severity success

normal		<Process FLA succeeded.>

.base 2000
.severity error

cfgdirclose	<Error closing direct config files. Status !SL>/fao=1
cfgdirget	<Error getting direct config files. Status !SL>/fao=1
cfgdiropen	<Error opening direct config files. Status !SL>/fao=1
cfgseqclose	<Error closing sequential config files. Status !SL>/fao=1
cfgseqget	<Error getting sequential config files. Status !SL>/fao=1
cfgseqopen	<Error opening sequential config files. Status !SL>/fao=1

ctinit		<Error from Cobetrieve initialization.  Status !SW>/fao=1

boltemp		<No bolometer temperature specified for /TEMP qualifier.>
bolvolt		<No bolometer voltage specified for /VOLT qualifier.>
cmdbias		<No commanded bias specified for /BIAS qualifier.>
failure		<Process FLA failed.>
nomodel		<No model solution specified for /MODEL_EXT qualifier.>
noreftime	<Neither /FLIGHT nor /INT reference time qualifier set.>
readetf		<Error reading electronics transfer function !SL.  Status !SL>/fao=1

repclose	<Error closing processing report !AS.  Status !SL>/fao=2
repopen		<Error opening processing report !AS.  Status !SL>/fao=2
repwrite	<Error writing to processing report !AS.  Status !SL>/fao=2

rmsclose	<Error closing RMS file !AS.  Status !SL>/fao=2
rmsopen		<Error opening RMS file !AS.  Status !SL>/fao=2
rmsread		<Error reading from RMS file !AS.  Status !SL>/fao=2
rmswrite	<Error writing to RMS file !AS.  Status !SL>/fao=2

.end
