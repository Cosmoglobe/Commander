!-------------------------------------------------------------------------------
!
!	FFI_MSG.MSG
!
!	Message definition file for FIRAS facility FFI_Fish_Input.
!
!	Author:   Gene Eplee
!		  General Sciences Corp.
!		  513-7768
!		  9 March 1993
!		  SER 10763
!
!-------------------------------------------------------------------------------

.title		FFI_MSG Message Definition File
.ident		'Version 10.6'
.facility	FFI,716/Prefix=FFI_

.base 0
.severity success

normal		<Process FFI succeeded.>

.base 1000
.severity informational

chanerr		<No good coadds records read for this channel and scan mode.>
eof		<End-of-file encountered while reading input.>

.base 2000
.severity error

cfgdirclose	<Error closing direct config files. Status !SL>/fao=1
cfgdirget	<Error getting direct config files. Status !SL>/fao=1
cfgdiropen	<Error opening direct config files. Status !SL>/fao=1
cfgseqclose	<Error closing sequential config files. Status !SL>/fao=1
cfgseqget	<Error getting sequential config files. Status !SL>/fao=1
cfgseqopen	<Error opening sequential config files. Status !SL>/fao=1

ctifgclose	<Error closing coadd file !AS.  Status !SL>/fao=2
ctifgopen	<Error opening coadd file !AS.  Status !SL>/fao=2
ctifgread	<Error reading from coadd file !AS.  Status !SL>/fao=2

ctinit		<Error from Cobetrieve initialization.  Status !SW>/fao=1

failure		<Process FFI failed.>
ftemplist	<FUT_TEMP_LIST failure for temperature !SL, spectrum !AS, pixel !SL.>/fao=3
invaltime	<Jstop time !AS is less than Jstart time !AS.>/fao=2
nooutfile	<No output file name extension specified.>

readapod	<Error reading apodization function !SL.  Status !SL>/fao=2
readetf		<Error reading electronics transfer function !SL.  Status !SL>/fao=1

repclose	<Error closing processing report !AS.  Status !SL>/fao=2
repopen		<Error opening processing report !AS.  Status !SL>/fao=2
repwrite	<Error writing to processing report !AS.  Status !SL>/fao=2

rmsclose	<Error closing RMS spectrum file !AS.  Status !SL>/fao=2
rmsopen		<Error opening RMS spectrum file !AS.  Status !SL>/fao=2
rmswrite	<Error writing to RMS spectrum file !AS.  Status !SL>/fao=2

.end
