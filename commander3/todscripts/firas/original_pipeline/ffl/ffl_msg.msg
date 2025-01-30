!-------------------------------------------------------------------------------
!
!	FFL_MSG.MSG
!
!	Message definition file for FIRAS facility FFL_Fishinput_Long.
!
!	Author:   Fred Shuman, HSTX, 1995 Jul 11
!	Derived from FFI_MSG.MSG by Gene Eplee, GSC, 1993 March 09, SER 10763
!
!-------------------------------------------------------------------------------

.title		FFL_MSG Message Definition File
.ident		'Version 12.0'
.facility	FFL,714/Prefix=FFL_

.base 0
.severity success

normal		<Process FFL succeeded.>

.base 1000
.severity informational

chanerr		<No good coadd records read for this channel and scan mode.>
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

etfclose	<Error closing electronics transfer function !AS.  Status !SL>/fao=1
etfopen		<Error opening electronics transfer function !AS.  Status !SL>/fao=1
etfread		<Error reading electronics transfer function !AS.  Status !SL>/fao=1

failure		<Process FFL failed.>
fapodrecnuml	<FUT_Apod_RecNumL failure; return status=!SL.>/fao=1
apodread	<Error reading apodization function !SL.  Status !SL>/fao=2
fapodrotl	<FUT_Apod_RotL failure; return status=!SL.>/fao=1
ftemplistgrt	<FUT_TEMP_LIST failure for temperature !SL, coadd GMT !AS.>/fao=2
ftemplistcall	<FUT_TEMP_LIST was called improperly; return status=!SL.>/fao=1
invaltime	<Jstop time !AS is less than Jstart time !AS.>/fao=2
nooutfile	<No output file name extension specified.>

repclose	<Error closing processing report !AS.  Status !SL>/fao=2
repopen		<Error opening processing report !AS.  Status !SL>/fao=2
repwrite	<Error writing to processing report !AS.  Status !SL>/fao=2

rmsclose	<Error closing RMS spectrum file !AS.  Status !SL>/fao=2
rmsopen		<Error opening RMS spectrum file !AS.  Status !SL>/fao=2
rmswrite	<Error writing to RMS spectrum file !AS.  Status !SL>/fao=2

invscnmode	<Invalid scan mode--FS and FL are valid only for low channels.>

.end
