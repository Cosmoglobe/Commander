!-------------------------------------------------------------------------------
!
!	FSL_MSG.MSG
!
!	Message definition file for FIRAS facility FSL_Spectra_Long
!
!	Author:
!                 FCF_Msg.Msg   
!                 Gene Eplee
!		  General Sciences Corporation
!		  11 May 1993
!
!                 FSL_Msg.Msg
!                 Shirley M. Read
!                 Hughes STX Corporation
!                 7 July 1995
!
!-------------------------------------------------------------------------------
!
!	Changes:
!
!	Added FEX_VAR error messages.
!	Gene Eplee, 25 October 1993
!       SER 11397
!
!       Modified FCF_MSG.MSG to FSL_MSG.MSG for new pipeline to process 
!       long spectra. Added command line input error messages for required 
!       qualifiers which were not present. Modified facility number for FSL.
!       Added messages for I/O for FEX_FLV reference file and invalid data
!       in the record.
!       Shirley M. Read, 7 July 1995
!
!-------------------------------------------------------------------------------

.title		FSL_MSG Message Definition File
.ident		'Version 13.9'
.facility	FSL,712/Prefix=FSL_

.base 0
.severity success

normal		<Process FSL succeeded.>

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

csafldoffset	<Error setting pixel offsets for skymap !AS.  Status !SL>/fao=2
csaifgclose	<Error closing coadd skymap !AS.  Status !SL>/fao=2
csaifgopen	<Error opening coadd skymap !AS.  Status !SL>/fao=2
csaifgread	<Error reading from coadd skymap !AS.  Status !SL>/fao=2
csapixnoinval	<Invalid pixel number !SL read from coadd skymap !AS.>/fao=2
csaspecclose	<Error closing spectrum skymap !AS.  Status !SL>/fao=2
csaspecnum	<Only !SL of !SL spectra written to the skymap !AS.>/fao=3
csaspecopen	<Error opening spectrum skymap !AS.  Status !SL>/fao=2
csaspecwrite	<Error writing to spectrum skymap !AS.  Status !SL>/fao=2

ctifgclose	<Error closing coadd file !AS.  Status !SL>/fao=2
ctifgopen	<Error opening coadd file !AS.  Status !SL>/fao=2
ctifgread	<Error reading from coadd file !AS.  Status !SL>/fao=2
ctspecclose	<Error closing spectrum file !AS.  Status !SL>/fao=2
ctspecnum	<Only !SL of !SL spectra written to the file !AS.>/fao=3
ctspecopen	<Error opening spectrum file !AS.  Status !SL>/fao=2
ctspecwrite	<Error writing to spectrum file !AS.  Status !SL>/fao=2

ctgetmodext	<Error reading file extension for model !AS.  Status !SL>/fao=2
ctinit		<Error from Cobetrieve initialization.  Status !SW>/fao=1
ctnocatrecs	<No catalog records found in archive !AS for dataset !AS.>/fao=2
ctquerycat	<Error querying archive !AS for file !AS.  Status !SL>/fao=3

dettempconv	<Bolometer temperature did not converge for spectrum !AS, pixel !SL.>/fao=2
failure		<Process FSL failed.>

flvclose	<Error closing FIL variances file !AS.  Status !SL>/fao=2
flvinval	<Incorrect variance in FIL variances file !AS>/fao=1
flvopen		<Error opening FIL variances file !AS.  Status !SL>/fao=2
flvread		<Error reading from FIL variances file !AS.  Status !SL>/fao=2

ftemplist	<FUT_TEMP_LIST failure for temperature !SL, spectrum !AS, pixel !SL.>/fao=3
invalbaryvel	<Invalid barycentric velocity for spectrum !AS, pixel !SL.>/fao=2
invaltime	<Jstop time !AS is less than Jstart time !AS.>/fao=2
maxcoaddexc	<Coadd buffer overflow.  !SL coadds read for pixel !SL.>/fao=2
nochan          <No channel specified for /CHANNEL qualifier.>
noinput         <No input data type (SKY or CAL) specified for /INPUT qualifier.>
nomodfile	<No model solution specified for /MODEL_EXT qualifier.>
nopixlist	<No pixels specified for /PIXEL qualifier.>
noscan          <No scan mode specified for /SCAN_MODE qualifier.>
modelclose	<Error closing calibration model solution file !AS.  Status !SL>/fao=2
modelinval	<Incorrect model solution in model file !AS>/fao=1
modelopen	<Error opening calibration model solution file !AS.  Status !SL>/fao=2
modelread	<Error reading from model calibration solution file !AS.  Status !SL>/fao=2

readapod	<Error reading apodization function !SL.  Status !SL>/fao=2
readetf		<Error reading electronics transfer function !SL.  Status !SL>/fao=1

repclose	<Error closing processing report !AS.  Status !SL>/fao=2
repopen		<Error opening processing report !AS.  Status !SL>/fao=2
repwrite	<Error writing to processing report !AS.  Status !SL>/fao=2

varclose	<Error closing variance file !AS.  Status !SL>/fao=2
varinval	<Incorrect variance in variance file !AS>/fao=1
varopen		<Error opening variance file !AS.  Status !SL>/fao=2
varread		<Error reading from variance file !AS.  Status !SL>/fao=2

.end
