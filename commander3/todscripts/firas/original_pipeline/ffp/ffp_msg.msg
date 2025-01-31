!-------------------------------------------------------------------------------
!
!	FFP_MSG.MSG
!
!	Message definition file for FIRAS facility FFP_FINAL_PRODUCT.
!
!	Author:  S. Brodd, HSTX, 3/21/96
!
!-------------------------------------------------------------------------------

.title		FFP_MSG Message Definition File
.ident		'Version 13.8'
.facility	FFP,704/Prefix=FFP_

.base 0
.severity success

normal		<Process FFP succeeded.>

.base 1000
.severity informational

eof		<End of pixel numbers encountered.>
abort		<Process aborted.>

.base 2000
.severity error

cfgdirclose	<Error closing direct config files. Status !SL>/fao=1
cfgdirget	<Error getting direct config files. Status !SL>/fao=1
cfgdiropen	<Error opening direct config files. Status !SL>/fao=1
cfgseqclose	<Error closing sequential config files. Status !SL>/fao=1
cfgseqget	<Error getting sequential config files. Status !SL>/fao=1
cfgseqopen	<Error opening sequential config files. Status !SL>/fao=1

csaclose	<Error closing skymap !AS.  Status !SL>/fao=2
csafld		<Error setting pixel offsets for skymap !AS.  Status !SL>/fao=2
csaopen		<Error opening skymap !AS.  Status !SL>/fao=2
csaread		<Error reading from skymap !AS.  Status !SL>/fao=2
csawrite	<Error writing to skymap !AS.  Status !SL>/fao=2

ctinit		<Error from Cobetrieve initialization.  Status !SW>/fao=1

adttime		<ADT time conversion error for spectrum !AS, pixel !SL.>/fao=2
boltemp		<No bolometer temperature specified for /TEMP qualifier.>
bolvolt		<No bolometer voltage specified for /VOLT qualifier.>
closerr		<Error closing file !AS; status !SL.>/fao=2
cmdbias		<No commanded bias specified for /BIAS qualifier.>
failure		<Process FFP failed.>
freqerr		<Low frequency cut !AS should be less than the high !AS.>/fao=2
galexc          <Invalid galactic latitude exclusion angle. Status !SL>/fao=1
gettimerr	<Error from Sys$GetTim; status = !SL.>/fao=1
invalfile       <Invalid input file specified for /FILE qualifier:  !AS.>/fao=1
lunerr		<Error getting logic unit number; status !SL.>/fao=1
maxrec		<Spectrum buffer overflow.  !SL spectra read for pixel !SL.>/fao=2
nocvs           <No C vector file specified for /FILE_EXT qualifier.>
nodust          <No dust spectrum file specified for /FILE_EXT qualifier.>
noejv           <No correction file specified for /FILE_EXT qualifier.>
noerr           <No error file specified for /FILE_EXT qualifier.>
nofef           <No FEF file specified for /FILE_EXT qualifier.>
nofile		<No input filename specified for /FILE qualifier.>
nofilext	<No output file extension specified in script file.>
noinput		<No input files specified in script file.>
nolines         <No lines file specified for /LINES_EXT qualifier.>
nomodel		<No model solution specified for /MODEL_EXT qualifier.>
noreftime	<Neither /FLIGHT nor /INT reference time qualifier set.>
numfreqerr	<Frequency range too large; number of frequency bins: !SL.>/fao=1
openerr		<Error opening file !AS; status !SL.>/fao=2
parserr		<Error command line parsing of: !AS>/fao=1
readapod	<Error reading apodization function !SL.  Status !SL>/fao=2
readerr		<Error reading file !AS; status !SL.>/fao=2
readetf		<Error reading electronics transfer function !SL.  Status !SL>/fao=1

repclose	<Error closing processing report !AS.  Status !SL>/fao=2
repopen		<Error opening processing report !AS.  Status !SL>/fao=2
repwrite	<Error writing to processing report !AS.  Status !SL>/fao=2

rmsclose	<Error closing RMS file !AS.  Status !SL>/fao=2
rmsopen		<Error opening RMS file !AS.  Status !SL>/fao=2
rmsread		<Error reading from RMS file !AS.  Status !SL>/fao=2
rmswrite	<Error writing to RMS file !AS.  Status !SL>/fao=2
writerr		<Error writing file !AS; status !SL.>/fao=2
ctnocatrecs	<No catalog records found in archive !AS for dataset !AS.>/fao=2
ctquerycat	<Error querying archive !AS for file !AS.  Status !SL>/fao=3

.end
