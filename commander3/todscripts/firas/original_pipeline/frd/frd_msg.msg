.Title       FRD_Msg
.Facility    FRD,710/Prefix=FRD_
.Base 0
.Severity    Success
NORMAL       <Successful completion.>
.Base 1000
.Severity    Informational
CTHKPEOF     <End of Housekeeping data found.>
eof	     <End of pixel numbers encountered.>
.Base 2000
.Severity    Warning
ABERR        <Process aborts due to non-recoverable error.>
CTINIT       <Error initializing cobetrieve. Status=!SL>/FAO=1
CTREADERR    <Error reading archive record. Status=!SL, Data set=!AS>/FAO=2
CTWRITERR    <Error writing archive record. Status=!SL, Data set=!AS>/FAO=2
CTOPENERR    <Error opening archive. Status=!SL, Data set=!AS>/FAO=2
CTCLOSERR    <Error closing archive. Status=!SL, Data set=!AS>/FAO=2
GETTIMERR    <Error return from SYS$GETTIM. Status=!SL>/FAO=1
GETLUNERR    <Error return from SYS$Get_Lun. Status=!SL>/FAO=1
RMSOPEN      <Error opening file !AS. Status=!SL>/FAO=2
RMSREAD      <Error writing file !AS. Status=!SL>/FAO=2
RMSWRITE     <Error writing file !AS. Status=!SL>/FAO=2
RMSCLOSE     <Error closing file !AS. Status=!SL>/FAO=2
GETLUERR     <Failed to get logical unit. status=!sl>/FAO=1
FILEMISS     <FEX TEXT FILE is not in the directory.>
BADTYPE      <Neither TOLS nor FLAG selected.>
ITEMNOFND    <Item field name not found from the TOLS' name table.>
RADBADLON    <Bad longitude in file !AS. Lonstep,Lonprev,Lon=!SL,!SL,!SL>/FAO=4
ABORT        <FEX failed on this status. Retstat=!sl>/FAO=1
closeconfig	<Error closing sequential config files. Status !SL>/fao=1
failure		<Process FRD failed.>
getconfig	<Error getting sequential config files. Status !SL>/fao=1
openconfig	<Error opening sequential config files. Status !SL>/fao=1
noinfile	<No input file specified for /FILE_EXT qualifier.>
rmsrew		<Error rewinding RMS file !AS.  Status !SL>/fao=2
rmssize		<Error sizing RMS file !AS.  Status !SL>/fao=2
csaclose	<Error closing skymap !AS.  Status !SL>/fao=2
csaopen	        <Error opening skymap !AS.  Status !SL>/fao=2
csaread	        <Error reading from skymap !AS.  Status !SL>/fao=2
csawrite	<Error writing to spectrum skymap !AS.  Status !SL>/fao=2
galexc          <Invalid galactic latitude exclusion angle. Status !SL>/fao=1
nofilext        <No file extension is specified in the script file !AS.>/fao=1
noscript        <No script file is specified for /SCRIPT qualifier.>
.END
