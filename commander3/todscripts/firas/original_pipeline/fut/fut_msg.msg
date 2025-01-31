.Title       FUT_Msg
.Facility    FUT,755/Prefix=FUT_
.Base 0
.Severity    Success
NORMAL       <Successful completion.>
ABERR        <Nonrecoverable error causes routine to abort.>
EOF          <End of file reached.>
.Base 1000
.Severity    Informational
NOSYMCTR     <IFG too noisy to locate central peak.>
SIG2SM       <FUT_SWG_MAKE_TEMPLATE: SIGMAS(!SL) IS < .5 . It has been set to .5 .>/FAO=1
WRTPLOTERR   <Bad return detected from FUT_Write_Engplot. Status=!SL>/FAO=1
TRIGPLOTERR  <Bad return detected from FUT_Trigger_Plot. Status=!SL>/FAO=1
SUMQUALBAD   <Alert value !SL of IFG Summary Quality Flag, channel !SL.>/FAO=2
LIBCHECKSM   <Bad return detected from Library Checksum. Status=!SL>/FAO=1
CHECKSUMER   <Bad return detected from FUT_Checksum. Status=!SL>/FAO=1
MJFCHANGER   <Bad return detected from FUT_MJF_Change. Status=!SL>/FAO=1
SUMFLAGER    <Bad return detected from FUT_Sum_Flg. Status=!SL>/FAO=1
QUALSUMER    <Bad return detected from FUT_Qual_Summary. Status=!SL>/FAO=1
DEFPEAK      <Using FUT_DEFAULT_PEAK to determine IFG peak position.>
TOOMANYPTS   <Too many points in plot vector for FUT_Vecplt. >
.Base 2000
.Severity    Warning
OFLWRNA      <FUT_AVE_BIN_TIME Warning: number of records exceeds 100. There is a>
OFLWRNB      <danger of integer overflow. See instructions within the subroutine.>
INVTIM       <JSTOP time is less than JSTART time.>
ICTRBDA      <Index of IFG peak is out of bounds.>
ICTRBDB      <Index= !SL. Nominal values are 11 .LT. Index .LT. NPTS - 10.>/FAO=1
RSEOPEN      <RSE file RMS open failure. Status=!SW>/FAO=1
RSEREAD      <RSE file RMS read failure. Status=!SW>/FAO=1
RSECLOSE     <RSE file RMS close failure. Status=!SW>/FAO=1
.Base 3000
.Severity    Error
CTOPEN       <Cobetrieve CT_OPEN_ARCV failure. Archive=!SW, Status=!SW>/FAO=2
PLOTSCALE    <Unknown plot scale requested.>
OPENCPY      <Display Copy File open error. Status=!SL>/FAO=1
CLOSECPY     <Display Copy File close error. Status=!SL>/FAO=1
UNKPLTDEV    <Unrecognized plot device requested.>
BADTEMPS     <Incompatable temperatures. GRT=!SL, Side !AS low= !ASK, high= !ASK>/FAO=4
ETFRECNO     <Unable to find electronics transfer function. Record number=!SL>/FAO=1
ILLADDS      <Illegal adds per group specified.>
ILLSMP       <Illegal sampling rate specified.>
ILLUNITS     <Illegal units specified. Hz or wavenumbers only.>
ILLMTMSPD    <Illegal MTM speed specified.>
ILLCHAN      <Illegal channel specified.>
ILLMODE      <Illegal MTM mode specified.>
FAKEITHZ     <Fake-it mode spectra can only be in Hz.>
CTQGET_ERR   <Archive query get failure. Status=!SW>/FAO=1
CFROPNFLD    <Field Retriever FR_OPEN_FIELD error. Field=!AS>/FAO=1
CFRATTR      <Field Retriever FR_ATTRIBUTES error. Field=!AS>/FAO=1
NoLUN        <No logical unit numbers available.>
InvLUN       <Requested logical unit number out of range.>
LUNFree      <Logical unit number already free.>
UnkSmth      <Unknown smoothing fucntion requested.>
InvAngle     <Averaged angle does not lie between -PI and PI.>
LUNGETERR    <Error obtaining logical unit number. Status=!SL>/FAO=1
OPENERR      <File open error. Status=!SL, File=!AS>/FAO=2
WRITERR      <File write error. Status=!SL, File=!AS>/FAO=2
CLOSERR      <File close error. Status=!SL, File=!AS>/FAO=2
INVALPLOT    <Invalid type of limits selected for Engplots. Type= !SW>/FAO=1 
INVDEFPEAK   <Invalid default IFG peak calculated by FUT_DEFAULT_PEAK.>
NOCATRECS    <No RSE catalog records found in timerange specified.>
CTQUERYCAT   <CCT_QUERY_TTG_CATALOG failure. Status=!SL>/FAO=1
CTOPENRSE    <RSE archive open failure. Status=!SL>/FAO=1
CTREADRSE    <RSE archive read failure. Status=!SW>/FAO=1
CTCLOSERSE   <RSE archive close failure. Status=!SW>/FAO=1
OPENORBHDR   <Orbit Header archive open failure. Status=!SL>/FAO=1
READORBHDR   <Orbit Header archive read failure. Status=!SL>/FAO=1
CLOSEORBHDR  <Orbit Header archive close failure. Status=!SL>/FAO=1
ATTOPNERR    <Attitude archive open failure.>
ORBOPNERR    <Orbit archive open failure.>
ATTORBOPNERR <Attitude and Orbit archive open failure.>
CombSwErr    <Erroneous combine-switch in FUT_Temp_List call. Should be 0,1,2.>
.END
