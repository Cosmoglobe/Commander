!=========================================================================
!  Modification History:
!
!  Jan 25, 1990 D. Bouler STX 
!      SPR 5533 "FIT output statement overflows record"
!      FIT_BADTEMPS was causing FIT to abort. The high and low temps were
!      converted to ascii strings in order to call lib$signal, but the
!      fortran internal write was failing. FIT_BADTEMPS had an element 
!      parameter that was not used, so it was removed. 
!
!
!
!
!=========================================================================
.Title       FIT_Msg
.Facility    FIT,785/Prefix=FIT_
.Base 0
.Severity    Success
NORMAL       <Successful completion.>
.Base 1000
.Severity    Informational
EXITMENU     <Exit main menu and make plots.>
QUITMENU     <Quit main menu and make plots.>
TRUNCDATA    <Data truncated after 6000 records.>
NODATA       <No data found in the current timerange.>
TBLEDGE      <Temperature lies outside or at edge of GRT table.>
NOGOODAT     <All orbit records in timerange were invalid.>
.Base 2000
.Severity    Warning
MENUELERR    <Menu element too large for menu buffer.>
.Base 3000
.Severity    Error
CTINIT       <Cobetrieve CT_INIT failure. Status=!SW>/FAO=1
CTOPENSTATS  <Trends data archive open failure. Status=!SW>/FAO=1
CTREADSTATS  <Trends data archive read failure. Status=!SW>/FAO=1
CTCLOSESTATS <Trends data archive close failure. Status=!SW>/FAO=1
INVTIM       <JSTOP time is less than JSTART time.>
INVALBIN     <Bin size is out of range.>
INVDATALEN   <Invalid data field length.>
GRTTBLOPEN   <GRT reference table RMS open failure. Status=!SW>/FAO=1
GRTTBLREAD   <GRT reference table RMS read failure. Status=!SW>/FAO=1
GRTTBLCLOSE  <GRT reference table RMS close failure. Status=!SW>/FAO=1
CFRINIT      <Field Retriever FR_INIT error.>
CFROPNFLD    <Field Retriever FR_OPEN_FIELD error. Field=!AS>/FAO=1
CFRATTR      <Field Retriever FR_ATTRIBUTES error. Field=!AS>/FAO=1
CFRCLOSE     <Field Retriever FR_CLOSE error.>
LLSQFERR     <Error returned by IMSL LLSQF function. Status=!SL>/FAO=1
CRESOPEN     <Calibrator resistor file open failure. Status=!SW>/FAO=1
CRESREAD     <Calibrator resistor file read failure. Status=!SW>/FAO=1
CRESCLOSE    <Calibrator resistor file close failure. Status=!SW>/FAO=1
OPENREP      <Report file open failure. Status=!SW>/FAO=1
WRITEREP     <Report file write failure. Status=!SW>/FAO=1
CLOSEREP     <Report file close failure. Status=!SW>/FAO=1
INVTMAT      <Matrix singular. Unable to invert matrix.>
BADTEMPS     <Incompatible temperatures. Side !AS low= !ASK, high= !ASK>/FAO=3
OPNCONFIGERR <FIT detects bad return from CCT_Open_Config. Status=!SL>/FAO=1
GETCONFIGERR <FIT detects bad return from CCT_Get_Config. Status=!SL>/FAO=1
CLSCONFIGERR <FIT detects bad return from CCT_Close_Config. Status=!SL>/FAO=1
.END
