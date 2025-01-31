.Title       FEP_Msg
.Facility    FEP,775/Prefix=FEP_
.Base 0
.Severity    Success
NORMAL       <Successful completion.>
ABERR        <Aborted due to error.>
.Base 1000
.Severity    Informational
EXITMENU     <Exit main menu and make plots.>
QUITMENU     <Quit main menu and make plots.>
TRUNCDATA    <Data truncated after 3000 records.>
NODATA       <No data found in the current timerange.>
MONREALTIME  <MONITOR option works only in real-time.>
TBLEDGE      <Temperature lies outside or at edge of GRT table.>
POLYSLNERR   <Estimated error in GRT polynomial solution greater than tolerance.>
MAXFLDSEL    <Maximum field selection limit of !SW was reached.>/FAO=1
.Base 2000
.Severity    Warning
MENUELERR    <Menu element too large for menu buffer.>
.Base 3000
.Severity    Error
CTINIT       <Cobetrieve CT_INIT failure. Status=!SW>/FAO=1
CTOPENHSKP   <Housekeeping data archive open failure. Status=!SW>/FAO=1
CTREADHSKP   <Housekeeping data archive read failure. Status=!SW>/FAO=1
CTCLOSEHSKP  <Housekeeping data archive close failure. Status=!SW>/FAO=1
INVTIM       <JSTOP time is less than JSTART time.>
INVDATALEN   <Invalid data field length.>
GRTTBLOPEN   <GRT reference table RMS open failure. Status=!SW>/FAO=1
GRTTBLREAD   <GRT reference table RMS read failure. Status=!SW>/FAO=1
GRTTBLCLOSE  <GRT reference table RMS close failure. Status=!SW>/FAO=1
GRTREFOPEN   <GRT conversion table RMS open failure. File=!AS, Status=!SW>/FAO=2
GRTREFREAD   <GRT conversion table RMS read failure. Status=!SL,File=!AS>/FAO=2
GRTREFCLOSE  <GRT conversion table RMS close failure. File=!AS, Status=!SW>/FAO=2
CFRINIT      <Field Retriever FR_INIT error.>
CFROPNFLD    <Field Retriever FR_OPEN_FIELD error. Field=!AS>/FAO=1
CFRATTR      <Field Retriever FR_ATTRIBUTES error. Field=!AS>/FAO=1
CFRGETFLD    <Field Retriever FR_GET_FIELDS error.>
CFRCLOSE     <Field Retriever FR_CLOSE error.>
LLSQFERR     <Error returned by IMSL LLSQF function. Status=!SL>/FAO=1
CRESOPEN     <Calibrator resisitor file open failure. Status=!SW>/FAO=1
CRESREAD     <Calibrator resisitor file read failure. Status=!SL. File=!as>/FAO=2
CRESCLOSE    <Calibrator resisitor file close failure. Status=!SW>/FAO=1
GETCONFIGERR <FEP detects bad return from CCT_Get_Config. Status=!SL>/FAO=1
OPNCONFIGERR <FEP detects bad return from CCT_Open_Config. Status=!SL>/FAO=1
CLSCONFIGERR <FEP detects bad return from CCT_Close_Config. Status=!SL>/FAO=1
LIMREAD      <Poly Coeffs. file read failure. Status=!SL. File=!AS>/FAO=2
NOPOLYDEF    <Database word !AS has no poly coefficients. Values set to 0>/FAO=1
DECODEREAD   <STOL database name decode read failure. Status=!SW>/FAO=1
GETPOLYERR   <Error retrieving polynomial conversion coefficients. Status=!SW>/FAO=1
OPENREP      <Report file open failure. Status=!SW>/FAO=1
WRITEREP     <Report file write failure. Status=!SW>/FAO=1
CLOSEREP     <Report file close failure. Status=!SW>/FAO=1
INVTMAT      <Matrix singular. Unable to invert matrix.>
WINDOW       <Window size must be an odd non-zero number.>
RMSOPENXREF  <Housekeeping/STOLDB Xref RMS open failure. File=!AS, Status=!SW>/FAO=2
RMSREADXREF  <Housekeeping/STOLDB Xref RMS read failure. File=!AS, Status=!SW>/FAO=2
RMSCLOSEXREF <Housekeeping/STOLDB Xref RMS close failure. File=!AS, Status=!SW>/FAO=2
.END
