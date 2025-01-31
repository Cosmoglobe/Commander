.Title       FTB_Msg
.Facility    FTB,770/Prefix=FTB_
.Base 0
.Severity    Success
NORMAL       <Successful completion.>
.Base 1000
.Severity    Informational
ABERR        <Run aborted.>
BAD_STEPSIZE <Inappropriate step size for FTB_Word31. >
STOP_LT_STAR <Error: Stop is less than start time for FTB_Word31. >
INV_CHOI1    <Choice of dataset not valid for this application.  >
INV_CHOI2    <Choose Raw Science, Noise Spectra>
WRITE        <Dump file write failure. Status=!SW>/FAO=1
CLOSE        <Dump file close failure. Status=!SW>/FAO=1
MAXRECCOU    <Maximum number of SCI records on READ.  First !SL taken.>/FAO=1
TOOMANYPTS   <Too many points in plot vector for FTB_Vecplt. >
CTQARC       <Cobetrieve CT_QUERY_ARCV failure. Status=!SW>/FAO=1
CTQGET       <Cobetrieve CT_QUERY_GET failure. Status=!SW>/FAO=1
CTQARCEOF    <No records found for Rse in CT_Query_Arcv. EOF. Status=!SW>/FAO=1
GET_CHAN     <No input channel equal : RH, RL, LH, LL > 
GET_CHAN_FILE < input fdq file not consistent with the input channel > 
GET_FILE     <no input fdq file >
FPPOVRLAP    <!SL overlaps encountered in candidate FPP run times.>/FAO=1
.Base 2000
.Severity    Warning
CHKWARN     <WARNING: Sience record header time is less than Eng header time>
.Base 3000
.Severity    Error
CTINIT       <Cobetrieve CT_INIT failure. Status=!SW>/FAO=1
CTQBLD       <Cobetrieve CT_QUERY_BLD failure. Status=!SW>/FAO=1
CTOPEN       <Cobetrieve CT_OPEN_ARCV failure. Status=!SW>/FAO=1
CTREAD       <Cobetrieve CT_READ_ARCV failure. Status=!SW>/FAO=1
CTCLOS       <Cobetrieve CT_CLOSE_ARCV failure. Status=!SW>/FAO=1
OPEN         <Dump file open failure. Status=!SW>/FAO=1
OPNFAIL_ANC  <Error: Open failure on ancillary data. Status= !SL>/FAO=1
OPNFAIL_HKP  <Error: Open failure on housekeeping data. Status= !SL>/FAO=1
RSEOPEN      <RSE file RMS open failure. Status=!SW>/FAO=1
RSEWRITE     <RSE file RMS write failure. Status=!SW>/FAO=1
RSEREAD      <RSE file RMS read failure. Status=!SW>/FAO=1
RSECLOSE     <RSE file RMS close failure. Status=!SW>/FAO=1
RMSOPEN      <Error opening file !AS. Status = !SL>/FAO=2
RMSREAD      <Error reading file !AS. Status = !SL>/FAO=2
RMSWRITE     <Error writing file !AS. Status = !SL>/FAO=2
RMSOPENPLT   <Plots data file RMS open failure. Status=!SW>/FAO=1
RMSCLOSEPLT  <Plots data file RMS close failure. Status=!SW>/FAO=1
INVTIM       <JSTART is greater than JSTOP time.>
DAFCREADERR  <Error reading DAFS Record with CSDR Key. Status=!SW>/FAO=1
DAFIREADERR  <Error reading DAFS Record with IGSE Key. Status=!SW>/FAO=1
DSNAMLEN     <Dataset name is too long. Max= 12 characters. Name= !AS>/FAO=1 
MISSDSNAME   <Missing dataset name. Name is required for batch mode.>
OPNCONFIGERR <FTB detects bad return from CCT_Open_Config. Status=!SL>/FAO=1
GETCONFIGERR <FTB detects bad return from CCT_Get_Config. Status=!SL>/FAO=1
CLSCONFIGERR <FTB detects bad return from CCT_Close_Config. Status=!SL>/FAO=1
NODATA       <No data in requested time interval.>
RMSOPENCTL   <RMS Segment Catalog open failure. Chan=!SL, Status=!SL>/FAO=2
RMSMODREADCTL  <RMS Segment Catalog modify read failure. Chan=!SL, Status=!SL>/FAO=2
RMSMODWRITECTL <RMS Segment Catalog modify write failure. Chan=!SL, Status=!SL>/FAO=2
RMSCLOSECTL  <RMS Segment Catalog close failure. Chan=!SL, Status=!SL>/FAO=2
RMSOPENREP   <Report file RMS open failure. Status=!SL, File=!AS>/FAO=2
RMSWRITEREP  <Report file RMS read failure. Status=!SL, File=!AS>/FAO=2
RMSCLOSEREP  <Report file RMS close failure. Status=!SL, File=!AS>/FAO=2
NOCATRECS    <No Raw Science catalog records found.>
CTQUERYCAT   <Raw Science catalog query failure. Chan=!SL, Status=!SW>/FAO=2
OPENSDF      <Raw science segment open failure. Segment=!AS, Status=!SL.>/FAO=2
READSDF      <Raw science segment read failure. Segment=!AS, Status=!SL.>/FAO=2
CLOSESDF     <Raw science segment close failure. Segment=!AS, Status=!SL.>/FAO=2
TOOLONG      <A time interval longer than 60 minutes was requested.>
NFSNAMLEN    <NFS Dataset name is too long. Max = 10. Current length=!SL>/FAO=1
NUMBINS      <Maximum number of analog bins exceeded. Max = 5.>
NOMATCH      <Selected field does not match any plot field name.>
BADFLD       <Invalid engineering field name selected.>
NOBINMATCH   <Engineering value does not fall in any selected bin.>
NOSCI        <No science data is associated with Engineering record.>
.Base 4000
.severity    Fatal 
CTRARC_ERR  <CT_KEYEDREAD fatal error. GMT = !AS .>/FAO=1
.END
