.Title       FNT_Msg
.Facility    FNT,780/Prefix=FNT_
.Base 0
.Severity    Success
NORMAL       <Successful completion.>
.Base 1000
.Severity    Informational
MAXRECCOU    <Maximum number of SCI records on READ.  First !SL taken.>/FAO=1
CON_VOLTS    <Error indicated by FNT_convert_to_volts RETURN status>
POW_DENSITY  <Error indicated by FNT_power_density RETURN status>
R_NOISE      <Error indicated by FNT_read_noise RETURN status>
RNOISE_NODATA    <FNT_READ_NOISE: no data to READ.>
RNOISE_INCONSIST <FNT_READ_NOISE: data failed consistency check.>
.Base 2000
.Severity    Error
P_NOISE      <Error indicated by FNT_Parse_NOISE RETURN status>
CTINIT       <Cobetrieve CT_INIT failure.  Status=!SW>/FAO=1
INVTIM       <JSTOP time is less than JSTART time.>
OPNCONFIGERR <Error indicated by CCT_OPEN_CONFIG, !SL>/FAO=1 
GETCONFIGERR <Error indicated by CCT_GET_CONFIG_IDX_TOD, !SL>/FAO=1 
READERR     <Error indicated by FNT_GET_TRANSFCN, failed to read elect. tran.> 
.Base 3000
.Severity    Severe
ABERR        <FNT aborts due to non-recoverable error.>
.END
