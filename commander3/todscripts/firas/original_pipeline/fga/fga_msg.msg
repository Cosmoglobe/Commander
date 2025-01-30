.Title       FGA_Msg
.Facility    FGA,790/Prefix=FGA_
.Base 0
.Severity    Success
NORMAL       <Successful completion.>
EOF          <End of file reached.>
.Base 1000
.Severity    Informational
WRITE_ERR    <Dump file write failure. Status=!SW>/FAO=1
CLOSE_ERR    <Dump file close failure. Status=!SW>/FAO=1
CTCLOSE_ERR  <Archive close failure. Status=!SW>/FAO=1
DSNAMLEN     <Dataset name is too long. Max= 12 characters. Name= !AS>/FAO=1 
MISSDSNAME   <Missing dataset name. Name is required for batch mode.>
RSEOPEN      <RSE file RMS open failure. Status=!SW>/FAO=1
RSEREAD      <RSE file RMS read failure. Status=!SW>/FAO=1
RSECLOSE     <RSE file RMS close failure. Status=!SW>/FAO=1
NODATA       <No data found in the archive.>
.Base 2000
.Severity    Warning
.Base 3000
.Severity    Error
OPEN_ERR     <Dump file open failure. Status=!SW>/FAO=1
CTINIT_ERR   <Cobetrieve initialization failure. Status=!SW>/FAO=1
CTQBLD_ERR   <Archive query build failure. Status=!SW>/FAO=1
CTOPEN_ERR   <Archive open failure. Status=!SW>/FAO=1
.Base 4000
.severity    Fatal 
CTRARC_ERR  <CT_KEYEDREAD fatal error. GMT = !AS .>/FAO=1
.END
