.TITLE     FSD_MSG
.FACILITY  FSD,760/PREFIX=FSD_

.BASE 0
.SEVERITY  SUCCESS
NORMAL     < Normal successful completion. >
ABORT      < Process run aborts. >

.BASE 1000
.SEVERITY  INFORMATIONAL
FREND      < End of file from IUnit=!SL Status=!SW>/FAO=2

.BASE 2000
.SEVERITY  WARNING
CNODATA    < No data from rms to open. Status =!SW>/FAO=1
TNODATA    < No data from this time range. Time_Range=!SW>/FAO=1
FAERR      < Error from Fut_Attitude. Status=!SW>/FAO=1
MAXIFG     < Maximum number (!SL) of IFGs Reached.>/FAO=1
MAXREC     < Maximum number (!SL) of input records reached.>/FAO=1
ENDTIME    < End Of Time Range Reached. Channel=!SL>/FAO=1
LPIXTER    < Limit Reached for Pix_terr position data discarded. NF=!SL>/FAO=1
ERRSIFG    < Error on smoothing IFGS. Status=!SW>/FAO=1
ERRCIFG    < Error comparing IFG to template. Status=!SW>/FAO=1
ILLORDFIL  < Input file !AS not properly ordered.>/FAO=1

.BASE 3000
.SEVERITY  ERROR
FOPEN      < Open failure from IUnit=!SL IOS =!SW>/FAO=2
FRERR      < Read error from IUnit=!SL   IOS =!SW>/FAO=2
CTINIT     < Cobetrieve CT_INIT failure. Status=!SW>/FAO=1
CTOPEN     < Cobetrieve CT_CONNECT_READ failure. Channel=!SL, Status=!SW>/FAO=2
CTQBLD     < Cobetrieve CT_QUERY_BLD failure. Channel=!SL, Status=!SW>/FAO=2
CTREAD     < Cobetrieve CT_READ_ARCV failure. Status=!SW>/FAO=1
CTREOF     < Cobetrieve CT_READ_ARCV end of file. Status=!SW>/FAO=1
CTQGET     < Cobetrieve CT_QUERY_GET failure. Status=!SW>/FAO=1
CTQEOF     < Cobetrieve CT_QUERY_GET end of file. Status=!SW>/FAO=1
CTCLOS     < Cobetrieve CT_CLOSE_ARCV failure. Status=!SW>/FAO=1
OPNCONFIGERR <FSD detects bad return from CCT_Open_Config. Status=!SL>/FAO=1
GETCONFIGERR <FSD detects bad return from CCT_Get_Config. Status=!SL>/FAO=1
NOCATREC   < No Cobetrieve catalog records. File=!AS>/FAO=1
QCATERR    < Cobetrieve catalog query error. Status=!SW>/FAO=1
INVTIM     < Invalid time selection. >

.BASE 4000
.SEVERITY  SEVERE
ABERR1     < FSD_ASTROPLOTS aborts due to a non-recoverable error.>
ABERR2     < FSD_POSERR aborts due to a non-recoverable error.>
ABERR3     < FSD_GLITCHMAP aborts due to a non-recoverable error.>
.END
