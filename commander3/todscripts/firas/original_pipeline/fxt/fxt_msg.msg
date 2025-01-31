.Title       FXT_Msg
.Facility    FXT,702/Prefix=FXT_
.Base 0
.Severity    Success
Normal       <Successful completion.>
.Base 1000
.Severity    Informational
Aberr        <Run aborted.>
MissEngSel   <Missing FDQ_Eng file selection. Filename or time range required.>
InvTim       <JSTART is greater than JSTOP time.>
InvPlot      <Invalid plot type selected. Valid types: 0=none, 1=line, 2=laser>
NoCatRec     <No catalog entry found for requested dataset file.>
NoCatRecTOD  <No catalog entries found for TOD dataset in requested time range.>
TooManyCat   <Too many catalog records in selected time range. Number=!SL>/FAO=1
AllCatRecDel <All catalog records in requested time range are deleted.>
.Base 2000
.Severity    Warning
.Base 3000
.Severity    Error
CTInit       <Cobetrieve CT_INIT failure. Status=!SW>/FAO=1
CTOpenErr    <Cobetrieve CT_OPEN_ARCV failure. Status=!SW>/FAO=1
CTReadErr    <Cobetrieve CT_READ_ARCV failure. Status=!SW>/FAO=1
CTWritErr    <Cobetrieve CT_WRITE_ARCV failure. Status=!SW>/FAO=1
CTClosErr    <Cobetrieve CT_CLOSE_ARCV failure. Status=!SW>/FAO=1
CCTSetTtgErr <Error return from CCT_Set_Ttg_Time_Range. Status=!SW>/FAO=1
CCTOpnConfig <Error return from CCT_Open_Config. Status=!SW>/FAO=1
CCTGetConfig <Error return from CCT_Get_Config. Status=!SW>/FAO=1
CCTCloConfig <Error return from CCT_Close_Config. Status=!SW>/FAO=1
OpenErr      <Report file open failure. Status=!SW>/FAO=1
WritErr      <Report file write failure. Status=!SW>/FAO=1
ClosErr      <Report file close failure. Status=!SW>/FAO=1
RMSOpenPlt   <Plots data file RMS open failure. Status=!SW>/FAO=1
RMSWritePlt  <Plots data file RMS write failure. Status=!SW>/FAO=1
RMSClosePlt  <Plots data file RMS close failure. Status=!SW>/FAO=1
QCatErr      <Error return from CCT_Query_Catalog. Status=!SW>/FAO=1
QTODCatErr   <Error return from CCT_Query_TOD_Catalog. Status=!SW>/FAO=1
GetLunErr    <Error return from Lib$Get_Lun. Status=!SW>/FAO=1
GettimErr    <Error return from Sys$Gettim. Status=!SW>/FAO=1
.End
