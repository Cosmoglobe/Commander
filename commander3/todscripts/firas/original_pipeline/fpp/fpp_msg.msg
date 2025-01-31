.Title       FPP_Msg
.Facility    FPP,703/Prefix=FPP_
.Base 0
.Severity    Success
Normal       <Successful completion.>
.Base 1000
.Severity    Informational
Aberr        <Run aborted.>
MissSciSel   <Missing NFS_SDF file selection. Filename or time range required.>
InvTim       <JSTART is greater than JSTOP time.>
NoCatRec     <No catalog entry found for segment !AS>/FAO=1
NoCatRecTOD  <No catalog entries found for dataset !AS in requested time range.>/FAO=1
SegDelt      <Raw Science !AS has been deleted. Continuing to search for data to process.>
AncNoMat     <Ancillary file !AS without raw science. Continuing to search for data to process.>
SciNoMat     <*** Raw science !AS without ancillary data noted. Continuing. ***>
ColBefAnc    <Consistency fails - Collection starts before Anc record 1. Start=!AS>/FAO=1
ColAftAnc    <Consistency fails - Collection starts after Anc record 2. Start=!AS>/FAO=1
Overlap      <FPP_SDF segment overlaps found in channel !AS.>/FAO=1
Clsconfigerr <Bad return from CCT_Close_Config. Status=!SL>/FAO=1
.Base 2000
.Severity    Warning
EmulErr      <Error return from Lib$Emul. Status = !SL>/FAO=1
AddxErr      <Error return from Lib$Addx. Status = !SL>/FAO=1
SubxErr      <Error return from Lib$Subx. Status = !SL>/FAO=1
EdivErr      <Error return from Lib$Ediv. Status = !SL>/FAO=1
.Base 3000
.Severity    Error
CTInit       <Cobetrieve CT_INIT failure. Status=!SL>/FAO=1
CTOpenErr    <Cobetrieve OPEN ARCHIVE failure. Status=!SL>/FAO=1
CTReadErr    <Cobetrieve CT_READ_ARCV failure. Status=!SL. File =!AS>/FAO=2
CTWritErr    <Cobetrieve CT_WRITE_ARCV failure. Status=!SL>/FAO=1
CTClosErr    <Cobetrieve CT_CLOSE_ARCV failure. Status=!SL>/FAO=1
CTAncEof     <End of file detected on read of first NFS_ANC record in file.>
CTRefEOF     <End of file detected on read of FEX reference file.>
LateRef      <Collect time was before first reference record.>
OpenErr      <Report file open failure. Status=!SL>/FAO=1
WritErr      <Report file write failure. Status=!SL>/FAO=1
ClosErr      <Report file close failure. Status=!SL>/FAO=1
QCatErr      <Error return from CCT_Query_Catalog. Status=!SL>/FAO=1
QTODCatErr   <Error return from CCT_Query_TOD_Catalog. Status=!SL>/FAO=1
GetLunErr    <Error return from Lib$Get_Lun. Status=!SL>/FAO=1
GettimErr    <Error return from Sys$Gettim. Status=!SL>/FAO=1
Fileseg      <Bad File Segment >
Brefchan     <Bad Reference Channel Selection Ref_Chan=!SL >/FAO=1
Bchannel     <Bad Channel Selection Chan=!AS>/FAO=1
Fnumchan     <Final Processing Channel Number EQ 0. Process Aborted!!>
CTGetFSpec   <Cobetrieve CT_Get_FileSpec_Fields failure. Status=!SL>/FAO=1
ClnCat       <FUT_Clean_Catalog failure. Status=!SL>/FAO=1
CtQueryCat   <FPP Raw Science catalog query failure. Chan=!SW, Status=!SL>/FAO=2
Getconfigerr <Bad return from CCT_Get_Config. Status=!SL>/FAO=1
Opnconfigerr <Bad return from CCT_Open_Config. Status=!SL>/FAO=1
.End
