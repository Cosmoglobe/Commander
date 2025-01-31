.Title       FDQ_Msg
.Facility    FDQ,730/Prefix=FDQ_
.Base 0
.Severity    Success
NORMAL       <Successful completion.>
.Base 1000
.Severity    Informational
SELFILE      <Segment selected by filename. File= !AS>/FAO=1
REFCHAN	     <Selected reference channel = !AS>/FAO=1
SELSEGNO     <Segment selected by segment catalog number. Catno=!SL>/FAO=1
SELTIMINT    <Data selected by time interval. JSTART=!AS, JSTOP=!AS>/FAO=2
SELASCFUL    <Full Ascii file printout of conversion coeffs selected.>
SELASCBRF    <Brief Ascii file printout of conversion coeffs selected.>
SELSEC       <Selected number of seconds for set grouping=!SL>/FAO=1
SELBST       <User selected binary data start time range for processing.>
SELDDT       <User selected data delivery time range for processing.>
DDTNOTIMP    <Selection of data by data delivery time is not implemented yet>
SELEXCL      <User excludes segments from processing.>
SELNUMEXC    <Selected number of segments to be excluded = !SL>/FAO=1
NUMSEGEXC    <Count of number of segments to be excluded = !SL>/FAO=1
SELNTOBIG    <Num. segments,!SL, excluded too big. List will be counted>/FAO=1
TYPATTSEL    <Selected attitude solution type = !AS>/FAO=1
SCSIDESEL    <Selected spacecraft side = !SL>/FAO=1
PLOTSEL      <Selected Engplots type for exceeding limits = !AS>/FAO=1 
NOSCIMATCH   <No science files found to match housekeeping file: !AS>/FAO=1
SEGDELETE    <Segment deleted from archive catalog. HKP file: !AS>/FAO=1
FPPRUN       <FPP not yet run on raw science. HKP file: !AS>/FAO=1
FDQRUN       <FDQ already run on FPP raw science. HKP file: !AS>/FAO=1
CTHKPFIR     <First failure on NFS_HKP archive open. Status=!SL, Retry>/FAO=1
CTHKPSUC     <Successful open of NFS_HKP archive on second try. DDT=!AS>/FAO=1
CURSEG       <Current segment being processed. START=!SL, STOP=!SL>/FAO=2 
MARKSEGMOD   <Marking excluded segment !SL as modified.>/FAO=1
NHSKREAD     <Number of successful reads for HSK records = !SL>/FAO=1
NCALREAD     <No Cal.resistor record found for Fex_AV_calrs>
NEWHKPSEG    <New housekeeping segment begins. Catno.= !SL>/FAO=1
NHKPREAD     <Number of successful reads for housekeeping record= !SL>/FAO=1
BADTIMINV    <JSTOP time is less than JSTART time.>
MISSHKP      <Corres. housekeeping data not found, Science IFG not processed.>
EOFHKP       <EOF in housekeeping file, Science IFG not processed.>
HKPRFAIL     <Failed to read housekeeping record, Science IFG not processed.> 
CATNOENTRY   <Catalog entry does not exist. File= !AS>/FAO=1
NOSCIENCE    <No raw science channel data files found for selected time range.>
WRITEIDX     <New index record to be written.>
NOIDXRPT     <Due to I/O errors, no index file will be written for segment.>
IDXFLGTRUE   <Due to error in getting index flags, all flags are set to true.>
IDXTOLONE    <Due to error getting index tolerances, all tolerances set to 1.> 
WITHINTRUE   <FDQ Within range is True.>
WITHINFALS   <FDQ Within range is False.>
FLGBADSCI    <Red flag set in segment for data ready bits. Count=!SL>/FAO=1
FLGBADHKP    <Red flag set in segment for HKP telem quality. Count=!SL>/FAO=1
FLGCAL       <Red flag set in segment: non-science-cal data. Count=!SL>/FAO=1
FLGMICRO     <Red flag set in segment for micro status word. Count=!SL>/FAO=1
FLGATTRED    <Red flag set in segment for attitude alert. Count=!SL>/FAO=1
FLGATTYEL    <Yellow flag set in segment for attitude alert. Count=!SL>/FAO=1
FLGGLTCNT    <Red flag set in segment for high glitch count. Count=!SL>/FAO=1
FLGSAT       <Red flag set in segment for high sat. sample. Count=!SL>/FAO=1
FLGCKSUM     <Red flag set in segment for checksum error. Count=!SL>/FAO=1
FLGGRTRED    <Red flag set in segment for GRTs and Calres. Count=!SL>/FAO=1
FLGGRTYEL    <Yellow flag set in segment for GRTs and Calres. Count=!SL>/FAO=1
FLGTCRED     <Red flag set in segment for TCs . Count=!SL>/FAO=1
FLGTCYEL     <Yellow flag set in segment for TCs. Count=!SL>/FAO=1
FLGIPDUR1    <Red flag set in segment IPDU temp. same side. Count=!SL>/FAO=1
FLGIPDUY1    <Yellow flag set in segment IPDU temp. same side. Count=!SL>/FAO=1
FLGIPDUY2    <Yellow flag set in segment IPDU temp. opp. side. Count=!SL>/FAO=1
FLGCHNTPR    <Red flag set in segment for channel temp. Count=!SL>/FAO=1
FLGCHNTPY    <Yellow flag set in segment for channel temp. Count=!SL>/FAO=1
FLGDBOXR1    <Red flag set in segment for drive tmp. same side. Count=!SL>/FAO=1
FLGDBOXY1    <Yel flag set in segment for drive tmp. same side. Count=!SL>/FAO=1
FLGDBOXY2    <Yel flag set in segment for drive tmp. opp. side. Count=!SL>/FAO=1
FLGSMTR      <Red flag set in segment for status monitor tmp. Count=!SL>/FAO=1
FLGSMTY      <Yel flag set in segment for status monitor tmp. Count=!SL>/FAO=1
FLGOPATR     <Red flag set in segment for opt. pre-amp tmp. Count=!SL>/FAO=1
FLGOPATY     <Yel flag set in segment for opt. pre-amp tmp. Count=!SL>/FAO=1
FLGCPATR     <Red flag set in segment for chan. pre-amp tmp. Count=!SL>/FAO=1
FLGCPATY     <Yel flag set in segment for chan. pre-amp tmp. Count=!SL>/FAO=1
FLGHSHAR     <Red flag set in segment for hot spot heater A. Count=!SL>/FAO=1
FLGHSHAY     <Yel flag set in segment for hot spot heater A. Count=!SL>/FAO=1
FLGHSHBR     <Red flag set in segment for hot spot heater B. Count=!SL>/FAO=1
FLGHSHBY     <Yel flag set in segment for hot spot heater B. Count=!SL>/FAO=1
FLGMTMMTR1   <Red flag set in segment for MTM tmp. act. side. Count=!SL>/FAO=1
FLGMTMMTY1   <Yel flag set in segment for MTM tmp. act. side. Count=!SL>/FAO=1
FLGMTMMTY2   <Yel flag set in segment for MTM tmp. opp. side. Count=!SL>/FAO=1
FLGBOLVR     <Red flag set in segment for bolometer voltage. Count=!SL>/FAO=1
FLGBOLVY     <Yel flag set in segment for bolometer voltage. Count=!SL>/FAO=1
FLGIPDUVR    <Red flag set in segment for IPDU voltages. Count=!SL>/FAO=1
FLGIPDUVY    <Yellow flag set in segment for IPDU voltages. Count=!SL>/FAO=1
FLGIPDUCR    <Red flag set in segment for IPDU currents. Count=!SL>/FAO=1
FLGIPDUCY    <Yellow flag set in segment for IPDU currents. Count=!SL>/FAO=1
FLGSCMJF     <Red flag set in segment for stat. chng. acr. MJF. Count=!SL>/FAO=1
FLGDAHY      <Yel flag set in segment for Dirbe anneal. htrs. Count=!SL>/FAO=1
FLGDAHR      <Red flag set in segment for Dirbe anneal. htrs. Count=!SL>/FAO=1
FLGDIRSY     <Yel flag set in segment for Dirbe IRS sources. Count=!SL>/FAO=1
FLGDIRSR     <Red flag set in segment for Dirbe IRS sources. Count=!SL>/FAO=1
FLGDCHPY     <Yel flag set in segment for Dirbe chopper. Count=!SL>/FAO=1
FLGDCHPR     <Red flag set in segment for Dirbe chopper. Count=!SL>/FAO=1
FLGDSHTY     <Yel flag set in segment for Dirbe shutter. Count=!SL>/FAO=1
FLGDSHTR     <Red flag set in segment for Dirbe shutter. Count=!SL>/FAO=1
FLGDLOY      <Yel flag set in segment for DMR local osc. Count=!SL>/FAO=1
FLGDLOR      <Red flag set in segment for DMR local osc. Count=!SL>/FAO=1
FLGDNSY      <Yel flag set in segment for DMR noise sources. Count=!SL>/FAO=1
FLGDNSR      <Red flag set in segment for DMR noise sources. Count=!SL>/FAO=1
FLGSCSY      <Yel flag set in segment for spacecraft status. Count=!SL>/FAO=1
FLGSCSR      <Red flag set in segment for spacecraft status. Count=!SL>/FAO=1
FLGSUM       <Alert count for overall quality summary in segment=!SL>/FAO=1
SUMQUALBAD   <Alert value !SL of IFG Summary Quality Flag, channel !SL.>/FAO=2
CKENTRYSEG   <FDQ_Check_Entry prevents segment from being processed: !AS>/FAO=1
GETOPTER     <FDQ detects bad return from FDQ_Get_Options. Status=!SL>/FAO=1
LOADCONVER   <FDQ detects bad return from FDQ_Load_Conv. Status=!SL>/FAO=1
LOADLIMER    <FDQ detects bad return from FDQ_Load_Lim. Status=!SL>/FAO=1
GETNXSEGER   <FDQ detects bad return from FDQ_Get_Nxt_Seg. Status=!SL>/FAO=1
GETDDTER     <FDQ detects bad return from FDQ_Get_DDT. Status=!SL>/FAO=1
GETVALSEG    <FDQ detects bad return from FDQ_Get_Valid_Segments.>
GETRANGE     <FDQ detects bad return from FDQ_Get_Time_Range. Status=!SL>/FAO=1
WITHINER     <FDQ detects bad return from FDQ_Within. Status=!SL>/FAO=1
GETNXSETER   <FDQ detects bad return from FDQ_Get_Nxt_Set. Status=!SL>/FAO=1
GETBRKHKER   <FDQ detects bad return from FDQ_Get_Brk_Hsk. Status=!SL>/FAO=1
OPENARCVER   <FDQ detects bad return from FDQ_Open_Arcv. Status=!SL>/FAO=1
PROCURSETR   <FDQ detects bad return from FDQ_Proc_Cur_Set. Status=!SL>/FAO=1
CONVERTER    <FDQ detects bad return from FDQ_Convert. Status=!SL>/FAO=1
OLPOLYCVER   <FDQ detects bad return from FDQ_Ol_Poly_Cvt. Status=!SL>/FAO=1
OLGRTCVTER   <FDQ detects bad return from FDQ_Ol_Grt_Cvt. Status=!SL>/FAO=1
GRTLOOKUP    <Bad return from GRT_Lookup. Status=!SL, MJF=!SL, GRT=!AS>/FAO=3
GRTCORRER    <FDQ detects bad return from GRT_Correction. Status=!SL>/FAO=1
INTERPOLER   <FDQ detects bad return from FDQ_Interpolate. Status=!SL>/FAO=1
MAKEIDXER    <FDQ detects bad return from FDQ_Make_Idx. Status=!SL>/FAO=1
GETENGERR    <FDQ detects bad return from FDQ_Get_Eng_Flds. Status=!SL>/FAO=1
GETHSKERR    <FDQ detects bad return from FDQ_Get_Hsk_Flds. Status=!SL>/FAO=1
GETSCIERR    <FDQ detects bad return from FDQ_Get_Sci_Flds. Status=!SL>/FAO=1
NEWIDXER     <FDQ detects bad return from FDQ_New_Idx. Status=!SL>/FAO=1
ATTITUDER    <FDQ detects bad return from FUT_Attitude. Status=!SL>/FAO=1
SIMQUATER    <FDQ detects bad return from FDQ_Sim_Quat. Status=!SL>/FAO=1
SIMPIXNOER   <FDQ detects bad return from FDQ_Sim_Pixno. Status=!SL>/FAO=1
MJFCHANGER   <FDQ detects bad return from FDQ_MJF_Change. Status=!SL>/FAO=1
SUMFLAGER    <FDQ detects bad return from FDQ_Sum_Flg. Status=!SL>/FAO=1
QUALSUMER    <FDQ detects bad return from FDQ_Qual_Summary. Status=!SL>/FAO=1
PRINTXYB     <FDQ detects bad return from FDQ_Print_XYB. Status=!SL>/FAO=1
SAVCATINF    <FDQ detects bad return from FDQ_Sav_Cat_Info. Status=!SL>/FAO=1
CHECKENTRY   <FDQ detects bad return from FDQ_Check_Entry. Status=!SL>/FAO=1
GETQUALFLG   <FDQ detects bad return from FDQ_Get_Qualflags. Status=!SL>/FAO=1
WRITEMATRIX  <FDQ detects bad return from FDQ_Write_Matrix. Status=!SL>/FAO=1
CHECKSUMER   <FDQ detects bad return from FDQ_Checksum. Status=!SL>/FAO=1
LIBCHECKSM   <FDQ detects bad return from Library Checksum. Status=!SL>/FAO=1
GENSTATER    <FDQ detects bad return from FDQ_Gen_Stat. Status=!SL>/FAO=1
FRSTSTATER   <FDQ detects bad return from first FDQ_Stat_Init. Status=!SL>/FAO=1
STATINITER   <FDQ detects bad return from FDQ_Stat_Init. Status=!SL>/FAO=1
FORMSTATER   <FDQ detects bad return from FDQ_Form_Stat. Status=!SL>/FAO=1
CLOSARCVER   <FDQ detects bad return from FDQ_Close_Arcv. Status=!SL>/FAO=1
WRTPLOTERR   <FDQ detects bad return from FDQ_Write_Engplot. Status=!SL>/FAO=1
TRIGPLOTERR  <FDQ detects bad return from FDQ_Trigger_Plot. Status=!SL>/FAO=1
GETSCCONFERR <FDQ detects bad return from FDQ_Get_SC_Config. Status=!SL>/FAO=1
GETCONFIGERR <FDQ detects bad return from CCT_Get_Config. Status=!SL>/FAO=1
OPNCONFIGERR <FDQ detects bad return from CCT_Open_Config. Status=!SL>/FAO=1
CLSCONFIGERR <FDQ detects bad return from CCT_Close_Config. Status=!SL>/FAO=1
SETLIMFLGERR <FDQ detects bad return from FDQ_Set_Limflags. Status=!SL>/FAO=1
ABERR        <FDQ aborts due to non-recoverable error.>
NUMCATINFO   <FDQ only processes!SL Cat entries. DID=!SL, Seg. trunc=!SL>/FAO=3
MODNEZERO    <FDQ_Check_Entry indicates requested segment already modified.>
GRTZERO      <Warning from Grt_Lookup. Result of GRT lookup is 0.0 value.>
NOCATRECS    <No Raw Science catalog records found.>
GRTEXCTBL    <GRT !AS Res. for MJF !SL exceeds table bounds. GRT val.=!AS>/FAO=3
GRTOLNOSAT   <GRT !AS for MJF !SL temp. soln. !AS may not be good. DR=!AS>/FAO=4
NOCALRESEX   <No calibrator resistors in range for GRT: !AS, MJF=!SL>/FAO=2
hgetconfgerr <hskp_conv: bad return from CCT_Get_Config. Status=!SL>/FAO=1
hopnconfgerr <hskp_conv: bad return from CCT_Open_Config. Status=!SL>/FAO=1
hclsconfgerr <hskp_conv: bad return from CCT_Close_Config. Status=!SL>/FAO=1
hskldconverr <hskp_conv: bad return from h_load_conv. Status=!SL>/FAO=1
hskrdhskerr  <hskp_conv: bad return from h_read_hsk. Status=!SL>/FAO=1
hskconverr   <hskp_conv: bad return from h_conv. Status=!SL>/FAO=1
hskpktemperr <hskp_conv: bad return from h_pack_temp. Status=!SL>/FAO=1
hgetqlflgerr <hskp_conv: bad return from h_get_qualflag. Status=!SL>/FAO=1
hfuttmrngerr <hskp_conv: can't get interval from futtimerange.>
hopenarcverr <hskp_conv: fails to open archives.>
hctwrtac1err <hskp_conv: fails to write eng. rec.#1. ct_stat(1)=!SL>/FAO=1
hctwrtac2err <hskp_conv: fails to write eng. rec.#2. ct_stat(1)=!SL>/FAO=1
hwrtidxrcerr <hskp_conv: fails to write last idx record. ios=!SL>/FAO=1
hidxlunclerr <hskp_conv: fails to close idx. lun. ct_stat(1)=!SL>/FAO=1
henglunclerr <hskp_conv: fails to close eng. lun. ct_stat(1)=!SL>/FAO=1
hsklunclserr <hskp_conv: fails to close hsk. lun. ct_stat(1)=!SL>/FAO=1
hreplunclerr <hskp_conv: fails to close report lun. status=!SL>/FAO=1
hskmkidx1err <hskp_conv: fails to make idx rec.#1 status=!SL>/FAO=1
hskmkidx2err <hskp_conv: fails to make idx rec.#2 status=!SL>/FAO=1
hidxwrterr   <h_make_idx: fails to write idx file. ios=!SL>/FAO=1
hopctiniterr <h_open_arcv: bad return from ct_init. ct_stat(1)=!SL>/FAO=1
hopctopnerr  <h_open_arcv: bad ret. fr. ct_open_arcv. ct_stat(1)=!SL>/FAO=1
hopengopnerr <h_open_arcv: error in open eng. file. ios=!SL>/FAO=1
hopidxopnerr <h_open_arcv: error in open idx. file. ios=!SL>/FAO=1
hoprepopnerr <h_open_arcv: error in open report file. ios=!SL>/FAO=1
hrdctrdacerr <h_read_hsk: bad ret. from ct_read_arcv. ct_stat(1)=!SL>/FAO=1
hoplmcuvkerr <h_load_conv: fails to open limcuvkey file. ios=!SL>/FAO=1
hpolynames   <h_load_conv: database name !AS has no poly.coef. defined>/FAO=1
hpolynameerr <h_load_conv: fails to read limcuvkey for !AS. ios=!SL>/FAO=2
hgrtnames    <h_load_conv: database name !AS has no grt lookup defined>/FAO=1
hgrtnameerr  <h_load_conv: fails to read grt db file for !AS. ios=!SL>/FAO=2
hopgrtlunerr <h_load_conv: fails to open grt db file. ios=!SL>/FAO=1
haberr       <hskp_conv aborts due to non-recoverable error.>
hfilexterr   <hskp_conv detects bad return from h_open_filext. Status=!SL>/FAO=1
hgetcnidxerr <hskp_conv: bad return from CCT_Get_Config_IDX. Status=!SL>/FAO=1
hinvalplot   <Invalid type of limits selected for Engplots. Type= !SW>/FAO=1 
hqualsumerr  <h_get_qualflags: bad return from h_qual_summary. Status=!SL>/FAO=1
htrigploterr <h_get_qualflags: bad return from h_trigger_plot. Status=!SL>/FAO=1
hqualwreperr <h_qual_summary: error in writing report file. ios=!SL>/FAO=1
hxtengopnerr <h_open_filext: error in open ancillary eng. file. ios=!SL>/FAO=1
hxtidxopnerr <h_open_filext: error in open ancillary eng. file. ios=!SL>/FAO=1
BADENG       < Eng. Record GMT Time = !AS less or equal than the prev. Eng. Record. It has been discarded>/FAO=1
.Base 2000
.Severity    Warning
STATLOOPER   <FDQ_Form_Stat detects error in loop index. # loops=!SL>/FAO=1
STATMODER    <FDQ_Form_Stat detects error in modulus =!SL>/FAO=1
FREEVMFAIL   <Failure in Lib$Free_VM. Status=!SL, Number Cat. Entries=!SL>/FAO=2
.Base 3000
.Severity    Error
ERSELFIL     <Error in obtaining user file selection. Status=!SL>/FAO=1
ERSELSEG     <Error in obtaining user seg. catno. selection. Status=!SL>/FAO=1
ERTIMINT     <Error in obtaining user selected time interval. Status=!SL>/FAO=1
SELNOTIM     <Error: No time range selected and no file-segment specified.>
FILEBLANK    <Error: Input filename is blank.>
BADFILE      <Error: Input file !AS does not conform to name convention.>/FAO=1
BADREFCHAN   <Error: Input file !AS does not contain valid channel ID.>/FAO=1
ERSELNEX     <Error obtaining number of segments to exclude. Status=!SL>/FAO=1
EXCSEGMIS    <User requested segments to be excluded but list is missing.>
ERSELEXC     <Error in obtaining selected segments to exclude. Status=!SL>/FAO=1
NEXZERO      <Count of number of segments from exclude list is zero.>
ERSELSEC     <Error in obtaining user defined seconds for set. Status=!SL>/FAO=1
SECTOOBIG    <User requested !SL seconds for IFG grouping out of range.>/FAO=1  
SELORBPD     <Error getting orbital period from command line. Status=!SL>/FAO=1
GETORBPD     <Error reading orbital period from NOR_ORB_HDR. Status=!SL>/FAO=1
BADFILNAM    <FDQ_Open_Arcv detects invalid filename: !AS>/FAO=1
NOFILEX      <FDQ_Open_Arcv detects no file extension for file: !AS>/FAO=1 
CTINITERR    <Error in Cobetrieve initialization. Status=!SL>/FAO=1
CTOPENERR    <Error opening archive. Status=!SL, Data set=!AS>/FAO=2
CTCLOSERR    <Error closing archive. Status=!SL, Data set=!AS>/FAO=2
CTREADERR    <Error reading archive record. Status=!SL, Data set=!AS>/FAO=2
CTWRITERR    <Error writing archive record. Status=!SL, Data set=!AS>/FAO=2
LUNGETERR    <Error obtaining logical unit number. Status=!SL>/FAO=1
LUNFREER     <Error in attempt to free logical unit number. Status=!SL>/FAO=1
OPENERR      <File open error. Status=!SL, File=!AS>/FAO=2
CLOSERR      <File close error. Status=!SL, File=!AS>/FAO=2
READERR      <File read error. Status=!SL, File=!AS>/FAO=2
WRITERR      <File write error. Status=!SL, File=!AS>/FAO=2
ERSUBX       <Error in lib. routine to subtract quad numbers. Status=!SL>/FAO=1
ERADDX       <Error in lib. routine to add quad numbers. Status=!SL>/FAO=1
MJFERR       <FDQ_Convert detects first major frame out of range. MJF=!SL>/FAO=1
AVBINTIMER   <FDQ_Form_Stat: bad return from FUT_Ave_Bin_Times. Stat=!SL>/FAO=1
NEGRADER     <FDQ_Form_Stat detects negative radicand for sigma. Rad=!SL>/FAO=1
GETCATENTF   <Failed to get catalog entry. Status= !SL, File= !AS>/FAO=2
BADATIDX     <Bad dataset array index !SL passed to FDQ_SAVE_CAT_INFO>/FAO=1
CHKENTRYER   <Segment entered is not processable due to: !AS>/FAO=1
CTCATINFOM   <FDQ detects bad return from CT_Cat_Info. Msg= !AS>/FAO=1
GETCATINFO   <Bad status from CT_Get_Cat_Info. Msg=!AS>/FAO=1
DSIDNERS1    <FDQ_Check_Entry indicates data set ID is not ID of Raw Science 1>
GETATTYPERR  <Bad return getting user selected attitude type. Status= !SL>/FAO=1
INVALATTSEL  <Invalid type of attitude solution selected: !AS>/FAO=1
GETSCSIDERR  <Bad return getting user chosen spacecraft side. Status=!SL>/FAO=1
INVALSCSIDE  <Invalid spacecraft side selected by user. Side= !SL>/FAO=1
GETPLOTERR   <Bad return getting type of limits for Engplot. Status= !SL>/FAO=1
INVALPLOT    <Invalid type of limits selected for Engplots. Type= !SW>/FAO=1 
NOPOLYDEF    <Database word !AS has no poly coefficients. Values set to 0>/FAO=1
NOGRTLOOK    <Database word !AS has no GRT lookup defined.>/FAO=1
POLYCVTER    <Bad return from Lib$Polyf. Status=!SL, MJF=!SL, DB NAME=!AS>/FAO=3
ENGANLGDIF   <FDQ detects bad eng. analog diff. Status= !SL>/FAO=1
ERBINTIM     <FDQ detects bad return from Sys$Bintim. Status= !SL>/FAO=1
CONVTULER    <FDQ detects bad return from OTS$CVT_TU_L. Status = !SL>/FAO=1
CONVTFER     <FDQ detects bad return from OTS$CVT_T_F. Status = !SL>/FAO=1
ASCTIMER     <Error: sys conversion from binary to Ascii time. Status=!SL>/FAO=1
ERGETTIM     <FDQ detects bad return from Sys$Gettim. Status= !SL>/FAO=1 
FDQEXIT      <FDQ_Exit invoked. Exit Status= !SL>/FAO=1
CTINIT       <CT_INIT failure. Status=!SW>/FAO=1
CTQUERYCAT   <Raw Science catalog query failure. Chan=!SL, Status=!SW>/FAO=2
CTOPENCTL    <Raw Science Segment Catalog archive open failure. Chan=!SL, Status=!SW>/FAO=2
CTWRITECTL   <Raw Science Segment Catalog archive write failure. Chan=!SL, Status=!SW>/FAO=2
CTCLOSECTL   <Raw Science Segment Catalog archive close failure. Chan=!SL, Status=!SW>/FAO=2
RMSOPENCTL   <RMS Segment Catalog open failure. Chan=!SL, Status=!SL>/FAO=2
RMSREADCTL   <RMS Segment Catalog read failure. Chan=!SL, Status=!SL>/FAO=2
RMSWRITECTL  <RMS Segment Catalog write failure. Chan=!SL, Status=!SL>/FAO=2
RMSMODREADCTL  <RMS Segment Catalog modify read failure. Chan=!SL, Status=!SL>/FAO=2
RMSMODWRITECTL <RMS Segment Catalog modify write failure. Chan=!SL, Status=!SL>/FAO=2
RMSCLOSECTL  <RMS Segment Catalog close failure. Chan=!SL, Status=!SL>/FAO=2
WRTSCIMTR    <Science Matrix CT_WRITE_ARCV failure. Chan=!SL, Status=!SW>/FAO=2
RMSOPENREP   <Report file RMS open failure. Status=!SL, File=!AS>/FAO=2
RMSWRITEREP  <Report file RMS read failure. Status=!SL, File=!AS>/FAO=2
RMSCLOSEREP  <Report file RMS close failure. Status=!SL, File=!AS>/FAO=2
CLNCAT       <FUT_Clean_Catalog failure. Status=!SL>/FAO=1
.END
