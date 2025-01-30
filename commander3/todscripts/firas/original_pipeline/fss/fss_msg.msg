! Message definition source file for fss facility
!
.Title    fss_MSG Message Definition File 
.Ident    'Version 9.0'
.Facility fss,721

.Severity SUCCESS 

   NORMAL      "Normal successful completion"

.Severity WARNING

   ALLSCIREJ   "All science records were rejected."
   ALLSKYREJ   "All skymap records were rejected."
   ALLENGREJ   "All engineering records were rejected."
   ALLMAPSREJ1 "All input skymaps rejected. Check JSTART."
   ALLMAPSREJ2 "JSTART must be later than JSTOP of input skymap."
   DUPSKYMAP   "Multiple versions of input skymap. Most recent used."
   DUPTIME     "Duplicate time in skymap. Time=!AS."                   /fao=1
   NOSKYRECS   "No records to write to skymap."
   CSAWRITERR  "Error writing skymap record with CSA.  status = !SL" /FAO=1

.Severity ERROR

   BADFSOPEN   "Error opening full science archive. FORTRAN STAT = !SL"/FAO=1
   BADFSCLOSE  "Error closing full science archive. CT_STATUS = !SL"   /FAO=1
   BADENGOPEN  "Error opening engineering archive. FORTRAN STAT = !SL" /FAO=1
   BADENGCLOSE "Error closing engineering archive.  CT_STATUS = !SL"   /FAO=1
   BADQUERY    "QUERY_ARCV error for eng archive.   CT_STATUS = !SL"   /fao=1
   BADQGET     "QUERY_GET error  for eng archive.   CT_STATUS = !SL"   /fao=1
   BADKEYREAD  "Key read error. status = !SL chan = !SL time = !AS"    /fao=3
   BADSEQRD    "Seqential read error on eng file.   CT_Status = !SL"   /fao=1
   BADREPOPEN  "Error opening report file.       FORTRAN STAT = !SL"   /FAO=1
   BIN_OVERLAP "The temperature bins on the command line overlap."
   DIHED_SORT  "Dihedral temp bins on the command line are out of order."
   CTQUERYCAT  "Catalog query failure.              CT_STATUS = !SL"   /fao=1
   CTINITERR   "Error on CT_INIT."
   NOCATRECS   "No files found in catalog.          File = !AS"        /fao=1
   NOENGRECS   "No engineering records found in catalog."
   NOINFILE    "No input skymap found in catalog."
   PIXELIZATION "Pixelization of new data and input skymap disagree."
   MAXDATA     "Maximum number of input records exceeded."

.Severity FATAL  

   ABERR      "FSS aborts due to non-recoverable error."

.End
