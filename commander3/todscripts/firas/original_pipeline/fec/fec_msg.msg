! Message definition source file for FEC facility
!
.Title    FEC_MSG Message Definition File 
.Ident    'Version 9.2'
.Facility FEC,715/prefix=FEC_

.Severity SUCCESS         !Success messages

   NORMAL "Normal successful completion"
   CONTINUE  "No errors in module !AS" /FAO_COUNT=1

.Severity INFORMATIONAL   !Informative messages
   NODATA "No engineering data found"

.Severity WARNING         !Warning messages
   TOOMANYTOL "Too many tolerances supplied in command line - using first !AS" /FAO_COUNT=1
   DATAOVERFLOW "Data array overflow - tried to read more records than allowable"

.Severity ERROR           !Error messages

.Severity FATAL           !Fatal messages
   NORSEFILE "RSE Filename not supplied"
   MISSINGTOL "Not enough tolerances provided - need !AS, found !AS" /FAO_COUNT=2
   INVTOL "Invalid tolerance supplied"
   PLATOVERFLOW "Plateau array overflow - tried to add more plateaus than allowed"
   HOTMERGE "Error in merging A side and B side Hot Spot Plateaus"
   USERUNDFND "COBETRIEVE error: User ID not defined"
   NOOPENLUN "COBETRIEVE error:  No open unit number available"
   BADLUN    "COBETRIEVE error:  Bad unit number !AS" /FAO_COUNT=1
   ABORT  "Program terminated due to fatal error in !AS" /FAO_COUNT=1
   ERROR  "Program terminated due to do unexpected error during processing"
   CTQBLDRERROR "Error in CT_QUERY_BLD"
   CTQGETERROR "Error in CT_QUERY_GET. Status=!SW"/FAO_COUNT=1
   CTWRITEERROR "Error in CT_WRITE_ARCV. Status=!SW"/FAO_COUNT=1
   QARCVERR "Error in CT_QUERY_ARCV"
   OPENARCVERR "Error in CT_OPEN_ARCV"
   CTCLOSEENG "Error closing engineering archive"
   CTOPENSS "Error opening short science archive. Channel=!SL" /FAO_COUNT=1
   CTCLOSESS "Error closing short science archive. Channel=!SL" /FAO_COUNT=1
   CLOSEREP "Error closing report file"
   OPENREP "Error opening report file"
   OPNCONFIGERR "FEC detects bad return from CCT_Open_Config. Status=!SL"/FAO=1
   GETCONFIGERR "FEC detects bad return from CCT_Get_Config. Status=!SL"/FAO=1
   CLSCONFIGERR "FEC detects bad return from CCT_Close_Config. Status=!SL"/FAO=1
.End
