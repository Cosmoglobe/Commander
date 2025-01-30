.title            fad_msg
.facility         fad, 718 /prefix=fad_
.base 0
.severity         success
 normal           <Successful completion.>
.base 1000
.severity         informational
 abort            <Process aborted.>
 model_mismatch   <FEX_EJ has different model solution than FCF data.>
.base 2000
.severity         warning
 parserr          <Error parsing command line.>
 lunerr           <Error getting logic unit number; status !SL.>/fao=1
 openerr          <Error opening reference file !AS; status !SL.>/fao=2
 readerr          <Error reading reference file !AS; status !SL.>/fao=2
 closerr          <Error closing reference file !AS; status !SL.>/fao=2
 ctinit           <Cobetrieve ct_init failure; status !SW.>/fao=1
 openrep          <Error opening report file.; status !SL.>/fao=1
 closerep         <Error closing report file.; status !SL.>/fao=1
 csaopen          <Error opening skymap !AS; status !SL.>/fao=2
 csaoffset        <Error calling csa_field_offset_values; status !SL.>/fao=1
 csaread          <Error calling csa_read_pixels; status !SL.>/fao=1
 csawrite         <Error calling csa_write_pixels; status !SL.>/fao=1
 csaclose         <Error closing skymap !AS; status !SL.>/fao=2
 query_cat	  <Error querying archive catalog; status !SL.>/fao=1
.end
