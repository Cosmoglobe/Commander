.title              fms_msg
.facility           fms, 799 /prefix=fms_
.base 0
.severity           success
 normal             <Successful completion.>
.base 1000
.severity           informational
 abort              <Process aborted.>
.base 2000
.severity           warning
 parserr            <Error parsing: !AS>/fao=1
 lunerr             <Error getting logic unit number; status !SL.>/fao=1
 openerr            <Error opening file !AS; status !SL.>/fao=2
 readerr            <Error reading file !AS; status !SL.>/fao=2
 writerr            <Error writing to file !AS; status !SL.>/fao=2
 closerr            <Error closing file !AS; status !SL.>/fao=2
 renamerr           <Error renaming file !AS to !AS; status !SL.>/fao=3
 num_maps           <Number of skymaps named in script file must be 2.>
 query_cat          <Error querying catalog; status !SL.>/fao=1
 TTGseterr          <Error setting time tags for file !AS>/fao=1
 ctinit             <Cobetrieve ct_init failure; status !SW.>/fao=1
 illegal_merge      <Illegal merge  !AS + !AS. Consult documentation. >/fao=2
 wrong_order        <Command line ORDER does not match input skymaps + 1.>
 wrong_combo        <Command line COMBINATION does not match input skymaps.>
 zeroth_weight      <Warning: Zeroth order combo should use D-Vector weight.>
 csaopen            <Error opening skymap !AS; status !SL.>/fao=2
 csaoffset          <Error calling csa_field_offset_values; status !SL.>/fao=1
 csaread            <Error calling csa_read_pixels; status !SL.>/fao=1
 csawrite           <Error calling csa_write_pixels; status !SL.>/fao=1
 csaclose           <Error closing skymap !AS; status !SL.>/fao=2
.end
