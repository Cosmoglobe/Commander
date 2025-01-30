.title		fcs_msg
.facility	fcs, 705 /prefix=fcs_
.base 0
.severity	success
 normal		<Successful completion.>
.base 1000
.severity	informational
 abort		<Process aborted.>
 zerospec	<Zero spectra averaged.>
 wrong_input	<Skymap file is wrong data type: !AS>/fao=1
 wrong_chan	<Skymap file has wrong channel !AS>/fao=1
 wrong_scan	<Skymap file has wrong scan mode !AS>/fao=1
 multmodel	<Multiple calibration solutions used among input skymaps.>
 galcut         <Pixel !SL has different galaxy cut flags: !SL & !SL>/fao=3
.base 2000
.severity	warning
 parserr	<Error parsing: !AS>/fao=1
 lunerr		<Error getting logic unit number; status !SL.>/fao=1
 openerr	<Error opening script file; status !SL.>/fao=1
 readerr	<Error reading script file; status !SL.>/fao=1
 closerr	<Error closing script file; status !SL.>/fao=1
 noinput	<No skymaps named in script file.>
 query_cat	<Error querying catalog; status !SL.>/fao=1
 invtime	<Jstop time is less than jstart time for skymap !AS.>/fao=1
 ctinit		<Cobetrieve ct_init failure; status !SW.>/fao=1
 zdatatype	<None of the skymaps have the specified datatype: !AS>/fao=1
 zchannel	<None of the skymaps have the specified channel.>
 zscan		<None of the skymaps have the specified scan modes.>
 openrep	<Error opening report file.; status !SL.>/fao=1
 closerep	<Error closing report file.; status !SL.>/fao=1
 renerr		<Error renaming report file.; status !SL.>/fao=1
 csaopen	<Error opening skymap !AS; status !SL.>/fao=2
 csaoffset	<Error calling csa_field_offset_values; status !SL.>/fao=1
 csaread	<Error calling csa_read_pixels; status !SL.>/fao=1
 csawrite	<Error calling csa_write_pixels; status !SL.>/fao=1
 csaclose	<Error closing skymap !AS; status !SL.>/fao=2
.end
