.facility fic,745 /prefix=fic_

.base 0
.severity success

 normal    <Process FIC succeeded>

.base 1000
.severity warning

 failure   <Process FIC failed>
 ctinit    <Error from Cobetrieve initialize; Status !SW>/fao=1
 nocatrecs <No catalog records found in archive !AS for dataset !AS>/fao=2
 invtime   <Jstop time !AS is less than jstart time !AS>/fao=2
 tempzero  <Primary template is zero; Channel !AS; Coadd GMT !AS>/fao=2
 nosectemp <Bad peak, sec template not formed; Channel !AS; Coadd GMT !AS>/fao=2
 primrms   <Primary template rms is zero; Channel !AS; Coadd GMT !AS>/fao=2
 secrms    <Secondary template rms is zero; Channel !AS; Coadd GMT !AS>/fao=2
 cfseqopen <Error opening sequential config files; Status !SL>/fao=1
 cfseqget  <Error getting sequential config files; Jstart !AS; Status !SL>/fao=2
 cfseqclose<Error closing sequential config files; Status !SL>/fao=1
 cfdiropen <Error opening direct config files; Status !SL>/fao=1
 cfdirget  <Error getting direct config files; Jstart !AS; Status !SL>/fao=2
 cfdirread <Error reading direct config file !AS; Status !SL>/fao=2
 cfdirclose<Error closing direct config files; Status !SL>/fao=1
 csaopen   <Error opening skymap !AS; Status !SL>/fao=2
 csaread   <Error reading skymap !AS; Status !SL>/fao=2
 csawrite  <Error writing skymap !AS; Status !SL>/fao=2
 csaclose  <Error closing skymap !AS; Status !SL>/fao=2
 ctopen    <Error opening file !AS; Status !SL>/fao=2
 ctread    <Error reading file !AS; Status !SW>/fao=2
 ctquery   <Error querying file !AS; Status !SW>/fao=2
 ctkeyread <Error key reading file !AS at time !AS; Status !SW>/fao=3
 ctwrite   <Error writing file !AS; Status !SW>/fao=2
 ctclose   <Error closing file !AS; Status !SW>/fao=2
 rseopen   <Error opening RSE file !AS; Status !SL>/fao=2
 rseread   <Error reading RSE file !AS; Status !SL>/fao=2
 rseclose  <Error closing RSE file !AS; Status !SL>/fao=2
 repopen   <Error opening report file !AS; Status !SL>/fao=2
 repwrite  <Error writing report file !AS; Status !SL>/fao=2
 repclose  <Error closing report file !AS; Status !SL>/fao=2

.end
