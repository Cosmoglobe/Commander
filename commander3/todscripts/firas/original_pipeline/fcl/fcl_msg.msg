.facility fcl,747 /prefix=fcl_

.base 0
.severity success

 normal    <Process FCL succeeded>

.base 1000
.severity warning

 failure   <Process FCL failed>
 parseerr  <Error parsing command line>
 ctinit    <Error from Cobetrieve initialize; Status !SW>/fao=1
 cfseqopen <Error opening sequential config files; Status !SL>/fao=1
 cfseqget  <Error getting sequential config files; Jstart !AS; Status !SL>/fao=2
 cfseqclose<Error closing sequential config files; Status !SL>/fao=1
 cfdiropen <Error opening direct config files; Status !SL>/fao=1
 cfdirget  <Error getting direct config files; Jstart !AS; Status !SL>/fao=2
 cfdirread <Error reading direct config file !AS; Status !SL>/fao=2
 cfdirclose<Error closing direct config files; Status !SL>/fao=1
 ctopen    <Error opening file !AS; Status !SL>/fao=2
 ctread    <Error reading file !AS; Status !SW>/fao=2
 ctwrite   <Error writing file !AS; Status !SW>/fao=2
 ctclose   <Error closing file !AS; Status !SW>/fao=2
 repopen   <Error opening report file !AS; Status !SL>/fao=2
 repwrite  <Error writing report file !AS; Status !SL>/fao=2

.end
