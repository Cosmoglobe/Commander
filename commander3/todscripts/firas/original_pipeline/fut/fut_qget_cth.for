	Integer*4 Function Fut_Qget_cth ( Ct_lun, cth_rec)
c----------------------------------------------------------------------
c Program Name : Fut_Qget_cth
c
c Programmer : Quoc C Chung/STX
c              NOV 15, 1989 SPR 4750 Formatted dumped out Firas Reference files.
c              
c Program description: This program read Fex_cth record and formatted
c                      write out the record until end of file.
c
c Include files : ct$library:Ctuser.inc
c
c--------------------------------------------------------------------------

	Implicit None
	Include  'CT$library:ctuser.inc'

        Integer *2 CT_lun
        Integer *2 CT_stat(20)
        Integer *4 istat

	external Fut_eof
	external Fut_normal
        external fut_ctqget_err

        Dictionary 'Fex_cth'
   	Record     /Fex_cth/ cth_rec


        Call ct_query_get (, ct_lun,cth_rec,ct_stat)
          If (ct_stat(1) .ne. ctp_normal) then
            If (ct_stat(1) .eq. ctp_endoffile) then
	      fut_qget_cth = %loc(fut_eof)
	    Else
              Call lib$signal (fut_ctqget_err, %val(1), %val(ct_stat(1)))
              fut_qget_cth = %loc(fut_ctqget_err)
            endif
            return
           endif
	Fut_qget_cth = %loc(fut_normal)
        return
        end
