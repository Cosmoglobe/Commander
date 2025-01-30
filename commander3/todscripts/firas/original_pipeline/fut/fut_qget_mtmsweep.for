	integer*4 Function Fut_Qget_mtmsweep ( Ct_lun, mtmsweep_rec)
c----------------------------------------------------------------------
c Program Name : Fut_Qget_Facr
c
c Programmer : N. Gonzales/STX, September 23, 1991
c              
c Program description: This program read Fex_MtmSweep record and formatted
c                      write out the record until end of file.
c
c Include files : ct$libraryt:Ctuser.inc
c
c--------------------------------------------------------------------------

	Implicit None
	Include  'CT$library:ctuser.inc'

        Integer *2 CT_lun
        Integer *2 CT_stat(20)
        Integer *4 istat

        Dictionary 'Fex_MtmSweep'
   	Record     /Fex_MtmSweep/ mtmsweep_rec

	external Fut_eof
	external Fut_normal
        external fut_ctqget_err


        Call ct_query_get (, ct_lun,mtmsweep_rec,ct_stat)
          If (ct_stat(1) .ne. ctp_normal) then
            If (ct_stat(1) .eq. ctp_endoffile) then
	      fut_qget_mtmsweep = %loc(fut_eof)
	    Else
              Call lib$signal (fut_ctqget_err, %val(1), %val(ct_stat(1)))
              fut_qget_mtmsweep = %loc(fut_ctqget_err)
            endif
            return
           endif
	Fut_qget_mtmsweep = %loc(fut_normal)
        return
        end
