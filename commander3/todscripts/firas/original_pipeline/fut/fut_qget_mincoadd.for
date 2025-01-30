	integer*4 Function Fut_Qget_mincoadd ( Ct_lun, mincoadd_rec)
c----------------------------------------------------------------------
c Program Name : Fut_Qget_mincoadd
c
c Programmer : N. Gonzales/HughesSTX, December 9, 1991, SER 9099
c              
c Program description: This program read Fex_Mincoadd record and formatted
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

        Dictionary 'Fex_Mincoadd'
   	Record     /Fex_Mincoadd/ mincoadd_rec

	external Fut_eof
	external Fut_normal
        external fut_ctqget_err


        Call ct_query_get (, ct_lun,mincoadd_rec,ct_stat)
          If (ct_stat(1) .ne. ctp_normal) then
            If (ct_stat(1) .eq. ctp_endoffile) then
	      fut_qget_mincoadd = %loc(fut_eof)
	    Else
              Call lib$signal (fut_ctqget_err, %val(1), %val(ct_stat(1)))
              fut_qget_mincoadd = %loc(fut_ctqget_err)
            endif
            return
           endif
	Fut_qget_mincoadd = %loc(fut_normal)
        return
        end
