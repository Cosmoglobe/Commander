	integer*4 Function Fut_Qget_facr ( Ct_lun, facr_rec)
c----------------------------------------------------------------------
c Program Name : Fut_Qget_Facr
c
c Programmer : N. Gonzales/STX, August 30, 1991
c              
c Program description: This program read Fex_av_calrs record and formatted
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

        Dictionary 'Fex_av_calrs'
   	Record     /Fex_av_calrs/ facr_rec

	external Fut_eof
	external Fut_normal
        external fut_ctqget_err


        Call ct_query_get (, ct_lun,facr_rec,ct_stat)
          If (ct_stat(1) .ne. ctp_normal) then
            If (ct_stat(1) .eq. ctp_endoffile) then
	      fut_qget_facr = %loc(fut_eof)
	    Else
              Call lib$signal (fut_ctqget_err, %val(1), %val(ct_stat(1)))
              fut_qget_facr = %loc(fut_ctqget_err)
            endif
            return
           endif
	Fut_qget_facr = %loc(fut_normal)
        return
        end
