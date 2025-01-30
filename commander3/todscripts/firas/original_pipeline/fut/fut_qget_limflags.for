	Integer*4 Function Fut_Qget_Limflags ( Ct_lun, limflags)
c----------------------------------------------------------------------
c Program Name : Fut_Qget_limflags
c
c Programmer : Quoc C Chung/STX
c              Oct 1989 SPR 4750 Formatted dumped out Firas Reference files.
c              
c Program description: This program read Fex_limflags record and formatted
c                      write out the record until end of file.
c
c Include files : Ct$library:Ctuser.inc
c
c--------------------------------------------------------------------------

	Implicit None
	Include  'CT$library:ctuser.inc'

        Integer *2 CT_lun
        Integer *2 CT_stat(20)
        Integer *4 istat

        Dictionary 'Fex_Limflags'
   	Record     /Fex_Limflags/ Limflags

	external Fut_eof
	external Fut_normal
        external fut_ctqget_err

        Call ct_query_get (, ct_lun,limflags,ct_stat)
          If (ct_stat(1) .ne. ctp_normal) then
            If (ct_stat(1) .eq. ctp_endoffile) then
	      fut_qget_limflags = %loc(fut_eof)
	    Else
              Call lib$signal (fut_ctqget_err, %val(1), %val(ct_stat(1)))
              fut_qget_limflags = %loc(fut_ctqget_err)
            endif
            return
           endif
	Fut_qget_limflags = %loc(fut_normal)
        return
        end
