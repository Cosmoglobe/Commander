
	subroutine fep_grtcoeffs (polycoefs,config_Lun,Con_lun,Btime,
	1	                  rcode )
c
c  This version editted by Rich Isaacman 19 May 1989 to use double
c  matrix inversion
c--------------------------------------------------------------------
C  Change Log:
C  
C       Spr 3253, Dwellplot has to fetch the reference data from the 
C                 reference data archive.
C                 H. Wang, STX, Mar. 2, 1990
C--------------------------------------------------------------------
        Implicit         None
        Dictionary       'FEX_CALRES'
        Record/FEX_CALRES/ CALRES
        Include  '(CCT_GET_CONFIG)'
        Record/config_status/ stat
	real*4 polycoefs(3,4), resis_cal(4,4), tol
	real*8 rhs(4), q(4,3), coeffs(3)
	integer*4 kal_counts(4,4), Config_lun, Con_lun
        Integer*4 Btime(2), status, iset, ncal
        Integer*4 index1, index2, rcode
        Integer*4 size/256/, tlm, kbasis, m, n, ier
        logical*1 new_cnts, new_cal
        Integer*4 CCT_GET_CONFIG_TOD, CCT_GET_CONFIG_IDX_TOD
        Integer*4 Fep_Get_Calres
        External  Fep_getconfigerr
C
C  
        rcode = 0
        Status = CCT_GET_CONFIG_TOD(BTIME,1,size,CON_LUN,Index1,calres,
	1	     NEW_CNTS,STAT)
        If (.not. status) then
          Call lib$signal(Fep_getconfigErr,%val(1),%val(status))
          rcode = 1
          return
        Endif  
        Status = CCT_GET_CONFIG_Idx_tod(btime,1,config_lun,index2,new_cal,stat)
        If (.not. status) then
          Call lib$signal(Fep_getconfigErr,%val(1),%val(status))
          rcode = 1
          return
        Endif  
        If (new_cal) then
          Status =  fep_get_calres (config_lun,resis_cal)
          If (.not. Status) Then
            rcode = 1
            return
          Endif  
        Endif
	do iset=1,4

	   do ncal=1,4
	      tlm = calres.kal_counts(ncal,iset)
	      q(ncal,1) = tlm
	      q(ncal,2) = 1.
	      q(ncal,3) = tlm * tlm
	      rhs(ncal) = resis_cal(ncal,iset)
	   enddo

	   tol = 0.
	   kbasis = 3
	   call dlsarg (kbasis, q, 4, rhs, 1, coeffs)

	   if (ier .ne. 0) then
	      do m=1,3
	         coeffs(m) = 0.
	      enddo
	   else
	      do n=1,3
	         polycoefs(n,iset) = coeffs(n)
	      enddo
	   endif

	enddo
	return
	end
