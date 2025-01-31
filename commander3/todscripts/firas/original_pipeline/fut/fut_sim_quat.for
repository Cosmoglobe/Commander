       integer*4 function fut_sim_quat(Time_Tag,Q,A,Time)
c
c      Module to simulate the COBE attitude determination facility
c
c 	  Shirley M. Read
c 	  STX
c 	  January 20, 1988
c 	  Modified: 	Modified to return Time, the time since JD 2440001.0,
c 			and A, the rotation matrix which describes the S/C
c 			orientation. These are needed for computations in
c 			FUT_Attitude. Also, modified for interface with the
c 		        Fut_Error condition handler. Added error checking
c 		        and calls to Lib$Signal.
c 
       IMPLICIT 	NONE

	include 	'($SSDef)'

	EXTERNAL   	FUT_NORMAL

        integer*4        SYS$AscTim
        integer*4        OTS$Cvt_TU_L
        integer*4        OTS$Cvt_T_F

c	Input/Output Parameters

        integer*4        time_tag(2)  !binary ifg time tag
        real*4           q(4)
        real*4           a(3,3)
        real*8           time         !Time format for attitude

c	Local Variables

        integer*4       status
  
       integer*2   	i
       integer*2        time_len     !Time conversion items

       integer*4        year
       integer*4        mon
       integer*4        day
       integer*4        hr
       integer*4        min

       real*4           sec
       real*4           prev(3)
       real*4           next(3), mag, dot


       character*3      months(12)
       character*32     ctime

       DATA months/'JAN','FEB','MAR','APR','MAY','JUN','JUL','AUG',
     &             'SEP','OCT','NOV','DEC'/

c      Extract year, day, month, etc from the time tag

        status = SYS$AscTim(time_len,ctime,time_tag,0)
	if ( status .ne. SS$_Normal ) then
	  call lib$signal(%val(status))
	  FUT_Sim_Quat = status
	  return
	endif

          mon = 0
          do while(ctime(4:6) .ne. months(mon))
            mon = mon + 1
          end do

          status = OTS$Cvt_TU_L(ctime(8:11),year,%val(4),%val(1))
	  if ( status .ne. SS$_Normal ) then
	    call lib$signal(%val(status))
	    FUT_Sim_Quat = status
	    return
	  endif

          status = OTS$Cvt_TU_L(ctime(1:2),day,%val(4),%val(1))
	  if ( status .ne. SS$_Normal ) then
	    call lib$signal(%val(status))
	    FUT_Sim_Quat = status
	    return
	  endif

          status = OTS$Cvt_TU_L(ctime(13:14),hr,%val(4),%val(1))
	  if ( status .ne. SS$_Normal ) then
	    call lib$signal(%val(status))
	    FUT_Sim_Quat = status
	    return
	  endif

          status = OTS$Cvt_TU_L(ctime(16:17),min,%val(4),%val(1))
	  if ( status .ne. SS$_Normal ) then
	    call lib$signal(%val(status))
	    FUT_Sim_Quat = status
	    return
	  endif

          status = OTS$Cvt_T_F(ctime(19:23),sec, , ,%val(1), )
	  if ( status .ne. SS$_Normal ) then
	    call lib$signal(%val(status))
	    FUT_Sim_Quat = status
	    return
	  endif

c      Convert the time tag to number of seconds since 5/24/68 UT noon

          call ttag(year,mon,day,hr,min,sec,time)

          call aspect(time,a)

          call m2q(a,q)

c      Check direction

          mag = sqrt(q(1)**2 + q(2)**2 + q(3)**2)

          do i = 1,3
            next(i) = q(i)/mag
          end do

          dot = prev(1)*next(1) + prev(2)*next(2) + prev(3)*next(3)
          if(dot .lt. 0.0)then
            do i = 1,4
              q(i) = -q(i)
            end do
            do i = 1,3
              next(i) = -next(i)
            end do
          endif

          do i = 1,3
            prev(i) = next(i)
          end do


	FUT_Sim_Quat = %Loc(FUT_Normal)

       return
       end
