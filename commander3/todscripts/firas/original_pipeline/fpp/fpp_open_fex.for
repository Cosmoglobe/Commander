C******************************************************************************
	Integer*4 function fpp_open_fex ( arch_ref, gmt_start, gmt_stop,
	1    lun_rpt, report, lun_fake, lun_gain, chan_num)
C
C  Opens reference data files of gain and fakeit data: FEX_FAKEIT.DAT and
C  FEX_GAIN.DAT.  Two opens are necessary for each since we need the 1st
C  record preceding time gmt_start.
C
C  Author:  Larry P. Rosen, STX, 7 March 1991
C******************************************************************************
	implicit none

C  Passed Parameters
	character*15	arch_ref
	character*14	gmt_start
	character*14	gmt_stop
	integer*4	lun_rpt
	logical*1	report
	integer*4	lun_fake
	integer*4	lun_gain
	integer*2	chan_num

C  Include files
	include      'ct$library:ctuser.inc'
	include      '(fut_params)'
	include      '(fpp_msg)'
	include      '($ssdef)'

C  Externals
	integer*4    ct_connect_read
	external     ct_connect_read	

C  Functions
	integer*4    lib$get_lun

C  Local variables
	dictionary 'fex_fakeit'
	record / fex_fakeit / fake_rec

	dictionary 'fex_gain'
	record / fex_gain / gain_rec

	character*13	fake_r_fname / 'FEX_FAKEIT_R/' /
	character*11	gain_r_fname / 'FEX_GAIN_R/' /
	character*13	fake_l_fname / 'FEX_FAKEIT_L/' /
	character*11	gain_l_fname / 'FEX_GAIN_L/' /
	character*64	fake_R_file
	character*62	gain_R_file
	character*64	fake_L_file
	character*62	gain_L_file
	character*64	fake_file
	character*62	gain_file
	character*13	fname
	character*11	gname
	integer*4	lun_f_r, lun_f_l, lun_g_r, lun_g_l
	integer*4	status
	integer*4	iostatus
	integer*2	ct_stat(20)
	character*14	start2			! new time to open fex record
	Logical*1	First_Time /.True./	! if first time, get logicals

C  Set the function status to Normal.

	fpp_open_fex = %loc(fpp_normal)

	if (First_Time) Then

C  Get logical unit numbers for the FOUR files

	   status = lib$get_lun (lun_f_r)
	   if ( status .ne. ss$_normal ) then
	      fpp_open_fex = %loc(fpp_aberr)
	      call lib$signal(fpp_getlunerr, %val(1), %val(status))
	   endif
	   status = lib$get_lun (lun_f_l)
	   if ( status .ne. ss$_normal ) then
	      fpp_open_fex = %loc(fpp_aberr)
	      call lib$signal(fpp_getlunerr, %val(1), %val(status))
	   endif
	   status = lib$get_lun (lun_g_r)
	   if ( status .ne. ss$_normal ) then
	      fpp_open_fex = %loc(fpp_aberr)
	      call lib$signal(fpp_getlunerr, %val(1), %val(status))
	   endif
	   status = lib$get_lun (lun_g_l)
	   if ( status .ne. ss$_normal ) then
	      fpp_open_fex = %loc(fpp_aberr)
	      call lib$signal(fpp_getlunerr, %val(1), %val(status))
	   endif

C  Create fakeit and gain file names.

	   fake_r_file = arch_ref // fake_r_fname // gmt_start // ';' //
	1     gmt_stop // ';'
	   fake_l_file = arch_ref // fake_l_fname // gmt_start // ';' //
	1     gmt_stop // ';'
	   gain_r_file = arch_ref // gain_r_fname // gmt_start // ';' //
	1     gmt_stop // ';'
	   gain_l_file = arch_ref // gain_l_fname // gmt_start // ';' //
	1     gmt_stop // ';'
	   first_time = .False.
  
	endif		! First time

C  If channel number = 1 or 2 then use right side.

	if (chan_num .lt. 3) then
	   lun_fake = lun_f_r
	   lun_gain = lun_g_r
	   fname = fake_r_fname
	   gname = gain_r_fname
	   fake_file = fake_r_file
	   gain_file = gain_r_file
	else
	   lun_fake = lun_f_l
	   lun_gain = lun_g_l
	   fname = fake_l_fname
	   gname = gain_l_fname
	   fake_file = fake_l_file
	   gain_file = gain_l_file
	endif
	if ( fpp_open_fex .eq. %loc(fpp_normal) ) then
	   open (unit=lun_fake, file=fake_file, status='old', iostat=iostatus,
	1     useropen=ct_connect_read)
	   if ( iostatus .ne. 0 ) then
	      fpp_open_fex = %loc(fpp_aberr)
	      call lib$signal (fpp_ctopenerr, %val(1), %val(iostatus))
	   endif
	endif
	if (fpp_open_fex .eq. %loc(fpp_normal)) then
	   open (unit=lun_gain, file=gain_file, status='old', iostat=iostatus,
	1     useropen=ct_connect_read)
	   if ( iostatus .ne. 0 ) then
	      fpp_open_fex = %loc(fpp_aberr)
	      call lib$signal (fpp_ctopenerr, %val(1), %val(iostatus))
	   endif
	endif

C  Read 1st record of fakeit file, extract previous record's gmt, close file,
C  open new file with previous record's gmt.

	if (fpp_open_fex .eq. %loc(fpp_normal)) then
	   call ct_read_arcv (,lun_fake, fake_rec, ct_stat)
	   if (ct_stat(1) .NE. ctp_normal .AND. ct_stat(1) .NE. ctp_Endoffile)
	1  Then
	      call lib$signal (fpp_ctreaderr, %val(1), %val(ct_stat(1)),
	1        'Fex_Fakeit.Dat')
	      fpp_open_fex = %loc(fpp_aberr)
	   else
	      start2 = fake_rec.prev_gmt
	      call ct_close_arcv (,lun_fake, ct_stat)
	      if (ct_stat(1) .ne. ctp_normal) then
	         call lib$signal (fpp_ctcloserr, %val(1), %val(ct_stat(1)))
	         fpp_open_fex = %loc(fpp_aberr)
	      else
	         fake_file = arch_ref // fname // start2 // ';' //
	1           gmt_stop // ';'
	         open (unit=lun_fake, file=fake_file, status='old',
	1           iostat=iostatus, useropen=ct_connect_read)
	         if ( iostatus .ne. 0 ) then
	            fpp_open_fex = %loc(fpp_aberr)
	            call lib$signal (fpp_ctopenerr, %val(1), %val(iostatus))
	         else
	            if (report) write (lun_rpt,fmt=10) fake_file
  10	            format(5X,'Archive file successfully opened: ',/,15x,A)
	         endif
	      endif
	   endif
	endif

C  Read 1st record of gain file, extract previous record's gmt, close file,
C  open new file with previous record's gmt.

	if (fpp_open_fex .eq. %loc(fpp_normal)) then
	   call ct_read_arcv (,lun_gain, gain_rec, ct_stat)
	   if (ct_stat(1) .NE. ctp_normal .AND. ct_stat(1) .NE. ctp_Endoffile)
	1  then
	      call lib$signal (fpp_ctreaderr, %val(1), %val(ct_stat(1)),
	1        'Fex_Gain.Dat')
	      fpp_open_fex = %loc(fpp_aberr)
	   else
	      start2 = gain_rec.prev_gmt
	      call ct_close_arcv (,lun_gain, ct_stat)
	      if (ct_stat(1) .ne. ctp_normal) then
	         call lib$signal (fpp_ctcloserr, %val(1), %val(ct_stat(1)))
	         fpp_open_fex = %loc(fpp_aberr)
	      else
	         gain_file = arch_ref // gname // start2 // ';' //
	1           gmt_stop // ';'
	         open (unit=lun_gain,file=gain_file,status='old',
	1           iostat=iostatus,useropen=ct_connect_read)
	         if ( iostatus .ne. 0 ) then
	            fpp_open_fex = %loc(fpp_aberr)
	            call lib$signal (fpp_ctopenerr, %val(1), %val(iostatus))
	         else
	            if (report) write (lun_rpt,fmt=10) gain_file
	         endif
	      endif
	   endif
	endif
	return
	end
