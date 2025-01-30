	integer*4 function fga_list_mincoadd (write_lun, mincoadd_rec, mincoadd_file,
	1			 cline, nrec)

C---------------------------------------------------------------------------
C
C	PROGRAM DESCRIPTION:
C	  This program dumps a formatted reference dataset FEX_MINCOADD.
C
C	AUTHOR:
C	  Nilo G. Gonzales/STX, December 9, 1991, Ref. SPR 9099.
C	  
C       
C	CALLING SEQUENCE:
C	  call fga_list_mincoadd (write_lun, mincoadd_rec, mincoadd_file, nrec)
C
C	INPUT PARAMETERS:
C	  write_lun	   i*4      List file logical unit number.
C			            if not supplied assume 6.  Assume the
C			            calling routine has opened this list file
C	  mincoadd_rec     record   mincoadd_rec record.
C	  mincoadd_file	   ch*64    fex_grtrawwt file specification.
C	  nrec		    I*2     The number of records processed in this run.
C
C	OUTPUT PARAMETERS:
C	  none
C
C	INPUT FILES:
C	  none
C
C	OUTPUT FILES:
C	  Listing file or FOR006.
C
C	INCLUDE FILES USED:
C	  none
C
C	SUBROUTINES CALLED:
C	  ct_binary_to_gmt
C
C	ERROR HANDLING:
C	  none
C
C---------------------------------------------------------------------------
C Changes:
C---------------------------------------------------------------------------

	implicit	none

	dictionary		'fex_mincoadd'
	record  /fex_mincoadd/	mincoadd_rec
	integer         *4      min_ifg_coadd(4)
	integer		*2	nrec, i
	integer		*4	write_lun
	integer		*4	iostatus
        integer         *4      status
        integer         *4      success/1/, error/2/
        logical         *1	first/.true./
	character       *64     mincoadd_file		
	character       *80	cline
	character	*1	line(130)/ 130 * '-'/


c
c Begin the dump.
c
        fga_list_mincoadd = success

	if (write_lun .eq. 0) write_lun = 6
	if ( first) then
	  write (write_lun, 20, iostat=iostatus) mincoadd_file,cline,line
          first = .false.
        endif
 20     format (// 1x, 'FEX_MINCOADD FILE : ',a/,x,a80/,130a/)

	if (iostatus .ne. 0) then
	  print *,' write error!!! unit = ', write_lun 
          fga_list_mincoadd = error          
	end if
c
c Header information.
c
	write (write_lun, 40, iostat=iostatus)
	1     (mincoadd_rec.ct_head.gmt(i:i),i=1,14),
	1     mincoadd_rec.ct_head.time(2),
	1     mincoadd_rec.ct_head.time(1),
	1     mincoadd_rec.ct_head.space_time,
	1     mincoadd_rec.ct_head.mjr_frm_no,
	1     mincoadd_rec.ct_head.orbit,
	1     mincoadd_rec.ct_head.dataset_id

40	format (/' <Start GMT> : ', 2a1, '-', 3a1, '-', 3(2a1, '-'), 3a1, 7x,
	1	' <Binary TIME> : ', z8.8, 1x, z8.8, 7x, /,
	1       ' <PB5> : ', 6(z2.2, 1x) /
	1       ' Major Frame Number : ', i6, 8X,
	1       ' Orbit Number : ', i12, 12X,
	1	' Data Set Id  : ', i10/,130a/)

	if (iostatus .ne. 0) then
	  print *,' write error!!! unit = ', write_lun 
          fga_list_mincoadd = error          
	end if

	do i=1,4
	   min_ifg_coadd(i) = mincoadd_rec.min_ifg_coadd(i)
	enddo

	write (write_lun, 50, iostat=iostatus) min_ifg_coadd

50 	Format (1X,'**** Minimum_IFG_Coadd **** '//,
	1       5x,' Channel RH:',i3,10x,'RL:',i3 /
	1       5x,'         LH:',i3,10x,'LL:',i3)

	if (iostatus .ne. 0) then
	  print *,' write error!!! unit = ', write_lun 
          fga_list_mincoadd = error          
	end if

        return
	end
