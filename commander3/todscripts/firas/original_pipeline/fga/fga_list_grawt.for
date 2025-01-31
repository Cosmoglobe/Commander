	integer*4 function fga_list_grawt (write_lun, grawt_rec, grawt_file,
	1			 cline, nrec)

C---------------------------------------------------------------------------
C
C	PROGRAM DESCRIPTION:
C	  This program dumps a formatted reference dataset FEX_GRTRAWWT.
C
C	AUTHOR:
C	  Nilo G. Gonzales/STX, September 4, 1991
C	  
C       
C	CALLING SEQUENCE:
C	  call fga_list_grawt (write_lun, grawt_rec, grawt_file, nrec)
C
C	INPUT PARAMETERS:
C	  write_lun	i*4     List file logical unit number.
C			        if not supplied assume 6.  Assume the
C			        calling routine has opened this list file
C	  grawt_rec     record	grawt_rec record.
C	  grawt_file	ch*64	fex_grtrawwt file specification.
C	  nrec		I*2     The number of records processed in this run.
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

	dictionary		'fex_grtrawwt'
	record  /fex_grtrawwt/	grawt_rec

	integer		*2	nrec, i
	integer		*4	write_lun
	integer		*4	iostatus
        integer         *4      status
        integer         *4      success/1/, error/2/
        logical         *1	first/.true./
	character       *64     grawt_file		
	character       *80	cline
	character	*1	line(130)/ 130 * '-'/

	character       *30     text(3) /'NUM_GRT',
	1                       'GRT_A_WEIGHT','GRT_B_WEIGHT'/

c
c Begin the dump.
c
        fga_list_grawt = success

	if (write_lun .eq. 0) write_lun = 6
	if ( first) then
	  write (write_lun, 20, iostat=iostatus) grawt_file,cline,line
          first = .false.
        endif
 20     format (// 1x, 'fex_grtrawwt file : ',a/,x,a80/,130a/)

	if (iostatus .ne. 0) then
	  print *,' write error!!! unit = ', write_lun 
          fga_list_grawt = error          
	end if

	write ( write_lun,30) nrec
  30    format(/,t32,' Rec # : ',i4/)
c
c Header information.
c
	write (write_lun, 40, iostat=iostatus)
	1     (grawt_rec.ct_head.gmt(i:i),i=1,14),
	1     grawt_rec.ct_head.time(2),
	1     grawt_rec.ct_head.time(1),
	1     grawt_rec.ct_head.space_time,
	1     grawt_rec.ct_head.mjr_frm_no,
	1     grawt_rec.ct_head.orbit,
	1     grawt_rec.ct_head.dataset_id

40	format (/' <Start GMT> : ', 2a1, '-', 3a1, '-', 3(2a1, '-'), 3a1, 7x,
	1	' <Binary TIME> : ', z8.8, 1x, z8.8, 7x, /,
	1       ' <PB5> : ', 6(z2.2, 1x) /
	1       ' Major Frame Number : ', i6, 8X,
	1       ' Orbit Number : ', i12, 12X,
	1	' Data Set Id  : ', i10/,130a/)

	write (write_lun, 50, iostat=iostatus)
	1       text(1),grawt_rec.num_grt,
	1       text(2),(grawt_rec.grt_a_weight(i),i=1,16),
	1       text(3),(grawt_rec.grt_b_weight(i),i=1,16)

50 	Format (1X,'FEX_GRTRAWWT Data Record : '//,
	1       1x,a,' : ',i12/,
	1       1x,a,' : ',8f11.2/,34x,8f11.2//,
	1       1x,a,' : ',8f11.2/,34x,8f11.2/)

	if (iostatus .ne. 0) then
	  print *,' write error!!! unit = ', write_lun 
          fga_list_grawt = error          
	end if

        return
	end
