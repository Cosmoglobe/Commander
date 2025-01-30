	integer*4 function fga_list_facr (write_lun, facr_rec, facr_file, cline,
	1				 nrec)

C---------------------------------------------------------------------------
C
C	PROGRAM DESCRIPTION:
C	  This program dumps a formatted reference dataset FEX_AV_CALRS.
C
C	AUTHOR:
C	  Nilo G. Gonzales/STX, September 4, 1991
C	  
C       
C	CALLING SEQUENCE:
C	  CALL FGA_LIST_FACR (WRITE_LUN, FACR_REC, FACR_FILE, NREC)
C
C	INPUT PARAMETERS:
C	  WRITE_LUN	I*4     List file logical unit number.
C			        if not supplied assume 6.  Assume the
C			        calling routine has opened this list file
C	  FACR_REC      RECORD	FACR_REC record.
C	  FACR_FILE	CH*64	FEX_AV_CALRS file specification.
C	  NREC			I*2	The number of records processed in this run.
C
C	OUTPUT PARAMETERS:
C	  NONE
C
C	INPUT FILES:
C	  NONE
C
C	OUTPUT FILES:
C	  Listing file or FOR006.
C
C	INCLUDE FILES USED:
C	  NONE
C
C	SUBROUTINES CALLED:
C	  CT_BINARY_TO_GMT
C
C	ERROR HANDLING:
C	  NONE
C
C---------------------------------------------------------------------------
C Changes:
C---------------------------------------------------------------------------

	implicit	none

	integer		*2	nrec
	integer		*4	write_lun
	integer		*4	iostatus
        integer         *4      status
        integer         *4      success/1/, error/2/
        logical         *1	first/.true./

	character       *64     facr_file		
	character       *80	cline
	character	*1	line(130)/ 130 * '-'/

	dictionary		'fex_av_calrs'
	record /fex_av_calrs/	facr_rec


C
C Begin the dump.
C
        fga_list_facr = success

	if (write_lun .eq. 0) write_lun = 6
	if ( first) then
	  write (write_lun, 20, iostat=iostatus) facr_file,cline,line
          first = .false.
        endif
 20     format (// 1x, 'fex_av_calrs file : ',a/,x,a80/,130a/)

	if (iostatus .ne. 0) then
	  print *,' write error!!! unit = ', write_lun 
          fga_list_facr = error          
	end if

	write ( write_lun,30) nrec
  30   format(/,t32,' Rec # : ',i4/)
C		
C Calibrator Resistors dumps information.

	Write (write_lun,40, iostat=iostatus) facr_rec.ct_head.gmt,
	1       FACR_REC.data_stop, 
	1       FACR_REC.prev_data_start,
	1	FACR_REC.Ave_Period, 
	1	FACR_REC.Num_good_record, 
	1       FACR_REC.Num_bad_record

  40	Format (1x,'START_TIME = ',A14,20X,'STOP_TIME = ',A14,/
	1	10x, 'PREV_GMT =',A14,/
	1	15X,'Average period = ',F6.3,/
	1	15X,'Number of Good Records = ',I5,/
	1       15x,'Number of Bad  Records = ',I5,/)

	Write (write_lun,50,iostat=iostatus) FACR_REC.Calres_Ave_A_Lo_dev(1), 
	1       FACR_REC.Calres_Ave_A_Lo_dev(2),
	1       FACR_REC.Calres_Ave_A_Lo_dev(3),
	1       FACR_REC.Calres_Ave_A_Lo_dev(4),
	1       FACR_REC.Calres_Ave_A_Hi_dev(1),
	1       FACR_REC.Calres_Ave_A_Hi_dev(2),
	1       FACR_REC.Calres_Ave_A_Hi_dev(3),
	1       FACR_REC.Calres_Ave_A_Hi_dev(4),
	1       FACR_REC.Calres_Ave_B_Lo_dev(1),
	1       FACR_REC.Calres_Ave_B_Lo_dev(2),
	1       FACR_REC.Calres_Ave_B_Lo_dev(3),
	1       FACR_REC.Calres_Ave_B_Lo_dev(4),
	1       FACR_REC.Calres_Ave_B_Hi_dev(1),
	1       FACR_REC.Calres_Ave_B_Hi_dev(2),
	1       FACR_REC.Calres_Ave_B_Hi_dev(3),
	1       FACR_REC.Calres_Ave_B_Hi_dev(4)

  50	Format (/10x,'Calres_Ave_A_Lo_Dev_1  = ',F6.3,/ 
	1       10x,'                Dev_2  = ',F6.3,/
	1       10x,'                Dev_3  = ',F6.3,/
	1       10x,'                Dev_4  = ',F6.3,/
	1       10x,'Calres_Ave_A_Hi_Dev_1  = ',F6.3,/
	1       10x,'                Dev_2  = ',F6.3,/
	1       10x,'                Dev_3  = ',F6.3,/
	1       10x,'                Dev_4  = ',F6.3,/
	1       10x,'Calres_Ave_B_Lo_Dev_1  = ',F6.3,/
	1       10x,'                Dev_2  = ',F6.3,/
	1       10x,'                Dev_3  = ',F6.3,/
	1       10x,'                Dev_4  = ',F6.3,/
	1       10x,'Calres_Ave_B_Hi_Dev_1  = ',F6.3,/
	1       10x,'                Dev_2  = ',F6.3,/
	1       10x,'                Dev_3  = ',F6.3,/
	1       10x,'                Dev_4  = ',F6.3,/)

        Write (write_lun,60,iostat=iostatus) FACR_REC.Calres_Ave_A_Lo(1), 
	1	FACR_REC.calres_bad_Pnts_a_lo(1),
	1       FACR_REC.Calres_Ave_A_Lo(2),
	1	FACR_REC.calres_bad_Pnts_a_lo(2),
	1       FACR_REC.Calres_Ave_A_Lo(3),
	1	FACR_REC.calres_bad_Pnts_a_lo(3),
	1       FACR_REC.Calres_Ave_A_Lo(4),
	1	FACR_REC.calres_bad_Pnts_a_lo(4)

  60	Format (10x,'Calres_Ave_A_Lo_1 = ',f12.3, ' Number Bad Points = ', I6,/ 
	1       10x,'           A_Lo_2 = ',f12.3,' Number Bad Points = ', I6,/
	1       10x,'           A_Lo_3 = ',f12.3,' Number Bad Points = ', I6,/
	1       10x,'           A_Lo_4 = ',f12.3,' Number Bad Points = ', I6,/)

	Write (write_lun,70,iostat=iostatus) FACR_REC.Calres_Ave_A_Hi(1), 
	1	FACR_REC.calres_bad_Pnts_a_Hi(1),
	1       FACR_REC.Calres_Ave_A_Hi(2),
	1	FACR_REC.calres_bad_Pnts_a_Hi(2),
	1       FACR_REC.Calres_Ave_A_Hi(3),
	1	FACR_REC.calres_bad_Pnts_a_hi(3),
	1       FACR_REC.Calres_Ave_A_Hi(4),
	1	FACR_REC.calres_bad_Pnts_a_Hi(4)

  70	Format (10x,'Calres_Ave_A_H1_1 = ',f12.3,' Number Bad Points = ', I6,/ 
	1       10x,'           A_Hi_2 = ',f12.3,' Number Bad Points = ', I6,/
	1       10x,'           A_Hi_3 = ',f12.3,' Number Bad Points = ', I6,/
	1       10x,'           A_Hi_4 = ',f12.3,' Number Bad Points = ', I6,/)

	Write (write_lun,80,iostat=iostatus) FACR_REC.Calres_Ave_B_Lo(1), 
	1	FACR_REC.calres_bad_Pnts_b_lo(1),
	1       FACR_REC.Calres_Ave_B_Lo(2),
	1	FACR_REC.calres_bad_Pnts_b_lo(2),
	1       FACR_REC.Calres_Ave_B_Lo(3),
	1	FACR_REC.calres_bad_Pnts_b_lo(3),
	1       FACR_REC.Calres_Ave_B_Lo(4),
	1	FACR_REC.calres_bad_Pnts_b_lo(4)

  80	Format (10x,'Calres_Ave_B_Lo_1 = ',f12.3,' Number Bad Points = ', I6,/
	1       10x,'           B_Lo_2 = ',f12.3,' Number Bad Points = ', I6,/
	1       10x,'           B_Lo_3 = ',f12.3,' Number Bad Points = ', I6,/
	1       10x,'           B_Lo_4 = ',f12.3,' Number Bad Points = ', I6,/)

	Write (write_lun,90,iostat=iostatus) FACR_REC.Calres_Ave_B_Hi(1),
	1	FACR_REC.calres_bad_Pnts_b_hi(1),
	1       FACR_REC.Calres_Ave_B_Hi(2),
	1	FACR_REC.calres_bad_Pnts_b_hi(2),
	1       FACR_REC.Calres_Ave_B_Hi(3),
	1	FACR_REC.calres_bad_Pnts_b_hi(3),
	1       FACR_REC.Calres_Ave_B_Hi(4),
	1	FACR_REC.calres_bad_Pnts_b_Hi(4)

  90	Format (10x,'Calres_Ave_B_Hi_1 = ',f12.3,' Number Bad Points = ', I6,/
	1       10x,'           B_Hi_2 = ',f12.3,' Number Bad Points = ', I6,/
	1       10x,'           B_Hi_3 = ',f12.3,' Number Bad Points = ', I6,/
	1       10x,'           B_Hi_4 = ',f12.3,' Number Bad Points = ', I6)

	if (iostatus .ne. 0) then
	  print *,' write error!!! unit = ', write_lun 
          fga_list_facr = error          
	end if

        return
	end
