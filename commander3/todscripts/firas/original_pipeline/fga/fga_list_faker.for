	integer*4 function fga_list_faker (write_lun, faker_rec, faker_file,
	1				    cline, nrec)

C---------------------------------------------------------------------------
C
C	PROGRAM DESCRIPTION:
C	  This program dumps a formatted reference dataset FEX_FAKEIT_R.
C
C	AUTHOR:
C	  Nilo G. Gonzales/STX, September 3, 1991.
C	  
C       
C	CALLING SEQUENCE:
C	  call fga_list_faker (write_lun, faker_rec, faker_file, nrec)
C
C	INPUT PARAMETERS:
C	  write_lun	i*4     List file logical unit number.
C			        if not supplied assume 6.  Assume the
C			        calling routine has opened this list file
C	  faker_rec     record	faker_rec record.
C	  faker_file	ch*64	fex_fakeit_r file specification.
C	  nrec			i*2 The number of records processed in this run.
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

	integer		*2	nrec
	integer		*4	write_lun
	integer		*4	iostatus
        integer         *4      status
        integer         *4      success/1/, error/2/
        logical         *1	first/.true./

	character       *64     faker_file	
	character       *80	cline
	character	*1	line(130)/ 130 * '-'/

	dictionary		'fex_fakeit'
	record /fex_fakeit/	faker_rec


C
C Begin the dump.
C
        fga_list_faker = success

	if (write_lun .eq. 0) write_lun = 6
	if ( first) then
	  write (write_lun, 10, iostat=iostatus) faker_file,cline,line
          first = .false.
        endif
 10   format(	// 1x, 'FEX_Fakeit_R file : ',a/,x,a80/,130a/)

	if (iostatus .ne. 0) then
	  print *,' write error!!! unit = ', write_lun 
          fga_list_faker = error          
	end if

	write ( write_lun,20) nrec
  20   format(/,t32,' Rec # : ',i4/)
c	
c Dump reference file fex_fakeit_r.
c
	write (write_lun,40,iostat=iostatus) Faker_rec.CT_head.GMT,
	1      Faker_rec.Prev_Gmt,
	1      Faker_rec.Start_GMT,
	1      Faker_rec.Stop_GMT

  40	Format (1x,'CT_HEAD.GMT = ',a14,20x,'Previous GMT = ',a14,/,
	1	5x,'Start_GMT = ',a14,10x,'Stop_GMT = ',a14,/)

	write (write_lun,50,iostat=iostatus) Faker_rec.Fakeit(1),
	1      Faker_rec.Fakeit(2),
	1      Faker_rec.Fakeit_Change(1),
	1      Faker_rec.Fakeit_Change(2),
	1      Faker_rec.TLM_Bad_Quality,
	1      Faker_rec.Fakeit_Data_Gap,
	1      Faker_rec.End_of_Data

  50 	Format (5x,'FakeIt Statuses :',2x,'RH:',i3,4x,'RL:',i3,/
	1       5x,'FakeIt Change   :',2x,'RH:',i3,4x,'RL:',i3,/
	1       5x,'TLM Quality: ',i2,3x,'GAP: ',i2,3x,'EOF: ',i2,/)

	if (iostatus .ne. 0) then
	  print *,' write error!!! unit = ', write_lun 
          fga_list_faker = error          
	end if
        return
	end
