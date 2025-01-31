	Integer*4  Function  FIP_READ_COV  ( arch_in, lun_rpt, report, filext,
	1                                    chan, scan, related_data,
	1                                    bin_total, wtd_bin_total, eff_wt,
	1                                    rdisp, idisp, temp, fcc_covar )

C Open FCC_COV_ccss.<fext>, read the FCC covariance matrix, and close the file.
C Convert three vectors of data into original covariance matrix (522x522) form.
C Unpacks all data.
C Unpacking the covariance matrix goes as follows:
C    Records 2-66 contain the Real-Real upper left triangle of matrix 
C       (1-265 x 1-265 = 35245 values)
C    Records 67-190 contain the Real-Imag upper right rectangle of matrix
C       (266-522 x 1-265 = 68105 values)
C    Records 191-251 contain the Imag-Imag lower right triangle of matrix
C       (266-522 x 266-522 = 33153 values)
C
C Returns the full covariance matrix for completeness, even though the lower
C left half is not going to be used.
C
C Larry P. Rosen, Hughes STX, 19 July 1993.
C
C    Input:
C	Character*(*) arch_in              ! input data archive
C	Integer*4     lun_rpt              ! Logical unit number - report file
C	Logical*1     report               ! Flag whether to write report
C	Character*20  filext               ! Input file name extension
C	Character*2   chan                 ! Channel to process (RH,RL,LH,LL)
C	Character*2   scan                 ! Scan mode to process (SS,SF,LS,LF)
C    Output:
C	Record / fcc_cov / related_data    ! Data in the first FCC record
C	Integer*4     bin_total (256)      ! Number of voltage spectra
C	Real*4        wtd_bin_total (256)  ! Weighted total voltage spectra
C	Real*8        eff_wt (256)         ! Effective weight of each volt spec
C	Real*8        rdisp (257,256)      ! Total real spectra dispersion
C	Real*8        idisp (257,256)      ! Total imag spectra dispersion
C	Real*8        temp (8,256)         ! Temp and glitch rate disp. total
C	Real*8        fcc_covar (522,522)  ! Covariance Matrix
C
	Implicit None

C  Include Files

	Include		'($ssdef)'
	Include		'Csdr$Library:CTParams.Inc'

C Passed Parameters:

	Character*(*)	arch_in            ! input data archive
	Integer*4       lun_rpt            ! Logical unit number - report file
	Logical*1       report             ! Flag whether to write report
	Character*20    filext             ! Input file name extension
	Character*2     chan               ! Channel to process (RH,RL,LH,LL)
	Character*2     scan               ! Scan mode to process (SS,SF,LS,LF)
	Dictionary      'fcc_cov'          ! FCC data structure.
	Record / fcc_cov / related_data    ! Data in the first FCC record
	Record / fcc_cov / cov_rec         ! Each record of fcc_cov
	Integer*4       bin_total (256)    ! Number of voltage spectra
	Real*4          wtd_bin_total (256)  ! Weighted total voltage spectra
	Real*8          eff_wt (256)       ! Effective weight of each volt spec
	Real*8          rdisp (257,256)    ! Total real spectra dispersion
	Real*8          idisp (257,256)    ! Total imag spectra dispersion
	Real*8          temp (8,256)       ! Temp and glitch rate disp. total
	Real*8          fcc_covar (522,522)    ! Covariance Matrix

C Function

	Integer*4	Lib$Get_Lun

C  Externals

	Integer*4	CT_Connect_Read
	External	CT_Connect_Read	
	External	fip_normal
	External	fip_abort
	External	fip_lunerr
	External	fip_openerr
	External	fip_readerr
	External	fip_closerr

C Local

	Integer*4	lun_in
	Character*47	infile
	Integer*4	rstat
	Integer*2	ct_stat (20)
	Integer*4	bin_number
	Integer*2	i
	Integer*2	col, row        ! column and row of covariance matrix
	Integer*4	rr_element      ! number of elements put in matrix R-R
	Integer*4	ri_element      ! number of elements put in matrix R-I
	Integer*4	ii_element      ! number of elements put in matrix I-I

C Begin

	FIP_READ_COV = %loc (fip_normal)

C Get a logical unit number

	rstat = Lib$Get_Lun (lun_in)
	If ( rstat .NE. SS$_Normal ) Then
	   FIP_READ_COV = %loc (fip_abort)
	   Call Lib$Signal (fip_lunerr, %val(1), %val(rstat))
	Else

C Open the input covariance matrix data file.

	   infile = arch_in // 'FCC_COV_' // chan // scan // '.' // filext
	   Open ( Unit=lun_in, File=infile, Status='OLD', Iostat=rstat,
	1         Useropen=CT_Connect_Read )

	   If (rstat .NE. 0) Then
	      FIP_READ_COV = %loc (fip_abort)
	      Call Lib$Signal (fip_openerr, %val(2), infile, %val(rstat))
	   Else
	      If (report) Write (lun_rpt, 10) 'opened', infile
   10	      Format (1X, 'Successfully ', A, ' ', A)
	   Endif
	Endif

C If everything ok, read the first record which contains the related data.

	If (FIP_READ_COV .EQ. %loc (fip_normal)) Then
	   Call CT_Read_Arcv (, lun_in, related_data, ct_stat )
	   If (ct_stat(1) .NE. CTP_Normal) Then
	      FIP_READ_COV = %loc (fip_abort)
	      Call Lib$Signal (fip_readerr, %val(2), infile, %val(rstat))
	   Endif
	Endif

C Next, unpack records 2-257 into the dispersion vectors and covariance matrix

	bin_number = 1                         ! = record number - 1
	rr_element = 0
	ri_element = 0
	ii_element = 0
	col = 0
	row = 1
	Do While ( (bin_number .LT. 257) .AND.
	1          (FIP_READ_COV .EQ. %loc (fip_normal) ))

	   Call CT_Read_Arcv (, lun_in, cov_rec, ct_stat )
	   If (ct_stat(1) .NE. CTP_Normal) Then
	      FIP_READ_COV = %loc (fip_abort)
	      Call Lib$Signal (fip_readerr, %val(2), infile, %val(rstat))
	   Else
	      bin_total (bin_number) = cov_rec.bin_total
	      wtd_bin_total (bin_number) = cov_rec.wtd_bin_total
	      eff_wt (bin_number) = cov_rec.eff_wt
	      Do i=1, 257
	         rdisp (i, bin_number) = cov_rec.rdisp (i)
	         idisp (i, bin_number) = cov_rec.idisp (i)
	      Enddo
	      Do i=1,8
	         temp (i, bin_number) = cov_rec.temp (i)
	      Enddo

C Unpack the covariance matrix.
C    Records 2-66 contain the Real-Real upper left triangle of matrix 
C       (1-265 x 1-265 = 35245 values)
C    Records 67-190 contain the Real-Imag upper right rectangle of matrix
C       (266-522 x 1-265 = 68105 values)
C    Records 191-251 contain the Imag-Imag lower right triangle of matrix
C       (266-522 x 266-522 = 33153 values)
C    bin_number = record number - 1

	      If (bin_number .LE. 65) Then
	         Do i=1,550
	            rr_element = rr_element + 1
	            If (rr_element .LE. 35245) Then
	               col = col + 1
	               If (col .GT. 265) Then
	                  row = row + 1
	                  col = row
	               Endif
	               fcc_covar (row, col) = cov_rec.covar (i)
	            Else
	               row = 1      ! Done with RR. set row for RI.  col = 265
	            Endif
	         Enddo
	      Elseif (bin_number .LE. 189) Then
	         Do i=1,550
	            ri_element = ri_element + 1
	            If (ri_element .LE. 68105) Then
	               col = col + 1
	               If (col .GT. 522) Then
	                  row = row + 1
	                  col = 266
	               Endif
	               fcc_covar (row, col) = cov_rec.covar (i)
	            Else
	               col = 265       ! Done with RI. set col = 265, row = 266
	               row = 266
	            Endif
	         Enddo
	      Elseif (bin_number .LE. 250) Then
	         Do i=1,550
	            ii_element = ii_element + 1
	            If (ii_element .LE. 33153) Then
	               col = col + 1
	               If (col .GT. 522) Then
	                  row = row + 1
	                  col = row
	               Endif
	               fcc_covar (row, col) = cov_rec.covar (i)
	            Endif
	         Enddo
	      Endif
	   Endif
	   bin_number = bin_number + 1
	Enddo

C Copy the top right half of the matrix to the bottom left, for completeness.

	If (FIP_READ_COV .EQ. %loc (fip_normal)) Then
	   If (report) Write (lun_rpt, 10) 'read  ', infile
	   Do col = 1, 522
	      Do row = 1, col
	         fcc_covar (col, row) = fcc_covar (row, col)
	      Enddo
	   Enddo
	Endif

C Close the input covariance matrix file.

	If (FIP_READ_COV .EQ. %loc (fip_normal)) Then
	   Call CT_Close_Arcv (, lun_in, ct_stat)
	   If (ct_stat(1) .NE. ctp_normal) Then
	      FIP_READ_COV = %loc (fip_abort)
	      Call Lib$Signal (fip_closerr, %val(2), infile, %val(ct_stat(1)))
	   Elseif (report) Then
	      Write (lun_rpt, 10) 'closed', infile
	      Write (lun_rpt, *)
	   Endif
	EndIf
	Return
	End
