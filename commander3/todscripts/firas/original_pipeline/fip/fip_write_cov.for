	Integer*4  Function  FIP_WRITE_COV  ( lun_rpt, report, related_data,
	1                                     bin_total, wtd_bin_total, eff_wt,
	1                                     rdisp, idisp, temp, fcc_covar,
	1                                     index_low, index_high, arch_out,
	1                                     chan, scan, filext, delta_nu,
	1                                     nu_zero )

C Write out the converted covariance matrix and related data in the FIP_COV
C RDL structure.
C
C Larry P. Rosen, Hughes STX, 26 July 1993.
C
C    Input:
C	Integer*4     lun_rpt              ! Logical unit number - report file
C	Logical*1     report               ! Flag whether to write report
C	Dictionary    'fcc_cov'            ! FCC data structure.
C	Record / fcc_cov / related_data    ! Data in the first FCC record
C	Integer*4     bin_total (256)      ! Number of voltage spectra
C	Real*4        wtd_bin_total (256)  ! Weighted total voltage spectra
C	Real*8        eff_wt (256)         ! Effective weight of each volt spec
C	Real*8        rdisp (257,256)      ! Total real spectra dispersion
C	Real*8        idisp (257,256)      ! Total imag spectra dispersion
C	Real*8        temp (8,256)         ! Temp and glitch rate disp. total
C	Real*8        fcc_covar (522,522)  ! FCC Covariance Matrix
C	Integer*2     index_low            ! Low frequency cut-off array index
C	Integer*2     index_high           ! High frequency cut-off array index
C	Character*15  arch_out             ! Output data archive
C	Character*2   chan                 ! Channel to process (RH,RL,LH,LL)
C	Character*2   scan                 ! Scan mode to process (SS,SF,LS,LF)
C	Character*20  filext               ! Input file name extension
C	Real*4        delta_nu             ! Optical frequency interval, GHz.
C	Real*4        nu_zero              ! Optical frequency of initial data
C                                          !    point, in GHz.
C    Output: writes data to output file FIP_COV_ccss.filext
C
C
C PDL:
C    1. Create the 368 x 368 real*8 matrix COVAR and transcribe converted
C    elements of FCC_COVAR into it.
C
C    2. Pack the COVAR matrix into the 3 vector form used for storage in
C    records.
C
C    3. Open output file FIP_COV_ccss.<fext>.
C
C    4. Copy the following quantities from the input FCC_COV file to the output
C    records with no changes:  MODEL_LABEL, COV_LABEL, CHANSCAN, NUM_IFGS,
C    NUM_COADDS, DEG_FREEDOM, BOL_AVG, VOLT_AVG, S0, TAU, TBOL, TEMP,
C    BIN_TOTALS, WTD_BIN_TOTALS.
C
C    5. Convert CBIAS_AVG from counts to volts by divide by 25.5.
C
C    6. COVAR is written to the output file according to RDL.
C
C    7.  Close the FIP_COV_ccss output file.
C

	Implicit None

C  Include Files

	Include		'($ssdef)'
	Include		'CT$Library:CTUser.Inc'

C Passed parameters

	Integer*4     lun_rpt              ! Logical unit number - report file
	Logical*1     report               ! Flag whether to write report
	Dictionary    'fcc_cov'            ! FCC data structure.
	Record / fcc_cov / related_data    ! Data in the first FCC record
	Integer*4     bin_total (256)      ! Number of voltage spectra
	Real*4        wtd_bin_total (256)  ! Weighted total voltage spectra
	Real*8        eff_wt (256)         ! Effective weight of each volt spec
	Real*8        rdisp (257,256)      ! Total real spectra dispersion
	Real*8        idisp (257,256)      ! Total imag spectra dispersion
	Real*8        temp (8,256)         ! Temp and glitch rate disp. total
	Real*8        fcc_covar (522,522)  ! FCC Covariance Matrix
	Integer*2     index_low            ! Low frequency cut-off array index
	Integer*2     index_high           ! High frequency cut-off array index
	Character*15  arch_out             ! Output data archive
	Character*2   chan                 ! Channel to process (RH,RL,LH,LL)
	Character*2   scan                 ! Scan mode to process (SS,SF,LS,LF)
	Character*20  filext               ! Input file name extension
	Real*4        delta_nu             ! Optical frequency interval, GHz.
	Real*4        nu_zero              ! Optical frequency of initial data
                                           !    point, in GHz.

C Functions

	Integer*4	Lib$Get_Lun
	Integer*4	FIP_TRANSCRIBE_COV
	Integer*4	FIP_PACK_COV

C  Externals

	Integer*4	CT_Connect_Write
	External	CT_Connect_Write
	External	fip_normal
	External	fip_abort
	External	fip_lunerr
	External	fip_openerr
	External	fip_writerr
	External	fip_closerr

C Local

	Real*4		cbias_volts	! Cbias_avg converted to volts
	Integer*4	lun_out
	Character*47	outfile
	Integer*4	rstat
	Integer*2	ct_stat (20)
	Dictionary      'fip_cov'       ! FIP_COV data structure.
	Record / fip_cov / fip_rec      ! One FIP_COV record
	Integer*4	rec             ! record number counter
	Integer*2	freq            ! counter for frequencies
	Integer*2	i               ! counter
	Real*8		covar (368,368) ! FIP converted covariance matrix
	Real*8		real_real (17766)	! Vector of covar matrix
	Real*8		real_imag (33840)	! Vector of covar matrix
	Real*8		imag_imag (16290)	! Vector of covar matrix
	Integer*4	rr_element      ! number of elements put in matrix R-R
	Integer*4	ri_element      ! number of elements put in matrix R-I
	Integer*4	ii_element      ! number of elements put in matrix I-I

C Begin

	FIP_WRITE_COV = %loc (fip_normal)

	cbias_volts = related_data.cbias_avg / 25.5

C Write info about the covariance matrix to the report file.

	If (report) Then
	   Write (lun_rpt, 10) nu_zero, delta_nu, related_data.bol_avg,
	1                      related_data.cbias_avg, cbias_volts,
	1                      related_data.volt_avg,
	1                      Real (related_data.s0), Real (related_data.tau),
	1                      related_data.tbol

   10	   Format (1X, 'Optical Frequency of First Point:', 7X, F8.4, ' GHz',/,
	1          1X, 'Optical Frequency Interval:', 13X, F8.4, ' GHz',/,
	1          1X, 'Average Bolometer Temp:', 18X, F7.4, ' deg K',/,
	1          1X, 'Average Commanded Bias:', 18X, F6.1,'  counts', /,
	1          1X, 'Average Commanded Bias:', 18X, F6.2,'  volts', /,
	1          1X, 'Average Readout Voltage:', 17X, F7.3, ' volts', /,
	1          1X, 'Detector Responsivity (S0):', 14X, 1PG10.3, /,
	1          1X, 'Detector Time Constant:', 18X, 1PG10.3, /,
	1          1X, 'Derived Bolometer Temp.:', 17X, F7.4, ' deg K', / )

	Endif

C Transcribe the Covariance matrix from fcc_covar (522,522) to covar (368,368).

	rstat = FIP_TRANSCRIBE_COV ( fcc_covar, covar, index_low, index_high )

C Pack covariance matrix COVAR into the three vector form for storage in
C the FIP record.

	rstat = FIP_PACK_COV ( covar, real_real, real_imag, imag_imag )

C Get a logical unit number

	rstat = Lib$Get_Lun (lun_out)
	If ( rstat .NE. SS$_Normal ) Then
	   FIP_WRITE_COV = %loc (fip_abort)
	   Call Lib$Signal (fip_lunerr, %val(1), %val(rstat))
	Else

C Open the output covariance matrix data file.

	   outfile = arch_out // 'FIP_COV_' // chan // scan // '.' // filext
	   Open ( Unit=lun_out, File=outfile, Status='NEW', Iostat=rstat,
	1         Useropen=CT_Connect_Write )

	   If (rstat .NE. 0) Then
	      FIP_WRITE_COV = %loc (fip_abort)
	      Call Lib$Signal (fip_openerr, %val(2), outfile, %val(rstat))
	   Else
	      If (report) Write (lun_rpt, 50) 'opened', outfile
   50	      Format (1X, 'Successfully ', A, ' ', A)
	   Endif
	Endif

C If everything ok, write the data.
C  Stored in 256 records: each record contains 
C   a.  General information - - the same for each record
C   b.  Information for one bin of the 256 categories of binned dispersion vectors
C   c.  Information for one frequency bin of the real and imaginary
C          parts of the average calibrated spectra (up to record 180)
C   d.  Up to 270 covariance matrix values.

	rr_element = 0
	ri_element = 0
	ii_element = 0
	rec = 1
	Do While ((FIP_WRITE_COV .EQ. %loc (fip_normal)) .AND. (rec .LE. 256))

C Initialize fip_cov record (5280 bytes of 0).

	   Call Lib$Movc5 ( 0, , 0, 5280, fip_rec )

C Put calibration model solution identification into the record.

	   fip_rec.Model_Label = related_data.Model_Label
	   fip_rec.Cov_Label = related_data.Cov_Label
	   fip_rec.ChanScan = chan // scan

C Put frequency range and info of model solution into the record.

	   fip_rec.Nu_Zero     = nu_zero
	   fip_rec.Delta_Nu    = delta_nu
	   fip_rec.Num_Freq    = index_high - index_low + 1
	   fip_rec.Num_Ifgs    = related_data.Num_Ifgs
	   fip_rec.Num_Coadds  = related_data.Num_Coadds
	   fip_rec.Deg_Freedom = related_data.Deg_Freedom

C Put average bolometer temperature, commanded bias, readout voltage,
C detector responsivity, time constant, and derived bolometer temperature
C used in calibration of the covariance matrix into the record.

	   fip_rec.Bol_Avg     = related_data.Bol_Avg
	   fip_rec.Cbias_Avg   = cbias_volts
	   fip_rec.Volt_Avg    = related_data.Volt_Avg
	   fip_rec.S0          = related_data.S0
	   fip_rec.Tau         = related_data.Tau
	   fip_rec.Tbol        = related_data.Tbol

C For the following 2 fields: record 1 contains info for frequency Nu_zero,
C record 2 contains info for frequency Nu_zero + delta_nu, and so forth up
C through num_freq.

	   i = rec + index_low - 1
	   If ( i .LE. index_high ) Then
	      fip_rec.Avg_Rcalspec  = related_data.Avg_Rcalspec (i)
	      fip_rec.Avg_Icalspec  = related_data.Avg_Icalspec (i)
	   Endif

C The following is bin information for the record's bin.

	   fip_rec.Bin_Total     = bin_total (rec)
	   fip_rec.Wtd_Bin_Total = wtd_bin_total (rec)
	   fip_rec.Eff_Wt        = eff_wt (rec)

C Dispersion vectors are stored, up to 180 frequencies for each bin.  The
C first data point which is at index_low gets put into the first element
C of the output array.

	   Do freq = index_low, index_high
	      fip_rec.Rdisp (freq-index_low+1) = rdisp (freq, rec)
	      fip_rec.Idisp (freq-index_low+1) = idisp (freq, rec)
	   Enddo

	   Do i=1,8
	      fip_rec.Temp (i) = temp (i, rec)
	   Enddo

C Store the covariance matrix vectors in records 1 - 253, each with up to
C 270 values:
C       Real - real (17766 values) stored in records 1 - 66;
C       Real - imag (33840 values) stored in records 67 - 192;
C       Imag - imag (16290 values) stored in records 193 - 253;

	   If ( rec .LE. 66 ) Then
	      Do i=1,270
	         rr_element = rr_element + 1
	         If (rr_element .LE. 17766) Then
	            fip_rec.covar (i) = real_real (rr_element)
	         Endif
	      Enddo
	   Elseif ( rec .LE. 192 ) Then
	      Do i=1,270
	         ri_element = ri_element + 1
	         If (ri_element .LE. 33840) Then
	            fip_rec.covar (i) = real_imag (ri_element)
	         Endif
	      Enddo
	   Elseif ( rec .LE. 253 ) Then
	      Do i=1,270
	         ii_element = ii_element + 1
	         If (ii_element .LE. 16290) Then
	            fip_rec.covar (i) = imag_imag (ii_element)
	         Endif
	      Enddo
	   Endif

C Write record to output file.

	   Call Ct_Write_Arcv ( , lun_out, fip_rec, ct_stat )
	   If (ct_stat(1) .NE. ctp_normal) Then
	      FIP_WRITE_COV = %loc (fip_abort)
	      Call Lib$Signal (fip_writerr, %Val(2), outfile,%loc (ct_stat(1)))
	   Endif                                               ! write error

	   rec = rec + 1

	Enddo                                           ! for each record (256)

C Close output file

	If (FIP_WRITE_COV .EQ. %loc (fip_normal)) Then
	   Call CT_Close_Arcv (, lun_out, ct_stat)
	   If (ct_stat(1) .NE. ctp_normal) Then
	      FIP_WRITE_COV = %loc (fip_abort)
	      Call Lib$Signal (fip_closerr, %val(2), outfile, %val(ct_stat(1)))
	   Else
	      If (report) Write (lun_rpt, 50) 'closed', outfile
	   EndIf	
	Endif
	Return
	End
