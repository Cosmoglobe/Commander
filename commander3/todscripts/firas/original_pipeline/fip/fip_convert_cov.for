	Integer*4  Function  FIP_CONVERT_COV  ( lun_rpt, report, chan, scan,
	1                                       related_data, bin_total,
	1                                       wtd_bin_total, eff_wt, rdisp,
	1                                       idisp, temp, fcc_covar,
	1                                       freq_low, freq_high, index_low,
	1                                       index_high, arch_ref, delta_nu,
	1                                       nu_zero )

C Convert the values in the covariance matrix from eplees^2 to MJy^2/sr^2,
C and the associated quantities from eplees to MJy/sr.  Convert the
C detector responsivity to V/W.
C
C Larry P. Rosen, Hughes STX, 20 July 1993.
C
C Input:
C	Integer*4     lun_rpt              ! Logical unit number - report file
C	Logical*1     report               ! Flag whether to write report
C	Character*2   chan                 ! Channel
C	Character*2   scan                 ! Scan mode
C	Record / fcc_cov / related_data    ! Data in the first FCC record
C	Integer*4     bin_total (256)      ! Number of voltage spectra
C	Real*4        wtd_bin_total (256)  ! Weighted total voltage spectra
C	Real*8        eff_wt (256)         ! Effective weight of each volt spec
C	Real*8        rdisp (257,256)      ! Total real spectra dispersion
C	Real*8        idisp (257,256)      ! Total imag spectra dispersion
C	Real*8        temp (8,256)         ! Temp and glitch rate disp. total
C	Real*8        fcc_covar (522,522)  ! Covariance Matrix
C	Real*4        freq_low             ! Low frequency cut-off (icm)
C	Real*4        freq_high            ! High frequency cut-off (icm)
C	Character*15  arch_ref             ! Archive containing reference file.
C Output:
C	Integer*2     index_low            ! Index for Low frequency cut-off
C	Integer*2     index_high           ! Index for High frequency cut-off
C	Real*4        delta_nu             ! Optical frequency interval, GHz.
C	Real*4        nu_zero              ! Optical frequency of initial data
C                                          !    point, in GHz.
C
C    PDL:
C
C Call FIP_GET_NYQUIST to get the nyquist frequency from the FEX_NYQUIST file.
C Call FIP_FREQUENCY_CUT to determine the index range of the frequency bins
C for the delivered product.
C
C 1.  For the following vectors, convert each quantity from units of Eplees to
C     MJy per sr for the designated frequencies.  Let BIN_1 designate the first
C     frequency bin, LAST_BIN designate the last, and NF =
C     (LAST_BIN - BIN_1 + 1) for convenience in the following. These bins will
C     be written in FIP_WRITE_COV to bins 1 through NF of the real*8 length
C     180 array of the same name in the FIP_COV output file.
C     a.  AVG_RCALSPEC:  bins BIN_1 through LAST_BIN of the real*8 length 257
C         array containing the overall average of the real parts of calibrated
C         spectra formed by the coadds.
C     b.  AVG_ICALSPEC:  bins BIN_1 through LAST_BIN of the real*8 length 257
C         array containing the overall average of the imaginary parts of
C         calibrated spectra formed by the coadds.
C     c.  RDISP:  for rows 1 - 256, bins BIN_1 through LAST_BIN of the real*8
C         256 x 257 array containing the real part of binned dispersion vector
C         totals.
C     d.  IDISP:  for rows 1 - 256, bins BIN_1 through LAST_BIN of the real*8
C         256 x 257 array containing the imaginary part of binned dispersion
C         vector totals.
C 2.  FIP_READ_COV reconstructed the real*8 522 x 522 covariance matrix,
C     henceforth denoted FCC_COVAR.  Convert it as follows:
C     a.  For the following indices, the units in the covariance matrix
C         FCC_COVAR(i,j) shall be converted from Eplees2 to MJy2 per sr2:
C         i.   1 <= i <= 257 and 1 <= j <= 257.
C         ii.  1 <= i <= 257 and 266 <= j <= 522.
C         iii. 266 <= i <= 522 and 1 <= j <= 257 (conversion of this piece is
C              not technically necessary since only the upper triangular
C              portion of the matrix will be stored).
C         iv.  266 <= i <= 522 and 266 <= j <= 522.
C     b.  For the following indices, the units in the covariance matrix are
C         Eplees times various quantities (temperature, glitch rate, etc.);
C         hence FCC_COVAR(i,j) should be multiplied only the factor which
C         converts from Eplees to Mjy per sr:
C         i.   1 <= i <= 257 and 258 <= j <= 265.
C         ii.  258 <= i <= 265 and 1 <= j <= 257.
C         iii. 266 <= i <= 522 and 258 <= j <= 265 (conversion of this piece is
C              not technically necessary since only the upper triangular
C              portion of the matrix will be stored).
C         iv.  258 <= i <= 265 and 266 <= j <= 522 (conversion of this piece is
C              also not technically necessary).
C     c.  For 258 <= i <= 265 and 258 <= j <= 265, no conversion of
C         FCC_COVAR(i,j) is necessary.
C
C  Convert the detector responsivity to V/W.
C------------------------------------------------------------------------------
C
C  Changes:
C
C	Converted FIP_CONFIG_LINES.TXT to FIP_CONFIG_FREQ.TXT and removed
C	FIP_INVOC_LINES.TXT
C	Gene Eplee, GSC, 10 February 1994.
C
C------------------------------------------------------------------------------

	Implicit None

C  Include Files

	Include		'($ssdef)'
	Include		'(fip_frequency)'
	Include 	'(fut_params)'       ! for conversion Eplees to MJy/sr
	Include 	'(fip_config_freq)'
	Include		'(cct_query_catalog_record)'

C Passed Parameters:

	Integer*4     lun_rpt              ! Logical unit number - report file
	Logical*1     report               ! Flag whether to write report
	Character*2   chan                 ! Channel
	Character*2   scan                 ! Scan mode
	Dictionary    'fcc_cov'            ! FCC data structure.
	Record / fcc_cov / related_data    ! Data in the first FCC record
	Integer*4     bin_total (256)      ! Number of voltage spectra
	Real*4        wtd_bin_total (256)  ! Weighted total voltage spectra
	Real*8        eff_wt (256)         ! Effective weight of each volt spec
	Real*8        rdisp (257,256)      ! Total real spectra dispersion
	Real*8        idisp (257,256)      ! Total imag spectra dispersion
	Real*8        temp (8,256)         ! Temp and glitch rate disp. total
	Real*8        fcc_covar (522,522)  ! Covariance Matrix
	Real*4        freq_low             ! Low frequency cut-off (icm)
	Real*4        freq_high            ! High frequency cut-off (icm)
	Integer*2     index_low            ! Index for Low frequency cut-off
	Integer*2     index_high           ! Index for High frequency cut-off
	Character*15  arch_ref             ! Archive containing reference file.
	Real*4        delta_nu             ! Optical frequency interval, GHz.
	Real*4        nu_zero              ! Optical frequency of initial data
                                           !    point, in GHz.

C Functions

	Integer*4     FIP_FREQUENCY_CUT
	Integer*4     FIP_READ_NYQUIST

C Local

	Integer*4     rstat
	Integer*2     findex               ! frequency
	Integer*2     bin                  ! counter for bin number
	Integer*2     i, j                 ! counters for covar matrix
	Real*8        Esqrd_2_MJ2          ! conversion factor squared
                                           !    for Eplees^2 to (MJy/sr)^2
	integer*4	cct_query_catalog
	record /query_catalog/ query_cat
	dictionary 'ccm_cme_catalog_entry'
	record /ccm_cme_catalog_entry/ cat
	Character*15	ref_filename /'FEX_NYQUIST.DAT'/

C External error messages

	External fip_normal
	External fip_abort
	External fip_numfreqerr

 
C Begin

	FIP_CONVERT_COV = %loc (fip_normal)

C The purpose of the following sections is to set the variables for the
C call to FIP_READ_NYQUIST.  References to variables for which you can't find
C the declarations or passes are in Common blocks (yecch!) in the FIP_*_LINES
C include files.

	Call CT_gmt_to_binary (ref_gmt_start, ref_start)
	Call CT_gmt_to_binary (ref_gmt_stop, ref_stop)

C     Query the catalog for the start and stop times so that they may be
C     copied to the output covariance matrix.

	query_cat.archive_id = arch_ref
	query_cat.filename = ref_filename
	rstat = cct_query_catalog (query_cat, cat)
	ref_time (1) = cat.initial_time (1)
	ref_time (2) = cat.initial_time (2)
	If (chan .EQ. 'RH') Then
	   fcc_fchan = 1
	ElseIf (chan .EQ. 'RL') Then
	   fcc_fchan = 2
	ElseIf (chan .EQ. 'LH') Then
	   fcc_fchan = 3
	Else
	   fcc_fchan = 4
	EndIf
	If (scan .EQ. 'SS') Then
	   fcc_fsmode = 1
	ElseIf (scan .EQ. 'SF') Then
	   fcc_fsmode = 2
	ElseIf (scan .EQ. 'LS') Then
	   fcc_fsmode = 3
	Else
	   fcc_fsmode = 4
	EndIf

C Call FIP_READ_NYQUIST to get the nyquist frequency from the FEX_NYQUIST file.

	rstat = FIP_READ_NYQUIST ()

C Call FIP_FREQUENCY_CUT to determine the index range of the frequency bins.
C Note that include file fip_frequency contains the fcc parameters used here.

	fcc_freq = fac_present
	fcc_lofreq = freq_low
	fcc_hifreq = freq_high

	rstat = FIP_FREQUENCY_CUT ( fnyq_icm, fnyq_hz)

	index_low = fcc_jlo
	index_high = fcc_jhi
	delta_nu = fcc_dnu
	nu_zero = fcc_nu0

C There must be less than 181 frequencies for fip_covar to work.  Write error
C if there are too many, else write info to report file.

	If (fcc_nfreq .GT. 180) Then
	   FIP_CONVERT_COV = %loc (fip_abort)
	   Call Lib$Signal (fip_numfreqerr, %val(1), %val(fcc_nfreq))
	Elseif (report) Then

C If specified frequencies in icm are integer or almost integer, write as
C whole number, else as real with 2 decimal places.

	   If ((MOD (freq_low, 1.) .LT. 0.01) .AND.
	1      (MOD (freq_high, 1.) .LT. 0.01)) Then
	      Write (lun_rpt, 10) freq_low, freq_high
   10	      Format (1X, 'Frequency Range:', 25X, F4.0, ' - ', F4.0, '  icm')
	   Else
	      Write (lun_rpt, 20) freq_low, freq_high, index_low, index_high,
	1                         fcc_nfreq
   20	      Format (1X, 'Frequency Range:', 25X, F6.2, ' - ', F6.2, '  icm')
	   Endif
	   Write (lun_rpt, 30) index_low, index_high
   30	   Format (1X, 'Frequency Indices:', 23X, I4, ' - ', I4)
	   Write (lun_rpt, 40) fcc_nfreq
   40      Format (1X, 'Number of Frequency Points:', 14X, I4)
	Endif

C 1.  Convert each quantity from units of Eplees to MJy per sr:
C	related_data.avg_rcalspec(257), related_data.avg_icalspec(257)
C	rdisp (257,256), idisp (257,256).  Eplee = ergs/sec/cm^2/sr/icm.

	If (FIP_CONVERT_COV .EQ. %loc (fip_normal)) Then
	   Do findex = index_low, index_high
	      related_data.avg_rcalspec (findex) =
	1        related_data.avg_rcalspec (findex) * fac_erg_to_mjy
	      related_data.avg_icalspec (findex) =
	1        related_data.avg_icalspec (findex) * fac_erg_to_mjy
	      Do bin = 1, 256
	         rdisp (findex, bin) = rdisp (findex, bin) * fac_erg_to_mjy
	         idisp (findex, bin) = idisp (findex, bin) * fac_erg_to_mjy
	      Enddo
	   Enddo

C 2. Convert the real*8 522 x 522 covariance matrix fcc_covar:
C     a.  For fcc_covar (i,j) convert from Eplees2 to MJy2 per sr2:
C           1 <= i <= 257 & 1 <= j <= 257;    1 <= i <= 257 & 266 <= j <= 522.
C         266 <= i <= 522 & 1 <= j <= 257;  266 <= i <= 522 & 266 <= j <= 522.
C     b.  For fcc_covar (i,j) convert from Eplees to Mjy per sr:
C           1 <= i <= 257 & 258 <= j <= 265;  258 <= i <= 265 & 1 <= j <= 257.
C         266 <= i <= 522 & 258 <= j <= 265;  258 <= i <= 265 & 266 <= j <= 522

	   Esqrd_2_MJ2 = fac_erg_to_mjy * fac_erg_to_mjy

	   Do i=1,257
	      Do j=1,257
	         fcc_covar (i,j) = fcc_covar (i,j) * Esqrd_2_MJ2
	      Enddo
	      Do j=266,522
	         fcc_covar (i,j) = fcc_covar (i,j) * Esqrd_2_MJ2
	      Enddo
	   Enddo
	   Do i=266,522
	      Do j=1,257
	         fcc_covar (i,j) = fcc_covar (i,j) * Esqrd_2_MJ2
	      Enddo
	      Do j=266,522
	         fcc_covar (i,j) = fcc_covar (i,j) * Esqrd_2_MJ2
	      Enddo
	   Enddo
C
	   Do j=258,265
	      Do i=1,257
	         fcc_covar (i,j) = fcc_covar (i,j) * fac_erg_to_mjy
	      Enddo
	      Do i=266,522
	         fcc_covar (i,j) = fcc_covar (i,j) * fac_erg_to_mjy
	      Enddo
	   Enddo
	   Do i=258,265
	      Do j=1,257
	         fcc_covar (i,j) = fcc_covar (i,j) * fac_erg_to_mjy
	      Enddo
	      Do j=266,522
	         fcc_covar (i,j) = fcc_covar (i,j) * fac_erg_to_mjy
	      Enddo
	   Enddo

C Convert bolometer responsivity from volts/ (ergs/sec) to volts / watts

	   related_data.s0 = 1.d7 * related_data.s0


	Endif

	Return
	End
