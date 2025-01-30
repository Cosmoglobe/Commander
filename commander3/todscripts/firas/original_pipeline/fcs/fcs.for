	Program FCS
c-------------------------------------------------------------------------------
c	FCS - FIRAS Combine Spectra is the final facility in the FIRAS
c	pipeline.  Its purpose is to combine spectra, averaging data and
c	calculating total variances.
c
c	Author: Larry P. Rosen, March 1992, Hughes STX, SER 10960
c
c	Modifications:  February 1992 - add qualifiers EXC_GAL and FREQ.
c	July 1994 - weighting for spectra and other quantities change from
c	just number of ifgs to   nifgs / (i + g * s)  where g is the glitch
c	rate, i and s are parameters from glitchcorr reference data dependent
c	on channel and scan mode.
c
c	Methods: (See requirements document and design document.  Excerpt:)
c
c	CALCULATION OF THE MEANS AND VARIANCES
c
c	1.  Some Notation
c		X(n,i) = n-th complex spectrum for the i-th pixel.
c		M(n,i) = number of IFGs which were averaged to obtain X(n,i).
c		S(n,i) = variance vector for the n-th spectrum, i-th pixel
c		N = number of spectra from the i-th pixel
c		K(n,i) = 1 if the secondary template was not subtracted from
c		coadd n and 2 if the secondary  template was subtracted.
c		GR = glitch rate during spectrum collection
c		GI = glitch model intercept
c		GS = glitch model slope
c
c	2.  The Means and Variances
c		Wt(n,i) = M(n,i) / (GI + GR * GS) is the weighting for FCS.
c		The weighted mean and variance of the mean of the N spectra
c		selected for the i-th pixel are given by:
c
c		<X(i)> = [Sum(n) Wt(n,i)*X(n,i)]/[Sum(n) Wt(n,i)]
c		<S(i)> = [Sum(n) M(n,i)*(M(n,i)-K(n,i))*S(n,i)] /
c		         [Sum(n) M(n,i)]*[Sum(n) M(n,i) - Sum(n) K(n,i)]
c		where Sum(n) means summation over n.
c
c	3.  Means and Variance of other Coadd Data
c		Some other non-engineering data associated with the output
c		coadd records should also be averaged using the sweeps times
c		number of ifgs weighting scheme.  For data values {Y(n,i)} the
c		mean for the i-th pixel is given by:
c
c		<Y(i)> = [Sum(n) Wt(n,i)*Y(n,i)]/[Sum(n) Wt(n,i)]
c
c		The associated variance must be calculated from the available
c		values (which are actually means) and variances:
c
c		V(i) = [ {Sum(n) V(n,i)*[Wt(n,i)-1] + Wt(n,i)*Y(n,i)^2} -
c		         {Sum(n) Wt(n,i)}*<Y(i)>^2 ] / [Sum(n) Wt(n,i)-1]
c		NOTE:  the quantity stored in the RDL is SQRT(V(i)).
c
c		The three sigmas PC_SIGMA, QRAD_SIGMA, and IR_POWER_SIGMA are
c		calculated directly:  the variance of the measurements {Y(i)}:
c
c		V(i)={Sum(n) Wt(n,i)*[Y(n,i)-<Y(i)>]^2}/[Sum(n) Wt(n,i)-1]
c
c		The variances of the "B" and "C" coefficients relating to the
c		subtraction of the secondary template and digital transient
c		respectively (COAD_SPEC_DATA.SEC_TEMPLATE.B_VARIANCE and
c		COAD_SPEC_DATA.TRANSIENT.C_VARIANCE) are calculated as:
c
c		V(i) = {Sum(n) V(n,i)*[M(n,i)-1]}/[Sum(n) Wt(n,i)-1]
c		     + {Sum(n) Wt(n,i)*[Y(n,i)-<Y(i)>]^2}/[Sum(n) Wt(n,i)-1]
c
c		The secondary template variance stored in the coadd record
c		(COAD_SPEC_DATA.SEC_TEMPLATE.VARIANCE) is a special case
c		because it is not a per-ifg variance; instead is is the
c		variance per bin.  The variance of the mean is calculated as:
c
c		V(i) = {Sum(n) Wt(n,i)*V(n,i)*511}/[(512*Sum(n) Wt(n,i))-1]
c
c	4.  Means of the Orbital and Engineering Data
c		Each combined spectrum record should contain averages of the
c		engineering and orbital data.  These will be weighted means and
c		total expected variance. Let Y(n,i) be a particular engineering
c		quantity, and V(n,i) its associated variance (standard
c		deviation SQUARED), if applicable.  The mean values are
c		calculated in the same manner as the mean vectors above.
c		The mean for the i-th pixel is given by
c
c		<Y(i)> = [Sum(n) Y(n,i)*Wt(n,i)]/[Sum(n) Wt(n,i)]
c
c		The associated variance is given by:
c
c		V(i) = [ {Sum(n) V(n,i)*[M(n,i)-1] +M(n,i)*Y(n,i)^2} 
c		         -[Sum(n) M(n,i)]*<Y(i)>^2 ] / [Sum(n) M(n,i)-1]
c
c		NOTE:  the quantity stored in the RDL is SQRT(V(i)).
c
c	The attitude unit vector will have to be renormalized after being
c	combined.  Other engineering and attitude parameters which refer to
c	maximum and minimum values will have to be updated, in order to reflect
c	the contents of the combined spectrum.
c
c	5.  The Bin Vectors
c		A total bin vector will also be saved for each mean spectrum.
c		It will simply be the sum of the bin vectors.
c
c		BIN(i,k) = Sum(n) BIN(n,i,k)
c		where BIN(i,k) refers to the k-th bin of the i-th pixel.
c
c  Modifications:
c               L. Rosen, HSTX, 11/93, SER 11417, Remove multiple scan modes.
c               L. Rosen, HSTX, 4/94, SER 11689, Add FL scan mode.
c               L. Rosen, HSTX, 7/94, ser 11840, Change weight from Nifgs to
c                  Nifgs/(I+G*S)  where G is glitchrate, I and S are parameters
c                  from a reference data set.  This new coadd weight will be
c                  saved in the field COAD_SPEC_HEAD.COMB_NUM_IFGS for use by
c                  FMS.
c               L. Rosen, HSTX, 8/94, SER 11861, Fix logic error for check
c                  of multiple model solutions among input data when using
c                  NOREPORT qualifier.
c-------------------------------------------------------------------------------
	Implicit None

c Include files

	Include		'($ssdef)'
	Include		'(fut_params)'
	Include		'(fut_error)'
	Include		'(fcs_msg)'
	Include		'(fut_fcs_include)'
	Include		'csdr$library:ctparams.inc'
	Include		'(upm_stat_msg)'
	Include		'(csa_pixel_input_rec)'
	Include		'(csa_pixel_output_rec)'

c External

	External	fut_error
	External	csa_normal

c Functions

	Integer*4	CUT_Register_Version, CUT_Display_Banner
	Integer*4	FCS_Parse, FCS_Report_Check, FCS_Open, FUT_FCS_Non_Avg
	Integer*4	FUT_FCS_Minmax, FUT_FCS_Sum, FCS_Spec, FCS_Avg
	Integer*4	FUT_FCS_Reference, FCS_Close
	Logical*1	Time_LT, Time_LE
	Integer*4	CSA_Read_Pixels, CSA_Write_Pixels

c Local

	Integer*4	rstat				! Return status
	Integer*4	status /0/, err /2/, ok /0/	! Processing status
	Character*6	version
	Parameter	(version = '12.3  ')
	Character*80	sfile			! Script file name
	Character*3	input			! Input data type; FAD or FCF
	Character*2	chan, scan		! Channel, scan mode
	Integer*2	nscan			! Number scan modes
	Logical*1	report
	Character*36	reportf			! Filename for report
	Character*79	cmdline(3)		! Command line with defaults
	Character*56	in_skymap (fac_max_num)		! Input skymaps
	Character*14	tstart (fac_max_num)		! Start times
	Character*14	tstop (fac_max_num)		! Stop times
	Character*20	fext				! Output file extension
	Integer*2	snum				! # of input skymaps
	Logical*1	repdef		! True = use default report name
	Integer*2	clin, clen	! command line number and length
	Integer*4	lun_rpt			! Report logical unit #
	Character*14	arch_in	/'CSDR$FIRAS_IN:'/	! input data archive
	Character*15	arch_out /'CSDR$FIRAS_OUT:'/	! output data archive
	Character*15	arch_ref /'CSDR$FIRAS_REF:'/	! reference archive
	Character*60	fileout
	Logical*1	goodmap (fac_max_num)	! map has right chan & scan
	Real*8		g_intercept, g_slope	! glitch rate model
	Integer*2	ct_stat(20)			! Cobetrieve status
	Integer*4	lun_out, lun_in (fac_max_num)
	Integer*2	npix			! Number of FIRAS pixels (6144)
	Integer*2	pix			! Pixel number
	Integer*2	numpix / 0 /		! Number of pixels with data.
	Record /pixel_input_list/  inlist
	Record /pixel_output_list/ outlist
	Integer*4	maxrex			! Max # of records in a pixel
	Parameter	(maxrex=fac_max_skymap_recs)
	Dictionary	'FCF_SKY'
	Record / fcf_sky / in_rec (maxrex)
	Dictionary	'FCS_SKY'
	Record / fcs_sky / fcs_rec
	Byte		blank_rec (fac_coad_spec_size) / fac_coad_spec_size*0 /
	Record / double_fcf / sum_rec, rec1	! Record to store sums.
	Byte		blank_sum (sum_rec_size) / sum_rec_size * 0 /
	Integer*4	anum			! Counter for angle array
	Logical*1	first			! Flag first record in pixel
	Logical*1	first_model		! Flag first model written.
	Integer*2	fil				! Input file number
	Integer*4	adtstart(2), adtstop(2)		! Check times
	Integer*4	dummy				! Junk data
	Integer*2	rec				! Record number
	Integer*2	numrecs			! number of records to average
	Integer*4	nifgs				! # of IFG's in coadd
	Integer*4	total			! Total # of IFGs combined
	Character*40	multmod	/'Multiple models used among input skymaps'/
	Character*40	model(32)
	Integer*2	num_mod /1/
	Real*8		weight		! nifgs / (i + g * s) see modif. notes
	Real*8		sum_weight		! sum of weights
	Real*4		moon_phase, orb_phase, scan_angle	! some angles
	Integer*4	maxsum			! max # for sum of nifgs
	Parameter	(maxsum=maxrex*100)
	Real*4		mphase(maxsum), ophase(maxsum), sangle(maxsum)
	Integer*2	two
	Parameter	(two=2)
	Integer*4	time_array (two, maxsum)	! time array to average
	Integer*4	time_weights (maxsum)
	Integer*4	n, m				! counters
	Real*8		dof	! degrees of freedom = weight - k  where
				! k = 2 if 2nd template subtracted, 1 if NOT.
	Real*8		sum_dof		! sum of dof
c------------------------------------------------------------------------------
c Begin

	rstat = CUT_Register_Version (version)
	rstat = CUT_Display_Banner (6, 80, 'FIRAS  FCS_Combine_Spectra')

c Call CT_INIT to initialize Cobetrieve.

	Call CT_Init ( ct_stat )
	If (ct_stat(1) .NE. ctp_normal) Then
	   Call Lib$Signal (fcs_ctinit, %val(1), %val(ct_stat(1)))
           status = err
	Endif

c Call FCS_PARSE to get the command line qualifiers, open, read, and close the
c script file.

	rstat = FCS_Parse ( sfile, input, chan, scan, report, reportf,
	1                   cmdline, in_skymap, tstart, tstop, fext,
	2                   snum, repdef, clin, clen )

	If ( rstat .NE. %Loc (fcs_normal) ) status = err

c Call FCS_REPORT_CHECK to open and initialize processing report, and to check
c the script file data.

	If (status .EQ. ok) Then
	   rstat = FCS_Report_Check ( report, reportf, repdef, version,
	1                             cmdline, clin, clen, input, chan, scan,
	2                             sfile, in_skymap, tstart, tstop,
	3                             fext, snum, arch_in, arch_out, lun_rpt,
	4                             goodmap )

	   If ( rstat .NE. %Loc (fcs_normal) ) status = err
	Endif
	Call Lib$Establish ( fut_error )

c Call FCS_OPEN to open the input and output skymap files.

	If (status .EQ. ok) Then
	   rstat = FCS_Open ( in_skymap, arch_in, fileout, arch_out,
	1                     chan, scan, fext, snum, tstart, tstop,
	2                     report, lun_rpt, goodmap, lun_in, lun_out )

	   If ( rstat .NE. %Loc (fcs_normal) ) status = err
	Endif

c Call FUT_FCS_REFERENCE to open, read, and close the Glitch rate model file.

	If (status .EQ. ok) Then
	   rstat = FUT_FCS_Reference ( arch_ref, g_intercept, g_slope, chan,
	1                              scan, report, lun_rpt )

	   If ( rstat .NE. 0 ) status = err
	Endif

c Pixel Loop: do for each FIRAS pixel.

	first_model = .TRUE.
	npix = ((4 ** fac_skymap_level) * 6) - 1
	inlist.level_no = fac_skymap_level
	pix = 0
	Do While ((pix .LE. npix) .AND. (status .EQ. ok))
	   						! zeroes out fcs_rec
	   Call Lib$MovC3 (fac_coad_spec_size, blank_rec, fcs_rec)
	   fcs_rec.spec_data.combine_script = sfile
	   inlist.pixel_no = pix
	   total = 0
	   sum_weight = 0.
	   anum = 0
	   numrecs = 0
	   sum_dof = 0.
	   first = .true.
						! Blank out 1st data record
	   Call Lib$MovC3 (sum_rec_size, blank_sum, rec1)
						! Blank out sum record
	   Call Lib$MovC3 (sum_rec_size, blank_sum, sum_rec)

c File Loop: do for each input skymap file.

	   Do fil = 1, snum
	      Call CT_GMT_To_Binary (tstart (fil), adtstart)
	      Call CT_GMT_To_Binary (tstop (fil), adtstop)

c Read the records in file (fil) for this pixel number (pix).

	      rstat = CSA_Read_Pixels ( lun_in (fil), inlist, 1, in_rec,
	1                               maxrex, outlist, dummy, 0 )

	      If ( rstat .NE. %Loc (csa_normal) ) Then
	         status = err
	         Call Lib$Signal (fcs_csaread, %val(1), %val(rstat))
	      Endif

c Record Loop: do for each input record.

	      rec = 1
	      Do While ((rec .LE. outlist.no_records) .AND. (status .EQ. ok))

c Check that the record has a coadd time within the time range specified
c for this file in the script file.  If not, get the next record.

	         If ( Time_LE (in_rec (rec).ct_head.time, adtstop) .AND.
	1             Time_LE (adtstart, in_rec (rec).ct_head.time) ) Then

c If first record, call FUT_FCS_NON_AVG to store non-averaged data in output
c record. Else, check model solution for consistency.  Flag and signal if not.

	            nifgs = in_rec (rec).coad_spec_head.num_ifgs
	            numrecs = numrecs + 1
	            If (nifgs .GT. 0) Then
	               If (first) Then
	                  rstat = FUT_FCS_Non_Avg ( fcs_rec, in_rec (rec),
	1                                           scan, scan, chan, chan,
	1                                           input )

	                  If (first_model) Then
	                     first_model = .FALSE.
	                     model(1) = in_rec (rec).spec_data.model_label
	                     If (report) Then
	                        Write (lun_rpt, 10)
	1                          in_rec (rec).spec_data.model_label
   10	                        Format (2X,'Model solution: ',A40)
	                     Endif
	                  Endif
	               Else
	                  m = 0
	                  Do n = 1, num_mod
	                    If (in_rec (rec).spec_data.model_label .NE.
	1                       model(n))  m = m + 1
	                  Enddo
	                  If (m .EQ. num_mod) Then
	                     fcs_rec.spec_data.model_label = multmod
	                     Call Lib$Signal ( fcs_multmodel )
	                     num_mod = num_mod + 1
	                     model (num_mod) =
	1                       in_rec (rec).spec_data.model_label
	                     If (report) Then
	                        Write (lun_rpt, 20) model (num_mod)
   20	                        Format (2X,'Other model solution used: ',A40)
	                     Endif
	                  Elseif (num_mod .GT. 1) Then
	                     fcs_rec.spec_data.model_label = multmod
	                  Endif
	               Endif

c Extract the number of IFGs coadded in the input record.
c Increase a counter of the total number of IFGs in this pixel by the number
c of IFGs coadded in the input spectrum.
c Increase the sum of the spectrum weights for this pixel.

	               total = total + nifgs
	               weight = dble (nifgs) / ( g_intercept + g_slope *
	1                 in_rec(rec).coad_spec_data.glitch_rate )
	               sum_weight = sum_weight + weight

c Store min and max of time and engineering quantities

	               rstat = FUT_FCS_Minmax ( first, fcs_rec, in_rec (rec) )

c Form sums of engineering quantities and attitude.

	               rstat = FUT_FCS_Sum ( weight, fcs_rec,
	1                                    in_rec (rec), moon_phase,
	2                                    orb_phase, scan_angle, sum_rec,
	3                                    rec1, first )

c Form sums of spectra, variance, and bin vector.

	               rstat = FCS_Spec ( weight, fcs_rec, in_rec (rec),
	1                                 sum_rec, rec1, first, dof )

	               sum_dof = sum_dof + dof

c Store the attitude angles for averaging.

	               Do n = 1, nifgs
	                  mphase(anum+n) = moon_phase
	                  ophase(anum+n) = orb_phase
	                  sangle(anum+n) = scan_angle
	               Enddo
	               time_array (1,numrecs) = in_rec (rec).ct_head.time(1)
	               time_array (2,numrecs) = in_rec (rec).ct_head.time(2)
	               time_weights (numrecs) = nifgs
	               anum = anum + nifgs
	               If (first) first = .FALSE.
	            Endif
	         Endif
	         rec = rec + 1
	      Enddo			! for each record
	   Enddo			! for each file

c Store the sum of spectrum weights in the output FCS record in the
c coad_spec_head.comb_num_ifgs for use by FMS.

	   fcs_rec.coad_spec_head.comb_num_ifgs = sum_weight

c Compute averages and variances.

	   If ((total .NE. 0) .AND. (status .EQ. ok)) Then
	      numpix = numpix + 1
	      rstat = FCS_Avg ( total, sum_weight, fcs_rec, maxsum,
	1                       mphase, ophase, sangle, sum_rec,
	2                       rec1, time_array, time_weights, numrecs,
	3                       sum_dof, scan )

c Write combined record for current pixel to output file.

	      rstat = CSA_Write_Pixels (lun_out, fcs_rec, 1, )

	      If ( rstat .NE. %Loc (csa_normal) ) Then
	         status = err
	         Call Lib$Signal (fcs_csawrite, %val(1), %val(rstat))
	      Endif
	   Endif
	   pix = pix + 1
	Enddo			! for each pixel

c Call FCS_CLOSE to close input data sets, close output dataset,
c write summary to and close processing report.  Signal completion status.

	rstat = FCS_Close (lun_in, snum, report, lun_rpt, lun_out,
	1                  in_skymap, fileout, numpix, status, ok)

	If ( rstat .NE. %Loc (fcs_normal) ) status = err

c Exit

	If (status .EQ. ok) Then
	   Call Exit (ss$_normal)
	Else
	   Call Exit (ss$_abort)
	Endif
	Stop
	End
