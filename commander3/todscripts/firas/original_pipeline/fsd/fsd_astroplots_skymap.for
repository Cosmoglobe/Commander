	Integer*4 Function FSD_Astroplots_Skymap(skydata,numdata,skyname)
C------------------------------------------------------------------------------
C
C	Purpose:  Take Skydata array with several quantities per pixel
C	          and produce a skymap file.
C
C	AUTHOR:   Fred Shuman, STX, 1989 Jan 29.  (SER 2373)
C
C	INVOCATION:   status = FSD_Astroplots_Skymap(skydata,numdata,skyname)
C
C	INPUT:
C	    R * 4	skydata(6144,3)
C	    I * 4	numdata(6144)
C	    Ch*64	skyname
C
C	OUTPUT:
C	    FSD_SDF_xx.<file_extension>		skymap file
C
C	SUBROUTINES CALLED:
C	    LIB$Get_Lun
C	    CSA_Open_Skymap
C	    CSA_Field_Offset_Values
C	    CSA_Write_Pixels
C	    CSA_Close_Skymap
C
C	COMMON VARIABLES USED:    None
C
C	INCLUDE FILES:            None
C
C------------------------------------------------------------------------------

	Implicit None

	Include		'($SSDEF)'

	Real      * 4	skydata(6144,3)
	Integer   * 4	numdata(6144)

	Integer   * 4	save_init(2,5)          ! Initial data times
	Integer   * 4	save_fin(2,5)           ! Final data times
	Logical   * 1   more_segments(4)        ! More channel data available
	Character *16   time_range              ! Data time range

C  Functions

	Logical   * 1	Time_Lt		        ! Ct time less than
	Logical   * 1	Time_Gt		        ! Ct time greater than

C  Local Declarations

	Integer   * 2	min, max		! Min and max index pointers
	Integer   * 2   ix, jx			! Indices
	Logical   * 1   found                   ! Found initial min and max
	Character *14   chartime		! Character time string

	Integer   * 2	ct_stat(20)
	Integer   * 4	status
	Integer   * 4	i
	Integer   * 4	lun
	Character * 64	skyname
	Integer   * 4	reclen
	Integer   * 4	numrecs
	Integer   * 2	multiblock_count /127/
	Integer   * 2	filesize
	Integer   * 2	indxlvls

	Integer   * 4	LIB$Get_Lun
	Integer   * 4	CSA_Open_Skymap
	Integer   * 4	CSA_Field_Offset_Values
	Integer   * 4	CSA_Write_Pixels
	Integer   * 4	CSA_Close_Skymap

	External	FSD_Normal
	External	FSD_AbErr1
	External	CSA_Open_Skymap

	Dictionary 'FSD_SKY'

	Record /FSD_SKY/ skyrec

	FSD_Astroplots_Skymap = %loc(FSD_Normal)

	reclen = 20
	indxlvls = 6

	numrecs = 0
	Do i=1,6144
	   If (numdata(i) .Gt. 0) Then
	      numrecs = numrecs + 1
	   End If
	End Do

	filesize = numrecs*reclen/512

	status = LIB$Get_Lun(lun)

	Open (UNIT = lun,
	2     FILE = skyname,
	3     STATUS = 'new',
	4     FORM = 'unformatted',
	5     INITIALSIZE = filesize,
	6     RECL = reclen/4,
	7     RECORDTYPE = 'fixed',
	8     USEROPEN = CSA_Open_Skymap)

	status = CSA_Field_Offset_Values (0,-1,-1,lun)

	Do i=1,6144
	   If (numdata(i) .Gt. 0) Then
	      skyrec.pixel_no    = i - 1
	      skyrec.peak_ht     = skydata(i,1)
	      skyrec.glitch_rate = skydata(i,2)
	      skyrec.degl_noise  = skydata(i,3)
	      skyrec.num_pts     = numdata(i)
	      status = CSA_Write_Pixels (lun,skyrec,1,multiblock_count)
	   End If
	End Do

	status = CSA_Close_Skymap (lun,indxlvls)

	End
