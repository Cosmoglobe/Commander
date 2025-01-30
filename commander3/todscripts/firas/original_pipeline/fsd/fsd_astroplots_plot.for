	Integer*4 Function FSD_Astroplots_Plot(chan,plt_device,
	1				       plt_com,plt_com_file,
	2				       xlabel,
	3                                      att_soln)

C-----------------------------------------------------------------------------
C
C       Configured in CSDR environment by Reid Wilson,  STX Inc., 31-JAN-1988
C
C-----------------------------------------------------------------------------
C   Changes:
C	SPR 3774.  Astroplots uses only simulated attitude & orbit.  If real
C	    ones exist, use them instead.  Fred Shuman, STX.  1989 May 22.
C
C	SPR 3966.  Suppress zooming (and prompting for it) in batch.
C	    Fred Shuman, STX.  1989 Jun 15.
C
C	SPR 5697, Remove dependency on MUT_YD_TO_YMD; reorganize time
C	    convertion.  R. Kummerer / STX, February 1, 1990.
C
C       SER 4569, Convert FSD_ASTROPLOTS from TEMPLATE to PLT graphics.
C	    R. Kummerer, STX / 1990 April 26
C-----------------------------------------------------------------------------

	Implicit None

	Include '(FSD_Astroplots)'

	Integer   *  4      chan
	Character * 32      plt_device
	Integer   *  4      plt_com
	Character * 64      plt_com_file
	Character * 32      xlabel
	Character * 32      ylabel
	Character * 20      att_soln

	Integer   *  4      status
	Character *100      title(3)
	Character * 60      label
	Character * 32      plot_label
	Character *  4      ans
	Integer   *  2      istart
	Integer   *  2      istop
	Character * 14      gmt_start
	Character * 14      gmt_stop

	External	FSD_Normal

	FSD_Astroplots_Plot = %loc(FSD_Normal)

	ans = 'Y'
	istart = 1
	istop = num

	plot_label = 'ASTROPLOTS'
	Call CT_Binary_To_GMT(gmt(1,istart),gmt_start)
	Call CT_Binary_To_GMT(gmt(1,istop),gmt_stop)

	Write (label,50) gmt_start, gmt_stop, istop-istart+1
50	Format(1x, ' JStart: ', a14, ' JStop: ', a14, ' IFGs: ', i4)

	Call FSD_Astroplots_Label(label,chan,att_soln,plot_label,title)
	Call FSD_Astroplots_Display(title,xlabel,plt_device,
	2			    plt_com,plt_com_file)

	Return
	End
