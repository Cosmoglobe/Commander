	Integer*4 Function FSD_Astroplots_Param(tstart, tend, deltat, t_unit)

C-----------------------------------------------------------------------------
C
C       Configured in CSDR environment by Reid Wilson,  STX Inc., 31-JAN-1988
C
C
C	INPUT/OUTPUT ARGUMENTS:
C
C	IN:
C	    R*8      tstart
C	    R*8      tend
C	    R*8      deltat
C	    R*8      t_unit
C	  Common /astro/                [ from  Include '(FSD_Astroplots)' ]
C	    R*4      timeaxis(3000)
C	    I*2      amoon(200)
C	  Common /astro_data/           [ from  Include '(FSD_Astroplots)' ]
C	    I*4      s_num
C	    I*4      vs_num
C	    I*4      vn_num
C	    R*4      saa_lon(100)
C	    R*4      saa_lmin(100)
C	    R*4      saa_lmax(100)
C	    R*4      vabs_lon(100)
C	    R*4      vabs_lmin(100)
C	    R*4      vabs_lmax(100)
C	    R*4      vabn_lon(100)
C	    R*4      vabn_lmin(100)
C	    R*4      vabn_lmax(100)
C	  Common /astro_num/            [ from  Include '(FSD_Astroplots)' ]
C	    I*4      num
C	    I*4      nmoon
C
C	OUT:
C	  Common /astro/                [ from  Include '(FSD_Astroplots)' ]
C	    R*4      saa(200)
C	    R*4      vab(200)
C	    R*4      galactic_cross(200)
C	    I*2      amoon(200)
C	    R*4      tmoon(200)
C	  Common /astro_num/            [ from  Include '(FSD_Astroplots)' ]
C	    I*4      nsaa
C	    I*4      nvab
C	    I*4      ngalactic
C	    I*4      nmoon
C
C-----------------------------------------------------------------------------
C   Changes:
C
C	SPR 3774.  Astroplots uses only simulated attitude & orbit.  If real
C	    ones exist, use them instead.  Fred Shuman, STX.  1989 May 22.
C
C	SPR 4142.  Accept galactic latitude information parsed in
C	    FSD_ASTROPLOTS_DATA and do not use science record (which was not
C	    passed anyway as source of this information.  R. Kummerer, STX.
C	    1989 Oct 24.
C
C	SPR 4820.  Correct checking of radiation fields.  R. Kummerer, STX.
C	    1989 Oct 25.
C
C       SER 4569, Convert FSD_ASTROPLOTS from TEMPLATE to PLT graphics.
C	    R. Kummerer, STX / 1990 April 26
C
C-----------------------------------------------------------------------------

	Implicit None

C   Note:  FSD_Astroplots.txt 'Includes' FUT_Params.txt

	Include '(FSD_Astroplots)'

	Real            *8      tstart
	Real            *8      tend
	Real            *8      deltat
	Real            *8      t_unit

	Integer         *4      i
	Integer         *4      j
	Integer         *4      k
	Integer         *4      status 
	Real            *8      t
	Real            *4      time
	Real            *4      dt
	Real            *4      rdot
	Real            *4      moon_cos(4)
	Integer         *2      moona(4) /5, 15, 30, 90/
	Logical         *4      m_flag
	Logical         *4      a_flag
	Logical         *4      saa_flag
	Logical         *4      vab_flag
	Real            *4      prev_gal
	Real            *4      curr_gal
	Real            *4      m_ang
	Real            *4      t_lat
	Real            *4      t_lon

	Logical		*4	FSD_Astroplots_Check_Radfld

	External        FSD_Normal

C   Determine galactic plane crossings, moon angle, and passages through
C   the South Atlantic Anomaly and Van Allen Belts.

	Do i=1,4
	   moon_cos(i) = Cosd(Real(moona(i)))
	End Do

	saa_flag = .False.
	vab_flag = .False.

	t = tstart
	dt = (tend - tstart) / t_unit
	i = 1

C   Break the time range into one second or 60 second intervals depending
C   on the duration of the time range.

	Do While (t .Le. tend)

	   time = Sngl((t-tstart)/t_unit)

	   Call FSD_Astroplots_IntAng(time, curr_gal, m_ang, t_lat, t_lon)

C   If previous and current values of galactic latitude differ in sign,
C   signal a galactic plane crossing.

	   If (i .Eq. 1) Then
	      prev_gal = curr_gal
	   End If

	   If (prev_gal*curr_gal .Le. 0.) Then
	      prev_gal = curr_gal
	      If (t_unit .Eq. 60.) Then
	         ngalactic = ngalactic + 2
	         galactic_cross(ngalactic-1) = time - 1.4 ! Assuming 2.8m spent
	         galactic_cross(ngalactic) = time + 1.4   ! crossing gal pln
	      Else
	         ngalactic = ngalactic + 2
	         galactic_cross(ngalactic-1) = time - 1.4 / 60.
	         galactic_cross(ngalactic) = time + 1.4 / 60.
	      End If
	      If (galactic_cross(ngalactic-1) .Lt. 0.)
	1		galactic_cross(ngalactic-1) = 0.
	      If (galactic_cross(ngalactic) .Gt. dt)
	1		galactic_cross(ngalactic) = dt
	   End If

C   Check the proximity of the Moon.

	   j = 1
	   m_flag = .True.

	   rdot = Cosd( m_ang )

	   Do While (m_flag)
	      If (rdot .Gt. moon_cos(j)) Then
	         m_flag = .False.

C   If NMOON = 0 then we mustn't even attempt to evaluate AMOON(NMOON)

	         If (nmoon .Eq. 0) Then
	            nmoon = nmoon + 1
	            tmoon(nmoon) = time
	            amoon(nmoon) = moona(j)
	         Else If ( amoon(nmoon) .Ne. moona(j) ) Then
	            nmoon = nmoon + 1
	            tmoon(nmoon) = time
	            amoon(nmoon) = moona(j)
	         End If
	      Else
	         j = j + 1
	         If (j .Gt. 4) Then
	            m_flag = .False.
	         End If
	      End If
	   End Do

C   Check if COBE is in the SAA.

	   a_flag = FSD_Astroplots_Check_Radfld(radfld.saa1_lat,
	1					radfld.saa2_lat,
	2					radfld.saa1_lon,
	3					radfld.saa2_lon,
	4					radfld.saa1_num,
	5					radfld.saa2_num,
	6					t_lat, t_lon)

	   If (saa_flag .Ne. a_flag) Then
	      nsaa = nsaa + 1
	      saa(nsaa) = time
	      saa_flag = a_flag 
	   End If

C   Check if COBE is in a VAB.

	   If (t_lat .Le. 0) Then
	      a_flag = FSD_Astroplots_Check_Radfld(radfld.vabs1_lat,
	1					   radfld.vabs2_lat,
	2					   radfld.vabs1_lon,
	3					   radfld.vabs2_lon,
	4					   radfld.vabs1_num,
	5					   radfld.vabs2_num,
	6					   t_lat, t_lon)
	   Else
	      a_flag = FSD_Astroplots_Check_Radfld(radfld.vabn1_lat,
	1					   radfld.vabn2_lat,
	2					   radfld.vabn1_lon,
	3					   radfld.vabn2_lon,
	4					   radfld.vabn1_num,
	5					   radfld.vabn2_num,
	6					   t_lat, t_lon)
	   End If

	   If (vab_flag .Ne. a_flag) Then
	      nvab = nvab + 1
	      vab(nvab) = time
	      vab_flag = a_flag
	   End If

	   i = i + 1
	   t = tstart + deltat*i

	End Do

	If (saa_flag) Then
	   nsaa = nsaa + 1
	   saa(nsaa) = timeaxis(num)
	End If

	If (vab_flag) Then
	   nvab = nvab + 1
	   vab(nvab) = timeaxis(num)
	End If


	FSD_Astroplots_Param = %loc(FSD_Normal)

	Return
	End
