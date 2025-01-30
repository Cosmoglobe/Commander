	subroutine fep_dmakeplot (config_lun,nrec, realval, dwell_loc, 
	1	                  timetags, jside)

C------------------------------------------------------------------------------
C Changes:
C	SPR 2910, Show all dwell data in the timerange.  Fred Shuman, STX,
C		1989 Feb 17.
C
C	SPR 2738, Add the ability to plot cal resistors.  Plot them in counts.
C		  Fred Shuman, STX,  1989 Feb 20.
C       
C       Spr 3253  Dwellplot has to fetch the reference data from the 
C                 reference data archive.
C                 H. Wang, STX, Mar. 2, 1990.  
C------------------------------------------------------------------------------
	Implicit None

	integer*2 nrec, dwell_loc(1)
	integer*4 jside, ifail, ngrt, grtnums(9), icurrent, config_lun
	real*4 realval(1)
	character*14 timetags(1000)
c
c  convert resistance values into temperatures, then plot everything
c
	call fep_dtemperatures (nrec, timetags,realval, dwell_loc, jside,
	2                       config_lun,ngrt, grtnums, icurrent, ifail)
	if (ifail .ne. 0) return
	call fep_dplotter (nrec, realval, dwell_loc, timetags, jside,
	2                  ngrt, grtnums, icurrent)
	return
	end
