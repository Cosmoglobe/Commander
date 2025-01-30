	integer*4 function fut_apod_recnum (mtm_speed,adds,fakeit,length,
     *     chan,u_mode,linearized,irec)
c
c  Given the MTM speed and data compression (adds_per_group), this 
c  routine returns the record number for the appropriate apodization 
c  function from the direct-access apodization file.
c
c  Written by Rich Isaacman
c	      Applied Research Corp.
c	      286-5758
c	      6 Jan 1989
c
c  Revision:		Alice Trenholme
c			General Sciences Corp.
c			9 February 1990
c			286-5569
c
c	SER 6795. Increased variables used in finding peak by including
c	channel, microprocessor mode, scan length, and linearization factor
c	in input variable list.
c
	implicit none
	integer*4 	mtm_speed		! mtm speed
	integer*4	length			! length of scan
	integer*4 	adds			! adds per group
	integer*4	chan		 	! channel number
	integer*4	ch			! channel (0 = Low, 1=High)
	integer*4	u_mode			! science mode 
	integer*4	sc_m			! sci. mode (0=1 or 3;1=2 or 4)
	integer*4	irec			! record number
	integer*4	fakeit			! fakeit bit off/on
	integer*4	linearized		! linearization switch (0=not
c						   linearized, 1=linearized)

	external	fut_normal

	if ((chan .eq. 1) .or. (chan .eq. 3)) then 
		ch = 1 
	else 
		ch = 0
	endif

	if ((u_mode .eq. 1) .or. (u_mode .eq. 3)) then 
		sc_m=0
	else 
		sc_m=1
	endif

	if (fakeit .eq. 1) then
	   irec = 385
	else
	   irec = 32*(adds-1) + 16*mtm_speed + 
     .		8*length + 4*ch + 2*sc_m + linearized + 1
	endif

	fut_apod_recnum = %loc(fut_normal)
	return
	end
