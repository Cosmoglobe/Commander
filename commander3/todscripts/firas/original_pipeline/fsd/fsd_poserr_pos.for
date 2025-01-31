	SUBROUTINE FSD_POSERR_POS(FARR,NPT,MEAN,SD,DMEAN,DSD)
C
C ***   COMPUTE MEAN,STANDARD DEVIATION OF GIVEN DATA FARR AND
C       DERIVATIVE OF FARR
C
C	INPUTS    FARR = INPUT DATA (NPT,NUM)
C		  NPT  = NUMBER OF DATA POINTS 
C
C	OUTPUTS   MEAN  = MEAN OF N VALUES AT EACH POINT
C		  SD    = RMSD VALUE AT EACH POINT
C		  DMEAN = MEAN VALUE OF DERIVATIVE AT EACH POINT
C		  DSD   = RMSD VALUE OF DERIVATIVE AT EACH POINT
C
	include '(fsd_poserr)'
	
	real*4 farr(npt,100),sd(npt),dsd(npt),mean(npt),dmean(npt)
	real*4 fn(100)
c
c ***  initialize ouput arrays
c
	call lib$movc5(0,,0,npt*4,mean)
	call lib$movc5(0,,0,npt*4,dmean)
	call lib$movc5(0,,0,npt*4,sd)
	call lib$movc5(0,,0,npt*4,dsd)
c
C ***  COMPUTE MEAN,SD
c
	do i=1,npt

	   mean(i)=0.0
	   call lib$movc5(0,,0,400,fn)
	   do j=1,num
	      fn(j)=farr(i,j)
	      mean(i)=mean(i)+fn(j)
	   enddo
	   mean(i)=mean(i)/num

	   s=0.0
	   do j=1,num
	      s=s+(fn(j)-mean(i))**2
	   enddo
	   sd(i)=sqrt(s/(num-1))

c ***	COMPUTE MEAN STANDARD DEVIATION OF DERIVATIVE

	  if (i .ne. npt)then
	   
	   call lib$movc5(0,,0,400,fn)
	   do j=1,num
	      fn(j)=farr(i+1,j)-farr(i,j)
	      dmean(i)=dmean(i)+fn(j)
	   enddo
	   dmean(i)=dmean(i)/num

	   s=0.0
	   do j=1,num
	      s=s+(fn(j)-dmean(i))**2
	   enddo
	   dsd(i)=sqrt(s/(num-1))

	endif

	enddo
c
	return
	end
