MODULE iotools
  !
  ! v1.0 : Olivier Dore et Eric Hivon, June-July 2002
  ! v1.1     2002-07-08 : bugs correction by E.H. 
  ! v1.2     2002-08-22 : pixel number starts at 0 for these routine,
  !                         but starts at 1 for cftisio routines ...
  
  USE healpix_types
  USE fitstools , ONLY : printerror

  IMPLICIT NONE
  
!   INTEGER, PARAMETER :: i8b = SELECTED_INT_KIND(18) ! already declared in Healpix 1.2

!   INTEGER(I8b) , PARAMETER :: nchunk_max = HUGE(1_i4b)/2
  INTEGER(I4B) , PARAMETER :: nchunk_max = 12000 ! faster

  INTERFACE input_tod_fits
     MODULE PROCEDURE input_tod_s,input_tod_d
  END INTERFACE

  PRIVATE
  
  PUBLIC input_tod_fits
  PUBLIC write_bintabh

CONTAINS
  ! *************************************************************************
  SUBROUTINE input_tod_s(filename,tod,ntods,npixtot,firstpix,fmissval)
  ! *************************************************************************

    ! ===================================================================================
    !     reads fits file
    !     filename = fits file (input)
    !     tod      = data rad from the file (ouput) = real*4 array of size (npixtot,ntods)
    !     npixtot  = length of the tod (input)
    !     ntods    = number of tods
    !     fmissval = OPTIONAL argument (input) with the value to be given to missing
    !             pixels, its default value is 0
    !     firstpix = first pixel to be read (starting at 0)
    !                                               O.D. & E.H. 06/02 @ IAP
    !      2002-07-08 : bugs correction by E.H. 
    !       (consistent use of fmiss_effct)
    ! =======================================================================
    
    IMPLICIT NONE

    INTEGER(I8B),     INTENT(IN)           :: npixtot,firstpix
    INTEGER(I4B),     INTENT(IN)           :: ntods
    REAL(SP),         INTENT(OUT)          :: tod(0:,1:)
    CHARACTER(LEN=*), INTENT(IN)           :: filename
    REAL(SP),         INTENT(IN), OPTIONAL :: fmissval
    
    INTEGER(I4B) :: unit,i,itod
    REAL(SP)     :: fmissing, fmiss_effct
    INTEGER(I4B) :: imissing

    LOGICAL(LGT) :: ok,anynull
    
    !-----------------------------------------------------------------------
    fmiss_effct = 0.
    IF (PRESENT(fmissval)) fmiss_effct = fmissval

    CALL read_bintod_s(filename,tod,npixtot,ntods,firstpix,fmissing,anynull)

    DO itod = 1, ntods
       anynull = .TRUE.
       IF (anynull) THEN
          imissing = 0
          DO i=0,npixtot-1
             IF ( ABS(tod(i,itod)/fmissing -1.) .LT. 1.e-5 ) THEN
                tod(i,itod) = fmiss_effct
                imissing = imissing + 1
             ENDIF
          ENDDO
          IF (imissing .GT. 0) THEN
             WRITE(*,'(a,1pe11.4)') 'blank value : ' ,fmissing
             WRITE(*,'(i7,a,f7.3,a,1pe11.4)') &
                  &           imissing,' missing pixels (', &
                  &           (100.*imissing)/npixtot,' %),'// &
                  &           ' have been set to : ',fmiss_effct
          ENDIF
       ENDIF
    ENDDO
    RETURN

  END SUBROUTINE input_tod_s
  !=======================================================================

  !**************************************************************************************
  SUBROUTINE input_tod_d(filename,tod,ntods,npixtot,firstpix,fmissval)
  !**************************************************************************************

    !=======================================================================
    !     reads fits file
    !     filename = fits file (input)
    !     tod      = data rad from the file (ouput) = real*4 array of size (npixtot,ntods)
    !     npixtot  = length of the tod (input)
    !     ntods     = number of tods
    !     fmissval  = OPTIONAL argument (input) with the value to be given to missing
    !             pixels, its default value is 0
    !     firstpix  = first pixel to be read (starting at 0)
    !                                                   O.D. & E.H. @ IAP
    !      2002-07-08 : bugs correction by E.H. 
    !       (consistent use of fmiss_effct)
    !=======================================================================
    
    IMPLICIT NONE

    INTEGER(I8B),     INTENT(IN)           :: npixtot,firstpix
    INTEGER(I4B),     INTENT(IN)           :: ntods
    REAL(DP),         INTENT(OUT)          :: tod(0:npixtot-1,1:ntods)
    CHARACTER(LEN=*), INTENT(IN)           :: filename
    REAL(DP),         INTENT(IN), OPTIONAL :: fmissval
    
    INTEGER(I4B) :: unit,i,itod
    REAL(DP)     :: fmissing, fmiss_effct
    INTEGER(I4B) :: imissing

    LOGICAL(LGT) :: ok,anynull
    
    !-----------------------------------------------------------------------
    
    fmiss_effct = 0.
    IF (PRESENT(fmissval)) fmiss_effct = fmissval
    
    CALL read_bintod_d(filename,tod,npixtot,ntods,firstpix,fmissing,anynull)
    
    DO itod = 1, ntods
       anynull = .TRUE.
       IF (anynull) THEN
          imissing = 0
          DO i=0,npixtot-1
             IF ( ABS(tod(i,itod)/fmissing -1.) .LT. 1.e-5 ) THEN
                tod(i,itod) = fmiss_effct
                imissing = imissing + 1
             ENDIF
          ENDDO
          IF (imissing .GT. 0) THEN
             WRITE(*,'(a,1pe11.4)') 'blank value : ' ,fmissing
             WRITE(*,'(i7,a,f7.3,a,1pe11.4)') &
                  &           imissing,' missing pixels (', &
                  &           (100.*imissing)/npixtot,' %),'// &
                  &           ' have been set to : ',fmiss_effct
          ENDIF
       ENDIF
    ENDDO
    RETURN

  END SUBROUTINE input_tod_d
  !=======================================================================
  
  !**************************************************************************
  SUBROUTINE read_bintod_s(filename,tod,npixtot,ntods,firstpix,nullval,anynull)
  !**************************************************************************
    
    !=======================================================================
    !     Read a FITS file
    !
    !     slightly modified to deal with vector column (ie TFORMi = '1024E')
    !     in binary table       EH/IAP/Jan-98
    !
    !     This routine is used for reading TODS by anafast.
    !     Modified to start at a given pix numb OD & RT 02/02
    !     Modified to handle huge array (npix_tot > 2^32) OD & EH 07/02
    !      2002-07-08 : bugs correction by E.H. 
    !=======================================================================
    
    IMPLICIT NONE
    
    CHARACTER(LEN=*),               INTENT(IN)  :: filename
    INTEGER(I8B)   ,                INTENT(IN)  :: npixtot,firstpix
    INTEGER(I4B),                   INTENT(IN)  :: ntods
    REAL(SP), DIMENSION(0:,1:),     INTENT(OUT) :: tod
    REAL(SP),                       INTENT(OUT) :: nullval
    LOGICAL(LGT),                   INTENT(OUT) :: anynull
    
    INTEGER(I4B) :: status,unit,readwrite,blocksize,naxes(2),nfound, naxis
    INTEGER(I4B) :: group,nbuffer,i,npix_32,firstpix_32
    REAL(SP)     :: blank, testval
    REAL(DP)     ::  bscale,bzero
    CHARACTER(LEN=80) :: comment
    LOGICAL(LGT) :: extend
    INTEGER(I4B) :: nmove, hdutype
    INTEGER(I4B) :: column, frow, itod
    INTEGER(I4B) :: datacode, repeat, width
    
    INTEGER(I4B), PARAMETER :: maxdim=20 !number of columns in the extension
    INTEGER(I4B) :: nrows, tfields, varidat,felem
    CHARACTER(LEN=20) :: ttype(maxdim), tform(maxdim),tunit(maxdim), extname
 
    INTEGER(I8B) :: q,iq,npix_tmp,firstpix_tmp
    
    !-----------------------------------------------------------------------
    status=0
    
    unit = 150
    naxes(1) = 1
    naxes(2) = 1
    nfound = -1
    anynull = .FALSE.
    bscale = 1.0d0
    bzero = 0.0d0
    blank = -2.e25
    nullval = bscale*blank + bzero
    
    readwrite=0
    CALL ftopen(unit,filename,readwrite,blocksize,status)
    IF (status .GT. 0) CALL printerror(status)
    !     -----------------------------------------
    
    !     determines the presence of image
    CALL ftgkyj(unit,'NAXIS', naxis, comment, status)
    IF (status .GT. 0) CALL printerror(status)
    
    !     determines the presence of an extension
    CALL ftgkyl(unit,'EXTEND', extend, comment, status)
    IF (status .GT. 0) status = 0 ! no extension : 
    !     to be compatible with first version of the code
    
    IF (naxis .GT. 0) THEN ! there is an image
       !        determine the size of the image (look naxis1 and naxis2)
       CALL ftgknj(unit,'NAXIS',1,2,naxes,nfound,status)
       
       !        check that it found only NAXIS1
       IF (nfound .EQ. 2 .AND. naxes(2) .GT. 1) THEN
          PRINT *,'multi-dimensional image'
          PRINT *,'expected 1-D data.'
          STOP
       END IF
       
       IF (nfound .LT. 1) THEN
          CALL printerror(status)
          PRINT *,'can not find NAXIS1.'
          STOP
       ENDIF
       
       npix_32=naxes(1)
       
       CALL ftgkyd(unit,'BSCALE',bscale,comment,status)
       IF (status .EQ. 202) THEN ! BSCALE not found
          bscale = 1.0d0
          status = 0
       ENDIF
       CALL ftgkyd(unit,'BZERO', bzero, comment,status)
       IF (status .EQ. 202) THEN ! BZERO not found
          bzero = 0.0d0
          status = 0
       ENDIF
       CALL ftgkye(unit,'BLANK', blank, comment,status)
       IF (status .EQ. 202) THEN ! BLANK not found 
          ! (according to fitsio BLANK is integer)
          blank = -2.e25
          status = 0
       ENDIF
       nullval = bscale*blank + bzero
       
       !        -----------------------------------------
       
       group=1
       !firstpix = 1
       ! COMMENTED  
       
       npix_32 = npixtot
       ! in fits everything (pixel, row, ...) is counted from 1
       CALL ftgpve(unit,group,firstpix+1,npix_32,nullval,tod,anynull,status)
         
       ! if there are any NaN pixels, (real data)
       ! or BLANK pixels (integer data) they will take nullval value
       ! and anynull will switch to .true.
       ! otherwise, switch it by hand if necessary
       testval = 1.e-6 * ABS(nullval)
       DO i=0, npix_32-1
          IF (ABS(tod(i,1)-nullval) .LT. testval) THEN
             anynull = .TRUE.
             GOTO 111
          ENDIF
       ENDDO
111    CONTINUE
       
    ELSE IF (extend) THEN ! there is an extension

         nmove = +1
         CALL ftmrhd(unit, nmove, hdutype, status)
         !cc         write(*,*) hdutype
         
         IF (hdutype .NE. 2) THEN ! not a binary table
            STOP 'this is not a binary table'
         ENDIF
         
         ! reads all the keywords
         CALL ftghbn(unit, maxdim, &
              &        nrows, tfields, ttype, tform, tunit, extname, varidat, &
              &        status)
         
         IF (tfields .LT. ntods) THEN
            PRINT *,'found ',tfields,' tods in the file'
            PRINT *,'expected ',ntods
            STOP
         ENDIF
         
         ! finds the bad data value
         CALL ftgkye(unit,'BAD_DATA',nullval,comment,status)
         IF (status .EQ. 202) THEN ! bad_data not found
            nullval = -1.6375e30 ! default value
            status = 0
         ENDIF

         IF (npixtot .LT. nchunk_max) THEN
            
            DO itod = 1, ntods
               
               !parse TFORM keyword to find out the length of the column vector (repeat)
               CALL ftbnfm(tform(itod), datacode, repeat, width, status)
               frow = (firstpix)/repeat+1          ! 1 based 
               felem = firstpix-(frow-1)*repeat+1  ! 1 based 
               npix_32 = npixtot 
            
               !reads the columns
               column = itod
               CALL ftgcve(unit,column,frow,felem,npix_32,nullval, &
                    &        tod(0,itod),anynull,status)
            END DO
            
         ELSE
            
            q = (npixtot-1)/nchunk_max
            DO iq = 0,q
               IF (iq .LT. q) THEN
                  npix_tmp = nchunk_max
               ELSE
                  npix_tmp = npixtot - iq*nchunk_max
               ENDIF
               firstpix_tmp = firstpix + iq*nchunk_max
               npix_32 = npix_tmp
                              
               DO itod = 1, ntods
                  ! parse TFORM keyword to find out the length of the column vector
                  CALL ftbnfm(tform(itod), datacode, repeat, width, status)
                  frow = (firstpix_tmp)/repeat+1          ! 1 based 
                  felem = firstpix_tmp-(frow-1)*repeat+1  ! 1 based 

                  CALL ftgcve(unit,itod,frow,felem,npix_32,nullval, &
                       &      tod(firstpix_tmp-firstpix,itod),anynull,status)
               END DO
                              
            ENDDO
            
         ENDIF
         
      ELSE ! no image no extension, you are dead, man
         STOP ' No image, no extension'
      ENDIF
      
      ! close the file
      CALL ftclos(unit, status)
      
      ! check for any error, and if so print out error messages
      IF (status .GT. 0) CALL printerror(status)
      
      RETURN
      
    END SUBROUTINE read_bintod_s
    !=======================================================================
    
    ! ********************************************************************************
    SUBROUTINE read_bintod_d(filename,tod,npixtot,ntods,firstpix,nullval,anynull)
    ! ********************************************************************************

      !=======================================================================
      !     Read a FITS file
      !
      !     slightly modified to deal with vector column 
      !     in binary table       EH/IAP/Jan-98
      !
      !     Reads a double-precision array with precomputed plms, used by syn/anafast
      !                FKH/Apr-99
      !     Modified to start at a given pix numb OD & RT 02/02
      !     Modified to handle huge array (npix_tot > 2^32) OD & EH 07/02
      !      2002-07-08 : bugs correction by E.H. 
      !=======================================================================

      USE healpix_types

      IMPLICIT NONE

      CHARACTER(LEN=*),               INTENT(IN)  :: filename
      INTEGER(I8B)    ,               INTENT(IN)  :: npixtot,firstpix
      INTEGER(I4B)    ,               INTENT(IN)  :: ntods
      REAL(DP), DIMENSION(0:,1:),     INTENT(OUT) :: tod
      REAL(DP)        ,               INTENT(OUT) :: nullval
      LOGICAL(LGT)    ,               INTENT(OUT) :: anynull

      INTEGER(I4B) :: status,unit,readwrite,blocksize,naxes(2),nfound, naxis
      INTEGER(I4B) :: group,nbuffer,npix,i
      REAL(SP)     :: blank, testval
      REAL(DP)     ::  bscale,bzero
      CHARACTER(LEN=80) :: comment
      LOGICAL(LGT) :: extend
      INTEGER(I4B) :: nmove, hdutype
      INTEGER(I4B) :: column, frow, itod
      INTEGER(I4B) :: datacode, repeat, width

      INTEGER(I4B), PARAMETER :: maxdim=20 !number of columns in the extension
      INTEGER(I4B) :: nrows, tfields, varidat,felem,npix_32
      CHARACTER(LEN=20) :: ttype(maxdim), tform(maxdim), tunit(maxdim), extname
   
      INTEGER(I8B) :: q,iq,npix_tmp,firstpix_tmp

      !-----------------------------------------------------------------------
      status=0

      unit = 150
      naxes(1) = 1
      naxes(2) = 1
      nfound = -1
      anynull = .FALSE.
      bscale = 1.0d0
      bzero = 0.0d0
      blank = -2.e25
      nullval = bscale*blank + bzero

      readwrite=0
      CALL ftopen(unit,filename,readwrite,blocksize,status)
      IF (status .GT. 0) CALL printerror(status)
      !     -----------------------------------------

      !     determines the presence of image
      CALL ftgkyj(unit,'NAXIS', naxis, comment, status)
      IF (status .GT. 0) CALL printerror(status)

      !     determines the presence of an extension
      CALL ftgkyl(unit,'EXTEND', extend, comment, status)
      IF (status .GT. 0) status = 0 ! no extension : 
      !     to be compatible with first version of the code

      IF (naxis .GT. 0) THEN ! there is an image
         !        determine the size of the image (look naxis1 and naxis2)
         CALL ftgknj(unit,'NAXIS',1,2,naxes,nfound,status)

         !        check that it found only NAXIS1
         IF (nfound .EQ. 2 .AND. naxes(2) .GT. 1) THEN
            PRINT *,'multi-dimensional image'
            PRINT *,'expected 1-D data.'
            STOP
         END IF

         IF (nfound .LT. 1) THEN
            CALL printerror(status)
            PRINT *,'can not find NAXIS1.'
            STOP
         ENDIF

         npix=naxes(1)
!!$         if (npix .ne. npixtot) then
!!$            print *,'found ',npix,' plms'
!!$            print *,'expected ',npixtot
!!$            stop
!!$         endif

         CALL ftgkyd(unit,'BSCALE',bscale,comment,status)
         IF (status .EQ. 202) THEN ! BSCALE not found
            bscale = 1.0d0
            status = 0
         ENDIF
         CALL ftgkyd(unit,'BZERO', bzero, comment,status)
         IF (status .EQ. 202) THEN ! BZERO not found
            bzero = 0.0d0
            status = 0
         ENDIF
         CALL ftgkye(unit,'BLANK', blank, comment,status)
         IF (status .EQ. 202) THEN ! BLANK not found 
            ! (according to fitsio BLANK is integer)
            blank = -2.e25
            status = 0
         ENDIF
         nullval = bscale*blank + bzero

         !        -----------------------------------------

         group=1
         !firstpix = 1
         npix = npixtot
         ! in fits everything (pixel, row, ...) is counted from 1
         CALL ftgpvd(unit,group,firstpix+1,npix,nullval,tod,anynull,status)
         ! if there are any NaN pixels, (real data)
         ! or BLANK pixels (integer data) they will take nullval value
         ! and anynull will switch to .true.
         ! otherwise, switch it by hand if necessary
         testval = 1.e-6 * ABS(nullval)
         DO i=0, npix-1
            IF (ABS(tod(i,1)-nullval) .LT. testval) THEN
               anynull = .TRUE.
               GOTO 111
            ENDIF
         ENDDO
111      CONTINUE

      ELSE IF (extend) THEN ! there is an extension
         nmove = +1
         CALL ftmrhd(unit, nmove, hdutype, status)
         !cc         write(*,*) hdutype

         IF (hdutype .NE. 2) THEN ! not a binary table
            STOP 'this is not a binary table'
         ENDIF

         ! reads all the keywords
         CALL ftghbn(unit, maxdim, &
              &        nrows, tfields, ttype, tform, tunit, extname, varidat, &
              &        status)

         IF (tfields .LT. ntods) THEN
            PRINT *,'found ',tfields,' tods in the file'
            PRINT *,'expected ',ntods
            STOP
         ENDIF

         ! finds the bad data value
         CALL ftgkyd(unit,'BAD_DATA',nullval,comment,status)
         IF (status .EQ. 202) THEN ! bad_data not found
            nullval = -1.6375e30 ! default value
            status = 0
         ENDIF
         
         ! Read repeat value
         CALL ftbnfm(tform(1),datacode,repeat,width,status)

         IF (npixtot .LT. nchunk_max) THEN
            
            DO itod = 1, ntods

               !parse TFORM keyword to find out the length of the column vector
               CALL ftbnfm(tform(itod), datacode, repeat, width, status)
               frow = (firstpix)/repeat+1
               felem = firstpix-(frow-1)*repeat+1
               npix_32 = npixtot 
               
               !reads the columns
               column = itod
               CALL ftgcvd(unit,column,frow,felem,npix_32,nullval, &
                    &        tod(0,itod),anynull,status)
            END DO

         ELSE
            
            q = (npixtot-1)/nchunk_max
            DO iq = 0,q
               IF (iq .LT. q) THEN
                  npix_tmp = nchunk_max
               ELSE
                  npix_tmp = npixtot - iq*nchunk_max
               ENDIF
               firstpix_tmp = firstpix + iq*nchunk_max
               npix_32 = npix_tmp
               
               DO itod = 1, ntods
                  ! parse TFORM keyword to find out the length of the column vector
                  CALL ftbnfm(tform(itod), datacode, repeat, width, status)
                  frow = (firstpix_tmp)/repeat+1
                  felem = firstpix_tmp-(frow-1)*repeat+1 

                  CALL ftgcvd(unit,itod,frow,felem,npix_32,nullval, &
                       &        tod(firstpix_tmp-firstpix,itod),anynull,status)
               END DO
               
            ENDDO
            
         ENDIF
         
      ELSE ! no image no extension, you are dead, man
         STOP ' No image, no extension'
      ENDIF
      
      ! close the file
      CALL ftclos(unit, status)
      
      ! check for any error, and if so print out error messages
      IF (status .GT. 0) CALL printerror(status)
      
      RETURN
   
    END SUBROUTINE read_bintod_d
    ! ================================================================================
 
    !======================================================================================
    SUBROUTINE write_bintabh(tod, npix, ntod, header, nlheader,filename,firstpix,repeat)
    !======================================================================================

      ! =================================================================================
      !     Create a FITS file containing a binary table extension in the first extension
      !
      !     Designed to deal with Huge file, (n_elements > 2^31)
      !
      !     OPTIONNAL NEW PARAMETERS:
      !     firstpix : position in the file of the first element to be written (start at 0) 
      !                default value =0
      !                8 bytes integer
      !                if NE 0 then suppose that the file already exists
      !
      !     repeat   : lenght of vector per unit rows and colomns of the first binary extension
      !                default value = 12000 (\equiv 1 mn of PLANCK/HFI data)
      !                4 byte integer
      ! 
      !     OTHER PARAMETERS
      !     unchanged as compare to the standard write_bintab of the HEALPIX package except 
      !     npix which is a 8 bytes integer
      !
      !     Adapted from write_bintab
      !                                           E.H. & O.D. @ IAP 07/02
      !
      !     Requires a compilation in 64 bits of the CFITSIO 
      !     Note that the flag -D_FILE_OFFSETS_BITS=64 has to be added 
      !         (cf page CFITIO 2.2 User's guide  Chap 4, section 4-13)
      ! 
      ! 2002-07-08 : bugs correction by E.H. 
      !    (uniform use of firstpix_tmp, introduction of firstpix_chunk)
      !==========================================================================================

      USE healpix_types
      IMPLICIT NONE

      INTEGER(I8B)     , INTENT(IN)           :: npix
      INTEGER(I8B)     , INTENT(IN), OPTIONAL :: firstpix
      INTEGER(I4B)     , INTENT(IN), OPTIONAL :: repeat
      INTEGER(I4B)     , INTENT(IN)           :: ntod,nlheader
      REAL(SP)         , INTENT(IN), DIMENSION(0:npix-1,1:ntod) :: tod
      CHARACTER(LEN=80), INTENT(IN), DIMENSION(1:nlheader)      :: header
      CHARACTER(LEN=*),  INTENT(IN)           :: filename
      
      INTEGER(I4B) :: status,unit,blocksize,bitpix,naxis,naxes(1),repeat_tmp,repeat_fits
      INTEGER(I4B) :: group,fpixel,nelements,i,npix_32
      LOGICAL(LGT) :: simple,extend
      CHARACTER(LEN=80) :: svalue, comment, ch
      REAL(DP)     :: bscale,bzero

      INTEGER(I4B), PARAMETER :: maxdim = 20 !number of columns in the extension
      INTEGER(I4B)      :: nrows,tfields,varidat,np32
      INTEGER(I4B)      :: frow,felem,colnum,new,readwrite,width,datacode
      CHARACTER(LEN=20) :: ttype(maxdim), tform(maxdim), tunit(maxdim), extname
      CHARACTER(LEN=8)  :: date
      CHARACTER(LEN=10) :: fulldate
      CHARACTER(LEN=10) :: card
      CHARACTER(LEN=2)  :: stn
      INTEGER(I4B)      :: itn  

      INTEGER(I8B) :: q,iq,npix_tmp,firstpix_tmp,firstpix_chunk
      
      !-----------------------------------------------------------------------
      

      IF (.NOT. PRESENT(repeat) ) THEN 
         repeat_tmp = 12000 
      ELSE 
         repeat_tmp = repeat
      ENDIF
      IF (.NOT. PRESENT(firstpix) ) THEN 
         firstpix_tmp = 0 
      ELSE 
         firstpix_tmp = firstpix
      ENDIF
      
      status=0
      unit = 100
      blocksize=1
      

      ! create the new empty FITS file

      IF (firstpix_tmp .EQ. 0) THEN
         
         CALL ftinit(unit,filename,blocksize,status)
         
         ! -----------------------------------------------------
         ! Initialize parameters about the FITS image
         simple=.TRUE.
         bitpix=32     ! integer*4
         naxis=0       ! no image
         naxes(1)=0
         extend=.TRUE. ! there is an extension
         
         !     ----------------------
         !     primary header
         !     ----------------------
         !     write the required header keywords
         CALL ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)
         
         !     writes supplementary keywords : none
         
         !     write the current date
         CALL ftpdat(unit,status) ! format (dd/mm/yy)
         
         !     update the date (format ccyy-mm-dd)
         CALL DATE_AND_TIME(date)
         fulldate = date(1:4)//'-'//date(5:6)//'-'//date(7:8)
         comment = 'FITS file creation date ccyy-mm-dd'
         CALL ftukys(unit,'DATE',fulldate,comment,status)
         
         !     ----------------------
         !     image : none
         !     ----------------------
         
         !     ----------------------
         !     extension
         !     ----------------------
         
         !     creates an extension
         CALL ftcrhd(unit, status)

         !     writes required keywords
         nrows    = npix / repeat_tmp ! naxis1
         tfields  = ntod
         WRITE(ch,'(i8)') repeat_tmp
         tform(1:ntod) = TRIM(ADJUSTL(ch))//'E'
         
 !      tform(1:ntod) = '1024E'
         IF (npix .LT. repeat_tmp) THEN
            nrows = npix
            tform(1:ntod) = '1E'
         ENDIF
         ttype(1:ntod) = 'simulation'   ! will be updated
         tunit(1:ntod) = ''      ! optional, will not appear
         extname  = ''      ! optional, will not appear
         varidat  = 0

         CALL ftphbn(unit, nrows, tfields, ttype, tform, tunit, &
              &     extname, varidat, status)
         
         !     write the header literally, putting TFORM1 at the desired place
         DO i=1,nlheader
            card = header(i)
            IF (card(1:5) == 'TTYPE') THEN ! if TTYPE1 is explicitely given
               stn = card(6:6)
               READ(stn,'(i1)') itn
               ! discard at their original location:
               CALL ftmcrd(unit,'TTYPE'//stn,'COMMENT',status)  ! old TTYPEi and 
               CALL ftmcrd(unit,'TFORM'//stn,'COMMENT',status)  !     TFORMi
               CALL ftprec(unit,header(i), status)           ! write new TTYPE1
               comment = 'data format of field: 4-byte REAL'
               CALL ftpkys(unit,'TFORM'//stn,tform(1),comment,status) ! and write new TFORM1 right after
            ELSEIF (header(i).NE.' ') THEN
               CALL ftprec(unit,header(i), status)
            ENDIF
10          CONTINUE
         ENDDO
      
      ELSE
         ! The file already exists
         readwrite=1
         CALL ftopen(unit,filename,readwrite,blocksize,status)
         CALL ftmahd(unit,2,2,status)  ; ! 2 for first extension ; 2 for binary table

         CALL ftgkys(unit,'TFORM1',tform(1),comment,status)
         CALL ftbnfm(tform(1),datacode,repeat_fits,width,status)

         IF (repeat_tmp .NE. repeat_fits) THEN
            WRITE(*,*) 'WARNING routine write_bintabh'
            WRITE(*,*) 'Inexact repeat value. Use the one read in the file'
         ENDIF

      ENDIF 

      
      IF (npix .LT. nchunk_max) THEN ! data is small enough to be written in one chunk

         frow = (firstpix_tmp)/repeat_tmp + 1
         felem = firstpix_tmp-(frow-1)*repeat_tmp+1
         npix_32 = npix 
         
         DO colnum = 1, ntod
            CALL ftpcle(unit,colnum,frow,felem,npix_32,tod(0,colnum),status)
         END DO
         
      ELSE ! data has to be written in several chunks
        
         q = (npix-1)/nchunk_max
         DO iq = 0,q
            IF (iq .LT. q) THEN
               npix_tmp = nchunk_max
            ELSE
               npix_tmp = npix - iq*nchunk_max
            ENDIF
            firstpix_chunk = firstpix_tmp + iq*nchunk_max
            frow  = (firstpix_chunk)/repeat_tmp+1
            felem =  firstpix_chunk-(frow-1)*repeat_tmp+1 
            npix_32 = npix_tmp
            
            DO colnum = 1, ntod
               CALL ftpcle(unit,colnum,frow,felem,npix_32,tod(firstpix_chunk-firstpix_tmp,colnum),status)
            END DO
         ENDDO
            
      ENDIF

      ! ----------------------
      ! close and exit
      ! ----------------------
      
      ! close the file and free the unit number
      CALL ftclos(unit, status)
      
      ! check for any error, and if so print out error messages
      IF (status .GT. 0) CALL printerror(status)
      
      RETURN
    
    END SUBROUTINE write_bintabh
    ! ==============================================================================
    

END MODULE iotools













