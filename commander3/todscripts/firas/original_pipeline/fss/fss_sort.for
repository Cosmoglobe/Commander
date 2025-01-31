      integer*4 function fss_sort(sky_buff, sci_recs, sky_recs,
     .                            tot_recs, chan)

C------------------------------------------------------------------------
C    PURPOSE: Sort the ifgs from the raw science and the input 
C             skymap. Sort the data on pixel, scan mode, ICAL temp, dihedral
C             temp, and time. Return sorted data and total records.
C
C    AUTHOR: D. Bouler, STX, Jan , 1991
C
C    INVOCATION:  status = fss_sort(sky_buff, sci_recs, sky_recs,
C                                   tot_recs, chan
C
C    INPUT PARAMETERS:     ARRAY sky_buff    records from science & skymap
C                          I*4   sci_recs    number of records in sci_buff
C                          I*4   sky_recs    number of records in sky_buff
C                          I*4   chan        channel
C
C    OUTPUT PARAMETERS:    ARRAY sky_buff    sorted array of FDQ and sky recs
C                          I*4   tot_recs    Total number of ifgs per channel
C
C    SUBROUTINES CALLED:    sor$begin_sort
C                           sor$release_rec
C                           sor$merge_sort
C                           sor$return_rec
C                           sor$end_sort
C                           lib$movc3
C
C    COMMON VARIABLES USED: None.
C 
C    INCLUDE FILES:         fss_include
C
C----------------------------------------------------------------------
C                      PDL for fss_sort
c                      D. BOULER, Jan , 1991
c
c  get total number of records
c
c  sort records in sky_buff array by pixel, scan mode, ical temp, time
c
c  fss_sort = status
c
c  return 
c  end (pdl)
c
c-----------------------------------------------------------------------------
C Changes:
C
C  4-20-95  L.Rosen - Sort on ICAL bin number (which is temporarily
C           being stored in the ORIGINAL_PIXEL field) instead of the
C           ICAL temperature.  Temperature is a float and not exact
C           for the same bin, and so time ordering was lost. This will
C           now produce time sorted records within each group.
C
C  SPR 12197 - Modifications to implement /DBINS qualifier in place
C              of /MAX_DIHED,  K.Jensen, HSTX, 23-May-1995
C              Sort on ICAL/Dihedral combined bin number. Combined bin
C              number = 100 * ICAL bin number + Dihedral bin number.
C
C-----------------------------------------------------------------------------

      implicit none
c
c Include files
c
      include  '(fss_include)'
c
c Function declarations
c
      integer*4 sor$begin_sort
      integer*4 sor$release_rec
      integer*4 sor$return_rec
      integer*4 sor$end_sort
      integer*4 sor$sort_merge
      integer*4 lib$movc3
C
C Variable declarations
C
      integer*2	sort_keys(21)           !Provide information on sorting 

      integer*4 tot_recs                !total recs from science and skymap
      integer*4 chan                    !channel
      integer*4 i,j,k,l                 !loop counters
      integer*4 status                  !return status
      integer*4 sci_recs                !number of recs input from full science
      integer*4 sky_recs                !number of recs from input skymap

c Declare externals
c
      external fss_normal
      external DSC$K_DType_B               !SOR type definition (byte)
      external DSC$K_DType_F               !SOR type definition (real*4)
      external DSC$K_DType_W               !SOR type definition (word)
      external DSC$K_DType_L               !SOR type definition (longword)
      external DSC$K_DType_Q               !SOR type definition (quadword)
c
c Dictionary and record declarations
c                                             
      dictionary 'fss_sssky'

      record /fss_sssky/     sky_buff(fss_max_buff)
c*
c******  Begin code *******************************************
c*
      fss_sort = %loc(fss_normal)         
c
c Get total number of ifgs.
c
      tot_recs = sci_recs + sky_recs
c*
c*    Intialize Sort_Keys for the SOR$ routines
c*
      Sort_Keys(01) = 5                     !  Total of five keys
      Sort_Keys(02) = %loc(DSC$K_Dtype_L)   !Key 1 = Pixel
      Sort_Keys(03) = 0                     !  Ascending sort
      Sort_Keys(04) = 9                     !  Start byte of key
      Sort_Keys(05) = 4                     !  Length in bytes of key
      Sort_Keys(06) = %loc(DSC$K_DType_B)   !Key 2 = MTM_Scan_Speed
      Sort_Keys(07) = 0                     
      Sort_Keys(08) = 13
      Sort_Keys(09) = 1
      Sort_Keys(10) = %loc(DSC$K_Dtype_B)  !Key 3 = MTM_Scan_Length
      Sort_Keys(11) = 0                    
      Sort_Keys(12) = 14
      Sort_Keys(13) = 1
      Sort_Keys(14) = %loc(DSC$K_DType_W)  !Key 4 = ICAL/Dihedral Temperature Bin
      Sort_Keys(15) = 0                    
      Sort_Keys(16) = 62
      Sort_Keys(17) = 2
      Sort_Keys(18) = %loc(DSC$K_DType_Q)  !Key 5 = Time
      Sort_Keys(19) = 0                    
      Sort_Keys(20) = 0
      Sort_Keys(21) = 8
c
c Initialize the sort process
c
      status = SOR$Begin_Sort(Sort_Keys, fss_byte_recsize)
      if (.not. status) call lib$signal(%val(status))
c
c Release records of sky_buff array to sort process
c
      i = 1
      do while (status .and. (i .le. (tot_recs)))
         call LIB$Movc3 ( fss_byte_recsize, sky_buff(i), %ref(sort_rec))
         status = SOR$Release_Rec ( sort_rec )
         if (.not. status) call lib$signal(%val(status))
         i = i + 1
      end do
c
c Now sort the combined array
c      
      if (status) then
          status = SOR$Sort_Merge ()
          if (.not. status) call lib$signal(%val(status))
      endif
c
c Place the combined array back into the sky_buff array.
c Simultaneously copy each records pixel number into its ORIGINAL_PIXEL
c field (used in neighbors).
c
      i = 1
      do while (status .and. (i .le. (tot_recs)))
         status = SOR$Return_Rec (sort_rec)
         if (.not. status) call lib$signal(%val(status))
         call lib$movc3 (fss_byte_recsize, %ref(sort_rec), sky_buff(i))
         sky_buff(i).original_pixel = sky_buff(i).pixel_no
         i = i + 1
      end do
c
c Close out the sort process
c
      if (status) then
          status = SOR$End_Sort ()
          if (.not. status) call lib$signal(%val(status))
      endif

      fss_sort = status

      return
      end
