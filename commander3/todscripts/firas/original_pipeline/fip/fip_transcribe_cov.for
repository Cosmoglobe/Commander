	Integer*4  Function  FIP_TRANSCRIBE_COV  ( fcc_covar, covar, index_low,
	1                                          index_high )

C  Transcribe the converted FCC covariance matrix to the FIP covariance matrix,
C  in accordance with the software requirements reproduced below.
C
C   Larry P. Rosen, Hughes STX, 26 July 1993.
C
C    Input:
C	Real*8        fcc_covar (522,522)  ! converted FCC Covariance Matrix
C	Integer*2     index_low            ! Low frequency cut-off array index
C	Integer*2     index_high           ! High frequency cut-off array index
C    Output:
C	Real*8        covar (368,368)      ! FIP converted covariance matrix
C
C  PDL - software requirements
C
C    Let BIN_1 designate the first frequency bin and LAST_BIN designate the
C    last.
C       a.  For BIN_1 <= [i and j] <= LAST_BIN, FCC_COVAR (i,j) shall be
C       transcribed to COVAR (i + 1 - BIN_1, j + 1 - BIN_1) (actually necessary
C       only for bins where j >= i since the upper triangular part of the matrix
C       only will be stored).
C       b.  For BIN1 <= i <= LAST_BIN and BIN_1 + 265 <= j <= LAST_BIN + 265,
C       FCC_COVAR (i,j) shall be transcribed to
C       COVAR (i + 1 - BIN_1, j - (76 + BIN_1)).
C       c.  For 265 + BIN_1 <= [i and j] <= 265 + LAST_BIN;  FCC_COVAR (i,j)
C       shall be transcribed to COVAR (i- (76 + BIN_1), j - (76 + BIN_1))
C       (also actually necessary only for bins where j >= i).
C       d.  For BIN_1 <= i <= LAST_BIN and 258 <= j <= 265, FCC_COVAR (i,j)
C       shall be transcribed to COVAR (i + 1 - BIN_1, j - 77).
C       e.  For 258 <= i <= 265 and BIN_1 <= j <= LAST_BIN, FCC_COVAR (i,j)
C       shall be transcribed to COVAR (i - 77, j + 1 - BIN_1).
C       f.  For 258 <= i <= 265 and 258 <= j <= 265, FCC_COVAR (i,j) shall be
C       transcribed to COVAR (i  - 77,j - 77) (also actually necessary only for
C       bins where j >= i).
C       g.  For 265 + BIN_1 <= i <= 265 + LAST_BIN and BIN_1 <= j <= LAST_BIN,
C       FCC_COVAR (i,j) shall be transcribed to COVAR (i - (76 + BIN_1),
C       j + 1 - BIN_1) (transcription of this piece is not technically
C       necessary since only the upper triangular portion of the matrix will
C       be stored).
C       h.  For 265 + BIN_1 <= i <= 265 + LAST_BIN and 258 <= j <= 265,
C       FCC_COVAR (i,j) shall be transcribed to COVAR (i - (76 + BIN_1),j - 77)
C       (transcription of this piece is also not technically necessary).
C       i.  For 258 <= i <= 265 and 265 + BIN_1 <= j <= 265 + LAST_BIN,
C       FCC_COVAR (i,j) shall be transcribed to COVAR (i - 77,j - (76 + BIN_1))
C       (transcription of this piece is also not technically necessary).
C

	Implicit None

C Passed parameters

	Real*8        fcc_covar (522,522)  ! FCC Covariance Matrix
	Real*8        covar (368,368)      ! FIP converted covariance matrix
	Integer*2     index_low            ! Low frequency cut-off array index
	Integer*2     index_high           ! High frequency cut-off array index

C Local

	Integer*2	row             ! row of covariance matrix
	Integer*2	col             ! column of covariance matrix

C External

	External	fip_normal

C Begin

	FIP_TRANSCRIBE_COV = %loc (fip_normal)

	Do col = index_low, index_high
	   Do row = index_low, index_high
	      covar (row+1-index_low, col+1-index_low) = fcc_covar (row, col)
	   Enddo
	Enddo
	Do col = index_low+265, index_high+265
	   Do row = index_low, index_high
	      covar (row+1-index_low, col-(76+index_low)) = fcc_covar (row, col)
	   Enddo
	Enddo
	Do col = 265+index_low, 265+index_high
	   Do row = 265+index_low, 265+index_high
	      covar (row-(76+index_low), col-(76+index_low)) =
	1        fcc_covar (row, col)
	   Enddo
	Enddo
	Do col = 258, 265
	   Do row = index_low, index_high
	      covar (row+1-index_low, col-77) = fcc_covar (row, col)
	   Enddo
	Enddo
	Do col = index_low, index_high
	   Do row = 258, 265
	      covar (row-77, col+1-index_low) = fcc_covar (row, col)
	   Enddo
	Enddo
	Do col = 258, 265
	   Do row = 258, 265
	      covar (row-77, col-77) = fcc_covar (row, col)
	   Enddo
	Enddo
	Do col = index_low, index_high
	   Do row = 265+index_low, 265+index_high
	      covar (row-(76+index_low), col+1-index_low) = fcc_covar (row,col)
	   Enddo
	Enddo
	Do col = 258, 265
	   Do row = 265+index_low, 265+index_high
	      covar (row-(76+index_low), col-77) = fcc_covar (row, col)
	   Enddo
	Enddo
	Do col = 265+index_low, 265+index_high
	   Do row = 258, 265
	      covar (row-77, col-(76+index_low)) = fcc_covar (row, col)
	   Enddo
	Enddo
	Return
	End
