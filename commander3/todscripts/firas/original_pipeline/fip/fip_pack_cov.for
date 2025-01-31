	Integer*4  Function  FIP_PACK_COV  ( covar, real_real, real_imag,
	1                                    imag_imag )

C Pack covariance matrix COVAR into the three vector form for storage in
C the FIP record.
C
C Larry P. Rosen, Hughes STX, 26 July 1993.
C
C    Input:
C	Real*8		covar (368,368)       ! FIP converted covariance matrix
C    Output:
C	Real*8		real_real (17766)     ! Vector of covar matrix
C	Real*8		real_imag (33840)     ! Vector of covar matrix
C	Real*8		imag_imag (16290)     ! Vector of covar matrix
C
C    Upper triangular part of symmetric covariance matrix, stored in records
C    1 - 253, each with up to 270 values, row x col.
C       Real - real (17766 values) stored in records 1 - 66;
C          upper half of symmetric submatrix 1-188 x 1-188.
C       Real - imag (33840 values) stored in records 67 - 192;
C          entirety of rectangular submatrix 189-368 x 1-188.
C       Imag - imag (16290 values) stored in records 193 - 253;
C          upper half of symmetric submatrix 189-368 x 189-368.

	Implicit None

C Passed Parameters

	Real*8		covar (368,368)       ! FIP converted covariance matrix
	Real*8		real_real (17766)     ! Vector of covar matrix
	Real*8		real_imag (33840)     ! Vector of covar matrix
	Real*8		imag_imag (16290)     ! Vector of covar matrix

C Local

	Integer*2	col, row        ! column and row of covariance matrix
	Integer*4	rr_element      ! number of elements put in matrix R-R
	Integer*4	ri_element      ! number of elements put in matrix R-I
	Integer*4	ii_element      ! number of elements put in matrix I-I

C External

	External	fip_normal


C Begin

	FIP_PACK_COV = %loc (fip_normal)

	rr_element = 0
	ri_element = 0
	ii_element = 0
	Do row = 1, 188
	   Do col = row, 188
	      rr_element = rr_element + 1
	      real_real (rr_element) = covar (row, col)
	   Enddo
	Enddo
	Do row = 1, 188
	   Do col = 189, 368
	      ri_element = ri_element + 1
	      real_imag (ri_element) = covar (row, col)
	   Enddo
	Enddo
	Do row = 189, 368
	   Do col = row, 368
	      ii_element = ii_element + 1
	      imag_imag (ii_element) = covar (row, col)
	   Enddo
	Enddo
	Return
	End
