      Integer*4 Function FUT_Clean_Catalog ( cats, cats_num, purge )

C------------------------------------------------------------------------
C
C    PURPOSE: Remove "deleted" catalog records and commandably purge
C	      duplicates from the input catalog record list.
C
C       REQUIREMENTS REFERENCE NUMBERS:
C
C       INPUT DATA:
C
C       OUTPUT DATA:
C
C    AUTHOR: Rob Kummerer
C            ST Systems Corporation
C            August 2, 1989
C
C    INVOCATION: status = FUT_Clean_Catalog ( cats,
C					      cats_num,
C					      purge )	      	      
C
C    INPUT PARAMETERS:
C
C	CATS(50)	RECORD	List of catalog records to be purged.
C	CATS_NUM	I*4	Total number of records.
C	PURGE		I*2	Flag turning ON/OFF purge.
C
C    OUTPUT PARAMETERS: 
C	STATUS		I*4	Success status.
C	CATS(50)	RECORD	Purged list of catalog records.
C	CATS_NUM	I*4	Total number of records after purge.
C
C    SUBROUTINES CALLED: None
C
C    COMMON VARIABLES USED: None
C
C    INCLUDE FILES:
C	CUT_Params
C	FUT_Params
C
C----------------------------------------------------------------------

	Implicit None

	Include		'(CUT_Params)'
	Include		'(FUT_Params)'

	Dictionary 'CCM_CME_Catalog_Entry'

	Integer		*2 	cats_num
	Record /CCM_CME_CATALOG_ENTRY/ cats(cats_num)
	Integer		*2	purge

	Record /CCM_CME_CATALOG_ENTRY/ purge_cats(CUT$_Max_Segments)
	Integer		*2	purge_num

	Integer		*4	i
	Integer		*4	j
	Character	*32	iext
	Character	*32	jext
	Character	*8	ivers
	Character	*8	jvers
	Logical		*1	idel
	Logical		*1	jdel

	External	FUT_Normal

C
C Perform the "catalog purge" if requested.  Here we assume that the highest
C version is the most recent version.  Only this file is to be kept.
C
	If (purge .Eq. fac_present) Then

	   Do i=1,cats_num

	      idel = cats(i).deletion_flag

	      If (.Not. idel) Then

	         iext = cats(i).filename_extension
	         iext = iext(1:index(iext,' ')-1)
	         ivers = cats(i).version_number
	         ivers = ivers(1:index(ivers,' ')-1)

	         Do j=i+1,cats_num

		    jdel = cats(j).deletion_flag

		    If (.Not. jdel) Then

		       jext = cats(j).filename_extension
		       jext = jext(1:index(jext,' ')-1)
		       jvers = cats(j).version_number
		       jvers = jvers(1:index(jvers,' ')-1)

		       If (iext .Eq. jext) Then

			  If (ivers .Gt. jvers) Then
		             cats(j).deletion_flag = .True.
		          Else
		             cats(i).deletion_flag = .True.
		          End If

		       End If

		    End If

	         End Do

	      End If

	   End Do

	End If

C
C First remove catalog records that have been marked as deleted.
C This includes those really marked as deleted and those marked
C as deleted from the above purge.
C
	purge_num = 0

	Do i=1,cats_num

	   If (.Not. cats(i).deletion_flag) Then
	      purge_num = purge_num + 1
	      purge_cats(purge_num) = cats(i)
	   End If

	End Do

C
C Return the purged catalog list.
C
	Do i=1,purge_num
	   cats(i) = purge_cats(i) 
	End Do

	cats_num = purge_num


	FUT_Clean_Catalog = %Loc(FUT_Normal)

	Return
	End
