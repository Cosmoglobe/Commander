function DDL_ST, DDL_Lun, Lun_out, DDL_NAME=DDL_name, ABBREVIATE=abbrev, $
								NOFILE=nofile
;+
; NAME:
;	DDL_ST
; PURPOSE:
;	Convert a DDL record-structure definition file
;			into IDL structure(s) code, packaged as a function.
; CATEGORY:
;			Structures
; CALLING SEQUENCE:
;			x = DDL_st( DDL="DDL_name" )
;			x = DDL_st( DDL="DDL_name", /nof )
;			x = DDL_st( DDL="DDL_name", /abbrev )
;	example:
;			print, ddl_st( ddl="FDQ_SDF", /abbrev )
; INPUTS/OUTPUTS:
;   keywords:
;		DDL_NAME = string giving DDL file name (type .DDL is assumed )
;   optional:
;		/ABBREV performs Field name abbreviations
;			as defined by DDL_abbrev.pro
;		/NOFILE to print results only to terminal.
;
;   inputs used recursively:
;		DDL_Lun = Logical unit number for reading DDL file.
;		Lun_out = Logical unit number for writing .PRO code.
;		ABBREV = keyword to indicate use of abbreviations.
;
;   function returns a string which is not important to user
;		(it is just the final return from chain of recursion).
; SIDE EFFECTS:
;		"DDL_name".DDL file is read,
;		"DDL_name"_ST.PRO is created (IDL function code).
; PROCEDURE:
;		Strategy is to call DDL_st recursively.
;  The DDL file (given by "DDL_name".DDL) is read one line at a time,
;
;  if declaration STRUCTURE is encountered,
;  then DDL_ST is called recursively, which then builds the structure
;   until the declaration END STRUCTURE is encountered in DDL file,
;    then DDL_ST returns the complete IDL structure definition,
;	and the result is printed to file "DDL_name"_struct.pro ,
;		or to terminal screen if /NOFILE is specified.
;
;  otherwise the variable type and size is used to construct IDL definition,
;	(calling DDL_abbrev,Var_Name for name abbreviations if /ABBREV),
;	which is then concatenated with the other IDL definitions,
;   	and returned when declaration END STRUCTURE is encountered in DDL file.
;
; EXTERNAL CALLS:
;		pro DDL_ABBREV		(should be with this routine)
;		function GET_WORDS	(from [varosi.idl.lib]VAROSI.TLB) 
; MODIFICATION HISTORY:
;	written 1990 Frank Varosi STX @ NASA/GSFC as DDL_STRUCT.PRO
;       October 2, 1991 K. Jensen STX - Changed to DDL_ST.PRO. Effect of
;           the change is to change the output IDL structure procedure
;           name from "DDL_Name"_STRUCT.PRO to "DDL_Name_ST.PRO.
;       
	if ( N_elements( DDL_name ) EQ 1 ) AND $	;do this just once.
	   ( N_elements( DDL_Lun ) NE 1 ) then begin

		DDL_name = strupcase( DDL_name )
		DDL_file = DDL_name + ".DDL"

		openr, DDL_Lun, DDL_file, /get_Lun, /SHARE
		message, "reading "+DDL_file, /CONTIN

		func_struct = DDL_name + "_st"

		func_code =  [  "function " + func_struct + ", Nstruct"      ,$
				" "					     ,$
				"  common "  + func_struct + ", defined"     ,$
				" "					     ,$
				"  if N_elements( defined ) NE 1 then begin" ,$
				" " ]


		func_end =  [   "       defined = 1 "			     ,$
				"    endif"				     ,$
				" "					     ,$
			"  if N_elements( Nstruct ) NE 1 then Nstruct = 1"   ,$
			" "						     ,$
			"return, replicate( {" + DDL_name + "}, Nstruct )"   ,$
			"end" ]

		Nline = N_elements( func_code )

		if NOT keyword_set( nofile ) then begin

			file_struct = func_struct + ".PRO"
			openw, Lun_out, file_struct, /get_Lun
			message, "creating file " + file_struct, /CONTIN

			for i=0,Nline-1 do printf, Lun_out, func_code(i)

		  endif else begin

			print," "
			for i=0,Nline-1 do print, func_code(i)
		   endelse
	   endif

	DDL_rec = ""
	DDL_rec2 = ""

	readf, DDL_Lun, DDL_rec

	while (strpos( DDL_rec, "." ) LT 0) do begin
		readf, DDL_Lun, DDL_rec2
		DDL_rec = DDL_rec + DDL_rec2
	  endwhile

	DDL_rec = strupcase( DDL_rec )
	struct_Tags = ""

   while ( strpos( DDL_rec, "END STRUCTURE" ) LT 0 ) AND $
         ( strpos( DDL_rec, "END VARIANT" ) LT 0 ) do begin

	if ( strpos( DDL_rec, "STRUCTURE" ) GT 0 ) then begin

		DDL_rec = strmid( DDL_rec, 0, strlen(DDL_rec)-1 )
		DDL_words = get_words( DDL_rec )
		struct_name = DDL_words(0)
		struct_namLc = strlowcase( struct_name )

		structure = "       " + struct_namLc + " = { " + struct_name

		struct_def = DDL_ST( DDL_Lun, Lun_out, ABBREV=abbrev )

		blanks = replicate( 32B, 80 )
		Lpad = strpos( structure, "{" ) + 1
		struct_def = string( blanks(0:Lpad) ) + struct_def

		structure = [ structure, struct_def ]

		Lpad = strlen( structure )
		Lpad = (max( Lpad ) - Lpad) < 70
		Nstr = N_elements( structure )

		for i=0,Nstr-2 do structure(i) = structure(i) + $
					string( blanks(0:Lpad(i)) ) + ", $"

		structure(Nstr-1) = structure(Nstr-1) + $
					string( blanks(0:Lpad(i)+4) ) + "}"

		if N_elements( Lun_out ) EQ 1 then begin

			for i=0,Nstr-1 do printf, Lun_out, structure(i)
			printf, Lun_out, " "

		  endif else begin

			for i=0,Nstr-1 do print, structure(i)
			print, " "
		   endelse

		if (strpos( DDL_rec, "ARRAY" ) GT 0) then begin

			w = where( DDL_words EQ "ARRAY" )
			size = DDL_words(w(0)+1)

			struct_Tags = [ struct_Tags, struct_name + $
			  ": replicate( " + struct_namLc + ", " + size + " )" ]

		  endif else  struct_Tags = $
			      [ struct_Tags, struct_name + ":" + struct_namLc ]

	 endif else if ( strpos( DDL_rec, "VARIANTS") GT 0 ) then begin

		message,"encountered VARIANTS, using only first VARIANT",/CONT
		readf, DDL_Lun, DDL_rec

		struct_Tags = [ struct_Tags , $
				DDL_ST( DDL_Lun, Lun_out, ABBREV=abbrev ) ]

		while ( strpos( DDL_rec, "END VARIANTS" ) LT 0 ) do begin
			readf, DDL_Lun, DDL_rec
			DDL_rec = strupcase( DDL_rec )
		  endwhile

	  endif else begin

		DDL_rec = strmid( DDL_rec, 0, strlen(DDL_rec)-1 )
		DDL_words = get_words( DDL_rec )
		Var_Name = DDL_words(0)
		Var_Type = DDL_words( N_elements( DDL_words )-1 )

		if keyword_set( abbrev ) then  DDL_abbrev, Var_Name

		if (strpos( DDL_rec, "TEXT" ) GT 0) then begin

			w = where( DDL_words EQ "SIZE" )
			size_text = DDL_words(w(0)+2)

			if (size_text EQ '1') then		$
				variable = Var_Name + ": 0B"	$
			else $
			  variable = Var_Name + ": bytarr( " + size_text + " )"

		 endif else if (strpos( DDL_rec, "ARRAY" ) GT 0) then begin

			w = where( DDL_words EQ "ARRAY" )
			size_array = DDL_words(w(0)+1)

			CASE Var_Type OF

			"BYTE":		variable = Var_Name + $
					": bytarr( " + size_array + " )"

			"WORD":		variable = Var_Name + $
					": intarr( " + size_array + " )"

			"LONGWORD":	variable = Var_Name + $
					": Lonarr( " + size_array + " )"

			"F_FLOATING":	variable = Var_Name + $
					": fltarr( " + size_array + " )"

			"D_FLOATING":	variable = Var_Name + $
					": dblarr( " + size_array + " )"

			"COMPLEX":	variable = Var_Name + $
					": complexarr( " + size_array + " )"

			else:		variable = ""

			ENDCASE

		  endif else begin

			CASE Var_Type OF

			"BYTE":		variable = Var_Name + ": 0B"

			"WORD":		variable = Var_Name + ": 0"

			"LONGWORD":	variable = Var_Name + ": 0L"

			"F_FLOATING":	variable = Var_Name + ": 0.0"

			"D_FLOATING":	variable = Var_Name + ": 0.D0"

			"COMPLEX":	variable = Var_Name + ": complex(0)"

			else:		variable = ""

			ENDCASE
	           endelse

		struct_Tags = [ struct_Tags , variable ]

	   endelse

	readf, DDL_Lun, DDL_rec

	while (strpos( DDL_rec, "." ) LT 0) do begin
		readf, DDL_Lun, DDL_rec2
		DDL_rec = DDL_rec + DDL_rec2
	   endwhile

	DDL_rec = strupcase( DDL_rec )

	if strpos( DDL_rec, "END RECORD" ) GT 0 then begin

		free_Lun, DDL_Lun
		Nline = N_elements( func_end )

		if N_elements( Lun_out ) EQ 1 then begin

			for i=0,Nline-1 do printf, Lun_out, func_end(i)
			free_Lun, Lun_out

		  endif else for i=0,Nline-1 do print, func_end(i)

		message,"finished IDL structure code for "+DDL_name,/CONTIN
		Nstr = N_elements( struct_Tags )
		return, struct_Tags(Nstr-1)
	   endif

   endwhile

	w = where( strlen( struct_Tags ) )

return, struct_Tags(w)
end
