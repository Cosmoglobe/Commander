pro DDL_abbrev, Field_Name

   common DDL_abbrev, Field_abbrev

	if N_elements( Field_abbrev ) LT 2 then begin

	   Field_abbrev = $
	      [	["GALACTIC_LONGITUDE",	"GALACTIC_LONG"	]		,$
		["GALACTIC_LATITUDE",	"GALACTIC_LAT"	]		,$
		["ECLIPTIC_LONGITUDE",	"ECLIPTIC_LONG"	]		,$
		["ECLIPTIC_LATITUDE",	"ECLIPTIC_LAT"	]		,$
		["TERR_LATITUDE",	"TERR_LAT"	]		,$
		["TERR_LONGITUDE",	"TERR_LONG"	]		,$
		["PROJECTED_GEOCENTRIC_VELOCITY",   "PROJ_GEO_VELOC"	],$
		["PROJECTED_HELIOCENTRIC_VELOCITY", "PROJ_HELIO_VELOC"	] ]

		Field_abbrev = transpose( Field_abbrev )
	   endif

	w = where( Field_abbrev(*,0) EQ Field_Name, Nf )

	if (Nf GT 0) then  Field_Name = Field_abbrev(w(0),1)
return
end
