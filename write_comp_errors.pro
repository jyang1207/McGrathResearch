; procedure to filter the bad components from a summary.txt output of *sumgalfit.txt
; and write them to a *summary_all_errors.txt file in the same directory
;
; positional parameters:
; infile - the full path filename of the component summary info (*summary.txt)
;
; named parameters:
; DIST_FROM_CENTER (e.g. DIST = 20.0) the maximum distance from (300,300) of
;			the center of the component model before being called an error
; SERSIC_INDEX (e.g. SERS = 0.5) the minimum sersic index of component before being called an error
;			
pro write_comp_errors, infile, DIST_FROM_CENTER = distLimit, SERS_LOWER_LIMIT = sersLowLimit, SERS_UPPER_LIMIT = sersHighLimit

	; if an error is detected anywhere in the procedure
	; stores the detection in file_error variable
	catch, file_error

	; read the summary.txt file, full path name is in infile parameter
	; store arrays for each field in 14 named variables (id, ts, etc.)
	readcol,infile,id,ts,age,cam,fil,px,py,mag,rad,ser,ba,ang,$
		SKIPLINE=2,FORMAT="A,A,F,A,A"

	; stop execution if an error has been detected
	if file_error NE 0 then STOP

	; find the position in the string just after any leading ../
	; so that the position of the . preceding the file extension can be found
	start_index = 0
	while ($
		( strlen( strmid(infile,start_index,strpos(infile,'.',start_index)) ) EQ 0 ) $
		and ( start_index LT strlen(infile) ) $
	) $
	do start_index = start_index + 3
	if not( start_index LT strlen(infile) ) then print,"failed to find the file extension for naming the output file"

	; use above info to define output filename
	out_filename = (strmid(infile,0,strpos(infile,".",start_index))+'_all_errors.txt') 

	; array of all records that start with a *, which indicates sersic or component seperation error
	; errors = where(strmid(id,0,1) EQ "*", errorCount)

	; open a file for writing the output		
	get_lun,outfile
	openw,outfile,out_filename

	; write the number of components with leading * to the output file
	; printf,outfile,errorCount," errors detected by sersic index and component separation (components with leading *s)"
	; printf,outfile,''

	; based on named param, write number and list of components too far from (300,300)
	if not(n_elements(distLimit) eq 0) then begin
		dx = (300-px)
		dy = (300-py)
		distErrors = where( ( ((dx^2 + dy^2)^(0.5)) GT float(distLimit) ), distCount )
		printf,outfile,distCount," components with (x, y) position greater than ",distLimit," from (300, 300)"
		for i=0,distCount-1 do $
			printf,outfile,id[distErrors[i]]," time:",ts[distErrors[i]],$
					" camera:",cam[distErrors[i]]," type:",typ[distErrors[i]],$
					" with a position of (",px[distErrors[i]],",",py[distErrors[i]],")"
		printf,outfile,''
	endif
	
	; based on named param, write number and list of components with sersic index too low
	if not(n_elements(sersLowLimit) eq 0) then begin
		sersLowErrors = where( ( ser LT float(sersLowLimit) ), sersLowCount )
		printf,outfile,sersLowCount," components with sersic index less than ",sersLowLimit
		for i=0,sersLowCount-1 do $
			printf,outfile,id[sersLowErrors[i]]," time:",ts[sersLowErrors[i]],$
					" camera:",cam[sersLowErrors[i]]," type:",typ[sersLowErrors[i]],$
					" with a sersic index of ",ser[sersLowErrors[i]]
		printf,outfile,''
	endif

	; based on named param, write number and list of components with sersic index too high
	if not(n_elements(sersHighLimit) eq 0) then begin
		sersHighErrors = where( ( ser GT float(sersHighLimit) ), sersHighCount )
		printf,outfile,sersHighCount," components with sersic index less than ",sersHighLimit
		for i=0,sersHighCount-1 do $
			printf,outfile,id[sersHighErrors[i]]," time:",ts[sersHighErrors[i]],$
					" camera:",cam[sersHighErrors[i]]," type:",typ[sersHighErrors[i]],$
					" with a sersic index of ",ser[sersHighErrors[i]]
		printf,outfile,''
	endif
	
	; close output file
	free_lun,outfile

	; print the location of the written output file
	print,"file ",out_filename," written"

end
