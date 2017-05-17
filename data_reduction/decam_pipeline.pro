pro decam_pipeline, recipe, PROGRAM=program, FILTER=filter, OVERWRITE=overwrite, TILE=tile, DEBUG=debug, STANDARD=standard, AHEAD=ahead, FWHM=fwhm, DITHER=dither, DITHER_EXCLUDE=dither_exclude, ALIGN_FWHM=align_fwhm, ZP_PHOT=zp_phot

do_program= (n_elements(program) GT 0) ? strupcase(program) : '20130808'
do_tile= (n_elements(tile) GT 0) ? tile : '1'
do_filter= (n_elements(filter) GT 0) ? filter : 'g'
do_dither_exclude= (n_elements(dither_exclude) GT 0) ? dither_exclude : 'NONE'

do_standard= (n_elements(standard) GT 0) ? standard : 'SDSS'
do_n_chip= (do_program EQ '20130808') ? '61' : '60'
do_scamp_ast='Y'
do_scamp_phot='Y'
do_overwrite= (n_elements(overwrite) GT 0) ? keyword_set(overwrite) : 0
do_ahead = (n_elements(ahead) GT 0) ? keyword_set(ahead) : 0
do_zp_phot = (n_elements(zp_phot) GT 0) ? keyword_set(zp_phot) : 0
do_debug = (n_elements(debug) GT 0) ? keyword_set(debug) : 0
do_tile_orig=do_tile
do_filter_orig=do_filter
do_program_orig=do_program
do_n_cpu=4
do_fwhm=(n_elements(fwhm) GT 0) ? fwhm : 99.
do_dither= (n_elements(dither) GT 0) ? keyword_set(dither) : 0
do_align_fwhm = (n_elements(align_fwhm) GT 0) ? align_fwhm : ''

input_dir='/Volumes/RAID/Magellanic_stream/DECam'
output_dir='/Volumes/RAID/Magellanic_stream/DECam'


do_program_full=strsplit(do_program,',', /extract)
if do_program EQ 'ALL' then begin
	temp_program=strsplit(file_search(input_dir+'/*', /TEST_DIRECTORY),'/', /extract)
	do_program=strarr(n_elements(temp_program))
	for i=0L, n_elements(temp_program)-1 do do_program[i]=(temp_program[i])[-1]
	gv=where(do_program NE 'stacks' AND do_program NE 'pipeline', n_gv)
	if n_gv GT 0 then do_program=do_program[gv]
endif else $
if n_elements(do_program_full) GT 1 then begin
	do_program=string(do_program_full)
endif

input_im_dir=input_dir+'/'+do_program+'/processed'
input_calib_dir=input_dir+'/'+do_program+'/calib'
input_chip_fwhm=[28,35] ; Where to measure FWHM
input_chip_flux_radius=[27,28,29,34,35,36] ; Where to measure FWHM
input_chip_iq=[27,28,35,36] ; Where to measure FWHM

if stregex(do_fwhm,'-') GT 0 AND stregex(do_fwhm,',') GT 0 then begin
	do_fwhm_type='range'
	input_fwhm_full=strsplit(strsplit(do_fwhm,',', /extract),'-', /extract)
endif else $
if stregex(do_fwhm,'-') GT 0 then begin
	do_fwhm_type='range'
	input_fwhm_full=list()
	input_fwhm_full.add, strsplit(strsplit(do_fwhm,',', /extract),'-', /extract)
endif else $
if stregex(do_fwhm,',') GT 0 then begin
	do_fwhm_type='limit'
	input_fwhm_full=strsplit(do_fwhm,',', /extract)
endif $
else begin
	do_fwhm_type='limit'
	input_fwhm_full=do_fwhm
endelse

if stregex(do_align_fwhm,'-') GT 0 then begin
	input_align_fwhm_full=strsplit(strsplit(do_align_fwhm,',', /extract),'-', /extract)
endif $
else begin
	input_align_fwhm_full=strsplit(do_align_fwhm,',', /extract)
endelse

output_im_dir=output_dir+'/'+do_program+'/pipeline/images'
output_calib_dir=output_dir+'/'+do_program+'/pipeline/calib'
output_sex_dir=output_dir+'/'+do_program+'/pipeline/sextractor'
output_sex_check_dir=output_dir+'/'+do_program+'/pipeline/sextractor/check'
output_scamp_dir=output_dir+'/'+do_program+'/pipeline/scamp'
output_stack_swarp_dir=output_dir+'/stacks';output_dir+'/'+do_program+'/pipeline/swarp'
output_stack_sex_dir=output_dir+'/stacks';output_dir+'/'+do_program+'/pipeline/swarp'
output_stack_check_dir=output_dir+'/stacks/check';output_dir+'/'+do_program+'/pipeline/swarp'
output_stack_psfex_dir=output_dir+'/stacks' ;output_dir+'/'+do_program+'/pipeline/psfex'

survey_info = list( {tile:'1', coo:['23:26:07.46','-53:28:13.7'], filter:['g','i','z']}, $
	{tile:'2', coo:['23:53:46.15','-52:32:21.7'], filter:['g','i','z']}, $
	{tile:'3', coo:['00:22:34.67','-51:17:37.6'], filter:['g','i','z']}, $
	{tile:'4', coo:['00:53:43.00','-49:36:27.2'], filter:['g','i','z']}, $
	{tile:'5', coo:['01:03:22.85','-48:35:37.4'], filter:['g','i','z']}, $
	{tile:'6', coo:['01:12:56.00','-47:35:04.6'], filter:['g','i','z']}, $
	{tile:'7', coo:['01:37:49.01','-45:26:39.8'], filter:['g','i','z']}, $
	{tile:'8', coo:['01:59:33.94','-42:55:45.3'], filter:['g','i']}, $
	{tile:'9', coo:['02:18:48.79','-40:11:55.5'], filter:['g','i']} )

survey_sequence = { seq1:['1','1'], seq2:['2','2'], seq3:['3','3'], seq4:['4','4'], seq5:['5','5'], seq6:['6','6'], seq7:['7','7'], seq8:['8','8'], seq9:['9','9']  }

for i=0L, n_elements(do_program)-1 do begin

	if file_test('survey_target_'+do_program[i]+'.dat') AND recipe NE 'database' then begin
		readcol, 'survey_target_'+do_program[i]+'.dat', temp_im_orig_file, temp_filter, temp_tile, temp_dither, temp_weight_orig_file, temp_zp, temp_fwhm, temp_mjd, temp_exptime, temp_n_chip, FORMAT='A,A,A,A,A,F,F,D,F,I', COMMENT='#'
		n_survey_info=n_elements(temp_im_orig_file)
		create_struct, temp_input_target, '', ['im_orig_file','program','filter','tile','dither','weight_orig_file','mjd','mjd_floor','zp','im_file','weight_file','fwhm','sky','exptime','n_chip','sex_cat_file','sex_xml_file','sex_check_file','scamp_cat_file','scamp_head_file','scamp_ahead_file','swarp_im_file','swarp_head_file','sex_zp_file','stack_flag'], 'A,A,A,A,A,A,D,D,F,A,A,F,F,F,I,A,A,A,A,A,A,A,A,A,A', dim=n_survey_info
		temp_input_target.im_orig_file=input_im_dir[i]+'/'+repstr(temp_im_orig_file,'.fits.fz','.fits')
		temp_input_target.weight_orig_file=input_im_dir[i]+'/'+repstr(temp_weight_orig_file,'.fits.fz','.fits')
		temp_input_target.program=do_program[i]
		temp_input_target.filter=temp_filter
		temp_input_target.tile=strtrim(temp_tile,2)
		temp_input_target.dither=strtrim(temp_dither,2)
		temp_input_target.mjd=temp_mjd
		temp_input_target.mjd_floor=floor(temp_mjd+0.2)
		temp_input_target.zp=temp_zp
		temp_input_target.fwhm=99.
		temp_input_target.sky=65000.
		temp_input_target.stack_flag='F'

		temp_input_target.exptime=temp_exptime
		temp_input_target.n_chip=temp_n_chip
		temp_input_target.im_file=output_im_dir[i]+'/'+'ms_t'+strtrim(temp_tile,2)+'_d'+strtrim(temp_dither,2)+'_'+temp_filter+'.fits'
		temp_input_target.weight_file=output_im_dir[i]+'/'+'ms_t'+strtrim(temp_tile,2)+'_d'+strtrim(temp_dither,2)+'_'+temp_filter+'.WEIGHT.fits'
		temp_input_target.sex_cat_file=output_sex_dir[i]+'/'+'ms_t'+strtrim(temp_tile,2)+'_d'+strtrim(temp_dither,2)+'_'+temp_filter+'.ldac'
		temp_input_target.sex_xml_file=output_sex_dir[i]+'/'+'ms_t'+strtrim(temp_tile,2)+'_d'+strtrim(temp_dither,2)+'_'+temp_filter+'.xml'
		temp_input_target.sex_check_file=output_sex_check_dir[i]+'/'+'ms_t'+strtrim(temp_tile,2)+'_d'+strtrim(temp_dither,2)+'_'+temp_filter+'.BACKGROUND.fits'
		temp_input_target.sex_zp_file=output_sex_dir[i]+'/'+'ms_t'+strtrim(temp_tile,2)+'_d'+strtrim(temp_dither,2)+'_'+temp_filter+'_zp.dat'
		temp_input_target.scamp_cat_file=output_scamp_dir[i]+'/'+'ms_t'+strtrim(temp_tile,2)+'_d'+strtrim(temp_dither,2)+'_'+temp_filter+'_stars.ldac'
		temp_input_target.scamp_head_file=output_scamp_dir[i]+'/'+'ms_t'+strtrim(temp_tile,2)+'_d'+strtrim(temp_dither,2)+'_'+temp_filter+'_stars.head'
		temp_input_target.scamp_ahead_file=output_scamp_dir[i]+'/'+'ms_t'+strtrim(temp_tile,2)+'_d'+strtrim(temp_dither,2)+'_'+temp_filter+'_stars.ahead'
		temp_input_target.swarp_im_file=output_im_dir[i]+'/'+'ms_t'+strtrim(temp_tile,2)+'_d'+strtrim(temp_dither,2)+'_'+temp_filter+'.fits'
		temp_input_target.swarp_head_file=output_im_dir[i]+'/'+'ms_t'+strtrim(temp_tile,2)+'_d'+strtrim(temp_dither,2)+'_'+temp_filter+'.head'

		if i EQ 0 then input_target=temp_input_target $
		else input_target=[input_target,temp_input_target]

	endif
endfor
	
if n_elements(input_target) GT 0 then begin
	
	do_tile_full=strsplit(do_tile,',', /extract)
	do_filter_full=strsplit(do_filter,',', /extract)
	do_dither_exclude_full=strsplit(do_dither_exclude,';', /extract)

	if n_elements(do_tile_full) GT 1 then begin
		if n_elements(do_filter_full) GT 1 then begin
			temp_tile1=byte(input_target.tile)
			temp_tile2=byte(do_tile_full)
			temp_filter1=byte(input_target.filter)
			temp_filter2=byte(do_filter_full)

			gv=where( total( strcmp( string(rebin(temp_tile1,[(size(temp_tile1,/dim))[0],n_elements(input_target.tile),n_elements(do_tile_full)])), string(rebin(reform(temp_tile2,[(size(temp_tile2,/dim))[0],1,(size(temp_tile2,/dim))[1]]),[(size(temp_tile2,/dim))[0],n_elements(input_target.tile),n_elements(do_tile_full)])) ), 2) AND total( strcmp( string(rebin(temp_filter1,[(size(temp_filter1,/dim))[0],n_elements(input_target.filter),n_elements(do_filter_full)])), string(rebin(reform(temp_filter2,[(size(temp_filter2,/dim))[0],1,(size(temp_filter2,/dim))[1]]),[(size(temp_filter2,/dim))[0],n_elements(input_target.filter),n_elements(do_filter_full)])) ), 2) GT 0, n_gv)
			
		endif $
		else begin
 			temp_tile1=byte(input_target.tile)
      temp_tile2=byte(do_tile_full)

			if do_filter EQ 'all' then begin
				gv=where( total( strcmp( string(rebin(temp_tile1,[(size(temp_tile1,/dim))[0],n_elements(input_target.tile),n_elements(do_tile_full)])), string(rebin(reform(temp_tile2,[(size(temp_tile2,/dim))[0],1,(size(temp_tile2,/dim))[1]]),[(size(temp_tile2,/dim))[0],n_elements(input_target.tile),n_elements(do_tile_full)])) ), 2) GT 0, n_gv)
			endif $
			else begin
				gv=where( (total( strcmp( string(rebin(temp_tile1,[(size(temp_tile1,/dim))[0],n_elements(input_target.tile),n_elements(do_tile_full)])), string(rebin(reform(temp_tile2,[(size(temp_tile2,/dim))[0],1,(size(temp_tile2,/dim))[1]]),[(size(temp_tile2,/dim))[0],n_elements(input_target.tile),n_elements(do_tile_full)])) ), 2) GT 0) AND input_target.filter EQ do_filter, n_gv)
			endelse
	
		endelse

	endif $
	else begin

		if n_elements(do_filter_full) GT 1 then begin
			temp_filter1=byte(input_target.filter)
			temp_filter2=byte(do_filter_full)

			if do_tile EQ 'all' then begin
				gv=where( total( strcmp( string(rebin(temp_filter1,[(size(temp_filter1,/dim))[0],n_elements(input_target.filter),n_elements(do_filter_full)])), string(rebin(reform(temp_filter2,[(size(temp_filter2,/dim))[0],1,(size(temp_filter2,/dim))[1]]),[(size(temp_filter2,/dim))[0],n_elements(input_target.filter),n_elements(do_filter_full)])) ), 2) GT 0, n_gv)
			endif $
			else begin
				gv=where( input_target.tile EQ do_tile AND (total( strcmp( string(rebin(temp_filter1,[(size(temp_filter1,/dim))[0],n_elements(input_target.filter),n_elements(do_filter_full)])), string(rebin(reform(temp_filter2,[(size(temp_filter2,/dim))[0],1,(size(temp_filter2,/dim))[1]]),[(size(temp_filter2,/dim))[0],n_elements(input_target.filter),n_elements(do_filter_full)])) ), 2) GT 0), n_gv)
			endelse
		
		endif $
		else begin

			if do_tile EQ 'all' AND do_filter EQ 'all' then begin
				n_gv=n_elements(input_target)
				gv=indgen(n_elements(input_target))
			endif else $
			if do_tile EQ 'all' then begin
				gv=where(input_target.filter EQ do_filter, n_gv)
			endif else $
			if do_filter EQ 'all' then begin
				gv=where(input_target.tile EQ do_tile, n_gv)
			endif $
			else begin
				gv=where(input_target.tile EQ do_tile AND input_target.filter EQ do_filter, n_gv)
			endelse

		endelse
	endelse

	if do_dither_exclude NE 'NONE' then begin
		gv_exclude=list()
		for i=0L, n_elements(do_dither_exclude_full)-1 do begin
			dither_exclude_info=strsplit(do_dither_exclude_full[i],':', /extract)
			dither_exclude_program=dither_exclude_info[0]
			dither_exclude_tile=dither_exclude_info[1]
			dither_exclude_filter=dither_exclude_info[2]
			dither_exclude_id=strsplit(dither_exclude_info[3],',', /extract)
			for j=0L, n_elements(dither_exclude_id)-1 do begin
				gv_dither=where(input_target.program EQ dither_exclude_program AND input_target.tile EQ dither_exclude_tile AND input_target.filter EQ dither_exclude_filter AND input_target.dither EQ dither_exclude_id[j], n_gv_dither)
				if n_gv_dither EQ 1 then begin
					gv_exclude.add, gv_dither
				endif
			endfor
		endfor
		gv_exclude=gv_exclude.toarray(type='long')
		do_dither_exclude=repstr(do_dither_exclude,',','_')
	endif

	if n_gv GT 0 then begin
		if n_elements(gv_exclude) GT 0 then begin
			print, 'PIPELINE - Excluding the following files becuase of DITHER_EXCLUDE keyword'
			forprint, input_target[gv_exclude].im_file, input_target[gv_exclude].tile, input_target[gv_exclude].filter, FORMAT='A,2X,A,2X,A', textout=2 

			gv_include=where( total( gv#make_array(n_elements(gv_exclude), value=1) EQ make_array(n_elements(gv), value=1)#gv_exclude, 2) EQ 0, n_gv_include)
			if n_gv_include GT 0 then gv=gv[gv_include]
		endif
		input_target=input_target[gv]
	endif $
	else stop
	wait, 1

	gv=where(input_target.exptime GT 100, n_gv)
	if n_gv GT 0 then input_target=input_target[gv] $
	else stop

	print, 'PIPELINE - The following files will be used'
	forprint, input_target.im_file, input_target.tile, input_target.filter, input_target.exptime, FORMAT='A,2X,A,2X,A,2X,F0.1', textout=2 
	wait, 1

endif

for i=0L, n_elements(do_program)-1 do begin

	if file_test('survey_iq_'+do_program[i]+'.dat') AND recipe NE 'iq' then begin

		readcol, 'survey_iq_'+do_program[i]+'.dat', temp_im_file, temp_filter, temp_tile, temp_dither, temp_sky, temp_fwhm, FORMAT='A,A,A,A,F,F', COMMENT='#'
		n_survey=n_elements(temp_im_file)
		
		if n_survey GT 0 then begin
			create_struct, temp_input_iq, '', ['im_file','filter','tile','dither','sky','fwhm','stack_flag'], 'A,A,A,A,F,F,A', dim=n_survey
			temp_input_iq.im_file=temp_im_file
			temp_input_iq.filter=temp_filter
			temp_input_iq.tile=temp_tile
			temp_input_iq.dither=temp_dither
			temp_input_iq.sky=temp_sky
			temp_input_iq.fwhm=temp_fwhm

			if n_elements(input_iq) EQ 0 then input_iq=temp_input_iq $
			else input_iq=[input_iq,temp_input_iq]
		endif
	endif
endfor

if file_test('survey_iq_manual.dat') AND n_elements(input_iq) GT 0 then begin
	print, 'Reading FWHM from IQ manual file'
	readcol, 'survey_iq_manual.dat', temp_im_file, temp_filter, temp_tile, temp_dither, temp_sky, temp_exptime, temp_fwhm, FORMAT='A,A,A,A,F,F,F', COMMENT='#'

	for j=0L, n_elements(temp_im_file)-1 do begin
		gv=where(input_iq.im_file EQ temp_im_file[j], n_gv)
		if n_gv EQ 1 then begin
			input_iq[gv].fwhm=temp_fwhm[j]
;			input_iq[gv].stack_flag=temp_stack_flag[j]
		endif
	endfor

endif

if n_elements(input_iq) GT 0 then begin

	for j=0L, n_elements(input_iq.im_file)-1 do begin
		gv=where(input_target.im_file EQ input_iq[j].im_file, n_gv)
		if n_gv EQ 1 then begin
			input_target[gv].fwhm=input_iq[j].fwhm
			input_target[gv].sky=input_iq[j].sky
;			input_target[gv].stack_flag=input_iq[j].stack_flag
		endif
	endfor

;	gv_iq=where(input_target.sky LT 22000. AND input_target.fwhm GT 0. AND input_target.stack_flag EQ 'T', n_gv_iq, COMPLEMENT=bv_iq)
	gv_iq=where(input_target.sky LT 22000. AND input_target.fwhm GT 0., n_gv_iq, COMPLEMENT=bv_iq)

	print, 'IQ Removing the following files'
	print, 'Filename     Sky_level (ADU)     FWHM (arsec)'
	forprint, input_target[bv_iq].im_file, input_target[bv_iq].sky, input_target[bv_iq].fwhm*0.2637, FORMAT='A,4X,F0.1,4X,F0.1', text=2
	input_target=input_target[gv_iq]

;	gv_iq=where(input_iq.sky LT 25000. AND input_iq.fwhm GT 0., n_gv_iq, COMPLEMENT=bv_iq)
;	print, 'Removing the following files'
;	forprint, input_iq[bv_iq].im_file, input_iq[bv_iq].sky, FORMAT='A,4X,F0.1', text=2

;	temp_file1=byte(input_target.im_file)
;	temp_file2=byte(input_iq[gv_iq].im_file)
	
;	gv=where( total( strcmp( string(rebin(temp_file1,[(size(temp_file1,/dim))[0],n_elements(input_target.im_file),n_gv_iq])), string(rebin(reform(temp_file2,[(size(temp_file2,/dim))[0],1,(size(temp_file2,/dim))[1]]),[(size(temp_file2,/dim))[0],n_elements(input_target.im_file),n_gv_iq])) ), 2) GT 0, n_gv)
;	input_target=input_target[gv]

endif

for i=0L, n_elements(do_program)-1 do begin

	if file_test('survey_calib_'+do_program[i]+'.dat') AND recipe NE 'database' then begin
		readcol, 'survey_calib_'+do_program[i]+'.dat', temp_im_orig_file, temp_object, temp_filter, temp_weight_orig_file, temp_zp, temp_fwhm, temp_mjd, temp_exptime, temp_n_chip, temp_photflag, FORMAT='A,A,A,A,F,F,D,F,I,A', COMMENT='#'
		n_survey=n_elements(temp_im_orig_file)
		create_struct, temp_input_calib, '', ['im_orig_file','object','program','filter','tile','weight_orig_file','skip','mjd','mjd_floor','zp','airmass','im_file','weight_file','fwhm','exptime','n_chip','sex_dir','sex_cat_file','sex_xml_file','sex_check_file','sex_zp_file','sex_zp_full_file','photflag','scamp_cat_file','scamp_head_file','scamp_ahead_file','swarp_im_file','swarp_head_file'], 'A,A,A,A,A,A,A,D,D,F,F,A,A,F,F,I,A,A,A,A,A,A,A,A,A,A,A,A', dim=n_survey
		temp_input_calib.im_orig_file=input_calib_dir[i]+'/'+repstr(temp_im_orig_file,'.fits.fz','.fits')
		temp_input_calib.weight_orig_file=input_calib_dir[i]+'/'+repstr(temp_weight_orig_file,'.fits.fz','.fits')
		temp_input_calib.object=temp_object
		temp_input_calib.program=do_program[i]
		temp_input_calib.filter=temp_filter
		temp_input_calib.skip='F'
		temp_input_calib.mjd=temp_mjd
		temp_input_calib.mjd_floor=floor(temp_mjd+0.2)
		temp_input_calib.zp=temp_zp
		temp_input_calib.airmass=0.
		temp_input_calib.fwhm=temp_fwhm
		temp_input_calib.exptime=temp_exptime
		temp_input_calib.n_chip=temp_n_chip
		temp_input_calib.photflag=temp_photflag
		temp_mjd=string(temp_mjd, FORMAT='(F0.8)')
		temp_input_calib.im_file=output_im_dir[i]+'/'+'standard_'+temp_object+'_'+temp_filter+'_MJD'+temp_mjd+'.fits'
		temp_input_calib.weight_file=output_im_dir[i]+'/'+'standard_'+temp_object+'_'+temp_filter+'_MJD'+temp_mjd+'.WEIGHT.fits'
		temp_input_calib.sex_dir=output_sex_dir[i]
		temp_input_calib.sex_cat_file=output_sex_dir[i]+'/'+'standard_'+temp_object+'_'+temp_filter+'_MJD'+temp_mjd+'.ldac'
		temp_input_calib.sex_xml_file=output_sex_dir[i]+'/'+'standard_'+temp_object+'_'+temp_filter+'_MJD'+temp_mjd+'.xml'
		temp_input_calib.sex_check_file=output_sex_check_dir[i]+'/'+'standard_'+temp_object+'_'+temp_filter+'_MJD'+temp_mjd+'.BACKGROUND.fits'
		temp_input_calib.sex_zp_file=output_calib_dir[i]+'/'+'standard_'+temp_object+'_'+temp_filter+'_MJD'+temp_mjd+'_zp.dat'
		temp_input_calib.sex_zp_full_file=output_calib_dir[i]+'/'+'standard_'+temp_object+'_'+temp_filter+'_MJD'+temp_mjd+'_zp_full.dat'
		temp_input_calib.scamp_cat_file=output_scamp_dir[i]+'/'+'standard_'+temp_object+'_'+temp_filter+'_MJD'+temp_mjd+'_stars.ldac'
		temp_input_calib.scamp_head_file=output_scamp_dir[i]+'/'+'standard_'+temp_object+'_'+temp_filter+'_MJD'+temp_mjd+'_stars.head'
		temp_input_calib.scamp_ahead_file=output_scamp_dir[i]+'/'+'standard_'+temp_object+'_'+temp_filter+'_MJD'+temp_mjd+'_stars.ahead'
		temp_input_calib.swarp_im_file=output_im_dir[i]+'/'+'standard_'+temp_object+'_'+temp_filter+'_MJD'+temp_mjd+'.fits'
		temp_input_calib.swarp_head_file=output_im_dir[i]+'/'+'standard_'+temp_object+'_'+temp_filter+'_MJD'+temp_mjd+'.head'

		if i EQ 0 then input_calib=temp_input_calib $
		else input_calib=[input_calib,temp_input_calib]

	endif
endfor

if n_elements(input_calib) GT 0 then begin

	do_filter_split=strsplit(do_filter,',', /extract)
	if n_elements(do_filter_split) GT 1 then begin
		gv=where( total( rebin(byte(input_calib.filter),[n_elements(do_filter_split),n_elements(input_calib.filter)]) EQ rebin(transpose(byte(do_filter_split)),[n_elements(do_filter_split),n_elements(input_calib.filter)]),1) GT 0, n_gv)
		if n_gv GT 0 then input_calib=input_calib[gv] $
		else stop
	endif $
	else begin
		case do_filter of
			'all': begin
				n_gv=n_elements(input_calib)
				gv=indgen(n_elements(input_calib))
			end
			else: gv=where(input_calib.filter EQ do_filter, n_gv)
		endcase
		if n_gv GT 0 then input_calib=input_calib[gv] $
		else stop
	endelse

endif

for i=0L, n_elements(do_program)-1 do begin

	if file_test('survey_zp_'+do_program[i]+'.dat') then begin
		readcol, 'survey_zp_'+do_program[i]+'.dat', temp_mjd, temp_filter, temp_zp, temp_k, temp_zp_error, temp_k_error, temp_photflag, FORMAT='F,A,F,F,F,F,A', COMMENT='#'
		n_survey=n_elements(temp_mjd)
		create_struct, temp_input_zp, '', ['mjd','filter','zp','k','zp_error','k_error','photflag'], 'F,A,F,F,F,F,A', dim=n_survey
		temp_input_zp.mjd=temp_mjd
		temp_input_zp.filter=temp_filter
		temp_input_zp.zp=temp_zp
		temp_input_zp.k=temp_k
		temp_input_zp.zp_error=temp_zp_error
		temp_input_zp.k_error=temp_k_error
		temp_input_zp.photflag=temp_photflag

		if i EQ 0 then input_zp=temp_input_zp $
		else input_zp=[input_zp,temp_input_zp]

	endif
endfor

do_program_plain=repstr(do_program_orig,',','_')
do_tile_plain=repstr(do_tile_orig,',','_')
do_filter_plain=repstr(do_filter_orig,',','_')
loadct, 12

if recipe EQ 'database' then begin

	print, 'DATABASE - The following programs were found ', do_program

	for ii=0L, n_elements(do_program)-1 do begin

		; Processing Science images
		im_file=file_search(input_im_dir[ii], '*.fz', count=n_im_file)
	
		create_struct, temp_target, '', ['dir','im_file','weight_file','expnum','type','date','ra','dec','mjd','mjd_floor','exptime','tile','dither','filter','zp','fwhm','n_chip'], 'A,A,A,I,A,A,D,D,D,D,F,A,A,A,F,F,I', dim=1
	
		for i=0L, n_im_file-1 do begin
			print, 'DATABASE - Processing file ', (strsplit(im_file[i],'/',/extract))[-1]
	
			im_h=headfits(im_file[i])
			fits_info, im_file[i], n_ext=n_ext, /silent
			temp_target.dir=input_im_dir[ii]
			temp_target.im_file=(strsplit(im_file[i],'/',/extract))[-1]
			temp_target.expnum=fxpar(im_h, 'EXPNUM')
			temp_target.type=strtrim(fxpar(im_h, 'PRODTYPE'),2)
			temp_target.date=fxpar(im_h, 'DATE-OBS')
			temp_target.mjd=fxpar(im_h, 'MJD-OBS')
			temp_target.mjd_floor=floor(fxpar(im_h, 'MJD-OBS')+0.2)
			temp_target.ra=ten(fxpar(im_h, 'RA'))*360./24.
			temp_target.dec=ten(fxpar(im_h, 'DEC'))
			temp_target.exptime=float(fxpar(im_h, 'EXPTIME'))
			temp_target.filter=strmid(fxpar(im_h, 'FILTER'),0,1)
			temp=size(fxpar(im_h, 'MAGZERO'),/type)
			temp_target.zp= (temp EQ 4 OR temp EQ 5) ? fxpar(im_h, 'MAGZERO') : 31.5
			temp_target.n_chip=n_ext
	
			if i EQ 0 then input_target=temp_target else input_target=[input_target,temp_target]
		endfor
		gv_sort=sort(input_target.date+' '+input_target.type)
	
		if file_test(input_im_dir[ii]+'/pipeline_database.dat') EQ 0 then begin
			openw, lun, input_im_dir[ii]+'/pipeline_database.dat', /get_lun
			printf, lun, '# FILE	 DATE-OBS	 EXPNUM	 OBJECT	 FILTER	 EXPTIME	 SEQID	 OBSTYPE	 PROCTYPE		PRODTYPE	 MAGZERO'
			for i=0L, n_im_file-1 do begin
				command='dfits '+im_file[gv_sort[i]]+' | fitsort -d DATE-OBS EXPNUM OBJECT FILTER EXPTIME SEQID OBSTYPE PROCTYPE PRODTYPE MAGZERO'
				spawn, command, result
				result=repstr(result, input_im_dir[ii]+'/', '')
				printf, lun, result
			endfor
			free_lun, lun
		endif
	
		gv_im=where(input_target.type EQ 'image' AND input_target.exptime GT 30., n_gv_im)
		gv_weight=where(input_target.type EQ 'wtmap' AND input_target.exptime GT 30., n_gv_weight)
		gv_sort=sort(input_target[gv_im].date+' '+input_target[gv_im].type)
		mjd_uniq=input_target[gv_im[uniq(input_target[gv_im].mjd_floor, sort(input_target[gv_im].mjd_floor))]].mjd_floor
		filter_uniq=input_target[gv_im[uniq(input_target[gv_im].filter, sort(input_target[gv_im].filter))]].filter
	
		for i=0L, n_gv_im-1 do begin
			gv=where(input_target[gv_weight].expnum EQ input_target[gv_im[i]].expnum, n_gv)
			if n_gv EQ 1 then input_target[gv_im[i]].weight_file=input_target[gv_weight[gv]].im_file $
			else stop
		endfor
	
		input_target=input_target[gv_im[gv_sort]]
	
		for i=0L, n_elements(survey_info)-1 do begin
			for j=0L, n_elements((survey_info[i]).filter)-1 do begin
				gv=where( sqrt( (input_target.ra-ten((survey_info[i]).coo[0])*360./24)^2 + (input_target.dec-ten((survey_info[i]).coo[1]))^2 ) LE 0.1 AND input_target.filter EQ (survey_info[i]).filter[j], n_gv)
				input_target[gv].tile=(survey_info[i]).tile
				print, (survey_info[i]).filter[j], (survey_info[i]).tile
	;			input_target[gv_im[gv_sort[gv]]].dither=string(indgen(n_gv)+1,FORMAT='(I0)')
			endfor
		endfor
		gv_im=where(input_target.type EQ 'image' AND input_target.exptime GT 30., n_gv_im)
		tile_uniq=input_target[gv_im[uniq(input_target[gv_im].tile, sort(input_target[gv_im].tile))]].tile
	
		for i=0L, n_elements(filter_uniq)-1 do begin
			i_dither=1L
			for j=0L, n_elements(tile_uniq)-1 do begin
				gv=where( input_target.filter EQ filter_uniq[i] AND input_target.tile EQ tile_uniq[j], n_gv)
				if n_gv GT 0 then begin
					input_target[gv].dither=string(indgen(n_gv)+1, FORMAT='(I0)')
				endif
			endfor
		endfor
		forprint, input_target.im_file, input_target.mjd_floor, input_target.filter, input_target.tile, input_target.dither, format='A,4X,F0,4X,A,4X,A,4X,A', textout=2
	
		gv=where(input_target.tile NE ' ' AND input_target.dither NE ' ', n_gv)
		if n_gv GT 0 then begin
			openw, lun, 'survey_target_'+do_program[ii]+'.dat', /get_lun
			printf, lun, '#  im_file        filter  tile  dither  weight_file     zp      FWHM    MJD   EXPTIME   NCHIP'
			for i=0L, n_gv-1 do begin
				printf, lun, input_target[gv[i]].im_file, input_target[gv[i]].filter, input_target[gv[i]].tile, input_target[gv[i]].dither, input_target[gv[i]].weight_file, input_target[gv[i]].zp, input_target[gv[i]].fwhm, input_target[gv[i]].mjd, input_target[gv[i]].exptime, input_target[gv[i]].n_chip, FORMAT='(A,4X,A,4X,A,4X,A,4X,A,4X,F5.2,4X,F4.1,4X,F0.6,4X,F0.1,4X,I0)'
			endfor
			free_lun, lun
		endif
	
		; Processing Calibration images
		im_file=file_search(input_calib_dir[ii], '*.fz', count=n_im_file)

		create_struct, temp_target, '', ['dir','im_file','weight_file','expnum','type','date','ra','dec','mjd','mjd_floor','exptime','tile','dither','filter','zp','fwhm','n_chip'], 'A,A,A,I,A,A,D,D,D,D,F,A,A,A,F,F,I', dim=1
	
		for i=0L, n_im_file-1 do begin
			print, 'DATABASE - Processing file ', (strsplit(im_file[i],'/',/extract))[-1]
	
			im_h=headfits(im_file[i])
			fits_info, im_file[i], n_ext=n_ext, /silent
			temp_target.dir=input_im_dir[ii]
			temp_target.im_file=(strsplit(im_file[i],'/',/extract))[-1]
			temp_target.expnum=fxpar(im_h, 'EXPNUM')
			temp_target.type=strtrim(fxpar(im_h, 'PRODTYPE'),2)
			temp_target.date=fxpar(im_h, 'DATE-OBS')
			temp_target.mjd=fxpar(im_h, 'MJD-OBS')
			temp_target.mjd_floor=floor(fxpar(im_h, 'MJD-OBS')+0.2)
			temp_target.ra=ten(fxpar(im_h, 'RA'))*360./24.
			temp_target.dec=ten(fxpar(im_h, 'DEC'))
			temp_target.exptime=float(fxpar(im_h, 'EXPTIME'))
			temp_target.filter=strmid(fxpar(im_h, 'FILTER'),0,1)
			temp=size(fxpar(im_h, 'MAGZERO'),/type)
			temp_target.zp= (temp EQ 4 OR temp EQ 5) ? fxpar(im_h, 'MAGZERO') : 31.5
			temp_target.n_chip=n_ext
	
			if i EQ 0 then input_target=temp_target else input_target=[input_target,temp_target]
		endfor
		gv_sort=sort(input_target.date+' '+input_target.type)



		im_date=list()
		im_type=list()

		if n_im_file GT 0 then begin	
			for i=0L, n_im_file-1 do begin
				im_h=headfits(im_file[i])
				im_date.add, fxpar(im_h, 'DATE-OBS')
				im_type.add, fxpar(im_h, 'PRODTYPE')
			endfor
			im_date=im_date.toarray(type='string')
			im_type=im_type.toarray(type='string')
			gv=sort(im_date+' '+im_type)
	
			openw, lun, input_calib_dir[ii]+'/pipeline_database.dat', /get_lun
			printf, lun, '# FILE	 DATE-OBS   MJD-OBS	 EXPNUM	 OBJECT	 FILTER	 EXPTIME	 SEQID	 OBSTYPE	 PROCTYPE		PRODTYPE	 MAGZERO'
			for i=0L, n_im_file-1 do begin
				command='dfits '+im_file[gv[i]]+' | fitsort -d DATE-OBS MJD-OBS EXPNUM OBJECT FILTER EXPTIME SEQID OBSTYPE PROCTYPE PRODTYPE MAGZERO'
				spawn, command, result
				result=repstr(result, input_calib_dir[ii]+'/', '')
				printf, lun, result
			endfor
			free_lun, lun
		endif

	endfor

endif else $
if recipe EQ 'setup' then begin

	for i=0L, n_elements(do_program)-1 do begin
		if file_test(output_im_dir[i], /directory) EQ 0 then file_mkdir, output_im_dir[i]
		if file_test(output_calib_dir[i], /directory) EQ 0 then file_mkdir, output_calib_dir[i]
		if file_test(output_sex_dir[i], /directory) EQ 0 then file_mkdir, output_sex_dir[i]
		if file_test(output_sex_check_dir[i], /directory) EQ 0 then file_mkdir, output_sex_check_dir[i]
		if file_test(output_scamp_dir[i], /directory) EQ 0 then file_mkdir, output_scamp_dir[i]
	endfor
	if file_test(output_stack_swarp_dir, /directory) EQ 0 then file_mkdir, output_stack_swarp_dir
	if file_test(output_stack_sex_dir, /directory) EQ 0 then file_mkdir, output_stack_sex_dir
	if file_test(output_stack_check_dir, /directory) EQ 0 then file_mkdir, output_stack_check_dir
	if file_test(output_stack_psfex_dir, /directory) EQ 0 then file_mkdir, output_stack_psfex_dir

	if do_overwrite then begin
		file_delete, input_target.im_file, /noexpand, /allow_non, /quiet
		file_delete, input_target.weight_file, /noexpand, /allow_non, /quiet
	endif

	for i=0, n_elements(input_target)-1 do begin

		print, 'SETUP - Processing file ', input_target[i].im_orig_file

		if file_test(input_target[i].im_orig_file+'.fz', /regular) EQ 0 OR file_test(input_target[i].weight_orig_file+'.fz', /regular) EQ 0 then begin
			print, 'Some files are missing'
			stop
		endif

		if file_test(input_target[i].im_orig_file, /regular) EQ 0 then begin
			command='funpack '+input_target[i].im_orig_file+'.fz'
			print, command
			spawn, command
		endif

		if file_test(input_target[i].weight_orig_file, /regular) EQ 0 then begin
			command='funpack '+input_target[i].weight_orig_file+'.fz'
			print, command
			spawn, command
		endif

		if file_test(input_target[i].im_file, /regular) EQ 0 then begin
			command='ln -s '+input_target[i].im_orig_file+' '+input_target[i].im_file
			print, command
			spawn, command
		endif

		if file_test(input_target[i].weight_file, /regular) EQ 0 then begin
			command='ln -s '+input_target[i].weight_orig_file+' '+input_target[i].weight_file
			print, command
			spawn, command
		endif
	endfor

	for i=0, n_elements(input_calib)-1 do begin

		if file_test(input_calib[i].im_orig_file+'.fz', /regular) EQ 0 OR file_test(input_calib[i].weight_orig_file+'.fz', /regular) EQ 0 then begin
			print, 'Some files are missing'
			stop
		endif

		if file_test(input_calib[i].im_orig_file, /regular) EQ 0 then begin
			command='funpack '+input_calib[i].im_orig_file+'.fz'
			print, command
			spawn, command

			fix_header, input_calib[i].im_orig_file
		endif

		if file_test(input_calib[i].weight_orig_file, /regular) EQ 0 then begin
			command='funpack '+input_calib[i].weight_orig_file+'.fz'
			print, command
			spawn, command

			fix_header, input_calib[i].weight_orig_file
		endif

		if file_test(input_calib[i].im_file, /regular) EQ 0 then begin
			command='ln -s '+input_calib[i].im_orig_file+' '+input_calib[i].im_file
			print, command
			spawn, command
		endif

		if file_test(input_calib[i].weight_file, /regular) EQ 0 then begin
			command='ln -s '+input_calib[i].weight_orig_file+' '+input_calib[i].weight_file
			print, command
			spawn, command
		endif
	endfor

endif else $
if recipe EQ 'compute zp' then begin

	vig_diam=101
	sex_mag_bin=0.5
	sex_radius_bin=0.2

	if do_overwrite then file_delete, input_calib.sex_cat_file, /noexpand, /allow_non, /quiet

	for i=0L, n_elements(input_calib)-1 do begin

		im_h=headfits(input_calib[i].im_file)
		im_filter=fxpar(im_h, 'FILTER')
		im_exptime=fxpar(im_h, 'EXPTIME')
		im_zp=fxpar(im_h, 'MAGZERO')
		im_ra=ten(fxpar(im_h, 'RA'))*360./24.
		im_dec=ten(fxpar(im_h, 'DEC'))
		im_mjd=fxpar(im_h, 'MJD-OBS')
		im_airmass=tai2airmass(im_ra,im_dec,2000., mjd=im_mjd, lon=-70.8065, latitude=-30.1692, alt=2243.31)
		input_calib[i].airmass=im_airmass

		print, 'Iteration ', i
		print, input_calib[i].im_file, date_conv(im_mjd+2400000.5D,'S'), im_zp, im_exptime, im_airmass, im_filter, FORMAT='("Image filename ",A,2X,A,2X,F0.2,2X,F0.2,2X,F0.2,2X,A)'

		case input_calib[i].filter of
			'u': begin
				sex_satur_level='30000.'
				sex_mag_range=[7,15.]
				sex_flux_radius_min=2
				plot_mag_range=[18,5]
				plot_flux_radius_range=[2,12]
				end
			'g': begin
				sex_satur_level='30000.'
				sex_mag_range=[11,19.]
				sex_flux_radius_min=2
				plot_mag_range=[20,8]
				plot_flux_radius_range=[2,12]
				end
			'r': begin
				sex_satur_level='30000.'
				sex_mag_range=[12,18.]
				sex_flux_radius_min=2
				plot_mag_range=[20,8]
				plot_flux_radius_range=[2,12]
				end
			'i': begin
				sex_satur_level='30000.'
				sex_mag_range=[12,18.]
				sex_flux_radius_min=2
				plot_mag_range=[20,8]
				plot_flux_radius_range=[2,12]
				end
			'z': begin
				sex_satur_level='30000.'
				sex_mag_range=[12,20.]
				sex_flux_radius_min=2
				plot_mag_range=[20,8]
				plot_flux_radius_range=[2,12]
				end
			'Y': begin
				input_calib[i].skip = 'T'
				end
			else: stop
		endcase

		if input_calib[i].skip EQ 'T' then begin
			print, 'Y-band filter images are not supported'
			continue
		endif

		im_h=headfits(input_calib[i].im_file)
		im_ra=ten(fxpar(im_h, 'RA'))*360./24.
		im_dec=ten(fxpar(im_h, 'DEC'))
		im_dec_sign=im_dec GT 0. ? '+' : '-'
		im_exptime=fxpar(im_h, 'EXPTIME')
		im_zp=fxpar(im_h, 'MAGZERO')
		im_mjd=fxpar(im_h, 'MJD-OBS')
		im_airmass=tai2airmass(im_ra,im_dec,2000., mjd=im_mjd, lon=-70.8065, latitude=-30.1692, alt=2243.31)
		input_calib[i].airmass=im_airmass

		file_sex_cat=file_info(input_calib[i].sex_cat_file, /noexpand)
		file_sex_zp=file_info(input_calib[i].sex_zp_file, /noexpand)
		file_sex_zp_full=file_info(input_calib[i].sex_zp_full_file, /noexpand)

		if file_sex_cat.exists AND file_sex_cat.size GT 0 then fits_info, input_calib[i].sex_cat_file, n_ext=sex_n_ext, /silent $
		else sex_n_ext=0

		if sex_n_ext EQ input_calib[i].n_chip*2 AND file_sex_zp.exists AND file_sex_zp.size GT 100 AND file_sex_zp_full.exists AND file_sex_zp_full.size GT 100 then begin
			print, 'file_sex_zp EXISTS: ', input_calib[i].sex_zp_file
			print, 'file_sex_zp_full EXISTS: ', input_calib[i].sex_zp_full_file
			continue
		endif

		im_fwhm=list()
		openw, lun_zp, input_calib[i].sex_zp_file, /get_lun
		openw, lun_zp_full, input_calib[i].sex_zp_full_file, /get_lun
		printf, lun_zp, '# Ext  DECam_zp   airmass   zp   zp_err  zp_nstars   FWHM'
		printf, lun_zp_full, '# Ext  NUMBER   RA   DEC   FLUX_APER   FLUX_AUTO   mag_ref   mag_diff   airmass'

		if file_sex_cat.exists AND file_sex_cat.size GT 0 then fits_info, input_calib[i].sex_cat_file, n_ext=sex_n_ext, /silent $
		else sex_n_ext=0

		if sex_n_ext NE input_calib[i].n_chip*2 OR do_overwrite EQ 1 then begin

			for j=0L, n_elements(input_chip_fwhm)-1 do begin
				im_data=readfits(input_calib[i].im_file, im_h, ext=input_chip_fwhm[j])
				writefits, input_calib[i].sex_dir+'/im_sextractor.fits', im_data, im_h
				wim_data=readfits(input_calib[i].weight_file, wim_h, ext=input_chip_fwhm[j])
				writefits, input_calib[i].sex_dir+'/wim_sextractor.fits', wim_data, wim_h
		
				im_size=size(im_data, /dim)
				im_gain=fxpar(im_h, 'GAINA')
				im_ron=fxpar(im_h, 'RDNOISEA')
	
				command = 'sex '+input_calib[i].sex_dir+'/im_sextractor.fits' +' -c sex_config/ctio_decam.sex -CATALOG_NAME '+input_calib[i].sex_dir+'/im_sextractor.ldac'+' -WEIGHT_IMAGE '+input_calib[i].sex_dir+'/wim_sextractor.fits'+' -XML_NAME '+input_calib[i].sex_xml_file+' -SATUR_LEVEL '+sex_satur_level+' -MAG_ZEROPOINT '+string(input_calib[i].zp,FORMAT='(F0.2)')+' -CHECKIMAGE_TYPE BACKGROUND -CHECKIMAGE_NAME '+input_calib[i].sex_check_file
				print, command
				spawn, command
		
				cat_sex=mrdfits(input_calib[i].sex_dir+'/im_sextractor.ldac', 2, cat_sex_h, COLUMNS=['NUMBER','X_IMAGE','Y_IMAGE','FLUX_RADIUS','MAG_AUTO','FLUX_AUTO','FLAGS'], /silent)
				plothist, (cat_sex[where(cat_sex.mag_auto GT sex_mag_range[0] and cat_sex.mag_auto LT sex_mag_range[1] AND cat_sex.flux_radius GT sex_flux_radius_min AND cat_sex.flags LT 8, n_gv)]).flux_radius, temp_xhist, temp_yhist, bin=sex_radius_bin, /noplot
				temp=max(temp_yhist, gv) & sex_radius=temp_xhist[gv]
				sex_radius=median((cat_sex[where(cat_sex.mag_auto GT sex_mag_range[0] AND cat_sex.mag_auto LT sex_mag_range[1] AND cat_sex.flux_radius GT sex_radius*0.9 AND cat_sex.flux_radius LT sex_radius*1.1 , n_gv)]).flux_radius)

				gv_plot=where(cat_sex.flags LT 8, n_gv_plot)
				plot, cat_sex[gv_plot].flux_radius, cat_sex[gv_plot].mag_auto, psym=1, xrange=plot_flux_radius_range, yrange=plot_mag_range, /ystyle, /xstyle
				oplot, sex_radius*[1,1], [0,100], color=200
				oplot, sex_radius*[0.9,0.9], [0,100], line=2, color=200
				oplot, sex_radius*[1.1,1.1], [0,100], line=2, color=200
				oplot, [0,100], sex_mag_range[0]*[1,1], line=2, color=100
				oplot, [0,100], sex_mag_range[1]*[1,1], line=2, color=100
				wait, 0.2
	
				gv_stars=where(cat_sex.mag_auto GT sex_mag_range[0] AND cat_sex.mag_auto LT sex_mag_range[1]-1 AND cat_sex.flux_radius GT sex_radius*0.9 AND cat_sex.flux_radius LT sex_radius*1.1 AND cat_sex.flags LE 1 AND cat_sex.x_image GT vig_diam/2. AND cat_sex.x_image LT (im_size[0]-vig_diam/2.) AND cat_sex.y_image GT vig_diam/2. AND cat_sex.y_image LT (im_size[1]-vig_diam/2.), n_gv_stars)
	
		    if n_gv_stars GT 0 then begin
					gv_stars=gv_stars[sort(cat_sex[gv_stars].mag_auto)]
	
					temp_n=n_gv_stars<5
			  	    temp_x=cat_sex[gv_stars].x_image
			   	  	temp_y=cat_sex[gv_stars].y_image
			
			  	    for k=0L, temp_n-1 do begin
			  	    	if floor(temp_x[k]-(vig_diam-1.)/2) GT 0 AND ceil(temp_x[k]+(vig_diam-1.)/2) LT im_size[0] $
			  	    		AND floor(temp_y[k]-(vig_diam-1.)/2) GT 0 AND ceil(temp_y[k]+(vig_diam-1.)/2) LT im_size[1] then begin

					        x_range=[ floor(temp_x[k]-(vig_diam-1.)/2), ceil(temp_x[k]+(vig_diam-1.)/2) ]
				    	    y_range=[ floor(temp_y[k]-(vig_diam-1.)/2), ceil(temp_y[k]+(vig_diam-1.)/2) ]
							vig_data = im_data[x_range[0]:x_range[1],y_range[0]:y_range[1]]; - im_sky
							vig_size=size(vig_data, /dim)

							x_center_range=[ floor((vig_diam-1.)/4), ceil(-1-(vig_diam-1.)/4) ]
							y_center_range=[ floor((vig_diam-1.)/4), ceil(-1-(vig_diam-1.)/4) ]
							vig_center_data=vig_data[x_center_range[0]:x_center_range[1],y_center_range[0]:y_center_range[1]]
							vig_center_size=size(vig_center_data, /dim)

							im_max=max(vig_center_data, gv_max)
							im_c= [gv_max mod vig_center_size[0], gv_max/vig_center_size[0]]
							im_c += [x_center_range[0], y_center_range[0]]

							dist_circle, vig_mask, vig_size, im_c[0], im_c[1]
							vig_sky=median(vig_data[where(vig_mask GT 20., n_sky)])
							vig_mag=25.-2.5*alog10(total(vig_data[where(vig_mask LE 20., n_star)]) - n_star*vig_sky )

							vig_res=abs( vig_data-vig_sky - max(vig_center_data-vig_sky)/2. )
							gv=sort(vig_res*vig_mask^2)
							vig_fwhm=2*median(vig_mask[gv[1:5]])
							print, 'Radius for computing FWHM ', vig_mask[gv[1:5]]
							plot, vig_mask, vig_data-vig_sky, xrange=[0,30], psym=1
							oplot, [0,100], max(vig_data-vig_sky)/2.*[1,1], line=2, color=200
							oplot, vig_fwhm/2.*[1,1], [-1e5,1e5], line=2, color=200

							vig_skyrad=4*vig_fwhm < 40.
							vig_psfrad=3*vig_fwhm < 50.
							vig_fitrad=vig_fwhm < 50.

							gcntrd, vig_data, im_c[0], im_c[1], im_cx, im_cy, vig_fwhm
							dist_circle, vig_mask, vig_size, im_cx, im_cy
							vig_sky=median(vig_data[where(vig_mask GT vig_skyrad, n_sky)])
							vig_mag=25.-2.5*alog10(total(vig_data[where(vig_mask LE vig_skyrad, n_star)]) - n_star*vig_sky )
							plot, vig_mask, vig_data-vig_sky, xrange=[0,30], psym=1
							oplot, [0,100], max(vig_data-vig_sky)/2.*[1,1], line=2, color=200
							oplot, vig_fwhm/2.*[1,1], [-1e5,1e5], line=2, color=200

							psf_param=[]
							getpsf, vig_data, im_cx, im_cy, vig_mag, vig_sky, im_ron, im_gain, psf_param, psf_residuals, [0], vig_psfrad, vig_fitrad, input_calib[i].sex_dir+'/im_sextractor_psf.fits' 

							catch, error
							if error EQ 0 then begin

								im_fwhm.add, 2*sqrt(2*alog(2))*sqrt( (psf_param[3]^2 + psf_param[4]^2)/2. )
								print, 'FWHM ', im_fwhm[-1],' pixels'

								dist_circle, vig_mask, vig_size, im_cx, im_cy
								plot, vig_mask, vig_data-vig_sky, xrange=[0,30], psym=1
								x=findgen(100)/10.
								y=gaussian(x,[psf_param[0], 0., mean(psf_param[3:4])])
								oplot, x, y, line=2, color=200
								oplot, im_fwhm[-1]/2.*[1,1], [-1e5,1e5], line=2, color=200
							endif else begin
								catch, /cancel
								print, 'Erro code ', error
								print, 'Procedure GETPSF failed'
								stop
							endelse

						endif
					endfor
				endif	
	
			endfor
	
			im_fwhm=median(im_fwhm.toarray(type='float'))
			print, 'Median FWHM ', im_fwhm,' pixels'

;		endif


;		im_h=headfits(input_calib[i].im_file)
;		im_ra=ten(fxpar(im_h, 'RA'))*360./24.
;		im_dec=ten(fxpar(im_h, 'DEC'))
;		im_dec_sign=im_dec GT 0. ? '+' : '-'
;		im_exptime=fxpar(im_h, 'EXPTIME')
;		im_zp=fxpar(im_h, 'MAGZERO')
;		im_mjd=fxpar(im_h, 'MJD-OBS')
;		im_airmass=tai2airmass(im_ra,im_dec,2000., mjd=im_mjd, lon=-70.8065, latitude=-30.1692, alt=2243.31)
;		input_calib[i].airmass=im_airmass


;		if file_test(input_calib[i].sex_cat_file, /regular) then fits_info, input_calib[i].sex_cat_file, n_ext=sex_n_ext, /silent $
;		else sex_n_ext=0
;
;		if sex_n_ext NE input_n_chip*2 OR do_overwrite EQ 1 then begin

			command = 'sex '+input_calib[i].im_file +' -c sex_config/ctio_decam_stack.sex -CATALOG_NAME '+input_calib[i].sex_cat_file+' -WEIGHT_IMAGE '+input_calib[i].weight_file+' -XML_NAME '+input_calib[i].sex_xml_file+' -SATUR_LEVEL '+sex_satur_level+' -MAG_ZEROPOINT '+string(input_calib[i].zp,FORMAT='(F0.2)')+' -CHECKIMAGE_TYPE BACKGROUND -CHECKIMAGE_NAME '+input_calib[i].sex_check_file+' -PHOT_APERTURES '+strjoin(string(2*[1,2,3,4,5]*im_fwhm, FORMAT='(F0.2)'),',')
			print, command
			spawn, command

		endif

;		pipe_airmass=fltarr(input_n_chip)
;		for j=0L, input_n_chip-1 do begin
;			im_h=headfits(input_calib[i].im_file, ext=j+1)
;			extast, im_h, im_ast
;			xy2ad, im_ast.naxis[0]/2., im_ast.naxis[1]/2., im_ast, pipe_ra, pipe_dec
;			pipe_airmass[j]=tai2airmass(pipe_ra,pipe_dec,2000., mjd=im_mjd, lon=-70.8065, latitude=-30.1692, alt=2243.31)
;		endfor

		im_h=headfits(input_calib[i].im_file)
		im_ra=ten(fxpar(im_h, 'RA'))*360./24.
		im_dec=ten(fxpar(im_h, 'DEC'))
		im_dec_sign=im_dec GT 0. ? '+' : '-'

		if do_standard EQ 'southern' then begin

			cat_file=file_search('standard/southern_stars', '*.dat.trimmed')
			for j=0L, n_elements(cat_file)-1 do begin
				readcol, cat_file[j], id, ra, dec, mag_u, mag_u_err, mag_u_n, mag_g, mag_g_err, mag_g_n, mag_r, mag_r_err, mag_r_n, mag_i, mag_i_err, mag_i_n, mag_z, mag_z_err, mag_z_n, FORMAT='A,A,A,F,F,I,F,F,I,F,F,I,F,F,I,F,F,I,X,X', comment='#'
				create_struct, cat_ref_partial, '', ['ra','dec','mag','magerr'], 'F,F,5F,5F', dim=n_elements(id)

				cat_ref_partial.ra=tenv(ra)*360./24
				cat_ref_partial.dec=tenv(dec)
				cat_ref_partial.mag[0]=mag_u
				cat_ref_partial.mag[1]=mag_g
				cat_ref_partial.mag[2]=mag_r
				cat_ref_partial.mag[3]=mag_i
				cat_ref_partial.mag[4]=mag_z
				cat_ref_partial.magerr[0]=mag_u_err
				cat_ref_partial.magerr[1]=mag_g_err
				cat_ref_partial.magerr[2]=mag_r_err
				cat_ref_partial.magerr[3]=mag_i_err
				cat_ref_partial.magerr[4]=mag_z_err
			
				if j EQ 0 then cat_ref=cat_ref_partial $
				else cat_ref=[cat_ref,cat_ref_partial]
				
			endfor
			n_gv_ref=n_elements(cat_ref)
	
		endif else $
		if do_standard EQ 'SDSS' then begin
			command='aclient_cgi cocat1.u-strasbg.fr sdss9 -c '+string(im_ra, im_dec_sign, abs(im_dec), FORMAT='(F0.6,A,F0.6)')+' -r 80. -lmg 10.,20. -m 10000000'
			spawn, command, cat_ref_data

			cat_ref_data=repstr(cat_ref_data,'---','NaN')
			gv_ref=where( strmid(cat_ref_data,0,1) NE '#', n_gv_ref)
			cat_ref_data=cat_ref_data[gv_ref]
			create_struct, cat_ref, '', ['ra','dec','mag','magerr'], 'F,F,5F,5F', dim=n_gv_ref
			cat_ref.ra=float(strmid(cat_ref_data, 26, 10))
			cat_ref.dec=float(strmid(cat_ref_data, 36, 10))
			cat_ref.mag[0]=float(strmid(cat_ref_data, 71, 6))
			cat_ref.mag[1]=float(strmid(cat_ref_data, 84, 6))
			cat_ref.mag[2]=float(strmid(cat_ref_data, 97, 6))
			cat_ref.mag[3]=float(strmid(cat_ref_data, 110, 6))
			cat_ref.mag[4]=float(strmid(cat_ref_data, 123, 6))
			cat_ref.magerr[0]=float(strmid(cat_ref_data, 78, 5))
			cat_ref.magerr[1]=float(strmid(cat_ref_data, 91, 5))
			cat_ref.magerr[2]=float(strmid(cat_ref_data, 104, 5))
			cat_ref.magerr[3]=float(strmid(cat_ref_data, 117, 5))
			cat_ref.magerr[4]=float(strmid(cat_ref_data, 130, 5))
		endif else begin
			print, 'Check if standard stars images belong to the SDSS or Southern star catalog'
			stop
		endelse

		flux_radius_median=list()
		foreach j, input_chip_flux_radius do begin
			cat_sex=mrdfits(input_calib[i].sex_cat_file, 2*(j+1), cat_sex_h, COLUMNS=['NUMBER','ALPHA_J2000','DELTA_J2000','FLUX_RADIUS','MAG_AUTO','FLUX_AUTO','MAG_APER','FLUX_APER','FLAGS'], /silent)
			plothist, (cat_sex[where(cat_sex.mag_auto GT sex_mag_range[0] and cat_sex.mag_auto LT sex_mag_range[1] AND cat_sex.flux_radius GT sex_flux_radius_min AND cat_sex.flags LT 8, n_gv)]).flux_radius, temp_xhist, temp_yhist, bin=sex_radius_bin, /noplot
			temp=max(temp_yhist, gv) & sex_radius=temp_xhist[gv]
			sex_radius=median((cat_sex[where(cat_sex.mag_auto GT sex_mag_range[0] AND cat_sex.mag_auto LT sex_mag_range[1] AND cat_sex.flux_radius GT sex_radius*0.9 AND cat_sex.flux_radius LT sex_radius*1.1 , n_gv)]).flux_radius)
			print, 'Flux radius for chip ', j, sex_radius
			flux_radius_median.add, sex_radius
		endforeach
		flux_radius_median=median(flux_radius_median.toarray(type='FLOAT'))
		
		for j=0L, input_calib[i].n_chip-1 do begin
			im_h=headfits(input_calib[i].im_file, ext=j+1)
			extast, im_h, im_ast
			xy2ad, im_ast.naxis[0]/2., im_ast.naxis[1]/2., im_ast, pipe_ra, pipe_dec
			pipe_airmass=tai2airmass(pipe_ra,pipe_dec,2000., mjd=input_calib[i].mjd, lon=-70.8065, latitude=-30.1692, alt=2243.31)

			cat_sex=mrdfits(input_calib[i].sex_cat_file, 2*(j+1), cat_sex_h, COLUMNS=['NUMBER','ALPHA_J2000','DELTA_J2000','FLUX_RADIUS','MAG_AUTO','FLUX_AUTO','MAG_APER','FLUX_APER','FLAGS'], /silent)
;			plothist, (cat_sex[where(cat_sex.mag_auto GT sex_mag_range[0] and cat_sex.mag_auto LT sex_mag_range[1] AND cat_sex.flux_radius GT sex_flux_radius_min AND cat_sex.flags LT 8, n_gv)]).flux_radius, temp_xhist, temp_yhist, bin=sex_radius_bin, /noplot
;			temp=max(temp_yhist, gv) & sex_radius=temp_xhist[gv]
			sex_radius=flux_radius_median
			sex_radius=median((cat_sex[where(cat_sex.mag_auto GT sex_mag_range[0] AND cat_sex.mag_auto LT sex_mag_range[1] AND cat_sex.flux_radius GT sex_radius*0.9 AND cat_sex.flux_radius LT sex_radius*1.1 , n_gv)]).flux_radius)
			gv_sex=where(cat_sex.mag_auto GT sex_mag_range[0] AND cat_sex.mag_auto LT sex_mag_range[1] AND cat_sex.flux_radius GT sex_radius*0.9 AND cat_sex.flux_radius LT sex_radius*1.1 AND cat_sex.flags LE 1, n_gv_sex)

			gv_plot=where(cat_sex.flags LT 8, n_gv_plot)
			plot, [cat_sex[gv_plot].flux_radius], [cat_sex[gv_plot].mag_auto], psym=1, xrange=plot_flux_radius_range, yrange=plot_mag_range, /ystyle, /xstyle
			oplot, [cat_sex[gv_sex].flux_radius], [cat_sex[gv_sex].mag_auto], psym=1, color=200
			oplot, sex_radius*[1,1], [0,100], color=200
			oplot, sex_radius*[0.9,0.9], [0,100], line=2, color=200
			oplot, sex_radius*[1.1,1.1], [0,100], line=2, color=200
			oplot, [0,100], sex_mag_range[0]*[1,1], line=2, color=100
			oplot, [0,100], sex_mag_range[1]*[1,1], line=2, color=100
			wait, 0.2
			if do_debug then wait, 2

		 	gv_match= where( min( (cat_sex[gv_sex].alpha_j2000#make_array(n_gv_ref,value=1,/double)-make_array(n_gv_sex,value=1.,/double)#cat_ref[*].ra)^2 + (cat_sex[gv_sex].delta_j2000#make_array(n_gv_ref,value=1,/double)-make_array(n_gv_sex,value=1.,/double)#cat_ref[*].dec)^2, id_match, DIM=1) LT (0.6/3600.)^2, n_gv_match)
		 	print, 'Number of matched stars ', n_gv_match

			if n_gv_match GT 2 then begin
				print, 'Using matched stars for computing the zp'
				gv_sex_match=gv_sex[(id_match[gv_match] mod n_gv_sex)]
				gv_ref_match=(id_match[gv_match]/n_gv_sex)

				case input_calib[i].filter of
					'u': begin
						gv=where(cat_ref[gv_ref_match].mag[0] GT 0 and cat_sex[gv_sex_match].flux_aper[4] GT 0., n_gv)

						mag_diff=cat_ref[gv_ref_match[gv]].mag[0] + 2.5*alog10( cat_sex[gv_sex_match[gv]].flux_aper[4]/im_exptime)
						pipe_zp=biweight_mean(mag_diff, pipe_zperr)
						plot, cat_ref[gv_ref_match].mag[0], mag_diff, psym=1, xrange=sex_mag_range+[0,2], yrange=median(mag_diff)+[-1.,1.]
						oplot, [0,100], pipe_zp*[1,1], color=200, line=2
						print, 'Pipeline zero point ', string(pipe_zp, pipe_zperr, FORMAT='(F0.2,4X,F0.2)')

						for k=0L, n_gv-1 do begin
							printf, lun_zp_full, j+1, cat_sex[gv_sex_match[gv[k]]].number, cat_sex[gv_sex_match[gv[k]]].alpha_j2000, cat_sex[gv_sex_match[gv[k]]].delta_j2000, cat_sex[gv_sex_match[gv[k]]].flux_aper[4], cat_sex[gv_sex_match[gv[k]]].flux_auto, cat_ref[gv_ref_match[gv[k]]].mag[0], mag_diff[k], pipe_airmass, FORMAT='(I4,4X,I4,4X,F0.6,4X,F0.6,4X,F15.3,4X,F15.3,4X,F0.3,4X,F0.3,4X,F0.3)'
						endfor
						end
					'g': begin
						gv=where(cat_ref[gv_ref_match].mag[1] GT 0 and cat_sex[gv_sex_match].flux_aper[4] GT 0., n_gv)

						mag_diff=cat_ref[gv_ref_match[gv]].mag[1] + 2.5*alog10( cat_sex[gv_sex_match[gv]].flux_aper[4]/im_exptime)
						pipe_zp=biweight_mean(mag_diff, pipe_zperr)
						plot, cat_ref[gv_ref_match].mag[1], mag_diff, psym=1, xrange=sex_mag_range+[0,2], yrange=median(mag_diff)+[-1.,1.]
						oplot, [0,100], pipe_zp*[1,1], color=200, line=2
						print, 'Pipeline zero point ', string(pipe_zp, pipe_zperr, FORMAT='(F0.2,4X,F0.2)')

						for k=0L, n_gv-1 do begin
							printf, lun_zp_full, j+1, cat_sex[gv_sex_match[gv[k]]].number, cat_sex[gv_sex_match[gv[k]]].alpha_j2000, cat_sex[gv_sex_match[gv[k]]].delta_j2000, cat_sex[gv_sex_match[gv[k]]].flux_aper[4], cat_sex[gv_sex_match[gv[k]]].flux_auto, cat_ref[gv_ref_match[gv[k]]].mag[1], mag_diff[k], pipe_airmass, FORMAT='(I4,4X,I4,4X,F0.6,4X,F0.6,4X,F15.3,4X,F15.3,4X,F0.3,4X,F0.3,4X,F0.3)'
						endfor
						end
					'r': begin
						gv=where(cat_ref[gv_ref_match].mag[2] GT 0 and cat_sex[gv_sex_match].flux_aper[4] GT 0., n_gv)

						mag_diff=cat_ref[gv_ref_match[gv]].mag[2] + 2.5*alog10( cat_sex[gv_sex_match[gv]].flux_aper[4]/im_exptime)
						pipe_zp=biweight_mean(mag_diff, pipe_zperr)
						plot, cat_ref[gv_ref_match].mag[2], mag_diff, psym=1, xrange=sex_mag_range+[0,2], yrange=median(mag_diff)+[-1.,1.]
						oplot, [0,100], pipe_zp*[1,1], color=200, line=2
						print, 'Pipeline zero point ', string(pipe_zp, pipe_zperr, FORMAT='(F0.2,4X,F0.2)')

						for k=0L, n_gv-1 do begin
							printf, lun_zp_full, j+1, cat_sex[gv_sex_match[gv[k]]].number, cat_sex[gv_sex_match[gv[k]]].alpha_j2000, cat_sex[gv_sex_match[gv[k]]].delta_j2000, cat_sex[gv_sex_match[gv[k]]].flux_aper[4], cat_sex[gv_sex_match[gv[k]]].flux_auto, cat_ref[gv_ref_match[gv[k]]].mag[2], mag_diff[k], pipe_airmass, FORMAT='(I4,4X,I4,4X,F0.6,4X,F0.6,4X,F15.3,4X,F15.3,4X,F0.3,4X,F0.3,4X,F0.3)'
						endfor
						end
					'i': begin
						gv=where(cat_ref[gv_ref_match].mag[3] GT 0 and cat_sex[gv_sex_match].flux_aper[4] GT 0., n_gv)

						mag_diff=cat_ref[gv_ref_match[gv]].mag[3] + 2.5*alog10( cat_sex[gv_sex_match[gv]].flux_aper[4]/im_exptime)
						pipe_zp=biweight_mean(mag_diff, pipe_zperr)
						plot, cat_ref[gv_ref_match].mag[3], mag_diff, psym=1, xrange=sex_mag_range+[0,2], yrange=median(mag_diff)+[-1.,1.]
						oplot, [0,100], pipe_zp*[1,1], color=200, line=2
						print, 'Pipeline zero point ', string(pipe_zp, pipe_zperr, FORMAT='(F0.2,4X,F0.2)')

						for k=0L, n_gv-1 do begin
							printf, lun_zp_full, j+1, cat_sex[gv_sex_match[gv[k]]].number, cat_sex[gv_sex_match[gv[k]]].alpha_j2000, cat_sex[gv_sex_match[gv[k]]].delta_j2000, cat_sex[gv_sex_match[gv[k]]].flux_aper[4], cat_sex[gv_sex_match[gv[k]]].flux_auto, cat_ref[gv_ref_match[gv[k]]].mag[3], mag_diff[k], pipe_airmass, FORMAT='(I4,4X,I4,4X,F0.6,4X,F0.6,4X,F15.3,4X,F15.3,4X,F0.3,4X,F0.3,4X,F0.3)'
						endfor
						end
					'z': begin
						gv=where(cat_ref[gv_ref_match].mag[4] GT 0 and cat_sex[gv_sex_match].flux_aper[4] GT 0., n_gv)

						mag_diff=cat_ref[gv_ref_match[gv]].mag[4] + 2.5*alog10( cat_sex[gv_sex_match[gv]].flux_aper[4]/im_exptime)
						pipe_zp=biweight_mean(mag_diff, pipe_zperr)
						plot, cat_ref[gv_ref_match].mag[4], mag_diff, psym=1, xrange=sex_mag_range+[0,2], yrange=median(mag_diff)+[-1.,1.]
						oplot, [0,100], pipe_zp*[1,1], color=200, line=2
						print, 'Pipeline zero point ', string(pipe_zp, pipe_zperr, FORMAT='(F0.2,4X,F0.2)')

						for k=0L, n_gv-1 do begin
							printf, lun_zp_full, j+1, cat_sex[gv_sex_match[gv[k]]].number, cat_sex[gv_sex_match[gv[k]]].alpha_j2000, cat_sex[gv_sex_match[gv[k]]].delta_j2000, cat_sex[gv_sex_match[gv[k]]].flux_aper[4], cat_sex[gv_sex_match[gv[k]]].flux_auto, cat_ref[gv_ref_match[gv[k]]].mag[4], mag_diff[k], pipe_airmass, FORMAT='(I4,4X,I4,4X,F0.6,4X,F0.6,4X,F15.3,4X,F15.3,4X,F0.3,4X,F0.3,4X,F0.3)'
						endfor
						end
					else: stop
				endcase

				printf, lun_zp, j+1, im_zp-2.5*alog10(im_exptime), pipe_airmass, pipe_zp, pipe_zperr, n_gv, im_fwhm, FORMAT='(I4,4X,F0.3,4X,F0.3,4X,F0.3,4X,F0.3,4X,I4,4X,F0.2)'
      		endif
     
		endfor
		print, 'DECam zero point ', im_zp
		free_lun, lun_zp
		free_lun, lun_zp_full
		
	endfor

	if do_zp_phot then gv_phot=where(input_calib.airmass LT 2. AND input_calib.mjd GT 0. AND input_calib.skip EQ 'F' AND input_calib.photflag EQ 'T', n_gv_phot) $
	else gv_phot=where(input_calib.airmass LT 2. AND input_calib.mjd GT 0. AND input_calib.skip EQ 'F', n_gv_phot)

	delvarx, pipe_zp, coeff
	for i=0L, n_gv_phot-1 do begin
		print, 'Processing file ', input_calib[gv_phot[i]].sex_zp_file

		readcol, input_calib[gv_phot[i]].sex_zp_file, temp_ext, temp_decam_zp, temp_airmass, temp_zp, temp_zperr, temp_zp_nstars, temp_fwhm, FORMAT='I,F,F,F,F,I,F', COUNT=n_lines
		if n_lines GT 0 then begin
			create_struct, temp_struct, '', ['ext','zp','airmass','filter','mjd','photflag'], 'I,F,F,A,F,A', dim=n_elements(temp_ext)
			temp_struct.ext=temp_ext
			temp_struct.zp=temp_zp
			temp_struct.airmass=temp_airmass
			temp_struct.filter=input_calib[gv_phot[i]].filter
			temp_struct.mjd=input_calib[gv_phot[i]].mjd
			temp_struct.photflag=input_calib[gv_phot[i]].photflag
			if size(pipe_zp, /N_ELEMENTS) EQ 0L then pipe_zp=temp_struct $
			else pipe_zp=[pipe_zp,temp_struct]
		endif

		if n_lines GT 0 then begin
			readcol, input_calib[gv_phot[i]].sex_zp_full_file, temp_ext, temp_id, temp_ra, temp_dec, temp_flux_aper, temp_flux_auto, temp_mag_ref, temp_zp, temp_airmass, FORMAT='I,I,F,F,D,D,F,F,F', COUNT=n_lines
			create_struct, temp_struct, '', ['ext','zp','airmass','filter','mjd'], 'I,F,F,A,F', dim=n_elements(temp_ext)
			temp_struct.ext=temp_ext
			temp_struct.zp=temp_zp
			temp_struct.airmass=temp_airmass
			temp_struct.filter=input_calib[gv_phot[i]].filter
			temp_struct.mjd=input_calib[gv_phot[i]].mjd
			if size(pipe_zp_full, /N_ELEMENTS) EQ 0L then pipe_zp_full=temp_struct $
			else pipe_zp_full=[pipe_zp_full,temp_struct]
		endif

	endfor

	filter_uniq=pipe_zp[uniq(pipe_zp.filter, sort(pipe_zp.filter))].filter
	mjd_uniq=floor(pipe_zp[uniq(floor(pipe_zp.mjd+2./24), sort(floor(pipe_zp.mjd+2./24)))].mjd)

	openw, lun, 'survey_zp_'+do_program+'.dat', /get_lun
	printf, lun, '# MJD		filter		zp		k		zp_err		k_err		PHOTFLAG'

	for i=0L, n_elements(filter_uniq)-1 do begin
		for j=0L, n_elements(mjd_uniq)-1 do begin

			gv_pipe=where(pipe_zp.filter EQ filter_uniq[i] AND floor(pipe_zp.mjd+2./24) EQ mjd_uniq[j], n_gv_pipe)
			if do_zp_phot then begin
				gv_pipe_phot=where(pipe_zp.filter EQ filter_uniq[i] AND floor(pipe_zp.mjd+2./24) EQ mjd_uniq[j] AND pipe_zp.photflag EQ 'T', n_gv_pipe_phot)
				coeff=robust_linefit(pipe_zp[gv_pipe_phot].airmass, pipe_zp[gv_pipe_phot].zp, temp_yfit, temp_sigma, coeff_error)
				photflag = coeff_error[0] LT 0.015 ? 'T' : 'F'
				printf, lun, mjd_uniq[j], filter_uniq[i], coeff[0], coeff[1], coeff_error[0], coeff_error[1], photflag, FORMAT='(F0.1,4X,A,4X,F0.2,4X,F0.2,4X,F0.2,4X,F0.2,4X,A)'
			endif $
			else begin
				gv_pipe_phot=where(pipe_zp.filter EQ filter_uniq[i] AND floor(pipe_zp.mjd+2./24) EQ mjd_uniq[j], n_gv_pipe_phot)
				coeff=robust_linefit(pipe_zp[gv_pipe_phot].airmass, pipe_zp[gv_pipe_phot].zp, temp_yfit, temp_sigma, coeff_error)
				photflag = coeff_error[0] LT 0.05 ? 'T' : 'F'
				printf, lun, mjd_uniq[j], filter_uniq[i], coeff[0], coeff[1], coeff_error[0], coeff_error[1], photflag, FORMAT='(F0.1,4X,A,4X,F0.2,4X,F0.2,4X,F0.2,4X,F0.2,4X,A)'
			endelse
			
			cgloadct, 0
			cgwindow, wxsize=800, wysize=600
			cgplot, [0], [0], xrange=[1.,1.9], yrange=mean(pipe_zp[gv_pipe].zp)+[-1,1], /nodata, /window, xtitle='airmass', ytitle='mag_std - mag_inst', title='Filter '+filter_uniq[i]+' - Date '+date_conv(mjd_uniq[j]+2400000.5D,'S')
			cgplot, pipe_zp[gv_pipe].airmass, pipe_zp[gv_pipe].zp, psym=cgsymcat('filled circle'), symsize=0.4, color='blue', /over, /addcmd
			if do_zp_phot then cgplot, pipe_zp[gv_pipe_phot].airmass, pipe_zp[gv_pipe_phot].zp, psym=cgsymcat('filled circle'), symsize=0.4, color='red', /over, /addcmd
			cgplot, [0.,100.], coeff[0]+coeff[1]*[0.,100.], color='black', line=2, /over, /addcmd
			cglegend, color='red', align=0, location=[0.14,0.86], length=0.02, titles=filter_uniq[i]+'-band coeff='+string(coeff,format='(F0.3,",",F0.3)')+' coeff_error='+string(coeff_error,format='(F0.3,",",F0.3)'), vspace=2., /addcmd
			cgcontrol, output='results/standard_zp_airmass_uncorrected_'+filter_uniq[i]+'_'+string(mjd_uniq[j],FORMAT='(I0)')+'.pdf'

			print, 'Airmass term compute - Filter ', filter_uniq[i], ' Date ', date_conv(mjd_uniq[j]+2400000.5D,'S')
			print, coeff[0], coeff_error[0], coeff[1], coeff_error[1], FORMAT='("zp = ",F0.3," +- ",F0.3," / k = ", F0.3, " +- ", F0.3)'

			gv_pipe=where(input_calib.filter EQ filter_uniq[i] AND floor(input_calib.mjd+2./24) EQ mjd_uniq[j], n_gv_pipe)
			if do_zp_phot then gv_pipe_phot=where(input_calib.filter EQ filter_uniq[i] AND floor(input_calib.mjd+2./24) EQ mjd_uniq[j] AND input_calib.photflag EQ 'T', n_gv_pipe_phot)
			cgloadct, 0
			cgwindow, wxsize=800, wysize=600
			cgplot, [0], [0], xrange=mjd_uniq[j]+[-2,12]/24., yrange=[1.,2.0], /nodata, /window, xtitle='MJD', ytitle='airmass', title='Filter '+filter_uniq[i]+' - Date '+date_conv(mjd_uniq[j]+2400000.5D,'S'), XTICKFORMAT='(F0.1)', XCHARSIZE=1
			cgplot, input_calib[gv_pipe].mjd, input_calib[gv_pipe].airmass, psym=cgsymcat('filled circle'), symsize=1., color='blue', /over, /addcmd
			if do_zp_phot then cgplot, input_calib[gv_pipe_phot].mjd, input_calib[gv_pipe_phot].airmass, psym=cgsymcat('filled circle'), symsize=1., color='red', /over, /addcmd
			cgcontrol, output='results/standard_airmass_mjd_'+filter_uniq[i]+'_'+string(mjd_uniq[j],FORMAT='(I0)')+'.pdf'

			cgdelete, /all

		endfor
	endfor

	free_lun, lun

	cgdelete, /all

	create_struct, pipe_phot, '', ['filter','zp','k'], 'A,F,F', dim=n_elements(filter_uniq)

endif else $
if recipe EQ 'sextractor' then begin

	if do_overwrite then begin 
		file_delete, input_target.sex_cat_file, /noexpand, /allow_non, /quiet
	  file_delete, input_target.scamp_head_file, /noexpand, /allow_non, /quiet
	endif

	if n_elements(bridges) eq 0 then bridges = build_bridges(do_n_cpu)
	n_cpu = n_elements(bridges)

	for k=0L, n_cpu-1 do begin
		(bridges[k])->execute, '.r worker'
		(bridges[k])->execute, '.r callback'
		(bridges[k])->setproperty, callback='callback'
	endfor

	for i=0L, n_elements(input_target)-1 do begin
			
		print, 'Sexractor - Iteration '+strn(i+1)+' of '+strn(n_elements(input_target))
		print, 'Checking if catalog exists ', input_target[i].sex_cat_file
;		if input_target[i].sex_cat_file EQ '/Volumes/Data2/Magellanic_stream/DECam/20130808/pipeline/sextractor/ms_t9_d5_i.ldac' then stop

		if file_test(input_target[i].sex_cat_file, /regular) then fits_info, input_target[i].sex_cat_file, n_ext=sex_n_ext, /silent $
		else sex_n_ext=0

		if sex_n_ext NE input_target[i].n_chip*2 OR do_overwrite EQ 1 then begin

			case input_target[i].filter of
				'u': sex_satur_level='25000.'
				'g': sex_satur_level='25000.'
				'i': sex_satur_level='30000.'
				'z': sex_satur_level='30000.'
				else: stop
			endcase
		

			if n_elements(input_zp) GT 0 then begin
				gv=where(input_zp.filter EQ input_target[i].filter AND input_zp.mjd EQ input_target[i].mjd_floor, n_gv)
				if n_gv GT 0 then sex_zp=input_zp[gv].zp $
				else stop
			endif else begin
				sex_zp=input_target[i].zp
			endelse
			;sex_zp=input_target[i].zp

			command = 'sex '+input_target[i].im_file +' -c sex_config/ctio_decam.sex -CATALOG_NAME '+input_target[i].sex_cat_file+' -WEIGHT_IMAGE '+input_target[i].weight_file+' -XML_NAME '+input_target[i].sex_xml_file+' -SATUR_LEVEL '+sex_satur_level+' -MAG_ZEROPOINT '+string(sex_zp,FORMAT='(F0.2)')+' -CHECKIMAGE_TYPE BACKGROUND -CHECKIMAGE_NAME '+input_target[i].sex_check_file+' -PHOT_APERTURES '+ string(2*3*input_target[i].fwhm,FORMAT='(F0.2)')

;			ud = {xrange:temp_xrange[ii:ii+1], yrange:temp_yrange[jj:jj+1], pout:pout, im_size:im_size}
			bridge = get_idle_bridge(bridges)
;			bridge->setproperty, userdata=ud
			bridge->setvar, 'command', command

			print, 'Sexractor - Running command ', command
			bridge->execute, /nowait, 'worker, command, out'


;			command = 'sex '+input_target[i].im_file +' -c sex_config/ctio_decam.sex -CATALOG_NAME '+input_target[i].sex_cat_file+' -WEIGHT_IMAGE '+input_target[i].weight_file+' -XML_NAME '+input_target[i].sex_xml_file+' -SATUR_LEVEL '+sex_satur_level+' -MAG_ZEROPOINT '+string(sex_zp,FORMAT='(F0.2)')+' -CHECKIMAGE_TYPE BACKGROUND -CHECKIMAGE_NAME '+input_target[i].sex_check_file
;			print, command
;			spawn, command

		endif

	endfor

	barrier_bridges, bridges

	burn_bridges, bridges
	temp=temporary(bridges)

endif else $
if recipe EQ 'iq' then begin

	window, 0, XSIZE=800, YSIZE=600
	window, 1, XSIZE=800, YSIZE=600


	for i=0L, n_elements(do_program)-1 do begin
		log_file='survey_iq_'+do_program[i]+'.dat'
		if file_test(log_file) EQ 0 OR do_overwrite then begin
			openw, lun, log_file, /get_lun
			printf, lun, '#  im_file        filter  tile  dither  sky_level   FWHM'
			free_lun, lun
		endif
	endfor

	for i=0L, n_elements(input_target)-1 do begin
		print, 'IQ - Iteration '+strn(i+1)+' of '+strn(n_elements(input_target))
		print, 'IQ - Analyzing image quality of file '+input_target[i].im_file

		log_file='survey_iq_'+input_target[i].program+'.dat'
		log_file_copy=input_dir+'/'+input_target[i].program+'/survey_iq_'+input_target[i].program+'.dat'
		readcol, log_file, temp_im_file, temp_filter, temp_tile, temp_dither, temp_sky, temp_fwhm, FORMAT='A,A,A,A,F,F', /silent
		if n_elements(temp_im_file) GT 0 then begin		
			gv_im_file=where(temp_im_file EQ input_target[i].im_file, n_gv_im_file)
			if n_gv_im_file GT 0 then continue
		endif

		vig_diam=101
		sex_mag_bin=0.5
		sex_radius_bin=0.2
		sex_flux_radius_min=1.5
		sex_ellipticity_min=0.8
		sex_flag_max=1
	
		case input_target[i].filter of
			'u': begin
				sex_mag_range=[10.,18.]
				plot_mag_range=[21,10]
			end
			'g': begin
				sex_mag_range=[15.,21.]
				plot_mag_range=[25,15]
			end
			'i': begin
				sex_mag_range=[15.,21.]
				plot_mag_range=[25,15]
			end
			'z': begin
				sex_mag_range=[15.,21.]
				plot_mag_range=[25,15]
			end
			else: stop
		endcase

		if file_test(input_target[i].sex_cat_file) EQ 0 then begin
			print, "IQ error - Please run decam_pipeline, 'sextractor'"
			stop
		endif

		iq_sky=list()
		iq_fwhm=list()

		foreach j, input_chip_iq do begin

			print
			print, 'IQ - Analyzing chip ', j

			im_data=readfits(input_target[i].im_file, im_h, ext=j, /silent)
			im_size=size(im_data, /dim)
			im_gain=fxpar(im_h, 'GAINA')
			im_ron=fxpar(im_h, 'RDNOISEA')

			im_sky=median(im_data[im_size[0]/2-300:im_size[0]/2+300, im_size[1]/2-300:im_size[1]/2+300])
			cat_sex=mrdfits(input_target[i].sex_cat_file, 2*j, cat_sex_h, COLUMNS=['NUMBER','X_IMAGE','Y_IMAGE','FLUX_RADIUS','MAG_AUTO','MAGERR_AUTO','FLAGS','A_IMAGE','B_IMAGE'], /silent)
			if size(cat_sex, /type) EQ 8 then gv_stars=where(cat_sex.mag_auto GT sex_mag_range[0] and cat_sex.mag_auto LT sex_mag_range[1] AND cat_sex.flux_radius GT sex_flux_radius_min AND cat_sex.flags LE 3 AND cat_sex.b_image/cat_sex.a_image GT sex_ellipticity_min, n_gv_stars) $
			else n_gv_stars=0

			if n_elements(cat_sex) GT 100 AND n_gv_stars GT 10 AND im_sky LT 30000  then begin
				wset, 0
				plot, cat_sex.flux_radius, cat_sex.mag_auto, psym=1, xrange=[1,6], yrange=plot_mag_range
	
;				if n_gv_stars LT 10 then begin
;					print, 'IQ - Error, there is not enough number of stars'
;					stop
;				endif
	
				oplot, cat_sex[gv_stars].flux_radius, cat_sex[gv_stars].mag_auto, psym=1, color=100
				oplot, [0,100], sex_mag_range[0]*[1,1], line=2, color=100
				oplot, [0,100], sex_mag_range[1]*[1,1], line=2, color=100
	
				plothist, cat_sex[gv_stars].flux_radius, temp_xhist, temp_yhist, bin=sex_radius_bin, /noplot
				temp=max(temp_yhist, gv) & sex_radius=temp_xhist[gv]
				sex_radius=median((cat_sex[where(cat_sex.mag_auto GT sex_mag_range[0] AND cat_sex.mag_auto LT sex_mag_range[1] AND cat_sex.flux_radius GT sex_radius*0.9 AND cat_sex.flux_radius LT sex_radius*1.1 , n_gv)]).flux_radius)
				print, j+1, sex_radius, FORMAT='("Chip ",I2,"  Flux_radius ",F0.1)'
	
				oplot, sex_radius*[1,1], [0,100], color=200
				oplot, sex_radius*[0.9,0.9], [0,100], line=2, color=200
				oplot, sex_radius*[1.1,1.1], [0,100], line=2, color=200
	
				gv_stars=where(cat_sex.mag_auto GT sex_mag_range[0] AND cat_sex.mag_auto LT sex_mag_range[1]-1 AND cat_sex.flux_radius GT sex_radius*0.9 AND cat_sex.flux_radius LT sex_radius*1.1 AND cat_sex.flags LE sex_flag_max AND cat_sex.x_image GT vig_diam/2. AND cat_sex.x_image LT (im_size[0]-vig_diam/2.) AND cat_sex.y_image GT vig_diam/2. AND cat_sex.y_image LT (im_size[1]-vig_diam/2.), n_gv_stars)
				oplot, cat_sex[gv_stars].flux_radius, cat_sex[gv_stars].mag_auto, psym=1, color=200
		
			  if n_gv_stars GT 0 then begin
					gv_stars=gv_stars[sort(cat_sex[gv_stars].mag_auto)]
		
					temp_n=n_gv_stars<40
		  	  temp_x=cat_sex[gv_stars].x_image
		   		temp_y=cat_sex[gv_stars].y_image
	
					im_fwhm=list()	
					im_sky=list()	
	
		  	  for k=0L, temp_n-1 do begin
			      x_range=[ floor(temp_x[k]-(vig_diam-1.)/2)>0, ceil(temp_x[k]+(vig_diam-1.)/2)<(im_size[0]-1) ]
		        y_range=[ floor(temp_y[k]-(vig_diam-1.)/2)>0, ceil(temp_y[k]+(vig_diam-1.)/2)<(im_size[1]-1) ]
		     	  vig_data = im_data[x_range[0]:x_range[1],y_range[0]:y_range[1]]; - im_sky
						vig_size=size(vig_data, /dim)
		
			      x_center_range=[ floor((vig_diam-1.)/4), ceil(-1-(vig_diam-1.)/4) ]
		    	  y_center_range=[ floor((vig_diam-1.)/4), ceil(-1-(vig_diam-1.)/4) ]
						vig_center_data=vig_data[x_center_range[0]:x_center_range[1],y_center_range[0]:y_center_range[1]]
						vig_center_size=size(vig_center_data, /dim)
		
			      im_max=max(vig_center_data, gv_max)
		       	im_c= [gv_max mod vig_center_size[0], gv_max/vig_center_size[0]]
						im_c += [x_center_range[0], y_center_range[0]]
	
						dist_circle, vig_mask, vig_size, im_c[0], im_c[1]
						vig_sky=median(vig_data[where(vig_mask GT 20., n_sky)])
						vig_mag=25.-2.5*alog10(total(vig_data[where(vig_mask LE 20., n_star)]) - n_star*vig_sky )
	
						vig_res=abs( vig_data-vig_sky - max(vig_center_data-vig_sky)/2. )
						gv=sort(vig_res*vig_mask^2)
						vig_fwhm=2*median(vig_mask[gv[1:5]])
;						print, 'Radius for computing FWHM ', vig_mask[gv[1:5]]
						wset, 1
						plot, vig_mask, vig_data-vig_sky, xrange=[0,30], psym=1
						oplot, [0,100], max(vig_data-vig_sky)/2.*[1,1], line=2, color=200
						oplot, vig_fwhm/2.*[1,1], [-1e5,1e5], line=2, color=200
	
						vig_skyrad=4*vig_fwhm < 40.
						vig_psfrad=3*vig_fwhm < 50.
						vig_fitrad=vig_fwhm < 50.
		
						gcntrd, vig_data, im_c[0], im_c[1], im_cx, im_cy, vig_fwhm
						dist_circle, vig_mask, vig_size, im_cx, im_cy
						vig_sky=median(vig_data[where(vig_mask GT vig_skyrad, n_sky)])
						vig_mag=25.-2.5*alog10(total(vig_data[where(vig_mask LE vig_skyrad, n_star)]) - n_star*vig_sky )
						plot, vig_mask, vig_data-vig_sky, xrange=[0,30], psym=1
						oplot, [0,100], max(vig_data-vig_sky)/2.*[1,1], line=2, color=200
						oplot, vig_fwhm/2.*[1,1], [-1e5,1e5], line=2, color=200
	
			      fit_x_range=[ floor(im_cx- 2*vig_fwhm)>0, ceil(im_cx+ 2*vig_fwhm)<(vig_size[0]-1)  ]
		    	  fit_y_range=[ floor(im_cy- 2*vig_fwhm)>0, ceil(im_cy+ 2*vig_fwhm)<(vig_size[1]-1)  ]
						fit_vig_data=vig_data[fit_x_range[0]:fit_x_range[1],fit_y_range[0]:fit_y_range[1]]
						fit_vig_size=size(fit_vig_data, /dim)
			      fit_peak=max(fit_vig_data-vig_sky, gv_max)
		       	im_c= [gv_max mod fit_vig_size[0], gv_max/fit_vig_size[0]]
						gcntrd, fit_vig_data, im_c[0], im_c[1], fit_im_cx, fit_im_cy, vig_fwhm
	
	          fit_im_cx=fit_im_cx + (floor(im_cx- 2*vig_fwhm)-floor(im_cx- 1.5*vig_fwhm))
	          fit_im_cy=fit_im_cy + (floor(im_cy- 2*vig_fwhm)-floor(im_cy- 1.5*vig_fwhm))
	          fit_x_range=[ floor(im_cx- 1.5*vig_fwhm)>0, ceil(im_cx+ 1.5*vig_fwhm)<(vig_size[0]-1) ]
	          fit_y_range=[ floor(im_cy- 1.5*vig_fwhm)>0, ceil(im_cy+ 1.5*vig_fwhm)<(vig_size[1]-1) ]
	          fit_vig_data=vig_data[fit_x_range[0]:fit_x_range[1],fit_y_range[0]:fit_y_range[1]]
	
	;					fit_sigma_range=vig_fwhm/2.35482*[0.8,1.4]
	;					fit_ellipticity_range=[0.95,1.05]
	;					fit_peak_range=[0.9,1.1]
	;					fit_niter=30
	;					fit_res=fltarr(fit_niter*[1,1,1])
	;
	;					fit_peak=(fit_peak_range[0] + findgen(fit_niter)/(fit_niter-1.)*(fit_peak_range[1]-fit_peak_range[0]))*max(fit_vig_data-vig_sky)
	;					fit_sigma=fit_sigma_range[0] + findgen(fit_niter)/(fit_niter-1.)*(fit_sigma_range[1]-fit_sigma_range[0])
	;					fit_ellipticity=fit_ellipticity_range[0] + findgen(fit_niter)/(fit_niter-1.)*(fit_ellipticity_range[1]-fit_ellipticity_range[0])
	;
	;					for ii=0L, fit_niter-1 do begin
	;						for jj=0L, fit_niter-1 do begin
	;							for kk=0L, fit_niter-1 do begin
	;								psf_param=[ [fit_peak[ii]*[1,1]],[[fit_im_cx,fit_im_cy]],[fit_sigma[jj]*[fit_ellipticity[kk],1]] ]
	;								fit_res[ii,jj,kk]=total( abs(fit_vig_data-vig_sky - psf_gaussian(psf_param, npix=fit_vig_size)) )
	;							endfor
	;						endfor
	;					endfor
	;					temp=min(fit_res, gv_min)
	;					gv_min=array_indices(fit_res, gv_min)
	;					psf_param1=[ fit_peak[gv_min[0]], fit_im_cx, fit_im_cy, fit_sigma[gv_min[1]]*[fit_ellipticity[gv_min[2]],1] ]
	
;	          parinfo = replicate({value:0., fixed:0, limited:[0,0], limits:[0.,0.]}, 7)
;	          parinfo[*].value = [0.,fit_peak, vig_fwhm, vig_fwhm, fit_im_cx, fit_im_cy, 0.]
	          parinfo = replicate({value:0., fixed:0, limited:[0,0], limits:[0.,0.]}, 8)
	          parinfo[*].value = [0.,fit_peak, vig_fwhm, vig_fwhm, fit_im_cx, fit_im_cy, 0., 2.]
	          parinfo[0].fixed = 1
;	          temp=mpfit2dpeak(fit_vig_data-vig_sky, psf_param1, /GAUSSIAN, /TILT, ESTIMATES=parinfo.value, PARINFO=parinfo)
	          temp=mpfit2dpeak(fit_vig_data-vig_sky, psf_param2, /MOFFAT, /TILT, ESTIMATES=parinfo.value, PARINFO=parinfo)
						temp_sigma=sqrt( (psf_param2[2]^2 + psf_param2[3]^2)/2 )

	
;			      im_fwhm.add, 2*sqrt(2*alog(2))*sqrt( (psf_param1[2]^2 + psf_param1[3]^2)/2. )
;			      im_fwhm.add, 2*sqrt(2*alog(2)) * temp_sigma
			      im_fwhm.add, temp_sigma
						im_sky.add, vig_sky
						print
						print, k+1, temp_x[k], temp_y[k], FORMAT='("Star: ",I0," - Coo: ",I0," ",I0)'
						print, 'FWHM and SKY: ', im_fwhm[-1], ' pixels', im_sky[-1], ' ADU'

;						im_cx=fit_x_range[0]+psf_param2[4]
;						im_cy=fit_y_range[0]+psf_param2[5]
;						dist_circle, vig_mask, vig_size, im_cx, im_cy
;						gv=where(vig_mask LE 6., n_gv)

;	          parinfo = replicate({value:0., fixed:0, limited:[0,0], limits:[0.,0.]}, 5)
;	          parinfo[*].value = [fit_peak, 0., vig_fwhm, 2., 0.]
;	          parinfo[4].fixed = 1
;	          temp=mpfitpeak(vig_mask[gv], vig_data[gv]-vig_sky, psf_param3, /MOFFAT, ESTIMATES=parinfo.value, PARINFO=parinfo, NTERMS=5)


						dist_circle, vig_mask, vig_size, im_cx, im_cy
						plot, vig_mask, vig_data-vig_sky, xrange=[0,30], psym=1
						x=findgen(100)/10.
;						y1=gaussian(x,[psf_param1[1], 0., mean(psf_param1[2:3])])
						y2=psf_param2[0] + psf_param2[1]/( (x/temp_sigma)^2 + 1. )^psf_param2[7]
;						oplot, x, y1, line=2, color=100
						oplot, x, y2, line=2, color=100
						oplot, im_fwhm[-1]/2.*[1,1], [-1e5,1e5], line=2, color=200
						wait, 0.2

							
		    	endfor
					if do_debug then begin
						stop
					endif	

					im_sky=median(im_sky.toarray(type='float'))
					im_fwhm=median(im_fwhm.toarray(type='float'))

	
					iq_sky.add, im_sky
					iq_fwhm.add, im_fwhm
	
				endif	

			endif $
			else begin
				mmm, im_data, temp_sky, temp_sigma, readnoise=im_ron
				iq_sky.add, temp_sky
				iq_fwhm.add, 0.
			endelse

		endforeach

		print, 'Summary of SKY and FWHM'
		print, iq_sky
		print, iq_fwhm
		iq_sky=mean(iq_sky.toarray(type='float'))
		iq_fwhm=mean(iq_fwhm.toarray(type='float'))

		openw, lun, log_file, /get_lun, /append
		printf, lun, input_target[i].im_file, input_target[i].filter, input_target[i].tile, input_target[i].dither, iq_sky, iq_fwhm, FORMAT='(A,4X,A,4X,A,4X,A,4X,F0.1,4X,F0.1)'
		free_lun, lun

		file_copy, log_file, log_file_copy, /overwrite

	endfor

endif else $
if recipe EQ 'scamp' then begin


	scamp_pos_error='1.2'
	scamp_scale_error='1.02'
	scamp_angle_error='0.02'
	scamp_sn_thresholds='30.,60.0'
	scamp_fwhm_thresholds='0.,100.'
	scamp_crossid_radius='4.'
	scamp_distort_degrees='4'
	scamp_astref_catalog = '2MASS' ;'USNO-B1' ;'USNO-B1' ;'2MASS'
	scamp_astref_band = 'Ks';'Rf' ;'DEFAULT' ;'Ks'
	scamp_astrefmag_limits='6.,20.' ;'-99.,18.' 
	scamp_match_resol='0.'

	scamp_refcat_dir = output_scamp_dir[0]+'/refcat'
	scamp_dir_out=output_scamp_dir[0]+'/PROGRAM_'+do_program_plain+'_TILE_'+do_tile_plain+'_FILTER_'+do_filter_plain+'_ORDER_'+scamp_distort_degrees+'_REFCAT_'+scamp_astref_catalog
	scamp_list_file=scamp_dir_out+'/scamp_decam_'+do_program_plain+'_'+do_filter_plain+'.dat'
	scamp_cat_file_out = scamp_dir_out+'/scamp_decam_'+do_program_plain+'_'+do_filter_plain+'.ldac'
	scamp_cat_type_out = 'FITS_LDAC'
	scamp_xml_file = scamp_dir_out+'/scamp_decam_'+do_program_plain+'_'+do_filter_plain+'.xml'
	scamp_check_type = 'SKY_ALL,FGROUPS,DISTORTION,ASTR_INTERROR2D,ASTR_INTERROR1D,ASTR_REFERROR2D,ASTR_REFERROR1D,ASTR_CHI2,PHOT_ERROR'
	scamp_check_file = strjoin(scamp_dir_out+'/'+['sky_all','fgroups','distort','astr_interror2d','astr_interror1d','astr_referror2d','astr_referror1d','astr_chi2','psphot_error'],',')

	scamp_mag_bin=0.5
	scamp_radius_bin=0.2

	if do_overwrite then begin
		file_delete, input_target.scamp_cat_file, /noexpand, /allow_non, /quiet
  		file_delete, input_target.scamp_head_file, /noexpand, /allow_non, /quiet
	endif

	if not file_test(scamp_refcat_dir, /directory) then file_mkdir, scamp_refcat_dir, /noexpand_path
	if not file_test(scamp_dir_out, /directory) then file_mkdir, scamp_dir_out, /noexpand_path
	loadct, 12

	for i=0L, n_elements(input_target)-1 do begin

		case input_target[i].filter of
			'u': begin
				scamp_mag_range=[10.,18.]
				plot_mag_range=[21,10]
  				scamp_sn_thresholds='20.,40.0'
				scamp_distort_degrees='4'
			end
			'g': begin
				scamp_mag_range=[15.,21.]
				plot_mag_range=[25,15]
  			scamp_sn_thresholds='20.,40.0'
				scamp_distort_degrees='4'
			end
			'i': begin
				scamp_mag_range=[15.,21.]
				plot_mag_range=[25,15]
  			scamp_sn_thresholds='20.,40.0'
				scamp_distort_degrees='4'
			end
			'z': begin
				scamp_mag_range=[9.,14.]
				plot_mag_range=[18,7]
  				scamp_sn_thresholds='20.,40.0'
				scamp_distort_degrees='4'
			end
			else: stop
		endcase
	
		if file_test(input_target[i].sex_cat_file) EQ 0 then stop

		if file_test(input_target[i].scamp_cat_file, /regular) then fits_info, input_target[i].scamp_cat_file, n_ext=scamp_n_ext, /silent $
		else scamp_n_ext=0
		if scamp_n_ext EQ input_target[i].n_chip*2 AND do_overwrite EQ 0 then continue

		print, 'SCAMP - Reading file '+input_target[i].sex_cat_file

		fits_open, input_target[i].sex_cat_file, fcb_in  ; Lee la tabla fits sin intervenerla           
		fits_read, fcb_in, cat_data0, cat_h0, exten=0 
;		fxaddpar, cat_h0, 'FILTER', input_target[i].filter
		writefits, input_target[i].scamp_cat_file, cat_data0, cat_h0

		for j=0L, input_target[i].n_chip-1 do begin
			cat_sex=mrdfits(input_target[i].sex_cat_file, 2*(j+1), cat_sex_h, COLUMNS=['NUMBER','FLUX_RADIUS','MAG_AUTO','MAGERR_AUTO','FLAGS','A_IMAGE','B_IMAGE','CLASS_STAR'], /silent)
			gv_stars=where(cat_sex.mag_auto GT scamp_mag_range[0] and cat_sex.mag_auto LT scamp_mag_range[1] AND cat_sex.flux_radius GT 1.5 AND cat_sex.flags LE 3 AND cat_sex.b_image/cat_sex.a_image GT 0.8, n_gv_stars)
			if n_gv_stars LT 10 then begin
				print, 'SCAMP - Error, there is not enough number of stars'
				stop
			endif

			plot, cat_sex.flux_radius, cat_sex.mag_auto, psym=1, xrange=[1,10], yrange=plot_mag_range
			oplot, cat_sex[gv_stars].flux_radius, cat_sex[gv_stars].mag_auto, psym=1, color=100
			oplot, [0,100], scamp_mag_range[0]*[1,1], line=2, color=100
			oplot, [0,100], scamp_mag_range[1]*[1,1], line=2, color=100

			plothist, cat_sex[gv_stars].flux_radius, temp_xhist, temp_yhist, bin=scamp_radius_bin, /noplot
			temp=max(temp_yhist, gv) & scamp_radius=temp_xhist[gv]
			scamp_radius=median((cat_sex[where(cat_sex.mag_auto GT scamp_mag_range[0] AND cat_sex.mag_auto LT scamp_mag_range[1] AND cat_sex.flux_radius GT scamp_radius*0.9 AND cat_sex.flux_radius LT scamp_radius*1.1 , n_gv)]).flux_radius)
			print, j+1, scamp_radius, FORMAT='("Chip ",I2,"  Flux_radius ",F0.1)'

			oplot, scamp_radius*[1,1], [0,100], color=200
			oplot, scamp_radius*[0.9,0.9], [0,100], line=2, color=200
			oplot, scamp_radius*[1.1,1.1], [0,100], line=2, color=200
			wait, 0.2
			if j EQ 0 AND do_debug then begin
				forprint, temp_xhist, temp_yhist, text=2
				print, 'DEBUG - CHECK mag range and flux radius are ok'
				stop
			endif

			gv_stars=where(cat_sex.mag_auto GT scamp_mag_range[0] AND cat_sex.mag_auto LT scamp_mag_range[1] AND cat_sex.flux_radius GT scamp_radius*0.9 AND cat_sex.flux_radius LT scamp_radius*1.1 AND cat_sex.flags LE 1 AND cat_sex.b_image/cat_sex.a_image GT 0.8, n_gv_stars)
			oplot, cat_sex[gv_stars].flux_radius, cat_sex[gv_stars].mag_auto, psym=1, color=200
			wait, 0.2

  			fits_read, fcb_in, cat_data1, cat_h1, exten=2*j+1                                   
			fits_read, fcb_in, cat_data2, cat_h2, exten=2*j+2

;			temp_filter=bytarr(80)
;			temp_filter[0:31] = byte('FILTER  =                    '+input_target[i].filter+' /')
;			cat_data1=reform(cat_data1, [80,n_elements(cat_data1)/80])
;			cat_data1=[reform(cat_data1[*,0:-2], n_elements(cat_data1[*,0:-2])), temp_filter, cat_data1[*,-1]]
			cat_data1=reform(cat_data1, [n_elements(cat_data1),1])
;			fxaddpar, cat_h1, 'NAXIS1', n_elements(cat_data1), 'BYTES PER ROW'
;			fxaddpar, cat_h1, 'TFORM1', string(n_elements(cat_data1),FORMAT='(I0)')+'A  '
;			fxaddpar, cat_h1, 'TDIM1', string(80,n_elements(cat_data1)/80,FORMAT='("(",I0,", ",I0,")")')

			cat_data2=cat_data2[*,[gv_stars]]
			fxaddpar, cat_h2, 'NAXIS2', (size(cat_data2, /dim))[1]

			print, 'SCAMP - Creating file '+input_target[i].scamp_cat_file
  			fits_open, input_target[i].scamp_cat_file, fcb_out, /append
		  	fits_write, fcb_out, cat_data1, cat_h1
			fits_write, fcb_out, cat_data2, cat_h2
  			fits_close, fcb_out

		endfor
  		fits_close, fcb_in
	endfor

	forprint, input_target.scamp_cat_file, textout=scamp_list_file, FORMAT='(A)', /NOCOMMENT

	for i=0L, n_elements(input_target)-1 do begin

		if file_test(input_target[i].scamp_ahead_file) EQ 0 OR do_overwrite OR do_ahead then begin

			im_h=headfits(input_target[i].im_file)
			im_mjd=fxpar(im_h, 'MJD-OBS')
			input_target[i].mjd=im_mjd
			print, input_target[i].im_file, input_target[i].filter, im_mjd, FORMAT='(A,2X,A,2X,F0.4)'

			openw, lun, input_target[i].scamp_ahead_file, /get_lun
			for j=0L, input_target[i].n_chip-1 do begin
				im_h=headfits(input_target[i].im_file, ext=j+1)
				extast, im_h, im_ast
				xy2ad, im_ast.naxis[0]/2., im_ast.naxis[1]/2., im_ast, im_ra, im_dec
				im_airmass=tai2airmass(im_ra, im_dec, 2000., mjd=im_mjd, lon=-70.8065, latitude=-30.1692, alt=2243.31)

				printf, lun, "PROGRAM = '"+input_target[i].program+"'"
				printf, lun, "FILTER  = '"+input_target[i].filter+"'"
				printf, lun, 'AIRMASS =     '+string(im_airmass,FORMAT='(F0.3)')
				printf, lun, 'EXPTIME =     '+string(input_target[i].exptime,FORMAT='(F0.2)')  ;string(600.,FORMAT='(F0.2)')
				if n_elements(input_zp) GT 0 then begin
					gv=where(input_zp.filter EQ input_target[i].filter AND input_zp.mjd EQ input_target[i].mjd_floor, n_gv)
					if n_gv EQ 1 then begin
						printf, lun, 'PHOTFLAG=     '+input_zp[gv].photflag
						printf, lun, 'PHOT_ZP =     '+string(input_zp[gv].zp,FORMAT='(F0.2)')
						printf, lun, 'PHOT_K  =     '+string(input_zp[gv].k,FORMAT='(F0.2)')
					endif	else begin
						gv=where(input_zp.filter EQ input_target[i].filter, n_gv)
						printf, lun, 'PHOTFLAG=     F'
						printf, lun, 'PHOT_ZP =     '+string(input_target[i].zp,FORMAT='(F0.3)')
						printf, lun, 'PHOT_K  =     '+string(input_zp[gv[0]].k,FORMAT='(F0.2)')
					endelse
				endif else begin
					printf, lun, 'PHOTFLAG=     T'
					printf, lun, 'PHOT_ZP =     '+string(input_target[i].zp,FORMAT='(F0.3)')
					printf, lun, 'PHOT_K  =      '+string(0.,FORMAT='(F0.3)')
				endelse
				printf, lun, 'END'
			endfor
			free_lun, lun

		endif
	endfor

	command='scamp @'+scamp_list_file+' -c scamp_config/ctio_decam.scamp'+' -MERGEDOUTCAT_TYPE '+scamp_cat_type_out+' -MERGEDOUTCAT_NAME '+scamp_cat_file_out+' -MATCH Y -WRITE_XML Y -XML_NAME '+scamp_xml_file+' -SAVE_REFCATALOG Y -REFOUT_CATPATH '+scamp_refcat_dir+' -CHECKPLOT_DEV PSC -CHECKPLOT_ANTIALIAS Y -CHECKPLOT_TYPE '+scamp_check_type + ' -CHECKPLOT_NAME '+scamp_check_file + ' -ASTREF_CATALOG '+scamp_astref_catalog+' -ASTREF_BAND '+scamp_astref_band+' -ASTREFMAG_LIMITS '+scamp_astrefmag_limits+' -DISTORT_DEGREES '+scamp_distort_degrees+' -PHOTCLIP_NSIGMA 2. -SOLVE_ASTROM Y -SOLVE_PHOTOM Y -POSITION_MAXERR '+scamp_pos_error+' -PIXSCALE_MAXERR '+scamp_scale_error+' -POSANGLE_MAXERR '+scamp_angle_error+' -SN_THRESHOLDS '+scamp_sn_thresholds+' -FWHM_THRESHOLDS '+scamp_fwhm_thresholds+' -CROSSID_RADIUS '+scamp_crossid_radius+' -MATCH_RESOL '+scamp_match_resol

	print, command
	spawn, command

endif else $
if recipe EQ 'swarp' then begin

	swarp_combine_type='CLIPPED'
	swarp_resampling_type='LANCZOS2'
	swarp_weight_suffix='.WEIGHT.fits'
	swarp_verbose='NORMAL'
	swarp_resample_do='Y'
	swarp_resample_dir='/Volumes/Data1/temp/resample'

	if not file_test(swarp_resample_dir, /directory) then file_mkdir, swarp_resample_dir, /noexpand_path

	tile_uniq=input_target[uniq(input_target.tile, sort(input_target.tile))].tile
	filter_uniq=input_target[uniq(input_target.filter, sort(input_target.filter))].filter

	for i=0L, n_elements(tile_uniq)-1 do begin
		for j=0L, n_elements(filter_uniq)-1 do begin
			gv=where(input_target.tile EQ tile_uniq[i] AND input_target.filter EQ filter_uniq[j], n_gv)
			if do_dither then input_dither_full=string(input_target[gv].dither) $
			else input_dither_full='ALL'

	     	if do_align_fwhm NE '' then begin
		
				swarp_alignref_h=list()	
				for k=0L, n_elements(input_align_fwhm_full)-1 do begin

					swarp_list_file=output_stack_swarp_dir+'/swarp_header_tile'+tile_uniq[i]+'_'+filter_uniq[j]+'_fwhm'+strjoin(string(input_align_fwhm_full[k],FORMAT='(F0.1)'),'-')+'.lst'
					swarp_im_out=output_stack_swarp_dir+'/swarp_header_tile'+tile_uniq[i]+'_'+filter_uniq[j]+'_fwhm'+strjoin(string(input_align_fwhm_full[k],FORMAT='(F0.1)'),'-')+'.fits'
					swarp_weight_out=output_stack_swarp_dir+'/swarp_header_tile'+tile_uniq[i]+'_'+filter_uniq[j]+'_fwhm'+strjoin(string(input_align_fwhm_full[k],FORMAT='(F0.1)'),'-')+'.WEIGHT.fits'
					swarp_xml_file=output_stack_swarp_dir+'/swarp_header_tile'+tile_uniq[i]+'_'+filter_uniq[j]+'_fwhm'+strjoin(string(input_align_fwhm_full[k],FORMAT='(F0.1)'),'-')+'.xml'

					gv=where(input_target.tile EQ tile_uniq[i] AND input_target.filter EQ filter_uniq[j] AND input_target.fwhm GE (input_align_fwhm_full[k])[0]/0.2637 AND input_target.fwhm LT (input_align_fwhm_full[k])[1]/0.2637, n_gv)

			      	if n_gv GT 0 AND ( file_test(swarp_im_out) EQ 0 OR do_overwrite ) then begin

			        	for l=0L, n_gv-1 do begin
			          		file_copy, input_target[gv[l]].scamp_head_file, input_target[gv[l]].swarp_head_file, /OVERWRITE
			          	endfor

			          	print, 'Processing the following image files'
			          	print, input_target[gv].swarp_im_file
			          	print

         				forprint, input_target[gv].swarp_im_file, textout=swarp_list_file, FORMAT='(A)', /NOCOMMENT
						command='swarp @'+swarp_list_file+' -c swarp_config/ctio_decam.swarp'+' -IMAGEOUT_NAME "'+swarp_im_out+'" -WEIGHTOUT_NAME "'+swarp_weight_out+'" -COMBINE_TYPE '+swarp_combine_type+' -RESAMPLE '+swarp_resample_do+' -RESAMPLE_DIR '+swarp_resample_dir+' -SATLEV_DEFAULT 35000 -WEIGHT_TYPE MAP_WEIGHT -WEIGHT_SUFFIX '+swarp_weight_suffix+' -WEIGHT_THRESH 0. -RESCALE_WEIGHTS N -BLANK_BADPIXELS Y -WRITE_XML Y -XML_NAME '+swarp_xml_file+' -VERBOSE_TYPE ' +swarp_verbose+' -RESAMPLING_TYPE '+swarp_resampling_type+' -SUBTRACT_BACK Y -BACK_SIZE 384 -HEADER_ONLY Y'
			          	print, command
			          	spawn, command

		          	swarp_alignref_h.add, headfits(swarp_im_out)
		       		endif else $
		      		if n_gv GT 0 then begin
		      	    	print, 'SWARP - Reading header from align filter image'
		        	 	swarp_alignref_h.add, headfits(swarp_im_out)
		        	endif else $
		         		print, 'SWARP - There are no images to create align filter image'

				endfor
 			endif

			for k=0L, n_elements(input_fwhm_full)-1 do begin

				print
				print, 'SWARP - Filtering FWHM of ', input_fwhm_full[k]
				swarp_resample_suffix='.resamp_fwhm'+( total(input_fwhm_full[k] EQ 99.) GT 0 ? 'ALL':strjoin(string(input_fwhm_full[k],FORMAT='(F0.1)'),'-'))+'.fits'

				for l=0L, n_elements(input_dither_full)-1 do begin
					print
					print, 'SWARP - Filtering dither number of ', input_dither_full[l]
	
					swarp_list_file=output_stack_swarp_dir+'/swarp_decam_tile'+tile_uniq[i]+'_'+filter_uniq[j]+'_fwhm'+( total(input_fwhm_full[k] EQ 99.) GT 0 ? 'ALL':strjoin(string(input_fwhm_full[k],FORMAT='(F0.1)'),'-'))+'_dither'+input_dither_full[l]+'_exclude'+do_dither_exclude+'_program'+do_program_plain+'.dat'
		
					swarp_im_out=output_stack_swarp_dir+'/ms_tile'+tile_uniq[i]+'_'+filter_uniq[j]+'_fwhm'+( total(input_fwhm_full[k] EQ 99.) GT 0 ? 'ALL':strjoin(string(input_fwhm_full[k],FORMAT='(F0.1)'),'-'))+'_dither'+input_dither_full[l]+'_exclude'+do_dither_exclude+'_program'+do_program_plain+'.fits'
					swarp_weight_out=output_stack_swarp_dir+'/ms_tile'+tile_uniq[i]+'_'+filter_uniq[j]+'_fwhm'+( total(input_fwhm_full[k] EQ 99.) GT 0  ? 'ALL':strjoin(string(input_fwhm_full[k],FORMAT='(F0.1)'),'-'))+'_dither'+input_dither_full[l]+'_exclude'+do_dither_exclude+'_program'+do_program_plain+'.WEIGHT.fits'
					swarp_xml_file=output_stack_swarp_dir+'/swarp_decam_'+tile_uniq[i]+'_'+filter_uniq[j]+'_fwhm'+( total(input_fwhm_full[k] EQ 99.) GT 0 ? 'ALL':strjoin(string(input_fwhm_full[k],FORMAT='(F0.1)'),'-'))+'_dither'+input_dither_full[l]+'_exclude'+do_dither_exclude+'_program'+do_program_plain+'.xml'
			
					sex_stack_im_file=swarp_im_out
					sex_stack_weight_file=swarp_weight_out
					sex_stack_cat_file=output_stack_sex_dir+'/ms_tile'+tile_uniq[i]+'_'+filter_uniq[j]+'_fwhm'+( total(input_fwhm_full[k] EQ 99.) GT 0  ? 'ALL':strjoin(string(input_fwhm_full[k],FORMAT='(F0.1)'),'-'))+'_dither'+input_dither_full[l]+'_exclude'+do_dither_exclude+'_program'+do_program_plain+'.ldac'
					sex_stack_xml_file=output_stack_sex_dir+'/ms_tile'+tile_uniq[i]+'_'+filter_uniq[j]+'_fwhm'+( total(input_fwhm_full[k] EQ 99.) GT 0  ? 'ALL':strjoin(string(input_fwhm_full[k],FORMAT='(F0.1)'),'-'))+'_dither'+input_dither_full[l]+'_exclude'+do_dither_exclude+'_program'+do_program_plain+'.xml'
					sex_stack_checkimage_type='BACKGROUND,SEGMENTATION'
					sex_stack_checkimage_file=strjoin(output_stack_check_dir+'/ms_tile'+tile_uniq[i]+'_'+filter_uniq[j]+'_fwhm'+( total(input_fwhm_full[k] EQ 99.) GT 0 ? 'ALL':strjoin(string(input_fwhm_full[k],FORMAT='(F0.1)'),'-'))+'_dither'+input_dither_full[l]+'_exclude'+do_dither_exclude+'_program'+do_program_plain+['.CHECK_BACK.fits','.CHECK_SEGMENTATION.fits'],',')

					if input_dither_full[l] EQ 'ALL' then begin
						case do_fwhm_type of
							'limit': gv=where(input_target.tile EQ tile_uniq[i] AND input_target.filter EQ filter_uniq[j] AND input_target.fwhm LE input_fwhm_full[k]/0.2637, n_gv)
							'range': gv=where(input_target.tile EQ tile_uniq[i] AND input_target.filter EQ filter_uniq[j] AND input_target.fwhm GE (input_fwhm_full[k])[0]/0.2637 AND input_target.fwhm LT (input_fwhm_full[k])[1]/0.2637, n_gv)
							else: stop
						endcase
					endif else begin
						case do_fwhm_type of
							'limit': gv=where(input_target.tile EQ tile_uniq[i] AND input_target.filter EQ filter_uniq[j] AND input_target.fwhm LE input_fwhm_full[k]/0.2637 AND input_target.dither EQ input_dither_full[l], n_gv)
							'range': gv=where(input_target.tile EQ tile_uniq[i] AND input_target.filter EQ filter_uniq[j] AND input_target.fwhm GE (input_fwhm_full[k])[0]/0.2637 AND input_target.fwhm LT (input_fwhm_full[k])[1]/0.2637 AND input_target.dither EQ input_dither_full[l], n_gv)
							else: stop
						endcase
					endelse

					if n_gv GT 0 then sex_stack_fwhm=mean(input_target[gv].fwhm) $
					else stop
	
					if n_gv GT 0 then begin
	
						for m=0L, n_gv-1 do begin
							file_copy, input_target[gv[m]].scamp_head_file, input_target[gv[m]].swarp_head_file, /OVERWRITE
						endfor
	
						print, 'Processing the following image files'
			          	print, input_target[gv].swarp_im_file
			          	print
		
						forprint, input_target[gv].swarp_im_file, textout=swarp_list_file, FORMAT='(A)', /NOCOMMENT
						if file_test(swarp_im_out) EQ 0 OR do_overwrite then begin
							command='swarp @'+swarp_list_file+' -c swarp_config/ctio_decam.swarp'+' -IMAGEOUT_NAME "'+swarp_im_out+'" -WEIGHTOUT_NAME "'+swarp_weight_out+'" -COMBINE_TYPE '+swarp_combine_type+' -RESAMPLE '+swarp_resample_do+' -RESAMPLE_DIR '+swarp_resample_dir+' -RESAMPLE_SUFFIX '+swarp_resample_suffix+' -SATLEV_DEFAULT 35000 -WEIGHT_TYPE MAP_WEIGHT -WEIGHT_SUFFIX '+swarp_weight_suffix+' -WEIGHT_THRESH 0. -RESCALE_WEIGHTS N -BLANK_BADPIXELS Y -WRITE_XML Y -XML_NAME '+swarp_xml_file+' -VERBOSE_TYPE ' +swarp_verbose+' -RESAMPLING_TYPE '+swarp_resampling_type+' -SUBTRACT_BACK Y -BACK_SIZE 384'
							print, command
							spawn, command
						endif

						print, 'FILE: '	, sex_stack_cat_file	
						if file_test(sex_stack_cat_file) EQ 0 OR do_overwrite then begin
							;command='sex "'+sex_stack_im_file+'" -c sex_config/ctio_decam_stack.sex -CATALOG_NAME "'+sex_stack_cat_file+'" -WEIGHT_IMAGE "'+sex_stack_weight_file+'" -MAG_ZEROPOINT 30. -XML_NAME '+sex_stack_xml_file+' -CHECKIMAGE_TYPE '+sex_stack_checkimage_type+' -CHECKIMAGE_NAME '+sex_stack_checkimage_file+' -BACK_SIZE 384'
							command='sex '+sex_stack_im_file+' -c sex_config/ctio_decam_stack.sex -PARAMETERS_NAME sex_config/ctio_decam_stack.param -CATALOG_NAME '+sex_stack_cat_file+' -WEIGHT_IMAGE '+sex_stack_weight_file+' -MAG_ZEROPOINT 30. -XML_NAME '+sex_stack_xml_file+' -CHECKIMAGE_TYPE '+sex_stack_checkimage_type+' -CHECKIMAGE_NAME '+sex_stack_checkimage_file + ' -PHOT_APERTURES '+strjoin(string(2*[1,2,3,4,5]*sex_stack_fwhm, FORMAT='(F0.2)'),',')
							print, command
							spawn, command
						endif
	
					endif
	
	      	if do_align_fwhm NE '' AND n_gv GT 0 then begin
						for m=0L, n_elements(input_align_fwhm_full)-1 do begin

							if total(input_fwhm_full[k] EQ input_align_fwhm_full[m]) LT 2 then begin

								swarp_align_list_file=output_stack_swarp_dir+'/ms_tile'+tile_uniq[i]+'_'+filter_uniq[j]+'_fwhm'+( total(input_fwhm_full[k] EQ 99.) GT 0 ? 'ALL':strjoin(string(input_fwhm_full[k],FORMAT='(F0.1)'),'-'))+'_dither'+input_dither_full[l]+'_exclude'+do_dither_exclude+'_program'+do_program_plain+'_ALIGNFWHM'+strjoin(string(input_align_fwhm_full[m],FORMAT='(F0.1)'),'-')+'.lst'
								swarp_align_im_out=output_stack_swarp_dir+'/ms_tile'+tile_uniq[i]+'_'+filter_uniq[j]+'_fwhm'+( total(input_fwhm_full[k] EQ 99.) GT 0 ? 'ALL':strjoin(string(input_fwhm_full[k],FORMAT='(F0.1)'),'-'))+'_dither'+input_dither_full[l]+'_exclude'+do_dither_exclude+'_program'+do_program_plain+'_ALIGNFWHM'+strjoin(string(input_align_fwhm_full[m],FORMAT='(F0.1)'),'-')+'.fits'
								swarp_align_weight_out=output_stack_swarp_dir+'/ms_tile'+tile_uniq[i]+'_'+filter_uniq[j]+'_fwhm'+( total(input_fwhm_full[k] EQ 99.) GT 0  ? 'ALL':strjoin(string(input_fwhm_full[k],FORMAT='(F0.1)'),'-'))+'_dither'+input_dither_full[l]+'_exclude'+do_dither_exclude+'_program'+do_program_plain+'_ALIGNFWHM'+strjoin(string(input_align_fwhm_full[m],FORMAT='(F0.1)'),'-')+'.WEIGHT.fits'
								swarp_align_xml_file=output_stack_swarp_dir+'/ms_tile_'+tile_uniq[i]+'_'+filter_uniq[j]+'_fwhm'+( total(input_fwhm_full[k] EQ 99.) GT 0 ? 'ALL':strjoin(string(input_fwhm_full[k],FORMAT='(F0.1)'),'-'))+'_dither'+input_dither_full[l]+'_exclude'+do_dither_exclude+'_program'+do_program_plain+'_ALIGNFWHM'+strjoin(string(input_align_fwhm_full[m],FORMAT='(F0.1)'),'-')+'_swarp.xml'
								swarp_align_im_head=output_stack_swarp_dir+'/ms_tile'+tile_uniq[i]+'_'+filter_uniq[j]+'_fwhm'+( total(input_fwhm_full[k] EQ 99.) GT 0 ? 'ALL':strjoin(string(input_fwhm_full[k],FORMAT='(F0.1)'),'-'))+'_dither'+input_dither_full[l]+'_exclude'+do_dither_exclude+'_program'+do_program_plain+'_ALIGNFWHM'+strjoin(string(input_align_fwhm_full[m],FORMAT='(F0.1)'),'-')+'.head'

								sex_alignref_stack_im_file=output_stack_swarp_dir+'/ms_tile'+tile_uniq[i]+'_'+filter_uniq[j]+'_fwhm'+( total(input_align_fwhm_full[m] EQ 99.) GT 0 ? 'ALL':strjoin(string(input_align_fwhm_full[m],FORMAT='(F0.1)'),'-'))+'_dither'+input_dither_full[l]+'_exclude'+do_dither_exclude+'_program'+do_program_plain+'.fits'
								sex_alignref_stack_weight_file=output_stack_swarp_dir+'/ms_tile'+tile_uniq[i]+'_'+filter_uniq[j]+'_fwhm'+( total(input_align_fwhm_full[m] EQ 99.) GT 0  ? 'ALL':strjoin(string(input_align_fwhm_full[m],FORMAT='(F0.1)'),'-'))+'_dither'+input_dither_full[l]+'_exclude'+do_dither_exclude+'_program'+do_program_plain+'.WEIGHT.fits'

								sex_align_stack_im_file=swarp_align_im_out
								sex_align_stack_weight_file=swarp_align_weight_out
								sex_align_stack_cat_file=repstr(repstr(swarp_align_im_out,output_stack_swarp_dir,output_stack_sex_dir), '.fits', '.ldac')
								sex_align_stack_xml_file=repstr(repstr(swarp_align_im_out,output_stack_swarp_dir,output_stack_sex_dir), '.fits', '_sex.xml')
								sex_align_stack_checkimage_type=sex_stack_checkimage_type
								sex_align_stack_checkimage_file=sex_stack_checkimage_file

								if file_test(swarp_align_im_out) EQ 0 OR do_overwrite then begin
	
									for n=0L, n_gv-1 do begin
										file_copy, input_target[gv[n]].scamp_head_file, input_target[gv[n]].swarp_head_file, /OVERWRITE
									endfor
	
									openw, lun, swarp_align_im_head, /get_lun
            			printf, lun, swarp_alignref_h[m], FORMAT='(A)'
            			close, lun, /all
	
									forprint, input_target[gv].swarp_im_file, textout=swarp_align_list_file, FORMAT='(A)', /NOCOMMENT
		
									command='swarp @'+swarp_align_list_file+' -c swarp_config/ctio_decam.swarp'+' -IMAGEOUT_NAME '+swarp_align_im_out+' -WEIGHTOUT_NAME '+swarp_align_weight_out+' -COMBINE_TYPE '+swarp_combine_type+' -RESAMPLE '+swarp_resample_do+' -RESAMPLE_DIR '+swarp_resample_dir+' -SATLEV_DEFAULT 35000 -WEIGHT_TYPE MAP_WEIGHT -WEIGHT_SUFFIX '+swarp_weight_suffix+' -WEIGHT_THRESH 0. -RESCALE_WEIGHTS N -BLANK_BADPIXELS Y -WRITE_XML Y -XML_NAME '+swarp_align_xml_file+' -VERBOSE_TYPE ' +swarp_verbose+' -RESAMPLING_TYPE '+swarp_resampling_type+' -SUBTRACT_BACK Y -BACK_SIZE 384'
									print, command
									spawn, command

								endif

								if file_test(sex_align_stack_cat_file) EQ 0 OR do_overwrite then begin
									command='sex '+sex_alignref_stack_im_file+','+sex_align_stack_im_file+' -c sex_config/ctio_decam_stack.sex -PARAMETERS_NAME sex_config/ctio_decam_stack.param -CATALOG_NAME '+sex_align_stack_cat_file+' -WEIGHT_IMAGE '+sex_alignref_stack_weight_file+','+sex_align_stack_weight_file+' -MAG_ZEROPOINT 30. -XML_NAME '+sex_align_stack_xml_file+' -CHECKIMAGE_TYPE '+sex_align_stack_checkimage_type+' -CHECKIMAGE_NAME '+sex_align_stack_checkimage_file + ' -PHOT_APERTURES '+strjoin(string(2*[1,2,3,4,5]*sex_stack_fwhm, FORMAT='(F0.2)'),',')
									print, command
									spawn, command
								endif

							endif
							print
		
						endfor				

					endif

				endfor
			endfor
		endfor
	endfor

endif else $
if recipe EQ 'psfex' then begin

	vig_diam=101
	sex_mag_bin=0.5
	sex_radius_bin=0.2
	sex_flux_radius_min=1.5
	sex_ellipticity_min=0.8
	sex_flag_max=1

	tile_uniq=input_target[uniq(input_target.tile, sort(input_target.tile))].tile
	filter_uniq=input_target[uniq(input_target.filter, sort(input_target.filter))].filter

	for i=0L, n_elements(tile_uniq)-1 do begin
		for j=0L, n_elements(filter_uniq)-1 do begin
			gv=where(input_target.tile EQ tile_uniq[i] AND input_target.filter EQ filter_uniq[j], n_gv)
			if do_dither then input_dither_full=string(input_target[gv].dither) $
			else input_dither_full='ALL'

			for k=0L, n_elements(input_fwhm_full)-1 do begin
				print
				print, 'PSFEX - Filtering FWHM of ', input_fwhm_full[k]

				for l=0L, n_elements(input_dither_full)-1 do begin
					print
					print, 'PSFEX - Filtering dither number of ', input_dither_full[l]

					case filter_uniq[j] of
						'u': begin
							sex_mag_range=[10.,18.]
							plot_mag_range=[21,10]
							psfex_mag_range=[15.,18.]+2.5*alog10(600.)
							plot_mag_range=[21,10]
						end
						'g': begin
							sex_mag_range=[15.,21.]
							plot_mag_range=[25,15]
							psfex_mag_range=[16.7,19.5]
							plot_mag_range=[25,15]
						end
						'i': begin
							sex_mag_range=[15.,21.]
							plot_mag_range=[25,15]
							psfex_mag_range=[15.,20.]
							plot_mag_range=[25,15]
						end
						'z': begin
							sex_mag_range=[15.,21.]
							plot_mag_range=[25,15]
							psfex_mag_range=[15.,20.]
						end
						else: stop
					endcase

					sex_stack_im_file=output_stack_swarp_dir+'/ms_tile'+tile_uniq[i]+'_'+filter_uniq[j]+'_fwhm'+( total(input_fwhm_full[k] EQ 99.) GT 0 ? 'ALL':strjoin(string(input_fwhm_full[k],FORMAT='(F0.1)'),'-'))+'_dither'+input_dither_full[l]+'_exclude'+do_dither_exclude+'_program'+do_program_plain+'.fits'
					sex_stack_weight_file=output_stack_swarp_dir+'/ms_tile'+tile_uniq[i]+'_'+filter_uniq[j]+'_fwhm'+( total(input_fwhm_full[k] EQ 99.) GT 0 ? 'ALL':strjoin(string(input_fwhm_full[k],FORMAT='(F0.1)'),'-'))+'_dither'+input_dither_full[l]+'_exclude'+do_dither_exclude+'_program'+do_program_plain+'.WEIGHT.fits'
					sex_stack_cat_psfex_file=output_stack_sex_dir+'/ms_tile'+tile_uniq[i]+'_'+filter_uniq[j]+'_fwhm'+( total(input_fwhm_full[k] EQ 99.) GT 0 ? 'ALL':strjoin(string(input_fwhm_full[k],FORMAT='(F0.1)'),'-'))+'_dither'+input_dither_full[l]+'_exclude'+do_dither_exclude+'_program'+do_program_plain+'_psfex.ldac'
					sex_stack_psfex_xml_file=output_stack_sex_dir+'/ms_tile'+tile_uniq[i]+'_'+filter_uniq[j]+'_fwhm'+( total(input_fwhm_full[k] EQ 99.) GT 0 ? 'ALL':strjoin(string(input_fwhm_full[k],FORMAT='(F0.1)'),'-'))+'_dither'+input_dither_full[l]+'_exclude'+do_dither_exclude+'_program'+do_program_plain+'_psfex.xml'

					sex_stack_cat_file=output_stack_sex_dir+'/ms_tile'+tile_uniq[i]+'_'+filter_uniq[j]+'_fwhm'+( total(input_fwhm_full[k] EQ 99.) GT 0 ? 'ALL':strjoin(string(input_fwhm_full[k],FORMAT='(F0.1)'),'-'))+'_dither'+input_dither_full[l]+'_exclude'+do_dither_exclude+'_program'+do_program_plain+'_psf.ldac'
					sex_stack_xml_file=output_stack_sex_dir+'/ms_tile'+tile_uniq[i]+'_'+filter_uniq[j]+'_fwhm'+( total(input_fwhm_full[k] EQ 99.) GT 0 ? 'ALL':strjoin(string(input_fwhm_full[k],FORMAT='(F0.1)'),'-'))+'_dither'+input_dither_full[l]+'_exclude'+do_dither_exclude+'_program'+do_program_plain+'_psf.xml'

					;gv=where(input_target.tile EQ tile_uniq[i] AND input_target.filter EQ filter_uniq[j] AND input_target.fwhm GE (input_fwhm_full[k])[0]/0.2637 AND input_target.fwhm LT (input_fwhm_full[k])[1]/0.2637, n_gv)

					if input_dither_full[l] EQ 'ALL' then begin
						case do_fwhm_type of
							'limit': gv=where(input_target.tile EQ tile_uniq[i] AND input_target.filter EQ filter_uniq[j] AND input_target.fwhm LE input_fwhm_full[k]/0.2637, n_gv)
							'range': gv=where(input_target.tile EQ tile_uniq[i] AND input_target.filter EQ filter_uniq[j] AND input_target.fwhm GE (input_fwhm_full[k])[0]/0.2637 AND input_target.fwhm LT (input_fwhm_full[k])[1]/0.2637, n_gv)
							else: stop
						endcase
					endif else begin
						case do_fwhm_type of
							'limit': gv=where(input_target.tile EQ tile_uniq[i] AND input_target.filter EQ filter_uniq[j] AND input_target.fwhm LE input_fwhm_full[k]/0.2637 AND input_target.dither EQ input_dither_full[l], n_gv)
							'range': gv=where(input_target.tile EQ tile_uniq[i] AND input_target.filter EQ filter_uniq[j] AND input_target.fwhm GE (input_fwhm_full[k])[0]/0.2637 AND input_target.fwhm LT (input_fwhm_full[k])[1]/0.2637 AND input_target.dither EQ input_dither_full[l], n_gv)
							else: stop
						endcase
					endelse
				
					if n_gv GT 0 then sex_stack_fwhm=mean(input_target[gv].fwhm) $
					else stop

					sex_stack_im_low_file=output_stack_swarp_dir+'/ms_tile'+tile_uniq[i]+'_'+filter_uniq[j]+'_fwhm'+strjoin(string(input_fwhm_full[0],FORMAT='(F0.1)'),'-')+'_dither'+input_dither_full[l]+'_exclude'+do_dither_exclude+'_program'+do_program_plain+'.fits'
					sex_stack_weight_low_file=output_stack_swarp_dir+'/ms_tile'+tile_uniq[i]+'_'+filter_uniq[j]+'_fwhm'+strjoin(string(input_fwhm_full[k],FORMAT='(F0.1)'),'-')+'_dither'+input_dither_full[l]+'_exclude'+do_dither_exclude+'_program'+do_program_plain+'.WEIGHT.fits'
					sex_stack_cat_match_low_file=output_stack_sex_dir+'/ms_tile'+tile_uniq[i]+'_'+filter_uniq[j]+'_fwhm'+strjoin(string(input_fwhm_full[k],FORMAT='(F0.1)'),'-')+'_dither'+input_dither_full[l]+'_exclude'+do_dither_exclude+'_program'+do_program_plain+'_match_low_fwhm_'+strjoin(string(input_fwhm_full[0],FORMAT='(F0.1)'),'-')+'.ldac'
					sex_stack_xml_match_low_file=output_stack_sex_dir+'/ms_tile'+tile_uniq[i]+'_'+filter_uniq[j]+'_fwhm'+strjoin(string(input_fwhm_full[k],FORMAT='(F0.1)'),'-')+'_dither'+input_dither_full[l]+'_exclude'+do_dither_exclude+'_program'+do_program_plain+'_match_low_fwhm_'+strjoin(string(input_fwhm_full[0],FORMAT='(F0.1)'),'-')+'.xml'

					sex_stack_im_high_file=output_stack_swarp_dir+'/ms_tile'+tile_uniq[i]+'_'+filter_uniq[j]+'_fwhm'+strjoin(string(input_fwhm_full[-1],FORMAT='(F0.1)'),'-')+'_dither'+input_dither_full[l]+'_exclude'+do_dither_exclude+'_program'+do_program_plain+'.fits'
					sex_stack_weight_high_file=output_stack_swarp_dir+'/ms_tile'+tile_uniq[i]+'_'+filter_uniq[j]+'_fwhm'+strjoin(string(input_fwhm_full[-1],FORMAT='(F0.1)'),'-')+'_dither'+input_dither_full[l]+'_exclude'+do_dither_exclude+'_program'+do_program_plain+'.WEIGHT.fits'
					sex_stack_cat_match_high_file=output_stack_sex_dir+'/ms_tile'+tile_uniq[i]+'_'+filter_uniq[j]+'_fwhm'+strjoin(string(input_fwhm_full[k],FORMAT='(F0.1)'),'-')+'_dither'+input_dither_full[l]+'_exclude'+do_dither_exclude+'_program'+do_program_plain+'_match_high_fwhm_'+strjoin(string(input_fwhm_full[-1],FORMAT='(F0.1)'),'-')+'.ldac'
					sex_stack_xml_match_high_file=output_stack_sex_dir+'/ms_tile'+tile_uniq[i]+'_'+filter_uniq[j]+'_fwhm'+strjoin(string(input_fwhm_full[k],FORMAT='(F0.1)'),'-')+'_dither'+input_dither_full[l]+'_exclude'+do_dither_exclude+'_program'+do_program_plain+'_match_high_fwhm_'+strjoin(string(input_fwhm_full[-1],FORMAT='(F0.1)'),'-')+'.xml'

	;				sex_stack_cat_psf_file=output_stack_sex_dir+'/ms_tile'+tile_uniq[i]+'_'+filter_uniq[j]+'_fwhm'+(input_fwhm_full[k] EQ 99. ? 'ALL':string(input_fwhm_full[k],FORMAT='(F0.1)'))+'_psfex.ldac'
	;				sex_stack_psf_xml_file=output_stack_sex_dir+'/ms_tile'+tile_uniq[i]+'_'+filter_uniq[j]+'_fwhm'+(input_fwhm_full[k] EQ 99. ? 'ALL':string(input_fwhm_full[k],FORMAT='(F0.1)'))+'_psfex.xml'

					sex_stack_checkimage_type='BACKGROUND,SEGMENTATION'
					sex_stack_checkimage_file=strjoin(output_stack_check_dir+'/ms_tile'+tile_uniq[i]+'_'+filter_uniq[j]+'_fwhm'+( total(input_fwhm_full[k] EQ 99.) GT 0 ? 'ALL': strjoin(string(input_fwhm_full[k],FORMAT='(F0.1)'),'-'))+'_dither'+input_dither_full[l]+'_exclude'+do_dither_exclude+'_program'+do_program_plain+['.CHECK_BACK.fits','.CHECK_SEGMENTATION.fits'],',')
				
					psfex_stack_cat_file=output_stack_sex_dir+'/ms_tile'+tile_uniq[i]+'_'+filter_uniq[j]+'_fwhm'+( total(input_fwhm_full[k] EQ 99.) GT 0 ? 'ALL':strjoin(string(input_fwhm_full[k],FORMAT='(F0.1)'),'-'))+'_dither'+input_dither_full[l]+'_exclude'+do_dither_exclude+'_program'+do_program_plain+'_sex_psfex_stars.ldac'
					psfex_stack_xml_file=output_stack_sex_dir+'/ms_tile'+tile_uniq[i]+'_'+filter_uniq[j]+'_fwhm'+( total(input_fwhm_full[k] EQ 99.) GT 0 ? 'ALL':strjoin(string(input_fwhm_full[k],FORMAT='(F0.1)'),'-'))+'_dither'+input_dither_full[l]+'_exclude'+do_dither_exclude+'_program'+do_program_plain+'_sex_psfex_stars.xml'
					psfex_stack_psf_file=repstr(psfex_stack_cat_file, '.ldac', '.psf')
					psfex_check_type = 'RESIDUALS,SNAPSHOTS,SAMPLES'
					psfex_check_file = strjoin(output_stack_psfex_dir+'/psfex_checkimage'+['.CHECK_RESIDUALS.fits','.CHECK_SNAPSHOTS.fits','.CHECK_SAMPLES.fits'],',')
					psfex_checkplot_type = 'FWHM,ELLIPTICITY,COUNTS'
					psfex_checkplot_file = strjoin(output_stack_psfex_dir+'/psfex_checkplot'+['.CHECK_FWHM','.CHECK_ELLIPTICITY','.CHECK_COUNTS'],',')
					psfex_psfvar_nsnap='14'
					psfex_psfvar_degrees='3'
					psfex_basis_type='PIXEL'
					psfex_basis_number='20'
					psfex_psf_size='47,47'
					psfex_psf_sampling='0.8'
				
				;	cat_stars_data=mrdfits('catalogs/psfex_stars_v2.fits', 1, cat_stars_h, COLUMNS=['RA','DEC','MAG_AUTO','FLAGS'])
				;	gv_stars=where(cat_stars_data.flags[0] LE 3 AND cat_stars_data.flags[1] LE 3 AND cat_stars_data.flags[2] LE 3, n_gv_stars)
	;				print, 'Filename is ', psfex_stack_cat_file
	;				print, file_test(psfex_stack_cat_file)
	;				print, file_test(sex_stack_im_file) AND (file_test(psfex_stack_cat_file) EQ 0 OR do_overwrite) 

					print, sex_stack_im_file, file_test(sex_stack_im_file)
					print, psfex_stack_cat_file, file_test(psfex_stack_cat_file)
					print

					if file_test(sex_stack_im_file) AND (file_test(psfex_stack_cat_file) EQ 0 OR do_overwrite) then begin

						im_h=headfits(sex_stack_im_file)
						im_size=long([fxpar(im_h,'NAXIS1'),fxpar(im_h,'NAXIS2')])
					
						if file_test(sex_stack_cat_psfex_file) EQ 0 OR do_overwrite then begin
							command='sex "'+sex_stack_im_file+'" -c sex_config/ctio_decam_psfex.sex -PARAMETERS_NAME sex_config/ctio_decam_psfex.param -CATALOG_NAME "'+sex_stack_cat_psfex_file+'" -WEIGHT_IMAGE "'+sex_stack_weight_file+'" -MAG_ZEROPOINT 30. -XML_NAME '+sex_stack_psfex_xml_file+' -CHECKIMAGE_TYPE '+sex_stack_checkimage_type+' -CHECKIMAGE_NAME "'+sex_stack_checkimage_file+'"'
							print, command
							spawn, command
						endif

						filter_sex_cat, sex_stack_cat_psfex_file, psfex_stack_cat_file, class='stars'

						command='psfex '+psfex_stack_cat_file+' -c psfex_config/ctio_decam.psfex'+' -PSF_DIR '+output_stack_psfex_dir+' -WRITE_XML Y -XML_NAME '+psfex_stack_xml_file+' -CHECKIMAGE_TYPE NONE '+' -CHECKIMAGE_NAME '+psfex_check_file+' -CHECKPLOT_DEV PNG -CHECKPLOT_TYPE '+psfex_checkplot_type+' -CHECKPLOT_NAME '+psfex_checkplot_file+' -PSFVAR_NSNAP '+psfex_psfvar_nsnap+' -PSFVAR_DEGREES '+psfex_psfvar_degrees+' -BASIS_TYPE '+psfex_basis_type+' -BASIS_NUMBER '+psfex_basis_number+' -SAMPLE_VARIABILITY 1.,1. -SAMPLE_MAXELLIP 1. -NEWBASIS_TYPE NONE -NEWBASIS_NUMBER 10 -SAMPLEVAR_TYPE NONE -STABILITY_TYPE EXPOSURE -SAMPLE_MINSN 1. -SAMPLE_FWHMRANGE 0.1,10. -SAMPLE_AUTOSELECT N -PSF_SIZE '+psfex_psf_size+' -PSF_SAMPLING '+psfex_psf_sampling+' -PSF_ACCURACY 0.01 -BADPIXEL_FILTER Y -BADPIXEL_NMAX 80'
						print, command
						spawn, command
					endif

					if file_test(psfex_stack_psf_file) AND ((file_test(sex_stack_cat_file) EQ 0 OR do_overwrite)) then begin
						command='sex '+sex_stack_im_file+' -c sex_config/ctio_decam_stack.sex -PARAMETERS_NAME sex_config/ctio_decam_stack_psf.param -CATALOG_NAME '+sex_stack_cat_file+' -WEIGHT_IMAGE '+sex_stack_weight_file+' -MAG_ZEROPOINT 30. -XML_NAME '+sex_stack_xml_file+' -CHECKIMAGE_TYPE '+sex_stack_checkimage_type+' -CHECKIMAGE_NAME '+sex_stack_checkimage_file + ' -PSF_NAME '+psfex_stack_psf_file+' -PSF_NMAX 1'
						print, command
						spawn, command
					endif

					case k of
						0: begin
							command='sex '+sex_stack_im_file+' -c sex_config/ctio_decam_stack.sex -PARAMETERS_NAME sex_config/ctio_decam_stack.param -CATALOG_NAME '+sex_stack_cat_match_low_file+' -WEIGHT_IMAGE '+sex_stack_weight_file+' -MAG_ZEROPOINT 30. -XML_NAME '+sex_stack_xml_match_low_file+' -CHECKIMAGE_TYPE '+sex_stack_checkimage_type+' -CHECKIMAGE_NAME '+sex_stack_checkimage_file + ' -PHOT_APERTURES '+strjoin(string(2*[1,2,3,4,5]*sex_stack_fwhm, FORMAT='(F0.2)'),',')
							print, command
							spawn, command
						end
						else: begin
							command='sex '+sex_stack_im_low_file+','+sex_stack_im_file+' -c sex_config/ctio_decam_stack.sex -PARAMETERS_NAME sex_config/ctio_decam_stack.param -CATALOG_NAME '+sex_stack_cat_match_low_file+' -WEIGHT_IMAGE '+sex_stack_weight_low_file+','+sex_stack_weight_file+' -MAG_ZEROPOINT 30. -XML_NAME '+sex_stack_xml_match_low_file+' -CHECKIMAGE_TYPE '+sex_stack_checkimage_type+' -CHECKIMAGE_NAME '+sex_stack_checkimage_file + ' -PHOT_APERTURES '+strjoin(string(2*[1,2,3,4,5]*sex_stack_fwhm, FORMAT='(F0.2)'),',')
							print, command
							spawn, command

							command='sex '+sex_stack_im_high_file+','+sex_stack_im_file+' -c sex_config/ctio_decam_stack.sex -PARAMETERS_NAME sex_config/ctio_decam_stack.param -CATALOG_NAME '+sex_stack_cat_match_high_file+' -WEIGHT_IMAGE '+sex_stack_weight_high_file+','+sex_stack_weight_file+' -MAG_ZEROPOINT 30. -XML_NAME '+sex_stack_xml_match_high_file+' -CHECKIMAGE_TYPE '+sex_stack_checkimage_type+' -CHECKIMAGE_NAME '+sex_stack_checkimage_file + ' -PHOT_APERTURES '+strjoin(string(2*[1,2,3,4,5]*sex_stack_fwhm, FORMAT='(F0.2)'),',')
							print, command
							spawn, command
						end
					endcase

				endfor
			endfor
		endfor
	endfor

endif else $
if recipe EQ 'sex psf' then begin

	sex_stack_im_file=swarp_dir+'/ms_tile'+do_tile+'_'+do_filter+'.fits'
	sex_stack_weight_file=swarp_dir+'/ms_tile'+do_tile+'_'+do_filter+'.WEIGHT.fits'
	sex_stack_cat_file=sextractor_dir+'/ms_tile'+do_tile+'_'+do_filter+'_psf.ldac'
	sex_stack_xml_file=sextractor_dir+'/ms_tile'+do_tile+'_'+do_filter+'_psf.xml'
	sex_stack_checkimage_type='BACKGROUND,SEGMENTATION'
	sex_stack_checkimage_file=strjoin(sextractor_check_dir+'/ms_tile'+do_tile+'_'+do_filter+['.CHECK_BACK.fits','.CHECK_SEGMENTATION.fits'],',')
	psfex_psf_file=psfex_dir+'/ms_tile'+do_tile+'_'+do_filter+'_psfex.psf'

	command='sex '+sex_stack_im_file+' -c sex_config/ctio_decam_stack.sex -PARAMETERS_NAME sex_config/ctio_decam_stack_psf.param -CATALOG_NAME '+sex_stack_cat_file+' -WEIGHT_IMAGE '+sex_stack_weight_file+' -MAG_ZEROPOINT 30. -XML_NAME '+sex_stack_xml_file+' -CHECKIMAGE_TYPE '+sex_stack_checkimage_type+' -CHECKIMAGE_NAME '+sex_stack_checkimage_file + ' -PSF_NAME '+psfex_psf_file+' -PSF_NMAX 1'
	print, command
	spawn, command

endif else $
if recipe EQ 'report' then begin

	tile_uniq=input_target[uniq(input_target.tile, sort(input_target.tile))].tile
	filter_uniq=input_target[uniq(input_target.filter, sort(input_target.filter))].filter
	
	for i=0L, n_elements(tile_uniq)-1 do begin
		for j=0L, n_elements(filter_uniq)-1 do begin
			for k=0L, n_elements(input_fwhm_full)-1 do begin

				case do_fwhm_type of
					'limit': gv=where(input_target.tile EQ tile_uniq[i] AND input_target.filter EQ filter_uniq[j] AND input_target.fwhm LE input_fwhm_full[k]/0.2637, n_gv)
					'range': gv=where(input_target.tile EQ tile_uniq[i] AND input_target.filter EQ filter_uniq[j] AND input_target.fwhm GT (input_fwhm_full[k])[0]/0.2637 AND input_target.fwhm LE (input_fwhm_full[k])[1]/0.2637, n_gv)
					else: stop
				endcase
				
				gv=gv[sort(input_target[gv].fwhm)]
				;print, 'Filename		Tile		Filter		FWHM(arsec)'
				;forprint, input_target[gv].im_file, input_target[gv].tile , input_target[gv].filter, input_target[gv].fwhm*0.2637, input_target[gv].exptime, textout=2, format='A,2X,A,2X,A,2X,F0.2,2X,F0.1'

				exptime_total=total(input_target[gv].exptime)
				exptime_uniq=input_target[gv[uniq(input_target[gv].exptime, sort(input_target[gv].exptime))]].exptime

				case do_fwhm_type of
					'limit': print, tile_uniq[i], filter_uniq[j], exptime_total, n_gv, input_fwhm_full[k], FORMAT='("Tile: ",I2," - Filter: ",A," - Total exptime: ",I5," s - Number of images: ", I3," - FWHM < ", A, " arcsec")'
					'range': print, tile_uniq[i], filter_uniq[j], exptime_total, n_gv, (input_fwhm_full[k])[0], (input_fwhm_full[k])[1], FORMAT='("Tile: ",I2," - Filter: ",A," - Total exptime: ",I5," s - Number of images: ", I3," - FWHM = ", F0.1, "-", F0.1, " arcsec")'
					else: stop
				endcase

				for l=0L, n_elements(exptime_uniq)-1 do begin

					case do_fwhm_type of
						'limit': gv_exptime=where(input_target.tile EQ tile_uniq[i] AND input_target.filter EQ filter_uniq[j] AND input_target.fwhm LE input_fwhm_full[k]/0.2637 AND input_target.exptime EQ exptime_uniq[l], n_gv_exptime)					
						'range': gv_exptime=where(input_target.tile EQ tile_uniq[i] AND input_target.filter EQ filter_uniq[j] AND input_target.fwhm GT (input_fwhm_full[k])[0]/0.2637 AND input_target.fwhm LE (input_fwhm_full[k])[1]/0.2637 AND input_target.exptime EQ exptime_uniq[l], n_gv_exptime)
						else: stop
					endcase

;					gv_exptime=where(input_target.tile EQ tile_uniq[i] AND input_target.filter EQ filter_uniq[j] AND input_target.fwhm LE fwhm_max/0.2637 AND input_target.exptime EQ exptime_uniq[k], n_gv_exptime)
					print, n_gv_exptime, exptime_uniq[l], FORMAT='("- ", I2," images with exposure time of ", I3, " seconds")'
					
				endfor
			endfor
		endfor
		print
	endfor

stop

if do_fwhm EQ 99. then plot_yrange=[0,24] else plot_yrange=[0,16]
gv=where(input_iq.sky LT 25000. AND input_iq.fwhm GT 0. AND input_iq.FWHM*0.2637 LT do_fwhm, n_gv)
input_iq=input_iq[gv]

im_exptime=list()
for i=0L, n_elements(input_iq)-1 do begin
	im_h=headfits(input_iq[i].im_file)
	im_exptime.add, fxpar(im_h, 'EXPTIME')
endfor

im_exptime=im_exptime.toarray(type='float')

forprint, input_iq.im_file, im_exptime, format='A,4X,F0.1', textout=2

openw, lun, 'summary_survey_iq.dat', /get_lun
printf, lun, '# FILE  FILTER  TILE  DITHER  SKY_LEVEL  FWHM  EXPTIME'
for i=0L, n_elements(input_iq)-1 do begin
	printf, lun, input_iq[i].im_file, input_iq[i].filter, input_iq[i].tile, input_iq[i].dither, input_iq[i].sky, input_iq[i].fwhm, im_exptime[i], format='(A,4X,A,4X,A,4X,A,4X,F0.1,4X,F0.1,4X,F0.1)'
endfor
free_lun, lun

tile_uniq=input_iq[uniq(input_iq.tile, sort(input_iq.tile))].tile
filter_uniq=input_iq[uniq(input_iq.filter, sort(input_iq.filter))].filter
exptime_uniq=im_exptime[uniq(im_exptime, sort(im_exptime))]

;for i=0L, n_elements(tile_uniq)-1 do begin
	for j=0L, n_elements(filter_uniq)-1 do begin
				
;		gv=where(input_iq.tile EQ tile_uniq[i] AND input_iq.filter EQ filter_uniq[j], n_gv)
		gv=where(input_iq.filter EQ filter_uniq[j], n_gv)

		if n_gv GT 0 then begin			
			exptime_uniq=im_exptime[gv[uniq(im_exptime[gv], sort(im_exptime[gv]))]]
			print, filter_uniq[j], exptime_uniq

      cgwindow, wxsize=800, wysize=300*n_elements(exptime_uniq);, WMulti=[0,1,n_elements(plot_extname)], WOXMargin=[2,6], WOYMargin=[2,2]
      plot_pos = cgLayout([1,n_elements(exptime_uniq)], OXMargin=[8,6], OYMargin=[2,2], XGap=1, YGap=2)

			for k=0L, n_elements(exptime_uniq)-1 do begin

;				gv=where(input_iq.tile EQ tile_uniq[i] AND input_iq.filter EQ filter_uniq[j] AND im_exptime EQ exptime_uniq[k], n_gv)
				gv=where(input_iq.filter EQ filter_uniq[j] AND im_exptime EQ exptime_uniq[k], n_gv)

				if n_gv GT 0 then begin			
					cghistoplot, fix(input_iq[gv].tile)-0.5, bin=1., position=plot_pos[*,k], /fill, POLYCOLOR='sky blue', MININPUT=0., /addcmd, /noerase, xtickvalues=0.5+indgen(9), xticknames=(k EQ n_elements(exptime_uniq)-1) ? string(1+indgen(9),FORMAT='(I0)') : replicate(' ',9), xticks=8, xrange=[-1,10], yrange=plot_yrange, title='Filter '+filter_uniq[j]+' EXPTIME '+string(exptime_uniq[k],FORMAT='(F0.1)'), charsize=1., yticklen=1., ygridstyle=1
				endif $
				else begin
;					cghistoplot, [0], /addcmd, /noerase, xrange=[-1,10], yrange=plot_yrange
				endelse

			endfor
	    cgcontrol, output='results/magellanic_histogram_filter'+filter_uniq[j]+'_fwhm'+(do_fwhm EQ 99. ? 'ALL':string(do_fwhm,FORMAT='(F0.1)'))+'.pdf'


		endif

	endfor
;endfor

cgdelete, /all

endif

end
