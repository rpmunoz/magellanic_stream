pro decam_pipeline, recipe, PROGRAM=program, FILTER=filter, OVERWRITE=overwrite, TILE=tile, DEBUG=debug, STANDARD=standard

do_program= (n_elements(program) GT 0) ? program : '20130808'
do_filter= (n_elements(filter) GT 0) ? filter : 'g'
do_tile= (n_elements(tile) GT 0) ? tile : '1'
do_standard= (n_elements(standard) GT 0) ? standard : 'SDSS'
do_n_chip= (do_program EQ '20130808' OR do_program EQ '20131115') ? '61' : '60'
do_scamp_ast='Y'
do_scamp_phot='Y'
do_overwrite= (n_elements(overwrite) GT 0) ? keyword_set(overwrite) : 0
do_debug= (n_elements(debug) GT 0) ? keyword_set(debug) : 0
o_filter_orig=do_filter
do_program_orig=do_program

input_dir='/Volumes/Data2/Magellanic_stream_2/DECam'
output_dir='/Volumes/Data2/Magellanic_stream_2/DECam'

if do_program EQ 'all' then begin
  temp_program=strsplit(file_search(input_dir+'/*', /TEST_DIRECTORY),'/',/extract)
  do_program=strarr(n_elements(temp_program))
  for i=0L, n_elements(temp_program)-1 do do_program[i]=(temp_program[i])[-1]
  gv=where(do_program NE 'stacks', n_gv)
  if n_gv GT 0 then do_program=do_program[gv]
endif

input_im_dir=input_dir+'/'+do_program+'/processed'
input_calib_dir=input_dir+'/'+do_program+'/calib'
;input_im_orig_dir='/Volumes/Data2/Magellanic_stream_2/DECam/data/'+do_program+'/processed'
;input_calib_dir='/Volumes/Data2/Magellanic_stream_2/DECam/data/'+do_program+'/calib'
;input_im_dir=output_dir+'/pipeline/images'

input_n_chip=do_n_chip
input_chip_fwhm=[28,35] ; Where to measure FWHM
input_chip_flux_radius=[27,28,29,34,35,36] ; Where to measure FWHM
output_im_dir=output_dir+'/'+do_program+'/pipeline/images'
output_calib_dir=output_dir+'/'+do_program+'/pipeline/calib'
output_sex_dir=output_dir+'/'+do_program+'/pipeline/sextractor'
output_sex_check_dir=output_dir+'/'+do_program+'/pipeline/sextractor/check'
output_scamp_dir=output_dir+'/'+do_program+'/pipeline/scamp'
output_stack_swarp_dir=output_dir+'/stacks';output_dir+'/'+do_program+'/pipeline/swarp'
output_stack_sex_dir=output_dir+'/stacks';output_dir+'/'+do_program+'/pipeline/swarp'
output_stack_check_dir=output_dir+'/stacks/check';output_dir+'/'+do_program+'/pipeline/swarp'
output_psfex_dir=output_dir+'/stacks' ;output_dir+'/'+do_program+'/pipeline/psfex'

survey_info = list( {tile:'1', coo:['23:26:07.46','-53:28:13.7'], filter:['g','i','z']}, $
	{tile:'2', coo:['23:53:46.15','-52:32:21.7'], filter:['g','i','z']}, $
	{tile:'3', coo:['00:22:34.67','-51:17:37.6'], filter:['g','i','z']}, $
	{tile:'4', coo:['00:53:43.00','-49:36:27.2'], filter:['g','i','z']}, $
	{tile:'5', coo:['01:03:22.85','-48:35:37.4'], filter:['g','i','z']}, $
	{tile:'6', coo:['01:12:56.00','-47:35:04.6'], filter:['g','i','z']}, $
	{tile:'7', coo:['01:37:49.01','-45:26:39.8'], filter:['g','i','z']}, $
	{tile:'8', coo:['01:59:33.94','-42:55:45.3'], filter:['g','i']}, $
	{tile:'9', coo:['02:18:48.79','-40:11:55.5'], filter:['g','i']} )

survey_sequence = { seq1:['1','1'], seq2:['2','2'], seq3:['3','3'], seq4:['4','4'], seq5:['5','5'], seq6:['6','6'], seq7=['7','7'], seq8=['8','8'], seq9=['9','9']  }

if file_test('survey_target_'+do_program+'.dat') then begin
	readcol, 'survey_target_'+do_program+'.dat', temp_im_orig_file, temp_filter, temp_tile, temp_dither, temp_weight_orig_file, temp_zp, temp_fwhm, temp_mjd, temp_exptime, FORMAT='A,A,A,A,A,F,F,D,F', COMMENT='#'
	n_survey_info=n_elements(temp_im_orig_file)
	create_struct, input_target, '', ['im_orig_file','filter','tile','dither','weight_orig_file','mjd','zp','im_file','weight_file','exptime','fwhm','sex_cat_file','sex_xml_file','sex_check_file','scamp_cat_file','scamp_head_file','scamp_ahead_file','swarp_im_file','swarp_head_file','sex_zp_file'], 'A,A,A,A,A,D,F,A,A,F,F,A,A,A,A,A,A,A,A,A', dim=n_survey_info
	input_target.im_orig_file=input_im_orig_dir+'/'+repstr(temp_im_orig_file,'.fits.fz','.fits')
	input_target.filter=temp_filter
	input_target.tile=strtrim(temp_tile,2)
	input_target.dither=strtrim(temp_dither,2)
	input_target.weight_orig_file=input_im_orig_dir+'/'+repstr(temp_weight_orig_file,'.fits.fz','.fits')
	input_target.mjd=temp_mjd
	input_target.zp=temp_zp
	input_target.exptime=temp_exptime
	input_target.im_file=input_im_dir+'/'+'fornax_t'+strtrim(temp_tile,2)+'_d'+strtrim(temp_dither,2)+'_'+temp_filter+'.fits'
	input_target.weight_file=input_im_dir+'/'+'fornax_t'+strtrim(temp_tile,2)+'_d'+strtrim(temp_dither,2)+'_'+temp_filter+'.WEIGHT.fits'
	input_target.fwhm=temp_fwhm
	input_target.sex_cat_file=sextractor_dir+'/'+'fornax_t'+strtrim(temp_tile,2)+'_d'+strtrim(temp_dither,2)+'_'+temp_filter+'.ldac'
	input_target.sex_xml_file=sextractor_dir+'/'+'fornax_t'+strtrim(temp_tile,2)+'_d'+strtrim(temp_dither,2)+'_'+temp_filter+'.xml'
	input_target.sex_check_file=sextractor_check_dir+'/'+'fornax_t'+strtrim(temp_tile,2)+'_d'+strtrim(temp_dither,2)+'_'+temp_filter+'.BACKGROUND.fits'
	input_target.sex_zp_file=sextractor_dir+'/'+'fornax_t'+strtrim(temp_tile,2)+'_d'+strtrim(temp_dither,2)+'_'+temp_filter+'_zp.dat'
	input_target.scamp_cat_file=scamp_dir+'/'+'fornax_t'+strtrim(temp_tile,2)+'_d'+strtrim(temp_dither,2)+'_'+temp_filter+'_stars.ldac'
	input_target.scamp_head_file=scamp_dir+'/'+'fornax_t'+strtrim(temp_tile,2)+'_d'+strtrim(temp_dither,2)+'_'+temp_filter+'_stars.head'
	input_target.scamp_ahead_file=scamp_dir+'/'+'fornax_t'+strtrim(temp_tile,2)+'_d'+strtrim(temp_dither,2)+'_'+temp_filter+'_stars.ahead'
	input_target.swarp_im_file=input_im_dir+'/'+'fornax_t'+strtrim(temp_tile,2)+'_d'+strtrim(temp_dither,2)+'_'+temp_filter+'.fits'
	input_target.swarp_head_file=input_im_dir+'/'+'fornax_t'+strtrim(temp_tile,2)+'_d'+strtrim(temp_dither,2)+'_'+temp_filter+'.head'

	case do_tile of
		'all': begin
			case do_filter of
				'all': gv=indgen(n_elements(input_target))
				else: gv=where(input_target.filter EQ do_filter, n_gv)
			endcase
		end
		else: gv=where(input_target.filter EQ do_filter AND input_target.tile EQ do_tile, n_gv)
	endcase
	input_target=input_target[gv]
endif

if file_test('survey_calib_'+do_program+'.dat') then begin
	readcol, 'survey_calib_'+do_program+'.dat', temp_im_orig_file, temp_mjd, temp_object, temp_filter, temp_weight_orig_file, temp_zp, temp_fwhm, temp_photflag, FORMAT='A,A,A,A,A,F,F,A', COMMENT='#'
	n_survey=n_elements(temp_im_orig_file)
	create_struct, input_calib, '', ['im_orig_file','object','filter','tile','mjd','weight_orig_file','zp','airmass','im_file','weight_file','fwhm','sex_cat_file','sex_xml_file','sex_check_file','sex_zp_file','sex_zp_full_file','photflag','scamp_cat_file','scamp_head_file','scamp_ahead_file','swarp_im_file','swarp_head_file'], 'A,A,A,A,A,A,F,F,A,A,F,A,A,A,A,A,A,A,A,A,A,A', dim=n_survey
	input_calib.im_orig_file=input_calib_dir+'/'+repstr(temp_im_orig_file,'.fits.fz','.fits')
	input_calib.object=temp_object
	input_calib.filter=temp_filter
	input_calib.mjd=temp_mjd
	input_calib.weight_orig_file=input_calib_dir+'/'+repstr(temp_weight_orig_file,'.fits.fz','.fits')
	input_calib.zp=temp_zp
	input_calib.airmass=0.
	input_calib.photflag=temp_photflag
	input_calib.im_file=input_im_dir+'/'+'standard_'+temp_object+'_'+temp_filter+'_MJD'+temp_mjd+'.fits'
	input_calib.weight_file=input_im_dir+'/'+'standard_'+temp_object+'_'+temp_filter+'_MJD'+temp_mjd+'.WEIGHT.fits'
	input_calib.fwhm=temp_fwhm
	input_calib.sex_cat_file=sextractor_dir+'/'+'standard_'+temp_object+'_'+temp_filter+'_MJD'+temp_mjd+'.ldac'
	input_calib.sex_xml_file=sextractor_dir+'/'+'standard_'+temp_object+'_'+temp_filter+'_MJD'+temp_mjd+'.xml'
	input_calib.sex_check_file=sextractor_check_dir+'/'+'standard_'+temp_object+'_'+temp_filter+'_MJD'+temp_mjd+'.BACKGROUND.fits'
	input_calib.sex_zp_file=calib_dir+'/'+'standard_'+temp_object+'_'+temp_filter+'_MJD'+temp_mjd+'_zp.dat'
	input_calib.sex_zp_full_file=calib_dir+'/'+'standard_'+temp_object+'_'+temp_filter+'_MJD'+temp_mjd+'_zp_full.dat'
	input_calib.scamp_cat_file=scamp_dir+'/'+'standard_'+temp_object+'_'+temp_filter+'_MJD'+temp_mjd+'_stars.ldac'
	input_calib.scamp_head_file=scamp_dir+'/'+'standard_'+temp_object+'_'+temp_filter+'_MJD'+temp_mjd+'_stars.head'
	input_calib.scamp_ahead_file=scamp_dir+'/'+'standard_'+temp_object+'_'+temp_filter+'_MJD'+temp_mjd+'_stars.ahead'
	input_calib.swarp_im_file=input_im_dir+'/'+'standard_'+temp_object+'_'+temp_filter+'_MJD'+temp_mjd+'.fits'
	input_calib.swarp_head_file=input_im_dir+'/'+'standard_'+temp_object+'_'+temp_filter+'_MJD'+temp_mjd+'.head'

	gv=where(input_calib.filter EQ do_filter, n_gv)
	input_calib=input_calib[gv]
;	case do_tile of
;		'all': gv=where(input_calib.filter EQ do_filter, n_gv)
;		else: gv=where(input_calib.filter EQ do_filter AND input_calib.tile EQ do_tile, n_gv)
;	endcase
;	input_calib=input_calib[gv]
endif

if file_test('survey_zp_'+do_program+'.dat') then begin
	readcol, 'survey_zp_'+do_program+'.dat', temp_mjd, temp_filter, temp_zp, temp_k, temp_zp_error, temp_k_error, temp_photflag, FORMAT='F,A,F,F,F,F,A', COMMENT='#'
	n_survey=n_elements(temp_mjd)
	create_struct, input_zp, '', ['mjd','filter','zp','k','zp_error','k_error','photflag'], 'F,A,F,F,F,F,A', dim=n_survey
	input_zp.mjd=temp_mjd
	input_zp.filter=temp_filter
	input_zp.zp=temp_zp
	input_zp.k=temp_k
	input_zp.zp_error=temp_zp_error
	input_zp.k_error=temp_k_error
	input_zp.photflag=temp_photflag
endif

loadct, 12

if recipe EQ 'database' then begin

	im_file=file_search(input_im_orig_dir, '*.fz')
	n_im_file=n_elements(im_file)

	create_struct, temp_target, '', ['dir','im_file','weight_file','expnum','type','date','ra','dec','mjd','mjd_floor','exptime','tile','dither','filter','zp','fwhm'], 'A,A,A,I,A,A,D,D,D,D,F,A,A,A,F,F', dim=1

	for i=0L, n_im_file-1 do begin
		im_h=headfits(im_file[i])
		temp_target.dir=input_im_orig_dir
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
		temp_target.zp= (temp EQ 4 OR temp EQ 5) ? fxpar(im_h, 'MAGZERO') : 0.
		if i EQ 0 then input_target=temp_target else input_target=[input_target,temp_target]
	endfor
	gv_sort=sort(input_target.date+' '+input_target.type)

	if file_test(input_im_orig_dir+'/pipeline_database.dat') EQ 0 then begin
		openw, lun, input_im_orig_dir+'/pipeline_database.dat', /get_lun
		printf, lun, '# FILE	 DATE-OBS	 EXPNUM	 OBJECT	 FILTER	 EXPTIME	 SEQID	 OBSTYPE	 PROCTYPE		PRODTYPE	 MAGZERO'
		for i=0L, n_im_file-1 do begin
			command='dfits '+im_file[gv_sort[i]]+' | fitsort -d DATE-OBS EXPNUM OBJECT FILTER EXPTIME SEQID OBSTYPE PROCTYPE PRODTYPE MAGZERO'
			spawn, command, result
			result=repstr(result, input_im_orig_dir+'/', '')
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

	for i=0L, n_elements(filter_uniq)-1 do begin
		i_dither=1L
		for j=0L, n_elements(mjd_uniq)-1 do begin
			gv=where( input_target.filter EQ filter_uniq[i] AND input_target.mjd_floor EQ mjd_uniq[j], n_gv)
;			print, filter_uniq[i], mjd_uniq[j], n_gv
			if n_gv GT 0 then begin
				for k=0L, n_tags(survey_sequence)-1 do begin
;					print, filter_uniq[i], mjd_uniq[j], survey_sequence.(k)
					gv_dither=gv[where(input_target[gv].dither EQ ' ', n_gv_dither)]
;					print, n_gv_dither
					if n_gv_dither GT 0 then begin
						temp_n_seq=n_elements(survey_sequence.(k))
						temp_tile=[input_target[gv_dither].tile,strarr(temp_n_seq-1)]
						temp_match=intarr(n_gv_dither)
						for l=0L, n_gv_dither-1 do begin
							temp_data=(shift(temp_tile,-l))[0:temp_n_seq-1]
							temp_match[l]=total(temp_data EQ survey_sequence.(k))
;						print, l, temp_match[l]
;						print, temp_data
						endfor
						gv_match=where(temp_match GE n_elements(survey_sequence.(k))/2, n_gv_match)
						gv_match=[gv_dither[gv_match[where( gv_match - [-n_elements(survey_sequence.(k))/2,gv_match] GE n_elements(survey_sequence.(k))/2, n_gv_match)]],gv_dither[-1]+1]
;						print, n_gv_match
;						print, gv_match
						for l=0L, n_gv_match-1 do begin
							m=gv_match[l]
							temp_tile=''
							while m LT gv_match[l+1] do begin
								gv_tile=where(input_target[m].tile EQ temp_tile, n_gv_tile)
								if n_gv_tile EQ 0 then input_target[m].dither=strn(i_dither) $
								else input_target[m].dither=strn(i_dither)+'b'
								temp_tile=[temp_tile,input_target[m].tile]
								m++
							endwhile
							i_dither++
						endfor
;						forprint, input_target[gv_dither].tile, input_target[gv_dither].dither, text=2, format='A,2X,A'
					endif
				endfor
			endif
			
		endfor
	endfor
	forprint, input_target.im_file, input_target.mjd_floor, input_target.filter, input_target.tile, input_target.dither, format='A,4X,F0,4X,A,4X,A,4X,A', textout=2

	gv=where(input_target.tile NE ' ' AND input_target.dither NE ' ', n_gv)
	if n_gv GT 0 then begin
		openw, lun, 'survey_target_'+do_program+'.dat', /get_lun
		printf, lun, '#  im_file        filter  tile  dither  weight_file     zp      FWHM    MJD'
		for i=0L, n_gv-1 do begin
			printf, lun, input_target[gv[i]].im_file, input_target[gv[i]].filter, input_target[gv[i]].tile, input_target[gv[i]].dither, input_target[gv[i]].weight_file, input_target[gv[i]].zp, input_target[gv[i]].fwhm, input_target[gv[i]].mjd, input_target[gv[i]].exptime , FORMAT='(A,4X,A,4X,A,4X,A,4X,A,4X,F5.2,4X,F4.1,4X,F0.6,4X,F0.1)'
		endfor
		free_lun, lun
	endif

	im_file=file_search(input_calib_dir, '*.fz')
	n_im_file=n_elements(im_file)
	im_date=list()
	im_type=list()

	for i=0L, n_im_file-1 do begin
		im_h=headfits(im_file[i])
		im_date.add, fxpar(im_h, 'DATE-OBS')
		im_type.add, fxpar(im_h, 'PRODTYPE')
	endfor
	im_date=im_date.toarray(type='string')
	im_type=im_type.toarray(type='string')
	gv=sort(im_date+' '+im_type)

	openw, lun, input_calib_dir+'/pipeline_database.dat', /get_lun
	printf, lun, '# FILE	 DATE-OBS   MJD-OBS	 EXPNUM	 OBJECT	 FILTER	 EXPTIME	 SEQID	 OBSTYPE	 PROCTYPE		PRODTYPE	 MAGZERO'
	for i=0L, n_im_file-1 do begin
		command='dfits '+im_file[gv[i]]+' | fitsort -d DATE-OBS MJD-OBS EXPNUM OBJECT FILTER EXPTIME SEQID OBSTYPE PROCTYPE PRODTYPE MAGZERO'
		spawn, command, result
		result=repstr(result, input_calib_dir+'/', '')
		printf, lun, result
	endfor
	free_lun, lun

endif else $
if recipe EQ 'setup' then begin

	if file_test(input_im_dir, /directory) EQ 0 then file_mkdir, input_im_dir
	if file_test(sextractor_dir, /directory) EQ 0 then file_mkdir, sextractor_dir
	if file_test(sextractor_check_dir, /directory) EQ 0 then file_mkdir, sextractor_check_dir
	if file_test(scamp_dir, /directory) EQ 0 then file_mkdir, scamp_dir
	if file_test(psfex_dir, /directory) EQ 0 then file_mkdir, psfex_dir
	if file_test(swarp_dir, /directory) EQ 0 then file_mkdir, swarp_dir
	if file_test(calib_dir, /directory) EQ 0 then file_mkdir, calib_dir
	if do_overwrite then begin
		file_delete, input_target.im_file, /noexpand, /allow_non, /quiet
		file_delete, input_target.weight_file, /noexpand, /allow_non, /quiet
	endif

	for i=0, n_elements(input_target)-1 do begin

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
		endif

		if file_test(input_calib[i].weight_orig_file, /regular) EQ 0 then begin
			command='funpack '+input_calib[i].weight_orig_file+'.fz'
			print, command
			spawn, command
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

		print, input_calib[i].im_file, date_conv(im_mjd+2400000.5D,'S'), im_zp, im_exptime, im_airmass, im_filter, FORMAT='("Image filename ",A,2X,A,2X,F0.2,2X,F0.2,2X,F0.2,2X,A)'

		case input_calib[i].filter of
			'u': begin
				sex_satur_level='30000.'
				sex_mag_range=[7,15.]
				sex_flux_radius_min=2
				plot_mag_range=[18,5]
				plot_flux_radius_range=[2,8]
				end
			'g': begin
				sex_satur_level='30000.'
				sex_mag_range=[11,19.]
				end
			'i': begin
				sex_satur_level='30000.'
				sex_mag_range=[12,20.]
				end
			'z': begin
				sex_satur_level='30000.'
				sex_mag_range=[12,20.]
				end
			else: stop
		endcase

		im_h=headfits(input_calib[i].im_file)
		im_ra=ten(fxpar(im_h, 'RA'))*360./24.
		im_dec=ten(fxpar(im_h, 'DEC'))
		im_dec_sign=im_dec GT 0. ? '+' : '-'
		im_exptime=fxpar(im_h, 'EXPTIME')
		im_zp=fxpar(im_h, 'MAGZERO')
		im_mjd=fxpar(im_h, 'MJD-OBS')
		im_airmass=tai2airmass(im_ra,im_dec,2000., mjd=im_mjd, lon=-70.8065, latitude=-30.1692, alt=2243.31)
		input_calib[i].airmass=im_airmass

		if file_test(input_calib[i].sex_cat_file, /regular) then fits_info, input_calib[i].sex_cat_file, n_ext=sex_n_ext, /silent $
		else sex_n_ext=0

		if sex_n_ext EQ input_n_chip*2 AND file_test(input_calib[i].sex_zp_file) AND file_test(input_calib[i].sex_zp_full_file) then continue

		im_fwhm=list()
		openw, lun_zp, input_calib[i].sex_zp_file, /get_lun
		openw, lun_zp_full, input_calib[i].sex_zp_full_file, /get_lun
		printf, lun_zp, '# Ext  DECam_zp   airmass   zp   zp_err  zp_nstars   FWHM'
		printf, lun_zp_full, '# Ext  NUMBER   RA   DEC   FLUX_APER   FLUX_AUTO   mag_ref   mag_diff   airmass'

		if file_test(input_calib[i].sex_cat_file, /regular) then fits_info, input_calib[i].sex_cat_file, n_ext=sex_n_ext, /silent $
		else sex_n_ext=0

		if sex_n_ext NE input_n_chip*2 OR do_overwrite EQ 1 then begin

			for j=0L, n_elements(input_chip_fwhm)-1 do begin
				im_data=readfits(input_calib[i].im_file, im_h, ext=input_chip_fwhm[j])
				writefits, sextractor_dir+'/im_sextractor.fits', im_data, im_h
				wim_data=readfits(input_calib[i].weight_file, wim_h, ext=input_chip_fwhm[j])
				writefits, sextractor_dir+'/wim_sextractor.fits', wim_data, wim_h
		
				im_size=size(im_data, /dim)
				im_gain=fxpar(im_h, 'GAINA')
				im_ron=fxpar(im_h, 'RDNOISEA')
	
				command = 'sex '+sextractor_dir+'/im_sextractor.fits' +' -c sex_config/ctio_decam.sex -CATALOG_NAME '+sextractor_dir+'/im_sextractor.ldac'+' -WEIGHT_IMAGE '+sextractor_dir+'/wim_sextractor.fits'+' -XML_NAME '+input_calib[i].sex_xml_file+' -SATUR_LEVEL '+sex_satur_level+' -MAG_ZEROPOINT '+string(input_calib[i].zp,FORMAT='(F0.2)')+' -CHECKIMAGE_TYPE BACKGROUND -CHECKIMAGE_NAME '+input_calib[i].sex_check_file
				print, command
				spawn, command
		
				cat_sex=mrdfits(sextractor_dir+'/im_sextractor.ldac', 2, cat_sex_h, COLUMNS=['NUMBER','X_IMAGE','Y_IMAGE','FLUX_RADIUS','MAG_AUTO','FLUX_AUTO','FLAGS'], /silent)
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

						getpsf, vig_data, im_cx, im_cy, vig_mag, vig_sky, im_ron, im_gain, psf_param, psf_residuals, [0], vig_psfrad, vig_fitrad, sextractor_dir+'/im_sextractor_psf.fits' 
		        im_fwhm.add, 2*sqrt(2*alog(2))*sqrt( (psf_param[3]^2 + psf_param[4]^2)/2. )
						print, 'FWHM ', im_fwhm[-1],' pixels'
	
						dist_circle, vig_mask, vig_size, im_cx, im_cy
						plot, vig_mask, vig_data-vig_sky, xrange=[0,30], psym=1
						x=findgen(100)/10.
						y=gaussian(x,[psf_param[0], 0., mean(psf_param[3:4])])
						oplot, x, y, line=2, color=200
						oplot, im_fwhm[-1]/2.*[1,1], [-1e5,1e5], line=2, color=200
						
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

			command = 'sex '+input_calib[i].im_file +' -c sex_config/ctio_decam.sex -CATALOG_NAME '+input_calib[i].sex_cat_file+' -WEIGHT_IMAGE '+input_calib[i].weight_file+' -XML_NAME '+input_calib[i].sex_xml_file+' -SATUR_LEVEL '+sex_satur_level+' -MAG_ZEROPOINT '+string(input_calib[i].zp,FORMAT='(F0.2)')+' -CHECKIMAGE_TYPE BACKGROUND -CHECKIMAGE_NAME '+input_calib[i].sex_check_file+' -PHOT_APERTURES '+strjoin(string(2*[1,2,3,4,5]*im_fwhm, FORMAT='(F0.2)'),',')
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
		if do_standard EQ 'sdss' then begin
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
		endif else $
			stop

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

		
		for j=0L, input_n_chip-1 do begin
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

		  if n_gv_match GT 2 then begin
				if do_debug then print, 'Running MATCH'
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

	gv_phot=where(input_calib.airmass LT 2. AND input_calib.mjd GT 0. AND input_calib.photflag EQ 'T', n_gv_phot)
	for i=0L, n_gv_phot-1 do begin
		print, 'Processing file ', input_calib[gv_phot[i]].sex_zp_file

		readcol, input_calib[gv_phot[i]].sex_zp_file, temp_ext, temp_decam_zp, temp_airmass, temp_zp, temp_zperr, temp_zp_nstars, temp_fwhm, FORMAT='I,F,F,F,F,I,F'
		create_struct, temp_struct, '', ['ext','zp','airmass','filter','mjd','photflag'], 'I,F,F,A,F,A', dim=n_elements(temp_ext)
		temp_struct.ext=temp_ext
		temp_struct.zp=temp_zp
		temp_struct.airmass=temp_airmass
		temp_struct.filter=input_calib[gv_phot[i]].filter
		temp_struct.mjd=input_calib[gv_phot[i]].mjd
		temp_struct.photflag=input_calib[gv_phot[i]].photflag
		if i EQ 0L then pipe_zp=temp_struct $
		else pipe_zp=[pipe_zp,temp_struct]

		readcol, input_calib[gv_phot[i]].sex_zp_full_file, temp_ext, temp_id, temp_ra, temp_dec, temp_flux_aper, temp_flux_auto, temp_mag_ref, temp_zp, temp_airmass, FORMAT='I,I,F,F,D,D,F,F,F'
		create_struct, temp_struct, '', ['ext','zp','airmass','filter','mjd'], 'I,F,F,A,F', dim=n_elements(temp_ext)
		temp_struct.ext=temp_ext
		temp_struct.zp=temp_zp
		temp_struct.airmass=temp_airmass
		temp_struct.filter=input_calib[gv_phot[i]].filter
		temp_struct.mjd=input_calib[gv_phot[i]].mjd
		if i EQ 0L then pipe_zp_full=temp_struct $
		else pipe_zp_full=[pipe_zp_full,temp_struct]

	endfor

	filter_uniq=pipe_zp[uniq(pipe_zp.filter, sort(pipe_zp.filter))].filter
	mjd_uniq=floor(pipe_zp[uniq(floor(pipe_zp.mjd+2./24), sort(floor(pipe_zp.mjd+2./24)))].mjd)

	openw, lun, 'survey_zp_'+do_program+'.dat', /get_lun
	printf, lun, '# MJD		filter		zp		k		zp_err		k_err		PHOTFLAG'
	for i=0L, n_elements(filter_uniq)-1 do begin
		for j=0L, n_elements(mjd_uniq)-1 do begin

			gv_pipe=where(pipe_zp.filter EQ filter_uniq[i] AND floor(pipe_zp.mjd+2./24) EQ mjd_uniq[j], n_gv_pipe)
			gv_pipe_phot=where(pipe_zp.filter EQ filter_uniq[i] AND floor(pipe_zp.mjd+2./24) EQ mjd_uniq[j] AND pipe_zp.photflag EQ 'T', n_gv_pipe_phot)
			coeff=robust_linefit(pipe_zp[gv_pipe_phot].airmass, pipe_zp[gv_pipe_phot].zp, temp_yfit, temp_sigma, coeff_error)
			photflag = coeff_error[0] LT 0.015 ? 'T' : 'F'
			printf, lun, mjd_uniq[j], filter_uniq[i], coeff[0], coeff[1], coeff_error[0], coeff_error[1], photflag, FORMAT='(F0.1,4X,A,4X,F0.2,4X,F0.2,4X,F0.2,4X,F0.2,4X,A)'

			cgloadct, 0
			cgwindow, wxsize=800, wysize=600
			cgplot, [0], [0], xrange=[1.,1.9], yrange=mean(pipe_zp[gv_pipe].zp)+[-1,1], /nodata, /window, xtitle='airmass', ytitle='mag_std - mag_inst', title='Filter '+filter_uniq[i]+' - Date '+date_conv(mjd_uniq[j]+2400000.5D,'S')
			cgplot, pipe_zp[gv_pipe].airmass, pipe_zp[gv_pipe].zp, psym=cgsymcat('filled circle'), symsize=0.4, color='blue', /over, /addcmd
			cgplot, pipe_zp[gv_pipe_phot].airmass, pipe_zp[gv_pipe_phot].zp, psym=cgsymcat('filled circle'), symsize=0.4, color='red', /over, /addcmd
			cgplot, [0.,100.], coeff[0]+coeff[1]*[0.,100.], color='black', line=2, /over, /addcmd
			cglegend, color='red', align=0, location=[0.14,0.86], length=0.02, titles=filter_uniq[i]+'-band coeff='+string(coeff,format='(F0.3,",",F0.3)')+' coeff_error='+string(coeff_error,format='(F0.3,",",F0.3)'), vspace=2., /addcmd
			cgcontrol, output=calib_dir+'/standard_zp_airmass_uncorrected_'+filter_uniq[i]+'_'+string(mjd_uniq[j],FORMAT='(I0)')+'.pdf'

			print, 'Airmass term compute - Filter ', filter_uniq[i], ' Date ', date_conv(mjd_uniq[j]+2400000.5D,'S')
			print, coeff[0], coeff_error[0], coeff[1], coeff_error[1], FORMAT='("zp = ",F0.3," +- ",F0.3," / k = ", F0.3, " +- ", F0.3)'

			gv_pipe=where(input_calib.filter EQ filter_uniq[i] AND floor(input_calib.mjd+2./24) EQ mjd_uniq[j], n_gv_pipe)
			gv_pipe_phot=where(input_calib.filter EQ filter_uniq[i] AND floor(input_calib.mjd+2./24) EQ mjd_uniq[j] AND input_calib.photflag EQ 'T', n_gv_pipe_phot)
			cgloadct, 0
			cgwindow, wxsize=800, wysize=600
			cgplot, [0], [0], xrange=mjd_uniq[j]+[-2,12]/24., yrange=[1.,1.8], /nodata, /window, xtitle='MJD', ytitle='airmass', title='Filter '+filter_uniq[i]+' - Date '+date_conv(mjd_uniq[j]+2400000.5D,'S'), XTICKFORMAT='(F0.1)'
			cgplot, input_calib[gv_pipe].mjd, input_calib[gv_pipe].airmass, psym=cgsymcat('filled circle'), symsize=1., color='blue', /over, /addcmd
			cgplot, input_calib[gv_pipe_phot].mjd, input_calib[gv_pipe_phot].airmass, psym=cgsymcat('filled circle'), symsize=1., color='red', /over, /addcmd
			cgcontrol, output=calib_dir+'/standard_airmass_mjd_'+filter_uniq[i]+'_'+string(mjd_uniq[j],FORMAT='(I0)')+'.pdf'

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

	for i=0L, n_elements(input_target)-1 do begin
		case input_target[i].filter of
			'u': sex_satur_level='12000.'
			'g': sex_satur_level='12000.'
			'i': sex_satur_level='19000.'
			'z': sex_satur_level='19000.'
			else: stop
		endcase

		command = 'sex '+input_target[i].im_file +' -c sex_config/ctio_decam.sex -CATALOG_NAME '+input_target[i].sex_cat_file+' -WEIGHT_IMAGE '+input_target[i].weight_file+' -XML_NAME '+input_target[i].sex_xml_file+' -SATUR_LEVEL '+sex_satur_level+' -MAG_ZEROPOINT '+string(input_target[i].zp,FORMAT='(F0.2)')+' -CHECKIMAGE_TYPE BACKGROUND -CHECKIMAGE_NAME '+input_target[i].sex_check_file
		print, command
		spawn, command
	endfor

endif else $
if recipe EQ 'sextractor standard' then begin

	case do_filter of
		'g': sex_satur_level='30000.'
		'i': sex_satur_level='30000.'
		else: stop
	endcase
 
	file_delete, input_calib.sex_cat_file, /noexpand, /allow_non, /quiet
  file_delete, input_calib.scamp_head_file, /noexpand, /allow_non, /quiet

	for i=0L, n_elements(input_calib)-1 do begin
		command = 'sex '+input_calib[i].im_file +' -c sex_config/ctio_decam.sex -CATALOG_NAME '+input_calib[i].sex_cat_file+' -WEIGHT_IMAGE '+input_calib[i].weight_file+' -XML_NAME '+input_calib[i].sex_xml_file+' -SATUR_LEVEL '+sex_satur_level+' -MAG_ZEROPOINT '+string(input_calib[i].zp,FORMAT='(F0.2)')+' -CHECKIMAGE_TYPE BACKGROUND -CHECKIMAGE_NAME '+input_calib[i].sex_check_file
		print, command
		spawn, command
	endfor

endif else $
if recipe EQ 'scamp' then begin

  scamp_pos_error='1.4'
  scamp_scale_error='1.01'
  scamp_angle_error='0.04'
  scamp_sn_thresholds='40.,80.0'
	scamp_fwhm_thresholds='0.,100.'
	scamp_crossid_radius='6.'
	scamp_distort_degrees='4'
	scamp_astref_catalog = 'USNO-B1' ;'USNO-B1' ;'2MASS'
	scamp_astref_band = 'Rf' ;'DEFAULT' ;'Ks'
	scamp_astrefmag_limits='6.,20.' ;'-99.,18.' 
	scamp_match_resol='0.'

	scamp_dir_out=scamp_dir+'/order_'+scamp_distort_degrees+'_refcat_'+scamp_astref_catalog+'_filter_'+do_filter
	scamp_list_file=scamp_dir+'/scamp_decam_'+do_filter+'.dat'
	scamp_cat_file_out = scamp_dir_out+'/scamp_decam_'+do_filter+'.ldac'
  scamp_cat_type_out = 'FITS_LDAC'
  scamp_xml_file = scamp_dir_out+'/scamp_decam_'+do_filter+'.xml'
  scamp_refcat_dir = scamp_dir+'/refcat'
  scamp_check_type = 'SKY_ALL,FGROUPS,DISTORTION,ASTR_INTERROR2D,ASTR_INTERROR1D,ASTR_REFERROR2D,ASTR_REFERROR1D,ASTR_CHI2,PHOT_ERROR'
  scamp_check_file = strjoin(scamp_dir_out+'/'+['sky_all','fgroups','distort','astr_interror2d','astr_interror1d','astr_referror2d','astr_referror1d','astr_chi2','psphot_error'],',')

	scamp_mag_bin=0.5
	scamp_radius_bin=0.2

  file_delete, input_target.scamp_head_file, /noexpand, /allow_non, /quiet
  if do_overwrite then file_delete, input_target.scamp_cat_file, /noexpand, /allow_non, /quiet
	if not file_test(scamp_refcat_dir, /directory) then file_mkdir, scamp_refcat_dir, /noexpand_path
	if not file_test(scamp_dir_out, /directory) then file_mkdir, scamp_dir_out, /noexpand_path
	loadct, 12

	for i=0L, n_elements(input_target)-1 do begin
		print, 'Creating file '+input_target[i].scamp_cat_file

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
			end
			'i': begin
				scamp_mag_range=[15.,21.]
				plot_mag_range=[25,15]
			end
			else: stop
		endcase
	
		if file_test(input_target[i].sex_cat_file) EQ 0 then stop

		if file_test(input_target[i].scamp_cat_file, /regular) then fits_info, input_target[i].scamp_cat_file, n_ext=scamp_n_ext, /silent $
		else scamp_n_ext=0
		if scamp_n_ext EQ input_n_chip*2 AND do_overwrite EQ 0 then continue

		fits_open, input_target[i].sex_cat_file, fcb_in  ; Lee la tabla fits sin intervenerla           
	  fits_read, fcb_in, cat_data0, cat_h0, exten=0 
  	writefits, input_target[i].scamp_cat_file, cat_data0, cat_h0

		for j=0L, input_n_chip-1 do begin
			cat_sex=mrdfits(input_target[i].sex_cat_file, 2*(j+1), cat_sex_h, COLUMNS=['NUMBER','FLUX_RADIUS','MAG_AUTO','MAGERR_AUTO','FLAGS','A_IMAGE','B_IMAGE'], /silent)
			plothist, (cat_sex[where(cat_sex.mag_auto GT scamp_mag_range[0] and cat_sex.mag_auto LT scamp_mag_range[1] AND cat_sex.flux_radius GT 2 AND cat_sex.flags LE 3 AND cat_sex.b_image/cat_sex.a_image GT 0.8, n_gv)]).flux_radius, temp_xhist, temp_yhist, bin=scamp_radius_bin, /noplot
			temp=max(temp_yhist, gv) & scamp_radius=temp_xhist[gv]
			scamp_radius=median((cat_sex[where(cat_sex.mag_auto GT scamp_mag_range[0] AND cat_sex.mag_auto LT scamp_mag_range[1] AND cat_sex.flux_radius GT scamp_radius*0.9 AND cat_sex.flux_radius LT scamp_radius*1.1 , n_gv)]).flux_radius)
			print, j+1, scamp_radius, FORMAT='("Chip ",I2,"  Flux_radius ",F0.1)'

			plot, cat_sex.flux_radius, cat_sex.mag_auto, psym=1, xrange=[1,6], yrange=plot_mag_range
			oplot, scamp_radius*[1,1], [0,100], color=200
			oplot, scamp_radius*[0.9,0.9], [0,100], line=2, color=200
			oplot, scamp_radius*[1.1,1.1], [0,100], line=2, color=200
			oplot, [0,100], scamp_mag_range[0]*[1,1], line=2, color=100
			oplot, [0,100], scamp_mag_range[1]*[1,1], line=2, color=100
			wait, 0.2
			if j EQ 0 AND do_debug then begin
				forprint, temp_xhist, temp_yhist, text=2
				print, 'DEBUG - CHECK mag range and flux radius are ok'
				stop
			endif

			gv_stars=where(cat_sex.mag_auto GT scamp_mag_range[0] AND cat_sex.mag_auto LT scamp_mag_range[1] AND cat_sex.flux_radius GT scamp_radius*0.9 AND cat_sex.flux_radius LT scamp_radius*1.1 AND cat_sex.flags LE 1, n_gv_stars)

  		fits_read, fcb_in, cat_data1, cat_h1, exten=2*j+1                                   
			fits_read, fcb_in, cat_data2, cat_h2, exten=2*j+2

			cat_data1=reform(cat_data1, [n_elements(cat_data1),1])
			cat_data2=cat_data2[*,[gv_stars]]
			fxaddpar, cat_h2, 'NAXIS2', (size(cat_data2, /dim))[1]

  		fits_open, input_target[i].scamp_cat_file, fcb_out, /append
		  fits_write, fcb_out, cat_data1, cat_h1
			fits_write, fcb_out, cat_data2, cat_h2
  		fits_close, fcb_out

		endfor
  	fits_close, fcb_in
	endfor

	forprint, input_target.scamp_cat_file, textout=scamp_list_file, FORMAT='(A)', /NOCOMMENT

	for i=0L, n_elements(input_target)-1 do begin
		im_h=headfits(input_target[i].im_file)
		im_mjd=fxpar(im_h, 'MJD-OBS')
		input_target[i].mjd=im_mjd
		print, input_target[i].im_file, input_target[i].filter, im_mjd

		openw, lun, input_target[i].scamp_ahead_file, /get_lun
		for j=0L, input_n_chip-1 do begin
			im_h=headfits(input_target[i].im_file, ext=j+1)
			extast, im_h, im_ast
			xy2ad, im_ast.naxis[0]/2., im_ast.naxis[1]/2., im_ast, im_ra, im_dec
			im_airmass=tai2airmass(im_ra, im_dec, 2000., mjd=im_mjd, lon=-70.8065, latitude=-30.1692, alt=2243.31)

			printf, lun, 'FILTER  =     '+input_target[i].filter
			printf, lun, 'AIRMASS =     '+string(im_airmass,FORMAT='(F0.3)')
			printf, lun, 'EXPTIME =     '+string(input_target[i].exptime,FORMAT='(F0.2)')  ;string(600.,FORMAT='(F0.2)')
			if n_elements(input_zp) GT 0 then begin
				gv=where(input_zp.filter EQ input_target[i].filter AND input_zp.mjd EQ floor(input_target[i].mjd+2./24), n_gv)
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
	endfor

	command='scamp @'+scamp_list_file+' -c scamp_config/ctio_decam.scamp'+' -MERGEDOUTCAT_TYPE '+scamp_cat_type_out+' -MERGEDOUTCAT_NAME '+scamp_cat_file_out+' -MATCH Y -WRITE_XML Y -XML_NAME '+scamp_xml_file+' -SAVE_REFCATALOG Y -REFOUT_CATPATH '+scamp_refcat_dir+' -CHECKPLOT_DEV PSC -CHECKPLOT_ANTIALIAS Y -CHECKPLOT_TYPE '+scamp_check_type + ' -CHECKPLOT_NAME '+scamp_check_file + ' -ASTREF_CATALOG '+scamp_astref_catalog+' -ASTREF_BAND '+scamp_astref_band+' -ASTREFMAG_LIMITS '+scamp_astrefmag_limits+' -DISTORT_DEGREES '+scamp_distort_degrees+' -PHOTCLIP_NSIGMA 2. -SOLVE_ASTROM Y -SOLVE_PHOTOM Y -POSITION_MAXERR '+scamp_pos_error+' -PIXSCALE_MAXERR '+scamp_scale_error+' -POSANGLE_MAXERR '+scamp_angle_error+' -SN_THRESHOLDS '+scamp_sn_thresholds+' -FWHM_THRESHOLDS '+scamp_fwhm_thresholds+' -CROSSID_RADIUS '+scamp_crossid_radius+' -MATCH_RESOL '+scamp_match_resol

	print, command
	spawn, command

endif else $
if recipe EQ 'swarp' then begin

	swarp_list_file=swarp_dir+'/swarp_decam_'+do_filter+'.dat'

	swarp_combine_type='AVERAGE'
	swarp_resampling_type='LANCZOS2'
	swarp_resample_do='Y'
	swarp_resample_dir=swarp_dir+'/resample'
	swarp_weight_suffix='.WEIGHT.fits'
	swarp_xml_file=swarp_dir+'/swarp_decam_'+do_filter+'.xml'
	swarp_verbose='NORMAL'

	swarp_im_out=swarp_dir+'/fornax_tile'+do_tile+'_'+do_filter+'.fits'
	swarp_weight_out=swarp_dir+'/fornax_tile'+do_tile+'_'+do_filter+'.WEIGHT.fits'
	
	sex_stack_im_file=swarp_im_out
	sex_stack_weight_file=swarp_weight_out
	sex_stack_cat_file=sextractor_dir+'/fornax_tile'+do_tile+'_'+do_filter+'.ldac'
	sex_stack_xml_file=sextractor_dir+'/fornax_tile'+do_tile+'_'+do_filter+'.xml'
	sex_stack_checkimage_type='BACKGROUND,SEGMENTATION'
	sex_stack_checkimage_file=strjoin(sextractor_check_dir+'/fornax_tile'+do_tile+'_'+do_filter+['.CHECK_BACK.fits','.CHECK_SEGMENTATION.fits'],',')

	for i=0L, n_elements(input_target)-1 do begin
		file_copy, input_target[i].scamp_head_file, input_target[i].swarp_head_file, /OVERWRITE
	endfor
	if not file_test(swarp_resample_dir, /directory) then file_mkdir, swarp_resample_dir, /noexpand_path
	forprint, input_target.swarp_im_file, textout=swarp_list_file, FORMAT='(A)', /NOCOMMENT

	command='swarp @'+swarp_list_file+' -c swarp_config/ctio_decam.swarp'+' -IMAGEOUT_NAME '+swarp_im_out+' -WEIGHTOUT_NAME '+swarp_weight_out+' -COMBINE_TYPE '+swarp_combine_type+' -RESAMPLE '+swarp_resample_do+' -RESAMPLE_DIR '+swarp_resample_dir+' -SATLEV_DEFAULT 10000 -WEIGHT_TYPE MAP_WEIGHT -WEIGHT_SUFFIX '+swarp_weight_suffix+' -WEIGHT_THRESH 0. -RESCALE_WEIGHTS N -BLANK_BADPIXELS Y -WRITE_XML Y -XML_NAME '+swarp_xml_file+' -VERBOSE_TYPE ' +swarp_verbose+' -RESAMPLING_TYPE '+swarp_resampling_type+' -SUBTRACT_BACK Y -BACK_SIZE 384'
	print, command
	spawn, command

	command='sex '+sex_stack_im_file+' -c sex_config/ctio_decam_stack.sex -CATALOG_NAME '+sex_stack_cat_file+' -WEIGHT_IMAGE '+sex_stack_weight_file+' -MAG_ZEROPOINT 30. -XML_NAME '+sex_stack_xml_file+' -CHECKIMAGE_TYPE '+sex_stack_checkimage_type+' -CHECKIMAGE_NAME '+sex_stack_checkimage_file
	print, command
	spawn, command

endif else $
if recipe EQ 'psfex' then begin

	case do_filter of
		'u': begin
			psfex_mag_range=[15.,18.]+2.5*alog10(600.)
			plot_mag_range=[21,10]
		end
		'g': begin
			psfex_mag_range=[16.7,19.5]
			plot_mag_range=[25,15]
		end
		'i': begin
			psfex_mag_range=[15.,20.]
			plot_mag_range=[25,15]
		end
		else: stop
	endcase

	sex_stack_im_file=swarp_dir+'/fornax_tile'+do_tile+'_'+do_filter+'.fits'
	sex_stack_weight_file=swarp_dir+'/fornax_tile'+do_tile+'_'+do_filter+'.WEIGHT.fits'
	sex_stack_cat_file=sextractor_dir+'/fornax_tile'+do_tile+'_'+do_filter+'_sex_psfex.ldac'
	sex_stack_xml_file=sextractor_dir+'/fornax_tile'+do_tile+'_'+do_filter+'_sex_psfex.xml'
	sex_stack_checkimage_type='BACKGROUND,SEGMENTATION'
	sex_stack_checkimage_file=strjoin(sextractor_check_dir+'/fornax_tile'+do_tile+'_'+do_filter+['.CHECK_BACK.fits','.CHECK_SEGMENTATION.fits'],',')

	psfex_stack_cat_file=sextractor_dir+'/fornax_tile'+do_tile+'_'+do_filter+'_psfex.ldac'
	psfex_stack_xml_file=sextractor_dir+'/fornax_tile'+do_tile+'_'+do_filter+'_psfex.xml'
	psfex_check_type = 'RESIDUALS,SNAPSHOTS,SAMPLES'
	psfex_check_file = strjoin(psfex_dir+'/psfex_checkimage'+['.CHECK_RESIDUALS.fits','.CHECK_SNAPSHOTS.fits','.CHECK_SAMPLES.fits'],',')
	psfex_checkplot_type = 'FWHM,ELLIPTICITY,COUNTS'
	psfex_checkplot_file = strjoin(psfex_dir+'/psfex_checkplot'+['.CHECK_FWHM','.CHECK_ELLIPTICITY','.CHECK_COUNTS'],',')
	psfex_psfvar_nsnap='20'
	psfex_psfvar_degrees='3'
	psfex_basis_type='PIXEL'
	psfex_basis_number='20'
	psfex_psf_size='47,47'
	psfex_psf_sampling='0.8'

	cat_stars_data=mrdfits('catalogs/psfex_stars_v2.fits', 1, cat_stars_h, COLUMNS=['RA','DEC','MAG_AUTO','FLAGS'])
	gv_stars=where(cat_stars_data.flags[0] LE 3 AND cat_stars_data.flags[1] LE 3 AND cat_stars_data.flags[2] LE 3, n_gv_stars)

	if file_test(sex_stack_cat_file) EQ 0 OR do_overwrite then begin
		command='sex '+sex_stack_im_file+' -c sex_config/ctio_decam_psfex.sex -PARAMETERS_NAME sex_config/ctio_decam_psfex.param -CATALOG_NAME '+sex_stack_cat_file+' -WEIGHT_IMAGE '+sex_stack_weight_file+' -MAG_ZEROPOINT 30. -XML_NAME '+sex_stack_xml_file+' -CHECKIMAGE_TYPE '+sex_stack_checkimage_type+' -CHECKIMAGE_NAME '+sex_stack_checkimage_file
		print, command
		spawn, command
	endif

	cat_sex_data=mrdfits(sex_stack_cat_file, 2, cat_sex_h, COLUMNS=['X_IMAGE','Y_IMAGE','ALPHA_J2000','DELTA_J2000','FLUX_APER','FLAGS','FLUX_RADIUS','ELONGATION'])
	gv_sex=where(cat_sex_data.flags LE 3 AND (30.-2.5*alog10(cat_sex_data.flux_aper)) GE psfex_mag_range[0] AND (30.-2.5*alog10(cat_sex_data.flux_aper)) LT psfex_mag_range[1], n_gv_sex)

	match_nmax=10000L
	j=0L
	gv_sex_match=list()
	gv_stars_match=list()
	repeat begin
		j1=j*match_nmax
		j2=((j+1)*match_nmax-1) < (n_gv_sex-1)
		gv_sex_j=gv_sex[j1:j2]
		n_gv_sex_j=j2-j1+1

		gv_match= where( sqrt(min( (cat_sex_data[gv_sex_j].alpha_j2000#make_array(n_gv_stars,value=1.,/double) - make_array(n_gv_sex_j,value=1.,/double)#cat_stars_data[gv_stars].ra)^2 + (cat_sex_data[gv_sex_j].delta_j2000#make_array(n_gv_stars,value=1.,/double) - make_array(n_gv_sex_j,value=1.,/double)#cat_stars_data[gv_stars].dec)^2, id_match, dim=2)) LT 2./3600, n_gv_match)
		if n_gv_match GT 0 then begin
			gv_sex_match.add, gv_sex_j[(id_match[gv_match] mod n_gv_sex_j)], /extract
			gv_stars_match.add, gv_stars[(id_match[gv_match]/n_gv_sex_j)], /extract
		endif
		gv_match=0
		gv_sex_j=0
		j++
	endrep until long(j*match_nmax) GT n_gv_sex

	cat_sex_data=0
	cat_sex_h=0

	gv_sex_match=gv_sex_match.toarray(type='long')
	gv_stars_match=gv_stars_match.toarray(type='long')
	print, 'Number of mathced stars ', n_elements(gv_sex_match)

	fits_open, sex_stack_cat_file, cat_fcb
	fits_read, cat_fcb, cat_data0, cat_h0, exten=0
	fits_read, cat_fcb, cat_data1, cat_h1, exten=1

	cat_data2=make_array(cat_fcb.axis[0:1,2], value=0, /byte)	
	temp_lines=10000L
	temp_max=double(product(cat_fcb.axis[0:1,2])-1.)
	k=0L
	repeat begin
		print, 'Reading line ', strn(k*temp_lines), ' of ', strn(cat_fcb.axis[1,2])
		k1=double(k*temp_lines*cat_fcb.axis[0,2])
		k2=double((k+1)*temp_lines*cat_fcb.axis[0,2] - 1) < temp_max
		fits_read, cat_fcb, temp_data, cat_h2, first=k1, last=k2, exten_no=2
		cat_data2[*,k1/cat_fcb.axis[0,2]:(k2+1)/cat_fcb.axis[0,2]-1]=reform(temp_data, [cat_fcb.axis[0,2],(k2-k1+1)/cat_fcb.axis[0,2]])
		k++
	endrep until double(k*temp_lines*cat_fcb.axis[0,2]) GT temp_max
	temp_data=0
	fits_close, fcb
	fits_close, cat_fcb

	cat_data1=reform(cat_data1, [n_elements(cat_data1),1])
	cat_data2=cat_data2[*,[gv_sex_match]]
	fxaddpar, cat_h2, 'NAXIS2', (size(cat_data2, /dim))[1]


	if file_test(psfex_stack_cat_file, /regular, /noexpand) then  file_delete, psfex_stack_cat_file, /noexpand
	writefits, psfex_stack_cat_file, cat_data0, cat_h0
	fits_open, psfex_stack_cat_file, cat_fcb, /append
	fits_write, cat_fcb, cat_data1, cat_h1
	fits_write, cat_fcb, cat_data2, cat_h2
	fits_close, cat_fcb
	cat_data0=0 & cat_data1=0 & cat_data2=0

	command='psfex '+psfex_stack_cat_file+' -c psfex_config/ctio_decam.psfex'+' -PSF_DIR '+psfex_dir+' -WRITE_XML Y -XML_NAME '+psfex_stack_xml_file+' -CHECKIMAGE_TYPE NONE '+' -CHECKIMAGE_NAME '+psfex_check_file+' -CHECKPLOT_DEV PNG -CHECKPLOT_TYPE '+psfex_checkplot_type+' -CHECKPLOT_NAME '+psfex_checkplot_file+' -PSFVAR_NSNAP '+psfex_psfvar_nsnap+' -PSFVAR_DEGREES '+psfex_psfvar_degrees+' -BASIS_TYPE '+psfex_basis_type+' -BASIS_NUMBER '+psfex_basis_number+' -SAMPLE_VARIABILITY 1.,1. -SAMPLE_MAXELLIP 1. -NEWBASIS_TYPE NONE -NEWBASIS_NUMBER 10 -SAMPLEVAR_TYPE NONE -STABILITY_TYPE EXPOSURE -SAMPLE_MINSN 1. -SAMPLE_FWHMRANGE 0.1,10. -SAMPLE_AUTOSELECT N -PSF_SIZE '+psfex_psf_size+' -PSF_SAMPLING '+psfex_psf_sampling+' -PSF_ACCURACY 0.02 -BADPIXEL_FILTER Y -BADPIXEL_NMAX 80'
	print, command
	spawn, command

endif else $
if recipe EQ 'sex psf' then begin

	sex_stack_im_file=swarp_dir+'/fornax_tile'+do_tile+'_'+do_filter+'.fits'
	sex_stack_weight_file=swarp_dir+'/fornax_tile'+do_tile+'_'+do_filter+'.WEIGHT.fits'
	sex_stack_cat_file=sextractor_dir+'/fornax_tile'+do_tile+'_'+do_filter+'_psf.ldac'
	sex_stack_xml_file=sextractor_dir+'/fornax_tile'+do_tile+'_'+do_filter+'_psf.xml'
	sex_stack_checkimage_type='BACKGROUND,SEGMENTATION'
	sex_stack_checkimage_file=strjoin(sextractor_check_dir+'/fornax_tile'+do_tile+'_'+do_filter+['.CHECK_BACK.fits','.CHECK_SEGMENTATION.fits'],',')
	psfex_psf_file=psfex_dir+'/fornax_tile'+do_tile+'_'+do_filter+'_psfex.psf'

	command='sex '+sex_stack_im_file+' -c sex_config/ctio_decam_stack.sex -PARAMETERS_NAME sex_config/ctio_decam_stack_psf.param -CATALOG_NAME '+sex_stack_cat_file+' -WEIGHT_IMAGE '+sex_stack_weight_file+' -MAG_ZEROPOINT 30. -XML_NAME '+sex_stack_xml_file+' -CHECKIMAGE_TYPE '+sex_stack_checkimage_type+' -CHECKIMAGE_NAME '+sex_stack_checkimage_file + ' -PSF_NAME '+psfex_psf_file+' -PSF_NMAX 1'
	print, command
	spawn, command

endif

end
