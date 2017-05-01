pro fix_header, im_file 

	fits_info, im_file, /silent, n_ext=n_chip
	rdfits_struct, im_file, im_data_struct, /silent

	im_h=im_data_struct.(0)
	im_ra=ten(fxpar(im_h, 'RA'))*360D/24
	im_dec=ten(fxpar(im_h, 'DEC'))

	if im_ra GT 360. then begin
		print, 'Fixing header of image ', im_file
		im_ra -= 360.
		val=string(sixty(im_ra*24D/360), format='(I,":",I02,":",F05.2)')
		fxaddpar, im_h, 'RA', val
	endif

	mwrfits, 3, im_file, im_h, /create	
	for j=0L, n_chip-1 do begin
		im_h=im_data_struct.((j+1)*2)
		im_data=float(im_data_struct.((j+1)*2+1))

		im_val1=double(fxpar(im_h, 'CRVAL1', /NAN))
		im_val2=double(fxpar(im_h, 'CRVAL2', /NAN))

		if finite(im_val1) then begin
			print, 'CRVAL1 exists - CHIP ', j
			if im_val1 GT 360. then begin
				val = im_val1-360.
				fxaddpar, im_h, 'CRVAL1', val
			endif
		endif $
		else begin
			print, 'CRVAL1 did not exist - CHIP ', j
			val = im_ra
			fxaddpar, im_h, 'CRVAL1', val
		endelse

		if finite(im_val2) then begin
			print, 'CRVAL2 exists - CHIP ', j
		endif $
		else begin
			print, 'CRVAL2 did not exist - CHIP ', j
			val = im_dec
			fxaddpar, im_h, 'CRVAL2', val
		endelse
		
		mwrfits, im_data, im_file, im_h, /silent
	endfor

	return
end