im_file='pipeline/ss_fornax_p92_t1_v01_s3nx_d1a_Ks.fits'
nim_file='pipeline/ss_fornax_p92_t1_v01_s3nx_d1a_Ks_mod.fits'
input_im_size=[2048,2048]
input_n_chip=16

rdfits_struct, im_file, im_data_struct, /silent
mwrfits, 3, nim_file, im_data_struct.(0), /create	
for j=0L, input_n_chip-1 do begin
	im_h=im_data_struct.((j+1)*2)
	im_data=float(im_data_struct.((j+1)*2+1))
	print, size(im_h)
	print, size(im_data)

	extast, im_h, im_ast
	temp_x=[input_im_size[0]/2-1,0,0,0,input_im_size[0]-1,input_im_size[0]-1,input_im_size[0]-1]
	temp_y=[input_im_size[1]/2-1,0,input_im_size[1]/2-1,input_im_size[1]-1,0,input_im_size[1]/2-1,input_im_size[1]-1]
	xy2ad, temp_x, temp_y, im_ast, temp_ra, temp_dec
	nim_ast=solve_astro(temp_ra, temp_dec, temp_x, temp_y)

	putast, im_h, nim_ast
	fxaddpar, im_h, 'PV2_1', 0.
	fxaddpar, im_h, 'PV2_2', 0.
	fxaddpar, im_h, 'PV2_3', 0.
	fxaddpar, im_h, 'PV2_4', 0.
	fxaddpar, im_h, 'PV2_5', 0.
	
	mwrfits, im_data, nim_file, im_h, /silent
endfor

END
