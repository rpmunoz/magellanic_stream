pro filter_sex_cat, cat_in_file, cat_out_file, CLASS=class, HD=hd, CHECK_PLOT=check_plot

	do_class= (n_elements(class) GT 0) ? class : 'stars'
	do_hd= (n_elements(hd) GT 0) ? keyword_set(hd) : 0
	do_check_plot= (n_elements(check_plot) GT 0) ? keyword_set(check_plot) : 0

	if do_check_plot then begin
		if do_hd then begin
		  cgdelete, /all
  		cgwindow, wxsize=800, wysize=600, wxpos=0, wypos=0, wobject=win1
		  cgwindow, wxsize=800, wysize=600, wxpos=600, wypos=0, wobject=win2
		endif $
		else begin
		  cgdelete, /all
 			cgwindow, wxsize=800, wysize=600, wxpos=0, wypos=0, wobject=win1
			window, 1, XSIZE=800, YSIZE=600
		endelse
	endif
	
	cat_data=mrdfits(cat_in_file, 2, cat_h, COLUMNS=['ALPHA_J2000','DELTA_J2000','MAG_AUTO','FLUX_RADIUS','FLAGS','ELONGATION'])

	if do_class EQ 'stars' then begin
		
		gv_flags=where(cat_data.flags LE 3 AND cat_data.mag_auto LT 99., n_gv_flags)
		gv_sort=gv_flags[sort(cat_data[gv_flags].mag_auto)]
		gv_mag=gv_sort[n_gv_flags*0.001:n_gv_flags*0.5]

    n_gv_mag=n_elements(gv_mag)
    gv_sort=gv_mag[sort(cat_data[gv_mag].flux_radius)]
    gv_mag_fwhm=gv_sort[n_gv_mag*0.01:n_gv_mag*0.4]

		temp_mean=biweight_mean(cat_data[gv_mag_fwhm].flux_radius, temp_sigma)
		temp_bin= round( (3.5 * temp_sigma) / n_elements(gv_mag_fwhm)^0.3333 * 100)/100.
		temp_yhist=cghistogram(cat_data[gv_mag_fwhm].flux_radius, binsize=temp_bin, locations=temp_xhist)
    temp_xhist+=temp_bin/2.
    temp_max=max(temp_yhist, gv_max)
    fwhm_mean=temp_xhist[gv_max]

		if do_check_plot then begin
	    cgset, win1
  	  cghistoplot, cat_data[gv_mag_fwhm].flux_radius, binsize=temp_bin, color='red', xtitle='FWHM (pixels)', /window
    	cgplot, fwhm_mean*[1,1], [0,1e6], line=2, color='blue', /addcmd, /over
		endif

    gv_low=where( cat_data[gv_mag_fwhm].flux_radius LE fwhm_mean, n_gv_low)
	  gv_sort=gv_mag_fwhm[gv_low[sort(cat_data[gv_mag_fwhm[gv_low]].flux_radius)]]
    gv_mag_fwhm_low=gv_sort[n_gv_low*0.1:n_gv_low-1]
    fwhm_range=fwhm_mean+abs(fwhm_mean - min(cat_data[gv_mag_fwhm_low].flux_radius))*[-1,1]
    mag_range=min(cat_data[gv_mag_fwhm_low].mag_auto)+[0.2,1.5] ; it was [0.5,2.5]

    gv_low_high=where(cat_data[gv_mag_fwhm].flux_radius GE fwhm_range[0] AND cat_data[gv_mag_fwhm].flux_radius LE fwhm_range[1], n_gv_low_high)
    fwhm_mean=biweight_mean(cat_data[gv_mag_fwhm[gv_low_high]].flux_radius, fwhm_sigma)
    fwhm_sigma >= 0.2
    fwhm_range= fwhm_mean + 2*fwhm_sigma*[-1,1]

    gv_stars=where(cat_data.mag_auto GE mag_range[0] AND cat_data.mag_auto LE mag_range[1] AND cat_data.flux_radius GE fwhm_range[0] AND cat_data.flux_radius LE fwhm_range[1] AND cat_data.flags EQ 0, n_gv_stars)

    plot_xrange=[-0.5,4.*median(cat_data[gv_flags].flux_radius)]
    plot_yrange=[max(cat_data[gv_flags].mag_auto)-2, min(cat_data[gv_flags].mag_auto)+1]

		if do_check_plot then begin
			if do_hd then begin
		    cgset, win2
		    cgplot, cat_data.flux_radius, cat_data.mag_auto, psym=cgsymcat('filled circle'), color='black', symsize=0.5, xrange=plot_xrange, yrange=plot_yrange, /xstyle, /ystyle, /window, xtitle='FWHM (pixels)', ytitle='Mag', title=(strsplit(cat_in_file,'/', /extract))[-1]
	 			cgplot, cat_data[gv_flags].flux_radius, cat_data[gv_flags].mag_auto, psym=cgsymcat('filled circle'), color='blue', symsize=0.5, /over, /addcmd
	    	cgplot, cat_data[gv_stars].flux_radius, cat_data[gv_stars].mag_auto, psym=cgsymcat('filled circle'), color='red', symsize=0.5, /over, /addcmd
	    	cgplot, fwhm_mean*[1,1], [0,1e6], line=2, color='blue', /addcmd, /over
	    	cgplot, fwhm_range[0]*[1,1], [0,1e6], line=1, color='blue', /addcmd, /over
	    	cgplot, fwhm_range[1]*[1,1], [0,1e6], line=1, color='blue', /addcmd, /over
			endif $
			else begin
				loadct, 12
				wset, 1
	   		plot, cat_data.flux_radius, cat_data.mag_auto, psym=3, symsize=0.5, xrange=plot_xrange, yrange=plot_yrange, /xstyle, /ystyle, xtitle='FWHM (pixels)', ytitle='Mag', title=(strsplit(cat_in_file,'/', /extract))[-1]
	    	oplot, cat_data[gv_flags].flux_radius, cat_data[gv_flags].mag_auto, psym=3, color=100, symsize=0.5
	    	oplot, cat_data[gv_stars].flux_radius, cat_data[gv_stars].mag_auto, psym=3, color=200, symsize=0.5

	    	oplot, fwhm_mean*[1,1], [0,1e6], line=2, color=100
	    	oplot, fwhm_range[0]*[1,1], [0,1e6], line=1, color=100
	    	oplot, fwhm_range[1]*[1,1], [0,1e6], line=1, color=100

			endelse
		endif

		print, 'Number of selected stars ', n_gv_stars
		print, 'Check WHY the stellar locus is so distorded in tile 1 - z band'
			
		fits_open, cat_in_file, cat_fcb
		fits_read, cat_fcb, cat_data0, cat_h0, exten=0
		fits_read, cat_fcb, cat_data1, cat_h1, exten=1
			
		cat_data2=make_array(cat_fcb.axis[0:1,2], value=0, /byte)	
		temp_lines=10000L
		temp_max=double(product(cat_fcb.axis[0:1,2])-1.)
		l=0L
		repeat begin
			print, 'Reading line ', strn(l*temp_lines), ' of ', strn(cat_fcb.axis[1,2])
			l1=double(l*temp_lines*cat_fcb.axis[0,2])
			l2=double((l+1)*temp_lines*cat_fcb.axis[0,2] - 1) < temp_max
			fits_read, cat_fcb, temp_data, cat_h2, first=l1, last=l2, exten_no=2
			cat_data2[*,l1/cat_fcb.axis[0,2]:(l2+1)/cat_fcb.axis[0,2]-1]=reform(temp_data, [cat_fcb.axis[0,2],(l2-l1+1)/cat_fcb.axis[0,2]])
			l++
		endrep until double(l*temp_lines*cat_fcb.axis[0,2]) GT temp_max
		temp_data=0
		fits_close, cat_fcb
			
		cat_data1=reform(cat_data1, [n_elements(cat_data1),1])
		cat_data2=cat_data2[*,[gv_stars]]
		fxaddpar, cat_h2, 'NAXIS2', (size(cat_data2, /dim))[1]
			
		if file_test(cat_out_file, /regular, /noexpand) then  file_delete, cat_out_file, /noexpand
		writefits, cat_out_file, cat_data0, cat_h0
		fits_open, cat_out_file, cat_fcb, /append
		fits_write, cat_fcb, cat_data1, cat_h1
		fits_write, cat_fcb, cat_data2, cat_h2
		fits_close, cat_fcb
		cat_data0=0 & cat_data1=0 & cat_data2=0

	endif

end
