import os, fnmatch, glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.cluster import MeanShift, estimate_bandwidth
from astropy.io import fits

survey_file = 'data/MS_VISTA_P98.csv'
im_folder= '/Volumes/RAID/Magellanic_stream/VISTA/P98/raw_data'
im_pattern = '*st.fit'

header_keys = ['HIERARCH ESO OBS PROG ID', 'HIERARCH ESO OBS NAME', 'HIERARCH ESO INS FILT1 NAME', \
	'HIERARCH ESO DET DIT', 'HIERARCH ESO DET NDIT', \
	'HIERARCH ESO TPL ID', 'HIERARCH ESO TPL NEXP', 'HIERARCH ESO TPL EXPNO', \
	'JITTR_ID', 'NJITTER', 'JITTER_I', 'ESOGRADE', 'RA', 'DEC']


if not(os.path.exists(survey_file)):
	im_files = glob.glob( os.path.join(im_folder,  im_pattern))

	survey_list=[]
	for im_file in im_files:
		print 'Processing file ', os.path.split(im_file)[1]
		hdulist = fits.open(im_file)
		im_header = hdulist[0].header

		im_fwhm = hdulist[1].header['SEEING'] * 0.34 # FWHM in arcsec

		result = {}
		for key in header_keys:
			result[key] = im_header[key]
		result['FWHM'] = im_fwhm
		result['Filename'] = im_file

		survey_list.append(result)

	survey_df = pd.DataFrame(survey_list, columns=header_keys+['FWHM','Filename'])
	survey_df.to_csv(survey_file, index_label='Index')
else:
	survey_df = pd.read_csv(survey_file, header=0, index_col=0, comment='#')

print 'Columns in survey_df: ', list(survey_df.columns.values)

survey_good = survey_df[ ((survey_df['ESOGRADE'] == 'A') | (survey_df['ESOGRADE'] == 'B')) \
	& (survey_df['FWHM'] < 1.0) ].reset_index(drop=True)

# All OBs
fig, ax = plt.subplots(1, 1, figsize=(8,6))
plt.plot(survey_df['RA'], survey_df['DEC'], linestyle='None', marker='+', color='red')

ax.set_xlabel('R.A. (deg)')
ax.set_ylabel('Dec. (deg)')
fig.savefig('figures/MS_VISTA_P98_ALL.pdf')

ax.set_xlim([351.,352.])
ax.set_ylim([-53.9,-52.9])
fig.savefig('figures/MS_VISTA_P98_ALL_zoom.pdf')

X = survey_df[['RA','DEC']].values
ms = MeanShift(bandwidth=0.2, bin_seeding=True)
ms.fit(X)
labels = ms.labels_
cluster_centers = ms.cluster_centers_

gv_sort = np.argsort(cluster_centers[:,1])
for i in gv_sort:
	center=cluster_centers[i,:]
	gv = (labels == i)
	print 'Center: ', center
	print 'Number of pointings: ', np.sum(gv)
	print 'Total exposure time in seconds per pixel: {:.1f}'.format(np.sum(survey_df[gv]['HIERARCH ESO DET DIT']*survey_df[gv]['HIERARCH ESO DET NDIT'])/6*2)


# Only OBs with grade A and B, with FWHM better than 1.0 arsec
fig, ax = plt.subplots(1, 1, figsize=(8,6))
plt.plot(survey_good['RA'], survey_good['DEC'], linestyle='None', marker='+', color='red')

ax.set_xlabel('R.A. (deg)')
ax.set_ylabel('Dec. (deg)')
fig.savefig('figures/MS_VISTA_P98_gradeAB_FWHM1.0.pdf')

ax.set_xlim([351.,352.])
ax.set_ylim([-53.9,-52.9])
fig.savefig('figures/MS_VISTA_P98_gradeAB_FWHM1.0_zoom.pdf')

ax.set_xlim([351.36,351.40])
ax.set_ylim([-53.5,-53.45])
fig.savefig('figures/MS_VISTA_P98_gradeAB_FWHM1.0_zoom2.pdf')

X = survey_good[['RA','DEC']].values
ms = MeanShift(bandwidth=0.2, bin_seeding=True)
ms.fit(X)
labels = ms.labels_
cluster_centers = ms.cluster_centers_

gv_sort = np.argsort(cluster_centers[:,1])
for i in gv_sort:
	center=cluster_centers[i,:]
	gv = (labels == i)
	print 'Center: ', center
	print 'Number of pointings: ', np.sum(gv)
	print 'Total exposure time in seconds per pixel: {:.1f}'.format(np.sum(survey_good[gv]['HIERARCH ESO DET DIT']*survey_good[gv]['HIERARCH ESO DET NDIT'])/6*2)
