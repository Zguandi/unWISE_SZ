import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
import pymaster as nmt
from astropy.io import fits
import pandas as pd

sample1 = 1
sample2 = 1

if sample1 == 5 or sample1 == 6:
    print('Please make the CMB lensing map the second argument.')

#HEALPix map resolution
NSIDE = 2048

JOB = '/mnt/d/data_large/unwise_sz/'
PATHMAP = JOB + 'unWISE/' 
########################################
### READING RAW MAPS PATH
########################################

def read_path(sample):
    if sample == 5:
        map_name = PATHMAP + 'PLANCK_LENSING/COM_Lensing_4096_R3.00/MV/dat_klm.fits'
    elif sample == 6:
        map_name = PATHMAP + 'PLANCK_LENSING/COM_Lensing_Szdeproj_4096_R3.00/TT/dat_klm.fits' # tSZ deprojected 2018 lensing map (TT only)
    elif sample == 4:
        map_name = PATHMAP + 'current_maps/numcounts_map3_2048-r1_w2lt166flag.fits'    # Red but brighter W2 < 16.6
    elif sample == 7:
        map_name = PATHMAP + 'current_maps/numcounts_map3_2048-r1_w2lt162flag.fits'    # Red but brighter W2 < 16.2
    elif sample == 8:
        map_name = PATHMAP + 'current_maps/numcounts_map3_2048-r1_w2lt165flag.fits'    # Red but brighter W2 < 16.5
    elif sample == 9:
        map_name = PATHMAP + 'star_maps/gaiastarmap.fits'
    elif sample == 10:
        map_name = PATHMAP + 'blue_stars.fits'
    elif sample == 11:
        map_name = PATHMAP + 'green_stars.fits'
    elif sample == 12:
        map_name = PATHMAP + 'red_stars.fits'
    else:
        map_name = PATHMAP + 'blue/numcounts_map1_2048-r1-v2_flag.fits'
    return map_name


########################################
### READING MASK AND APODIZE
########################################

#mask_name = PATH + 'MASKS/sdss_mask.fits'     # Mask already apodized
lostmap = PATHMAP+"loss/unmaskedareafrac-flag.fits"
pl_mask_name = PATHMAP + 'mask/mask_unWISE_full_v10.fits'

print('Reading mask...')
mask = hp.read_map(PATHMAP+'mask/mask_unWISE_full_v10.fits')
#mask = hp.ud_grade(hp.read_map(mask_name, verbose=False),2048)*mask_default

print('Reading Planck mask...')
mask_pl = hp.read_map(pl_mask_name)
#hp.mollview(mask)
#pl.show()

lost = fits.open(lostmap)
mask_lost = lost[0].data

########################################
### MAKING MAPS
########################################

#   First field

PATHWEIGHTS = PATHMAP + 'weights/'
pathweight1 = PATHWEIGHTS + 'blue_w2_5sig_weights.fits'
pathweight2 = PATHWEIGHTS + 'blue_star_weights.fits'
weight1 = hp.read_map(pathweight1)
weight2 = hp.read_map(pathweight2)
weights = hp.ud_grade(weight1*weight2,2048)

def makemap(sample):
    if (sample == 5) or (sample == 6):
        print('Making lensing map')
        kappa_map_alm = hp.read_alm(read_path(sample))
        map = hp.alm2map(kappa_map_alm, NSIDE)
        # Note that namaster autmatically multiplies by the mask compute_full_master
        # I verify this in test_namaster.py (using Eiichiro Komatsu's MASTER code for mask deconvolution)
        # Can also check this by comparing nmt.compute_coupled_cell to hp.anafast
    elif sample >= 9:
        #numcounts_map = hp.read_map(read_path(sample))
        numcounts_map = fits.open(read_path(sample))[0].data
        masked_count = numcounts_map * mask
        mean_count = np.nansum(masked_count)/np.nansum(mask)
        masked_count_dn = numcounts_map / mean_count - 1.
        map = masked_count_dn
    else:
        # Converting the masked number counts to delta_n/n. Only consider unmasked regions!
        print('Making galaxy map ' + str(sample))
        numcounts_map = hp.read_map(read_path(sample), field=[0]) * weights
        # Correct for lower density in regions of high area lost due to stars or stellar masking
        numcounts_map = numcounts_map / mask_lost
        masked_count = numcounts_map * mask
        mean_count = np.nansum(masked_count) / np.nansum(mask)
        masked_count_dn = numcounts_map / mean_count - 1.
        
        map = masked_count_dn
        map[mask_lost == 0] = 0
        #std_map = np.sqrt( np.sum(map**2) / np.sum(mask) )
        #print std_map
        #hp.mollview(map)
        #pl.show()
    return map

print('Making map...')
galaxy_density = makemap(sample1)

########################################
### GALAXY NAMASTER FIELD
########################################

BIN = 50
print('Making galaxy field...')
b = nmt.NmtBin.from_nside_linear(NSIDE, BIN)
galaxy = nmt.NmtField(mask, [galaxy_density])

cl_gg = nmt.compute_full_master(galaxy, galaxy,b)

cl_gg_1 = cl_gg[0]
ll = b.get_effective_ells()


########################################
### SAVING
########################################

header = ['ell', 'cl_gg']
data = np.column_stack((ll, cl_gg_1))


PATHRESULTS = '/mnt/c/Users/gdzhao/projects/unwise_sz/ver3/result/nobeam/'

df = pd.DataFrame(data, columns=header)
df.to_csv(PATHRESULTS + 'namaster_comparison_gg.csv', index=False)
