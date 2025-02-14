import sys
import os

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
import pymaster as nmt
from astropy.io import fits

# cuntom package
from assets import deprojection_index
# from assets import make_galaxy_map

sample1 = 1
sample2 = 1

if sample1 == 5 or sample1 == 6:
    print('Please make the CMB lensing map the second argument.')

#HEALPix map resolution
NSIDE = 2048

DAT = '/mnt/d/data_large/unwise_sz/'
PATHMAP = DAT + 'unWISE/' 

PATHWEIGHTS = PATHMAP + 'weights/'

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
### Y_MAP HANDLING
########################################

NSIDE = 2048
lmax = 4096
bin = 50

y_mask = np.load(DAT+'Planck/mask/y_mask.npy')

comparison_group = 2

ymap_name_list, ymap_path_list = deprojection_index.get_ymap_index_planck(comparison_group)

########################################
### GALAXY NAMASTER FIELD
########################################

print('Making galaxy field...')
b = nmt.NmtBin.from_nside_linear(NSIDE, bin)
galaxy = nmt.NmtField(mask, [galaxy_density])

########################################
### Y NAMASTER FIELD
########################################

b = nmt.NmtBin.from_nside_linear(NSIDE, bin)
ell = np.arange(3*NSIDE)
FWHM = 10.0 # arcmin
sigma = FWHM/2.35
beam = np.exp(-ell*(ell+1)*sigma**2)
# galaxy = nmt.NmtField(mask, [galaxy_density])
# print('g initialized')

cols = [b.get_effective_ells()]
names = ['ell']

for i in range(len(ymap_name_list)):
    print('Initializing y for '+ymap_name_list[i]+'...',end='\r')
    y_CIB  = nmt.NmtField(y_mask, [hp.read_map(ymap_path_list[i])])
    print('y initialization complete for '+ymap_name_list[i])
    
    print('computing gy for '+ymap_name_list[i]+'...',end='\r')
    cl_gy = nmt.compute_full_master(galaxy, y_CIB, b)
    print('namaster gy compelete for '+ymap_name_list[i])

    print('computing yy for '+ymap_name_list[i]+'...',end='\r')
    cl_yy = nmt.compute_full_master(y_CIB, y_CIB, b)
    print('namaster yy compelete for '+ymap_name_list[i])

    names.append(ymap_name_list[i]+'_gy')
    cols.append(cl_gy[0])

    names.append(ymap_name_list[i]+'_yy')
    cols.append(cl_yy[0])


OUTPATH = '/mnt/c/Users/gdzhao/projects/unwise_sz/ver3/result/nobeam/'
print('Saving gy and yy...')
np.savetxt(OUTPATH + f'namaster_comparison{comparison_group}.txt',np.column_stack(cols),header=' '.join(names))
print('Done!')