
import numpy as np
import healpy as hp
from astropy.io import fits

# sample1 = 1
# sample2 = 1

# if sample1 == 5 or sample1 == 6:
#     print('Please make the CMB lensing map the second argument.')

#HEALPix map resolution
NSIDE = 2048

DAT = '/mnt/d/data_large/unwise_sz/'
PATHMAP = DAT + 'unWISE/'

PATHWEIGHTS = PATHMAP + 'weights/'
pathweight1 = PATHWEIGHTS + 'blue_w2_5sig_weights.fits'
pathweight2 = PATHWEIGHTS + 'blue_star_weights.fits'
pathmask = PATHMAP+'mask/mask_unWISE_full_v10.fits'
lostmap = PATHMAP+"loss/unmaskedareafrac-flag.fits"

########################################
### READING RAW MAPS PATH
########################################

def read_path(sample):
    # if sample == 5:
    #     map_name = PATHMAP + 'PLANCK_LENSING/COM_Lensing_4096_R3.00/MV/dat_klm.fits'
    # elif sample == 6:
    #     map_name = PATHMAP + 'PLANCK_LENSING/COM_Lensing_Szdeproj_4096_R3.00/TT/dat_klm.fits' # tSZ deprojected 2018 lensing map (TT only)
    # elif sample == 4:
    #     map_name = PATHMAP + 'current_maps/numcounts_map3_2048-r1_w2lt166flag.fits'    # Red but brighter W2 < 16.6
    # elif sample == 7:
    #     map_name = PATHMAP + 'current_maps/numcounts_map3_2048-r1_w2lt162flag.fits'    # Red but brighter W2 < 16.2
    # elif sample == 8:
    #     map_name = PATHMAP + 'current_maps/numcounts_map3_2048-r1_w2lt165flag.fits'    # Red but brighter W2 < 16.5
    # elif sample == 9:
    #     map_name = PATHMAP + 'star_maps/gaiastarmap.fits'
    # elif sample == 10:
    #     map_name = PATHMAP + 'blue_stars.fits'
    # elif sample == 11:
    #     map_name = PATHMAP + 'green_stars.fits'
    # elif sample == 12:
    #     map_name = PATHMAP + 'red_stars.fits'
    # else:
    map_name = PATHMAP + 'blue/numcounts_map1_2048-r1-v2_flag.fits'
    return map_name

def makemap(sample):
        # if (sample == 5) or (sample == 6):
        #     print('Making lensing map')
        #     kappa_map_alm = hp.read_alm(read_path(sample))
        #     map = hp.alm2map(kappa_map_alm, NSIDE)
        
        # Note that namaster autmatically multiplies by the mask compute_full_master
        # I verify this in test_namaster.py (using Eiichiro Komatsu's MASTER code for mask deconvolution)
        # Can also check this by comparing nmt.compute_coupled_cell to hp.anafast
        
        # elif sample >= 9:
        #     #numcounts_map = hp.read_map(read_path(sample))
        #     numcounts_map = fits.open(read_path(sample))[0].data
        #     masked_count = numcounts_map * mask
        #     mean_count = np.nansum(masked_count)/np.nansum(mask)
        #     masked_count_dn = numcounts_map / mean_count - 1.
        #     map = masked_count_dn
        # else:
        
    ########################################
    print('Reading weights...') 
    weight1 = hp.read_map(pathweight1)
    weight2 = hp.read_map(pathweight2)
    weights = hp.ud_grade(weight1*weight2,2048)
    ########################################
    
    
    ########################################
    print('Reading mask...')
    mask = hp.read_map(pathmask)
    lost = fits.open(lostmap)
    mask_lost = lost[0].data
    ########################################
    
    # omit divide by zero errors
    
    np.seterr(divide='ignore', invalid='ignore')
    ########################################
    # Converting the masked number counts to delta_n/n. Only consider unmasked regions!
    print('Making galaxy map ' + str(sample))
    numcounts_map = hp.read_map(read_path(sample), field=[0]) * weights
    # Correct for lower density in regions of high area lost due to stars or stellar masking
    ########################################
    
    
    ########################################
    numcounts_map = numcounts_map / mask_lost
    masked_count = numcounts_map * mask
    mean_count = np.nansum(masked_count) / np.nansum(mask)
    masked_count_dn = numcounts_map / mean_count - 1.
    ########################################
    
    
    map = masked_count_dn
    map[mask_lost == 0] = 0
    
    
    ########################################
    #std_map = np.sqrt( np.sum(map**2) / np.sum(mask) )
    #print std_map
    #hp.mollview(map)
    #pl.show()
    ########################################
        
    return map

def readmask():
    ########################################
    print('Reading mask...')
    mask = hp.read_map(pathmask)
    lost = fits.open(lostmap)
    mask_lost = lost[0].data
    invalid = mask_lost == 0
    
    galmask = mask * invalid
    ########################################
    return galmask
    
if __name__ == "__main__":
    readmask()