
import numpy as np
import healpy as hp
from astropy.io import fits
import matplotlib.pyplot as plt

NSIDE = 2048
DAT = '/mnt/d/data_large/unwise_sz/unWISE/'


########################################
### READING RAW MAPS PATH
########################################


def read_path(sample: str):
    match sample:
        case "blue": # blue sample
            map_path = DAT + 'blue/map_healpix/numcounts_map1_2048-r1-v2_flag.fits'
            mask_path = DAT + 'blue/mask_healpix/mask_unWISE_full_v10.fits'
            loss_path = DAT + 'blue/loss_healpix/unmaskedareafrac-flag.fits'
        # case "bright": # bright sample
        #     map_path = DAT + 'bright/map_healpix/numcounts_map0_2048-r1_w2lt166flag.fits'
        #     mask_path = DAT + 'bright/mask_healpix/mask_unWISE_full_v10.fits'
        case "midz": # midz sample
            map_path = DAT + 'midz/map_healpix/midz_fullsky_catalog_agglomerative_clustering_fix_healpix2048.fits'
            mask_path = DAT + 'midz/mask_healpix/mask_unWISE_full_v10.fits'
            loss_path = DAT + 'midz/loss_healpix/neo8-areamask-midz-2048.fits'
            
        case "lowz": # lowz sample
            map_path = DAT + 'lowz/map_healpix/lowz_fullsky_catalog_agglomerative_clustering_fix_healpix2048.fits'
            mask_path = DAT + 'lowz/mask_healpix/mask_unWISE_full_v10.fits'
            loss_path = DAT + 'lowz/loss_healpix/neo8-areamask-lowz-2048.fits'
    
    #"D:\data_large\unwise_sz\unWISE\weights\blue_star_weights.fits"
    path_weight1 = DAT + 'weights/blue_star_weights.fits'
    # "D:\data_large\unwise_sz\unWISE\weights\blue_w2_5sig_weights.fits"
    path_weight2 = DAT + 'weights/blue_w2_5sig_weights.fits'
    
    return map_path, mask_path, loss_path, path_weight1, path_weight2


def makemap_healpix(sample: str):

    print('Galaxy sample generating started for ' + sample)
    map_path, mask_path, loss_path, path_weight1, path_weight2 = read_path(sample)
    
    ########################################
    
    print('\t Reading weights...',end = "\r") 
    weight1 = hp.read_map(path_weight1)
    weight2 = hp.read_map(path_weight2)
    weights = hp.ud_grade(weight1*weight2,2048)

    ########################################
    
    # Converting the masked number counts to delta_n/n. Only consider unmasked regions!
    print('\t Making galaxy map ' + sample,end = "\r")
    numcounts_map = hp.read_map(map_path, field=[0])
    numcounts_map = numcounts_map * weights
    # Correct for lower density in regions of high area lost due to stars or stellar masking
    
    ########################################
    
    # omit divide by zero errors
    np.seterr(divide='ignore', invalid='ignore')
    
    ########################################
    
    print('\t Reading mask and loss...',end = "\r")
    mask = hp.read_map(mask_path)
    lost = fits.open(loss_path)
    
    if sample == "blue":
        loss_map = lost[0].data
        
        print("\t loss_map generated with shape:",loss_map.shape,end = "\r")
    elif sample == "midz" or sample == "lowz":
        loss_data = lost[1].data
        loss_column = loss_data.field(0)
        loss_map = np.concatenate(loss_column)
        
        # set nan values of loss to 1.0
        loss_map[np.isnan(loss_map)] = 1.0
        
        # print("\t loss_map generated with shape:",loss_map.shape,end = "\r")
    else:
        raise ValueError("Invalid sample")
    
    lost.close()
    
    ########################################
    
    print("\t compositing loss and normalizing...",end = "\r")
    numcounts_map = numcounts_map / loss_map
    masked_count = numcounts_map * mask
    mean_count = np.nansum(masked_count) / np.nansum(mask)
    masked_count_dn = numcounts_map / mean_count - 1.
    
    ########################################

    map = masked_count_dn
    map[loss_map == 0] = 0
        
    ########################################
    #std_map = np.sqrt( np.sum(map**2) / np.sum(mask) )
    #print std_map
    #hp.mollview(map)
    #pl.show()
    ########################################
    print('Galaxy sample generating finished for ' + sample)
    return map

def readmask(sample: str):
    
    map_path, mask_path, loss_path, path_weight1, path_weight2 = read_path(sample)
    
    print('Reading mask...')
    mask = hp.read_map(mask_path)
    
    # the mask is apodized already, and ready to use
    
    return mask


if __name__ == "__main__":
    # mask = readmask()
    # print(np.min(mask),np.max(mask))
    # print(np.sum(mask))
    # blue_counts_map = makemap()
    # print(np.min(blue_counts_map),np.max(blue_counts_map))
    # print(np.sum(blue_counts_map))
    
    sample = 'lowz'
    
    map = makemap_healpix(sample)
    # mask = readmask(sample)
    # print(np.min(mask),np.max(mask))
    # map = makemap_healpix(sample)
    # print(np.min(map),np.max(map))