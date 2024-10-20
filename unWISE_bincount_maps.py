import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy import units as u
import healpy as hp
import os

DAT = '/mnt/d/data_large/unwise_sz/unWISE/catalog/blue_fullsky/'
catalog_list = os.listdir(DAT)

OUT = '/mnt/d/data_large/unwise_sz/unWISE/blue_bincounts/'

COMBINED = '/mnt/d/data_large/unwise_sz/unWISE/blue_bincounts_combined/'

NSIDE = 2048

def bincount_catalog(catalogPath):
    try:
        catalog = fits.open(catalogPath)[1].data
    except:
        print('Error reading catalog: ', catalogPath)
        return None
    
    ra = catalog['ra']
    dec = catalog['dec']
    
    coords = SkyCoord(ra=ra*u.degree, dec=dec*u.degree, frame='icrs')
    l, b = coords.galactic.l.value, coords.galactic.b.value
    
    flux_w1 = catalog['flux_w1']
    flux_w2 = catalog['flux_w2']
    
    pix = hp.ang2pix(NSIDE, l, b, lonlat=True)
    l_w1_map = np.bincount(pix, weights=flux_w1, minlength=hp.nside2npix(NSIDE))
    l_w2_map = np.bincount(pix, weights=flux_w2, minlength=hp.nside2npix(NSIDE))
    
    # release memory
    return l_w1_map, l_w2_map

def map_gen():
    for catalog in catalog_list:
        print('Processing: ', catalog)
        l_w1_map, l_w2_map = bincount_catalog(DAT + catalog)
        
        if l_w1_map is not None and l_w2_map is not None:
            hp.write_map(OUT + catalog.replace('.fits', '_w1_map.fits'), l_w1_map, overwrite=True)
            hp.write_map(OUT + catalog.replace('.fits', '_w2_map.fits'), l_w2_map, overwrite=True)
    
    return None

def numcount_catalog(catalogPath):
    try:
        catalog = fits.open(catalogPath)[1].data
    except:
        print('Error reading catalog: ', catalogPath)
        return None
    ra = catalog['ra']
    dec = catalog['dec']
    
    coords = SkyCoord(ra=ra*u.degree, dec=dec*u.degree, frame='icrs')
    l, b = coords.galactic.l.value, coords.galactic.b.value
    
    pix = hp.ang2pix(NSIDE, l, b, lonlat=True)
    numcount_map = np.bincount(pix, minlength=hp.nside2npix(NSIDE))
    
    # release memory
    return numcount_map

def countmap_gen():
    for catalog in catalog_list:
        print('Processing: ', catalog)
        numcount_map = numcount_catalog(DAT + catalog)
        
        if numcount_map is not None:
            hp.write_map(OUT + catalog.replace('.fits', '_numcount_map.fits'), numcount_map, overwrite=True)
    
    return None

def bincount_map_read(mapPath):
    try:
        map_data = hp.read_map(mapPath)
    except:
        print('Error reading map: ', mapPath)
        return None
    if 'w1' in mapPath:
        map_data = np.nan_to_num(map_data)
        type_data = 'w1'
    elif 'w2' in mapPath:
        map_data = np.nan_to_num(map_data)
        type_data = 'w2'
    elif 'numcount' in mapPath:
        type_data = 'numcount'
    return map_data, type_data

def combine_maps(type_data):
    map_list = os.listdir(OUT)
    numcount_map = np.zeros(hp.nside2npix(NSIDE))
    for map in map_list:
        if type_data in map:
            numcount_map += hp.read_map(OUT + map)
    hp.write_map(COMBINED + type_data + '_map.fits', numcount_map, overwrite=True)
    return None

if __name__ == '__main__':
    combine_maps('w2')
    # map_list = os.listdir(OUT)
    
    # map_data, type_data = bincount_map_read(OUT + map_list[1])
    
    # if map_data is not None:
    #     print(map_data.shape)
    #     print(type_data)