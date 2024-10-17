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
    return map_data, type_data

if __name__ == '__main__':    
    map_list = os.listdir(OUT)
    
    map_data, type_data = bincount_map_read(OUT + map_list[1])
    
    if map_data is not None:
        print(map_data.shape)
        print(type_data)