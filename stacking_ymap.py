import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy import units as u
import healpy as hp
import os

NSIDE = 2048

UNWISE = '/mnt/d/data_large/unwise_sz/unWISE/catalog/blue_fullsky/'
catalog_list = os.listdir(UNWISE)

CMBY = '/mnt/d/data_large/unwise_sz/CMB_ymap/Planck/ymap/'
cmb_list = os.listdir(CMBY)

def read_catalog(catalogPath):
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
    
    base_catalog = {'ra': ra, 'dec': dec, 'l': l, 'b': b, 'flux_w1': flux_w1, 'flux_w2': flux_w2, 'pix': pix}
    return base_catalog