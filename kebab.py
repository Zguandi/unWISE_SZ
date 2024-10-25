# stacking of CMB y-maps around regions close to unWISE objects.

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

def test_catalog(catalog,ntest):
    
    indices = np.random.choice(len(catalog['ra']),ntest,replace=False)
    
    ra = catalog['ra'][indices]
    dec = catalog['dec'][indices]
    l = catalog['l'][indices]
    b = catalog['b'][indices]
    flux_w1 = catalog['flux_w1'][indices]
    flux_w2 = catalog['flux_w2'][indices]
    pix = catalog['pix'][indices]
    
    return {'ra': ra, 'dec': dec, 'l': l, 'b': b, 'flux_w1': flux_w1, 'flux_w2': flux_w2, 'pix': pix}

def get_pixel_around_centre(l,b,nside=NSIDE,radius_deg = 1.0):
    vec_centre = hp.ang2vec(l,b,lonlat=True)
    radius_rad = np.radians(radius_deg)
    pix_neighbors = hp.query_disc(nside, vec_centre, radius_rad)
    vec_neighbors = np.array(hp.pix2vec(nside, pix_neighbors)).T
    
    # define local 2D coordinate system
    x_axis = np.cross(vec_centre, [0,0,1])
    y_axis = np.cross(vec_centre, x_axis)
    x_axis /= np.linalg.norm(x_axis)
    y_axis /= np.linalg.norm(y_axis)

    # project the neighbors onto the local 2D coordinate system
    proj_x = np.dot(vec_neighbors, x_axis)
    proj_y = np.dot(vec_neighbors, y_axis)
    
    proj_catalog = {'x': proj_x, 'y': proj_y, 'pix': pix_neighbors}
    return proj_catalog

def get_yvalue(proj_catalogs,ymap):
    for proj_catalog in proj_catalogs:
        ymap_pix = proj_catalog['pix']
        ymap_values = ymap[ymap_pix]
        proj_catalog['ymap'] = ymap_values
    return proj_catalog
    
def stack_yvalue(proj_catalog_list,nbins = 100,boundary_angle = 1.0):
    counts = np.zeros((nbins,nbins),dtype=np.int64)
    ystack = np.zeros((nbins,nbins),dtype=np.float64)
    
    boundary_angle = np.radians(boundary_angle)
    x_min, x_max = -np.sin(boundary_angle), np.sin(boundary_angle)
    y_min, y_max = -np.sin(boundary_angle), np.sin(boundary_angle)
    
    x_bins = np.linspace(x_min,x_max,nbins+1)
    y_bins = np.linspace(y_min,y_max,nbins+1)
    
    for proj_catalog in proj_catalog_list:
        x = proj_catalog['x']
        y = proj_catalog['y']
        ymap = proj_catalog['ymap']
        
        x_idx = np.digitize(x,x_bins) - 1
        y_idx = np.digitize(y,y_bins) - 1
        
        for i in range(len(x)):
            counts[x_idx[i],y_idx[i]] += 1
            ystack[x_idx[i],y_idx[i]] += ymap[i]
    
    ystack /= counts
    return counts,ystack

if __name__ == '__main__':
    catalog_base = read_catalog(UNWISE + catalog_list[0])
    catalog_test = test_catalog(catalog_base,10)
    
    print(catalog_test['ra'], catalog_test['dec'])