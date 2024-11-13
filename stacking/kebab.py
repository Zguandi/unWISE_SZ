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

class catalog_object:
    def __init__(self,ra = None,dec = None,flux_w1 = None,flux_w2 = None,nside = NSIDE):
        
        self.ra = ra
        self.dec = dec
        self.fix_coords()
        
        self.flux_w1 = flux_w1
        self.flux_w2 = flux_w2
        self.pix = hp.ang2pix(nside, self.l, self.b, lonlat=True)
    
    def __repr__(self):
        return 'ra: {}, dec: {}, l: {}, b: {}, flux_w1: {}, flux_w2: {}, pix: {}'.format(self.ra, self.dec, self.l, self.b, self.flux_w1, self.flux_w2, self.pix)
    
    def fix_coords(self):
        self.coords = SkyCoord(ra=self.ra*u.degree, dec=self.dec*u.degree, frame='icrs')
        self.l, self.b = self.coords.galactic.l.value, self.coords.galactic.b.value
        return None
        
    def get_neighbours(self,angle_distance_deg, nside=NSIDE):
        angular_distance = np.radians(angle_distance_deg)
        neighbour_pix = hp.query_disc(nside, hp.ang2vec(self.l, self.b, lonlat=True), angular_distance)
        return neighbour_pix

    def get_ymap_vals(self,angle_distance_deg,ymap=None,nside=NSIDE):
        if ymap is None:
            return None
        else:
            neighbour_pixels = self.get_neighbours(angle_distance_deg, nside=nside)
            return ymap[neighbour_pixels]
    
    def get_pixel_vals_hist(self,angle_distance_deg,ymap=None,nside=NSIDE,nbins=100):
        pixels = self.get_neighbours(angle_distance_deg,nside=nside)
        values = self.get_ymap_vals(angle_distance_deg,ymap=ymap,nside=nside)
        centre_pixel = self.pix
        return self.cast_values_2d(pixels,centre_pixel,values,nbins=nbins)
    
    def cast_values_2d(self,pixels,centre_pixel,values,nbins = 100):
        vec_centre = hp.pix2vec(NSIDE,centre_pixel)
        vecs = hp.pix2vec(NSIDE,pixels)

        # set x axis perpendicular to the centre pixel
        x_axis = np.cross(vec_centre,[0,0,1])
        x_axis = x_axis/np.linalg.norm(x_axis)
        
        # set y axis perpendicular to the centre pixel and the x axis
        y_axis = np.cross(vec_centre,x_axis)
        y_axis = y_axis/np.linalg.norm(y_axis)
        
        # project the vectors onto the x and y axes
        x_projections = np.dot(x_axis,vecs)
        y_projections = np.dot(y_axis,vecs)
        
        # create a 2d histogram
        hist, xedges, yedges = np.histogram2d(x_projections, y_projections, bins=nbins, weights=values)
        
        # get the bin centers
        xcenters = (xedges[:-1] + xedges[1:]) / 2
        ycenters = (yedges[:-1] + yedges[1:]) / 2
        
        return xcenters, ycenters, hist.T
    


def read_catalog(catalogPath,ntest = None):
    try:
        catalog = fits.open(catalogPath)[1].data
    except:
        print('Error reading catalog: ', catalogPath)
        return None
    
    ra = catalog['ra']
    dec = catalog['dec']
    
    flux_w1 = catalog['flux_w1']
    flux_w2 = catalog['flux_w2']
    
    # yield a iterable object containing the catalog objects
    # if ntest is interval (a,b) then it will return objects from a to b
    # if ntest is a number it will return the first ntest objects
    
    if ntest is None:
        ntest = range(len(ra))
    elif isinstance(ntest,list):
        pass
    elif isinstance(ntest,tuple):
        a,b = ntest
        ntest = range(a,b)
    elif isinstance(ntest,int):
        ntest = range(ntest)
    
    for i in ntest:
        yield catalog_object(ra[i], dec[i], flux_w1[i], flux_w2[i])
        
        
if __name__ == '__main__':
    catalogPath = UNWISE + catalog_list[0]
    ymapPath = CMBY + cmb_list[0]
    
    catalog = read_catalog(catalogPath,ntest=(0,10))
    ymap = hp.read_map(ymapPath)
    
    for obj in catalog:
        print(obj)
        print(obj.get_pixel_vals_hist(1,ymap=ymap,nbins=100))
        print('\n\n')
        
    print('done')