import healpy as hp
from assets import galaxy_map_read as gmr
from pixell import enmap, utils, reproject, wcsutils
from assets import deprojection_index as di
from astropy.io import fits
import numpy as np
import os

# "D:\data_large\unwise_sz\unWISE\mask"
OUTDIR = '/mnt/d/data_large/unwise_sz/unWISE/mask/'

def default_wcs():
    codex = di.get_ymap_index_act()
    index = codex[0]
    
    with fits.open(index.path) as hdul:
        # data = np.array(hdul[0].data,dtype=np.float32)
        header = hdul[0].header
    wcs = wcsutils.WCS(header)
    # shape = data.shape
    return wcs

def mask_healpix_to_car(outfile:str):
    if not os.path.exists(OUTDIR):
        os.makedirs(OUTDIR)
    unwisesample = 'blue'
    mask = gmr.readmask(unwisesample)
    
    wcs = default_wcs()
    mask = reproject.healpix2map(mask, wcs = wcs, shape = (10320,43200),method='spline',verbose = True)
    
    outpath = OUTDIR + outfile
    hp.write_map(outpath,mask)

if __name__ == '__main__':
    mask_healpix_to_car('mask_unWISE_full_v10_car.fits')