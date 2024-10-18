import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
import pymaster as nmt
from astropy.io import fits

# custom package
from assets import deprojection_index
from assets import make_galaxy_map
#HEALPix map resolution
NSIDE = 2048

DAT = '/mnt/d/data_large/unwise_sz/'

pathmaskact = DAT + 'ACT/mask/wide_mask_GAL070_apod_1.50_deg_wExtended.fits'

galmask = make_galaxy_map.readmask()


with fits.open(pathmaskact) as hdul:
    hdul.info()
    data = hdul[0].data
    
    print(data.shape)

print(galmask.shape)
# print(actmask.shape)
OUTPATH = '/mnt/c/Users/gdzhao/projects/unwise_sz/ver3/result/nobeam/'

