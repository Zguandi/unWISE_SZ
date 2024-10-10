import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
import pymaster as nmt
from astropy.io import fits
import pandas as pd

#HEALPix map resolution
NSIDE = 2048

DAT = '/mnt/d/data_large/unwise_sz/'
PATHMAP = DAT + 'unWISE/' 

unWISE_mask_name = PATHMAP + 'mask/mask_unWISE_full_v10.fits'
pl_mask_name = PATHMAP + 'mask/mask_unWISE_full_v10.fits'
lostmap_name = PATHMAP+"loss/unmaskedareafrac-flag.fits"

print('Reading mask...')
mask = hp.read_map(unWISE_mask_name)

print('Reading Planck mask...')
mask_pl = hp.read_map(pl_mask_name)

lost = fits.open(lostmap_name)
mask_lost = (lost[0].data != 0)

complete_mask = mask * mask_pl * mask_lost

hp.mollview(complete_mask, title='Complete mask')
plt.savefig('/mnt/c/Users/gdzhao/projects/unwise_sz/ver3/result/images/complete_mask.png')
plt.close()

f_sky = np.sum(complete_mask) / len(complete_mask)

print(f'Sky fraction: {f_sky}')

# At default, the sky fraction is 0.572050134340922