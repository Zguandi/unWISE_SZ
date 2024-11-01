# This is the apodization for act SZ times galaxy mask galaxy mask

DAT = '/mnt/d/data_large/unwise_sz/'

import numpy as np
import healpy as hp
from astropy.io import fits
import pymaster as nmt

path_act_mask = DAT+'unwiseact/act_mask_composite/healpix_act_mask_nside2048.fits'
path_unwise_mask = DAT+'unwiseact/unwise_mask_composite/healpix_unwise_mask_nside2048.fits'

# mask_act = hp.read_map(path_act_mask, field=0).astype(np.float32)

mask_unwise = hp.read_map(path_unwise_mask, field=0).astype(np.float32)
# mask_composite = mask_act * mask_unwise

aposcale = 1.5
OUT = DAT+'unwiseact/unwise_mask_composite/healpix_unwise_mask_nside2048_apo1_5.fits'

mask_apodized = nmt.mask_apodization(mask_unwise, aposcale, apotype='C1')

hp.write_map(OUT, mask_apodized, overwrite=True)