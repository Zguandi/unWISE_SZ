# This is the apodization for planck galaxy mask for Temperature map

DAT = '/mnt/d/data_large/unwise_sz/'

pathMaskPlanck = DAT+ 'Planck/temp/COM_Mask_CMB-common-Mask-Int_2048_R3.00.fits'

import numpy as np
import healpy as hp
from astropy.io import fits
import pymaster as nmt

aposcale = 1.0
OUT = DAT + 'Planck/temp/COM_Mask_CMB-common-Mask-apo1_2048_R3.00.fits'

mask = hp.read_map(pathMaskPlanck, field=0).astype(np.float64)
mask_apodized = nmt.mask_apodization(mask, aposcale, apotype='C1')

hp.write_map(OUT, mask_apodized, overwrite=True)