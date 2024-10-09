import numpy as np
import pymaster as nmt
import os
import matplotlib.pyplot as plt
import healpy as hp
from astropy.io import fits

y_mask = hp.read_map('/home/zgd/cmb_y/ver1/Planck/mask/COM_CompMap_Compton-SZMap-masks_2048_R2.01.fits', field=1)
y_mask_LFI = hp.read_map('/home/zgd/cmb_y/ver1/Planck/mask/LFI_inpainting_bool.fits', field=0)
y_mask_HFI = hp.read_map('/home/zgd/cmb_y/ver1/Planck/mask/HFI_inpainting_bool.fits', field=0)

aposize = 1.
y_mask_HFI_apodized = nmt.mask_apodization(y_mask_HFI, aposize=aposize, apotype='Smooth')
y_mask_LFI_apodized = nmt.mask_apodization(y_mask_LFI, aposize=aposize, apotype='Smooth')

y_mask = y_mask_HFI_apodized*y_mask_LFI_apodized

im = hp.mollview(y_mask)
im.savefig('/home/zgd/cmb_y/ver2/output/ymask.png')


