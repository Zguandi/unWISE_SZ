import os
import matplotlib.pyplot as plt
import healpy as hp
import numpy as np
import pymaster as nmt
from scipy.interpolate import interp1d
from astropy.io import fits

apodized = True
UNWISE = '/mnt/d/data_large/unwise_sz/unWISE/'
PLANCK = '/mnt/d/data_large/unwise_sz/Planck/'
OUT = '/mnt/d/data_large/unwise_sz/unwiseplanck/combined_mask/'

y_mask = hp.read_map(PLANCK+'mask/COM_CompMap_Compton-SZMap-masks_2048_R2.01.fits', field=1)
y_mask_LFI = hp.read_map(PLANCK+'mask/LFI_inpainting_bool.fits', field=0)
y_mask_HFI = hp.read_map(PLANCK+'mask/HFI_inpainting_bool.fits', field=0)

unwise_mask = hp.read_map(UNWISE+'mask/mask_unWISE_full_v10.fits', field=0)
lost = fits.open(UNWISE+"loss/unmaskedareafrac-flag.fits")
map_lost = lost[0].data

unwise_mask[map_lost < 0.1] = 0.0

if apodized:
    print('read complete, apodizing...')

    aposcale = 1.0
    # apodize
    y_mask_LFI_apodized = nmt.mask_apodization(y_mask_LFI, aposcale, apotype='C1')
    y_mask_HFI_apodized = nmt.mask_apodization(y_mask_HFI, aposcale, apotype='C1')

    y_mask = y_mask*y_mask_LFI_apodized*y_mask_HFI_apodized
    
    g_mask = nmt.mask_apodization(unwise_mask, aposcale, apotype='C1')
    total_mask = y_mask*g_mask
    
    hp.write_map(OUT+'total_mask_apo.fits', total_mask, overwrite=True)
    
    mask_binary = (total_mask > 0.3).astype(bool)
    
    hp.write_map(OUT+'total_mask_apo_bool.fits', mask_binary, overwrite=True)
    
else:
    print('read complete, no apodization...')
    y_mask = y_mask*y_mask_LFI*y_mask_HFI
    
    total_mask = y_mask*unwise_mask
    
    hp.write_map(OUT+'total_mask_noapo.fits', total_mask, overwrite=True)
    
    mask_binary = (total_mask > 0.3).astype(bool)
    
    hp.write_map(OUT+'total_mask_noapo_bool.fits', mask_binary, overwrite=True)


