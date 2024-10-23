import numpy as np
import healpy as hp
from astropy.io import fits
import pymaster as nmt

NSIDE = 2048
binnum = 50

DAT = '/mnt/d/data_large/unwise_sz/'
OUT = '/mnt/c/Users/gdzhao/projects/unwise_sz/Planck_temperature/'

# pathMaskHFI = DAT + 'Planck/temp/HFI_Mask_GalPlane-apo1_2048_R2.00.fits'

pathMaskPlanck = DAT+ 'Planck/temp/COM_Mask_CMB-common-Mask-apo1_2048_R3.00.fits'

pathMapPlanck = DAT + 'Planck/temp/COM_CMB_IQU-smica_2048_R3.00_full.fits'


print('Reading mask and map...')
mask_apodized = hp.read_map(pathMaskPlanck, field=0)
mask_apodized = hp.ud_grade(mask_apodized, nside_out=NSIDE)

map = hp.read_map(pathMapPlanck, field=0)
map = hp.ud_grade(map, nside_out=NSIDE)

b = nmt.NmtBin.from_nside_linear(NSIDE, binnum)

cols = [b.get_effective_ells()]
names = ['ell']

print('Making temperature field...')
temp = nmt.NmtField(mask_apodized, [map])

print('Computing cltt...')
cl_tt = nmt.compute_full_master(temp, temp, b)

cols.append(cl_tt[0])
names.append('cltt')

print('Saving gy and yy...')
np.savetxt(OUT + f'namaster_cltt.txt',np.column_stack(cols),header=' '.join(names))
print('Done!')