import sys
import os

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from assets import act_car as act
from assets import make_galaxy_map as mgm
import healpy as hp
import numpy as np
from astropy.io import fits
import pymaster as nmt
import matplotlib.pyplot as plt
from pixell import enmap

import gc

NSIDE = 2048
binwidth = 50
outpath = '/mnt/c/Users/gdzhao/projects/unwise_sz/unwiseact_results/car/'

# deprotype = 'cib_cibdBeta_cibdT'
# deprotype = 'cib_cibdBeta'
# deprotype = 'cib'

selected_maps = act.get_ymap_index_act_selected(deprotype='cib_cibdBeta',beta_range=[1.21,1.4],T_range=[10.0,12.0],verbose=True)
# map_index = selected_maps[0]

b = nmt.NmtBin.from_nside_linear(NSIDE, binwidth)
ells = b.get_effective_ells()

# unwise_map = mgm.makemap()
# unwise_mask = mgm.readmask()

# g_field = nmt.NmtField(unwise_mask, [unwise_map])

# print('field created, with shape:',g_field.get_maps()[0].shape)
# del unwise_map, unwise_mask
# gc.collect()

# cl_gg = nmt.compute_full_master(g_field, g_field, b)
# print('namaster gg complete')

# np.savetxt(outpath+'gg_results.txt',np.array([ells,cl_gg[0]]).T,header='ell cl_gg')

# ymap, wcs = selected_maps[0].read_map_to_array()
# print('map read complete with shape:',ymap.shape)
# print('wcs:',wcs)

# mask, wcs_n = act.read_composite_mask_to_array()
# print('mask read complete with shape:',mask.shape)
# print('wcs:',wcs_n)

##########################
# DOWNGRADE for MEM LIMIT#
##########################

for map_index in selected_maps:
    print(map_index)
    print('downgrading ymaps')
    ymap_enmap = map_index.read_map_to_enmap()
    ymap_enmap = enmap.downgrade(ymap_enmap, 2)
    wcs = ymap_enmap.wcs
    ymap = np.array(ymap_enmap)

    del ymap_enmap
    gc.collect()

    print('downgrading mask')
    mask_enmap = act.read_composite_mask_to_enmap()
    mask_enmap = enmap.downgrade(mask_enmap, 2)
    wcs_n = mask_enmap.wcs
    mask = np.array(mask_enmap)

    del mask_enmap
    gc.collect()

    lls, beam = act.read_beam()
    y_field = nmt.NmtField(mask, [ymap], wcs = wcs, n_iter=0,beam = beam)

    print('downgraded y field created')
    del ymap, mask
    gc.collect()

    model_beta = map_index.beta
    model_t = map_index.T
    model_deprotype = map_index.deprotype

    yy_filename =  f"{outpath}yy_{model_deprotype}_T{model_t}_beta{model_beta}.txt"

    print('y field created, with shape:',y_field.get_maps()[0].shape)

    cl_yy = nmt.compute_full_master(y_field, y_field, b)
    print('namaster yy complete')

    np.savetxt(yy_filename,np.array([ells,cl_yy[0]]).T,header='ell cl_yy')

# cl_gy = nmt.compute_full_master(g_field, y_field, b)
# print('namaster gy complete')

# np.savetxt(outpath+'cibdbeta_gy_results.txt',np.array([ells,cl_gy[0]]).T,header='ell cl_gy')
# plt.figure()
# plt.title("Map from NmtField")
# plt.imshow(y_field.get_maps().reshape([ymap.shape[0], -1]),
#            interpolation='nearest', origin='lower')
# plt.savefig('y_field.png')