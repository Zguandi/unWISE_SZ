
DAT = '/mnt/d/data_large/unwise_sz/'
BEAM = '/mnt/d/data_large/unwise_sz/ACT/act_beam/effective_beam.txt'

import os
import numpy as np
from astropy.io import fits
from pixell import enmap, utils, reproject, wcsutils
import healpy as hp

def get_ymap_index_planck(comparison_group):
    if comparison_group == 0:
        #CIB deprojection
        ymap_name_list = ['no_deprojection',
                            'CIB+CMB_T=10.17beta=1.7',
                            'CIB+CMB_T=24beta=1.0',
                            'CIB+CMB_T=24beta=1.4',
                            'CIB+CMB_T=10.14beta=1.4',
                            'CIB+CMB_T=10.14beta=1.6']
        ymap_path_list = [
                        DAT+'CMB_ymap/Planck/ymap/no_deprojection_standard_full.fits',
                        DAT+'CMB_ymap/Planck/ymap/deproject_CMB_CIB_default_standard_full.fits',
                        DAT+'CMB_ymap/Planck/ymap/SED/deproject_CMB_CIB_beta1.0_T24_standard_full.fits',
                        DAT+'CMB_ymap/Planck/ymap/SED/deproject_CMB_CIB_beta1.4_T24_standard_full.fits',
                        DAT+'CMB_ymap/Planck/ymap/SED/deproject_CMB_CIB_beta1.4_T10.14_standard_full.fits',
                        DAT+'CMB_ymap/Planck/ymap/SED/deproject_CMB_CIB_beta1.6_T10.14_standard_full.fits']
    elif comparison_group == 1:
        #CIB + CIBdbeta deprojection
        ymap_name_list = [
                        'no_deprojection',
                        'CMB+CIB+CIBdbeta_T=10.17beta=1.7',
                        'CMB+CIB+CIBdbeta_T=24beta=1.0',
                        'CMB+CIB+CIBdbeta_T=24beta=1.4',
                        'CMB+CIB+CIBdbeta_T=10.14beta=1.4',
                        'CMB+CIB+CIBdbeta_T=10.14beta=1.6']
        ymap_path_list = [
                        DAT+'CMB_ymap/Planck/ymap/no_deprojection_standard_full.fits',
                        DAT+'CMB_ymap/Planck/ymap/deproject_CMB5_CIB_CIBdbeta_default_standard_full.fits',
                        DAT+'CMB_ymap/Planck/ymap/SED_dbeta/deproject_CMB5_CIB_CIBdbeta_beta1.0_T24_standard_full.fits',
                        DAT+'CMB_ymap/Planck/ymap/SED_dbeta/deproject_CMB5_CIB_CIBdbeta_beta1.4_T24_standard_full.fits',
                        DAT+'CMB_ymap/Planck/ymap/SED_dbeta/deproject_CMB5_CIB_CIBdbeta_beta1.4_T10.14_standard_full.fits',
                        DAT+'CMB_ymap/Planck/ymap/SED_dbeta/deproject_CMB5_CIB_CIBdbeta_beta1.6_T10.14_standard_full.fits']
    elif comparison_group == 2:
        ymap_name_list = [
                        'no_deprojection',
                        'CMB+CIB+CIBdbetadT_T=10.17beta=1.7',
                        'CMB+CIB+CIBdbetadT_T=24beta=1.0',
                        'CMB+CIB+CIBdbetadT_T=24beta=1.4',
                        'CMB+CIB+CIBdbetadT_T=10.14beta=1.4',
                        'CMB+CIB+CIBdbetadT_T=10.14beta=1.6']

        ymap_path_list = [
                        DAT+'CMB_ymap/Planck/ymap/no_deprojection_standard_full.fits',
                        DAT+'CMB_ymap/Planck/ymap/deproject_CMB5_CIB_CIBdbeta_CIBdT_default_standard_full.fits',
                        DAT+'CMB_ymap/Planck/ymap/SED_dbetadT/deproject_CMB5_CIB_CIBdbeta_CIBdT_beta1.0_T24_standard_full.fits',
                        DAT+'CMB_ymap/Planck/ymap/SED_dbetadT/deproject_CMB5_CIB_CIBdbeta_CIBdT_beta1.4_T24_standard_full.fits',
                        DAT+'CMB_ymap/Planck/ymap/SED_dbetadT/deproject_CMB5_CIB_CIBdbeta_CIBdT_beta1.4_T10.14_standard_full.fits',
                        DAT+'CMB_ymap/Planck/ymap/SED_dbetadT/deproject_CMB5_CIB_CIBdbeta_CIBdT_beta1.6_T10.14_standard_full.fits']
    
    return ymap_name_list, ymap_path_list


class deprojection_index_act:
    def __init__(self,directory:str,filename:str):
        self.dir = directory
        self.filename = filename
        self.deprotype = filename.split('_')[-3]
        self.nu = float(filename.split('_')[-2])
        self.T = float(filename.split('_')[-1][:-5])
        self.path = directory+filename
    
    def __repr__(self):
        return f'map object with {self.deprotype} {self.nu} {self.T}'

    def shape(self):
        with fits.open(self.path) as hdul:
            data = hdul[0].data
        return data.shape
    
    def read(self):
        with fits.open(self.path) as hdul:
            data = np.array(hdul[0].data,dtype=np.float32)
            header = hdul[0].header
        return data, header
    
    def generate_map(self,nside=2048):
        data,header = self.read()
        wcs = wcsutils.WCS(header)
        enmap_obj = enmap.ndmap(data,wcs)

        # print("Generating HEALPix map...")
        heapix_map = reproject.map2healpix(enmap_obj, nside=nside, lmax=3*nside-1, out=None, rot='cel,gal', spin=[0, 2], method='harm', order=1, extensive=False, bsize=100000, nside_mode='pow2', boundary='constant', verbose=True, niter=0)
        
        return heapix_map

class mask_index_act:
    def __init__(self,directory:str,filename:str):
        self.dir = directory
        self.filename = filename
        self.path = directory+filename
        self.name = None
    
    def __repr__(self):
        return f'mask_object at {self.filename}'
    
    def shape(self):
        with fits.open(self.path) as hdul:
            data = hdul[0].data
        return data.shape
    
    def read(self):
        with fits.open(self.path) as hdul:
            data = np.array(hdul[0].data,dtype=np.float32)
            header = hdul[0].header
        return data, header
    
    def generate_map(self,nside=2048):
        data,header = self.read()
        wcs = wcsutils.WCS(header)
        enmap_obj = enmap.ndmap(data,wcs)

        print("Generating HEALPix map...")
        heapix_mask = reproject.healpix_from_enmap(enmap_obj, nside=nside, lmax=3*nside-1)
        
        return heapix_mask

def get_ymap_index_act():
    filename_cib = os.listdir(DAT+'ACT/deprojections/cib/')
    filename_cibdbeta = os.listdir(DAT+'ACT/deprojections/cib_cibdbeta/')
    
    codex = []
    for i in range(len(filename_cib)):
        index =  deprojection_index_act(DAT+'ACT/deprojections/cib/',filename_cib[i])
        codex.append(index)
    for i in range(len(filename_cibdbeta)):
        index =  deprojection_index_act(DAT+'ACT/deprojections/cib_cibdbeta/',filename_cibdbeta[i])
        codex.append(index)
    return codex

def get_ymap_index_act_selected(deprotype,nu_range=[1.0,2.0],T_range=[10.7,24.0]):
    filename_cib = os.listdir(DAT+'ACT/deprojections/cib/')
    filename_cibdbeta = os.listdir(DAT+'ACT/deprojections/cib_cibdbeta/')
    
    codex = []
    for i in range(len(filename_cib)):
        index =  deprojection_index_act(DAT+'ACT/deprojections/cib/',filename_cib[i])
        
        if index.deprotype == deprotype:
            if index.nu<=nu_range[1] and index.nu>=nu_range[0]:
                if index.T<=T_range[1] and index.T>=T_range[0]:
                    codex.append(index)
                    
    for i in range(len(filename_cibdbeta)):
        index =  deprojection_index_act(DAT+'ACT/deprojections/cib_cibdbeta/',filename_cibdbeta[i])
        
        if index.deprotype == deprotype:
            if index.nu<=nu_range[1] and index.nu>=nu_range[0]:
                if index.T<=T_range[1] and index.T>=T_range[0]:
                    codex.append(index)
                    
    return codex

def get_ACT_mask_path():
    filename = os.listdir(DAT+'ACT/mask/')
    pathlist = []
    for i in range(len(filename)):
        pathlist.append(DAT+'ACT/mask/'+filename[i])
    return pathlist

def get_mask_index_act():
    filename = os.listdir(DAT+'ACT/mask/')
    codex = []
    for i in range(len(filename)):
        index = mask_index_act(DAT+'ACT/mask/',filename[i])
        codex.append(index)
    return codex

def read_beam():
    beam = np.loadtxt(BEAM,skiprows=1)
    ells = beam[:,0]
    beam = beam[:,1]
    return ells,beam

def read_composite_mask(apodize = True):
    if apodize:
        pathmask = DAT + 'unwiseact/act_mask_composite/healpix_act_mask_nside2048_apo1_5.fits'
    else:
        pathmask = DAT + 'unwiseact/act_mask_composite/healpix_act_mask_nside2048.fits'
    mask = hp.read_map(pathmask)
    return mask

if __name__ == '__main__':
    # codex = get_ymap_index_act()
    
    # for i in range(len(codex)):
    #     print(codex[i].deprotype,codex[i].nu,codex[i].T)
        
    # pathlist = get_ACT_mask_path()
    # for i in range(len(pathlist)):
    #     print(pathlist[i])
        
    ells,beam = read_beam()
    print(ells,beam)
    
    # codex = get_ymap_index_act_selected(deprotype = 'cibdBeta',nu_range=[1.0,1.4],T_range=[10.7,12.0])
    # for i in range(len(codex)):
    #     print(codex[i].deprotype,codex[i].nu,codex[i].T)
