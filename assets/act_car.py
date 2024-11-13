DAT = '/mnt/d/data_large/unwise_sz/'
BEAM = '/mnt/d/data_large/unwise_sz/ACT/act_beam/effective_beam.txt'

import os
import numpy as np
from astropy.io import fits
from pixell import enmap, utils, reproject, wcsutils
import healpy as hp

class act_car_file:
    def __init__(self,directory:str,filename:str):
        self.dir = directory
        self.filename = filename
        self.path = directory+filename

    def __repr__(self):
        return f'file object at {self.path}'
    
    def read_map_to_enmap(self):
        '''Read the map to an enmap object, return the enmap object'''
        path = self.path
        enmap_data = enmap.read_map(path)
        return enmap_data
    
    def read_map_to_array(self):
        '''Read the map to a numpy array, return a dictionary with data and wcs'''
        hdul = fits.open(self.path)
        data = hdul[0].data
        wcs = wcsutils.WCS(hdul[0].header)
        hdul.close()
        return {'data':data,'wcs':wcs}
    
class act_car_map(act_car_file):
    def __init__(self,directory:str,filename:str):
        super().__init__(directory,filename)
        self.deprotype = None
        self.beta = float(filename.split('_')[-2])
        self.T = float(filename.split('_')[-1][:-5])
    
    def __repr__(self):
        return f'map object of {self.deprotype} with beta = {self.beta}, T = {self.T}'
    
class act_car_mask(act_car_file):
    def __init__(self,directory:str,filename:str):
        super().__init__(directory,filename)
    
    def __repr__(self):
        return f'mask object at {self.path}'


def get_ymap_index_act(verbose=False):
    filename_cib = os.listdir(DAT+'ACT/deprojections/cib/')
    filename_cibdbeta = os.listdir(DAT+'ACT/deprojections/cib_cibdBeta/')
    filename_cibdbetadT = os.listdir(DAT+'ACT/deprojections/cib_cibdBeta_cibdT/')
    
    car_maps = []
    for i in range(len(filename_cib)):
        car_map =  act_car_map(DAT+'ACT/deprojections/cib/',filename_cib[i])
        car_map.deprotype = 'cib'
        car_maps.append(car_map)
        print(f'loaded {car_map}') if verbose else None
    for i in range(len(filename_cibdbeta)):
        car_map =  act_car_map(DAT+'ACT/deprojections/cib_cibdBeta/',filename_cibdbeta[i])
        car_map.deprotype = 'cib_cibdBeta'
        car_maps.append(car_map)
        print(f'loaded {car_map}') if verbose else None
    for i in range(len(filename_cibdbetadT)):
        car_map =  act_car_map(DAT+'ACT/deprojections/cib_cibdBeta_cibdT/',filename_cibdbetadT[i])
        car_map.deprotype = 'cib_cibdBeta_cibdT'
        car_maps.append(car_map)
        print(f'loaded {car_map}') if verbose else None
    return car_maps

def get_ymap_index_act_selected(deprotype,beta_range=[1.0,2.0],T_range=[10.7,24.0],verbose=False):
    car_maps = get_ymap_index_act()
    selected = []
    for car_map in car_maps:
        if car_map.deprotype == deprotype:
            if car_map.beta >= beta_range[0] and car_map.beta <= beta_range[1]:
                if car_map.T >= T_range[0] and car_map.T <= T_range[1]:
                    selected.append(car_map)
                    print(f'selected {car_map}') if verbose else None
    return selected

def get_mask_index_act(verbose=False):
    '''Get all the ACT masks, return a list of mask objects'''
    filename = os.listdir(DAT+'ACT/masks/')
    car_masks = []
    for i in range(len(filename)):
        car_mask =  act_car_mask(DAT+'ACT/masks/',filename[i])
        car_masks.append(car_mask)
        print(f'loaded {car_mask}') if verbose else None
    return car_masks


def read_beam():
    '''Read the ACT beam file, return ells and beam values in dictionary'''
    beam = np.loadtxt(BEAM,skiprows=1)
    ells = beam[:,0]
    beam = beam[:,1]
    return {'ells':ells,'beam':beam}