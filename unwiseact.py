import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
import pymaster as nmt
from astropy.io import fits
from scipy.interpolate import interp1d
from assets import make_galaxy_map as mgm
from assets import deprojection_index as di
import pandas as pd

DAT = di.DAT
NSIDE = 2048
BIN = 50
OUTPATH = '/mnt/c/Users/gdzhao/projects/unwise_sz/unwiseact_results/'

def unwiseact(deprotype:str,nu_range=[1.0,1.2],T_range=[10.7,12.0]):
    
    print("reading unwisemap")
    unwise_map = mgm.makemap()
    
    print("reading unwisemask")
    unwise_mask = mgm.readmask()
    
    print("reading act mask")
    act_mask = di.read_composite_mask(apodize=False)
    
    act_map_codex = di.get_ymap_index_act_selected(deprotype,nu_range,T_range)
    print("number of act samples: ",len(act_map_codex))
    
    print("initializing pseudocl...")
    b = nmt.NmtBin.from_nside_linear(NSIDE, 50)
    galaxy = nmt.NmtField(unwise_mask, [unwise_map])
    
    ells,beam = di.read_beam()
    
    cols = [b.get_effective_ells()]
    names = ['ell']
    
    print("main loop for act map...")
    for index in act_map_codex:
        print("  generating y for "+index.filename,end='\r')
        y_field = nmt.NmtField(act_mask, [index.generate_map()],beam=beam)
        
        print("  computing gy for "+index.filename,end='\r')
        cl_gy = nmt.compute_full_master(galaxy, y_field, b)
        print("  namaster gy complete for "+index.filename)
        
        print("  computing yy for "+index.filename,end='\r')
        cl_yy = nmt.compute_full_master(y_field, y_field, b)
        print("  namaster yy complete for "+index.filename)
        
        names.append('{}_{}_T_{}_gy'.format(index.deprotype,index.nu,index.T))
        cols.append(cl_gy[0])
        
        names.append('{}_{}_T_{}_yy'.format(index.deprotype,index.nu,index.T))
        cols.append(cl_yy[0])
    
    print("saving results...")
    np.savetxt(OUTPATH+'{}_results_first_sample.txt'.format(deprotype),np.array(cols).T,header=' '.join(names))

# def beam_manual_interp(ells,beam):
#     f = interp1d(ells,beam)
#     return f

def read_results(filename:str):
    path = OUTPATH+filename
    # ignore the '#' in the header
    
    data = np.loadtxt(path,skiprows=1)
    header = open(path).readline().strip().split()
    header.remove('#')
    table = {header[i]:data[:,i] for i in range(len(header))}
    return table
    
if __name__ == "__main__":
    # unwiseact('cibdBeta',nu_range=[1.0,1.1],T_range=[10.7,12.0])
    unwiseact('cibdBeta',nu_range=[1.0,1.2],T_range=[10.7,12.0])
    # table = read_results('cibdBeta')
    # print(table.keys())
    # print(table['cibdBeta_1.0_T_10.7_yy'])
    

