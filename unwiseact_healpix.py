import os
import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
import pymaster as nmt
from astropy.io import fits
from scipy.interpolate import interp1d
from assets import make_galaxy_map as mgm
from assets import galaxy_map_read as gmr
from assets import deprojection_index as di
import pandas as pd

DAT = di.DAT
NSIDE = 2048
BIN = 50
OUTPATH = '/mnt/c/Users/gdzhao/projects/unwise_sz/unwiseact_results/healpix_midz/'
wisesample = 'midz'
deprotype = 'cibdBeta'


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

def unwise_zsample_act(wisesample:str,deprotype:str,nu_range=[1.0,1.2],T_range=[10.7,12.0]):
    
    print("reading unwisemap")
    unwise_map = gmr.makemap_healpix(wisesample)
    
    print("reading unwisemask")
    unwise_mask = gmr.readmask(wisesample)
    
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

def unwise_zsample_act_jitsave(wisesample:str,deprotype:str,nu_range=[1.0,1.2],T_range=[10.7,12.0]):
    
    print("reading unwisemap")
    unwise_map = gmr.makemap_healpix(wisesample)
    
    print("reading unwisemask")
    unwise_mask = gmr.readmask(wisesample)
    
    print("reading act mask")
    act_mask = di.read_composite_mask(apodize=False)
    
    act_map_codex = di.get_ymap_index_act_selected(deprotype,nu_range,T_range)
    print("number of act samples: ",len(act_map_codex))
    
    print("initializing pseudocl...")
    b = nmt.NmtBin.from_nside_linear(NSIDE, 50)
    galaxy = nmt.NmtField(unwise_mask, [unwise_map])
    
    ells,beam = di.read_beam()
    
    print("main loop for act map...")
    total_depro = len(act_map_codex)
    
    for i in range(total_depro):
        index = act_map_codex[i]
        
        cols = [b.get_effective_ells()]
        names = ['ell','cl_gy','cl_yy']
        
        print("  generating y for "+index.filename,end='\r')
        y_field = nmt.NmtField(act_mask, [index.generate_map()],beam=beam)
        
        print("  computing gy for "+index.filename,end='\r')
        cl_gy = nmt.compute_full_master(galaxy, y_field, b)
        print("  namaster gy complete for "+index.filename,end='\r')
        
        print("  computing yy for "+index.filename,end='\r')
        cl_yy = nmt.compute_full_master(y_field, y_field, b)
        print("  namaster yy complete for "+index.filename,end='\r')
        
        cols.append(cl_gy[0])
        cols.append(cl_yy[0])
        
        outfile_path = OUTPATH + '{}_{}_{}_T_{}_results.txt'.format(wisesample,deprotype,index.nu,index.T)
        np.savetxt(outfile_path,np.array(cols).T,header=' '.join(names))
        print("  namaster gyyy results saved for " + index.filename,"{} of {}".format(i+1,total_depro))

def unwise_auto(wisesample:str):
    
    print("reading unwisemap")
    unwise_map = gmr.makemap_healpix(wisesample)
    
    print("reading unwisemask")
    unwise_mask = gmr.readmask(wisesample)
    
    print("initializing pseudocl...")
    b = nmt.NmtBin.from_nside_linear(NSIDE, 50)
    galaxy = nmt.NmtField(unwise_mask, [unwise_map])
    
    cols = [b.get_effective_ells()]
    names = ['ell']
    
    print("computing auto for unwise map...")
    cl_auto = nmt.compute_full_master(galaxy, galaxy, b)
    
    names.append('cl_auto')
    cols.append(cl_auto[0])
    
    print("saving results...")
    np.savetxt(OUTPATH+'{}_auto_results.txt'.format(wisesample),np.array(cols).T,header=' '.join(names))

def read_results(filename:str):
    path = OUTPATH+filename
    # ignore the '#' in the header
    
    data = np.loadtxt(path,skiprows=1)
    header = open(path).readline().strip().split()
    header.remove('#')
    table = {header[i]:data[:,i] for i in range(len(header))}
    return table

def read_results_jit(filename:str):
    path = OUTPATH+filename
    # ignore the '#' in the header
    
    unwisesample = filename.split('_')[0]
    deprotype = filename.split('_')[1]
    nu = filename.split('_')[2]
    T = filename.split('_')[4]
    
    print("reading sample: ",unwisesample,deprotype,nu,T)
    
    data = np.loadtxt(path,skiprows=1)
    header = open(path).readline().strip().split()
    header.remove('#')
    table = {header[i]:data[:,i] for i in range(len(header))}
    return table

if __name__ == "__main__":
    # unwiseact('cibdBeta',nu_range=[1.0,1.1],T_range=[10.7,12.0])
    # unwiseact('cibdBeta',nu_range=[1.0,1.2],T_range=[10.7,12.0])
    # table = read_results('cibdBeta')
    # print(table.keys())
    # print(table['cibdBeta_1.0_T_10.7_yy'])
    
    if not os.path.exists(OUTPATH):
        os.makedirs(OUTPATH)
        
    unwise_zsample_act_jitsave(wisesample,deprotype,nu_range=[1.0,2.0],T_range=[10.7,24.0])
    
    # unwise_auto('lowz')
