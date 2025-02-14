import numpy as np
import assets.deprojection_index as di
from astropy.io import fits
from pixell import enmap, utils, reproject, wcsutils
import healpy as hp

DAT = di.DAT
NSIDE = 2048

def main():
    # map_ACT_codex = di.get_ymap_index_act()
    mask_ACT_codex = di.get_mask_index_act()
    mask_act_enmap = np.ones((10320, 43200)).astype(np.float32)
    for index in mask_ACT_codex:
        mask_enmap,header = index.read()
        wcs = wcsutils.WCS(header)
        mask_act_enmap *= mask_enmap
        print("finished reading {}".format(index))

    enmap_obj = enmap.ndmap(mask_act_enmap, wcs)
    print("Generating HEALPix map...")
    healpix_mask = reproject.map2healpix(enmap_obj, nside=NSIDE, lmax=3*NSIDE-1, rot='equ,gal', verbose=True, niter=0)
    
    print("Writing HEALPix map...")
    hp.write_map(DAT+'unwiseact/act_mask_composite/healpix_act_mask_nside2048.fits', healpix_mask, overwrite=True)
    
    print("Finished writing HEALPix map.")

def readmask():
    mask = hp.read_map(DAT+'unwiseact/act_mask_composite/healpix_act_mask_nside2048.fits')
    return mask

def cutoffmask():
    mask = hp.read_map(DAT+'unwiseact/act_mask_composite/healpix_act_mask_nside2048.fits')
    mask[mask < 0.0] = 0.0
    mask[mask >= 1.0] = 1.0
    hp.write_map(DAT+'unwiseact/act_mask_composite/healpix_act_mask_nside2048_cutoff.fits', mask, overwrite=True)
    return mask
if __name__ == '__main__':
    cutoffmask()