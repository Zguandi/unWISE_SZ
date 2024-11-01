import numpy as np
import assets.make_galaxy_map as mgm
from astropy.io import fits
from pixell import enmap, utils, reproject, wcsutils
import healpy as hp

DAT = mgm.DAT
NSIDE = 2048

def main():
    # map_ACT_codex = di.get_ymap_index_act()
    mask_unWISE = mgm.readmask()
    
    print("Writing HEALPix map...")
    hp.write_map(DAT+'unwiseact/unwise_mask_composite/healpix_unwise_mask_nside2048.fits', mask_unWISE, overwrite=True,dtype=np.float32)
    
    print("Finished writing HEALPix map.")

if __name__ == '__main__':
    main()