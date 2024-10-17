from pixell import enmap, utils,enplot, reproject
import healpy as hp
import numpy as np
from matplotlib import pyplot as plt

def eshow(x,**kwargs):
    plots = enplot.get_plots(x, **kwargs)
    enplot.show(plots, method = "ipython")
    
