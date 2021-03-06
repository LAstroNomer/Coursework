#! /usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
from tractor import *
from tractor.psfex import PixelizedPsfEx
from astropy.io import fits

if __name__ == '__main__':

    
    PSF = PixelizedPsfEx(fn='default.fits', ext=1)
    shape = PSF.shape
    W = shape[0]
    H = shape[1]

    img = np.zeros((W, H))

    tim = Image(data=img, invvar=np.ones_like(img),
    psf=PSF,
    wcs=NullWCS(), photocal=NullPhotoCal(),
    sky=ConstantSky(0.))
    src = PointSource(PixPos(W//2, H//2), Flux(1))
    tr = Tractor([tim], [src])
    mod = tr.getModelImage(0)

    fits.writeto('PSF_IMAG.fits', mod , overwrite=True)
