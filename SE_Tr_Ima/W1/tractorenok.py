#! /usr/bin/env python3 

from tractor import *
from astropy.table import Table
from tractor import *
from tractor.sersic import *
import numpy as np

if __name__ == '__main__':
    cat = open('main.cat', 'r')   
    cat_tbl = None
    hdr = []
    for sr in cat:
        sp = sr.split()
        if sp[0] == '#':
            hdr.append(sp[2])
        else:
            if cat_tbl is None:
                cat_tbl = Table(names=hdr)
            cat_tbl.add_row(sp)
    
    cx = cat_tbl['X_IMAGE']
    cy = cat_tbl['Y_IMAGE']
    cgf = cat_tbl['FLUX_SPHEROID']
    ctp = cat_tbl['SPREAD_MODEL']
    ctp_err = cat_tbl['SPREADERR_MODEL']
    cn = cat_tbl['SPHEROID_SERSICN']
    cref = cat_tbl['SPHEROID_REFF_IMAGE']
    ce1 = cat_tbl['ELLIP1MODEL_IMAGE']
    ce2 = cat_tbl['ELLIP2MODEL_IMAGE']
    src = []
    eps = 5.0 * 10**(-3)
    k = 4
    switch = np.sqrt(eps**2 + (k * ctp_err)**2)
    gc = []
    pc = []
    
    for x, y, gf, tp, n, re, e1, e2, key in zip(cx, cy, cgf, ctp, cn, cref, ce1, ce2, switch):
#        if tp <0: 
#            pass
        if tp < 1:
            pc.append((x, y))
            src.append(PointSource(PixPos(x, y), Flux(gf)))
        else:
            print('gal x = %f5.3, y = %f5.3, mag = %f5.3, re = %f5.3, n = %f5.3' %(x, y, 20.5 - 2.5*np.log10(gf), re, n))
            gc.append((x, y))
            src.append(SersicGalaxy(PixPos(x, y), Flux(gf), EllipseE(5., 0.,0.), SersicIndex(3.)))
    
    import os
    home = os.getcwd()
    os.chdir('results')
    from astropy.io import fits 
    bkg = fits.getdata('sex_bkg.fits')
    print('bkg = %4.2f' %np.median(bkg))
    
    sigmN = fits.getdata('sex_rms.fits')
    sigmN = np.median(sigmN)
    print('Noice Sigm = %4.2f' %sigmN)
    os.chdir(home)
    
    from tractor.psfex import PixelizedPsfEx
    PSF = PixelizedPsfEx(fn='default.fits')
    img, hdr = fits.getdata('sources.fits', header=True)
    img = img
    #pixscale = hdr['PXSCAL1']
  
    tim = Image(data=img, invvar=np.ones_like(img) / (sigmN**2),
    psf=PSF,
    wcs=NullWCS(), photocal=NullPhotoCal(),
    sky=ConstantSky(np.median(bkg)))

    tractor = Tractor([tim], src)

    mod0 = tractor.getModelImage(0)  
    chi0 = tractor.getChiImage(0)
    fits.writeto('Trmod0.fits', mod0 , overwrite=True)
    fits.writeto('Trchi0.fits', chi0, overwrite=True)
    tractor.freezeParam('images')

    for _ in range(100):
        dlnp,X,alpha = tractor.optimize(damp = 0.1)
        print( 'dlnp', dlnp)
        if dlnp < 1e-2:
            break
    mod = tractor.getModelImage(0)
    chi = tractor.getChiImage(0)

    fits.writeto('Trmod.fits', mod , overwrite=True)
    fits.writeto('Trchi.fits', chi, overwrite=True)
    fits.writeto('TSub.fits', img - mod, overwrite=True)
    new_cat = tractor.getCatalog()
    

    from matplotlib import pyplot as plt
    from astropy.visualization import simple_norm

    norm = simple_norm(img, 'log', percent=99.)
    plt.imshow(img, norm=norm, origin='lower', cmap='viridis' )
    # plot init
    x_0 = cx
    y_0 = cy    
    plt.plot(x_0, y_0, 'r+', ms=3, label='init')

    # plot galaxies
    for x0, y0 in gc:
        plt.plot(x0, y0, 'o', mec='y', mfc='none', label='galaxies')


    # plot PSs
    for x0, y0 in pc:
        plt.plot(x0, y0, 's', mec='r', mfc='none', label='PSs')
    
    #plt.legend()
    plt.savefig('sources_res.png')

    p_file = open('tr_point.cat','w')
    g_file = open('tr_galaxies.cat', 'w')

    print('N_source', 'x_0', 'y_0', 'Flux', 're', 'n', 'e1', 'e2', sep=' | ' , end=' | \n', file=g_file)
    print('N_source', 'x_0', 'y_0', 'Flux', sep=' | ' , end=' | \n', file=p_file)

    p = g = 0
    for i in range(len(new_cat)):
        if len(new_cat[i]) <=3:
            x, y, flux = new_cat[i].getParams()
            print("%9d %5.3f %5.3f %6.3f " % (p, x, y, flux), file=p_file)
            p += 1
        else:
            x0, y0, fx, re, e1, e2, n = new_cat[i].getParams()
            print("%9d %5.3f %5.3f %6.3f %4.1f %3.1f %3.2f %3.2f" % (g, x0, y0, fx, re, n, e1, e2), file=g_file)
            g += 1
