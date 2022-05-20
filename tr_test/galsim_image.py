#! /usr/bin/env python3

import numpy  as np
import galsim
from astropy.stats import gaussian_sigma_to_fwhm
import os
 
if __name__ == '__main__':
    # Model params
    m = 2.3
    std = 0.08
    shape = (500, 500)
    psf_sigma = 1.0   
    PSF = galsim.Gaussian(flux=1., sigma=psf_sigma)
    out_filename = os.path.join('images','sources.fits') 
    file = open('init_sources.cat', 'w')
    print('N_source', 'x_0', 'y_0', 'Flux', 're', 'n', 'g1', 'g2', 'theta', sep=' | ' , end=' | \n', file=file)

    #Создание подстилающего изображения
    
    img = galsim.ImageF(500, 500)
    random_seed = 1234567
    rng = galsim.BaseDeviate(random_seed)
    

    # Создаем галактики

    n_sources = 20
    
    flux  = [100, 1000]
    r_eff = [10,50]
    n     = [10, 5]
    x_0   = [0, shape[1]]
    y_0   = [0, shape[0]]
    ellip = [0.3, 0.9]
    theta = [0, np.pi]

    for i in range(n_sources):
        fx = np.random.randint(flux[0], flux[1])
        re = np.random.randint(r_eff[0], r_eff[1])/10
        n = np.random.randint(10, 50)/10
        x0 = np.random.randint(x_0[0],  x_0[1])
        y0 = np.random.randint(y_0[0],  y_0[1])
        g1 = np.random.randint(1, 30)/100
        g2 = np.random.randint(1, 30)/100        
        print("%9d %5d %5d %6d %4.1f %3.1f %3.2f %3.2f" % (i, x0-1, y0-1, fx, re, n, g1, g2),  file=file)

        sers = galsim.Sersic(n=n, half_light_radius=re, flux=fx)
        sers = sers.shear(g1=g1, g2=g2) 
        sers = galsim.Convolve([sers, PSF])

        sers.drawImage(image=img, center=(x0,y0), add_to_image=True) 
            
    
    #Создаем звезды
    n_sources = 10
    
    x_mean   = [10, shape[0]]
    y_mean   = [10, shape[0]]
    x_stddev = [1, 1]
    y_stddev = [1, 1]
    theta    =  [0, 2 * np.pi]
    x = np.arange(shape[0])[np.newaxis,:]
    y = np.arange(shape[1])[:,np.newaxis]
    for i in range(n_sources):
        xm = np.random.randint(x_mean[0], x_mean[1])
        ym = np.random.randint(y_mean[0], y_mean[1])
        xd = np.random.randint(1,  x_stddev[0]+1)
        yd = np.random.randint(1,  y_stddev[0]+1)
        th = np.random.randint(theta[0], theta[1])
        print("%9d %5d %5d %6d " % (i, xm-1, ym-1, (i+1)*100),  file=file)

        gauss1 = galsim.Gaussian(fwhm=xd*gaussian_sigma_to_fwhm, flux=(i+1)*100)
        gauss1.drawImage(image=img, center=(xm,ym), add_to_image=True)

        # draw profile through LSST filters
    
    gaussian_noise = galsim.GaussianNoise(rng, sigma=std)
    
    
    
    images = []
    newImg = img.copy()
    newImg.addNoise(gaussian_noise)
    images.append(newImg)
    galsim.fits.writeCube(images, out_filename)
    

