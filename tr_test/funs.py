#! /usr/bin/env python3


from astropy.io import fits
from astropy.stats import sigma_clipped_stats, gaussian_fwhm_to_sigma
from astropy.nddata import NDData
from astropy.convolution import Gaussian2DKernel, convolve
from astropy.visualization import simple_norm

import matplotlib.pyplot as pl
import matplotlib.colors as colors
import numpy as np

from photutils.detection import find_peaks
from photutils.segmentation import detect_threshold, detect_sources, SourceCatalog
from photutils.psf import extract_stars

from tractor.psfex import PixelizedPsfEx
from tractor import *
from tractor.sersic import *

import os

 

def Bkg2D(data):
    
    from astropy.stats import SigmaClip
    from photutils.background import Background2D
   
    # Сырой поиск источников
    threshold = detect_threshold(data, nsigma=3.)
    sigma = 3.0 * gaussian_fwhm_to_sigma  # FWHM = 3.
    kernel = Gaussian2DKernel(sigma, x_size=3, y_size=3)
    convolved_data = convolve(data, kernel, normalize_kernel=True)

    # Создание маски для повышения качества оценки фона
    segm = detect_sources(convolved_data, threshold, npixels=5)
    mask = (segm.data > 0)
    sigma_clip = SigmaClip(sigma=3.)
 
    #Оценка фона
    Bkg = Background2D(data, (20,20), mask=mask,  filter_size=(5,5), sigma_clip=sigma_clip, fill_value=0.0) 
    
    return Bkg   


def Sources(data, bkg):
    
    # Фильтрация с учетом полученного фона
    threshold = bkg.background_median + 2.0 * bkg.background_rms

    #   Фильтрация изображения Гауссовым ядром (? может другие)
    sigma = 3.0 * gaussian_fwhm_to_sigma  # FWHM = 3.
    kernel = Gaussian2DKernel(sigma, x_size=3, y_size=3)
    convolved_data = convolve(data, kernel, normalize_kernel=True)
    
    # Поиск источников
    segm = detect_sources(convolved_data, threshold, npixels=5)
    
    # Маскирова больших (слишком) источников ???
    
    # Деблендирование
    from photutils.segmentation import deblend_sources
    segm_deblend = deblend_sources(convolved_data, segm, npixels=5, nlevels=32, contrast=0.005)

    return segm_deblend


def psfe(sub_data, tresh = 60., size = 25, dbg=False):

    
    # Поиск точечных источников
    peaks = find_peaks(sub_data, threshold=tresh)
    
    if len(peaks) == 0:
        print('Источников не обнаружено. Текущий threshold =', tresh)
        os._exit()    
    
    x = peaks['x_peak']  
    y = peaks['y_peak']  
    # Надо для get_stars
    nddata = NDData(data=sub_data)
      
    # Маскировка источников близких к краю

    hsize = (size - 1) // 2

    mask = ((x > hsize) & (x < (sub_data.shape[1] -1 - hsize)) &
        (y > hsize) & (y < (sub_data.shape[0] -1 - hsize))) 

    from astropy.table import Table
    stars_tbl = Table()

    stars_tbl['x'] = x[mask]
    stars_tbl['y'] = y[mask]

    if len(stars_tbl) == 0:
        print('Все источники перекрыты.', 'Текущий threshold =', tresh, 'Текущий барьер маски =', hsize)
        os._exit()

    print('Извлечение звезд')
    stars = extract_stars(nddata, stars_tbl, size=size)
    
    print('Подгонка PSF гауссианой')

    # Нормализация звездных изображений
    lis = [] 
    for i in range(len(stars)):
        lis.append(stars[i].data/np.sum(stars[i].data))

    # Создание источников и изображений для Трактора
    tims =[]
    pos = []
    for image in lis:
        tim = Image(data=image, invvar=np.ones_like(image),
        psf=NCircularGaussianPSF([1],[1]),
        wcs=NullWCS(), photocal=NullPhotoCal(),
        sky=ConstantSky(0.))
        tims.append(tim)
        pos.append(PointSource(PixPos(size//2,size//2) , Flux(1.)))

    # Создание модели Трактора
    tr = Tractor(tims, pos)

    # Рендер начальных изображений
    mod0 = [tr.getModelImage(i) for i in range(len(stars))]
    chi0 = [tr.getChiImage(i) for i in range(len(stars))]

    # Фиксация яркости и фона. (Свободны координаты и параметры PSF)
    for i in range(len(stars)):
        tr.images[i].freezeParam('sky')
        tr.catalog[i].freezeParam('brightness')           
    #print(tr.printThawedParams())
    
    # Подгонка
    for i in range(10):
        dlnp,X,alpha = tr.optimize()
        print( 'dlnp', dlnp)
        if dlnp < 1e-3:
            break    

    # Рендер финальных изображений
    mod0 = [tr.getModelImage(i) for i in range(len(stars))]
    chi0 = [tr.getChiImage(i) for i in range(len(stars))]
   

    # Проверка качества подгонки PSF
    if dbg:
        nrows = 2
        ncols = 2
        fig, ax = pl.subplots(nrows=nrows, ncols=ncols, figsize=(20, 20),
                       squeeze=True)
        ax = ax.ravel()
        for i in range(nrows * ncols):
            norm = simple_norm(chi0[i], 'log', percent=99.)
            ax[i].imshow(chi0[i], norm=norm, origin='lower', cmap='viridis')
        pl.show()




    # Усреднение получнных данных (медианное). Возможно стоит поступать иначе ???
    psf_sigm = []
    wh = []
    for i in range(len(stars)):
        psf_sigm.append(tr.images[i].psf.sigmas[0])
        wh.append(tr.images[i].psf.weights[0])    
    psf_sigm = np.median(psf_sigm)
    wh       = np.median(wh)
    
    # Сохранение полученного PSF в FITS
    fits.writeto('PSF.fits', mod0[0], overwrite=True)

    return psf_sigm, wh
    
    




def main(image, hdr, pathout=None):


    home = os.getcwd()
    # Сюда сохраняем изображения. Иначе в рабочей дирректории
    if pathout is not None:
        os.chdir(pathout)   
    
    print('Draw init image')    
    fits.writeto('Init_imag.fits', image, overwrite=True)
    
    print('Estimate Bkg')
    Bkg = Bkg2D(image)
    sky = Bkg.background_median
    noice_sigm = Bkg.background_rms
    noisesigma = Bkg.background_rms_median

    print('sky=', sky, '\nsigma=', noisesigma)
    sub = image - sky


    
    fits.writeto('Bakground.fits', Bkg.background, overwrite=True)
    
    fits.writeto('SUB.fits', sub, overwrite=True)


    print('Define Sources')
    segm = Sources(image, Bkg)
    
    fits.writeto('Segm.fits', segm.data, overwrite=True)
    

    
    print('Create catalog for Tractor')

    cat = SourceCatalog(sub, segm)
    tbl = cat.to_table()
    N_sources = len(tbl)
    print('N_sources', N_sources)

    print('Crea PSF')
    psf, wh = psfe(sub)
    
    print('Fiting by PSs')
    pixscale = hdr['GS_SCALE']
    # Create Tractor Image
    tim = Image(data=image, invvar=np.ones_like(image) / (noice_sigm**2),
    psf=NCircularGaussianPSF([psf],[1.]),
    wcs=NullWCS(pixscale=pixscale), photocal=NullPhotoCal(),
    sky=ConstantSky(sky))
    
    # Create Tractor source with approximate position and flux
    src_p = [PointSource(PixPos(cat.xcentroid[i], cat.ycentroid[i]), Flux(cat.kron_flux[i])) for i in range(N_sources)]

    # Create Tractor object itself
    tractor = Tractor([tim], src_p)

    # Render the model image
    mod0 = tractor.getModelImage(0)
    chi0 = tractor.getChiImage(0)

    # Plots
    fits.writeto('Mod0.fits', mod0, overwrite=True)
    fits.writeto('Chi0.fits', chi0, overwrite=True)
    


    # Freeze all image calibration params -- just fit source params
    tractor.freezeParam('images')


    # Take several linearized least squares steps    
    for i in range(100):
        dlnp,X,alpha = tractor.optimize(damp=0.1)
        print( 'dlnp', dlnp)
        if dlnp < 1e-3:
            break

    # Get the fit model and residual images for plotting
    mod = tractor.getModelImage(0)
    chi = tractor.getChiImage(0)
    
    fits.writeto('Mod.fits', mod+0.08*np.random.normal(size=mod.shape), overwrite=True)
    fits.writeto('Chi.fits', chi, overwrite=True)


    print('Поиск кандидатов в Галактики')
    tresh = 10.
    err = find_peaks(chi, threshold=tresh)
    x_0 = y_0 = -10 
    
    if len(err) == 0:
        print('Ничего не найдено.','\n'+'Текущий threshold =', tresh)

    print('psf=', psf) 

    # Отбор найденных пиков. Близкие пики -> один источник. Надо сделать красивее ???
    galaxies = []
    src = []    
    for x, y, f in err:
        if np.hypot(x-x_0, y-y_0) > 3*psf:
            galaxies.append((x,y))            
        x_0 = x; y_0 = y

    # Создание нового каталога
    
    for i in range(len(tractor.catalog)):
        x, y = tractor.catalog[i].pos
        fl   = tractor.catalog[i].brightness

        
        for xi, yi in galaxies:
            if np.hypot(x-xi,y-yi) < 3*psf:
                src.append(SersicGalaxy(PixPos(x, y), fl, EllipseE(2., 0.5, 0.5), SersicIndex(3.)))
                key = 0
                break
            else:
                key = 1
        if key == 1:
            src.append(PointSource(PixPos(x, y), fl))
    #src.append(SersicGalaxy(PixPos(155, 155), Flux(100.), EllipseE(2., 0.5, 0.5), SersicIndex(3.)))
    # Создание модели Трактора
    trac = Tractor([tim], src)
    
    # Рендер начального изображения
    mod01 = trac.getModelImage(0)
    chi01 = trac.getChiImage(0)
    
    # Заморозка изображения
    trac.freezeParam('images')

    # Подгонка
    for i in range(100):
        dlnp,X,alpha = trac.optimize(damp=0.1)
        print( 'dlnp', dlnp)
        if dlnp < 1e-8:
            break
    
    # Итоговое изображение
    mod1 = trac.getModelImage(0)
    chi1 = trac.getChiImage(0)

    fits.writeto('Mod1.fits', mod1, overwrite=True)
    fits.writeto('Chi1.fits', chi1, overwrite=True)

    os.chdir(home) 
    
    #for i  in zip(trac.getParamNames(),trac.getParams()):
    #    print(i)

    return trac.getCatalog()


if __name__ == '__main__':

    from photutils.datasets import load_simulated_hst_star_image

    trueimage = load_simulated_hst_star_image().data
    noisesigma = 0.08

    '''
    #Point_sources
    W,H = 500,500
    noisesigma = 0.08
    cx,cy = 150, 150
    flux = 1200.

    cx2, cy2 = 200, 200
    flux2 = 3000

    # PSF size
    psfsigma = 2.

    # Create synthetic Gaussian star image
    G = np.exp(((np.arange(W)-cx)[np.newaxis,:]**2 +
(np.arange(H)-cy)[:,np.newaxis]**2)/(-2.*psfsigma**2))
    trueimage = flux * G/G.sum()

    G2 = np.exp(((np.arange(W)-cx2)[np.newaxis,:]**2 +
(np.arange(H)-cy2)[:,np.newaxis]**2)/(-2.*psfsigma**2))
    trueimage += flux2 * G2/G2.sum() 
    '''

    image = trueimage +2.3 + noisesigma * np.random.normal(size=trueimage.shape)
    
    # Start program
    main(image)
