#! /usr/bin/env python3

import numpy  as np
from matplotlib import pyplot as plt
from astropy.io import fits
from astropy.table import Table
from astropy.visualization import simple_norm
import os

if __name__ == '__main__':

    
    init = open('init_sources.cat', 'r')
    
    key = 0 
    for lines in init:
        if key == 0:
            header = lines.split('|')[:-1]
            tbl_init = Table(names=header)
            key = 1    
        else:
            lines = list(map(float, lines.split()))
            while len(lines)<len(header):
                 lines.append(0)      
            tbl_init.add_row(lines)
    ind = len(tbl_init[' n '][tbl_init[' n '] > 0])
    tbl_init_g = tbl_init[0:ind]
    tbl_init_p = tbl_init[ind:]
    tbl_init_g.sort(' x_0 ')
    tbl_init_p.sort(' Flux ')
    
    gal = open('galaxies.cat', 'r')    
    key = 0
    for lines in gal:
        if key == 0:
            header = lines.split('|')[:-1]
            tbl_g = Table(names=header)
            key = 1    
        else:
            lines = list(map(float, lines.split()))
            tbl_g.add_row(lines)
    tbl_g.sort(' x_0 ')


    point = open('point.cat', 'r')

    key = 0
    for lines in point:
        if key == 0:
            header = lines.split('|')[:-1]
            tbl_p = Table(names=header)
            key = 1    
        else:
            lines = list(map(float, lines.split()))
            tbl_p.add_row(lines)
    tbl_p.sort(' Flux ')

    path = os.path.join('results')
    os.chdir(path)
    img = fits.getdata('Init_imag.fits')

    plt.subplot(1, 1, 1)
    norm = simple_norm(img, 'log', percent=99.)
    plt.imshow(img, norm=norm, origin='lower', cmap='viridis' )
    # plot init
    x_0 = tbl_init[' x_0 ']
    y_0 = tbl_init[' y_0 ']    
    plt.plot(x_0, y_0, 'r+', ms=3, label='init')

    # plot galaxies
    x_0 = tbl_g[' x_0 ']
    y_0 = tbl_g[' y_0 ']    
    plt.plot(x_0, y_0, 'o', mec='r', mfc='none', label='galaxies')


    # plot PSs
    x_0 = tbl_p[' x_0 ']
    y_0 = tbl_p[' y_0 ']    
    plt.plot(x_0, y_0, 's', mec='r', mfc='none', label='PSs')
    #plt.legend()
    plt.savefig('sources_res.png')

    tbl_g['x_err'] = tbl_init_g[' x_0 '] - tbl_g[' x_0 ']
    tbl_g['y_err'] = tbl_init_g[' y_0 '] - tbl_g[' y_0 ']
    tbl_g['f_err'] = tbl_init_g[' Flux '] - tbl_g[' Flux ']
    tbl_g['re_err'] = tbl_init_g[' re '] - tbl_g[' re ']    
    tbl_g['n_err'] = tbl_init_g[' n '] - tbl_g[' n ']
    tbl_g['g1_err'] = tbl_init_g[' g1 '] - tbl_g[' e1 ']
    tbl_g['g2_err'] = tbl_init_g[' g2 '] - tbl_g[' e2 ']
    tbl_init_g[' Flux '] = 20-2.5*np.log10(tbl_init_g[' Flux '])

    
    plt.subplot(3, 3, 1)
    plt.plot(tbl_init_g[' Flux '], tbl_g['x_err'], 'bo')
    plt.plot( [min(tbl_init_g[' Flux ']), max(tbl_init_g[' Flux '])], [0, 0], 'y--')    
    plt.xlabel('Mag')
    plt.ylabel('x_err')
    plt.title('X_err')

    plt.subplot(3, 3, 2)
    plt.plot(tbl_init_g[' Flux '], tbl_g['y_err'], 'bo')
    plt.plot( [min(tbl_init_g[' Flux ']), max(tbl_init_g[' Flux '])], [0, 0], 'y--')
    plt.xlabel('Mag')
    plt.ylabel('y_err')
    plt.title('Y_err')

    plt.subplot(3, 3, 3)
    plt.plot(tbl_init_g[' Flux '], tbl_g['f_err'], 'bo')
    plt.plot( [min(tbl_init_g[' Flux ']), max(tbl_init_g[' Flux '])], [0, 0], 'y--')
    plt.xlabel('Mag')
    plt.ylabel('f_err')
    plt.title('Flux_err')

    plt.subplot(3, 3, 4)
    plt.plot(tbl_init_g[' Flux '], tbl_g['re_err'], 'bo')
    plt.plot( [min(tbl_init_g[' Flux ']), max(tbl_init_g[' Flux '])], [0, 0], 'y--')
    plt.xlabel('Mag')
    plt.ylabel('re_err')
    plt.title('Re_err')

    plt.subplot(3, 3, 5)
    plt.plot(tbl_init_g[' Flux '], tbl_g['n_err'], 'bo')
    plt.plot( [min(tbl_init_g[' Flux ']), max(tbl_init_g[' Flux '])], [0, 0], 'y--')
    plt.xlabel('Mag')
    plt.ylabel('n_err')
    plt.title('Sersic_Index_err')
    
    plt.subplot(3, 3, 7)
    plt.plot(tbl_init_g[' Flux '], tbl_g['g1_err'], 'bo')
    plt.plot( [min(tbl_init_g[' Flux ']), max(tbl_init_g[' Flux '])], [0, 0], 'y--')
    plt.xlabel('Mag')
    plt.ylabel('g1_err')
    plt.title('g1_err')

    plt.subplot(3, 3, 8)
    plt.plot(tbl_init_g[' Flux '], tbl_g['g2_err'], 'bo')
    plt.plot( [min(tbl_init_g[' Flux ']), max(tbl_init_g[' Flux '])], [0, 0], 'y--')
    plt.xlabel('Mag')
    plt.ylabel('g2_err')
    plt.title('g2_err')

    plt.tight_layout()

    plt.savefig('g_res.png')
    
 


    tbl_p['x_err'] = tbl_init_p[' x_0 ']  - tbl_p[' x_0 ']
    tbl_p['y_err'] = tbl_init_p[' y_0 ']  - tbl_p[' y_0 ']
    tbl_p['f_err'] = tbl_init_p[' Flux '] - tbl_p[' Flux ']
    tbl_init_p[' Flux '] = 20-2.5*np.log10(tbl_init_p[' Flux '])

    plt.subplot(3, 1, 1)
    plt.plot(tbl_init_p[' Flux '], tbl_p['x_err'], 'bo')
    plt.plot( [min(tbl_init_p[' Flux ']), max(tbl_init_p[' Flux '])], [0, 0], 'y--')    
    plt.xlabel('Mag')
    plt.ylabel('x_err')
    plt.title('X_err')

    plt.subplot(3, 1, 2)
    plt.plot(tbl_init_p[' Flux '], tbl_p['y_err'], 'bo')
    plt.plot( [min(tbl_init_p[' Flux ']), max(tbl_init_p[' Flux '])], [0, 0], 'y--')
    plt.xlabel('Mag')
    plt.ylabel('y_err')
    plt.title('Y_err')

    plt.subplot(3, 1, 3)
    plt.plot(tbl_init_p[' Flux '], tbl_p['f_err'], 'bo')
    plt.plot( [min(tbl_init_p[' Flux ']), max(tbl_init_p[' Flux '])], [0, 0], 'y--')
    plt.xlabel('Mag')
    plt.ylabel('f_err')
    plt.title('Flux_err')

    plt.savefig('p_res.png')



    


   
