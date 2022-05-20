#! /usr/bin/env python3

from funs import *
from astropy.io import fits
import os


 
if __name__ == '__main__':
    
    path = os.path.join('images','sources.fits')
    pathout = os.path.join('results')
    data, hdr = fits.getdata(path, header=True)
    data = data[0]
    cat = main(data, hdr, pathout)

    p_file = open('point.cat','w')
    g_file = open('galaxies.cat', 'w')

    print('N_source', 'x_0', 'y_0', 'Flux', 're', 'n', 'e1', 'e2', sep=' | ' , end=' | \n', file=g_file)
    print('N_source', 'x_0', 'y_0', 'Flux', sep=' | ' , end=' | \n', file=p_file)



    p = g = 0
    for i in range(len(cat)):
        if len(cat[i]) <=3:
            x, y, flux = cat[i].getParams()
            print("%9d %5.3f %5.3f %6.3f " % (p, x, y, flux), file=p_file)
            p += 1
        else:
            x0, y0, fx, re, e1, e2, n = cat[i].getParams()
            print("%9d %5.3f %5.3f %6.3f %4.1f %3.1f %3.2f %3.2f" % (g, x0, y0, fx, re*0.8452134572561297, n, e1, e2), file=g_file)
            g += 1
    
