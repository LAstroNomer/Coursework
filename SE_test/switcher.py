#! /usr/bin/env python3

from astropy.table import Table
import numpy as np 
if __name__ == '__main__':
    cat = open('SE_test.cat', 'r')   
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
    #print(hdr)
    cx = cat_tbl['XMODEL_IMAGE']
    cy = cat_tbl['YMODEL_IMAGE']
    cgf = cat_tbl['FLUX_SPHEROID']
    ctp = cat_tbl['SPREAD_MODEL']
    ctp_err = cat_tbl['SPREADERR_MODEL']
    cn = cat_tbl['SPHEROID_SERSICN']
    cref = cat_tbl['SPHEROID_REFF_IMAGE']
    ce1 = cat_tbl['ELLIP1MODEL_IMAGE']
    ce2 = cat_tbl['ELLIP2MODEL_IMAGE']

    eps = 5.0 * 10**(-3)
    k = 4
    switch = np.sqrt(eps**2 + (k * ctp_err)**2)
   
    p_file = open('point.cat','w')
    g_file = open('galaxies.cat', 'w')

    print('N_source', 'x_0', 'y_0', 'Flux', 're', 'n', 'e1', 'e2', sep=' | ' , end=' | \n', file=g_file)
    print('N_source', 'x_0', 'y_0', 'Flux', sep=' | ' , end=' | \n', file=p_file)

    p = g = 0
    for x, y, gf, re, e1, e2, n, tp, sw in zip(cx-1, cy-1, cgf, cref, ce1, ce2, cn, ctp, switch):
        if tp < sw:
            print("%9d %5.3f %5.3f %6.3f " % (p, x, y, gf), file=p_file)
            p += 1
        else:
            print("%9d %5.3f %5.3f %6.3f %4.1f %3.1f %3.2f %3.2f" % (g, x, y, gf, re, n, e1, e2), file=g_file)
            g += 1
    





