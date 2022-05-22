res: extr
	mv psf_sources.psf default.psf; mv *_psf_sources.ps results; cp default.psf default.fits; ./psf_viwer.py

extr: cat
	psfex psf_sources.cat -c extr.sex; mv *_psf_sources.fits results


cat: results
	sex sources.fits -c psfex.sex; mv psfex_*.fits results

 
results:
	rm -fr results
	mkdir results
	
