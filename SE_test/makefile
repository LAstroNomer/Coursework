res: plot_tr
	mv T*.fits results; mv *.png results; echo "All good!"

plot_tr: trac
	./tr_resulter.py

trac: plot_sex
	./tractorenok.py

plot_sex: switch
	./resulter.py

switch: sex
	./switcher.py

sex: results
	sex sources.fits -c SE_test.sex; mv sex_*.fits results

results:
	rm -fr results
	mkdir results
        
clear:
	rm -f Tr*.fits, *.png 
