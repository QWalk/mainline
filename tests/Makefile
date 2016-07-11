all: h.results h2.results

clean: 
	rm -f *.results

h.results: 
	cd h; ./run_test | tee ../h.results

h2.results:
	cd h2; ./run_test | tee ../h2.results
