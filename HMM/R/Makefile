NB_REPETITIONS = 10
SCENARIO = 1
WITH_C = 0
# should be 0 for pure R

all: HMM.so

%.so: %.c
	R CMD SHLIB $<

clean: 
	rm -f *.o *.so *.Rout

benchmark:
	echo "source('HMM_discret.R'); benchmark($(NB_REPETITIONS), $(SCENARIO), $(WITH_C));" |R  --no-save --slave 
