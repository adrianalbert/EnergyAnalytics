NB_REPETITIONS = 10
SCENARIO = 2

clean:
	rm -f *~

benchmark:
	echo "exec('HMM_discret.sci'); benchmark($(NB_REPETITIONS), $(SCENARIO), 1);" > tmp.sci; scilab -nw -nb -f tmp.sci; rm tmp.sci
