NB_REPETITIONS = 10
SCENARIO = 1
WITH_C = 0
# should be 0 for pure matlab/octave

MYFLAGS = COPTIMFLAGS='-O2' -v

all: matlab octave

matlab: HMMfilter_C.mex HMMsmoother_C.mex

octave: HMMfilter_C.oct HMMsmoother_C.oct

%.oct: %.c
	mkoctfile --mex $<

%.mex: %.c
	mex $(MYFLAGS) $<

benchmark:
	matlab -nodisplay -nosplash -r "benchmark($(NB_REPETITIONS), $(SCENARIO), $(WITH_C), 1)"

octave_benchmark:
	octave --eval "benchmark($(NB_REPETITIONS), $(SCENARIO), $(WITH_C), 1)"

statbox_benchmark:
	matlab -nosplash -nodisplay < benchmarkStatBox.m

clean:
	rm -f *.mex* *.o *.oct *~
