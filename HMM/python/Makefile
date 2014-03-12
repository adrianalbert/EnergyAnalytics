NB_REPETITIONS = 10
SCENARIO = 1
WITH_C = 0 
# should be 0 for pure python

all: build install

build: HMM_C.c setup.py
	python setup.py build clean

install:
	cp build/lib*/HMM_C* .

benchmark: clean all
	python HMM_discret.py $(NB_REPETITIONS) $(SCENARIO) $(WITH_C) 

clean: 
	rm -rf build *.so *.dll *~
