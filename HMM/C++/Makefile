all: HMM_discret.cc
	g++ -O2 -I/home/garivier/lib/usr/include -L/home/garivier/lib/usr/lib64/ -larmadillo HMM_discret.cc -o HMM_discret.out
clean: 
	rm -f *.out
benchmark: clean all
	./HMM_discret.out

