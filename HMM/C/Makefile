all: HMM_discret.c
	gcc -lm -O2 HMM_discret.c -o HMM_discret.out
clean: 
	rm -f *.out

benchmark: clean all
	./HMM_discret.out

