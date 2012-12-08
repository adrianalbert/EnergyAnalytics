# @(#)makefile	2-99
# $Revision:$

.SUFFIXES: .o .C .c

# FLAGS=  -O2
FLAGS= -g

# BINDIR = /homes/daniel/bin
# SRCDIR = /homes/daniel/src
SRCDIR = .
BINDIR = .

HEADER = utils.h  obs.h obsSeq.h stateTrans.h plainStateTrans.h gammaProb.h explicitDurationTrans.h obsProb.h discreteObsProb.h  gaussianObsProb.h vectorObsProb.h initStateProb.h hmm.h

OBJ = utils.o readConfigFile.o obsSeq.o plainStateTrans.o gammaProb.o explicitDurationTrans.o initStateProb.o hmm.o discreteObsProb.o gaussianObsProb.o vectorObsProb.o 

OBJ3 = hmmFind.o

#===============================================================================


hmmFind : $(HEADER) $(OBJ) $(OBJ3)
	   g++ $(FLAGS) -o $(BINDIR)/hmmFind $(OBJ) $(OBJ3) -lm

hmmPur : $(HEADER) $(OBJ) $(OBJ3)
	   purify g++ $(FLAGS) -o $(BINDIR)/hmmPur $(OBJ) $(OBJ3) -lm


#===============================================================================

#===============================================================================

utils.o: utils.h

obsSeq.o: obs.h obsSeq.h

plainStateTrans.o: stateTrans.h plainStateTrans.h

gammaProb.o: gammaProb.h

explicitDurationTrans.o: stateTrans.h gammaProb.h explicitDurationTrans.h 

discreteObsProb.o: obs.h obsProb.h discreteObsProb.h

gaussianObsProb.o: obs.h obsProb.h gaussianObsProb.h

vectorObsProb.o: obs.h obsProb.h gaussianObsProb.h vectorObsProb.h

initStateProb.o: initStateProb.h

hmm.o: utils.h obs.h obsSeq.h stateTrans.h plainStateTrans.h gammaProb.h explicitDurationTrans.h obsProb.h discreteObsProb.h gaussianObsProb.h initStateProb.h hmm.h

hmmFind.o: utils.h obs.h obsSeq.h stateTrans.h plainStateTrans.h gammaProb.h explicitDurationTrans.h obsProb.h discreteObsProb.h gaussianObsProb.h initStateProb.h hmm.h

#===============================================================================

#===============================================================================

.C.o: $(HEADERI)
	g++ -c  $(FLAGS) $<

#checkout:
#	co -l $(REV) $(ALL)

#checkin:
#	ci $(REV)  $(ALL)

clean:
	rm -f *.o core tags

tags: $(SRC)
	ctags *.c

lint:	rm -f lint.out
	lint *.c > lint.out

tar.all: $(SRC) $(DOC) $(HEADERS) makefile
	tar cvf tar.all $(SRC) $(DOC) $(HEADERS) makefile

#===============================================================================
#===============================================================================
