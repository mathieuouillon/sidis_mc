NAME = simulation

FFLAGS  =  -w -g

SRC=	simulation.f fermimotion.f pythia-6.4.19.f \
	leptoconfig.f book.f transfo.f
#	leptoconfig.f book.f transfo.f BOS.f
        
OBJ=	simulation.o fermimotion.o pythia-6.4.19.o \
	leptoconfig.o book.o transfo.o
#	leptoconfig.o book.o transfo.o BOS.o

.f.o:
	gfortran -c $(FFLAGS) -o $@ $*.f

go: ${OBJ}  
	gfortran  $(FFLAGS) -o $(NAME) $(OBJ)  -L/cern/pro/lib -lpawlib -lpacklib -lkernlib -lmathlib
#	gfortran  $(FFLAGS) -o $(NAME) $(OBJ)  -L$(CERN_ROOT)/lib -lpawlib -lpacklib -lkernlib -lmathlib \
                                               -L$(CLAS_LIB) -lbosio -lbos -lfpack -lc_bos_io -lrecutl
#	gfortran  $(FFLAGS) -o $(NAME) $(OBJ)  -L/usr/lib/cernlib/2006/lib -lpawlib -lpacklib -lkernlib -lmathlib

#	f77  $(FFLAGS) -o go $(OBJ)  `cernlib`
clean:
	rm -f *.o 
