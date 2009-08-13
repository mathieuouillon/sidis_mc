NAME = simulation

FFLAGS  =  -w -g

# add BOS.f for BOS file output
SRC=    simulation.f fermimotion.f pythia-6.4.19.f \
        leptoconfig.f book.f transfo.f \
        fermimotion2.f nucdens.f density.f qweight.f

OBJ=    simulation.o fermimotion.o pythia-6.4.19.o \
        leptoconfig.o book.o transfo.o \
        fermimotion2.o nucdens.o density.o qweight.o

.f.o:
	gfortran -c $(FFLAGS) -o $@ $*.f

go: ${OBJ}
	gfortran  $(FFLAGS) -o $(NAME) $(OBJ)  -L/cern/pro/lib -lpawlib -lpacklib -lkernlib -lmathlib
#	gfortran  $(FFLAGS) -o $(NAME) $(OBJ)  -L$(CERN_ROOT)/lib -lpawlib -lpacklib -lkernlib -lmathlib \
                                               -L$(CLAS_LIB) -lbosio -lbos -lfpack -lc_bos_io -lrecutl

