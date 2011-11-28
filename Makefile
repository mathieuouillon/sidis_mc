NAME = simulation

FFLAGS  =  -O3 -w -g
CPPFLAGS  =  -g -O3

ROOTLIBS = $(shell root-config --libs)
CPPFLAGS += $(shell root-config --cflags)

# add BOS.f for BOS file output
SRC=    simulation.f fermimotion.f pythia-6.4.22.f \
        pythiaconfig.f book.f transfo.f \
        fermimotion2.f nucdens.f density.f qweight.f \
        main.cc

OBJ=    simulation.o fermimotion.o pythia-6.4.22.o \
        pythiaconfig.o book.o transfo.o \
        fermimotion2.o nucdens.o density.o qweight.o \
        main.o

.f.o:
	gfortran -c $(FFLAGS) -o $@ $*.f
.c.o:
	g++ $(CPPFLAGS) -c $*.cc

go: ${OBJ}
	g++ $(CPPFLAGS) -o $(NAME) $(OBJ)  -lgfortran -L$(CERN_ROOT)/lib -lpawlib -lpacklib $(ROOTLIBS)
#	gfortran  $(FFLAGS) -o $(NAME) $(OBJ)  -L$(CERN_ROOT)/lib -lpawlib -lpacklib \
                                               -L$(CLAS_LIB) -lbosio -lbos -lfpack -lc_bos_io -lrecutl

clean:
	rm -f $(NAME) simulation.o fermimotion.o pythiaconfig.o book.o \
        fermimotion2.o nucdens.o density.o qweight.o transfo.o \
        fort.9 last.kumac last.kumacold paw.metafile \
	main.o root.o

