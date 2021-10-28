NAME = simulation

FFLAGS  =  -O3 -w -g
CPPFLAGS  =  -g -O3

ROOTLIBS = $(shell root-config --libs)
CPPFLAGS += $(shell root-config --cflags)

# add BOS.f for BOS file output
SRC=    simulation.f fermimotion.f pythia6428.f \
        pythiaconfig.f book.f transfo.f \
        fermimotion2.f nucdens.f density.f qweight.f \
        accep_fun.f clas_at12g.f clas12_accept.f read_par_clas12g.f \
        smear_fun.f ALERT_fastMC.cc main.cc

OBJ=    simulation.o fermimotion.o pythia6428.o \
        pythiaconfig.o book.o transfo.o \
        fermimotion2.o nucdens.o density.o qweight.o \
        accep_fun.o clas_at12g.o clas12_accept.o read_par_clas12g.o \
        smear_fun.o ALERT_fastMC.o main.o

.f.o:
	gfortran -c $(FFLAGS) -o $@ $*.f
.c.o:
	g++ $(CPPFLAGS) -c $*.cc

go: ${OBJ}
	g++ $(CPPFLAGS) -o $(NAME) $(OBJ)  -lgfortran -L$(CERN_ROOT)/lib -lmathlib -lpawlib -lpacklib -lkernlib $(ROOTLIBS)
#	gfortran  $(FFLAGS) -o $(NAME) $(OBJ)  -L$(CERN_ROOT)/lib -lpawlib -lpacklib -lkernlib \
                                               -L$(CLAS_LIB) -lbosio -lbos -lfpack -lc_bos_io -lrecutl

clean:
	rm -f $(NAME) simulation.o fermimotion.o pythiaconfig.o book.o \
        fermimotion2.o nucdens.o density.o qweight.o transfo.o \
        accep_fun.o clas_at12g.o clas12_accept.o read_par_clas12g.o \
        smear_fun.o eloss.o fort.9 last.kumac last.kumacold paw.metafile \
	ALERT_fastMC.o main.o root.o

