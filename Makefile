
BIN     =  pimcf90

OBJ     =  math.o \
           init.o \
           refbeadmodule.o \
           rtable.o \
           paircorspecies.o \
           ag_density.o \
	   observables.o \
	   energies.o \
	   action.o \
           exchangemodule.o \
	   moves.o \
           montecarlo.o \
	   main.o

LIBS = 

FFLAGS = -O2
#FFLAGS = -Wall -fbounds-check
#FFLAGS = -CB -warn all
FC = mpif90



$(BIN): $(OBJ)
	$(FC) $(FFLAGS) -o $(BIN) $(OBJ) $(LIBS)


%.o : %.f90
	$(FC) -c $(FFLAGS) $<

#.o: 
#	$(FC) $(FLAGS) -c $(.SOURCE)

clean: 
	rm -f core
	rm -f *.mod
	rm -f $(OBJ)
	rm -f $(BIN)
