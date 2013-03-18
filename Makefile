.SUFFIXES: .f90 .o
CC = ifort -openmp
CFLAGS = -c 
CAMBLIB = /u/sudeep/CAMBLATEST/CAMB
#LAPACKLIB = /u/sudeep/LAPACK
OBJECTS = miscUtilsForCAMB.o FisherModules.o FileIO.o LensingNoise.o GetSpectra.o FisherGenerate.o
LDFLAGS =  -L$(CAMBLIB) -lcamb #$(LAPACKLIB)/lapack_LINUX3.a $(LAPACKLIB)/blas_LINUX3.a 
EXECUTABLE = FisherCombine.x
.f90.o:
	$(CC) $(CFLAGS) $< -I$(CAMBLIB) $(LDFLAGS)
default: lib $(EXECUTABLE)

FisherCombine.x: FisherCombine.f90  $(OBJECTS)
	$(CC) $(OBJECTS) $< -I$(CAMBLIB) $(LDFLAGS) -o $@
lib:$(OBJECTS)
	ar rv libfisher.a *.o
doc:

clean:
	rm -rf *.o *mod *.x
