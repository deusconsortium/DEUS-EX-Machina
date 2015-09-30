# PRECISION should be DOUB or SING
PRECISION = SING

MAKESLICES = NOSLICES
#if MAKESLICES=SLICES then do not use C routine to write a single file, instead split all files slice by slice.  


#intel
#CC = icc
#CFLAGS = -D_LARGEFILE64_SOURCE  -D_FILE_OFFSET_BITS=64
#FC = mpif90
#FFLAGS = -I/usr/local/fftw-2.1.5-intel/include -O3 -g -traceback  -cpp  -D$(PRECISION) -D$(MAKESLICES) -convert big_endian
#ifeq ($(PRECISION),DOUB)
#	LIBS =-L/usr/local/fftw-2.1.5-intel/lib -ldrfftw_mpi -ldfftw_mpi -ldrfftw -ldfftw -limf -lm
#else
#	LIBS = -L/usr/local/fftw-2.1.5-intel/lib -lsrfftw_mpi -lsfftw_mpi -lsrfftw -lsfftw -limf -lm
#endif

#gnu
CC = gcc
CFLAGS = -D_LARGEFILE64_SOURCE  -D_FILE_OFFSET_BITS=64
FC = mpif90

# 

FFLAGS = -I/usr/local/fftw-2.1.5-gnu/include -frecord-marker=4 -O3 -ffree-line-length-none -g -fbacktrace -x f95-cpp-input -DLONGINT -D$(PRECISION) -D$(MAKESLICES) -fconvert=big-endian
ifeq ($(PRECISION),DOUB)
	LIBS =-L/usr/local/fftw-2.1.5-gnu/lib -ldrfftw_mpi -ldfftw_mpi -ldrfftw -ldfftw -limf -lm
else
	LIBS = -L/usr/local/fftw-2.1.5-gnu/lib -lsrfftw_mpi -lsfftw_mpi -lsrfftw -lsfftw -limf -lm
endif



ifeq ($(MAKESLICES),SLICES)
	OBJS = grafic_types.o common_var_ramses.o grafic_io.o modpart2extrema.o
else
	OBJS = grafic_types.o parallel_io.o common_var_ramses.o grafic_io.o modpart2extrema.o
endif
MODS = grafic_types.mod grafic_io.mod common_var_ramses.mod


all:	compute_extrema

compute_extrema:	$(OBJS) $(MODS) compute_extrema.o
	$(FC)  -limf  -o $@ compute_extrema.o $(OBJS) $(LIBS)

clean:
	rm -rf $(OBJS) $(MODS) powergrid.o *.mod powergrid compute_extrema.o compute_extrema

#.f90.o:
#	$(FC) $(FFLAGS) -c $< -o $@
%.o:	%.f90
	$(FC) $(FFLAGS) -c $< -o $@

