#
#-I/apps/plpa/plpa-1.0.3/include 
#-L/apps/plpa/plpa-1.0.3/lib 
#-Wl,-rpath,/apps/plpa/plpa-1.0.3/lib -lplpa

#
######################### general section ########################
#

#CC: Program for compiling C programs
RZ_CC ?= g++
RZ_MPICC ?= g++

#CXX: Program for compiling C++ programs
RZ_CXX ?= g++
RZ_MPICXX ?= g++

#CPP: Program for running the C preprocessor, with results to standard output
RZ_CPP ?= cpp

#FC: Program for compiling or preprocessing Fortran and Ratfor programs
RZ_FC ?= ifort
RZ_MPIFC ?= ifort

#CFLAGS: Extra flags to give to the C compiler.
RZ_CFLAGS ?= -Ofast -Wall -pedantic -ansi 

#CXXFLAGS: Extra flags to give to the C++ compiler.
RZ_CXXFLAGS ?= -Ofast -Wall -pedantic -ansi 

#CPPFLAGS: Extra flags to give to the C preprocessor and programs that use it (the C and Fortran compilers).
RZ_CPPFLAGS ?= -P

#FFLAGS: Extra flags to give to the Fortran compiler.
RZ_FFLAGS ?= -Ofast -Wall -pedantic -ansi 

#LDFLAGS: Extra flags to give to compilers when linking
RZ_LDFLAGS ?= -O3 -Ofast -IPF-fp-relaxed -ansi 

RZ_FNNOALIAS ?= -fno-fnalias

#CPPSYMBOL: how to force Fortran/C preprocessing
RZ_FPP ?= -cpp
RZ_CPP ?= UNKNOWN

#OMP: how to switch on OpenMP for C and Fortran
RZ_COMP ?= #-fopenmp
RZ_FOMP ?= #-fopenmp

#BLAS: how to link BLAS
RZ_BLAS ?= -L/opt/intel/mkl80_22/lib/64 -lmkl_ipf

#BLASLAPACK: how to link with LAPACK and BLAS
RZ_BLASLAPACK ?= -L/opt/intel/mkl80_22/lib/64 -lmkl_lapack -lmkl_ipf -lguide -lpthread

#MATH library
RZ_MATHLIB ?=

#FFTW: how to link with FFTW3
RZ_FFTW3 ?= -L/opt/fftw-3.1.1/lib -lfftw3
RZ_FFTW3INC ?= -I/opt/fftw-3.1.1/include

#MODDIR: how to deal with F90 modulefiles
RZ_MODDIR ?= 

#MPILIB: how to specify MPI libs (if required)
RZ_MPILIB ?= -lmpi
RZ_MPICXXLIB ?= -lmpi++

#MPIINC: how to specify MPI include files
RZ_MPIINC ?= 

RZ_FSYMBOL ?= -D
RZ_CSYMBOL ?= -D
RZ_CPPSYMBOL ?= -D

