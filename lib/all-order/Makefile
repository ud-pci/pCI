# Specify install directory
PREFIX=../../bin

# Specify Fortran compiler (override with FC=ifort on the command line)
FC ?= ifx

# Specify include directory
INCLUDE=-I./include
INCLUDE_CI=-I../include
INCLUDE_BASIS=-I./basis

# MKL flag differs between ifort (-mkl) and ifx (-qmkl)
ifeq ($(findstring ifort,$(notdir $(FC))),ifort)
    INCLUDE_LAPACK = -mkl=sequential
else
    INCLUDE_LAPACK = -qmkl=sequential
endif

# Specify Fortran compiler flags
FFLAGS=-O -r8

# Specify executables
TARGET=nspl bspl tdhf bdhf bas_wj bas_x allcore-ci sdvw-ci valsd-ci second-ci
 
all: $(TARGET)

nspl: ./basis/nspl/nspl.f
		$(FC) $(FFLAGS) $(INCLUDE_LAPACK) $(INCLUDE) $(INCLUDE_BASIS)/nspl -o nspl ./basis/nspl/nspl.f

bspl: ./basis/bspl/bspl.f
		$(FC) $(FFLAGS) $(INCLUDE_LAPACK) $(INCLUDE) $(INCLUDE_BASIS)/bspl -o bspl ./basis/bspl/bspl.f

tdhf: ./basis/tdhf/tdhf.f
		$(FC) $(FFLAGS) $(INCLUDE_LAPACK) $(INCLUDE) $(INCLUDE_BASIS)/tdhf -o tdhf ./basis/tdhf/tdhf.f

bdhf: ./basis/bdhf/bdhf.f
		$(FC) $(FFLAGS) $(INCLUDE_LAPACK) $(INCLUDE) $(INCLUDE_BASIS)/bdhf -o bdhf ./basis/bdhf/bdhf.f

bas_wj: ./bas_wj.for
		$(FC) $(FFLAGS) $(INCLUDE) $(INCLUDE_CI) -o bas_wj ./bas_wj.for

bas_x: ./bas_x.for
		$(FC) $(FFLAGS) $(INCLUDE) $(INCLUDE_CI) -o bas_x ./bas_x.for

allcore-ci: ./allcore-ci.f
		$(FC) $(FFLAGS) $(INCLUDE) -mcmodel=large -o allcore-ci ./allcore-ci.f

sdvw-ci: ./sdvw-ci.f
		$(FC) $(FFLAGS) $(INCLUDE) -mcmodel=large -o sdvw-ci ./sdvw-ci.f

valsd-ci: ./valsd-ci.f
		$(FC) $(FFLAGS) $(INCLUDE) -mcmodel=large -o valsd-ci ./valsd-ci.f

second-ci: ./second-ci.f
		$(FC) $(FFLAGS) $(INCLUDE) $(INCLUDE_CI) -mcmodel=large -o second-ci ./second-ci.f

install: all
	install -d $(PREFIX)
	install $(TARGET) $(PREFIX)

clean:
	$(RM) $(TARGET) *.mod *genmod.f90
