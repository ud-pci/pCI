# Specify install directory
PREFIX=../../bin

# Specify include directory
INCLUDE=-I./include
INCLUDE_CI=-I../include
INCLUDE_BASIS=-I./basis
INCLUDE_LAPACK=-mkl=sequential

# Specify Fortran compiler
FC=ifort

# Specify Fortran compiler flags
FFLAGS=-O -r8

# Specify executables
TARGET=nspl bspl tdhf bdhf bas_wj bas_x allcore-cis sdvw-cis valsd-cis second-cis
 
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

allcore-cis: ./allcore-cis.f
		$(FC) $(FFLAGS) $(INCLUDE) -mcmodel=large -o allcore-cis ./allcore-cis.f

sdvw-cis: ./sdvw-cis.f
		$(FC) $(FFLAGS) $(INCLUDE) -mcmodel=large -o sdvw-cis ./sdvw-cis.f

valsd-cis: ./valsd-cis.f
		$(FC) $(FFLAGS) $(INCLUDE) -mcmodel=large -o valsd-cis ./valsd-cis.f

second-cis: ./second-cis.f
		$(FC) $(FFLAGS) $(INCLUDE) $(INCLUDE_CI) -mcmodel=large -o second-cis ./second-cis.f

install: all
	install -d $(PREFIX)
	install $(TARGET) $(PREFIX)

clean:
	$(RM) $(TARGET) *.mod *genmod.f90
