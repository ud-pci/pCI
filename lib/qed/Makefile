# Specify install directory
PREFIX=../../bin

# Specify include directory
INCLUDE=../include

# Specify Fortran compiler
FC=ifort

# Specify Fortran compiler flags
FFLAGS=-O -r8 -I $(INCLUDE)

# Specify executables
TARGET=qed_rot qedpot_conf sgc0

all: $(TARGET)

qed_rot: ./qed_rot.for
		$(FC) $(FFLAGS) -o qed_rot ./qed_rot.for

qedpot_conf: ./qedpot_conf.for
		$(FC) $(FFLAGS) -o qedpot_conf ./qedpot_conf.for

sgc0: ./sgc0.f
		$(FC) $(FFLAGS) -o sgc0 ./sgc0.f

install: all
	install -d $(PREFIX)
	install $(TARGET) $(PREFIX)

clean:
	$(RM) $(TARGET)