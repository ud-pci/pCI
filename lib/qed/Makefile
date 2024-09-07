# Specify install directory
PREFIX=../bin

# Specify include directory
INCLUDE=../include

# Specify Fortran compiler
FC=ifort

# Specify Fortran compiler flags
FFLAGS=-O -r8 -I $(INCLUDE)

# Specify executables
TARGET=qed_rot qedpot_conf sgc

all: $(TARGET)

qed_rot: ./qed_rot.for
		$(FC) $(FFLAGS) -o qed_rot ./qed_rot.for

qedpot_conf: ./qedpot_conf.for
		$(FC) $(FFLAGS) -o qedpot_conf ./qedpot_conf.for

sgc: ./sgc.for
		$(FC) $(FFLAGS) -o sgc ./sgc.for 

install: all
	install -d $(PREFIX)
	install $(TARGET) $(PREFIX)

clean:
	$(RM) $(TARGET)