# Specify install directory
PREFIX=../bin

# Specify include directory
INCLUDE=../include

# Specify Fortran compiler
FC=ifort

# Specify Fortran compiler flags
FFLAGS=-O -r8 -I $(INCLUDE)

# Specify executables
TARGET=rpa rpa_dtm

all: $(TARGET)

rpa: ./rpa.for
		$(FC) $(FFLAGS) -o rpa ./rpa.for

rpa_dtm: ./rpa_dtm.for
		$(FC) $(FFLAGS) -o rpa_dtm ./rpa_dtm.for

install: all
	install -d $(PREFIX)
	install $(TARGET) $(PREFIX)

clean:
	$(RM) $(TARGET)