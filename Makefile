Q=rlwrap q
QSCRIPT=test.q

EXT = so
LONG_BIT := 32 #$(shell getconf LONG_BIT)
OPTFLAGS=-O2 -Wall -pedantic
INCLUDES=-I../../../../kx/kdb+/c/c
CFLAGS=$(OPTGLAGS) -fPIC -DKXVER=3 -Wno-strict-aliasing -std=c99 $(INCLUDES)
FFLAGS=$(OPTFLAGS) -fPIC -Wno-uninitialized -Wno-unused -ffloat-store -std=f95

CC=gcc -m$(LONG_BIT)
FORTRAN=gfortran -m$(LONG_BIT)
MINPACK_SRCDIR=minpack
OBJDIR := obj

all:shared
vpath %.f $(MINPACK_SRCDIR)

CSOURCES = $(wildcard *.c)
FSOURCES = dpmpar enorm fdjac2 lmdif lmpar qrfac qrsolv

COBJECTS := $(patsubst %.c,%.o,$(CSOURCES))
FOBJECTS := $(addsuffix .o,$(FSOURCES))
OBJECTS := $(addprefix $(OBJDIR)/,$(COBJECTS) $(FOBJECTS))
LINKFLAGS = -Wl,--version-script=nls.map

$(OBJECTS): | $(OBJDIR)

shared: $(OBJECTS)
	@echo "Making DLL"
	$(FORTRAN) $(LINKFLAGS) -shared $^ -o nls.$(EXT) -lm -lpthread

objects: $(OBJECTS)

$(OBJDIR):
	@echo "Creating obj directory"
	@mkdir -p $(OBJDIR)

$(OBJDIR)/%.o : %.f
	@echo "Making " $@ " from " $<
	@$(FORTRAN) $(FFLAGS) -c $< -o $@

$(OBJDIR)/%.o : %.c
	@echo "Compiling " $<
	@$(CC) $(CFLAGS) -c $< -o $@

clean:
	$(RM) $(OBJECTS)

test: shared
	$(Q) $(QSCRIPT)

