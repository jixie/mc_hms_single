## This makefile must be executed with gmake (gnu make).

## These will have to be modified when setting up your own simc.  They
## point to the software necessary to run simc.

## ARGONNE DEFAULT SETUP FLAGS:
simcdir = .
#Csoft = /disk1/users/reinhold/Csoft/05Dec1996
#CERN_ROOT = /disk1/lib/cern/97a
## CEBAF DEFAULT SETUP FLAGS: ALSO NEED TO UNCOMMENT CERN_ROOT DEFNS FOR HP/SUN/...
#simcdir = .
#CERN_ROOT = /u/site/cernlib/pc_linux/99
#CERN_ROOT = /u/site/cernlib/x86_64_rhel5/2005
#Csoft = ~christy/Csoft


## THE REST SHOULD BE OK WITHOUT MODIFICATION.

RM        = rm -f 
SHELL     = /bin/sh

my_objs =  check_dipole.o check_sieve.o locforunt.o lowcase.o stringlib.o \
          ran3.o musc.o musc_ext.o project.o rotate_haxis.o transp.o \
          mt19937.o gauss1.o mc_hms.o mc_hms_recon.o mc_hms_hut.o \
          radlength_tuna.o hms_mispoint.o musc_ld.o trgInit.o \
          track_from_tgt.o track_to_tgt.o trgTrackToPlane.o \
          trgDerive.o trgField.o trgRK4Bdl.o simc_hms_recon.o \


MYOS := $(subst -,,$(shell uname))
CERNLIBS = -lgeant$(GEANTVER) -lpawlib -lgraflib -lgrafX11 -lpacklib -lmathlib

ifeq ($(MYOS),HPUX)
  CERNROOT = /site/cernlib/hp_ux102/97a
  ifneq (,$(findstring 09,$(shell uname -r)))
    HPUXVERSION := 09
  else
    HPUXVERSION := 10
  endif
  LIBROOT = $(Csoft)/../$(MYOS)$(HPUXVERSION)/lib
else
  LIBROOT = $(Csoft)/../$(MYOS)/lib
endif

ifeq ($(MYOS),HPUX)
  CERN_ROOT = /site/cernlib/hp_ux102/97a
  FFLAGS=+U77 +ppu -C +e +es +FPVZOU -O +Onolimit -R8
  LDFLAGS=-Wl,-a archive
  OTHERLIBS = \
	-Wl,-L$(CERN_ROOT)/lib -lpacklib $(CERNLIBS) \
	-Wl,-L/usr/lib/X11R5 -lX11 -lm
endif


ifeq ($(MYOS),ULTRIX)
  FFLAGS=-check_bounds
  LDFLAGS=
  OTHERLIBS = -L$(CERN_ROOT)/lib -lpacklib $(CERNLIBS)
endif

ifeq ($(MYOS),Linux)

  FC=gfortran

#  FFLAGS = -O2 -Wall -Wsurprising -ffixed-line-length-132 -C -I$(Csoft)/Analyzer/INCLUDE
  FFLAGS = -O2 -Wall -Wsurprising -ffixed-line-length-132 -C 
#  OTHERLIBS = -L$(CERN_ROOT)/lib -lpacklib -lmathlib -lpawlib -lc -lm
  CERNLIBS = -L$(CERN_ROOT)/lib -lpacklib -lmathlib
  OTHERLIBS = -L$(CERN_ROOT)/lib $(CERNLIBS) -lpawlib -lc -lm 
endif

ifeq ($(MYOS),SunOS)
  CERN_ROOT = /site/cernlib/sun4_solaris2/97a
  FFLAGS=-e -O -I$(Csoft)/SRC/INCLUDE -ftrap=underflow
  ifeq ($(OSTYPE),SunOS4)
    OTHERLIBS = -L$(CERN_ROOT)/lib $(CERNLIBS) -lnsl -lX11
  else
    OTHERLIBS = -L$(CERN_ROOT)/lib $(CERNLIBS) -lnsl -lsocket -lX11
  endif
endif

ifeq ($(MYOS),AIX)
  F77=f77
  FFLAGS=-g -qfixed=132 -qextname -O -I$(Csoft)/SRC/INCLUDE
  OTHERLIBS = -L$(CERN_ROOT)/lib -lpacklib $(CERNLIBS) -lX11
endif

ifeq ($(MYOS),OSF1)
  F77=f77
  LIBROOT = $(Csoft)/OSF1/lib
  FFLAGS= -r8 -extend_source -Wl,-taso -I -warn argument_checking \
        -warn declarations -warn truncated_source -warn unused
  LDFLAGS= 
  OTHERLIBS = -Wl,-L$(CERN_ROOT)/lib \
        -lpacklib $(CERNLIBS) -Wl,-L/usr/lib/X11R5 -lX11 -lm 
endif

%.o: %.f
	$(FC) $(FFLAGS) -c $< -o $@

mc_hms_single: $(my_objs) Makefile mc_hms_single.o
	$(FC) -o $@ $(FFLAGS) mc_hms_single.o $(my_objs) $(OTHERLIBS)

clean:
	rm -f *.o mc_hms_single





