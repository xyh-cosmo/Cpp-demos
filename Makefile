#   ============================================================================
#   Youhua Xu Dec-13-2017
#   ============================================================================
#   Make a build dir for compilation
MDIR := $(shell pwd)
WRKDIR = $(MDIR)/build
BINDIR = $(MDIR)/bin
# DEBUG  = $(MDIR)/debug

#	make build & binary dirs
.base:
	if ! [ -e $(WRKDIR) ]; then mkdir $(WRKDIR) ; fi;
	touch $(WRKDIR)/.base;
	if ! [ -e $(BINDIR) ]; then mkdir $(BINDIR) ; fi;
	touch $(BINDIR)/.base;
	# if ! [ -e $(DEBUG) ]; then mkdir $(DEBUG) ; fi;
	# touch $(DEBUG)/.base;
#   ============================================================================
#   Set the source file path
vpath %.cpp main:src:test:area_stats
# vpath %.cc  main:src:test:area_stats
# vpath %.c   main:src:test:area_stats
vpath %.hpp include:area_stats
# vpath %.h   include:area_stats
vpath %.o build
vpath .base build

INCLUDES = -I $(MDIR)/include

# GCC			= gcc
CPP			= g++
# MCC         = mpicc
CCFLAG  	= -Wall
OPTFLAG		= #-O2 #-ffast-math #( not recognized by intel compiler )
# OPTFLAG += -pg

LIBS	= -lgsl -lgslcblas -lm

LDFLAG      = 

#   http://www.tuicool.com/articles/NBfeYj
ifeq ($(shell uname -s),Linux)
	LDFLAG	+= -Wl,--no-as-needed
endif

ifeq ($(shell uname -s),Darwin)
	LDFLAG	+= -framework Accelerate #(-framework Accelerate is for Mac OSX)
endif

%.o: %.cpp .base
	cd $(WRKDIR);$(CPP) $(OMPFLAG) $(OPTFLAG) $(CCFLAG) $(INCLUDES) -c ../$< -o $*.o

# %.o: %.cc .base
# 	cd $(WRKDIR);$(CC) $(OMPFLAG) $(OPTFLAG) $(CCFLAG) $(INCLUDES) -c ../$< -o $*.o

# %.o: %.c .base
# 	cd $(WRKDIR);$(CC) $(OMPFLAG) $(OPTFLAG) $(CCFLAG) $(INCLUDES) -c ../$< -o $*.o

VERBOSE = Verbose.o

MAIN 	= main_4th_RK.o
RK4TH 	= ODE_solver_4th_RK.o
MISC	= misc.o

all: demo_4th_RK

demo_4th_RK: ${MAIN} ${RK4TH} ${MISC}
	${CPP} ${OPTFLAG} ${LDFLAG} $(addprefix build/,$(notdir $^)) ${LIBS} -o $(BINDIR)/$@


.PHONY:clean tidy run
clean: .base
	rm -rf $(WRKDIR);
tidy:
	make clean; rm -rf $(BINDIR) $(WRKDIR)
