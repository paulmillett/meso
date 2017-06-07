
# -------------------------------------------------------------------------------
# define variables...
# -------------------------------------------------------------------------------

BIN = meso
MPICXX = mpicxx
CXXFLAGS = -O3
LIBS = -lfftw3_mpi -lfftw3 -lm
INCL =
HDRS =
OBJS = $(SRCS:.cpp=.o)
SRCS = $(wildcard src/cahn_hilliard/finite_difference/*.cpp \
                  src/cahn_hilliard/finite_difference/*/*.cpp \
						src/cahn_hilliard/spectral/*.cpp \
						src/cahn_hilliard/spectral/*/*.cpp \
						src/brownian_dynamics/*.cpp \
						src/brownian_dynamics/*/*.cpp \
						src/lattice_boltzmann/*.cpp \
						src/phase_field/*.cpp \
						src/phase_field/*/*.cpp \
						src/particle_dynamics/*.cpp \
						src/particle_dynamics/*/*.cpp \
                  src/base/*.cpp)

						#src/fipi/*.cpp \

# -------------------------------------------------------------------------------
# rules for compiling and linking...
# -------------------------------------------------------------------------------

$(BIN): Start CompileMessage $(OBJS)
	@echo
	@echo ------------------------------------------------------------------------
	@echo LINK OBJECT FILES TO PRODUCE AN EXECUTABLE.......
	@echo ------------------------------------------------------------------------
	@echo
	$(MPICXX) $(CXXFLAGS) $(OBJS) $(LIBPATH) $(LIBS) -o $(BIN)
	@echo
	@echo ------------------------------------------------------------------------
	@echo REMOVE OBJECT FILES.......
	@echo ------------------------------------------------------------------------
	@echo
	rm $(OBJS)
	@echo

# -------------------------------------------------------------------------------
# print message indicating compile start...
# -------------------------------------------------------------------------------

Start:
	@echo
	@echo ------------------------------------------------------------------------
	@echo LOAD NEEDED MODULES.......
	@echo ------------------------------------------------------------------------
	@echo
	@echo need to work on this...
#	module load /share/apps/modules/Modules/modulefiles/fftw/3.3.4

# -------------------------------------------------------------------------------
# print message indicating compile start...
# -------------------------------------------------------------------------------

CompileMessage:
	@echo
	@echo ------------------------------------------------------------------------
	@echo COMPILE SOURCE FILES.......
	@echo ------------------------------------------------------------------------
	@echo

# -------------------------------------------------------------------------------
# pattern rule for creating an object file from source file...
# -------------------------------------------------------------------------------

.SUFFIXES:            # delete the default suffixes
.SUFFIXES: .cpp .o    # define our suffix list
.cpp.o:
	$(MPICXX) $(CXXFLAGS) $(INCL) -c -o $@ $<

# -------------------------------------------------------------------------------
# other "generic" target rules...
# -------------------------------------------------------------------------------

all:
	@echo
	@echo ------------------------------------------------------------------------
	@echo TOUCH SOURCE FILES TO UPDATE TIMESTAMP.......
	@echo ------------------------------------------------------------------------
	@echo
	touch $(SRCS)
	make

clean:
	rm -f $(OBJS) $(LIB).a core

tar:
	tar -cvf $(TARNAME).tar $(SRCS) $(HDRS)  makefile
