# Main makefile

# Did the user tell us to use a particular machine? If so, use that.
ifdef MACHINE
     include Make.mach.$(MACHINE)
else
     # Machine not specified, so try to guess
     UNAME		= $(shell uname)
     UNAMEN		= $(shell uname -n)

     # Do we have a makefile that matches the machine name? If so, use
     # that. If not, use a generic makefile depending on the OS.
     ifeq ($(UNAMEN), hyades.ucsc.edu)
          include Make.mach.ucsc-hyades
     else ifeq ($(UNAMEN), cosma-e)
          include Make.mach.icc-cosma4
     else ifeq ($(UNAMEN), mosura)
          include Make.mach.coala
     else ifneq (,$(findstring raijin, $(UNAMEN)))
          include Make.mach.icc-raijin
     else ifeq ($(UNAME), Linux)
          include Make.mach.linux-gnu
     else ifeq ($(UNAME), Darwin)
          include Make.mach.darwin
     else
          $(info Cannot detect system type. Suggest you specify MACHINE= manually.)
          include Make.mach.generic
     endif
endif

# Set compiler
CXX		= $(MACH_CXX) $(MACH_C11FLAG)

# Set optimization mode flags
CXXOPTFLAGS	= $(MACH_CXXOPTFLAGS) -DNDEBUG \
	-DHAVE_INLINE -MMD -MP
LDOPTFLAGS	= $(MACH_LDOPTFLAGS) $(MACH_CXXFLAG)

# Set debug mode flags
CXXDEBFLAGS     = $(MACH_CXXDEBFLAGS) -MMD -MP
LDDEBFLAGS	= $(MACH_LDDEBFLAGS) $(MACH_CXXFLAG)

# Read any user overrides
-include Make.config.override

# Include flags
ifdef C_INCLUDE_PATH
     INCFLAGS += -I$(subst :, -I ,$(C_INCLUDE_PATH))
endif
ifdef CXX_INCLUDE_PATH
     INCFLAGS += -I$(subst :, -I ,$(CXX_INCLUDE_PATH))
endif
ifdef GSL_HDR_PATH
     INCFLAGS += -I$(GSL_HDR_PATH)
endif

# Link flags
LDLIBFLAGS      = -lgsl -lgslcblas
ifdef LD_LIBRARY_PATH
     LDLIBFLAGS += -L$(subst :, -L ,$(LD_LIBRARY_PATH))
endif
ifdef LIBRARY_PATH
     LDLIBFLAGS += -L$(subst :, -L ,$(LIBRARY_PATH))
endif
ifdef GSL_LIB_PATH
     LDLIBFLAGS += -L$(GSL_LIB_PATH)
endif

CXXFLAGS +=  $(INCFLAGS) $(DEFINES)
LDFLAGS  +=  $(LDLIBFLAGS)

# Name for compiled library, shared or static
LIBNAME		= libpwind$(LIB_EXTENSION)
LIBSTATNAME	= libpwind$(LIBSTAT_EXTENSION)

# Pointers to sources
SOURCES		= $(wildcard *.cpp)
OBJECTS		= $(SOURCES:%.cpp=%.o)
DEPS		= $(SOURCES:%.cpp=%.d)

# Default target
.PHONY: lib lib-debug libstatic libstatic-debug compile clean

lib: CXXFLAGS += $(CXXOPTFLAGS) $(CLIBFLAGS)
lib: LDFLAGS  += $(LDOPTFLAGS)
lib: $(LIBNAME)

lib-debug: CXXFLAGS += $(CXXDEBFLAGS) $(CLIBFLAGS)
lib-debug: LDFLAGS  += $(LDDEBFLAGS)
lib-debug: $(LIBNAME)

libstatic: CXXFLAGS += $(CXXOPTFLAGS) $(CLIBFLAGS)
libstatic: $(LIBSTATNAME)

libstatic-debug: CXXFLAGS += $(CXXDEBFLAGS) $(CLIBFLAGS)
libstatic-debug: $(LIBSTATNAME)

all: exe
all: lib

# Include dependencies
-include $(DEPS)

$(LIBNAME): $(OBJECTS)
	$(CXX) -o $(LIBNAME) $^ $(LDFLAGS) $(DYNLIBFLAG)

$(LIBSTATNAME): $(OBJECTS)
	ar rcs $(LIBSTATNAME) $^

compile:	$(OBJECTS)

clean:
	rm -f $(EXENAME) $(OBJECTS) $(DEPS) $(LIBNAME)
