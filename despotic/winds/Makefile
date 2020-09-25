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
     ifeq ($(UNAMEN), morticia.anu.edu.au)
          include Make.mach.morticia
     else ifeq ($(UNAMEN), avatar)
          include Make.mach.avatar
     else ifneq (,$(findstring gadi, $(UNAMEN)))
          include Make.mach.icc-gadi
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
CC		= $(MACH_CC)

# Set optimization mode flags
CXXOPTFLAGS	= $(MACH_CXXOPTFLAGS) $(MACH_OMPFLAGS) -DNDEBUG \
	-DHAVE_INLINE -MMD -MP
LDOPTFLAGS	= $(MACH_LDOPTFLAGS) $(MACH_OMPFLAGS) $(MACH_CXXFLAG)

# Set debug mode flags
CXXDEBFLAGS     = $(MACH_CXXDEBFLAGS) $(MACH_OMPFLAGS) -MMD -MP
LDDEBFLAGS	= $(MACH_LDDEBFLAGS) $(MACH_OMPFLAGS) $(MACH_CXXFLAG)

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

table:
	$(CC) -o hot_wind_tab hot_wind_tab.c $(CXXOPTFLAGS) $(LDOPTFLAGS) \
		$(INCFLAGS) $(LDLIBFLAGS)

table-debug:
	$(CC) -o hot_wind_tab hot_wind_tab.c $(CXXDEBFLAGS) $(LDDEBFLAGS) \
		$(INCFLAGS) $(LDLIBFLAGS)

