default: all

# define compiler
CC=g++

# flags to compiler
CCFLAGS=-c -Wall
LDFLAGS=

# folders with header files
INCLUDE+= -I./include

# include also these libraries
#LIBS+=${PETSC_SYS_LIB}

# define what to compile
SOURCES=pascinference.cpp $(wildcard src/*.cpp) $(wildcard src/*/*.cpp) $(wildcard src/*/*/*.cpp)

# directory where the .o files will be stored
BUILDDIR=bin/

# define final executable file
EXECUTABLE=pascinference

# the list of all objects created from sources
_OBJECTS=$(SOURCES:.cpp=.o) 
OBJECTS=$(patsubst %,$(BUILDDIR)%,$(notdir $(_OBJECTS)))

# include PETSc
include ${PETSC_DIR}/conf/variables
include ${PETSC_DIR}/conf/rules

# include Permon
include ${FLLOP_DIR}/conf/fllop_variables
LIB+=-lm ${FLLOP_LIB}

# make all
all: $(SOURCES) $(EXECUTABLE) chkopts
	mkdir -p bin
	mkdir -p output
    
# make main file
$(EXECUTABLE): $(_OBJECTS) 
	$(CC) $(LDFLAGS) $(OBJECTS) ${LIB} -o $@ ${PETSC_SYS_LIB}

%.o: %.cpp
	${PETSC_COMPILE} $(CCFLAGS) ${INCLUDE} $< -o $(BUILDDIR)$(notdir $@) 
	
# to clean the last build
clean::
	@rm -f $(BUILDDIR)*.o
	@rm -f $(EXECUTABLE)
