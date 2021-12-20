
CXX     ?= g++
OPT		?= -O2
LDFLAGS   +=  -lz 
CXXFLAGS  = ${OPT}  -std=c++17   -Wall -Wextra -Ilib/cgranges -Ilib/cxxopts -Ilib/fmtlib -Ilib/stattest 
ifdef CONDA_PREFIX
	CXXFLAGS += -I${CONDA_PREFIX}/include
	LDFLAGS += -L${CONDA_PREFIX}/lib
endif
#-fopenmp -lpthread
SRCDIR = src
SOURCES_NAMES   = main.cpp paf.cpp filters.cpp format.cpp  util.cpp chainchain.cpp cigar.cpp annotate.cpp
SOURCES = $(addprefix $(SRCDIR)/, $(SOURCES_NAMES))
HEADERS = $(SOURCES:.cpp=.h)
OBJDIR = obj
OBJECTS   = $(addprefix $(OBJDIR)/, $(SOURCES_NAMES:.cpp=.o))

genion: $(OBJECTS) directories
	$(CXX) $(OBJECTS) -o $@ ${LDFLAGS}


$(OBJDIR)/%.o : $(SRCDIR)/%.cpp directories
	$(CXX) $(CXXFLAGS) -c -o $@  $< 


directories: ${OBJDIR}

${OBJDIR}:
	mkdir -p ${OBJDIR}

.PHONY: clean
clean: clean-exe
	@rm -f $(OBJECTS)
.PHONY: clean-exe
clean-exe:
	@rm -f genion
