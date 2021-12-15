
CXX     ?= g++
OPT ?= -O2
LDFLAGS   +=  -lz 
CXXFLAGS  = ${OPT}  -std=c++17  -Wall -Wextra -Ilib/cgranges -Ilib/cxxopts -Ilib/fmtlib -Ilib/stattest 
ifdef CONDA_PREFIX
	CXXFLAGS += -I${CONDA_PREFIX}/include
	LDFLAGS += -L${CONDA_PREFIX}/lib
endif
#-fopenmp -lpthread
SOURCES   = main.cpp paf.cpp filters.cpp format.cpp reference.cpp util.cpp chainchain.cpp cigar.cpp annotate.cpp
HEADERS = $(SOURCES:.cpp=.h)
OBJDIR = obj
OBJECTS   = $(addprefix $(OBJDIR)/, $(SOURCES:.cpp=.o))

genion: $(OBJECTS) directories
	$(CXX) $(OBJECTS) -o $@ ${LDFLAGS}


$(OBJDIR)/%.o : %.cpp directories
	$(CXX) $(CXXFLAGS) -c -o $@  $< 


directories: ${OBJDIR}

${OBJDIR}:
	mkdir -p ${OBJDIR}
clean:
	@rm -f $(OBJECTS)
clean-exe:
	@rm -f genion
