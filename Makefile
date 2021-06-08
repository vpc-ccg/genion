
CXX      ?= g++
OPT ?= -O2
CXXFLAGS  = ${OPT}  -std=c++17 -lstdc++fs   -pg -g  -Wall -Wextra -Ilib/cgranges -Ilib/cxxopts -Ilib/fmtlib $(LDFLAGS)
LDFLAGS   +=  -lz -lstdc++fs 
#-fopenmp -lpthread
SOURCES   = main.cpp paf.cpp filters.cpp format.cpp cluster.cpp reference.cpp util.cpp chainchain.cpp cigar.cpp annotate.cpp
HEADERS = $(SOURCES:.cpp=.h)
OBJDIR = obj
OBJECTS   = $(addprefix $(OBJDIR)/, $(SOURCES:.cpp=.o))

fusion: $(OBJECTS) directories
	$(CXX) $(OBJECTS) -o $@ ${LDFLAGS} -pg


$(OBJDIR)/%.o : %.cpp directories
	$(CXX) $(CXXFLAGS) -c -o $@  $<


directories: ${OBJDIR}

${OBJDIR}:
	mkdir -p ${OBJDIR}
clean:
	@rm -f $(OBJECTS)
clean-exe:
	@rm -f fusion
