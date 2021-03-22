
CXX      ?= g++
OPT ?= -O2
CXXFLAGS  = ${OPT}  -std=c++17 -lstdc++fs   -pg -g  -Wall -Wextra -Ilib/cgranges -Ilib/cxxopts -Ilib/fmtlib $(LDFLAGS)
LDFLAGS   += -fopenmp -lpthread -lz -lstdc++fs
SOURCES   = main.cpp paf.cpp filters.cpp format.cpp cluster.cpp reference.cpp util.cpp chainchain.cpp cigar.cpp annotate.cpp
HEADERS = $(SOURCES:.cpp=.h)
OBJDIR = obj
OBJECTS   = $(addprefix $(OBJDIR)/, $(SOURCES:.cpp=.o))

fusion: $(OBJECTS) 
	$(CXX) $(OBJECTS) -o $@ ${LDFLAGS} -pg
$(OBJDIR)/%.o : %.cpp
	$(CXX) $(CXXFLAGS) -c -o $@  $<
clean:
	@rm -f $(OBJECTS)
clean-exe:
	@rm -f fusion
