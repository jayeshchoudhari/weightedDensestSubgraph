#DEBUG=-O2
DEBUG=-O4
CPP=g++
# CFLAGS=-Wall -fmessage-length=0 -g -D DEBUG -std=c++0x  -Wextra -pedantic -pedantic-errors $(DEBUG)
CFLAGS=-Wall -g -D DEBUG -std=c++11  -Wextra $(DEBUG)

LDFLAGS=

#

# SOURCES.
#SOURCES=DynGraph.cpp UDynGraph.cpp DynGraphUtils.cpp DSAlgs.cpp DynDSAlg.cpp DynDSAlgAddRem.cpp GraphScheduler.cpp Stats.cpp
#BINARY_SOURCES=ExecAdd.cpp ExecAddRem.cpp 
SOURCES=Hypergraph.cpp ApproximateDS.cpp GraphScheduler.cpp Utility.cpp Stats.cpp ApproximateDSFullyDyn.cpp Decompose.cpp
BINARY_SOURCES=Add.cpp AddRemove.cpp Approximate.cpp





# OBJECTS.
OBJECTS=$(SOURCES:.cpp=.o)
BINARY_OBJECTS=$(BINARY_SOURCES:.cpp=.o)
ALL_LOCAL_OBJECTS=$(OBJECTS) $(BINARY_OBJECTS)

# DEPENDENCIES.
DEPS=$(patsubst %.o,%.d,$(ALL_LOCAL_OBJECTS))

# BINARIES.
BINARIES=$(BINARY_SOURCES:.cpp=)

# RULES.

all: $(BINARIES) $(OBJECTS)

$(ALL_LOCAL_OBJECTS): %.o: %.cpp
	$(CPP) -MMD -MP $(CFLAGS) -c $< -o $@
	@sed -i -e '1s,\($*\)\.o[ :]*,\1.o $*.d: ,' $*.d

$(BINARIES): %: %.o $(OBJECTS)
	$(CPP) $^ $(LDFLAGS) -o $@

clean:
	rm -f $(DEPS) $(ALL_LOCAL_OBJECTS) $(BINARIES)

-include $(DEPS)
