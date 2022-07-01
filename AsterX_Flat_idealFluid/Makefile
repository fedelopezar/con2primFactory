

CXX=g++
INCDIR=..
LDIR=.
LIBS=#-lpthread -lm
CXXFLAGS=-std=c++11 -I$(INCDIR) -L$(LDIR) $(LIBS)

DEPS = AsterX_Flat_idealFluid.hxx  ../con2primFactory.hxx 

OBJ = ../2DNRNoble.o  AsterX_Flat_idealFluid.o main.o

%.o: %.cxx $(DEPS)
	$(CXX) -c -o $@ $< $(CXXFLAGS)

AsterX_Flat_idealFluid: $(OBJ)
	$(CXX) -o $@ $^ $(CXXFLAGS)

clean:
	rm *.o 
