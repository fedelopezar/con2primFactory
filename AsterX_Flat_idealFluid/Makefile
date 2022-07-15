

CXX=g++
INC=-I.. -I. -I/usr/local/hdf5/include
LDIR=-L. -L/usr/local/hdf5/lib
LIBS=-lhdf5 -lhdf5_cpp #-lpthread -lm
CXXFLAGS=-std=c++11 $(INC) -L$(LDIR) $(LIBS)


DEPS = AsterX_Flat_idealFluid.hxx  ../con2primFactory.hxx 

OBJ = ../2DNRNoble.o  AsterX_Flat_idealFluid.o main.o

%.o: %.cxx $(DEPS)
	$(CXX) -c -o $@ $< $(CXXFLAGS)

all: $(OBJ)
	$(CXX) -o AsterX_Flat_idealFluid $^ $(CXXFLAGS)

clean:
	rm *.o  ../*o AsterX_Flat_idealFluid
