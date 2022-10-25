CXX=g++
INC=-I. -I/usr/local/hdf5/include
LDIR=-L. -L/usr/local/hdf5/lib
LIBS=-lhdf5 -lhdf5_cpp -lm #-lpthread
CXXFLAGS=-std=c++11 $(INC) $(LDIR) $(LIBS)

DEPS = con2primFactory.hxx con2primIdealFluid.hxx

OBJ = 2DNRNoble.o main.o

%.o: %.cxx $(DEPS)
	$(CXX) -c -o $@ $< $(CXXFLAGS)

all: $(OBJ)
	$(CXX) -o idealFluid $^ $(CXXFLAGS)

clean:
	rm *.o idealFluid 2DNRNoble.h5 testResults.pdf
