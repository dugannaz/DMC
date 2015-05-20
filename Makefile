CXX = g++ -fopenmp -O2

CXXFLAGS =

OBJS1 =  dmc.o \
		HeliumPsi.o \
		HarmonicPsi.o \
		Walker.o \
		CoulombPotential.o \
		HarmonicPotential.o


LIBS = -lgsl -lgslcblas -lnewmat -lm
INC = /home/nazim/program/gsl-1.16/

TARGET =       dmc

all:    $(TARGET)

dmc:	$(OBJS1) $(INC)
	$(CXX) -o $@ $(OBJS1) $(LIBS)         

clean:
	rm -f $(OBJS1) $(OBJS2) $(TARGET)
