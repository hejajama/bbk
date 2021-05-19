CXXFLAGS = `gsl-config --cflags` -O2 -openmp -pedantic -I ../amplitudelib_v2/
LDFLAGS = `gsl-config --libs` -lm 

include filelist.m

all: rbk

rbk: $(OBJECTS) 
	g++ $(CXXFLAGS) $(LDFLAGS) $(OBJECTS) -o rbk 
.cpp.o:
	 g++ $(CXXFLAGS) $< -c -o $@
.c.o:
	gcc $(CXXFLAGS) $< -c -o $@

clean:
	rm -f $(OBJECTS) $(AMPLITUDELIBO)
	rm -f rbk	
