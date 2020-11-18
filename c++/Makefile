all: tracer

tracer: tracer.o
	g++ -Wall tracer.o -o tracer -lSDL2 -lpthread

tracer.o: tracer.cpp
	g++ -c -Wall -march=native -O2 -ffast-math tracer.cpp

clean:
	rm -rf *.o tracer
