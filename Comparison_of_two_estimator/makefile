Simulation: Simulation.o Newton_method.o sketch.o MurmurHash3.o
	g++ Simulation.o Newton_method.o sketch.o MurmurHash3.o -o Simulation -lmpfr -lgmp

Simulation.o: Simulation.cpp 
	g++ -c Simulation.cpp

Newton_method: Newton_method.o
	g++ Newton_method.o -o Newton_method -lmpfr -lgmp

Newton_method.o: Newton_method.cpp
	g++ -c Newton_method.cpp 

sketch: sketch.o MurmurHash3.o
	g++ sketch.o MurmurHash3.o -o sketch

sketch.o: sketch.cpp 
	g++ -c sketch.cpp

MurmurHash3: MurmurHash3.o
	g++ MurmurHash3.o -o MurmurHash3

MurmurHash3.o: MurmurHash3.cpp
	g++ -c MurmurHash3.cpp

clean:
	rm -f *.o

