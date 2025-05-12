Simulation: Simulation.o Newton_method.o Var.o
	g++ -std=c++17 Simulation.o Newton_method.o Var.o -o Simulation -lmpfr -lgmp

Simulation.o: Simulation.cpp 
	g++ -std=c++17 -c Simulation.cpp

Newton_method: Newton_method.o
	g++ -std=c++17 Newton_method.o -o Newton_method -lmpfr -lgmp

Newton_method.o: Newton_method.cpp
	g++ -std=c++17 -c Newton_method.cpp 

Var: Var.o
	g++ -std=c++17 Var.o -o Var -lmpfr -lgmp

Var.o: Var.cpp
	g++ -std=c++17 -c Var.cpp 

clean:
	rm -f *.o

