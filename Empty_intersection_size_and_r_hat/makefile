Simulation: Simulation.o Newton_method.o DP_Prob_Empty.o
	g++ Simulation.o Newton_method.o DP_Prob_Empty.o -o Simulation -lmpfr -lgmp

Simulation.o: Simulation.cpp 
	g++ -c Simulation.cpp

Newton_method: Newton_method.o
	g++ Newton_method.o -o Newton_method -lmpfr -lgmp

Newton_method.o: Newton_method.cpp
	g++ -c Newton_method.cpp 

DP_Prob_Empty: DP_Prob_Empty.o
	g++ DP_Prob_Empty.o -o DP_Prob_Empty -lmpfr -lgmp

DP_Prob_Empty.o: DP_Prob_Empty.cpp
	g++ -c DP_Prob_Empty.cpp 

clean:
	rm -f *.o

