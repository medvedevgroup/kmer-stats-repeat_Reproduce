P_empty: P_empty.o Newton_method.o DP_Prob_Empty.o
	g++ P_empty.o Newton_method.o DP_Prob_Empty.o -o P_empty -lmpfr -lgmp

P_empty.o: P_empty.cpp 
	g++ -c P_empty.cpp

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

