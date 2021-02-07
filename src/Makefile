
#NTL_FLAGS=-lntl
#NTL_FLAGS=-I/home/username/ntl/include/ -L/home/username/ntl/lib -lntl
NTL_FLAGS=-I/usr/local/include/  -L/usr/local/lib -lntl
GMP_FLAGS=-lgmp

all: tests_GAHE tests_ScalarNandHE tests_BootstrapperSingleNandHE

utils.o: utils.h utils.cpp
	g++ utils.h utils.cpp -c -g -O2 -std=c++11 -pthread -march=native $(NTL_FLAGS) $(GMP_FLAGS) -lm

GAHE.o: GAHE.h GAHE.cpp utils.o
	g++ GAHE.cpp GAHE.h -c -g -O2 -std=c++11 -pthread -march=native $(NTL_FLAGS) $(GMP_FLAGS) -lm

ScalarNandHE.o: ScalarNandHE.h ScalarNandHE.cpp utils.o
	g++ ScalarNandHE.cpp ScalarNandHE.h -c -g -O2 -std=c++11 -pthread -march=native $(NTL_FLAGS) $(GMP_FLAGS) -lm

BootstrapperSingleNandHE.o: BootstrapperSingleNandHE.h BootstrapperSingleNandHE.cpp GAHE.o ScalarNandHE.o utils.o
	g++ BootstrapperSingleNandHE.cpp BootstrapperSingleNandHE.h -c -g -O2 -std=c++11 -pthread -march=native $(NTL_FLAGS) $(GMP_FLAGS) -lm

tests_BootstrapperSingleNandHE: tests_BootstrapperSingleNandHE.cpp BootstrapperSingleNandHE.o utils.o PolynomialMultiplier.o
	g++ tests_BootstrapperSingleNandHE.cpp *.o -g -O2 -std=c++11 -pthread -march=native $(NTL_FLAGS) $(GMP_FLAGS) -lm -o tests_BootstrapperSingleNandHE

tests_GAHE: tests_GAHE.cpp GAHE.o utils.o
	g++ tests_GAHE.cpp *.o -g -O2 -std=c++11 -pthread -march=native $(NTL_FLAGS) $(GMP_FLAGS) -lm -o tests_GAHE

tests_ScalarNandHE: tests_ScalarNandHE.cpp ScalarNandHE.o utils.o
	g++ tests_ScalarNandHE.cpp *.o -g -O2 -std=c++11 -pthread -march=native $(NTL_FLAGS) $(GMP_FLAGS) -lm -o tests_ScalarNandHE

NTT.o: NTT.h NTT.cpp utils.o
	g++ -c NTT.cpp NTT.h -c -g -O2 -std=c++11 -pthread -march=native $(NTL_FLAGS) $(GMP_FLAGS) -lm

PolynomialMultiplier.o: utils.o NTT.o PolynomialMultiplier.h
	g++ PolynomialMultiplier.h PolynomialMultiplier.cpp -c -g -O2 -std=c++11 -pthread -march=native $(NTL_FLAGS) $(GMP_FLAGS) -lm

clean:
	rm *.o *.gch tests_GAHE tests_ScalarNandHE tests_BootstrapperSingleNandHE 
