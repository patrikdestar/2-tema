OMPDO:
	clear
	make -B OMP

OMP:
	g++ -fopenmp -g -v main.cpp -Wall -o something
	./something