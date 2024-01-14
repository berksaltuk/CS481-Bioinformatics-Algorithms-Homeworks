all: hw3

hw2: main.cpp
	g++ -std=c++11 -Wall -g -o hw3 main.cpp

run:
	./hw3 -i input.fasta -o output.phy -s 1:-1:-5 
	./hw3 -i input2.fasta -o output2.phy -s 1:-1:-5 

clean:
	rm -fr hw3 *~
