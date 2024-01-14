all: hw2

hw2: main.cpp
	g++ -std=c++11 -Wall -g -o hw2 main.cpp

run:
	./hw2 -g -p pattern.fasta -t text.fasta -s "+1|-1|-2" -o global.txt
	./hw2 -l -p pattern.fasta -t text.fasta -s "+1|-1|-2" -o local.txt
	./hw2 -g -p pattern2.fasta -t text2.fasta -s "+1|-1|-2" -o global2.txt
	./hw2 -l -p pattern2.fasta -t text2.fasta -s "+1|-1|-2" -o local2.txt
	
clean:
	rm -fr hw2 *~
