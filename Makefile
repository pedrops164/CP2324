CC = g++
SRC = src/
CFLAGS = -pg -Ofast -ftree-vectorize -fprefetch-loop-arrays -march=native -fno-math-errno -mavx -msse4.1 -funroll-all-loops -fno-inline # -fno-fast-math -ftree-vectorizer-verbose=6

.DEFAULT_GOAL = all

all: MDseq.exe MDpar.exe

#MD.exe: $(SRC)/MD.cpp
#	$(CC) $(CFLAGS) $(SRC)MD.cpp -lm -o MD.exe

MDseq.exe: $(SRC)/MDseq.cpp
	module load gcc/11.2.0;
	$(CC) $(CFLAGS) $(SRC)MDseq.cpp -lm -o MDseq.exe

MDpar.exe: $(SRC)/MDpar.cpp
	module load gcc/11.2.0;
	$(CC) $(CFLAGS) $(SRC)MDpar.cpp -lm -fopenmp -o MDpar.exe

clean:
	rm ./MD.exe

runseq:
	./MDseq.exe < inputdata.txt

runpar:
	./MDpar.exe < inputdata.txt
