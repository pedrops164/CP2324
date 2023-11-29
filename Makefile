CC = gcc
SRC = src/
CFLAGS = -pg -Ofast -ftree-vectorize -fprefetch-loop-arrays -mavx -msse4.1 -funroll-all-loops -march=ivybridge -flto
#CFLAGS += -flto -fgraphite-identity -floop-nest-optimize -mtune=native
.DEFAULT_GOAL = all

all: MDseq.exe MDpar.exe

#MD.exe: $(SRC)/MD.cpp
#	$(CC) $(CFLAGS) $(SRC)MD.cpp -lm -o MD.exe

MDseq.exe: $(SRC)/MDseq.cpp
	module load gcc/11.2.0; $(CC) $(CFLAGS) $(SRC)MDseq.cpp -lm -o MDseq.exe

MDpar.exe: $(SRC)/MDpar.cpp
	module load gcc/11.2.0; $(CC) $(CFLAGS) $(SRC)MDpar.cpp -lm -fopenmp -o MDpar.exe

clean:
	rm ./MDpar.exe
	rm ./MDseq.exe

runseq:
	./MDseq.exe < inputdata.txt

runpar:
	OMP_NUM_THREADS=40 ./MDpar.exe < inputdata.txt
