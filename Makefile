CC = g++
SRC = src/
CFLAGS = -pg -Ofast -ftree-vectorize -fprefetch-loop-arrays -march=native -fno-math-errno -mavx -msse4.1 -funroll-all-loops # -fno-fast-math -ftree-vectorizer-verbose=6

.DEFAULT_GOAL = MD.exe

MD.exe: $(SRC)/MD.cpp
	$(CC) $(CFLAGS) $(SRC)MD.cpp -lm -o MD.exe

clean:
	rm ./MD.exe

run:
	./MD.exe < inputdata.txt
