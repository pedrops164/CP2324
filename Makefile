CC = g++
SRC = src/
#CFLAGS = -pg -O2 -ftree-vectorize -fprefetch-loop-arrays -march=native -fno-omit-frame-pointer # -funroll-all-loops# -ftree-vectorizer-verbose=6
CFLAGS = -pg -Ofast -ftree-vectorize -fprefetch-loop-arrays -march=native -fno-math-errno -funroll-all-loops# -ftree-vectorizer-verbose=6

.DEFAULT_GOAL = MD.exe

MD.exe: $(SRC)/MD.cpp
	$(CC) $(CFLAGS) $(SRC)MD.cpp -lm -o MD.exe

clean:
	rm ./MD.exe

run:
	./MD.exe < inputdata.txt
