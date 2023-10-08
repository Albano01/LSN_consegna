CC = mpicxx
CFLAGS = -Wall -O3 --std=c++11

main.exe : main.o random.o file_classi.h
	$(CC) random.o main.o -o main.exe
main.o : main.cpp file_classi.h
	$(CC) -c main.cpp -o main.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)

run:
	mpiexec --hostfile hostfile -np 4 main.exe

clean :
	rm *.o main.exe *.out