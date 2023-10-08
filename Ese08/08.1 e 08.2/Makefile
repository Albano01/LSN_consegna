CC = g++
CFLAGS = -Wall -O3 --std=c++11

main.exe : MC.o random.o
	$(CC) random.o MC.o -o main.exe
MC.o : MC.cpp
	$(CC) -c MC.cpp -o MC.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o main.exe