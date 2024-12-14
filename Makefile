CC = g++
CFLAGS  = -O3 -Wall -pedantic -std=c++17
TARGET = lbm

all: main.cpp 
	$(CC) $(CFLAGS) main.cpp -o $(TARGET)