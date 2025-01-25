CC = g++
CFLAGS = 
LDFLAGS = -lstdc++

all: p1

p1: MSTgraph.o
	$(CC) $(CFLAGS) $(LDFLAGS) MSTgraph.o -o p1

MSTgraph.o: MSTgraph.cpp
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -f MSTgraph.o p1
