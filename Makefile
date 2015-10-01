CC=g++
CFLAGS=-I
DEPS=vector.o

main: main.cpp $(DEPS)
	$(CC) main.cpp $(DEPS) -o main $(FLAGS)

%.o: %.cpp
	$(CC) $< -o $@ $(FLAGS) -c

clean:
	rm -rf *.o main
