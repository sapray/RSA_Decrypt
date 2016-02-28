CFLAGS=-g -O2 -Wall 
CC = mpicc

PROGRAM_NAME= factor
OBJS = main.o

$(PROGRAM_NAME): $(OBJS)
	$(CC) -o $@ $? -lgmp -lm

clean:
	rm  *.o $(PROGRAM_NAME) *~
