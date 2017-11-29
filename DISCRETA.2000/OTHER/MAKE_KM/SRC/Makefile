CC = cc
LIBS= -lm

FLAGS= -O2

OBJECTS= pack.o bincoef.o setoper.o utils.o

all:	pres kmgen pack.bin

.c.o :
	$(CC) $(FLAGS) -c $<

pack.bin:	$(OBJECTS)
	$(CC) $(FLAGS) -o pack.bin $(OBJECTS) $(LIBS)

pres:	pres.o
	$(CC) $(FLAGS) -o pres pres.o $(LIBS)

kmgen:	kmgen.o
	$(CC) $(FLAGS) -o kmgen kmgen.o $(LIBS)

clean:
	rm -f *.o *~* core

