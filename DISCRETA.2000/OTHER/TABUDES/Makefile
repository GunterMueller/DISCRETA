OBJPATH = ../../obj/$(DISCRETA_ARCH)
OBJ = $(OBJPATH)/TABUDES

#CC = g++
CLIBS= -lm -static

#CFLAGS= -g -Wall
CFLAGS= -O3 -Wall
TABUDESOBJECTS= $(OBJ)/dmain.o $(OBJ)/discreta.o $(OBJ)/utils.o $(OBJ)/move.o $(OBJ)/dp.o \
$(OBJ)/matrix.o $(OBJ)/attr.o $(OBJ)/tenure.o $(OBJ)/tabulist.o \
$(OBJ)/neighborhood.o $(OBJ)/solution.o $(OBJ)/args.o $(OBJ)/checked.o \
$(OBJ)/penalty.o $(OBJ)/algorithms.o $(OBJ)/warning.o 

all: $(OBJPATH)/BIN/orbithint $(OBJPATH)/BIN/tabudes 

$(OBJPATH)/BIN/orbithint: orbithint.c
	$(CC) $(CFLAGS) -o $(OBJPATH)/BIN/orbithint orbithint.c

$(OBJPATH)/BIN/tabudes:	
	cp tabudes.linux $(OBJPATH)/BIN/tabudes 

#$(OBJPATH)/BIN/tabudes:	$(TABUDESOBJECTS)
#	$(CC) $(CFLAGS) -o $(OBJPATH)/BIN/tabudes $(TABUDESOBJECTS) $(CLIBS)

tabudes.tar:
	tar vcf tabudes.tar Makefile algorithms.C algorithms.H args.C args.H attr.C attr.H checked.C checked.H discreta.C discreta.H dp.H dp.C matrix.H matrix.C move.C move.H neighborhood.C neighborhood.H penalty.C penalty.H solution.C solution.H tabulist.H tabulist.C tenure.C tenure.H utils.C utils.H warning.C warning.H packcover.H limits.H main.H dmain.C dpackcover changepenalty tabudes.tex

clean:
	rm -f $(OBJ)/*.o *~* core

$(OBJ)/dmain.o:	dmain.C
	$(CC) $(CFLAGS) -O0 -c dmain.C -o $(OBJ)/dmain.o

$(OBJ)/%.o:	%.C
	$(CC) $(CFLAGS) -o $@ -c $<
