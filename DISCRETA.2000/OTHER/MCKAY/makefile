OBJPATH = ../../obj/$(DISCRETA_ARCH)
OBJ = $(OBJPATH)/MCKAY
CFLAGS = $(MCKAYOPTIONS)

all :  $(OBJPATH)/BIN/discreta_mckay
	
$(OBJPATH)/BIN/discreta_mckay : $(OBJ)/brute.o brute.h $(OBJ)/possolve.o
	$(ANSICC) -o $(OBJPATH)/BIN/discreta_mckay $(OBJ)/brute.o $(OBJ)/possolve.o  -lm

$(OBJ)/possolve.o : possolve.c brute.h
	$(ANSICC) $(CFLAGS) -c possolve.c -o $(OBJ)/possolve.o

$(OBJ)/brute.o : brute.c brute.h
	$(ANSICC) $(CFLAGS) -c brute.c -o $(OBJ)/brute.o

clean: ; - rm *.o *.out $(OBJ)/*.o $(OBJ)/*.out *.BAK $(OBJPATH)/BIN/discreta_mckay

