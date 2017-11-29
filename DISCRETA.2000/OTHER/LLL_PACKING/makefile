OBJPATH = ../../obj/$(DISCRETA_ARCH)
OPTIONS=$(LLLOPTIONS)
OBJ = $(OBJPATH)/LLL_PACKING


all: $(OBJPATH)/BIN/discreta_lllpacking \
     $(OBJPATH)/BIN/discreta_deletelines \
     $(OBJPATH)/BIN/discreta_lllpacking_shell 
	
$(OBJPATH)/BIN/discreta_lllpacking: $(OBJ)/designl.o llll.h $(OBJ)/lllfpl.o $(OBJ)/enumall_with.o 
	$(ANSICC) $(OPTIONS) $(OBJ)/designl.o $(OBJ)/lllfpl.o $(OBJ)/enumall_with.o -lm -o $(OBJPATH)/BIN/discreta_lllpacking

$(OBJPATH)/BIN/discreta_deletelines: deletelines.c 
	$(ANSICC) $(OPTIONS) deletelines.c -o $(OBJPATH)/BIN/discreta_deletelines

$(OBJPATH)/BIN/discreta_lllpacking_shell: discreta_lllpacking_shell.sh
	cp discreta_lllpacking_shell.sh $(OBJPATH)/BIN/discreta_lllpacking_shell
	chmod a+x $(OBJPATH)/BIN/discreta_lllpacking_shell

#$(OBJPATH)/BIN/discreta_lllpacking_shell: discreta_lllpacking_shell.c
#	$(ANSICC) $(OPTIONS) discreta_lllpacking_shell.c -o $(OBJPATH)/BIN/discreta_lllpacking_shell

$(OBJ)/lllfpl.o: lllfpl.c  llll.h
	$(ANSICC) $(OPTIONS) -DGIVENS -DDEEPINSERT -c lllfpl.c -o $(OBJ)/lllfpl.o
#       $(ANSICC) $(OPTIONS) -c lllfpl.c

$(OBJ)/designl.o: designl.c llll.h
	$(ANSICC) $(OPTIONS) -c designl.c -o $(OBJ)/designl.o

$(OBJ)/enumall.o: enumall.c llll.h
	$(ANSICC) $(OPTIONS) -DGIVENS -c enumall.c -o $(OBJ)/enumall.o
	
$(OBJ)/enumall_with.o: enumall.c llll.h
	cp enumall.c enumall_with.c
	$(ANSICC) $(OPTIONS) -DPRINTSOLUTION -DGIVENS -c enumall_with.c -o $(OBJ)/enumall_with.o
	
clean: ; - rm -f $(OBJ)/*.o 

