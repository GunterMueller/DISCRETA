
OBJPATH = ../../obj/$(DISCRETA_ARCH)

OPTIONS=$(LLLOPTIONS)

OBJ = $(OBJPATH)/LLL_ENUM

#####################################################
#  Digital Unix 4.0:
#OPTIONS= -newc -O4 -fast -inline speed
#
#####################################################
#  Standard (z.B Digital Unix 3.2):
#OPTIONS= -O2
#
#####################################################
#  GNU gcc + profile:
#OPTIONS= -O2 -Wall  -pg -funroll-loops -fstrength-reduce 
#
#####################################################
# GNU gcc:
#OPTIONS= -O9 -Wall -funroll-loops -funroll-all-loops -fomit-frame-pointer -finline-functions
#
#####################################################
#  CC=gcc
#CC=cc

all: $(OBJPATH)/BIN/discreta_lll $(OBJPATH)/BIN/discreta_lll_with 
	
$(OBJPATH)/BIN/discreta_lll: $(OBJ)/designl.o llll.h $(OBJ)/lllfpl.o $(OBJ)/enumall.o 
	$(ANSICC) $(OPTIONS) $(OBJ)/designl.o $(OBJ)/lllfpl.o $(OBJ)/enumall.o -lm -o $(OBJPATH)/BIN/discreta_lll 

$(OBJPATH)/BIN/discreta_lll_with: $(OBJ)/designl.o llll.h $(OBJ)/lllfpl.o $(OBJ)/enumall_with.o 
	$(ANSICC) $(OPTIONS) $(OBJ)/designl.o $(OBJ)/lllfpl.o $(OBJ)/enumall_with.o -lm -o $(OBJPATH)/BIN/discreta_lll_with 


$(OBJ)/lllfpl.o: lllfpl.c  llll.h
	$(ANSICC) $(OPTIONS) -DDEEPINSERT -c lllfpl.c -o $(OBJ)/lllfpl.o
#       $(ANSICC) $(OPTIONS) -c lllfpl.c

$(OBJ)/designl.o: designl.c llll.h
	$(ANSICC) $(OPTIONS) -c designl.c -o $(OBJ)/designl.o

$(OBJ)/enumall.o: enumall.c llll.h
	$(ANSICC) $(OPTIONS) -c enumall.c -o $(OBJ)/enumall.o
	
$(OBJ)/enumall_with.o: enumall.c llll.h
	cp enumall.c enumall_with.c
	$(ANSICC) $(OPTIONS) -DPRINTSOLUTION -c enumall_with.c -o $(OBJ)/enumall_with.o
	
#enumall_with_prob.o: enumall.c llll.h
#	cp enumall.c enumall_with_prob.c
#	$(ANSICC) $(OPTIONS) -DPRINTSOLUTION -DPROBABILISTIC -c enumall_with_prob.c
	
tar: designl.c lllfpl.c *.h KM.dat makefile enumall.c start
	tar -cvf lllcomp.tar designl.c lllfpl.c *.h KM.dat makefile enumall.c start
	compress lllcomp.tar

clean: ; - rm -f $(OBJ)/*.o lgs *.output *.out core lll.log

