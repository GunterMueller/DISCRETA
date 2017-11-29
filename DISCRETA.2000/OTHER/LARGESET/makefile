# $Header: /usr/local/cvsroot/LARGESET/makefile,v 1.1.1.1 1998/01/15 08:51:13 alfred Exp $ $
#####################################################
OBJPATH = ../../obj/$(DISCRETA_ARCH)
OBJ = $(OBJPATH)/LARGESET
OPTIONS=$(LLLOPTIONS)

all:  $(OBJPATH)/BIN/discreta_ls1 $(OBJPATH)/BIN/discreta_ls2
	
$(OBJPATH)/BIN/discreta_ls1: largeset.c 
	$(ANSICC) $(OPTIONS) -o $(OBJPATH)/BIN/discreta_ls1 largeset.c -lm

$(OBJPATH)/BIN/discreta_ls2: largeset_r.c 
	$(ANSICC) $(OPTIONS) -o $(OBJPATH)/BIN/discreta_ls2 largeset_r.c -lm

clean: ; - rm *.o *.out $(OBJ)/*.o $(OBJ)/*.out *.BAK 

