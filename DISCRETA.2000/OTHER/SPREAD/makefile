# $Header: /usr/local/cvsroot/spread/makefile,v 1.2 1998/01/15 10:20:06 alfred Exp $
#####################################################
OBJPATH = ../../obj/$(DISCRETA_ARCH)
OPTIONS=$(LLLOPTIONS)
#OBJ = $(OBJPATH)/SPREAD

all: $(OBJPATH)/BIN/discreta_spread \
     $(OBJPATH)/BIN/discreta_spreadpacking
	
$(OBJPATH)/BIN/discreta_spread: spread.c 
	$(ANSICC) $(OPTIONS) -o $(OBJPATH)/BIN/discreta_spread spread.c

$(OBJPATH)/BIN/discreta_spreadpacking: spreadpacking.c 
	$(ANSICC) $(OPTIONS) -o $(OBJPATH)/BIN/discreta_spreadpacking spreadpacking.c

clean: ; 
