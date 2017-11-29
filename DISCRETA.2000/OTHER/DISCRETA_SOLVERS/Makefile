OBJPATH = ../../obj/$(DISCRETA_ARCH)
OPTIONS=$(LLLOPTIONS)
OBJ = $(OBJPATH)/SPREAD
CC = $(ANSICC)

all: $(OBJPATH)/BIN/discreta_spread  \
	$(OBJPATH)/BIN/discreta_lll  \
	$(OBJPATH)/BIN/discreta_dance  \
	
################################################################################
#  spread
################################################################################

$(OBJPATH)/BIN/discreta_spread: spread.c 
	$(CC) $(OPTIONS) -o $(OBJPATH)/BIN/discreta_spread spread.c

################################################################################
#  LLL
################################################################################

$(OBJPATH)/BIN/discreta_lll: solvedesign.c $(OBJPATH)/BIN/stripcolumns $(OBJPATH)/BIN/solvediophant
	$(CC) $(OPTIONS) -o $(OBJPATH)/BIN/discreta_lll solvedesign.c 

$(OBJPATH)/BIN/stripcolumns: stripcolumns.c 
	$(CC) $(OPTIONS) stripcolumns.c -o $(OBJPATH)/BIN/stripcolumns

$(OBJPATH)/BIN/solvediophant: solvediophant.c $(OBJPATH)/OTHER/diophant.o diophant.h
	$(CC) $(OPTIONS) -o $(OBJPATH)/BIN/solvediophant solvediophant.c $(OBJPATH)/OTHER/diophant.o -lm 

$(OBJPATH)/OTHER/diophant.o: diophant.c
	$(CC) $(OPTIONS) -c diophant.c -o $(OBJPATH)/OTHER/diophant.o

###############################################################################
#  dance
################################################################################

$(OBJPATH)/BIN/discreta_dance: dance.c $(OBJPATH)/BIN/stripcolumnsdance
	$(CC) $(OPTIONS) -o $(OBJPATH)/BIN/discreta_dance dance.c

$(OBJPATH)/BIN/stripcolumnsdance: stripcolumnsdance.c
	$(CC)  $(LOCALOPT) -o $(OBJPATH)/BIN/stripcolumnsdance stripcolumnsdance.c 


clean: ; 
	- rm -f *.bck core solution* lll.lo* KMtmp km_fname  *.log spreadtmp.txt
	- rm $(OBJPATH)/OTHER/diophant.o $(OBJPATH)/BIN/discreta_spread $(OBJPATH)/BIN/discreta_lll $(OBJPATH)/BIN/discreta_dance

