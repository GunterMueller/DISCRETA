OBJPATH = ../../obj/$(DISCRETA_ARCH)
OPTIONS=$(LLLOPTIONS)
# OBJ = $(OBJPATH)/MAXCLIQUE

all: \
   $(OBJPATH)/BIN/wclique \
	$(OBJPATH)/BIN/discreta2wclique \
	$(OBJPATH)/BIN/wclique2discreta \
	$(OBJPATH)/BIN/discreta_clique \
	$(OBJPATH)/BIN/discreta_cliqued \
	$(OBJPATH)/BIN/group2graph 

$(OBJPATH)/BIN/wclique: wclique.c
	$(ANSICC)  $(OPTIONS) -o $(OBJPATH)/BIN/wclique wclique.c 
		
$(OBJPATH)/BIN/discreta2wclique: discreta2wclique.c
	$(ANSICC)  $(OPTIONS) -o $(OBJPATH)/BIN/discreta2wclique discreta2wclique.c 
		
$(OBJPATH)/BIN/wclique2discreta: wclique2discreta.c
	$(ANSICC)  $(OPTIONS) -o $(OBJPATH)/BIN/wclique2discreta wclique2discreta.c 

$(OBJPATH)/BIN/discreta_clique: discreta_clique
	cp discreta_clique $(OBJPATH)/BIN/discreta_clique
	chmod a+x $(OBJPATH)/BIN/discreta_clique

$(OBJPATH)/BIN/discreta_cliqued: discreta_cliqued
	cp discreta_cliqued $(OBJPATH)/BIN/discreta_cliqued
	chmod a+x $(OBJPATH)/BIN/discreta_cliqued

$(OBJPATH)/BIN/group2graph: group2graph.c
	$(ANSICC)  $(OPTIONS) -o $(OBJPATH)/BIN/group2graph group2graph.c 

clean: ; 
