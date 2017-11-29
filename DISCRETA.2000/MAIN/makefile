# DISCRETA makefile
# Anton Betten 1997
#
# MAIN

OBJPATH = $(DISCRETA_HOME)/obj/$(DISCRETA_ARCH)

OBJ = $(OBJPATH)/XDKM
EXE = $(OBJPATH)/BIN


INCLUDES = \
	-I/usr/include/X11 \
	-I/usr/include/Xm \
	-I/usr/include/X11/Xm \
	-I/usr/X11R6/include \
	-I.. \
	-I.


LIBS = $(OBJPATH)/libDISCRETA.a  -lm


MODULES = \
	$(OBJ)/discreta.o \
	$(OBJ)/interface.o \
	$(OBJ)/discreta_stubs.o 

CCFLAGS = $(DISCRETA_CFLAGS) $(DISCRETA_XFLAGS) $(INCLUDES) -DDISCRETA_ARCH=\"$(DISCRETA_ARCH)\"

all : \
	$(EXE)/discreta \
	$(EXE)/discreta_report \
	$(EXE)/discreta_batch \
	$(EXE)/t116.out \
	$(EXE)/t117.out \
	$(EXE)/t122.out \
	$(EXE)/t140e.out \
	$(EXE)/t143.out \
	$(EXE)/t146.out \
	$(EXE)/t147.out \
	$(EXE)/t148.out \
	$(EXE)/t149.out \
	$(EXE)/t150.out \
	$(EXE)/t151.out \
	$(EXE)/t153.out \
	;

$(EXE)/discreta : $(MODULES) $(OBJPATH)/libDISCRETA.a ; 
	$(DISCRETA_CC) $(MODULES)  $(LIBS) $(DISCRETA_XLIBS) -o $(EXE)/discreta

$(EXE)/discreta_report : $(OBJPATH)/libDISCRETA.a $(OBJ)/discreta_report.o ; 
	$(DISCRETA_CC) $(OBJ)/discreta_report.o $(LIBS) -o $(EXE)/discreta_report

$(EXE)/discreta_batch : $(OBJPATH)/libDISCRETA.a $(OBJ)/discreta_batch.o ; 
	$(DISCRETA_CC) $(OBJ)/discreta_batch.o $(LIBS) -o $(EXE)/discreta_batch

$(EXE)/t116.out : t116.C  ; 
	$(DISCRETA_CC) $(CCFLAGS) t116.C $(LIBS) -o $(EXE)/t116.out

$(EXE)/t117.out : t117.C  ; 
	$(DISCRETA_CC) $(CCFLAGS) t117.C $(LIBS) -o $(EXE)/t117.out

$(EXE)/t122.out : t122.C  ; 
	$(DISCRETA_CC) $(CCFLAGS) t122.C $(LIBS) -o $(EXE)/t122.out

$(EXE)/t140e.out : t140e.C  ; 
	$(DISCRETA_CC) $(CCFLAGS) t140e.C $(LIBS) -o $(EXE)/t140e.out

$(EXE)/t143.out : t143.C  ; 
	$(DISCRETA_CC) $(CCFLAGS) t143.C $(LIBS) -o $(EXE)/t143.out

$(EXE)/t146.out : t146.C  ; 
	$(DISCRETA_CC) $(CCFLAGS) t146.C $(LIBS) -o $(EXE)/t146.out

$(EXE)/t147.out : t147.C  ; 
	$(DISCRETA_CC) $(CCFLAGS) t147.C $(LIBS) -o $(EXE)/t147.out

$(EXE)/t148.out : t148.C  ; 
	$(DISCRETA_CC) $(CCFLAGS) t148.C $(LIBS) -o $(EXE)/t148.out

$(EXE)/t149.out : t149.C  ; 
	$(DISCRETA_CC) $(CCFLAGS) t149.C $(LIBS) -o $(EXE)/t149.out

$(EXE)/t150.out : t150.C  ; 
	$(DISCRETA_CC) $(CCFLAGS) t150.C $(LIBS) -o $(EXE)/t150.out

$(EXE)/t151.out : t151.C  ; 
	$(DISCRETA_CC) $(CCFLAGS) t151.C $(LIBS) -o $(EXE)/t151.out

$(EXE)/t153.out : t153.C  ; 
	$(DISCRETA_CC) $(CCFLAGS) t153.C $(LIBS) -o $(EXE)/t153.out

$(OBJ)/discreta.o : discreta.C disc.h ; $(DISCRETA_CC) -c $(CCFLAGS) discreta.C -o $(OBJ)/discreta.o
$(OBJ)/interface.o : interface.C disc.h ; $(DISCRETA_CC) -c $(CCFLAGS) interface.C -o $(OBJ)/interface.o
$(OBJ)/discreta_batch.o : discreta_batch.C ; $(DISCRETA_CC) -c $(CCFLAGS) discreta_batch.C -o $(OBJ)/discreta_batch.o
$(OBJ)/discreta_stubs.o : discreta_stubs.C disc.h dial_g.C dial_well.C dial_solvable.C dial_linear.C dial_sporadic.C dial_solid.C cb.C db.C tvt.C ; $(DISCRETA_CC) -c $(CCFLAGS) discreta_stubs.C -o $(OBJ)/discreta_stubs.o

$(OBJ)/discreta_report.o : discreta_report.C disc.h ; $(DISCRETA_CC) -c $(CCFLAGS) discreta_report.C -o $(OBJ)/discreta_report.o


gett: ;
	cp ../TEST/t116.C .
	cp ../TEST/t117.C .
	cp ../TEST/t122.C .
	cp ../TEST/t140e.C .
	cp ../TEST/t143.C .
	cp ../TEST/t146.C .
	cp ../TEST/t147.C .
	cp ../TEST/t148.C .
	cp ../TEST/t149.C .
	cp ../TEST/t150.C .
	cp ../TEST/t151.C .
	cp ../TEST/t153.C .
	cp ../TEST/t_graph.C .

clean : 
	- rm -f *.o *.out *.idx *.db *.log KM_*.bin *.BAK lll.output *.tbl *.tdf 
	- rm -f *.in.bin *.out.bin d?.html tables.dat despar_searchresult.txt 
	- rm -f labra_*.dsc 
	- rm -f KM_*.dsc KM_*.txt.lp
	- rm -f solutions dp.txt
	- rm -f $(OBJ)/*.o $(OBJPATH)/BIN/discreta $(OBJ)/*.idx \
			$(OBJ)/*.db $(OBJ)/*.log $(OBJ)/KM_*.txt $(OBJ)/KM_*.bin \
			$(OBJ)/*.BAK $(OBJ)/lll.output 

Clean :
	rm -f lowindex.dat permutation.dat permutation
	rm -f *.dsc KM_*.txt
	rm -f $(OBJ)/*.dsc


