# makefile
# for Juergen Mueller database

OBJPATH = ../../obj/$(DISCRETA_ARCH)

OBJ = $(OBJPATH)/DB



ALL: \
	$(OBJPATH)/BIN/discreta_dbcreate.out \
	$(OBJPATH)/BIN/discreta_dbinsert.out \
	$(OBJPATH)/BIN/discreta_dbremove.out \
	$(OBJPATH)/BIN/discreta_dbselect.out

#OPTS = -g -c
#CC = cc
#OPTS = -c
#CC = cxx
OPTS =

#LFLAGS = -g
LFLAGS = 

OBJCRE = \
	$(OBJ)/dbcreate.o \
	$(OBJ)/create.o \
	$(OBJ)/baserout.o \
	$(OBJ)/filerout.o

OBJINS = \
	$(OBJ)/dbinsert.o \
	$(OBJ)/init.o \
	$(OBJ)/fio.o \
	$(OBJ)/inserth.o \
	$(OBJ)/inserto.o \
	$(OBJ)/splitd.o \
	$(OBJ)/spliti.o \
	$(OBJ)/splitih.o \
	$(OBJ)/baserout.o \
	$(OBJ)/filerout.o \
	$(OBJ)/tdefini.o

OBJREM = \
	$(OBJ)/dbremove.o \
	$(OBJ)/delete.o \
	$(OBJ)/deleteh.o \
	$(OBJ)/tdefini.o \
	$(OBJ)/filerout.o \
	$(OBJ)/baserout.o \
	$(OBJ)/init.o \
	$(OBJ)/inserth.o \
	$(OBJ)/inserto.o \
	$(OBJ)/splitd.o \
	$(OBJ)/splitih.o \
	$(OBJ)/spliti.o

OBJSEL = \
	$(OBJ)/dbselect.o \
	$(OBJ)/baserout.o \
	$(OBJ)/filerout.o \
	$(OBJ)/init.o \
	$(OBJ)/tdefini.o \
	$(OBJ)/select.o


$(OBJPATH)/BIN/discreta_dbcreate.out : $(OBJCRE) ; 
	$(ANSICC) -o $@ $(OBJCRE) $(LFLAGS) 

$(OBJPATH)/BIN/discreta_dbinsert.out : $(OBJINS) ; 
	$(ANSICC) -o $@ $(OBJINS) $(LFLAGS)

$(OBJPATH)/BIN/discreta_dbremove.out : $(OBJREM);  
	$(ANSICC) -o $@ $(OBJREM) $(LFLAGS) 

$(OBJPATH)/BIN/discreta_dbselect.out : $(OBJSEL);  
	$(ANSICC) -o $@ $(OBJSEL) $(LFLAGS)


$(OBJ)/baserout.o : baserout.c ; $(ANSICC) -c -O $(OPTS) baserout.c -o $(OBJ)/baserout.o
$(OBJ)/create.o   : create.c   ; $(ANSICC) -c -O $(OPTS) create.c   -o $(OBJ)/create.o
$(OBJ)/dbcreate.o : dbcreate.c ; $(ANSICC) -c -O $(OPTS) dbcreate.c -o $(OBJ)/dbcreate.o
$(OBJ)/dbinsert.o : dbinsert.c ; $(ANSICC) -c -O $(OPTS) dbinsert.c -o $(OBJ)/dbinsert.o
$(OBJ)/dbremove.o : dbremove.c ; $(ANSICC) -c -O $(OPTS) dbremove.c -o $(OBJ)/dbremove.o
$(OBJ)/dbselect.o : dbselect.c ; $(ANSICC) -c -O $(OPTS) dbselect.c -o $(OBJ)/dbselect.o
$(OBJ)/delete.o   : delete.c   ; $(ANSICC) -c -O $(OPTS) delete.c   -o $(OBJ)/delete.o
$(OBJ)/deleteh.o  : deleteh.c  ; $(ANSICC) -c -O $(OPTS) deleteh.c  -o $(OBJ)/deleteh.o
$(OBJ)/filerout.o : filerout.c ; $(ANSICC) -c -O $(OPTS) filerout.c -o $(OBJ)/filerout.o
$(OBJ)/fio.o      : fio.c      ; $(ANSICC) -c -O $(OPTS) fio.c      -o $(OBJ)/fio.o
$(OBJ)/init.o     : init.c     ; $(ANSICC) -c -O $(OPTS) init.c     -o $(OBJ)/init.o
$(OBJ)/inserth.o  : inserth.c  ; $(ANSICC) -c -O $(OPTS) inserth.c  -o $(OBJ)/inserth.o
$(OBJ)/inserto.o  : inserto.c  ; $(ANSICC) -c -O $(OPTS) inserto.c  -o $(OBJ)/inserto.o
$(OBJ)/select.o   : select.c   ; $(ANSICC) -c -O $(OPTS) select.c   -o $(OBJ)/select.o
$(OBJ)/sio.o      : sio.c      ; $(ANSICC) -c -O $(OPTS) sio.c      -o $(OBJ)/sio.o
$(OBJ)/splitd.o   : splitd.c   ; $(ANSICC) -c -O $(OPTS) splitd.c   -o $(OBJ)/splitd.o
$(OBJ)/spliti.o   : spliti.c   ; $(ANSICC) -c -O $(OPTS) spliti.c   -o $(OBJ)/spliti.o
$(OBJ)/splitih.o  : splitih.c  ; $(ANSICC) -c -O $(OPTS) splitih.c  -o $(OBJ)/splitih.o
$(OBJ)/tdefini.o  : tdefini.c  ; $(ANSICC) -c -O $(OPTS) tdefini.c  -o $(OBJ)/tdefini.o

clean : ; - rm -f *.o *.out $(OBJ)/*.o $(OBJPATH)/BIN/discreta_db*.out

