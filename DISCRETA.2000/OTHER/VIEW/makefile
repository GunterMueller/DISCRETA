# makefile
# Anton Betten 1999
#
# MAIN

OBJPATH = ../../obj/$(DISCRETA_ARCH)

OBJ = $(OBJPATH)/VIEW
EXE = $(OBJPATH)/BIN


INCLUDES = \
	-I/usr/local/include \
	-I/usr/include/X11 \
	-I/usr/include/Xm \
	-I/usr/include/X11/Xm \
	-I/usr/X11R6/include \
	-I$(MOTIF_HOME)/include \
	-I$(MESA_HOME)/include \
	-I.. \
	-I.

LIBS = -L/usr/local/lib -L$MESA_HOME/lib -lGL -lGLU -lMesaGLwM 
#-lglut
#LIBS = $(OBJPATH)/libDISCRETA.a 


MODULES = \
	$(OBJ)/callbacks.o \
	$(OBJ)/file_select.o \
	$(OBJ)/graph.o \
	$(OBJ)/quat.o \
	$(OBJ)/main.o \
	$(OBJ)/render.o \
	$(OBJ)/rotate.o 

CCFLAGS = $(DISCRETA_CFLAGS) $(DISCRETA_XFLAGS) $(INCLUDES) -DDISCRETA_ARCH=\"$(DISCRETA_ARCH)\"

all : $(EXE)/view3d ;

$(EXE)/view3d : $(MODULES) ; 
	$(DISCRETA_CC) $(MODULES) $(LIBS) $(DISCRETA_XLIBS) -o $(EXE)/view3d


$(OBJ)/callbacks.o : callbacks.C ; $(DISCRETA_CC) -c $(CCFLAGS) callbacks.C -o $(OBJ)/callbacks.o
$(OBJ)/file_select.o : file_select.C ; $(DISCRETA_CC) -c $(CCFLAGS) file_select.C -o $(OBJ)/file_select.o
$(OBJ)/graph.o : graph.C ; $(DISCRETA_CC) -c $(CCFLAGS) graph.C -o $(OBJ)/graph.o
$(OBJ)/quat.o : quat.C ; $(DISCRETA_CC) -c $(CCFLAGS) quat.C -o $(OBJ)/quat.o
$(OBJ)/main.o : main.C ; $(DISCRETA_CC) -c $(CCFLAGS) main.C -o $(OBJ)/main.o
$(OBJ)/render.o : render.C ; $(DISCRETA_CC) -c $(CCFLAGS) render.C -o $(OBJ)/render.o
$(OBJ)/rotate.o : rotate.C ; $(DISCRETA_CC) -c $(CCFLAGS) rotate.C -o $(OBJ)/rotate.o



clean : 
	- rm -f *.o *.out 
	- rm -f $(OBJ)/*.o 
	- rm -f $(EXE)/view3d

