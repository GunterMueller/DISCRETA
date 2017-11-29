OBJPATH = ../../obj/$(DISCRETA_ARCH)

OBJ = $(OBJPATH)/BIN


CC=cc


LIBS= 
#LIBS= -L/usr/X11R6/lib -lXm -lXt -lX11 -lXext -lXpm -lICE -lSM

all: $(OBJPATH)/BIN/grepdocu ;



$(OBJPATH)/BIN/grepdocu : $(OBJ)/grepdocu.o ;
	$(CC)  $(OBJ)/grepdocu.o $(LIBS) -o $(OBJPATH)/BIN/grepdocu

$(OBJ)/grepdocu.o : grepdocu.c ; $(CC) grepdocu.c -c -o $(OBJ)/grepdocu.o

clean: ;
	rm -f $(OBJ)/grepdocu.o $(OBJPATH)/BIN/grepdocu


