CFLAGS = -Wno-deprecated -Wall -O -I../libdnumrec
LIBFLAGS = -L../libdnumrec
SRC = dcosmology.cc
OBJ = dcosmology.o

all: libdcos.a

libdcos.a:$(OBJ)
	ar rs libdcos.a $(OBJ)

dcosmology.o:dcosmology.cc
	g++ $(CFLAGS) -c dcosmology.cc

clean:
	rm -f $(OBJ) libdcos.a 