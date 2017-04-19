CFLAGS = -Wno-deprecated -Wall -O -Wno-write-strings
SRC = dnumrecipes.cc dcomplex.cc
OBJ = dnumrecipes.o dcomplex.o

all: libdnumrec.a

libdnumrec.a:$(OBJ)
	ar rs libdnumrec.a $(OBJ)

dnumrecipes.o:dnumrecipes.cc
	g++ $(CFLAGS) -c dnumrecipes.cc

dcomplex.o:dcomplex.cc
	g++ $(CFLAGS) -c dcomplex.cc

clean:
	rm -f $(OBJ) libdnumrec.a 
