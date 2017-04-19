CFLAGS = -Wall -I../libdnumrec/ -I../libdcosmo/
LIBFLAGS = -L../libdnumrec/ -L../libdcosmo/
OBJ = main.o
EXEC = ESMR.x

all: $(EXEC)

$(EXEC): $(OBJ)
	g++ $(CFLAGS) -o $(EXEC) $(OBJ) $(LIBFLAGS) -lm -ldcos -ldnumrec

main.o: main.cpp
	g++ $(CFLAGS) -c main.cpp

clean:
	rm -f $(OBJ) $(EXEC)