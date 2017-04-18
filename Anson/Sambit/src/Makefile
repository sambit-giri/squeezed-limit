CFLAGS = -Wall -I../libdnumrec/ -I../libdcosmo/
LIBFLAGS = -L../libdnumrec/ -L../libdcosmo/
OBJ = main.o
EXEC = ESMR.x
#OBJ = reionization.o
#EXEC = reionzation.x

all: $(EXEC)

$(EXEC): $(OBJ)
	g++ $(CFLAGS) -o $(EXEC) $(OBJ) $(LIBFLAGS) -lm -ldcos -ldnumrec

main.o: main.cpp
	g++ $(CFLAGS) -c main.cpp
	#g++ $(CFLAGS) -c reionization_history.cpp

clean:
	rm -f $(OBJ) $(EXEC)