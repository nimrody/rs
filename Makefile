CFLAGS=-g -Wall
CC=gcc

OBJ=gf.o rs.o print.o

rstest: $(OBJ)
	$(CC)  -o rstest $(OBJ)

clean:
	rm -rf $(OBJ) rstest	

