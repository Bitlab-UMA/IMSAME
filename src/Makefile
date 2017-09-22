CC=gcc
CXX=g++
CFLAGS=-O3 -D_FILE_OFFSET_BITS=64 -D_LARGEFILE64_SOURCE -Wall #-DVERBOSE
BIN=../bin

all: IMSAME revcomp combine covident



IMSAME: IMSAME.c
	$(CC) $(CFLAGS) alignmentFunctions.c -lpthread -lm commonFunctions.c -msse -lm IMSAME.c -lpthread -lm -o $(BIN)/IMSAME
	
revcomp: reverseComplement.c
	$(CC) $(CFLAGS) commonFunctions.c -lm  reverseComplement.c -o $(BIN)/revComp
	
combine: combine_reads.c
	$(CC) $(CFLAGS) commonFunctions.c -lm combine_reads.c -o $(BIN)/combine

covident: covident.c
	$(CC) $(CFLAGS) covident.c -o $(BIN)/combine

clean:
	rm -rf $(BIN)/IMSAME $(BIN)/revComp $(BIN)/combine
