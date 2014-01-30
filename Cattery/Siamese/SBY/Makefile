
CC = gcc
CFLAGS = -O3

default: lofar_beams_lib.so

%.o: %.c
	$(CC) $(CFLAGS) -c -fpic $<

lofar_beams_lib.so: lba_beams.o hba_beams.o
	ld  -shared -lm $^ -o $@
