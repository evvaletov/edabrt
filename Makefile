CC = gcc
CFLAGS = -Wall -Wextra -O2
LDLIBS = -lm

all: edabrt

edabrt: edabrt.c
	$(CC) $(CFLAGS) -o $@ $< $(LDLIBS)

clean:
	rm -f edabrt *.o *.d

.PHONY: all clean
