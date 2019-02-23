CC     = gcc
CFLAGS = -g -Wall -Wextra -O3 -DNDEBUG -march=native

HEADERS = api.h bch.h bin-lwe.h ecc.h lac_param.h rand.h fips202.h bch128.h bch192.h bch256.h
SOURCES = ake.c bch.c bin-lwe.c ecc.c encrypt.c ke.c kem.c rand.c fips202.c test_correctness.c test_cpucycles.c test_speed.c main.c
OBJECTS = ake.o ke.o kem.o encrypt.o ecc.o bch.o bin-lwe.o rand.o fips202.o test_correctness.o test_cpucycles.o test_speed.o main.o

all: lac $(OBJECTS)

lac : $(OBJECTS)
	$(CC) -o lac $(OBJECTS)

%.o: %.c $(HEADERS)
	$(CC) $(CFLAGS) -c -o $@ $<

.PHONY: clean
clean:
	-$(RM) lac
	-$(RM) *.o
