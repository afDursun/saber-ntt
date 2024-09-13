CC        = /usr/bin/gcc
CFLAGS    = -Wall -Wextra -Wmissing-prototypes -Wredundant-decls \
            -O3 -fomit-frame-pointer -march=native \
            -I/opt/homebrew/opt/openssl@3/include -L/opt/homebrew/opt/openssl@3/lib -lcrypto
NISTFLAGS = -Wno-unused-result -O3 -fomit-frame-pointer -march=native -std=c99 \
            -I/opt/homebrew/opt/openssl@3/include -L/opt/homebrew/opt/openssl@3/lib -lcrypto

RM        = /bin/rm


all: test/kem \


afd: test/PQCgenKAT_kem \
     test/test_kex \
     test/kem \

SOURCES = pack_unpack.c poly.c fips202.c verify.c cbd.c SABER_indcpa.c kem.c
HEADERS = SABER_params.h pack_unpack.h poly.h rng.h fips202.h verify.h cbd.h SABER_indcpa.h 

test/test_kex: $(SOURCES) $(HEADERS) rng.o test/test_kex.c
	$(CC) $(CFLAGS) -o $@ $(SOURCES) rng.o test/test_kex.c -lcrypto

test/PQCgenKAT_kem: $(SOURCES) $(HEADERS) rng.o test/PQCgenKAT_kem.c 
	$(CC) $(NISTFLAGS) -o $@ $(SOURCES) rng.o test/PQCgenKAT_kem.c -lcrypto

test/kem: $(SOURCES) $(HEADERS) rng.o test/kem.c
	$(CC) $(CFLAGS) -o $@ $(SOURCES) rng.o test/kem.c -lcrypto

rng.o: rng.c
	$(CC) $(NISTFLAGS) -c rng.c -lcrypto -o $@ 

# fips202.o: fips202.c
# 	$(CLANG) -c $^ -o $@

.PHONY: clean test

test: 
	./test/test_kex
	./test/PQCgenKAT_kem
	./test/kem

clean:
	-$(RM) -f *.o
	-$(RM) -rf test/test_kex
	-$(RM) -rf test/kem
	-$(RM) -rf test/PQCgenKAT_kem
	-$(RM) -f *.req
	-$(RM) -f *.rsp
