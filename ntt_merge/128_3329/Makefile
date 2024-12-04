CC =gcc
CFLAGS += -Wall -Wextra -Wpedantic -Wmissing-prototypes -Wredundant-decls \
  -Wshadow -Wpointer-arith -O3 -fomit-frame-pointer
NISTFLAGS += -Wno-unused-result -O3 -fomit-frame-pointer
RM = /bin/rm

LDFLAGS=-lcrypto

SOURCES= ntt.c reduce.c
HEADERS= params.h ntt.h reduce.h

.PHONY: all main clean

all: \
	main

main: $(HEADERS) $(SOURCES) main.c 
	$(CC) $(CFLAGS) -o $@ $(SOURCES) main.c

clean:
	-$(RM) main