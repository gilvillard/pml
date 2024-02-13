INC_DIR = ../../include
LIB_DIR = ../../lib
OBJS = $(wildcard $(LIB_DIR)/*)
SRC = $(wildcard test*c)
EXEC = $(patsubst %.c, %, $(SRC))
CHK = $(patsubst %.c, %.chk, $(patsubst test_%, %, $(SRC)))
COMPILE = $(patsubst %.c, %.exe, $(patsubst test_%, %, $(SRC)))

CC = gcc
CFLAGS = -Wall -O3 -march=native -I$(INC_DIR) -I$(INC_DIR) -I/usr/local/include -I/opt/homebrew/include
LIBS = -lflint -lmpfr -lgmp -lm
LDFLAGS = -L/opt/homebrew/lib -L/usr/local/lib

.PHONY: clean

clean:
	rm -f $(EXEC) $(CHK)

all:	$(COMPILE)

run:	$(CHK)

%.exe:
	if test -e test_$*.c; then $(CC) -o test_$* test_$*.c $(OBJS) $(CFLAGS) $(LDFLAGS) $(LIBS); fi

%.chk:
	if test -e test_$*; then ./test_$* > $*.chk; fi
