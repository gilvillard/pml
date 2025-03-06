INC_DIR = ../../include
LIB_DIR = ../../lib
SOURCES = $(wildcard *.c)
HEADERS = $(wildcard ../*.h)
OBJS = $(patsubst %.c, $(LIB_DIR)/%.o, $(SOURCES))
INSTALLED_HEADERS = $(patsubst ../%.h, $(INC_DIR)/%.h, $(HEADERS))

CC = gcc
#CC = clang
#CFLAGS = -fverbose-asm -S -Wall -Wextra -O3 -march=native -I$(INC_DIR)
<<<<<<< HEAD
CFLAGS = -Wall -Wextra -O3 -march=native -I$(INC_DIR) -I/opt/homebrew/include -I/usr/local/include
=======
CFLAGS = -Wall -Wextra -O2 -march=native -I$(INC_DIR)
>>>>>>> 3022ab666e09dcf166c23229f47c6a295ef76d0c

.PHONY: clean

all:    cp $(OBJS)

cp:
	cp -f $(HEADERS) $(INC_DIR)

clean:
	rm -f $(INSTALLED_HEADERS)
	rm -f $(OBJS)

$(LIB_DIR)/%.o:       %.c
	$(CC) -c $< $(CFLAGS) -o $@
