INC_DIR = include 
LIB_DIR = lib
OBJS = $(wildcard $(LIB_DIR)/*.o)

CC = gcc
CFLAGS = -Wall -O3 -march=native -I$(INC_DIR) -I/opt/homebrew/include -I/usr/local/include -L$(LIB_DIR) 
LIBS = -lflint -lmpfr -lgmp -lm -lpthread
LDFLAGS = -L/opt/homebrew/lib -L/usr/local/lib

INSTALL_DIR = /usr/local/lib

DIRS = nmod_extra nmod_vec_extra nmod_poly_extra fmpz_extra nmod_mat_extra fmpz_mat_extra nmod_mat_poly_extra nmod_poly_mat_extra sagemath_extra

all:    clean
	$(foreach dir, $(DIRS), if [ -d $(dir) ]; then cd $(dir)/src ; make ; cd ../.. ; fi ;)
	cd lib; make -f makelib; cd ..; 

clean:
	rm -f lib/*.o lib/*.dylib
	rm -f include/*.h
	$(foreach dir, $(DIRS), if [ -d $(dir)/src ]; then cd $(dir)/src ; make clean ; cd ../.. ; fi ;)
	$(foreach dir, $(DIRS), if [ -d $(dir)/test ]; then cd $(dir)/test ; make clean ; cd ../.. ; fi ;)
	$(foreach dir, $(DIRS), if [ -d $(dir)/timings ]; then cd $(dir)/timings ; make clean ; cd ../.. ; fi ;)


.PHONY: doc
doc:
	doxygen doxygen.conf

docclean:
	rm -r include/html
	rm -r include/latex


shared: 
	$(CC) -shared -fPIC -o libpmlflint.dylib $(OBJS) $(CFLAGS)  $(LDFLAGS) $(LIBS)

install:
	cp libpmlflint.dylib $(INSTALL_DIR)/

