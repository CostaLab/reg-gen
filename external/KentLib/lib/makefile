CC=gcc
CXX=g++
MACHTYPE=x86_64
LDFLAGS=-pthread
CFLAGS=-c -fPIC -Wall -I../inc  -DMACHTYPE_$(MACHTYPE)
LIBDIR=./
INCDIR=../inc/
SOURCES=$(shell find $(LIBDIR) -type f -name '*.c')
HEADERS=$(shell find $(INCDIR) -type f -name '*.h')
OBJECTS=$(patsubst %.c,%.o,$(SOURCES))
LIBOUT=jkweb.a

all: $(SOURCES) $(LIBOUT) $(HEADERS)


$(LIBOUT): $(OBJECTS)
	ar rcus $(LIBOUT) $(OBJECTS)

$(OBJECTS): $(HEADERS)

.c.o:
	$(CC) $(CFLAGS) $< -o $@
clean:
	rm -f $(OBJECTS) $(LIBOUT) 



