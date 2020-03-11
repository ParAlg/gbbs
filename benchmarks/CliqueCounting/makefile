include ../../makefile.variables

INCLUDE_DIRS = -I../../

ALL= Clique

all: $(ALL)

LIGRA_OBJS = $(wildcard ../../bin/ligra/*.o)
LIGRA_ENCODING_OBJS = $(wildcard ../../bin/ligra/encodings/*.o)
PBBSLIB_OBJS = $(wildcard ../../bin/pbbslib/*.o)
PBBSLIB_STRINGS_OBJS = $(wildcard ../../bin/pbbslib/strings/*.o)

ligra : 
	make -C ../../ligra/

pbbslib : 
	make -C ../../pbbslib/

ligra_encodings : 
	make -C ../../ligra/encodings/

pbbslib_strings :
	make -C ../../pbbslib/strings/

% : %.cc pbbslib ligra ligra_encodings pbbslib_strings
	$(CC) $(INCLUDE_DIRS) $(CFLAGS) $(PFLAGS) -o $@ $< $(LIGRA_OBJS) $(LIGRA_ENCODING_OBJS) $(PBBSLIB_OBJS) $(PBBSLIB_STRINGS_OBJS)

.PHONY : clean

clean :
	rm -f *.o $(ALL) $(LIGRA_OBJS) $(LIGRA_ENCODING_OBJS) $(PBBSLIB_OBJS) $(PBBSLIB_STRINGS_OBJS)
