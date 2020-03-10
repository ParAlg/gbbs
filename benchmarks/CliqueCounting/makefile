include ../../makefile.variables

INCLUDE_DIRS = -I../../

ALL= Clique

all: $(ALL)

LIGRA_OBJS = $(wildcard ../../ligra/*.o)
LIGRA_ENCODING_OBJS = $(wildcard ../../ligra/encodings/*.o)
PBBSLIB_OBJS = $(wildcard ../../pbbslib/*.o)
PBBSLIB_STRINGS_OBJS = $(wildcard ../../pbbslib/strings/*.o)

ligra : 
	make -C ../../ligra/

pbbslib : 
	make -C ../../pbbslib/

ligra_encodings : 
	make -C ../../ligra/encodings/

pbbslib_strings :
	make -C ../../pbbslib/strings/

% : %.cc ligra pbbslib ligra_encodings pbbslib_strings
	$(CC) $(INCLUDE_DIRS) $(CFLAGS) $(PFLAGS) -o $@ $< $(LIGRA_OBJS) $(LIGRA_ENCODING_OBJS) $(PBBSLIB_OBJS) $(PBBSLIB_STRINGS_OBJS)

.PHONY : clean

clean :
	rm -f *.o $(ALL)
