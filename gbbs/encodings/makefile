include ../../makefile.variables

INCLUDE_DIRS = -I../../ -I../../external/parlaylib/include/

OBJDIR = ../../bin/gbbs/encodings/

ALL_PRE = byte byte_pd byte_pd_amortized
ALL= $(addprefix $(OBJDIR), $(addsuffix .o, $(ALL_PRE))) $(addprefix $(OBJDIR), $(addsuffix .a, $(ALL_PRE)))

all: $(ALL)

$(OBJDIR)byte.o : byte.cc
	$(CC) $(INCLUDE_DIRS) $(CFLAGS) $(PFLAGS) -c $< -o $@

$(OBJDIR)byte_pd.o : byte_pd.cc
	$(CC) $(INCLUDE_DIRS) $(CFLAGS) $(PFLAGS) -c $< -o $@

$(OBJDIR)byte_pd_amortized.o : byte_pd_amortized.cc
	$(CC) $(INCLUDE_DIRS) $(CFLAGS) $(PFLAGS) -c $< -o $@

$(OBJDIR)%.a : $(OBJDIR)%.o
	ar -rcs $@ $<

.PHONY : clean

clean :
	rm -f *.o $(ALL)
