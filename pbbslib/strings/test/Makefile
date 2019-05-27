CC = g++
INCLUDE = -I../../
CFLAGS = -mcx16 -O3 -std=c++14 -DHOMEGROWN $(INCLUDE)
LFLAGS = -pthread
OBJS = suffix_array_test.o suffix_array.o

%.o : %.cpp
	$(CC) $(CFLAGS) -c $< -o $@

suffix_array_test : $(OBJS)
	$(CC) $(LFLAGS) -o suffix_array_test $(OBJS)

clean :
	rm *.o suffix_array_test
