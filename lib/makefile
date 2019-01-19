CFLAGS = -mcx16 -O2 -ldl -std=c++11 -march=native

OMPFLAGS = -DOPENMP -fopenmp
CILKFLAGS = -DCILK -fcilkplus
HGFLAGS = -DHOMEGROWN -pthread

ifdef CLANG
CC = clang++
PFLAGS = $(CILKFLAGS)
else ifdef CILK
CC = g++
PFLAGS = $(CILKFLAGS)
else ifdef OPENMP
CC = g++
PFLAGS = $(OMPFLAGS)
else ifdef HOMEGROWN
CC = g++
PFLAGS = $(HGFLAGS)
else ifdef SERIAL
CC = g++
PFLAGS =
else # default is cilk
CC = g++
PFLAGS = $(CILKFLAGS)
endif

AllFiles = allocator.h alloc.h bag.h binary_search.h block_allocator.h collect_reduce.h concurrent_stack.h counting_sort.h get_time.h hash_table.h histogram.h integer_sort.h list_allocator.h memory_size.h merge.h merge_sort.h monoid.h parallel.h parse_command_line.h par_string.h quicksort.h random.h random_shuffle.h reducer.h sample_sort.h seq.h sequence_ops.h sparse_mat_vec_mult.h time_operations.h transpose.h utilities.h scheduler.h

time_tests:	$(AllFiles) time_tests.cpp time_operations.h
	$(CC) $(CFLAGS) $(PFLAGS) time_tests.cpp -o time_tests

test_scheduler:	test_scheduler.cpp scheduler.h
	$(CC) $(CFLAGS) $(PFLAGS)

test_scheduler_%:	test_scheduler.cpp scheduler.h
	$(CC) $($(subst test_scheduler_,,$@)FLAGS) $(CFLAGS) test_scheduler.cpp -o $@

test_schedulers: test_scheduler_OMP test_scheduler_CILK test_scheduler_HG

all:	time_tests

clean:
	rm -f time_tests test
