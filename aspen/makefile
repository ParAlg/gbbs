all: test_cpam

test_cpam:		test.cpp
	g++ -O3 -mcx16 -march=native -DHOMEGROWN -pthread -std=c++17 -Wall -I../parlaylib/include -I../ -o test_cpam test.cpp -L/usr/local/lib -ljemalloc
#	g++ -DUSE_PAM -O0 -g -DDEBUG -mcx16 -march=native -DHOMEGROWN -pthread -std=c++17 -Wall -I../parlaylib/include -I../ -o test_pam test.cpp -L/usr/local/lib -ljemalloc

#test_pam:		test.cpp
#	g++ -DUSE_PAM -O3 -mcx16 -march=native -DHOMEGROWN -pthread -std=c++17 -Wall -I../parlaylib/include -I../ -o test_pam test.cpp -L/usr/local/lib -ljemalloc
#	g++ -O0 -g -DDEBUG -mcx16 -march=native -DHOMEGROWN -pthread -std=c++17 -Wall -I../parlaylib/include -I../ -o test_cpam test.cpp -L/usr/local/lib


#	g++ -O0 -g -DDEBUG -mcx16 -march=native -DHOMEGROWN -pthread -std=c++17 -Wall -I../parlaylib/include -I../ -o test test.cpp -L/usr/local/lib -ljemalloc
#	g++ -O3 -mcx16 -march=native -DHOMEGROWN -pthread -std=c++17 -Wall -I../parlaylib/include -o test test.cpp -L/usr/local/lib -ljemalloc

clean:
	rm -f test test_pam test_cpam
