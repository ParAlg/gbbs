index:		index.cpp
	g++ -O0 -g -DDEBUG -mcx16 -march=native -DHOMEGROWN -pthread -std=c++17 -Wall -I../parlaylib/include -I../ -o index index.cpp -L/usr/local/lib -ljemalloc
#	g++ -O3 -mcx16 -march=native -DHOMEGROWN -pthread -std=c++17 -Wall -I../parlaylib/include -o index index.cpp -L/usr/local/lib -ljemalloc

clean:
	rm -f test index
