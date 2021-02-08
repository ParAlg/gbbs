all: intervalTree_pam intervalTree_cpam

intervalTree_pam:		intervalTree.cpp
	g++ -DNDEBUG -DUSE_PAM -O3 -mcx16 -march=native -DHOMEGROWN -pthread -std=c++17 -Wall -I../parlaylib/include -I../ -o intervalTree_pam intervalTree.cpp -L/usr/local/lib -ljemalloc
#	g++ -DUSE_PAM -O0 -g -DDEBUG -mcx16 -march=native -DHOMEGROWN -pthread -std=c++17 -Wall -I../parlaylib/include -I../ -o intervalTree_pam intervalTree.cpp -L/usr/local/lib -ljemalloc

intervalTree_cpam:		intervalTree.cpp
	g++ -DNDEBUG -O3 -mcx16 -march=native -DHOMEGROWN -pthread -std=c++17 -Wall -I../parlaylib/include -I../ -o intervalTree_cpam intervalTree.cpp -L/usr/local/lib -ljemalloc
#	g++ -O0 -g -DDEBUG -mcx16 -march=native -DHOMEGROWN -pthread -std=c++17 -Wall -I../parlaylib/include -I../ -o intervalTree_cpam intervalTree.cpp -L/usr/local/lib


#	g++ -O0 -g -DDEBUG -mcx16 -march=native -DHOMEGROWN -pthread -std=c++17 -Wall -I../parlaylib/include -I../ -o intervalTree intervalTree.cpp -L/usr/local/lib -ljemalloc
#	g++ -O3 -mcx16 -march=native -DHOMEGROWN -pthread -std=c++17 -Wall -I../parlaylib/include -o intervalTree intervalTree.cpp -L/usr/local/lib -ljemalloc

clean:
	rm -f test intervalTree_cpam intervalTree_pam
