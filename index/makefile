all: index_cpam index_cpam_de index_pam

index_pam:		index.cpp
	g++ -DUSE_PAM -O3 -mcx16 -march=native -DHOMEGROWN -pthread -std=c++17 -Wall -I../parlaylib/include -I../ -o index_pam index.cpp -L/usr/local/lib -ljemalloc
#	g++ -DUSE_PAM -O0 -g -DDEBUG -mcx16 -march=native -DHOMEGROWN -pthread -std=c++17 -Wall -I../parlaylib/include -I../ -o index_pam index.cpp -L/usr/local/lib -ljemalloc

index_cpam:		index.cpp
	g++ -O3 -mcx16 -march=native -DHOMEGROWN -pthread -std=c++17 -Wall -I../parlaylib/include -I../ -o index_cpam index.cpp -L/usr/local/lib -ljemalloc
#	g++ -O0 -g -DDEBUG -mcx16 -march=native -DHOMEGROWN -pthread -std=c++17 -Wall -I../parlaylib/include -I../ -o index_cpam index.cpp -L/usr/local/lib


index_cpam_de:		index.cpp
	g++ -O3 -DUSE_DIFF_ENCODING -mcx16 -march=native -DHOMEGROWN -pthread -std=c++17 -Wall -I../parlaylib/include -I../ -o index_cpam_de index.cpp -L/usr/local/lib -ljemalloc


#	g++ -O0 -g -DDEBUG -mcx16 -march=native -DHOMEGROWN -pthread -std=c++17 -Wall -I../parlaylib/include -I../ -o index index.cpp -L/usr/local/lib -ljemalloc
#	g++ -O3 -mcx16 -march=native -DHOMEGROWN -pthread -std=c++17 -Wall -I../parlaylib/include -o index index.cpp -L/usr/local/lib -ljemalloc

clean:
	rm -f test index_cpam index_pam index_cpam_de
