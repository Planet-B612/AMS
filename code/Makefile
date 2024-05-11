all: ams run test Gene_cost

ams: *.cpp *.h
		g++ ams.cpp -Wall -std=c++11 -O3 ../dSFMT/dSFMT.c -o ams

# run: *.cpp *.h
# 		g++ run.cpp -Wall -std=c++11 -O3 ../dSFMT/dSFMT.c -o run

test: *.cpp *.h
		g++ test.cpp -Wall -std=c++11 -O3 ../dSFMT/dSFMT.c -o test

Gene_cost: *.cpp *.h
		g++ Gene_cost.cpp -Wall -std=c++11 -O3 ../dSFMT/dSFMT.c -o Gene_cost