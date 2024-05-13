all: ams run test Gene_userInfo

ams: *.cpp *.h
		g++ ams.cpp -Wall -std=c++11 -O3 ../../dSFMT/dSFMT.c -o ams

# run: *.cpp *.h
# 		g++ run.cpp -Wall -std=c++11 -O3 ../dSFMT/dSFMT.c -o run

test: *.cpp *.h
		g++ test.cpp -Wall -std=c++11 -O3 ../../dSFMT/dSFMT.c -o test

Gene_userInfo: *.cpp *.h
		g++ Gene_userInfo.cpp -Wall -std=c++11 -O3 ../../dSFMT/dSFMT.c -o Gene_userInfo