all: test

test : main.cpp GetParams.cpp mersenne.cpp 
       g++ -std=c++11 main.cpp GetParams.cpp mersenne.cpp -Wall -O3 -o run
