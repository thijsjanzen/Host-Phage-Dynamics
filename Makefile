all: test

test : main.cpp mersenne.cpp 
       g++ -std=c++11 main.cpp mersenne.cpp -Wall -O3 -o run
