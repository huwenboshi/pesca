CC = g++

EIGEN = ./eigen-3.3.3
COPTS =	-I$(EIGEN) -lblas -w -O3 -lboost_program_options -std=c++11 

BINS = pesca

all: $(BINS)

pesca: fit.o post.o main.o mcmc.o utils.o numeric.o
	$(CC) main.o numeric.o fit.o post.o mcmc.o utils.o -o pesca $(COPTS)

%.o: %.cpp
	$(CC) -c $< $(COPTS)

main.o: main.cpp
fit.o: fit.h fit.cpp
post.o: post.h post.cpp
utils.o: utils.h utils.cpp
mcmc.o: mcmc.h mcmc.cpp
numeric.o: numeric.h numeric.cpp

clean:
	rm -f *.o $(BINS)
