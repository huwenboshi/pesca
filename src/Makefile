CC = g++

# change this to the boost library on your computer / cluster
BOOST_INC = /n/app/boost/1.57.0/include
BOOST_LIB = /n/app/boost/1.57.0/lib

EIGEN = $(PWD)/eigen-3.3.3
BLAS = $(PWD)/BLAS-3.8.0
COPTS =	-I$(EIGEN) -L$(BLAS) -lblas -lgfortran \
		-I$(BOOST_INC) -L$(BOOST_LIB) -lboost_program_options \
		-w -O3 -std=c++11

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
