.phony: all
all: calc

calc: calc.cpp
	g++ -std=c++20 $^ -o $@

.phony: run
run:
	./calc
