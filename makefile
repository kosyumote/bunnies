.PHONY: run all

debug-flags = -Wall -g
flags = $(debug-flags) -Ofast -std=gnu++11 -fopenmp

run: all
	@echo "Running bunnies simulation:"
	@exec/bunnies

all: exec/bunnies

exec/bunnies: objects/bunnies.o | exec
	@echo "Creating executable..."
	@g++ $(flags) -o $@ $^

exec:
	@mkdir $@

objects/bunnies.o: bunnies.cpp | objects
	@echo "Compiling bunnies.cpp..."
	@g++ $(flags) -o $@ -c $^

objects:
	@mkdir $@
