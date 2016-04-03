.PHONY: run all

debug-flags = -Wall -g
release-flags = -Ofast

other-flags = -std=gnu++11 -fopenmp
flags = $(release-flags) $(other-flags)

run: all
	@echo "Running bunnies simulation:"
	@exec/bunnies

all: exec/bunnies

exec/bunnies: objects/bunnies.o | exec
	@echo "Creating executable..."
	@g++ $(flags) -o $@ $^ -lpng -lz

exec:
	@mkdir $@

objects/bunnies.o: bunnies.cpp | objects
	@echo "Compiling bunnies.cpp..." 
	@g++ $(flags) -o $@ -c $^ 

objects:
	@mkdir $@
