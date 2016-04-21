.PHONY: run iv all

debug-flags = -Wall -g
release-flags = -Ofast

other-flags = -std=gnu++11 -fopenmp
flags = $(release-flags) $(other-flags)

run: all
	@echo "Running bunnies simulation:"
	@exec/bunnies

iv: all
	@echo "Running bunnies intial values simulation:"
	@exec/ivbunnies

all: exec/bunnies exec/ivbunnies

exec/bunnies: objects/bunnies.o | exec
	@echo "Creating bunnies executable..."
	@g++ $(flags) -o $@ $^ -lpng -lz

objects/bunnies.o: bunnies.cpp | objects
	@echo "Compiling bunnies.cpp..." 

exec/ivbunnies: objects/ivbunnies.o | exec
	@echo "Creating ivbunnies executable..."
	@g++ $(flags) -o $@ $^ -lpng -lz

exec:
	@mkdir $@

objects/ivbunnies.o: bunnies_initial_value.cpp | objects
	@echo "Compiling bunnies_initial_value.cpp..." 
	@g++ $(flags) -o $@ -c $^ 

objects:
	@mkdir $@


