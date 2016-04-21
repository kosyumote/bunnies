.PHONY: run iv div all

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

div: all
	@echo "Running bunnies divergence simulation:"
	@exec/divbunnies

all: exec/bunnies exec/ivbunnies exec/divbunnies

exec/bunnies: objects/bunnies.o | exec
	@echo "Creating bunnies executable..."
	@g++ $(flags) -o $@ $^ -lpng -lz

objects/bunnies.o: bunnies.cpp | objects
	@echo "Compiling bunnies.cpp..." 
	@g++ $(flags) -o $@ -c $^ 

exec/divbunnies: objects/divbunnies.o | exec
	@echo "Creating divbunnies executable..."
	@g++ $(flags) -o $@ $^ -lpng -lz

objects/divbunnies.o: bunnies_diverge.cpp | objects
	@echo "Compiling bunnies_diverge.cpp..." 
	@g++ $(flags) -o $@ -c $^ 

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


