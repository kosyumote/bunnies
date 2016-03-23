#include <random> // for stochastic sampling
#include <iostream> // for output
#include <chrono> // for seeding the random generator with the time
#include <cstring> // for low-level memory manipulations
#include <cmath> // for ln(x)
#include <gif_lib.h> // for making GIF output

#define MAX_BUNNIES_PER_SQUARE 1000
#define MAX_PANTHERS_PER_SQUARE 100
#define XSIZE 500
#define YSIZE 300 // these are approximately the dimensions of PA, if each square has a size of 1 km x 1 km 
#define BUNNY_MIGRATION_PROBABILITY 0.06 // this is the probability that a bunny moves to a specific square... if this is 0.06, then the probability of moving is 0.48 and not moving is 0.52
#define PANTHER_MIGRATION_PROBABILITY 0.1
#define MAX_STEPS 100
#define INITIAL_BUNNY_POPULATION 100000
#define INITIAL_PANTHER_POPULATION 1000
#define HUNGRY_PANTHER_HUNTING_RATE 0.2
#define FED_PANTHER_HUNTING_RATE 0.1
#define HUNGRY_PANTHER_DEATH_RATE 0.1
#define FED_PANTHER_DEATH_RATE 0.03

int main()
{
    
    std::mt19937 rand(std::chrono::system_clock::now().time_since_epoch().count()); // random number generator
    std::chrono::high_resolution_clock::time_point startTime; // for timing how long one step takes
    std::chrono::high_resolution_clock::time_point endTime = std::chrono::high_resolution_clock::now();

    ColorMapObject* outputPalette = MakeMapObject(256, NULL); // palette for making GIF output
    

    //////////////////////////
    // VARIABLE DEFINITIONS //
    //////////////////////////

    int bunnyPopulation = 0;
    int pantherPopulation = 0;

    //////////////////////////////////////////////
    // ALLOCATE MEMORY FOR BUNNY/PANTHER ARRAYS //
    //////////////////////////////////////////////

    int** bunnies = new int*[XSIZE]; // in consideration of memory, maybe make these a smaller int type; however, they still have to be able to hold at least MAX_BUNNIES_PER_SQUARE * 9
    int** panthers = new int*[XSIZE];
    int** tempBunnies = new int*[XSIZE];
    int** tempPanthers = new int*[XSIZE];

    for(int i = 0; i < XSIZE; i++)
    {
	bunnies[i] = new int[YSIZE];
	panthers[i] = new int[YSIZE];
	tempBunnies[i] = new int[YSIZE];
	tempPanthers[i] = new int[YSIZE];
    }


    //////////////////////////////////////
    // SET INITIAL CONDITIONS SOMEHOW?? //
    //////////////////////////////////////

    int initializationMethod = 0; // change this to pick the initialization method

    switch(initializationMethod)
    {
    case 0:
    default:
    {
	//////////////////////////////////////////////////////////////////////////
	// This initialization method only approximates the initial populations //
	// but is quite fast. It usually underestimates the populations.        //
	//////////////////////////////////////////////////////////////////////////

	std::binomial_distribution<int> bunnyDistribution(MAX_BUNNIES_PER_SQUARE, (double) INITIAL_BUNNY_POPULATION / XSIZE / YSIZE / MAX_BUNNIES_PER_SQUARE);
	std::binomial_distribution<int> pantherDistribution(MAX_PANTHERS_PER_SQUARE, (double) INITIAL_PANTHER_POPULATION / XSIZE / YSIZE / MAX_PANTHERS_PER_SQUARE);
	for(int i = 0; i < XSIZE; i++)
	{
	    for(int j = 0; j < YSIZE; j++)
	    {
		bunnies[i][j] = bunnyDistribution(rand);
		bunnyPopulation += bunnies[i][j];
		panthers[i][j] = pantherDistribution(rand);
		pantherPopulation += panthers[i][j];
		
		tempBunnies[i][j] = 0;
		tempPanthers[i][j] = 0;
	    }
	}
    }
    }


    ////////////////////////
    // Run the simulation //
    ////////////////////////

    std::uniform_real_distribution<double> realDistribution(0.0, 1.0);

    for(int step = 0; step < MAX_STEPS; step++)
    {
	startTime = endTime;
	endTime = std::chrono::high_resolution_clock::now();
	
	std::chrono::duration<double> timeSpan = std::chrono::duration_cast<std::chrono::duration<double>>(endTime - startTime);

	std::cout << "Step " << step + 1 << " (" << timeSpan.count() << " seconds):\n    bunny population: " << bunnyPopulation << "\n    panther population: " << pantherPopulation << std::endl; 
	
///////// MIGRATION STEP (NONTEMP -> TEMP) //////////////

#pragma omp parallel for collapse(2)
	for(int i = 0; i < XSIZE; i++)
	{
	    for(int j = 0; j < YSIZE; j++)
	    {
		for(int k = 0; k < bunnies[i][j]; k++)
		{
		    double r = realDistribution(rand);
		    if((r -= BUNNY_MIGRATION_PROBABILITY) < 0) // up
		    {
#pragma omp atomic
			tempBunnies[i][(j + YSIZE - 1) % YSIZE]++;
		    }
		    else if((r -= BUNNY_MIGRATION_PROBABILITY) < 0) // up-right
		    {
#pragma omp atomic
			tempBunnies[(i + 1) % XSIZE][(j + YSIZE - 1) % YSIZE]++;
		    }
		    else if((r -= BUNNY_MIGRATION_PROBABILITY) < 0) // right
		    {
#pragma omp atomic
			tempBunnies[(i + 1) % XSIZE][j]++;
		    }
		    else if((r -= BUNNY_MIGRATION_PROBABILITY) < 0) // down-right
		    {
#pragma omp atomic
			tempBunnies[(i + 1) % XSIZE][(j + 1) % YSIZE]++;
		    }
		    else if((r -= BUNNY_MIGRATION_PROBABILITY) < 0) // down
		    {
#pragma omp atomic
			tempBunnies[i][(j + 1) % YSIZE]++;
		    }
		    else if((r -= BUNNY_MIGRATION_PROBABILITY) < 0) // down-left
		    {
#pragma omp atomic
			tempBunnies[(i + XSIZE - 1) % XSIZE][(j + 1) % YSIZE]++;
		    }
		    else if((r -= BUNNY_MIGRATION_PROBABILITY) < 0) // left
		    {
#pragma omp atomic
			tempBunnies[(i + XSIZE - 1) % XSIZE][j]++;
		    }
		    else if((r -= BUNNY_MIGRATION_PROBABILITY) < 0) // up-left
		    {
#pragma omp atomic
			tempBunnies[(i + XSIZE - 1) % XSIZE][(j + YSIZE - 1) % YSIZE]++;
		    }
		    else
		    {
			tempBunnies[i][j]++;
		    }
		}

		for(int k = 0; k < panthers[i][j]; k++)
		{
		    double r = realDistribution(rand);
		    if((r -= PANTHER_MIGRATION_PROBABILITY) < 0) // up
		    {
#pragma omp atomic
			tempPanthers[i][(j + YSIZE - 1) % YSIZE]++;
		    }
		    else if((r -= PANTHER_MIGRATION_PROBABILITY) < 0) // up-right
		    {
#pragma omp atomic
			tempPanthers[(i + 1) % XSIZE][(j + YSIZE - 1) % YSIZE]++;
		    }
		    else if((r -= PANTHER_MIGRATION_PROBABILITY) < 0) // right
		    {
#pragma omp atomic
			tempPanthers[(i + 1) % XSIZE][j]++;
		    }
		    else if((r -= PANTHER_MIGRATION_PROBABILITY) < 0) // down-right
		    {
#pragma omp atomic
			tempPanthers[(i + 1) % XSIZE][(j + 1) % YSIZE]++;
		    }
		    else if((r -= PANTHER_MIGRATION_PROBABILITY) < 0) // down
		    {
#pragma omp atomic
			tempPanthers[i][(j + 1) % YSIZE]++;
		    }
		    else if((r -= PANTHER_MIGRATION_PROBABILITY) < 0) // down-left
		    {
#pragma omp atomic
			tempPanthers[(i + XSIZE - 1) % XSIZE][(j + 1) % YSIZE]++;
		    }
		    else if((r -= PANTHER_MIGRATION_PROBABILITY) < 0) // left
		    {
#pragma omp atomic
			tempPanthers[(i + XSIZE - 1) % XSIZE][j]++;
		    }
		    else if((r -= PANTHER_MIGRATION_PROBABILITY) < 0) // up-left
		    {
#pragma omp atomic
			tempPanthers[(i + XSIZE - 1) % XSIZE][(j + YSIZE - 1) % YSIZE]++;
		    }
		    else // no migration
		    {
			tempPanthers[i][j]++;
		    }
		}
	    }
	}



////////////////// RUN A GILLESPIE SIMULATION ON EACH SQUARE AND UPDATE THE NUMBERS OF BUNNIES/PANTHERS (TEMP -> NONTEMP) //////////////////////

	bunnyPopulation = 0;
	pantherPopulation = 0;

#pragma omp parallel for collapse(2)
	for(int i = 0; i < XSIZE; i++)
	{
	    for(int j = 0; j < YSIZE; j++)
	    {
		double t = 0; // simulation time
		int u = tempPanthers[i][j]; // unfed panthers
		int f = 0; // fed panthers
		int b = tempBunnies[i][j]; // number of bunnies

		if(u != 0)
		{
		    while(t < 1 && u + f > 0) // for a month, as long as panthers still exist -- otherwise bunnies remain constant
		    {
			double tot = HUNGRY_PANTHER_HUNTING_RATE * u * b + FED_PANTHER_HUNTING_RATE * f * b + HUNGRY_PANTHER_DEATH_RATE * u + FED_PANTHER_DEATH_RATE * f;
			double r = tot * realDistribution(rand);
			if((r -= HUNGRY_PANTHER_HUNTING_RATE * u * b) < 0)
			{
			    b--;
			    u--;
			    f++;
			}
			else if((r -= FED_PANTHER_HUNTING_RATE * f * b) < 0)
			{
			    b--;
			}
			else if((r -= HUNGRY_PANTHER_DEATH_RATE * u) < 0)
			{
			    u--;
			}
			else // if((r -= FED_PANTHER_DEATH_RATE * f) < 0) -- should always be true; no need to check
			{
			    f--;
			}
			t -= log(realDistribution(rand)) / tot;

			
		    }
		}

		bunnies[i][j] = std::min(2 * b, MAX_BUNNIES_PER_SQUARE);
		panthers[i][j] = std::min(u + 2 * f, MAX_PANTHERS_PER_SQUARE);

#pragma omp atomic
		bunnyPopulation += bunnies[i][j];
#pragma omp atomic
		pantherPopulation += panthers[i][j];
	    }
	}

////////////////// CLEAR TEMP ARRAYS ////////////////////////////////////////////

	for(int i = 0; i < XSIZE; i++)
	{
	    memset(tempBunnies[i], 0, YSIZE * sizeof(int));
	    memset(tempPanthers[i], 0, YSIZE * sizeof(int));
	}

//////////////// DRAW THE PICTURE ////////////////////////////////////////////

	std::vector<GifByteType> r, g, b;
	
			    

    }



    ////////////////////////////////////////////////
    // DEALLOCATE MEMORY FOR BUNNY/PANTHER ARRAYS //
    ////////////////////////////////////////////////

    for(int i = 0; i < XSIZE; i++)
    {
	delete[] tempBunnies[i];
	delete[] tempPanthers[i];
	delete[] bunnies[i];
	delete[] panthers[i];
    }

    delete[] tempBunnies;
    delete[] tempPanthers;
    delete[] bunnies;
    delete[] panthers;
    
}
