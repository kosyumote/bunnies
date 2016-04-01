#include <random> // for stochastic sampling
#include <iostream> // for output
#include <chrono> // for seeding the random generator with the time
#include <cstring> // for low-level memory manipulations
#include <cmath> // for ln(x)
#include <png.h> // for making PNG output
#include <stdio.h> // for writing to files
#include <omp.h> // for getting thread id so different threads have different seeds

#define MAX_BUNNIES_PER_SQUARE 10
#define MAX_PANTHERS_PER_SQUARE 3
#define XSIZE 500 // 500 x 300
#define YSIZE 300 // are approximately the dimensions of PA, if each square has a size of 1 km x 1 km 
#define BUNNY_MIGRATION_PROBABILITY 0.006 // this is the probability that a bunny moves to a specific square... if this is 0.06, then the probability of moving is 0.48 and not moving is 0.52
#define PANTHER_MIGRATION_PROBABILITY 0.01
#define MAX_STEPS 997
#define INITIAL_BUNNY_POPULATION 100000
#define INITIAL_PANTHER_POPULATION 10000
#define HUNGRY_PANTHER_HUNTING_RATE 0.02
#define FED_PANTHER_HUNTING_RATE 0.01
#define HUNGRY_PANTHER_DEATH_RATE 0.01
#define FED_PANTHER_DEATH_RATE 0.003
#define TIME_STEP 1 // monthish

#define RANDOM_SEEDING 0
#define SEEDING_METHOD RANDOM_SEEDING

#define RANDOM_MOVEMENT 0
#define BUNNIES_RUN_AWAY 1
#define BUNNIES_NS_PANTHERS_EW 2 // not yet implemented
#define MIGRATION_METHOD BUNNIES_RUN_AWAY

#define BUNNIES_RUN_AWAY_DAMPENING 0.01 // a positive number -- as it goes to infinity, the model approaches RANDOM_MOVEMENT

#define OUTPUT_FILENAME "./pics/run.png"
#define FRAME_DELAY 5

void write_output(int** bunnies, int** panthers, int pic)
{
    png_byte **outputPNG = new png_byte*[YSIZE]; //[XSIZE][YSIZE];
    for(int i = 0; i < YSIZE; i++)
    {
	outputPNG[i] = new png_byte[3 * XSIZE];
    } 

    std::string fileName(OUTPUT_FILENAME);
    fileName = fileName.substr(0, fileName.length() - 4);
    std::string fileNameAppender = std::to_string(pic);
    fileNameAppender = fileNameAppender.insert(0, 3 - fileNameAppender.length(), '0');
    fileName += fileNameAppender;
    fileName += ".png";
	
    FILE *png = fopen(fileName.c_str(), "wb");

    png_structp png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, (png_voidp) NULL, NULL, NULL);

    png_infop info_ptr = png_create_info_struct(png_ptr);
    png_init_io(png_ptr, png);

    png_set_IHDR(png_ptr, info_ptr, XSIZE, YSIZE, 8, PNG_COLOR_TYPE_RGB, PNG_INTERLACE_NONE, PNG_COMPRESSION_TYPE_DEFAULT, PNG_FILTER_TYPE_DEFAULT);

    png_write_info(png_ptr, info_ptr);
	
    for(int i = 0; i < XSIZE; i++)
    {
	for(int j = 0; j < YSIZE; j++)
	{
	    int bunnyNumber = 256 * bunnies[i][j] / (MAX_BUNNIES_PER_SQUARE + 1);
	    int pantherNumber = 256 * panthers[i][j] / (MAX_PANTHERS_PER_SQUARE + 1); 
	    if(bunnyNumber > pantherNumber)
	    {
		outputPNG[j][3 * i] = 255;
		outputPNG[j][3 * i + 1] = 255 - bunnyNumber + pantherNumber;
		outputPNG[j][3 * i + 2] = 255 - bunnyNumber;
	    }
	    else
	    {
		outputPNG[j][3 * i] = 255 - pantherNumber + bunnyNumber;
		outputPNG[j][3 * i + 1] = 255;
		outputPNG[j][3 * i + 2] = 255 - pantherNumber;
	    }
	}
    }

    png_write_image(png_ptr, outputPNG);

    png_write_end(png_ptr, info_ptr);
    png_destroy_write_struct(&png_ptr, &info_ptr);
    fclose(png);

    for(int i = 0; i < YSIZE; i++)
    {
	delete[] outputPNG[i];
    }
    delete[] outputPNG;
}

int main()
{
    std::mt19937 rand(std::chrono::system_clock::now().time_since_epoch().count());

    std::chrono::high_resolution_clock::time_point startTime; // for timing how long one step takes
    std::chrono::high_resolution_clock::time_point endTime = std::chrono::high_resolution_clock::now();
    std::chrono::high_resolution_clock::time_point completeStartTime = std::chrono::high_resolution_clock::now();

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

    switch(SEEDING_METHOD)
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


    write_output(bunnies, panthers, 1);

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
	
///////// MIGRATION STEP (NONTEMP -> TEMP) ////////////// to be honest, this whole section should be in a switch (or #if) statement, so that separate migration methods can implement it as efficiently as possible
#if MIGRATION_METHOD == BUNNIES_RUN_AWAY
#pragma omp parallel for collapse(2)
	for(int i = 0; i < XSIZE; i++)
	{
	    for(int j = 0; j < YSIZE; j++)
	    {
		double migrationFactor = 0;
		for(int k = 0; k < 8; k++)
		{
		    migrationFactor += 1 / (panthers[(i + XSIZE + (1 - k / 4 * 2) * ((bool) (k % 4))) % XSIZE][(j + YSIZE + (1 - ((k + 2) % 8) / 4 * 2) * ((bool) ((k + 2) % 4))) % YSIZE] + BUNNIES_RUN_AWAY_DAMPENING);
		}
		
		for(int k = 0; k < 8; k++)
		{
		    double bunnyProbability = BUNNY_MIGRATION_PROBABILITY * 8 / (panthers[(i + XSIZE + (1 - k / 4 * 2) * ((bool) (k % 4))) % XSIZE][(j + YSIZE + (1 - ((k + 2) % 8) / 4 * 2) * ((bool) ((k + 2) % 4))) % YSIZE] + BUNNIES_RUN_AWAY_DAMPENING) / migrationFactor;
		    
		    std::binomial_distribution<int> binomialDistribution(bunnies[i][j], bunnyProbability);
		    int movers = binomialDistribution(rand);

#pragma omp atomic
		    tempBunnies[(i + XSIZE + (1 - k / 4 * 2) * ((bool) (k % 4))) % XSIZE][(j + YSIZE + (1 - ((k + 2) % 8) / 4 * 2) * ((bool) ((k + 2) % 4))) % YSIZE] += movers; // this line "magically" gets the right direction
		    bunnies[i][j] -= movers;
		}

#pragma omp atomic
		tempBunnies[i][j] += bunnies[i][j];
	    }
	}


#pragma omp parallel for collapse(2)
	for(int i = 0; i < XSIZE; i++)
	{
	    for(int j = 0; j < YSIZE; j++)
	    {
		for(int k = 0; k < 8; k++)
		{
		    std::binomial_distribution<int> binomialDistribution(panthers[i][j], PANTHER_MIGRATION_PROBABILITY);
		    int movers = binomialDistribution(rand);

#pragma omp atomic
		    tempPanthers[(i + XSIZE + (1 - k / 4 * 2) * ((bool) (k % 4))) % XSIZE][(j + YSIZE + (1 - ((k + 2) % 8) / 4 * 2) * ((bool) ((k + 2) % 4))) % YSIZE] += movers;
		    panthers[i][j] -= movers;
		}

#pragma omp atomic
		tempPanthers[i][j] += panthers[i][j];

	    }
	}

#else // this is the default, or when MIGRATION_METHOD == RANDOM_MOTION
#pragma omp parallel for collapse(2)
	for(int i = 0; i < XSIZE; i++)
	{
	    for(int j = 0; j < YSIZE; j++)
	    {
		for(int k = 0; k < 8; k++)
		{
		    std::binomial_distribution<int> binomialDistribution(bunnies[i][j], BUNNY_MIGRATION_PROBABILITY);
		    int movers = binomialDistribution(rand);

#pragma omp atomic
		    tempBunnies[(i + XSIZE + (1 - k / 4 * 2) * ((bool) (k % 4))) % XSIZE][(j + YSIZE + (1 - ((k + 2) % 8) / 4 * 2) * ((bool) ((k + 2) % 4))) % YSIZE] += movers; // this line "magically" gets the right direction
		    bunnies[i][j] -= movers;
		    
		    binomialDistribution = std::binomial_distribution<int>(panthers[i][j], PANTHER_MIGRATION_PROBABILITY);
		    movers = binomialDistribution(rand);

#pragma omp atomic
		    tempPanthers[(i + XSIZE + (1 - k / 4 * 2) * ((bool) (k % 4))) % XSIZE][(j + YSIZE + (1 - ((k + 2) % 8) / 4 * 2) * ((bool) ((k + 2) % 4))) % YSIZE] += movers;
		    panthers[i][j] -= movers;
		}

#pragma omp atomic
		tempBunnies[i][j] += bunnies[i][j];
#pragma omp atomic
		tempPanthers[i][j] += panthers[i][j];

	    }
	}
#endif

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
		    while(t < TIME_STEP && u + f > 0) // for a month, as long as panthers still exist -- otherwise bunnies remain constant
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

//////////////// STORE A FRAME ////////////////////////////////////////////

	write_output(bunnies, panthers, step + 2);

    }


    std::chrono::duration<double> simTime = std::chrono::duration_cast<std::chrono::duration<double>>(std::chrono::high_resolution_clock::now() - completeStartTime);
    std::cout << "\nTotal simulation time: " << simTime.count() << " seconds." << std::endl;

//////////////// SAVE THE ANIMATED PNG USING APNGASM ////////////////////////////////////////////
    
    std::cout << "Processing files..." << std::endl;
    system("utils/apngasm ./pics/run.apng ./pics/run001.png 1 20 > log.txt");
    system("rm ./pics/run*.png");


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
