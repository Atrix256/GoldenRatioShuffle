#define _CRT_SECURE_NO_WARNINGS // for stb

#include <stdio.h>
#include <cmath>
#include <random>
#include <vector>
#include <direct.h>

#include "LDShuffle.h"

typedef unsigned int uint;

// -1 means random seed
// Any other value is used as a deterministic seed
#define SEED() -1

#define DO_TEST_COPRIMES() false
static const uint c_coprimeTestStart = 2;
static const uint c_coprimeTestEnd = 1000;

#define DO_TEST_INVERSION() false
static const uint c_inversionTestItemCount = 65536;

#define DO_TEST_CONVERGENCE() false
static const uint c_convergenceTestItemCount = 10000;
static const uint c_convergenceTestTestCount = 1000;

#define DO_TEST() true
#define PRINT_NUMBERS() true
static const uint c_itemCounts[] = { 16 };// , 64, 37, 64, 100, 1000, 1337, 1546 };

template <typename T>
T Frac(T f)
{
	return std::fmod(f, T(1.0));
}

static const float c_goldenRatio = (1.0f + std::sqrt(5.0f)) / 2.0f;
static const float c_goldenRatioFract = Frac(c_goldenRatio); // Same value as 1.0f / golden ratio aka golden ratio conjugate

template <typename T>
T Clamp(T value, T themin, T themax)
{
	if (value <= themin)
		return themin;
	else if (value >= themax)
		return themax;
	else
		return value;
}

float Lerp(float A, float B, float t)
{
	return A * (1.0f - t) + B * t;
}

static inline uint ExtendedEuclidianAlgorithm(int smaller, int larger, int& s, int& t)
{
	return 0;
}

template <typename LAMBDA>
void DoTest_Single(uint itemCount, const LAMBDA& lambda)
{
	std::vector<int> visitCount(itemCount, 0);
	int errorCount = 0;
	for (uint index = 0; index < itemCount; ++index)
	{
		uint shuffleItem = lambda();

		if(PRINT_NUMBERS())
			printf("%u ", shuffleItem);
		visitCount[shuffleItem]++;
		if (visitCount[shuffleItem] > 1)
			errorCount++;
	}

	if (PRINT_NUMBERS())
	{
		printf("\n");
		if (errorCount > 0)
		{
			printf("Duplicates:");
			for (uint index = 0; index < itemCount; ++index)
			{
				if (visitCount[index] > 1)
					printf(" %u", index);
			}
			printf("\n");

			printf("Missing:");
			for (uint index = 0; index < itemCount; ++index)
			{
				if (visitCount[index] == 0)
					printf(" %u", index);
			}
			printf("\n");
		}
	}

	printf("%i Errors (%0.2f%%)\n\n", errorCount, 100.0f * float(errorCount) / float(itemCount));
}

void DoTest(uint itemCount, std::mt19937& rng)
{
	std::uniform_real_distribution<float> dist(0.0f, 1.0f);
	float startValue = dist(rng);

	printf("========== %i items (%0.2f startValue) ==========\n\n", itemCount, startValue);

	// Golden Ratio LDS
	{
		printf("GR LDS:");

		if (PRINT_NUMBERS())
			printf("\n");
		else
			printf(" ");

		float valueF = startValue;
		DoTest_Single(itemCount,
			[&valueF, itemCount]()
			{
				uint ret = Clamp<uint>((uint)(valueF * float(itemCount)), 0, itemCount - 1);
				valueF = Frac(valueF + c_goldenRatioFract);
				return ret;
			}
		);
	}

	// Integer solution
	{
		uint shuffleSeed = Clamp<uint>((uint)(startValue * float(itemCount)), 0, itemCount - 1);

		LDShuffle shuffle(itemCount, shuffleSeed);

		printf("GR INT (%u vs %0.2f):", shuffle.GetCoprime(), c_goldenRatioFract * float(itemCount));

		if (PRINT_NUMBERS())
			printf("\n");
		else
			printf(" ");

		DoTest_Single(itemCount,
			[&shuffle]()
			{
				return shuffle.Next();
			}
		);
	}
}

int main(int argc, char** argv)
{
	_mkdir("out");

	// initialize RNG
	std::random_device rd;
	uint seed = (SEED() == -1) ? rd() : SEED();
	printf("Seed = %u\n\n", seed);

	// Coprime test
	if (DO_TEST_COPRIMES())
	{
		FILE* file = nullptr;
		fopen_s(&file, "out/coprimes.csv", "wb");

		fprintf(file, "\"itemCount\",\"goldenIndexF\",\"goldenIndexI\",\"coprime\",\"diffF\",\"diffI\"\n");

		uint minDiff = 0;
		uint maxDiff = 0;

		for (uint itemCount = c_coprimeTestStart; itemCount < c_coprimeTestEnd; ++itemCount)
		{
			uint target = uint(c_goldenRatioFract * float(itemCount) + 0.5f);

			LDShuffle shuffle(itemCount, 0);
			uint coprime = shuffle.GetCoprime();

			uint diff = (target >= coprime) ? target - coprime : coprime - target;

			if (itemCount == c_coprimeTestStart)
			{
				minDiff = diff;
				maxDiff = diff;
			}
			else
			{
				minDiff = std::min(minDiff, diff);
				maxDiff = std::max(maxDiff, diff);
			}

			float diffF = std::abs(c_goldenRatioFract * float(itemCount) - float(coprime));

			fprintf(file, "\"%u\",\"%f\",\"%u\",\"%u\",\"%f\",\"%u\"\n", itemCount, c_goldenRatioFract * float(itemCount), uint(c_goldenRatioFract * float(itemCount) + 0.5f), coprime, diffF, diff);
		}

		printf(
			"When looking for coprimes near the golden ratio index, the diffs were:\n"
			"  min: %u\n"
			"  max: %u\n"
			"\n",
			minDiff, maxDiff
		);

		fclose(file);
	}

	// Inversion test
	if (DO_TEST_INVERSION())
	{
		std::mt19937 rng(seed);
		std::uniform_int_distribution<uint> dist(0, c_inversionTestItemCount - 1);
		LDShuffle shuffle(c_inversionTestItemCount, dist(rng));
		for (uint i = 0; i < c_inversionTestItemCount; ++i)
		{
			uint valueNext = shuffle.Next();
			uint value = shuffle.RandomAccess(i);

			if (value != valueNext)
			{
				printf("Inversion Test Error: iterative vs random access values mismatched: [%u] %u vs %u.\n", i, valueNext, value);
				return 1;
			}

			uint index = shuffle.GetIndex(value);
			if (index != i)
			{
				printf("Inversion Test Error: [%u] Inversion failed.\n", i);
				return 1;
			}

			//printf("[%u] %u\n", index, value);
		}
		printf("Inversion test passed!\n\n");
	}

	// Convergence test
	if (DO_TEST_CONVERGENCE())
	{
		std::mt19937 rng(seed);
		std::uniform_int_distribution<uint> dist(0, c_convergenceTestItemCount - 1);

		std::vector<uint> shuffled(c_convergenceTestItemCount);
		for (uint i = 0; i < c_convergenceTestItemCount; ++i)
			shuffled[i] = i;

		// Do multiple tests and average the results
		std::vector<float> convergenceLDS(c_convergenceTestItemCount, 0.0f);
		std::vector<float> convergenceWhite(c_convergenceTestItemCount, 0.0f);
		for (uint testIndex = 0; testIndex < c_convergenceTestTestCount; ++testIndex)
		{
			// Init the low discrepancy shuffler with a random starting seed
			LDShuffle shuffle(c_convergenceTestItemCount, dist(rng));

			// Make the white noise shuffled array
			std::shuffle(shuffled.begin(), shuffled.end(), rng);

			float convergedValueLDS = 0.0f;
			float convergedValueWhite = 0.0f;
			for (uint i = 0; i < c_convergenceTestItemCount; ++i)
			{
				float valueLDS = float(shuffle.Next()) / float(c_convergenceTestItemCount - 1) - 0.5f;
				float valueWhite = float(shuffled[i]) / float(c_convergenceTestItemCount - 1) - 0.5f;

				convergedValueLDS = Lerp(convergedValueLDS, valueLDS, 1.0f / float(i + 1));
				convergedValueWhite = Lerp(convergedValueWhite, valueWhite, 1.0f / float(i + 1));

				convergenceLDS[i] = Lerp(convergenceLDS[i], convergedValueLDS, 1.0f / float(testIndex + 1));
				convergenceWhite[i] = Lerp(convergenceWhite[i], convergedValueWhite, 1.0f / float(testIndex + 1));
			}
		}

		// Write the results out
		{
			FILE* file = nullptr;
			fopen_s(&file, "out/convergence.csv", "wb");
			fprintf(file, "\"itemCount\",\"LDShuffle\",\"Shuffle\"\n");

			for (uint i = 0; i < c_convergenceTestItemCount; ++i)
				fprintf(file, "\"%u\",\"%f\",\"%f\"\n", i, std::abs(convergenceLDS[i]), std::abs(convergenceWhite[i]));

			fclose(file);
		}
	}

	// Shuffle Test
	if (DO_TEST())
	{
		for (uint itemCount : c_itemCounts)
		{
			std::mt19937 rng(seed);
			DoTest(itemCount, rng);
		}
	}

	return 0;
}
