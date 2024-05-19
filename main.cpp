#define _CRT_SECURE_NO_WARNINGS // for stb

#include <stdio.h>
#include <cmath>
#include <random>
#include <vector>
#include <array>
#include <direct.h>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb/stb_image_write.h"

typedef unsigned int uint;

typedef std::array<uint, 2> uint2;

typedef std::array<float, 2> float2;
typedef std::array<float, 3> float3;

// -1 means non deterministic
#define SEED() 0

#define DO_TEST_COPRIMES() false
// When looking for coprimes near the golden ratio integer, this is the range looked at
static const uint c_coprimeTestStart = 2;
static const uint c_coprimeTestEnd = 1000;

#define DO_TEST() true
#define PRINT_NUMBERS() false
static const uint c_itemCounts[] = { 10, 37, 64, 100 , 1000, 1337, 1546 };

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

// From https://blog.demofox.org/2015/01/24/programmatically-calculating-gcd-and-lcm/
uint CalculateGCD(uint smaller, uint larger)
{
	// make sure A <= B before starting
	if (larger < smaller)
		std::swap(smaller, larger);

	// loop
	while (1)
	{
		// if the remainder of larger / smaller is 0, they are the same
		// so return smaller as the GCD
		uint remainder = larger % smaller;
		if (remainder == 0)
			return smaller;

		// otherwise, the new larger number is the old smaller number, and
		// the new smaller number is the remainder
		larger = smaller;
		smaller = remainder;
	}
}

// From https://blog.demofox.org/2015/09/10/modular-multiplicative-inverse/
unsigned int ExtendedEuclidianAlgorithm(int smaller, int larger, int& s, int& t)
{
	// make sure A <= B before starting
	bool swapped = false;
	if (larger < smaller)
	{
		swapped = true;
		std::swap(smaller, larger);
	}

	// set up our storage for the loop.  We only need the last two values so will
	// just use a 2 entry circular buffer for each data item
	std::array<int, 2> remainders = { larger, smaller };
	std::array<int, 2> ss = { 1, 0 };
	std::array<int, 2> ts = { 0, 1 };
	int indexNeg2 = 0;
	int indexNeg1 = 1;

	// loop
	while (1)
	{
		// calculate our new quotient and remainder
		int newQuotient = remainders[indexNeg2] / remainders[indexNeg1];
		int newRemainder = remainders[indexNeg2] - newQuotient * remainders[indexNeg1];

		// if our remainder is zero we are done.
		if (newRemainder == 0)
		{
			// return our s and t values as well as the quotient as the GCD
			s = ss[indexNeg1];
			t = ts[indexNeg1];
			if (swapped)
				std::swap(s, t);
			return remainders[indexNeg1];
		}

		// calculate this round's s and t
		int newS = ss[indexNeg2] - newQuotient * ss[indexNeg1];
		int newT = ts[indexNeg2] - newQuotient * ts[indexNeg1];

		// store our values for the next iteration
		remainders[indexNeg2] = newRemainder;
		ss[indexNeg2] = newS;
		ts[indexNeg2] = newT;

		// move to the next iteration
		std::swap(indexNeg1, indexNeg2);
	}
}

bool IsCoprime(uint A, uint B)
{
	return CalculateGCD(A, B) == 1;
}

// Find the coprime of n that is nearest to target
uint GetNearestCoprime(uint target, uint n)
{
	uint coprime = 0;

	uint offset = 0;
	while (1)
	{
		if (offset < target)
		{
			coprime = target - offset;
			if (IsCoprime(coprime, n))
				break;
		}

		coprime = target + offset + 1;
		if (coprime < n && IsCoprime(coprime, n))
			break;

		offset++;
	}

	return coprime;
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
		uint coprime = GetNearestCoprime(uint(c_goldenRatioFract * float(itemCount) + 0.5f), itemCount);

		printf("GR INT (%u vs %0.2f):", coprime, c_goldenRatioFract * float(itemCount));

		if (PRINT_NUMBERS())
			printf("\n");
		else
			printf(" ");

		uint valueI = Clamp<uint>((uint)(startValue * float(itemCount)), 0, itemCount - 1);

		DoTest_Single(itemCount,
			[&valueI, itemCount, coprime]()
			{
				uint ret = valueI;
				valueI = (valueI + coprime) % itemCount;
				return ret;
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
			uint coprime = GetNearestCoprime(target, itemCount);

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

// TODO: need to make a fully featured 1d shuffle iterator with random access and inversion
// TODO: show average of shuffle sequence, vs white noise, to show quality
