#define _CRT_SECURE_NO_WARNINGS // for stb

#include <stdio.h>
#include <cmath>
#include <random>
#include <vector>
#include <array>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb/stb_image_write.h"

typedef unsigned int uint;

typedef std::array<float, 2> float2;
typedef std::array<uint, 2> uint2;

// -1 means non deterministic
#define SEED() 0

#define DO_TEST_COPRIMES() false
// When looking for coprimes near the golden ratio integer, this is the range looked at
static const uint c_coprimeTestStart = 2;
static const uint c_coprimeTestEnd = 1000;

#define DO_TEST_1D() false
#define PRINT_NUMBERS_1D() true
static const uint c_itemCounts1D[] = { 10, 37, 64, 100 };// , 1000, 1337, 1546 };

#define DO_TEST_2D() true
#define PRINT_NUMBERS_2D() true
static const uint2 c_imageSizes2D[] = { {3, 3} };// { 64, 64 }, { 640, 480 } };

template <typename T>
T Frac(T f)
{
	return std::fmod(f, T(1.0));
}

static const float c_goldenRatio = (1.0f + std::sqrt(5.0f)) / 2.0f;
static const float c_goldenRatioFract = Frac(c_goldenRatio); // Same value as 1.0f / golden ratio aka golden ratio conjugate

// A generalization of the golden ratio to 2D, from http://extremelearning.com.au/unreasonable-effectiveness-of-quasirandom-sequences/
static const float c_R2_g = 1.32471795724474602596f;
static const float c_R2_a1 = 1 / c_R2_g;
static const float c_R2_a2 = 1 / (c_R2_g * c_R2_g);

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

// R2 is from http://extremelearning.com.au/unreasonable-effectiveness-of-quasirandom-sequences/
// R2 Low discrepancy sequence
float2 R2LDS(int index)
{
	return float2{ Frac(float(index) * c_R2_a1), Frac(float(index) * c_R2_a2) };
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

// Find the coprime of m and n that is nearest to target
uint GetNearestCoprime(uint target, uint m, uint n)
{
	uint coprime = 0;

	uint offset = 0;
	while (1)
	{
		if (offset < target)
		{
			coprime = target - offset;
			if (IsCoprime(coprime, m) && IsCoprime(coprime, n))
				break;
		}

		coprime = target + offset + 1;
		if (coprime < n && IsCoprime(coprime, m) && IsCoprime(coprime, n))
			break;

		offset++;
	}

	return coprime;
}

template <typename LAMBDA>
void DoTest1D_Single(uint itemCount, const LAMBDA& lambda)
{
	std::vector<int> visitCount(itemCount, 0);
	int errorCount = 0;
	for (uint index = 0; index < itemCount; ++index)
	{
		uint shuffleItem = lambda();

		if(PRINT_NUMBERS_1D())
			printf("%u ", shuffleItem);
		visitCount[shuffleItem]++;
		if (visitCount[shuffleItem] > 1)
			errorCount++;
	}

	if (PRINT_NUMBERS_1D())
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

void DoTest1D(uint itemCount, std::mt19937& rng)
{
	printf("========== 1D: %i items ==========\n\n", itemCount);

	std::uniform_real_distribution<float> dist(0.0f, 1.0f);
	float startValue = dist(rng);

	// Golden Ratio LDS
	{
		printf("GR LDS:");

		if (PRINT_NUMBERS_1D())
			printf("\n");
		else
			printf(" ");

		float valueF = startValue;
		DoTest1D_Single(itemCount,
			[&valueF, itemCount] ()
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

		if (PRINT_NUMBERS_1D())
			printf("\n");
		else
			printf(" ");

		uint valueI = Clamp<uint>((uint)(startValue * float(itemCount)), 0, itemCount - 1);

		DoTest1D_Single(itemCount,
			[&valueI, itemCount, coprime]()
			{
				uint ret = valueI;
				valueI = (valueI + coprime) % itemCount;
				return ret;
			}
		);
	}
}

template <typename LAMBDA>
void DoTest2D_Single(const uint2& dims, const LAMBDA& lambda)
{
	std::vector<int> visitCount(dims[0]*dims[1], 0);

	int errorCount = 0;
	for (uint index = 0; index < dims[0] * dims[1]; ++index)
	{
		uint2 shuffleItem = lambda(index);

		if (PRINT_NUMBERS_2D())
			printf("(%u,%u) ", shuffleItem[0], shuffleItem[1]);

		uint shuffleItemFlat = shuffleItem[1] * dims[0] + shuffleItem[0];

		visitCount[shuffleItemFlat]++;
		if (visitCount[shuffleItemFlat] > 1)
			errorCount++;
	}

	if (PRINT_NUMBERS_2D())
	{
		printf("\n");
		if (errorCount > 0)
		{
			printf("Duplicates:");
			for (uint index = 0; index < dims[0] * dims[1]; ++index)
			{
				if (visitCount[index] > 1)
					printf(" (%u,%u)", index / dims[0], index % dims[0]);
			}
			printf("\n");

			printf("Missing:");
			for (uint index = 0; index < dims[0] * dims[1]; ++index)
			{
				if (visitCount[index] == 0)
					printf(" (%u,%u)", index / dims[0], index % dims[0]);
			}
			printf("\n");
		}
	}

	printf("%i Errors (%0.2f%%)\n\n", errorCount, 100.0f * float(errorCount) / float(dims[0] * dims[1]));
}

void DoTest2D(const uint2& dims, std::mt19937& rng)
{
	printf("========== 2D: %u x %u ==========\n\n", dims[0], dims[1]);

	std::uniform_real_distribution<float> dist(0.0f, 1.0f);
	float2 startValue = { dist(rng), dist(rng) };

	// R2 LDS
	{
		printf("R2 LDS:");

		if (PRINT_NUMBERS_2D())
			printf("\n");
		else
			printf(" ");

		float2 valueF = startValue;
		DoTest2D_Single(dims,
			[&valueF, &dims](uint index)
			{
				uint x = index % dims[0];
				uint y = index / dims[0];

				uint2 ret;
				ret[0] = Clamp<uint>((uint)(valueF[0] * float(dims[0])), 0, dims[0] - 1);
				ret[1] = Clamp<uint>((uint)(valueF[1] * float(dims[1])), 0, dims[1] - 1);

				valueF[0] = Frac(valueF[0] + c_R2_a1);
				valueF[1] = Frac(valueF[1] + c_R2_a2);

				return ret;
			}
		);
	}

	// Integer solution
	{
		uint2 coprimes;
		coprimes[0] = GetNearestCoprime(uint(c_R2_a1 * float(dims[0]) + 0.5f), dims[0]);
		coprimes[1] = GetNearestCoprime(uint(c_R2_a2 * float(dims[1]) + 0.5f), dims[1], coprimes[0]);

		printf("R2 Int (%u, %u vs %0.2f, %0.2f):", coprimes[0], coprimes[1], c_R2_a1 * float(dims[0]), c_R2_a2 * float(dims[1]));

		if (PRINT_NUMBERS_2D())
			printf("\n");
		else
			printf(" ");

		uint2 valueI;
		valueI[0] = Clamp<uint>((uint)(startValue[0] * float(dims[0])), 0, dims[0] - 1);
		valueI[1] = Clamp<uint>((uint)(startValue[1] * float(dims[1])), 0, dims[1] - 1);

		DoTest2D_Single(dims,
			[&valueI, &dims, coprimes](uint index)
			{
				uint x = index % dims[0];
				uint y = index / dims[0];

				uint2 ret = valueI;
				valueI[0] = (valueI[0] + coprimes[0]) % dims[0];
				valueI[1] = (valueI[1] + coprimes[1]) % dims[1];
				return ret;
			}
		);
	}
}

int main(int argc, char** argv)
{
	// initialize RNG
	std::random_device rd;
	uint seed = (SEED() == -1) ? rd() : SEED();
	printf("Seed = %u\n\n", seed);

	// Coprime test
	if (DO_TEST_COPRIMES())
	{
		FILE* file = nullptr;
		fopen_s(&file, "coprimes.csv", "wb");

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

	// 1D Tests
	if (DO_TEST_1D())
	{
		for (uint itemCount : c_itemCounts1D)
		{
			std::mt19937 rng(seed);
			DoTest1D(itemCount, rng);
		}
	}

	// 2D Tests
	if (DO_TEST_2D())
	{
		std::mt19937 rng(seed);
		for (const uint2& size : c_imageSizes2D)
			DoTest2D(size, rng);
	}

	return 0;
}

/*
TODO:
- work out 2D stuff? should we show some kind of image demo for it? like a digital disolve?

Blog:
- talk about how to extend it to 2d with R2?
 - need to think through how to make this work without repeats though. maybe co-irrational to each other, as well as the edge lengths? idk.
- GCD: https://blog.demofox.org/2015/01/24/programmatically-calculating-gcd-and-lcm/
- the coprime version of the golden ratio is real close to the actual thing. but drift.
- show where the sequnces diverge! lil bit of drift due to rounding to integer, but we want that cause it's a coprime generator.
- mention that you are not going to be analzying the quality of the sequence, but you can see it's very, very similar to the pure golden ratio one

find marc reynolds post about golden ratio integer LDS thing.
and link this too: ttps://mastodon.gamedev.place/@demofox/112452866510781822

*/